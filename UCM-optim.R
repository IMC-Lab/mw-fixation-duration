##################################################################################################
##                 Log-Likelihood Functions
##################################################################################################

#' Turn a vector into a discrete distribution of binned values
#'
#' @param x the vector of values to discretize
#' @param min the minimum value of allowed values in x
#' @param max the maximum value of allowed values in x
#' @param binwidth the width of the distribution's bins
#' @return a dataframe with the following columns, where each row corresponds to a bin:
#'  - lower: the lower bound of the bin
#'  - upper: the upper bound of the bin
#'  - mid: the midpoint of the bin
#'  - count: the number of data points in this bin
#'  - p: the probability of finding a data point in this bin
discretize <- function(x, min=0, max=1.2, binwidth=.06) {
    breaks <- seq(min, max, by=binwidth)
    x <- x[x >= min & x <= max]
    h <- hist(x, breaks=breaks, plot=FALSE)
    d <- data.frame(lower=h$breaks[-length(h$breaks)], upper=h$breaks[-1], mid=h$mids,
                    count=h$counts, p=h$counts/length(x), density=h$density)

    ## set all probabilities to 0 if x is empty
    if (length(x) == 0)
        d$p <- 0

    d
}

#' Obtain the approximate log likelihood of x given
#' the model that produced xhat, where x and xhat are vectors
#' of countinuous values. This function creates a histogram
#' of x and of xhat by binning the values with the discretize function.
#' Then, the likelihood of observing x within any histogram bin is
#' simply the proportion of xhat that lies inside that bin.
#'
#' @param x a vector of raw data
#' @param xhat a vector of simulated data from a model
#' @param delta a small fraction (0 < delta << 1) added to the histogram
#' of xhat to ensure that xhat places non-zero probability at each bin.
#' This avoids taking log(0), creating infinite log likelihoods.
#' @param ... optional arguments to discretize (min, max, binwidth)
#' @return the approximate log likelihood of x given the model that produced xhat
LL_discrete <- function(x, xhat, delta=1e-250, ...) {
    x.disc <- discretize(x, ...)
    xhat.disc <- discretize(xhat, ...)
    sum(x.disc$count * log(xhat.disc$p + delta))
}

#' Obtain the approximate negative log likelihood of x given
#' the model that produced xhat, where x and xhat are vectors
#' of countinuous values. This function creates a histogram
#' of x and of xhat by binning the values with the discretize function.
#' Then, the likelihood of observing x within any histogram bin is
#' simply the proportion of xhat that lies inside that bin.
#'
#' @param x a vector of raw data
#' @param xhat a vector of simulated data from a model
#' @param ... optional arguments to discretize (min, max, binwidth)
#' @return the approximate negative log likelihood of x
#' given the model that produced xhat
NLL_discrete <- function(x, xhat, ...) {
    -LL_discrete(x, xhat, ...)
}


between <- function(x, xmin, xmax) {
    if (is.na(x) || is.null(x))
        stop('Error: x is null')
    if (is.na(xmin) || is.null(xmin))
        stop('Error: xmin is null')
    if (is.na(xmax) || is.null(xmax))
        stop('Error: xmax is null')
    
    x >= xmin & x <= xmax
}

in_bounds <- function(parameters, bounds) {
    all(sapply(names(bounds),
               function (p) {
                   if (is.na(parameters[[p]]) || is.null(parameters[[p]]))
                       stop(paste0('Error: parameter ', p, ' is null.'))
                   
                   between(parameters[[p]],
                           bounds[[p]][1],
                           bounds[[p]][2])
               }))
}

#' Get a function to calculate the LL of x given some parameterization of the UCM
#'
#' @param x a vector of fixation durations to model using the UCM
#' @param bounds a named list of UCM parameters with upper and lower bounds.
#' If empty, parameter bounds will not be checked.
#' @param n_trials the number of UCM simulations to run in calculating the LL
#' @param trial_dur the duration of each simulated trial in seconds
#' @param min, max, binwidth parameters used to bin fixation durations into histograms
#' @param delta the minimum probability allowed in each histogram bin (used to avoid infinite values)
#' @param default_params a named list of default parameter values provided to the returned log likelihood function.
#' @return a function that calculates the LL of x given some parameterization of the UCM
LL <- function(x, bounds, n_trials, trial_dur, min, max, binwidth, delta=1/n_trials, default_params=list()) {
    function(N_states=default_params$N_states, t_timer=default_params$t_timer,
             t_labile=default_params$t_labile, t_nonlabile=default_params$t_nonlabile,
             t_motor=default_params$t_motor, t_execution=default_params$t_execution,
             modulation=1, N_mod=1) {
        ## calculate N, initialize LL variables
        N <- round(N_states * N_mod)
        ll <- 0
        
        ## calculate LL, assigning minimum value for out-of-bound parameter settings
        if (!in_bounds(list(N_states=N, t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile,
                            t_motor=t_motor, t_execution=t_execution, modulation=modulation, N_mod=N_mod),
                       bounds)) {
            ll <- LL_discrete(x, numeric(), min=min, max=max, binwidth=binwidth, delta=delta)        
        } else {
            ll <- 1:n_trials %>%
                mclapply(function (i) {
                    UCM(N_timer=N, N_labile=N, N_nonlabile=N, N_motor=N_states, N_execution=N_states,  ## only modulate first three stages
                        t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile,
                        t_motor=t_motor, t_execution=t_execution, modulation=modulation, mon=0) %>%
                        run(trial_dur) %>%
                        wrap()
                }) %>%
                get_fixations() %>%
                pull(duration) %>%
                LL_discrete(x, ., min=min, max=max, binwidth=binwidth, delta=delta)
        }
        
        list(Score=ll, n_trials=n_trials, min=min, max=max, binwidth=binwidth, delta=delta)
    }
}


#' Get a function returning the log likelihood (LL) of the MW data given some parameterization of UCM
#'
#' @param N_states the number of state transitions for *all* stages
#' @param t_timer, t_labile, t_nonlabile the mean duration of the timer, labile, and nonlabile stages
#' @return a function taking the argument MW_mod (the rate modulation factor), returning the same list as LL.on_task
LL.joint <- function(x.on_task, x.mw, bounds, n_trials, trial_dur, min, max, binwidth, delta=1/n_trials,
                     default_params=list()) {
    ## Create separate helper LL functions for on_task & MW trials
    LL.on_task <- LL(x.on_task, bounds, n_trials, trial_dur, min, max, binwidth, delta)
    LL.mw <- LL(x.mw, bounds, n_trials, trial_dur, min, max, binwidth, delta)
    
    function(N_states=default_params$N_states, t_timer=default_params$t_timer,
             t_labile=default_params$t_labile, t_nonlabile=default_params$t_nonlabile,
             t_motor=default_params$t_motor, t_execution=default_params$t_execution,
             modulation=1, N_mod=1) {
        ll.on_task <- LL.on_task(N_states, t_timer, t_labile, t_nonlabile, t_motor, t_execution)
        ll.mw <- LL.mw(N_states, t_timer, t_labile, t_nonlabile, t_motor, t_execution,
                       modulation=modulation, N_mod=N_mod)
        
        list(Score=ll.on_task$Score+ll.mw$Score, LL.on_task=ll.on_task$Score, LL.mw=ll.mw$Score,
             n_trials=ll.on_task$n_trials, min=ll.on_task$min,
             max=ll.on_task$max, binwidth=ll.on_task$binwidth, delta=ll.on_task$delta)
    }
}



#' Fit the UCM using bayesOpt
#'
#' @param saveFile the file to save the model to
#' @param LL.fun the likelihood function being optimized
#' @param bounds the parameter bounds to explore (see bayesOpt)
#' @param initPoints the number of initial parameter settings to explore
#' @param kappa the weighting of uncertainty information
#' @param utilityThresh the convergence threshold (between 0 and 1)
optimize <- function(saveFile, LL.fun, bounds, initPoints=16, kappa=5.0, utilityThresh=0.001) {
    if (file.exists(saveFile)) {
        bopt <- readRDS(saveFile)   ## load the cached optimizer
    } else {
        ## start the optimizer from scratch
        bopt <- bayesOpt(LL.fun, bounds=bounds, saveFile=saveFile, initPoints=initPoints,
                         iters.n=5, kappa=kappa, plotProgress=TRUE)
    }
    
    ## run the optimizer until convergence
    while (bopt$scoreSummary$gpUtility[nrow(bopt$scoreSummary)] >= utilityThresh) {
        bopt <- addIterations(bopt, iters.n=1, plotProgress=TRUE)
    }

    bopt
}


## Test for differences in distributions using KS test
KS <- function(fix) {
    options(warn=-1)
    D <- ks.test(fix %>% filter(mw==0) %>% pull(fixdur),
                 fix %>% filter(mw==1) %>% pull(fixdur),
                 exact=FALSE)$statistic
    options(warn=0)
    return(D)
}

## Perform a KS test using within-participant permutations
KS.permutation <- function(fix, n=1000) {
    fix %>%
        group_by(ID) %>%
        nest() %>%
        expand_grid(permutation=1:n) %>%
        mutate(data=map(data, ~ mutate(., mw=sample(mw)))) %>%
        unnest(data) %>%
        group_by(permutation) %>%
        nest() %>%
        mutate(D.null=map_dbl(data, KS),
               D=KS(fix))
}

mw_plot <- function(title, data, binwidth) {
    p.1 <- ggplot(data) +
        aes(x=fixdur, fill=factor(mw)) +
        geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=binwidth, show.legend=FALSE) +
        ylab('Proportion') +
        scale_fill_jco() +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
    
    p.2 <- ggplot(data) +
        aes(x=fixdur, y=factor(mw), fill=factor(mw)) +
        stat_summary(fun.data=mean_cl_normal, geom='col') +
        stat_summary(fun.data=mean_cl_normal, geom='errorbar', width=.5) +
        scale_fill_jco(name='Probe Response', labels=c('Attentive Viewing', 'Mind Wandering')) +
        scale_x_continuous(name='Fixation Duration (ms)', labels=ms_format()) +
        theme_bw() +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major.y=element_blank())
    
    (p.1 / p.2 & coord_cartesian(xlim=c(0, 1))) +
        plot_layout(heights=c(1, .25), guides='collect') +
        plot_annotation(title=title)
}


mw_hist <- function(fix, legend=FALSE, x.axis=FALSE) {
    p <- fix %>% unnest(hist) %>%
        ggplot(aes(x=mid, y=p, color=factor(mw), linetype=type)) +
        geom_line(size=.75, show.legend=c(color=legend)) +
        scale_linetype_manual(name='', values=c('dotted', 'solid')) +
        scale_color_jco(name='Probe Response', labels=c('Attentive Viewing', 'Mind Wandering')) +
        scale_x_continuous(name='Fixation Duration (ms)', labels=ms_format()) +
        scale_y_continuous('Proportion', labels=no_leading_zeros) +
        theme_bw()

    if (!x.axis)
        p <- p +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())
    p
}

mw_means_plot <- function(fix) {
    fix %>% unnest(data) %>%
        ggplot() +
        aes(x=fixdur, y=type, fill=factor(mw)) +
        stat_summary(fun.data=mean_cl_normal, geom='col', position=position_dodge(.95)) +
        stat_summary(fun.data=mean_cl_normal, geom='errorbar', width=.5, position=position_dodge(.95)) +
        scale_fill_jco(name='Probe Response', labels=c('Attentive Viewing', 'Mind Wandering')) +
        scale_x_continuous(name='Fixation Duration (ms)', labels=ms_format()) +
        theme_bw() +
        theme(axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major.y=element_blank(),
              strip.background=element_blank(),
              strip.text=element_blank())
}

reciprobit_plot <- function(title, fix, min) {
    inverse_probit_trans <- trans_new('inverse_probit', qnorm, pnorm)
    
    fix %>% unnest(data) %>%
        mutate(fixrate=1/fixdur) %>%
        ggplot(aes(x=fixrate, color=type)) +
        stat_ecdf() +
        scale_color_brewer(name='', palette='Set1') +
        facet_grid( ~ mw, labeller=labeller(mw=c('0'='Attentive Viewing', '1'='Mind Wandering'))) +
        scale_x_continuous(name='Fixation Duration (s)',
                           trans='reverse',
                           breaks=seq(0, 1/min),
                           labels=round(1/seq(0, 1/min), 2)) +
        scale_y_continuous(name='Cumulative Probability',
                           trans=inverse_probit_trans,
                           limits=pnorm(c(-4,4)),
                           breaks=pnorm(-4:4),
                           labels=round(pnorm(-4:4), 3)) +
        theme_bw() +
        coord_cartesian(xlim=c(1/min, 0)) +
        ggtitle(title)    
}

mw_ecdf_plot <- function(fix, n=5000) {
    fix %>%
        unnest(data) %>%
        ggplot(aes(x=fixdur, group=interaction(mw, type), color=factor(mw), linetype=type)) +
        stat_ecdf(n=n, geom='line') +
        scale_color_jco(name='Probe Response', labels=c('Attentive Viewing', 'Mind Wandering')) +
        scale_linetype_manual(name='', values=c('dotted', 'solid')) +
        scale_x_continuous(name='Fixation Duration (ms)', labels=ms_format()) +
        scale_y_continuous(name='Cumulative Probability', expand=c(0, 0)) +
        coord_cartesian(xlim=c(0, 1)) +
        theme_bw()
}

mw_ecdf_difference_plot <- function(fix, n=5000) {
    fix %>%
        mutate(ecdf_fun=map(data, ~ ecdf(.$fixdur))) %>%
        expand_grid(fixdur=seq(0, 2, length.out=n)) %>%
        mutate(ecdf=map2_dbl(ecdf_fun, fixdur, ~ .x(.y))) %>%
        select(type, mw, fixdur, ecdf) %>%
        pivot_wider(names_from=mw, values_from=ecdf) %>%
        mutate(ecdf_diff=`1` - `0`) %>%
        ggplot(aes(x=fixdur, y=ecdf_diff, color=type, fill=type)) +
        geom_line() +
        geom_area(alpha=.33, position='identity') +
        scale_x_continuous(name='Fixation Duration (ms)', labels=ms_format()) +
        scale_y_continuous(name='Cumulative Probability Difference\n(Mind Wandering - Attentive Viewing)') +
        scale_color_brewer(name='', palette='Set1') +
        scale_fill_brewer(name='', palette='Set1') +
        coord_cartesian(xlim=c(0, 1)) +
        theme_bw()
}

## Fit LATER on the dataframe df
LATER.fit <- function(df, n_samples) {
    df %>%
        mutate(fixrate=1/fixdur) %>%
        summarize(m=mean(fixrate, na.rm=TRUE),
                  s=sd(fixrate, na.rm=TRUE)) %>%
        expand_grid(sample=1:n_samples) %>%
        mutate(type='LATER',
               fixrate=rnorm(n(), m, s),
               fixdur=1/fixrate)
}
