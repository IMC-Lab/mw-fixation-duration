library(httr)
library(simmer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggsci)
library(parallel)
library(ParBayesianOptimization)

## Load UCM source code
source('UCM.R')
options(mc.cores=parallel::detectCores())

## Set optimization parameters
T_MOTOR <- .03           ## duration of motor planning stage (fixed parameter)
T_EXECUTION <- .02       ## duration of saccade execution (fixed parameter)
N_TRIALS <- 2500         ## number of UCM trials to simulate per epoch to estimate LL
MIN_FIXDUR <- 0.08       ## only keep fixations over this limit
MAX_FIXDUR <- 2          ## only keep fixations under this limit
BINWIDTH <- 0.04         ## the width of histogram bins to calculate LL
INIT_POINTS <- 16        ## number of points to initialize optimization
KAPPA <- 3.0             ## importance of uncertainty in optimization
UTILITY_THRESH <- 0.05   ## stop optimization when utility is smaller than this value
on_task_bounds <- list(N_states=c(10L, 30L),
                       t_timer=c(.1, .3),
                       t_labile=c(.1, .25),
                       t_nonlabile=c(.05, .1))
mw_bounds <- list(MW_mod=c(0.25, 1))
bounds <- c(on_task_bounds, mw_bounds)

## Load data
GET('https://osf.io/fpqsj/download', write_disk(tf <- tempfile(fileext = ".csv")))
fix_data <- read.csv(tf, header=TRUE) %>%
    tibble %>%
    mutate(fixdur=duration/1000,
           proberesp=factor(Attention, levels=c('On-task', 'Mind-wandering')),
           mw=ifelse(proberesp == 'Mind-wandering', 1, 0)) %>%
    ## filter as per the manuscript
    drop_na() %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR) %>%
    filter(x_pos > 0 & x_pos < 1024) %>%
    filter(y_pos > 0 & y_pos < 768)

fix_data.on_task <- fix_data[fix_data$mw == 0,]
fix_data.mw <- fix_data[fix_data$mw == 1,]

## Plot data
p.1 <- ggplot(fix_data) +
    aes(x=fixdur, fill=factor(mw)) +
    geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=BINWIDTH, show.legend=FALSE) +
    ylab('Proportion') +
    scale_fill_jco() +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

p.2 <- ggplot(fix_data) +
    aes(x=fixdur, y=Attention, fill=factor(mw)) +
    stat_summary(fun.data=mean_cl_normal, geom='col') +
    stat_summary(fun.data=mean_cl_normal, geom='errorbar') +
    scale_fill_jco(name='Probe Response', labels=c('On-task', 'Mind-wandering')) +
    xlab('Fixation Duration (ms)') +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y=element_blank())

(p.1 / p.2 & xlim(0, 1)) +
    plot_layout(heights=c(1, .25), guides='collect') +
    plot_annotation(title='Zhang et al., (2021) Raw Data')
ggsave('plots/hist_zhang_2021.png', width=8, height=6)


##################################################################################################
##                 Log-Likelihood Functions
##################################################################################################

#' Get the log likelihood (LL) of the on-task data given some parameterization of UCM
#'
#' @param N_states the number of state transitions for *all* stages
#' @param t_timer, t_labile, t_nonlabile the mean duration of the timer, labile, and nonlabile stages
#' @return a list containing:
#'  Score- the LL of on-task trials
#'  t_motor- the mean duration of the motor planning stage (fixed parameter)
#'  t_execution- the mean duration of saccade execution (fixed parameter)
#'  N_TRIALS- the number of simulated UCM trials per condition
#'  min- the minimum fixation duration
#'  max- the maximum fixation duration
#'  binwidth- the width of the histogram bins used to approximate the LL
LL.on_task <- function(N_states, t_timer, t_labile, t_nonlabile) {
    LL <- 1:N_TRIALS %>%
        mclapply(function (i) {
            UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                N_motor=N_states, N_execution=N_states,
                t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile,
                t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
                run(10) %>%
                wrap()
        }) %>%
        get_fixations() %>%
        pull(duration) %>%
        LL_discrete(fix_data.on_task$fixdur, ., min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
    
    list(Score=LL, t_motor=T_MOTOR, t_execution=T_EXECUTION,
         N_TRIALS=N_TRIALS, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
}

#' Get a function returning the log likelihood (LL) of the MW data given some parameterization of UCM
#'
#' @param N_states the number of state transitions for *all* stages
#' @param t_timer, t_labile, t_nonlabile the mean duration of the timer, labile, and nonlabile stages
#' @return a function taking the argument MW_mod (the rate modulation factor), returning the same list as LL.on_task
LL.mw <- function(N_states, t_timer, t_labile, t_nonlabile) {
    function (MW_mod) {
        LL <- 1:N_TRIALS %>%
            mclapply(function (i) {
                UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                    N_motor=N_states, N_execution=N_states,
                    t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile,
                    t_motor=T_MOTOR, t_execution=T_EXECUTION, modulation=MW_mod) %>%
                    run(10) %>%
                    wrap()
            }) %>%
            get_fixations() %>%
            pull(duration) %>%
            LL_discrete(fix_data.mw$fixdur, ., min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
        
        list(Score=LL, t_motor=T_MOTOR, t_execution=T_EXECUTION,
             N_TRIALS=N_TRIALS, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
    }
}

#' Get the log likelihood (LL) of all of the data given some parameterization of UCM
#'
#' @param N_states the number of state transitions for *all* stages
#' @param t_timer, t_labile, t_nonlabile the mean duration of the timer, labile, and nonlabile stages
#' @param MW_mod the factor by which rates are modulated during mind-wandering
#' @return a list containing:
#'  Score- the LL of the data
#'  LL.on_task- the LL of on-task trials
#'  LL.mw- the LL of mind-wandering trials
#'  N_TRIALS- the number of simulated UCM trials per condition
#'  min- the minimum fixation duration
#'  max- the maximum fixation duration
#'  binwidth- the width of the histogram bins used to approximate the LL
LL <- function(N_states, t_timer, t_labile, t_nonlabile, MW_mod) {
    LL.on_task <- 1:N_TRIALS %>%
        mclapply(function (i) {
            UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                N_motor=N_states, N_execution=N_states,
                t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile,
                t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
                run(10) %>%
                wrap()
        }) %>%
        get_fixations() %>%
        pull(duration) %>%
        LL_discrete(fix_data.on_task$fixdur, ., min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)

    LL.mw <- 1:N_TRIALS %>%
        mclapply(function (i) {
            UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                N_motor=N_states, N_execution=N_states,
                t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile,
                t_motor=T_MOTOR, t_execution=T_EXECUTION,
                modulation=MW_mod) %>%
                run(10) %>%
                wrap()
        }) %>%
        get_fixations() %>%
        pull(duration) %>%
        LL_discrete(fix_data.mw$fixdur, ., min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)

    list(Score=LL.on_task + LL.mw, LL.on_task=LL.on_task, LL.mw=LL.mw,
         N_TRIALS=N_TRIALS, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
}




##################################################################################################
##                 Use this code to optimze on-task & mw trials simultaneously
##################################################################################################

## Optimize the UCM, loading a cached file if it exists
if (file.exists('zhang_2021.rds')) {
    bopt <- readRDS('zhang_2021.rds')   ## load the cached optimizer
} else {
    ## start the optimizer from scratch
    bopt <- bayesOpt(LL, bounds=bounds, saveFile='zhang_2021.rds', initPoints=INIT_POINTS,
                     iters.n=5, kappa=KAPPA, plotProgress=TRUE)
}

## run the optimizer until convergence
while (bopt$scoreSummary$gpUtility[nrow(bopt$scoreSummary)] >= UTILITY_THRESH) {
    bopt <- addIterations(bopt, iters.n=1, plotProgress=TRUE)
}


## Run the UCM with the best-fitting parameters
ucm.MLE <- getBestPars(bopt)
ucm.MLE



##################################################################################################
##                 Use this code to optimze on-task & mw trials separately
##################################################################################################

## Optimize the UCM on on-task trials, loading a cached file if it exists
if (file.exists('zhang_2021_on_task.rds')) {
    bopt.on_task <- readRDS('zhang_2021_on_task.rds')   ## load the cached optimizer
} else {
    ## start the optimizer from scratch
    bopt.on_task <- bayesOpt(LL.on_task, bounds=bounds.on_task, saveFile='zhang_2021_on_task.rds',
                             initPoints=INIT_POINTS, iters.n=5, kappa=KAPPA, plotProgress=TRUE)
}

## run the optimizer until convergence
while (bopt.on_task$scoreSummary$gpUtility[nrow(bopt.on_task$scoreSummary)] >= UTILITY_THRESH) {
    bopt.on_task <- addIterations(bopt.on_task, iters.n=1, plotProgress=TRUE)
}

ucm.MLE <- getBestPars(bopt.on_task)
ucm.MLE

## Optimize the UCM on MW trials, loading a cached file if it exists
if (file.exists('zhang_2021_mw.rds')) {
    bopt.mw <- readRDS('zhang_2021_mw.rds')   ## load the cached optimizer
} else {
    ## start the optimizer from scratch
    bopt.mw <- bayesOpt(do.call(LL.mw, ucm.MLE), bounds=bounds.mw, saveFile='zhang_2021_mw.rds',
                        initPoints=INIT_POINTS, iters.n=5, kappa=KAPPA, plotProgress=TRUE)
}

## run the optimizer until convergence
while (bopt.mw$scoreSummary$gpUtility[nrow(bopt.mw$scoreSummary)] >= UTILITY_THRESH) {
    bopt.mw <- addIterations(bopt.mw, iters.n=1, plotProgress=TRUE)
}

ucm.MLE <- c(getBestPars(bopt.on_task), getBestPars(bopt.mw))
ucm.MLE




##################################################################################################
##                 Print the best-fitting parameters & plot results
##################################################################################################
ucm.MLE


ucm.on_task <- mclapply(1:N_TRIALS, function(i) {
    UCM(N_timer=ucm.MLE$N_states, N_labile=ucm.MLE$N_states,
        N_nonlabile=ucm.MLE$N_states, N_motor=ucm.MLE$N_states,
        N_execution=ucm.MLE$N_states, t_timer=ucm.MLE$t_timer,
        t_labile=ucm.MLE$t_labile, t_nonlabile=ucm.MLE$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
        run(10) %>%
        wrap()
})

ucm.mw <- mclapply(1:N_TRIALS, function(i) {
    UCM(N_timer=ucm.MLE$N_states, N_labile=ucm.MLE$N_states,
        N_nonlabile=ucm.MLE$N_states, N_motor=ucm.MLE$N_states,
        N_execution=ucm.MLE$N_states, t_timer=ucm.MLE$t_timer,
        t_labile=ucm.MLE$t_labile, t_nonlabile=ucm.MLE$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION, modulation=ucm.MLE$MW_mod) %>%
        run(10) %>%
        wrap()
})

fix_ucm.on_task <- get_fixations(ucm.on_task)
fix_ucm.mw <- get_fixations(ucm.mw)


## Plot the simulated fixations against the raw data
rbind(mutate(discretize(fix_data.on_task$fix_dur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH),
             type='Data', mw=0),
      mutate(discretize(fix_data$fixdur[fix_data$mw==1], min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH),
             type='Data', mw=1),
      mutate(discretize(fix_ucm.on_task$duration, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH),
             type='UCM', mw=0),
      mutate(discretize(fix_ucm.mw$duration, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH),
             type='UCM', mw=1)) %>%
    ggplot(aes(x=mid, y=p, color=factor(mw), linetype=type)) +
    geom_line(size=1, show.legend=c(color=FALSE)) +
    scale_linetype_discrete(name='') +
    scale_color_jco(name='Probe Response', labels=c('On-task', 'Mind-wandering')) +
    facet_grid( ~ mw, labeller=labeller(mw=c('0'='On-task', '1'='Mind-wandering'))) +
    scale_x_continuous(limits=c(0, 2)) +
    xlab('Fixation Duration (s)') + ylab('Proportion') +
    theme_bw() + ggtitle('Zhang et al., (2021) UCM Fit')

ggsave('plots/fit_zhang_2021.png', width=8, height=4)
