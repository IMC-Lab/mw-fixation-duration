library(simmer)
library(dplyr)
library(tidyr)

#' Create a trajectory with a discrete state random walk from 0 to N.
#'
#' @param .trj The trajectory to add a random walk to. By default, creates a new blank trajectory.
#' @param env The simmer environment that will simulate this trajectory
#' @param N The number of states of the random walk (> 0)
#' @param t The mean duration of the random walk (> 0)
#' @param modulation Logical indicating whether the rate of the random walk should be multiplied by the global attribute 'modulation'
#' @return .trj with a random walk added to its trajectory
#' @examples
#' ## Simulate a random walk with 10 states for a mean duration of 1s
#' env <- simmer()
#' walk <- trajectory() %>%
#'     log_('starting random walk') %>%
#'     randomWalk(env, N=10, t=1.0) %>%
#'     log_('finished random walk')
#'
#' env %>%
#'     add_generator('walk', walk, at(0)) %>%
#'     run()
#' @export
randomWalk <- function(.trj=trajectory(), env, N, t, modulation=TRUE) {
    if (N <= 0)
        stop('Cannot create a random walk with N<0 states.')
    if (t <= 0)
        stop('Cannot create a random walk with mean duration t<0.')

    ## Define a function to randomly draw wait times.
    ## If modulation is true, this time is multiplied by the global
    ## attribute 'modulation'.
    wait_time <- ifelse(modulation,
                        function() rexp(1, get_global(env, 'modulation') * N/t),
                        function() rexp(1, N/t))
    
    .trj %>%
        set_attribute('state', 0) %>%
        timeout(wait_time) %>%
        set_attribute('state', 1, mod='+') %>%
        rollback(2, check=function() get_attribute(env, 'state') < N-1) %>%
        timeout(wait_time)
}

#' Create a simmer environment to simulate the UCM model of fixation durations.
#' This model is comprised of random walks over 5 stages (timer, labile, nonlabile,
#' motor, and saccade execution). Each stage has two parameters: N is the number
#' of states in the random walk for that stage, and t is the average duration of
#' that stage. The timer cycles continuously from 0 to N. Every time it reaches N,
#' it begins a new timer and enters the labile stage. Any previous walks in the labile
#' stage are cancelled at this point. From there, the walk enters the nonlabile stage,
#' the motor planning stage, and finally saccade execution. A fixation begins at
#' the end of a saccade execution and ends at the start of the next saccade.
#'
#' Random walks can be locally modified in two ways: rate modulation and cancellation.
#' Walks in the timer, labile, and nonlabile stages can be sped up or slowed down using
#' the global attribute "modulation", which is a multiplicative factor (e.g.
#' set_global("modulation", 2.0) doubles the speed of the first three stages).
#' Walks in the labile stage can be cancelled by sending the signal "cancel-labile"
#' (e.g. send('cancel-labile')).
#'
#' @param N_timer The number of states in the random timer (> 0).
#' @param N_labile The number of states in the labile stage (> 0).
#' @param N_nonlabile The number of states in the nonlabile stage (> 0).
#' @param N_motor The number of states in the motor programming stage (> 0).
#' @param N_execution The number of states in the saccade execution stage (> 0).
#' @param t_timer The mean duration of the random timer (> 0).
#' @param t_labile The mean duration of the labile stage (> 0).
#' @param t_nonlabile The mean duration of the nonlabile stage (> 0).
#' @param t_motor The mean duration of the motor programming stage (> 0).
#' @param t_execution The mean duration of the saccade execution stage (> 0).
#' @return A simmer environment simulating the UCM model.
#' @examples
#' ## Run the baseline model for 5 seconds
#' UCM() %>% run(until=5)
#'
#' ## Run UCM with custom parameters for 10 seconds
#' UCM(N_motor=5, N_execution=5, t_timer=.5) %>% run(until=10)
#' @export
UCM <- function(N_timer=14, N_labile=14, N_nonlabile=14, N_motor=14, N_execution=14,
                t_timer=.25, t_labile=.175, t_nonlabile=.07, t_motor=.03, t_execution=.02) {
    if (!all(c(N_timer, N_labile, N_nonlabile, N_motor, N_execution) > 0))
        stop('Cannot create UCM with N < 0')
    if (!all(c(t_timer, t_labile, t_nonlabile, t_motor, t_execution) > 0))
        stop('Cannot create UCM with t < 0')
    
    s <- simmer()
    
    ## create a trajectory for the timer
    timer <- trajectory() %>%
        set_attribute('id', function() get_n_generated(s, 'timer.')) %>%
        ## cancel any labile programs at the start of a new timer
        send('cancel-labile') %>%
        randomWalk(s, N_timer, t_timer, modulation=TRUE) %>%
        ## start a new timer and labile program after completion
        activate('timer.') %>%
        activate('labile.')
    
    ## create a trajectory for the labile stage
    labile <- trajectory() %>%
        set_attribute('id', function() get_n_generated(s, 'labile.')) %>%
        set_attribute('cancelled', 0) %>%
        ## cancel if a new timer has started
        renege_if('cancel-labile',
                  out=trajectory() %>%
                      set_attribute('cancelled', 1)) %>%
        randomWalk(s, N_labile, t_labile, modulation=TRUE) %>%
        set_global('labile_id', function() get_attribute(s, 'id')) %>%
        ## start a new nonlabile program after completion
        activate('nonlabile.')
    
    ## create a trajectory for the nonlabile stage
    nonlabile <- trajectory() %>%
        set_attribute('id', function() get_global(s, 'labile_id')) %>%
        randomWalk(s, N_nonlabile, t_nonlabile, modulation=TRUE) %>%
        set_global('nonlabile_id', function() get_attribute(s, 'id')) %>%
        ## start a motor program after completion
        activate('motor.')
    
    ## create a trajectory for the motor planning stage
    motor <- trajectory() %>%
        set_attribute('id', function() get_global(s, 'nonlabile_id')) %>%
        randomWalk(s, N_motor, t_motor, modulation=FALSE) %>%
        set_global('motor_id', function() get_attribute(s, 'id')) %>%
        ## start a new saccade execution after completion
        activate('execution.')
    
    ## create a trajectory for the saccade execution stage
    execution <- trajectory() %>%
        send('end-fixation') %>%
        set_attribute('id', function() get_global(s, 'motor_id')) %>%
        randomWalk(s, N_execution, t_execution, modulation=FALSE) %>%
        set_global('execution_id', function() get_attribute(s, 'id')) %>%
        activate('fixation.')

    ## create a trajectory for fixations
    fixation <- trajectory() %>%
        set_attribute('id', function() get_global(s, 'execution_id')) %>%
        set_attribute('state', 0) %>%
        trap('end-fixation') %>%
        wait()
    
    s %>%
        add_generator('init.', trajectory() %>%
                               set_global('modulation', 1) %>%
                               activate('timer.'), at(0)) %>%
        add_generator('timer.', timer, when_activated(), mon=2) %>%
        add_generator('labile.', labile, when_activated(), mon=2) %>%
        add_generator('nonlabile.', nonlabile, when_activated(), mon=2) %>%
        add_generator('motor.', motor, when_activated(), mon=2) %>%
        add_generator('execution.', execution, when_activated(), mon=2) %>%
        add_generator('fixation.', fixation, when_activated(), mon=2)
}

#' Get a dataframe containing all of the IDs of each trajectory.
#' IDs can be used to tie fixations to particular cycles of each stage
#'
#' @param ucm The UCM simmer environment to get trajectory IDs for.
#' @return a dataframe containing the ID of every trajectory so far.
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     get_ids()
#' @export
get_ids <- function(ucm) {
    ucm %>%
        get_mon_attributes() %>%
        tibble() %>%
        filter(key=='id') %>%
        select(-time) %>%
        separate(name, c('stage', 'n'), sep='\\.') %>%
        mutate(stage=factor(stage, levels=c('timer', 'labile', 'nonlabile', 'motor', 'execution', 'fixation'))) %>%
        pivot_wider(names_from=key)
}

#' Get a dataframe specifying which IDs were cancelled during the labile stage.
#' cancelled states whether the ith labile program was cancelled or not.
#' n_cancellations states how many cancellations are within the saccade generated
#' by the ith timer.
#' 
#' @param ucm The UCM simmer environment to get cancellation information for.
#' @return a dataframe detailing which IDs were cancelled, and how many cancellations were present before each uncancelled ID
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     get_cancellations()
#' @export
get_cancellations <- function(ucm) {
    ucm %>%
        get_mon_attributes() %>%
        tibble() %>%
        filter(key == 'cancelled') %>%
        separate(name, c('stage', 'n'), sep='\\.') %>%
        mutate(n=as.numeric(n)) %>%
        group_by(stage, n, replication, key) %>%
        summarize(value=max(value), .groups='drop') %>%
        pivot_wider(names_from=key) %>%
        arrange(replication, n) %>%
        group_by(replication) %>%
        ## count the number of cancellations
        mutate(s=sequence(rle(cancelled)$lengths),
               g=cumsum(s==1)) %>%
        group_by(replication, g) %>%
        mutate(count=n()) %>%
        group_by(replication) %>%
        mutate(n_cancellations=ifelse(cancelled==0 & s==count, lead(count, default=0), 0),
               n=as.character(n)) %>%
        select(-s, -g, -count) %>%
        ungroup() %>%
        left_join(filter(get_ids(ucm), stage=='labile')) %>%
        select(-stage, -n) %>%
        relocate(id)
}

#' Get a dataframe specifying every state transition of every random walk
#'
#' @param ucm The UCM simmer environment to get state information for.
#' @param add A logical indicating whether the state numbers of each stage should begin at the end of the previous stage (TRUE; default) or at 0 (FALSE).
#' @param N_timer The number of states in the random timer (> 0).
#' @param N_labile The number of states in the labile stage (> 0).
#' @param N_nonlabile The number of states in the nonlabile stage (> 0).
#' @param N_motor The number of states in the motor programming stage (> 0).
#' @param N_execution The number of states in the saccade execution stage (> 0).
#' @return a dataframe specifying the time of every stage transiton in UCM
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     get_states()
#' @export
get_states <- function(ucm) {
    s <- ucm %>%
        get_mon_attributes() %>%
        tibble() %>%
        filter(key == 'state') %>%
        separate(name, c('stage', 'n'), sep='\\.') %>%
        mutate(stage=factor(stage, levels=c('timer', 'labile', 'nonlabile', 'motor', 'execution', 'fixation'))) %>%
        pivot_wider(names_from=key)
    
    N <- s %>% group_by(stage) %>%
        summarize(N=length(unique(state))) %>%
        pivot_wider(names_from=stage, values_from=N)
    
    s %>%
        mutate(cum_state=ifelse(stage=='labile', N$timer+state,
                         ifelse(stage=='nonlabile', N$timer+N$labile+state,
                         ifelse(stage=='motor', N$timer+N$labile+N$nonlabile+state,
                         ifelse(stage=='execution', N$timer+N$labile+N$nonlabile+N$motor+state,
                         ifelse(stage=='fixation', N$timer+N$labile+N$nonlabile+N$motor+N$execution+state,
                                state)))))) %>%
        left_join(get_ids(ucm))
}

#' Get a dataframe specifying the start, end, and duration of every simulated fixation.
#'
#' @param ucm The UCM simmer environment to get fixations from
#' @return A dataframe specifying information about each fixation
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     get_fixations()
#' @export
get_fixations <- function(ucm) {
    ucm %>%
        get_mon_arrivals() %>%
        tibble() %>%
        separate(name, c('stage', 'n'), sep='\\.') %>%
        filter(stage == 'fixation') %>%
        mutate(duration=end_time-start_time,
               stage=factor(stage, levels=c('timer', 'labile', 'nonlabile',
                                            'motor', 'execution', 'fixation'))) %>%
        left_join(get_ids(ucm), by=c('stage', 'n', 'replication')) %>%
        select(-stage, -activity_time, -finished) %>%
        relocate(id, n, replication)
}

#' Get a traceplot of UCM's simulation from start to end.
#'
#' @param ucm The UCM simmer environment to plot traces for
#' @param start The start time of the plot
#' @param end The end time of the plot
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     trace_plot()
#' @export
trace_plot <- function(ucm, start=0, end=now(ucm)) {
    f <- ucm %>%
        get_fixations() %>%
        filter(end_time >= start & start_time <= end)
    s <- ucm %>%
        get_states() %>%
        filter(time >= start & time <= end)

    ## Get the starting state of each stage
    N <- s %>% group_by(stage) %>%
        summarize(N=min(cum_state)) %>%
        pivot_wider(names_from=stage, values_from=N) %>%
        as.numeric
    
    s %>%
        ggplot(aes(x=time, y=cum_state, group=id)) +
        geom_rect(aes(xmin=start_time, xmax=end_time, ymin=-Inf, ymax=Inf),
                  inherit.aes=FALSE, alpha=0.25, data=f) +
        geom_hline(yintercept=N) +
        geom_step() +
        scale_x_continuous(name='Time (s)', limits=c(start, end), expand=c(0, 0)) +
        scale_y_continuous(breaks=N[-length(N)] + diff(N)/2,
                           labels=c('timer', 'labile', 'non-labile', 'motor', 'execution'),
                           limits=c(0, max(s$cum_state)),
                           expand=c(0, 0)) +
        theme_bw() +
        theme(axis.text.y=element_text(angle=90, hjust=0.5, vjust=0.5),
              axis.title.y=element_blank())
}


cancellations_hist <- function(ucm, max_duration=1.2, binwidth=0.06) {
    c <- get_cancellations(envs)
    f <- left_join(get_fixations(envs), c)
    
    f %>%
        filter(activity_time <= max_diration) %>%
        mutate(n_cancellations=pmin(n_cancellations, 3)) %>%
        ggplot(aes(x=activity_time, group=n_cancellations, fill=factor(n_cancellations))) +
        scale_fill_viridis(name='Cancellations', discrete=TRUE) +
        geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=binwidth) +
        scale_x_continuous(limits=c(0, 1.2)) +
        xlab('Fixation Duration (s)') + ylab('Proportion') +
        theme_bw()
}


#' Turn a vector into a discrete distribution of binned values
#'
#' @param x the vector of values to discretize
#' @param min the minimum value of allowed values in x
#' @param max the maximum value of allowed values in x
#' @param binsize the width of the distribution's bins
#' @return a dataframe with the following columns, where each row corresponds to a bin:
#'  - lower: the lower bound of the bin
#'  - upper: the upper bound of the bin
#'  - mid: the midpoint of the bin
#'  - count: the number of data points in this bin
#'  - p: the probability of finding a data point in this bin
discretize <- function(x, min=0, max=1.2, binsize=.06) {
    breaks <- seq(min, max, by=binsize)
    x <- x[x >= min & x <= max]
    h <- hist(x, breaks=breaks, plot=FALSE)
    data.frame(lower=h$breaks[-length(h$breaks)], upper=h$breaks[-1], mid=h$mids,
               count=h$counts, p=h$counts/length(x))
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
#' @param ... optional arguments to discretize (min, max, binsize)
#' @return the approximate log likelihood of x given the model that produced xhat
LL_discrete <- function(x, xhat, delta=1e-10, ...) {
    x.disc <- do.call(discretize, list(x, ...))
    xhat.disc <- do.call(discretize, list(xhat, ...))
    sum(x.disc$count * log(xhat.disc$p+delta))
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
#' @param ... optional arguments to discretize (min, max, binsize)
#' @return the approximate negative log likelihood of x
#' given the model that produced xhat
NL_discrete <- function(x, xhat, ...) {
    -do.call(LL_discrete, x, xhat, ...)
}


#' Return a function to run the UCM and calculate its approximate
#' Log Likelihood (LL) with respect to the data (fixation durations) x.
#'
#' @param x the data for which to calculate the NLL
#' @param N_trials the number of UCM instances to run
#' @param trial_duration the length of time to run each instance of the UCM
#' @param parallel run each instance of the UCM as a separate process?
#' @param save if TRUE, the returned function writes the simulated fixations to a file
#' which can be used to amortize computational load
#' @param dir if save==TRUE, the directory in which to save files
#' @param ... optional arguments passed to LL_discrete and to discretize.
#' @return a function taking N_states, t_timer, t_labile, t_nonlabile,
#' t_motor, and t_execution as arguments (see UCM). Calling this function
#' will run N_trials instances of UCM for trial_duration seconds, and return
#' the LL of the data given the simulated data as a list under the name Score
#' (for use with ParBayesOptimization).
LL.UCM <- function (x, N_trials, trial_duration, parallel=TRUE, save=TRUE, dir='.', ...) {
    if (parallel) {
        function (N_states, t_timer, t_labile, t_nonlabile) {
            fname <- paste0(dir, '/ucm_sim_', N_trials, '_', trial_duration, '_',
                    N_states, '_', t_timer, '_', t_labile, '_',
                    t_nonlabile, '.csv')
            
            if (file.exists(fname)) {
                fix <- read.csv(fname)
            } else {
                fix <- 1:N_trials %>%
                    mclapply(function (i) {
                        UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                            N_motor=N_states, N_execution=N_states,
                            t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile) %>%
                            run(trial_duration) %>%
                            wrap()
                    }) %>%
                    get_fixations()
            }
            
            if (save)
                write.csv(fix, fname)

            list(Score=do.call(LL_discrete, list(x, fix$activity_time, ...)))
        }            
    } else {
        function (N_states, t_timer, t_labile, t_nonlabile) {
            fname <- paste0(dir, '/ucm_sim_', N_trials, '_', trial_duration, '_',
                    N_states, '_', t_timer, '_', t_labile, '_',
                    t_nonlabile, '.csv')
            
            if (file.exists(fname)) {
                fix <- read.csv(fname)
            } else {
                fix <- 1:N_trials %>%
                    lapply(function (i) {
                        UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                            N_motor=N_states, N_execution=N_states,
                            t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile) %>%
                            run(trial_duration) %>%
                            wrap()
                    }) %>%
                    get_fixations()
            }
            
            if (save)
                write.csv(fix, fname)
            
            list(Score=do.call(LL_discrete, list(x, fix$activity_time, ...)))
        }
    }
}
