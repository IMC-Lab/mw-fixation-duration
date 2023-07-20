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
randomWalk <- function(.trj=simmer::trajectory(), env, N, t, modulation=TRUE) {
    if (N <= 0)
        stop('Cannot create a random walk with N<0 states.')
    if (t <= 0)
        stop('Cannot create a random walk with mean duration t<0.')
    
    ## Define a function to randomly draw wait times.
    ## If modulation is true, this time is multiplied by the global
    ## attribute 'modulation'.
    wait_time <- ifelse(modulation,
                        function() rexp(1, simmer::get_global(env, 'modulation') * N/t),
                        function() rexp(1, N/t))
    
    .trj %>%
        simmer::set_attribute('state', 0) %>%
        simmer::timeout(wait_time) %>%
        simmer::set_attribute('state', 1, mod='+') %>%
        simmer::rollback(2, check=function() simmer::get_attribute(env, 'state') < N-1) %>%
        simmer::timeout(wait_time)
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
#' @param modulation The starting modulation of the timer, labile stage, and nonlabile stage.
#' @param mon An integer signal for whether to monitor individual state transitions in the UCM (0 = none; 1 = arrivals; 2 = attributes, default)
#' @return A simmer environment simulating the UCM model.
#' @examples
#' ## Run the baseline model for 5 seconds
#' UCM() %>% run(until=5)
#'
#' ## Run UCM with custom parameters for 10 seconds
#' UCM(N_motor=5, N_execution=5, t_timer=.5) %>% run(until=10)
#' @export
UCM <- function(N_timer=14, N_labile=14, N_nonlabile=14, N_motor=14, N_execution=14,
                t_timer=.25, t_labile=.175, t_nonlabile=.07, t_motor=.03, t_execution=.02,
                modulation=1.0, mon=2) {
    if (!all(c(N_timer, N_labile, N_nonlabile, N_motor, N_execution) > 0))
        stop('Cannot create UCM with N < 0')
    if (!all(c(t_timer, t_labile, t_nonlabile, t_motor, t_execution) > 0))
        stop('Cannot create UCM with t < 0')
    if (modulation <= 0 | modulation > 1)
        stop('Cannot create UCM with modulation <= 0 or > 1')
    
    s <- simmer::simmer()
    
    ## create a trajectory for the timer
    timer <- simmer::trajectory() %>%
        simmer::set_attribute('id', function() simmer::get_n_generated(s, 'timer.')) %>%
        ## cancel any labile programs at the start of a new timer
        simmer::send('cancel-labile') %>%
        randomWalk(s, N_timer, t_timer, modulation=TRUE) %>%
        ## start a new timer and labile program after completion
        simmer::activate('timer.') %>%
        simmer::activate('labile.')
    
    ## create a trajectory for the labile stage
    labile <- simmer::trajectory() %>%
        simmer::set_attribute('id', function() simmer::get_n_generated(s, 'labile.')) %>%
        simmer::set_attribute('cancelled', 0) %>%
        ## cancel if a new timer has started
        simmer::renege_if('cancel-labile',
                          out=simmer::trajectory() %>%
                              simmer::set_attribute('cancelled', 1)) %>%
        randomWalk(s, N_labile, t_labile, modulation=TRUE) %>%
        simmer::set_global('labile_id', function() simmer::get_attribute(s, 'id')) %>%
        ## start a new nonlabile program after completion
        simmer::activate('nonlabile.')
    
    ## create a trajectory for the nonlabile stage
    nonlabile <- simmer::trajectory() %>%
        simmer::set_attribute('id', function() simmer::get_global(s, 'labile_id')) %>%
        randomWalk(s, N_nonlabile, t_nonlabile, modulation=TRUE) %>%
        simmer::set_global('nonlabile_id', function() simmer::get_attribute(s, 'id')) %>%
        ## start a motor program after completion
        simmer::activate('motor.')
    
    ## create a trajectory for the motor planning stage
    motor <- simmer::trajectory() %>%
        simmer::set_attribute('id', function() simmer::get_global(s, 'nonlabile_id')) %>%
        randomWalk(s, N_motor, t_motor, modulation=FALSE) %>%
        simmer::set_global('motor_id', function() simmer::get_attribute(s, 'id')) %>%
        ## start a new saccade execution after completion
        simmer::activate('execution.')
    
    ## create a trajectory for the saccade execution stage
    execution <- simmer::trajectory() %>%
        simmer::send('end-fixation') %>%
        simmer::set_attribute('id', function() simmer::get_global(s, 'motor_id')) %>%
        randomWalk(s, N_execution, t_execution, modulation=FALSE) %>%
        simmer::set_global('execution_id', function() simmer::get_attribute(s, 'id')) %>%
        simmer::activate('fixation.')

    ## create a trajectory for fixations
    fixation <- simmer::trajectory() %>%
        simmer::set_attribute('id', function() simmer::get_global(s, 'execution_id')) %>%
        simmer::set_attribute('state', 0) %>%
        simmer::trap('end-fixation') %>%
        simmer::wait()
    
    s %>%
        simmer::add_generator('init.', simmer::trajectory() %>%
                                       simmer::set_global('modulation', modulation) %>%
                                       simmer::activate('timer.'), simmer::at(0)) %>%
        simmer::add_generator('timer.', timer, simmer::when_activated(), mon=mon) %>%
        simmer::add_generator('labile.', labile, simmer::when_activated(), mon=mon) %>%
        simmer::add_generator('nonlabile.', nonlabile, simmer::when_activated(), mon=mon) %>%
        simmer::add_generator('motor.', motor, simmer::when_activated(), mon=mon) %>%
        simmer::add_generator('execution.', execution, simmer::when_activated(), mon=mon) %>%
        simmer::add_generator('fixation.', fixation, simmer::when_activated(), mon=2)
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
        simmer::get_mon_attributes() %>%
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
#' @param ucm The UCM simmer environment to get cancellation information for
#' @param time if true (default to false), include the time of cancellation as a column
#' @return a dataframe detailing which IDs were cancelled, and how many cancellations were present before each uncancelled ID
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     get_cancellations()
#' @export
get_cancellations <- function(ucm, time=FALSE) {
    c <- ucm %>%
        simmer::get_mon_attributes() %>%
        tibble() %>%
        filter(key == 'cancelled') %>%
        separate(name, c('stage', 'n'), sep='\\.') %>%
        mutate(n=as.numeric(n)) %>%
        group_by(stage, n, replication, key)

    if (time==TRUE)
        c <- c %>%
            summarize(value=max(value), time=ifelse(value>0, max(time), NA), .groups='drop')
    else
        c <- c %>%
            summarize(value=max(value), .groups='drop')

    
    c %>%
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
        relocate(replication, id, cancelled)
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
        simmer::get_mon_attributes() %>%
        tibble() %>%
        filter(key == 'state') %>%
        separate(name, c('stage', 'n'), sep='\\.') %>%
        mutate(stage=factor(stage, levels=c('timer', 'labile', 'nonlabile', 'motor', 'execution', 'fixation'))) %>%
        pivot_wider(names_from=key)
    
    N <- s %>% group_by(stage, replication) %>%
        summarize(N=length(unique(state)))

    s %>% left_join(N, by=c('replication', 'stage')) %>%
        mutate(cum_state=ifelse(stage=='timer', state/N,
                         ifelse(stage=='labile', 1 + state/N,
                         ifelse(stage=='nonlabile', 2 + state/N,
                         ifelse(stage=='motor', 3 + state/N,
                         ifelse(stage=='execution', 4 + state/N,
                         ifelse(stage=='fixation', 5, state))))))) %>%
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
        simmer::get_mon_arrivals() %>%
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

#' Return a function to format numerics in seconds to milliseconds
ms_format <- function() {
    function(x) format(x*1000, digits=2) 
}

#' Format numbers by removing leading 0's
no_leading_zeros <- function (x, digits=2) {
    if (is.character(x))
        x <- as.numeric(x)
    
    str_replace(sprintf('%.2f', x), '0\\.', '\\.')
}

#' Get a traceplot of UCM's simulation from start to end.
#'
#' @param ucm The UCM simmer environment to plot traces for
#' @param start The start time of the plot (ms)
#' @param end The end time of the plot (ms)
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     trace_plot()
#' @export
trace_plot <- function(ucm, start=0, end=NULL, fix_areas=TRUE, cancellations=TRUE) {
    if (is.null(end)) {
        if (is.list(ucm))
            end <- simmer::now(ucm[[1]])
        else
            end <- simmer::now(ucm)
    }
    
    s <- ucm %>%
        get_states() %>%
        filter(time >= start-1 & time <= end+1)
    c <- ucm %>%
        get_cancellations(time=TRUE) %>%
        filter(cancelled == 1 & time >= start-1 & time <= end+1)
    
    ## add states for cancelled labile programs
    s.cancelled <- s %>%
        filter(stage=='labile' & id %in% c$id) %>%
        group_by(replication, id) %>%
        slice(n()) %>%
        ungroup() %>%
        mutate(time=c$time)
    s <- bind_rows(s, s.cancelled)
    
    p <- s %>%
        ggplot(aes(x=time, y=cum_state, group=interaction(replication, id)))

    if (fix_areas) {
        f <- ucm %>%
            get_fixations() %>%
            filter(end_time >= start & start_time <= end)

        if (nrow(f) > 0)
            p <- p + geom_rect(aes(xmin=start_time, xmax=end_time, ymin=-Inf, ymax=Inf),
                               inherit.aes=FALSE, alpha=0.25, data=f)
    }
    
    p <- p + geom_hline(yintercept=0:5) +
        geom_step() +
        scale_x_continuous(name='Time (ms)', expand=c(0, 0), labels=ms_format()) +
        scale_y_continuous(name='', breaks=0.5:4.5,
                           labels=c('timer', 'labile', 'non-labile', 'motor', 'execution'),
                           limits=c(0, 5), expand=c(0, 0)) +
        coord_cartesian(xlim=c(start, end), expand=FALSE) +
        theme_bw() +
        theme(axis.text.y=element_text(angle=90, hjust=0.5, vjust=0.5),
              axis.title.y=element_blank())

    if (cancellations & nrow(s.cancelled) > 0)
        p <- p + geom_point(data=s.cancelled, color='red', shape=4, size=2, stroke=1) +
            geom_segment(aes(xend=time, yend=1), data=s.cancelled, color='red', linetype='dashed')
        
    p
}


#' Get a dataframe containing all of the fixations from the UCM,
#' including two extra columns:
#'   - next_id: the id of the trace ending the current fixation
#'   - prev_id: the id of the previous fixation, which the current fixation ended
#'
#' @param ucm the UCM simmer environment to get fixations from
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     get_aligned_fixations()
#' @export
get_aligned_fixations <- function(ucm) {
    ucm %>% get_fixations() %>%
        full_join(get_cancellations(ucm)) %>%
        mutate(next_id=id+n_cancellations+1) %>%
        group_by(replication) %>%
        complete(id=1:max(id)) %>%  ## fill in cancelled fixations
        mutate(c=cumsum(!cancelled) + cancelled) %>%
        group_by(replication, c=cumsum(!cancelled)) %>%
        mutate(next_id=ifelse(is.na(n), next_id[1], next_id)) %>%
        group_by(replication) %>%
        mutate(next_start_time=map_dbl(next_id, ~ ifelse(!is.na(.), start_time[id == .], NA)),
               next_end_time=map_dbl(next_id, ~ ifelse(!is.na(.), end_time[id == .], NA)),
               next_duration=map_dbl(next_id, ~ ifelse(!is.na(.), duration[id == .], NA)),
               prev_id=map_dbl(id, ~ ifelse(. %in% next_id, id[next_id==.], NA))) %>%
        group_by(replication, c=cumsum(!cancelled)+cancelled) %>%
        mutate(prev_id=ifelse(is.na(n), prev_id[n()], prev_id)) %>%
        group_by(replication) %>%
        mutate(prev_start_time=map_dbl(prev_id, ~ ifelse(!is.na(.), start_time[id == .], NA)),
               prev_end_time=map_dbl(prev_id, ~ ifelse(!is.na(.), end_time[id == .], NA)),
               prev_duration=map_dbl(prev_id, ~ ifelse(!is.na(.), duration[id == .], NA))) %>%
        select(-c)
}

#' Get a trace of UCM's simulation, where each time is aligned
#' to the start of the previous fixation. Most useful to show the
#' timecourse of model behavior leading to the distribution of
#' fixation durations.
#'
#' @param ucm the UCM simmer environment to plot traces for
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     get_aligned_states()
#' @export
get_aligned_states <- function(ucm) {
    get_states(ucm) %>%
        filter((stage != 'execution' | state == 0) & stage != 'fixation') %>%
        left_join(ucm %>% get_aligned_fixations() %>% select(-n)) %>%
        mutate(time=time - prev_start_time)
}

#' Get a trace of the UCM's simulation, where each time is aligned
#' to the start of the timer on each iteration. Most useful to show the
#' duration of each stage of processing
get_timer_aligned_states <- function(ucm) {
    ## get the state traces aligned by the timer's start time
    get_states(ucm) %>%
        group_by(replication, id) %>%
        mutate(time=time-time[1]) %>%
        left_join(get_cancellations(ucm)) %>%
        group_by(replication)
}

#' Get a traceplot of UCM's simulation, where each trace is aligned
#' to the start of the previous fixation. Most useful to show the
#' timecourse of model behavior leading to the distribution of
#' fixation durations.
#'
#' @param ucm the UCM simmer environment to plot traces for
#' @param aligned_states optional argument to use cached state data
#' @param n the number of traces to plot (by default, plot all traces)
#' @param ids a list of individual fixation ids to plot (by default, plot all traces)
#' @param cancelled if true, plot traces of cancelled saccade programs. if
#' false, only plot the saccade programs which exited the labile stage
#' @examples
#' UCM() %>%
#'     run(until=5) %>%
#'     aligned_trace_plot()
#' @export
aligned_trace_plot <- function(ucm, aligned_states=NULL, n=NULL, fix_ids=NULL, cancelled=FALSE) {
    s <- aligned_states
    if (is.null(s))
        s <- get_aligned_states(ucm)
    
    if (!cancelled)
        s <- s %>% filter(!cancelled)

    if (!is.null(fix_ids)) {
        if (is.numeric(fix_ids)) {
            s <- s %>% group_by(replication) %>%
                nest %>%
                mutate(data=map(data, ~ filter(., prev_id %in% fix_ids))) %>%
                unnest(data)
        } else if (is.data.frame(fix_ids)) {
            s <- s %>% group_by(replication) %>%
                nest() %>%
                mutate(data=map(data,
                                ~ filter(., prev_id %in%
                                            fix_ids$id[fix_ids$replication == replication[1]]))) %>%
                unnest(data)
        }
    } else if (!is.null(n)) {
        if (n <= 0)
            stop('Error: N must be an integer above 0')

        ## select the fixations to keep
        s <- s %>% mutate(ID=paste0(replication, '_', id))
        fix_ids <- s %>% group_by(replication) %>%
            distinct(ID) %>%
            slice_sample(n=n) %>%
            pull(ID)
        
        s <- s %>% group_by(replication) %>%
            nest %>%
            mutate(data=map(data, ~filter(., ID %in% fix_ids))) %>%
            unnest(data) %>%
            select(-ID)
    }
    
    s %>%
        ggplot(aes(x=time, y=cum_state, color=interaction(id, replication))) +
        geom_hline(yintercept=0:4) +
        geom_step() +
        scale_x_continuous(name='Time (ms)', labels=ms_format()) +
        scale_y_continuous(name='', breaks=0.5:3.5,
                           labels=c('timer', 'labile', 'non-labile', 'motor'),
                           limits=c(0, 4), expand=c(0, 0)) +
        theme_bw() +
        theme(axis.text.y=element_text(angle=90, hjust=0.5, vjust=0.5),
              axis.title.y=element_blank())
}


cancellations_hist <- function(ucm, max_duration=1.2, binwidth=0.06, max_cancellations=10) {
    c <- get_cancellations(ucm)
    f <- left_join(get_fixations(ucm), c)
    
    f %>%
        filter(duration <= max_duration) %>%
        mutate(n_cancellations=factor(pmin(n_cancellations, max_cancellations),
                                      levels=0:max_cancellations,
                                      labels=c(0:(max_cancellations-1), paste0(max_cancellations, '+')))) %>%
        ggplot(aes(x=duration, group=n_cancellations, fill=n_cancellations)) +
        scale_fill_viridis(name='Cancellations', discrete=TRUE) +
        geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=binwidth) +
        coord_cartesian(xlim=c(0, max_duration)) +
        xlab('Fixation Duration (s)') + ylab('Proportion') +
        theme_bw()
}

#' For a given set of UCM parameters, calculate the cancellation rate
#'
#' @param N_timer the number of states in the timer
#' @param N_labile the number of states in the labile stage
#' @param t_timer the mean duration of the timer
#' @param t_labile the mean duration of the labile stage
cancellation_prob <- function(N_timer, N_labile, t_timer, t_labile) {
    r_timer <- N_timer/t_timer
    r_labile <- N_labile/t_labile

    return(pbeta(r_timer / (r_timer+r_labile), N_timer, N_labile))
}
