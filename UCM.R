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
#' @param modulation The starting modulation of the timer, labile stage, and nonlabile stage.
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
                modulation=1.0) {
    if (!all(c(N_timer, N_labile, N_nonlabile, N_motor, N_execution) > 0))
        stop('Cannot create UCM with N < 0')
    if (!all(c(t_timer, t_labile, t_nonlabile, t_motor, t_execution) > 0))
        stop('Cannot create UCM with t < 0')
    if (modulation <= 0 | modulation > 1)
        stop('Cannot create UCM with modulation <= 0 or > 1')
    
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
                               set_global('modulation', modulation) %>%
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
        get_mon_attributes() %>%
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
        get_mon_attributes() %>%
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
trace_plot <- function(ucm, start=0, end=NULL, fix_areas=TRUE, cancellations=TRUE) {
    if (is.null(end)) {
        if (is.list(ucm))
            end <- now(ucm[[1]])
        else
            end <- now(ucm)
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
        scale_x_continuous(name='Time (s)', expand=c(0, 0)) +
        scale_y_continuous(name='', breaks=0.5:4.5,
                           labels=c('timer', 'labile', 'non-labile', 'motor', 'execution'),
                           limits=c(0, 5), expand=c(0, 0)) +
        coord_cartesian(xlim=c(start, end), expand=FALSE) +
        theme_bw() +
        theme(axis.text.y=element_text(angle=90, hjust=0.5, vjust=0.5),
              axis.title.y=element_blank())

    if (cancellations & nrow(s.cancelled > 0))
        p <- p + geom_point(data=s.cancelled, color='red', shape=4, size=2, stroke=1) +
            geom_segment(aes(xend=time, yend=1), data=s.cancelled, color='red', linetype='dashed')
        
    p
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
    ## Get the model fixations, and ID them by the ID
    ## of the saccade program that *ended* that fixation
    fix <- get_fixations(ucm) %>%
        group_by(replication) %>%
        mutate(end_id=lead(id),
               start_id=id,
               id=end_id) %>%
        relocate(id, start_id, end_id)

    ## fill in blank end_ids
    for (r in 1:max(fix$replication)) {
        i <- max(which(fix$replication==r))
        fix$end_id[i] <- get_states(ucm) %>%
            filter(id > fix$start_id[i] & state==0 & stage=='execution' & time == fix$end_time[i]) %>%
            pull(id)
        fix$id[i] <- fix$end_id[i]
    }
    
    ## Get fixation cancellation information
    f <- get_cancellations(ucm) %>%
        left_join(fix) %>%
        select(-n_cancellations, -n)

    for (r in 1:max(f$replication)) {
        fr <- f %>% filter(replication==r)
        
        ## remove unstarted/unfinished fixations
        while (is.na(fr$start_id[1]) & !fr$cancelled[1])
            fr <- fr[-1,]
        while (is.na(fr$start_id[nrow(fr)]))
            fr <- fr[-nrow(fr),]
        
        ## find the start_id and end_id of cancelled saccade programs
        for (i in 1:nrow(fr)) {
            if (fr$cancelled[i]) {
                j <- i+1
                while (fr$cancelled[j] & j <= nrow(fr))
                    j <- j+1
                if (j <= nrow(fr)) {
                    fr$start_id[i] <- fr$start_id[j]
                    fr$end_id[i] <- fr$end_id[j]
                    fr$start_time[i] <- fr$start_time[j]
                    fr$end_time[i] <- fr$end_time[j]
                    fr$duration[i] <- fr$duration[j]            
                }
            }
        }

        f <- f %>% filter(replication != r) %>%
            bind_rows(fr)
    }
    
    ## count number of cancellations
    f <- f %>%
        group_by(replication, end_id) %>%
        mutate(n_cancellations=n()-1)
    
    
    ## get the model trace and align by start_time
    get_states(ucm) %>%
        left_join(f) %>%
        filter(!is.na(start_id)) %>%
        mutate(time=time-start_time)
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
    s <- s %>% filter(stage != 'fixation')
    
    if (!cancelled)
        s <- s %>% filter(cancelled==0)

    if (!is.null(fix_ids)) {
        s <- s %>% group_by(replication) %>%
            nest %>%
            mutate(data=map(data,
                            ~filter(., end_id %in% fix_ids))) %>%
            unnest(data)
    } else if (!is.null(n)) {
        if (n <= 0)
            stop('Error: N must be an integer above 0')

        ## select the fixations to keep
        s <- s %>% mutate(END_ID=paste0(replication, '_', end_id))
        fix_ids <- s %>% group_by(replication) %>%
            distinct(END_ID) %>%
            slice_sample(n=n) %>%
            pull(END_ID)
        
        s <- s %>% group_by(replication) %>%
            nest %>%
            mutate(data=map(data, ~filter(.x, END_ID %in% fix_ids))) %>%
            unnest(data) %>%
            select(-END_ID)
    }
    
    s %>%
        ggplot(aes(x=time, y=cum_state, group=interaction(id, replication))) +
        geom_hline(yintercept=0:5) +
        geom_step() +
        scale_x_continuous(name='Time (s)', expand=c(0, 0)) +
        scale_y_continuous(name='', breaks=0.5:3.5,
                           labels=c('timer', 'labile', 'non-labile', 'motor'),
                           limits=c(0, 4), expand=c(0, 0)) +
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
