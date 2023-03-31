library(simmer)
library(tidyverse)
library(parallel)
library(viridis)
library(patchwork)
library(latex2exp)

## Load the UCM source code
source('UCM.R')
options(mc.cores=parallel::detectCores())

## Make a sample trace plot
set.seed(12345678)
ucm <- UCM() %>%
    run(until=5)
trace_plot(ucm, start=1, end=3) +
    scale_x_continuous(name='Time (s)', labels=c('0', '.5', '1', '1.5', '2'))
ggsave('plots/UCM-trace2.png', width=7, height=3.5)
ggsave('plots/UCM-trace2.pdf', width=7, height=3.5, dpi=1000)

## get cancellation times
ucm %>% get_cancellations(time=TRUE) %>%
    mutate(time=time - 1) %>%
    filter(cancelled == 1, time > 0, time < 3)

## split up fixations by quantiles
add_quantiles <- function(fix, n_quantiles=6) {
    fix %>%
        group_by(replication) %>%
        mutate(quantile=map_dbl(duration,
                                function(d) detect_index(quantile(duration, probs=seq(0, 1, length.out=n_quantiles), na.rm=TRUE), ~ d <= .)) - 1)
}

## sample n_samples ids in each of n_quantiles of fixation durations
sample_quantiles <- function(ucm, n_quantiles=6, n_samples=1) {
    ucm %>%
        get_fixations() %>%
        add_quantiles() %>%
        filter(quantile > 0) %>%
        group_by(replication, quantile) %>%
        sample_n(n_samples) %>%
        select(quantile, replication, id)
}

## helper function to run/plot UCM with various parameter settings
fixdur_traceplot <- function(ucm, n=NULL, fix_ids=NULL, xmin=-.5, xmax=1, base_size=10,
                             bw='nrd0', title='Fixation Duration (s)', rep_name='Replication',
                             rep_labels=as.character(seq(1, ifelse(is.list(ucm), length(ucm), 1)))) {
    ## plot the hist of fixation durations
    dens <- ucm %>% get_fixations %>%
        mutate(replication=factor(replication)) %>%
        ggplot(aes(x=duration, group=replication,
                   color=replication, fill=replication)) +
        geom_density(alpha=0.33, bw=bw) +
        scale_x_continuous(limits=c(0, xmax), expand=c(0,0)) +
        scale_color_viridis(discrete=TRUE, name=rep_name, labels=rep_labels) +
        scale_fill_viridis(discrete=TRUE, name=rep_name, labels=rep_labels) +
        ggtitle(title) +
        theme_void(base_size=base_size) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())

    ## plot some of the model traces
    trace <- aligned_trace_plot(ucm, n=n, fix_ids=fix_ids) +
        ##geom_vline(xintercept=0) +
        aes(color=factor(replication), group=id) +
        scale_x_continuous(name='Time relative to previous fixation onset (ms)',
                           minor_breaks=c(0), expand=c(0,0), labels=ms_format()) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size=base_size) +
        theme(legend.position='none',
              panel.grid.major.y=element_blank())
    
    ## overlay the plots
    ((dens/trace) &
     coord_cartesian(xlim=c(xmin, xmax))) +
        plot_layout(heights=c(.25, .75), guides='collect')
}

tracedur_traceplot <- function(ucm, n=NULL, fix_ids=NULL, xmin=0, xmax=1.25, base_size=10,
                               bw='nrd0', title='Fixation Duration (s)', rep_name='Replication',
                               rep_labels=as.character(seq(1, ifelse(is.list(ucm), length(ucm), 1)))) {
    s <- get_timer_aligned_states(ucm) %>%
        mutate(prev_id=id) %>%
        ungroup()
    
    ## plot the hist of fixation durations
    dens <- s %>% filter(stage=='fixation') %>%
        mutate(replication=factor(replication)) %>%
        ggplot(aes(x=time, group=replication,
                   color=replication, fill=replication)) +
        geom_density(alpha=0.33, bw=bw) +
        scale_x_continuous(limits=c(0, xmax), expand=c(0,0)) +
        scale_color_viridis(discrete=TRUE, name=rep_name, labels=rep_labels) +
        scale_fill_viridis(discrete=TRUE, name=rep_name, labels=rep_labels) +
        ggtitle(title) +
        theme_void(base_size=base_size) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
    
    ## plot some of the model traces
    trace <- aligned_trace_plot(ucm, aligned_states=s, n=n, fix_ids=fix_ids) +
        aes(group=id, color=factor(replication)) +
        scale_x_continuous(name='Time relative to start of the timer (ms)',
                           minor_breaks=c(0), expand=c(0,0), labels=ms_format()) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size=base_size) +
        theme(legend.position='none',
              panel.grid.major.y=element_blank())
    
    ## overlay the plots
    ((dens/trace) &
     coord_cartesian(xlim=c(xmin, xmax))) +
        plot_layout(heights=c(.25, .75), guides='collect')
}

set.seed(1)
ucm <- UCM() %>%
    run(1000)
ucm.rate_mod <- UCM(modulation=.7) %>%
    run(1000)
ucm.n_mod <- UCM(N_timer=5, N_labile=5, N_nonlabile=5) %>%
    run(1000)


## plot rate mod
ids.rate_mod <- sample_quantiles(list(ucm, ucm.rate_mod))
plot.rate_mod <- fixdur_traceplot(list(ucm, ucm.rate_mod), fix_ids=ids.rate_mod,
                                  title='Fixation Duration (ms)',
                                  rep_name='Rate\nModulation\nParameter',
                                  rep_labels=TeX(c(expression('$\\beta_{rate} = \\,1$'),
                                                   expression('$\\beta_{rate} = .7$'))))
plot.rate_mod

ggsave('plots/UCM-rate-mod.png', width=6, height=2)

## plot n mod
ids.n_mod <- ids.rate_mod %>% filter(replication==1) %>%
    bind_rows(sample_quantiles(ucm.n_mod) %>% mutate(replication=2))
plot.n_mod <- fixdur_traceplot(list(ucm, ucm.n_mod), fix_ids=ids.n_mod,
                               title='Fixation Duration (ms)',
                               rep_name='Threshold\nModulation\nParameter',
                               rep_labels=TeX(c(expression('$\\beta_{threshold} = \\,1\\;$'),
                                                expression('$\\beta_{threshold} = .36$'))))
plot.n_mod

ggsave('plots/UCM-N-mod.png', width=6, height=2)

## combine the plots
wrap_plots(plot.rate_mod, plot.n_mod, ncol=1, tag_level='new') +
    plot_annotation(tag_levels=list(c('A', '', 'B', '')))
ggsave('plots/UCM-mod2.png', width=6, height=5)



## make plots aligned by timer
plot.rate_mod.2 <- tracedur_traceplot(list(ucm, ucm.rate_mod), fix_ids=ids.rate_mod,
                                      title='Visual Processing Duration (ms)',
                                      rep_name='Rate\nModulation\nParameter',
                                      rep_labels=TeX(c(expression('$\\beta_{rate} = \\,1$'),
                                                   expression('$\\beta_{rate} = .7$'))))
plot.n_mod.2 <- tracedur_traceplot(list(ucm, ucm.n_mod), fix_ids=ids.n_mod,
                                   title='Visual Processing Duration (ms)',
                                   rep_name='Threshold\nModulation\nParameter',
                                   rep_labels=TeX(c(expression('$\\beta_{threshold} = \\,1\\;$'),
                                                    expression('$\\beta_{threshold} = .36$'))))
wrap_plots(plot.rate_mod.2, plot.n_mod.2, ncol=1, tag_level='new') +
    plot_annotation(tag_levels=list(c('A', '', 'B', '')))
ggsave('plots/UCM-mod3.png', width=6, height=5)


