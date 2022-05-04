library(simmer)
library(tidyverse)
library(parallel)
library(viridis)
library(patchwork)

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

## get cancellation times
ucm %>% get_cancellations(time=TRUE) %>%
    mutate(time=time - 1) %>%
    filter(cancelled == 1, time > 0, time < 3)


## helper function to run/plot UCM with various parameter settings
fixdur_traceplot <- function(ucm, aligned_states=NULL, n=NULL, xmin=-0.5, xmax=1.0, base_size=12,
                             bw='nrd0', rep_name='Replication',
                             rep_labels=as.character(seq(1, ifelse(is.list(ucm), length(ucm), 1)))) {
    ## plot the hist of fixation durations
    dens <- ucm %>% get_fixations %>%
        mutate(replication=factor(replication, labels=rep_labels)) %>%
        ggplot(aes(x=duration, group=replication,
                   color=replication, fill=replication)) +
        geom_density(alpha=0.33, bw=bw) +
        scale_x_continuous(limits=c(0, xmax), expand=c(0,0)) +
        scale_color_viridis(discrete=TRUE, name=rep_name) +
        scale_fill_viridis(discrete=TRUE, name=rep_name) +
        ggtitle('Fixation Duration (s)') +
        theme_void(base_size=base_size) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())

    ## plot some of the model traces
    trace <- aligned_trace_plot(ucm, aligned_states=aligned_states, n=n) +
        ##geom_vline(xintercept=0) +
        aes(color=factor(replication)) +
        scale_x_continuous(name='Time relative to previous fixation onset (s)',
                           minor_breaks=c(0), limits=c(xmin, xmax), expand=c(0,0)) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size=base_size) +
        theme(legend.position='none',
              panel.grid.major.y=element_blank())
    
    ## overlay the plots
    ((dens/trace) &
     coord_cartesian(xlim=c(xmin, xmax))) +
        plot_layout(heights=c(.25, .75), guides='collect')
}


ucmNormal <- UCM() %>%
    run(1000)
ucmMod <- UCM(modulation=.7) %>%
    run(1000)
sRateMod <- get_aligned_states(list(ucmNormal, ucmMod))

plot.rate_mod <- fixdur_traceplot(list(ucmNormal, ucmMod),
                                  aligned_states=sRateMod,
                                  n=5, base_size=10,
                                  rep_name='Modulation',
                                  rep_labels=c('M = 1', 'M = .7'))
plot.rate_mod
ggsave('plots/UCM-rate-mod.png', width=6, height=2)



ucm5 <- UCM(N_timer=5, N_labile=5, N_nonlabile=5, N_motor=5, N_execution=5) %>%
    run(1000)
sNMod <- get_aligned_states(list(ucmNormal, ucm5))

plot.N_mod <- fixdur_traceplot(list(ucmNormal, ucm5),
                               aligned_states=sNMod,
                               n=5, base_size=10, rep_name='Threshold',
                               rep_labels=c('N = 14', 'N = 5'))
plot.N_mod
ggsave('plots/UCM-N-mod.png', width=6, height=2)



wrap_plots(plot.rate_mod, plot.N_mod, ncol=1, tag_level='new') +
    plot_annotation(tag_levels=list(c('A', '', 'B', '')))
ggsave('plots/UCM-mod2.png', width=6, height=5)
