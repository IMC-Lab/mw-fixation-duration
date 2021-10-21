library(simmer)
library(tidyverse)
library(parallel)
library(viridis)
library(patchwork)

## Load the UCM source code
source('UCM.R')
options(mc.cores=parallel::detectCores())
set.seed(1234)


## helper function to run/plot UCM with various parameter settings
fixdur_traceplot <- function(ucm, n=NULL, xmin=-0.5, xmax=1.0, legend.x=.033,
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
        theme_void(base_size=16) +
        theme(legend.position=c(legend.x, .5), 
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())

    ## plot some of the model traces
    trace <- aligned_trace_plot(ucm, n=n) +
        ##geom_vline(xintercept=0) +
        aes(color=factor(replication)) +
        scale_x_continuous(name='Time relative to previous fixation onset (s)',
                           minor_breaks=c(0), limits=c(xmin, xmax), expand=c(0,0)) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size=16) +
        theme(legend.position='none',
              panel.grid.major.y=element_blank())
    
    ## overlay the plots
    ((dens/trace) &
     coord_cartesian(xlim=c(xmin, xmax))) +
        plot_layout(heights=c(.25, .75))
}


ucm5 <- UCM(N_timer=5, N_labile=5, N_nonlabile=5, N_motor=5, N_execution=5) %>%
    run(1000)
ucm25 <- UCM(N_timer=25, N_labile=25, N_nonlabile=25, N_motor=25, N_execution=25) %>%
    run(1000)
ucm <- list(ucm5, ucm25)

fixdur_traceplot(ucm, rep_name='Threshold', rep_labels=c('N = 5', 'N = 25'), n=15)
ggsave('plots/UCM-N-mod.png', width=16, height=8)



ucmNormal <- UCM() %>%
    run(1000)
ucmMod <- UCM(modulation=.7) %>%
    run(1000)
ucm <- list(ucmNormal, ucmMod)

fixdur_traceplot(ucm, rep_name='Modulation', rep_labels=c('Modulation = 1', 'Modulation = .7'), n=15, legend.x=.05)
ggsave('plots/UCM-rate-mod.png', width=16, height=8)
