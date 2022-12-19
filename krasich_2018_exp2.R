library(simmer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(parallel)
library(ParBayesianOptimization)

## Load UCM source code
source('UCM.R')
options(mc.cores=parallel::detectCores())

## Load data
fix_data <- read_excel('data/BE13Triad_datat.xlsx') %>%
    mutate(fixdur=fixdur/1000,
           proberesp=factor(proberesp, levels=c('NotMW', 'MW')),
           mw=ifelse(proberesp == 'MW', 1, 0))


m.out <- lmer(fixdur ~ TimeIntoImage + proberesp + (1|ID), fix_data)
m.med <- lmer(TimeIntoImage ~ proberesp + (1|ID), fix_data)
med <- mediate(m.med, m.out, treat='proberesp', mediator='TimeIntoImage')
summary(med)


## Set optimization parameters
bounds <- list(N_states=c(5L, 30L),
               t_timer=c(.1, .3),
               t_labile=c(.1, .25),
               t_nonlabile=c(.05, .1))


duration.1 <- function () runif(1, 45, 75)

function LL(N_states, t_timer, t_labile, t_nonlabile) {
    fix <- 1:N_trials %>%
        mclapply(function (i) {
            UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                N_motor=N_states, N_execution=N_states,
                t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile) %>%
                run(duration.1()) %>%
                wrap()
        }) %>%
        get_fixations()
    
    list(Score=do.call(LL_discrete, list(x, fix$activity_time, ...)))
}

bopt <- bayesOpt(LL, bounds, 'ucm.rds', initPoints=16,
                 iters.n=10, kappa=5.0, plotProgress=TRUE)

ucm.MLE <- getBestPars(bopt)

ucm.best <- mclapply(1:1000, function(i) {
    UCM(N_timer=ucm.MLE$N_states, N_labile=ucm.MLE$N_states,
        N_nonlabile=ucm.MLE$N_states, N_motor=ucm.MLE$N_states,
        N_execution=ucm.MLE$N_states, t_timer=ucm.MLE$t_timer,
        t_labile=ucm.MLE$t_labile, t_nonlabile=ucm.MLE$t_nonlabile) %>%
        run(30) %>%
        wrap()
})


f <- get_fixations(ucm.best)

after_stat(count / sum(count))

rbind(data.frame(type='Data', fixdur=df.1$fixdur),
      data.frame(type='UCM', fixdur=f$activity_time)) %>%
    ggplot(aes(x=fixdur)) +
    geom_histogram(aes(y=stat(density)), binwidth=0.05) +
    facet_grid(type ~ .) +
    scale_x_continuous(limits=c(0, 3)) +
    xlab('Fixation Duration (s)') + ylab('Proportion') +
    theme_bw()







fix_data %>%
    ggplot(aes(x=fixdur, color=factor(mw))) +
    stat_ecdf() +
    scale_linetype_manual(name='', values=c('solid', 'dotted', 'dashed')) +
    scale_color_jco(name='Probe Response', labels=c('On-task', 'Mind-wandering')) +
    xlab('Fixation Duration (s)') + ylab('Cumulative Probability') +
    coord_cartesian(xlim=c(0, 2)) +
    theme_bw() + ggtitle('Krasich et al. (2018) Experiment 2')
ggsave('plots/reciprobit_krasich_2018_exp2.png', width=8, height=4)
