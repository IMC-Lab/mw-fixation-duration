library(readr)
library(simmer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)
library(ggsci)
library(parallel)
library(ParBayesianOptimization)
library(scales)
library(ggallin)


## Load UCM source code
source('UCM.R')
source('UCM-optim.R')
options(mc.cores=parallel::detectCores())

## Set optimization parameters
z.T_MOTOR <- .03            ## duration of motor planning stage (fixed parameter)
z.T_EXECUTION <- .04        ## duration of saccade execution (fixed parameter)
z.N_TRIALS <- 2500          ## number of UCM trials to simulate per epoch to estimate LL
z.MIN_FIXDUR <- 0.08        ## only keep fixations over this limit
z.MAX_FIXDUR <- 2           ## only keep fixations under this limit
z.BINWIDTH <- 0.04          ## the width of histogram bins to calculate LL
z.TRIAL_DUR <- 10           ## trial length (in seconds)
z.bounds.separate <- list(N_states=c(2L, 30L),
                          t_timer=c(.15, .375),
                          t_labile=c(.1, .225),
                          t_nonlabile=c(.025, .08))
z.bounds.rate_mod <- c(z.bounds.separate, list(modulation=c(0.25, 1)))
z.bounds.N_mod <- c(z.bounds.rate_mod, list(N_mod=c(0.25, 1)))

## Load data
z.fix_data <- read_csv('https://osf.io/fpqsj/download') %>%
    tibble %>%
    mutate(fixdur=duration/1000,
           attention=factor(attention, levels=c('On-task', 'MW: Unintentional', 'MW: Intentional')),
           Attention=factor(Attention, levels=c('On-task', 'Mind-wandering')),
           mw=ifelse(Attention == 'Mind-wandering', 1, 0),
           type='Data') %>%
    rename(proberesp=Attention,
           ID=subject_nr) %>%
    ## filter as per the manuscript
    drop_na() %>%
    filter(fixdur >= z.MIN_FIXDUR & fixdur <= z.MAX_FIXDUR) %>%
    filter(x_pos > 0 & x_pos < 1024) %>%
    filter(y_pos > 0 & y_pos < 768) %>%
    filter(attention != 'MW: Intentional')

## count number of trials with/without MW
z.fix_data %>% group_by(ID, s_name, proberesp) %>%
    summarize %>% group_by(proberesp) %>% count

z.fix_data.on_task <- z.fix_data %>% filter(mw == 0)
z.fix_data.mw <- z.fix_data %>% filter(mw == 1)

## Plot data
mw_plot('Zhang et al. (2021) Raw Data', z.fix_data, z.BINWIDTH)
ggsave('plots/hist_zhang_2021.png', width=8, height=6)


## Define log likelihood functions
z.LL.on_task <- LL(z.fix_data.on_task$fixdur, z.bounds.separate,
                 n_trials=z.N_TRIALS, trial_dur=z.TRIAL_DUR, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR,
                 binwidth=z.BINWIDTH, default_params=list(t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION))
z.LL.mw <- LL(z.fix_data.mw$fixdur, z.bounds.separate,
            n_trials=z.N_TRIALS, trial_dur=z.TRIAL_DUR, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR,
            binwidth=z.BINWIDTH, default_params=list(t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION))
z.LL.rate_mod <- LL.joint(z.fix_data.on_task$fixdur, z.fix_data.mw$fixdur, z.bounds.rate_mod,
                        n_trials=z.N_TRIALS, trial_dur=z.TRIAL_DUR, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR,
                        binwidth=z.BINWIDTH, default_params=list(t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION))
z.LL.N_mod <- LL.joint(z.fix_data.on_task$fixdur, z.fix_data.mw$fixdur, z.bounds.N_mod,
                     n_trials=z.N_TRIALS, trial_dur=z.TRIAL_DUR, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR,
                     binwidth=z.BINWIDTH, default_params=list(t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION))



## Test for differences between MW/Attentive Viewing using permuted KS test
z.KS <- KS.permutation(z.fix_data, n=10000)
z.KS$D[1]
mean(z.KS$D.null >= z.KS$D)
ggplot(z.KS, aes(x=D.null)) +
    geom_histogram(bins=100) +
    geom_vline(xintercept=z.KS$D) +
    theme_classic()



##################################################################################################
##                      Optimze using rate modulation & N modulation
##################################################################################################
z.bopt.N_mod <- optimize('zhang_2021_N_mod_constrained.rds', z.LL.N_mod, z.bounds.N_mod)
z.params.N_mod <- getBestPars(z.bopt.N_mod)
z.params.N_mod

## Simulate lots of fixations from the best parameter values
z.ucm.on_task.N_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=z.params.N_mod$N_states, N_labile=z.params.N_mod$N_states,
        N_nonlabile=z.params.N_mod$N_states, N_motor=z.params.N_mod$N_states,
        N_execution=z.params.N_mod$N_states, t_timer=z.params.N_mod$t_timer,
        t_labile=z.params.N_mod$t_labile, t_nonlabile=z.params.N_mod$t_nonlabile,
        t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION) %>%
        run(z.TRIAL_DUR) %>%
        wrap()
})
z.fix_ucm.on_task.N_mod <- z.ucm.on_task.N_mod %>%
    get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= z.MIN_FIXDUR & fixdur <= z.MAX_FIXDUR)
z.ucm.mw.N_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=round(z.params.N_mod$N_states * z.params.N_mod$N_mod),
        N_labile=round(z.params.N_mod$N_states * z.params.N_mod$N_mod),
        N_nonlabile=round(z.params.N_mod$N_states * z.params.N_mod$N_mod),
        N_motor=z.params.N_mod$N_states, N_execution=z.params.N_mod$N_states,
        t_timer=z.params.N_mod$t_timer,
        t_labile=z.params.N_mod$t_labile, t_nonlabile=z.params.N_mod$t_nonlabile,
        t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION, modulation=z.params.N_mod$modulation) %>%
        run(z.TRIAL_DUR) %>%
        wrap()
})
z.fix_ucm.mw.N_mod <- z.ucm.mw.N_mod %>%
    get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= z.MIN_FIXDUR & fixdur <= z.MAX_FIXDUR)

## Extract the model fixations
z.fix_ucm.N_mod <- bind_rows(z.fix_ucm.on_task.N_mod, z.fix_ucm.mw.N_mod)
z.fix.N_mod <- bind_rows(z.fix_data %>% select(type, mw, fixdur),
                         z.fix_ucm.N_mod %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, z.MIN_FIXDUR, z.MAX_FIXDUR, z.BINWIDTH)))

## calculate descriptive statistics
z.fix.N_mod %>%
    mutate(mean=map_dbl(data, ~mean(.x$fixdur*1000)),
           sd=map_dbl(data, ~sd(.x$fixdur*1000)),
           lower=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymin),
           upper=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymax)) %>%
    select(-data, -hist) %>%
    as.data.frame

## plot model fit (histograms/means)
z.hist.N_mod <- mw_hist(z.fix.N_mod)
z.means.N_mod <- mw_means_plot(z.fix.N_mod)
(z.hist.N_mod / z.means.N_mod & coord_cartesian(xlim=c(0, 1))) +
        plot_layout(heights=c(1, .25), guides='collect') +
        plot_annotation(title='Zhang et al. (2021)')
ggsave('plots/fit_zhang_2021_N_mod.png', width=8, height=4)

## plot model fit (ecdfs/ecdf_differences)
z.ecdf.N_mod <- mw_ecdf_plot(z.fix.N_mod)
z.ecdf_diff.N_mod <- mw_ecdf_difference_plot(z.fix.N_mod, n=500)
(z.ecdf.N_mod / z.ecdf_diff.N_mod) +
    plot_annotation(title='Zhang et al. (2021)')
ggsave('plots/ecdf_zhang_2021_N_mod.png', width=8, height=8)

## Calculate log-likelihood
z.ll.N_mod <- LL_discrete(z.fix_data.on_task$fixdur, z.fix_ucm.on_task.N_mod$fixdur, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR, binwidth=z.BINWIDTH, delta=1/z.N_TRIALS) +
    LL_discrete(z.fix_data.mw$fixdur, z.fix_ucm.mw.N_mod$fixdur, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR, binwidth=z.BINWIDTH, delta=1/z.N_TRIALS)
z.ll.N_mod




(wrap_plots(k.hist.N_mod+ggtitle('Krasich et al. (2018)')+coord_cartesian(xlim=c(0, 1), ylim=c(0, .18)),
            z.hist.N_mod+ggtitle('Zhang et al. (2021)')+coord_cartesian(xlim=c(0, 1), ylim=c(0, .18)),
            k.means.N_mod+coord_cartesian(xlim=c(0, 1)),
            z.means.N_mod+coord_cartesian(xlim=c(0, 1))) &
 theme(legend.position='bottom')) +
    plot_annotation(tag_levels=list(c('A', 'B', '', '')), title='Full Model') +
    plot_layout(heights=c(1, .33), guides='collect')
ggsave('plots/fit_Nmod.png', width=9, height=4.5)


((k.ecdf.N_mod+ggtitle('Krasich et al. (2018)')+coord_cartesian(xlim=c(0, 1)) |
  z.ecdf.N_mod+ggtitle('Zhang et al. (2021)')+coord_cartesian(xlim=c(0, 1))) +
 plot_layout(guides='collect')) /
    ((k.ecdf_diff.N_mod+coord_cartesian(xlim=c(0, 1), ylim=c(-.075, .03)) |
      z.ecdf_diff.N_mod+coord_cartesian(xlim=c(0, 1), ylim=c(-.075, .03))) +
     plot_layout(guides='collect')) +
    plot_annotation(tag_levels=list(c('A', 'B', '', '')), title='Full Model')
ggsave('plots/ecdf_Nmod.png', width=10, height=6)


## get cancellation rates
z.ucm.on_task.N_mod %>%
    get_cancellations() %>%
    summarize(M=mean(cancelled))

cancellation_prob(z.params.N_mod$N_states, z.params.N_mod$N_states, ## derived value
                  z.params.N_mod$t_timer, z.params.N_mod$t_labile)


z.ucm.mw.N_mod %>%
    get_cancellations() %>%
    summarize(M=mean(cancelled))

cancellation_prob(round(z.params.N_mod$N_states * z.params.N_mod$N_mod),  ## derived value
                  round(z.params.N_mod$N_states * z.params.N_mod$N_mod),
                  z.params.N_mod$t_timer, z.params.N_mod$t_labile)



##################################################################################################
##                                   Optimze w/ rate modulation only
##################################################################################################
z.bopt.rate_mod <- optimize('zhang_2021_rate_mod.rds', z.LL.rate_mod, z.bounds.rate_mod)
z.params.rate_mod <- getBestPars(z.bopt.rate_mod)
z.params.rate_mod

z.fix_ucm.on_task.rate_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=z.params.rate_mod$N_states, N_labile=z.params.rate_mod$N_states,
        N_nonlabile=z.params.rate_mod$N_states, N_motor=z.params.rate_mod$N_states,
        N_execution=z.params.rate_mod$N_states, t_timer=z.params.rate_mod$t_timer,
        t_labile=z.params.rate_mod$t_labile, t_nonlabile=z.params.rate_mod$t_nonlabile,
        t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION) %>%
        run(z.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= z.MIN_FIXDUR & fixdur <= z.MAX_FIXDUR)
z.fix_ucm.mw.rate_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=z.params.rate_mod$N_states, N_labile=z.params.rate_mod$N_states,
        N_nonlabile=z.params.rate_mod$N_states, N_motor=z.params.rate_mod$N_states,
        N_execution=z.params.rate_mod$N_states, t_timer=z.params.rate_mod$t_timer,
        t_labile=z.params.rate_mod$t_labile, t_nonlabile=z.params.rate_mod$t_nonlabile,
        t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION, modulation=z.params.rate_mod$modulation) %>%
        run(z.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= z.MIN_FIXDUR & fixdur <= z.MAX_FIXDUR)

## extract model fixations
z.fix_ucm.rate_mod <- bind_rows(z.fix_ucm.on_task.rate_mod, z.fix_ucm.mw.rate_mod)
z.fix.rate_mod <- bind_rows(z.fix_data %>% select(type, mw, fixdur),
                            z.fix_ucm.rate_mod %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, z.MIN_FIXDUR, z.MAX_FIXDUR, z.BINWIDTH)))

## calculate descriptive statistics
z.fix.rate_mod %>%
    mutate(mean=map_dbl(data, ~mean(.x$fixdur*1000)),
           sd=map_dbl(data, ~sd(.x$fixdur*1000)),
           lower=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymin),
           upper=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymax)) %>%
    select(-data, -hist) %>%
    as.data.frame

## plot model fit (histograms/means)
z.hist.rate_mod <- mw_hist(z.fix.rate_mod)
z.means.rate_mod <- mw_means_plot(z.fix.rate_mod)
(z.hist.rate_mod / z.means.rate_mod & coord_cartesian(xlim=c(0, 1))) +
        plot_layout(heights=c(1, .25), guides='collect') +
        plot_annotation(title='Zhang et al. (2021)')
ggsave('plots/fit_zhang_2021_rate_mod.png', width=8, height=4)

## plot model fit (ecdfs/ecdf_differences)
z.ecdf.rate_mod <- mw_ecdf_plot(z.fix.rate_mod)
z.ecdf_diff.rate_mod <- mw_ecdf_difference_plot(z.fix.rate_mod, n=500)
(z.ecdf.rate_mod / z.ecdf_diff.rate_mod) +
    plot_annotation(title='Zhang et al. (2021)')
ggsave('plots/ecdf_zhang_2021_rate_mod.png', width=8, height=8)

## Calculate log-likelihood
z.ll.rate_mod <- LL_discrete(z.fix_data.on_task$fixdur, z.fix_ucm.on_task.rate_mod$fixdur, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR, binwidth=z.BINWIDTH, delta=1/z.N_TRIALS) +
    LL_discrete(z.fix_data.mw$fixdur, z.fix_ucm.mw.rate_mod$fixdur, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR, binwidth=z.BINWIDTH, delta=1/z.N_TRIALS)
z.ll.rate_mod

## compare to N_mod model
z.likelihood.ratio <- -2*(z.ll.rate_mod - z.ll.N_mod)
z.likelihood.ratio
pchisq(z.likelihood.ratio, df=1, lower.tail=FALSE)


(wrap_plots(k.hist.rate_mod+ggtitle('Krasich et al. (2018)')+coord_cartesian(xlim=c(0, 1), ylim=c(0, .18)),
            z.hist.rate_mod+ggtitle('Zhang et al. (2021)')+coord_cartesian(xlim=c(0, 1), ylim=c(0, .18)),
            k.means.rate_mod+coord_cartesian(xlim=c(0, 1)),
            z.means.rate_mod+coord_cartesian(xlim=c(0, 1))) &
 theme(legend.position='bottom')) +
    plot_annotation(tag_levels=list(c('A', 'B', '', '')), title='Reduced Model') +
    plot_layout(heights=c(1, .33), guides='collect')
ggsave('plots/fit_ratemod.png', width=9, height=4.5)

((k.ecdf.rate_mod+ggtitle('Krasich et al. (2018)')+coord_cartesian(xlim=c(0, 1)) |
  z.ecdf.rate_mod+ggtitle('Zhang et al. (2021)')+coord_cartesian(xlim=c(0, 1))) +
 plot_layout(guides='collect')) /
    ((k.ecdf_diff.rate_mod+coord_cartesian(xlim=c(0, 1), ylim=c(-.075, .03)) |
      z.ecdf_diff.rate_mod+coord_cartesian(xlim=c(0, 1), ylim=c(-.075, .03))) +
     plot_layout(guides='collect')) +
    plot_annotation(tag_levels=list(c('A', 'B', '', '')), title='Reduced Model')
ggsave('plots/ecdf_ratemod.png', width=10, height=6)


##################################################################################################
##                           Optimze on-task & mw trials separately
##################################################################################################
z.bopt.on_task <- optimize('zhang_2021_on_task.rds', z.LL.on_task, z.bounds.separate)
z.params.on_task <- getBestPars(z.bopt.on_task)
z.params.on_task

z.bopt.mw <- optimize('zhang_2021_mw.rds', z.LL.mw, z.bounds.separate)
z.params.mw <- getBestPars(z.bopt.mw)
z.params.mw

z.fix_ucm.on_task.separate <- mclapply(1:10000, function(i) {
    UCM(N_timer=z.params.on_task$N_states, N_labile=z.params.on_task$N_states,
        N_nonlabile=z.params.on_task$N_states, N_motor=z.params.on_task$N_states,
        N_execution=z.params.on_task$N_states, t_timer=z.params.on_task$t_timer,
        t_labile=z.params.on_task$t_labile, t_nonlabile=z.params.on_task$t_nonlabile,
        t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION) %>%
        run(z.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= z.MIN_FIXDUR & fixdur <= z.MAX_FIXDUR)
z.fix_ucm.mw.separate <- mclapply(1:10000, function(i) {
    UCM(N_timer=z.params.mw$N_states, N_labile=z.params.mw$N_states,
        N_nonlabile=z.params.mw$N_states, N_motor=z.params.mw$N_states,
        N_execution=z.params.mw$N_states, t_timer=z.params.mw$t_timer,
        t_labile=z.params.mw$t_labile, t_nonlabile=z.params.mw$t_nonlabile,
        t_motor=z.T_MOTOR, t_execution=z.T_EXECUTION) %>%
        run(z.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= z.MIN_FIXDUR & fixdur <= z.MAX_FIXDUR)

## extract model fixations
z.fix_ucm.separate <- bind_rows(z.fix_ucm.on_task.separate, z.fix_ucm.mw.separate)
z.fix.separate <- bind_rows(z.fix_data %>% select(type, mw, fixdur),
                          z.fix_ucm.separate %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, z.MIN_FIXDUR, z.MAX_FIXDUR, z.BINWIDTH)))

## plot model fit (histograms/means)
z.hist.separate <- mw_hist(z.fix.separate)
z.means.separate <- mw_means_plot(z.fix.separate)
(z.hist.separate / z.means.separate & coord_cartesian(xlim=c(0, 1))) +
        plot_layout(heights=c(1, .25), guides='collect') +
        plot_annotation(title='Zhang et al. (2021)')
ggsave('plots/fit_zhang_2021_separate.png', width=8, height=4)

## plot model fit (ecdfs/ecdf_differences)
z.ecdf.separate <- mw_ecdf_plot(z.fix.separate)
z.ecdf_diff.separate <- mw_ecdf_difference_plot(z.fix.separate, n=500)
(z.ecdf.separate / z.ecdf_diff.separate) +
    plot_annotation(title='Zhang et al. (2021)')
ggsave('plots/ecdf_zhang_2021_separate.png', width=8, height=8)

## Calculate log-likelihood
LL_discrete(z.fix_data.on_task$fixdur, z.fix_ucm.on_task.separate$fixdur, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR, binwidth=z.BINWIDTH, delta=1/z.N_TRIALS) +
    LL_discrete(z.fix_data.mw$fixdur, z.fix_ucm.mw.separate$fixdur, min=z.MIN_FIXDUR, max=z.MAX_FIXDUR, binwidth=z.BINWIDTH, delta=1/z.N_TRIALS)

