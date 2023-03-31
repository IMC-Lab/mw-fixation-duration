library(readxl)
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
k.T_MOTOR <- .03            ## duration of motor planning stage (fixed parameter)
k.T_EXECUTION <- .04        ## duration of saccade execution (fixed parameter)
k.N_TRIALS <- 2500          ## number of UCM trials to simulate per epoch to estimate LL
k.MIN_FIXDUR <- 0.05        ## only keep fixations over this limit
k.MAX_FIXDUR <- 2           ## only keep fixations under this limit
k.BINWIDTH <- 0.05          ## the width of histogram bins to calculate LL
k.TRIAL_DUR <- 15           ## trial length (in seconds)
k.bounds.separate <- list(N_states=c(2L, 30L),
                          t_timer=c(.15, .375),
                          t_labile=c(.1, .225),
                          t_nonlabile=c(.025, .08))
k.bounds.rate_mod <- c(k.bounds.separate, list(modulation=c(0.25, 1)))
k.bounds.N_mod <- c(k.bounds.rate_mod, list(N_mod=c(0.25, 1)))

## Load data
k.fix_data <- read_excel('data/MWScenesTriad_datat.xlsx') %>%
    mutate(fixdur=fixdur/1000,
           Image_Duration=Image_Duration/1000,
           ID=factor(ID),
           proberesp=factor(proberesp, levels=c('notmw', 'mw')),
           mw=ifelse(proberesp == 'mw', 1, 0),
           type='Data') %>%
    filter(proberesp != 'noprobe' & `5sBins` <= 15 &
           fixdur >= k.MIN_FIXDUR & fixdur <= k.MAX_FIXDUR)

## count number of trials with/without MW
k.fix_data %>% group_by(ID, Image_Duration, proberesp) %>%
    summarize %>% group_by(proberesp) %>% count

k.fix_data.on_task <- k.fix_data %>% filter(mw==0)
k.fix_data.mw <- k.fix_data %>% filter(mw==1)

## Plot data
mw_plot('Krasich et al. (2018) Raw Data', k.fix_data, k.BINWIDTH)
ggsave('plots/hist_krasich_2018_exp1.png', width=8, height=6)


## Define log likelihood functions
k.LL.on_task <- LL(k.fix_data.on_task$fixdur, k.bounds.separate,
                   n_trials=k.N_TRIALS, trial_dur=k.TRIAL_DUR, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR,
                   binwidth=k.BINWIDTH, default_params=list(t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION))
k.LL.mw <- LL(k.fix_data.mw$fixdur, k.bounds.separate,
              n_trials=k.N_TRIALS, trial_dur=k.TRIAL_DUR, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR,
              binwidth=k.BINWIDTH, default_params=list(t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION))
k.LL.rate_mod <- LL.joint(k.fix_data.on_task$fixdur, k.fix_data.mw$fixdur, k.bounds.rate_mod,
                          n_trials=k.N_TRIALS, trial_dur=k.TRIAL_DUR, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR,
                          binwidth=k.BINWIDTH, default_params=list(t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION))
k.LL.N_mod <- LL.joint(k.fix_data.on_task$fixdur, k.fix_data.mw$fixdur, k.bounds.N_mod,
                       n_trials=k.N_TRIALS, trial_dur=k.TRIAL_DUR, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR,
                       binwidth=k.BINWIDTH, default_params=list(t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION))



## Test for differences between MW/Attentive Viewing using permuted KS test
k.KS <- KS.permutation(k.fix_data, n=10000)
mean(k.KS$D.null >= k.KS$D)
ggplot(k.KS, aes(x=D.null)) +
    geom_histogram(bins=100) +
    geom_vline(xintercept=k.KS$D) +
    theme_classic()


##################################################################################################
##                      Optimze using rate modulation & N modulation
##################################################################################################
k.bopt.N_mod <- optimize('krasich_2018_N_mod.rds', k.LL.N_mod, k.bounds.N_mod)
k.params.N_mod <- getBestPars(k.bopt.N_mod)
k.params.N_mod

## Simulate lots of fixations from the best parameter values
k.ucm.on_task.N_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=k.params.N_mod$N_states, N_labile=k.params.N_mod$N_states,
        N_nonlabile=k.params.N_mod$N_states, N_motor=k.params.N_mod$N_states,
        N_execution=k.params.N_mod$N_states, t_timer=k.params.N_mod$t_timer,
        t_labile=k.params.N_mod$t_labile, t_nonlabile=k.params.N_mod$t_nonlabile,
        t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION) %>%
        run(k.TRIAL_DUR) %>%
        wrap()
})
k.fix_ucm.on_task.N_mod <- k.ucm.on_task.N_mod %>%
    get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= k.MIN_FIXDUR & fixdur <= k.MAX_FIXDUR)
k.ucm.mw.N_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=round(k.params.N_mod$N_states * k.params.N_mod$N_mod),
        N_labile=round(k.params.N_mod$N_states * k.params.N_mod$N_mod),
        N_nonlabile=round(k.params.N_mod$N_states * k.params.N_mod$N_mod),
        N_motor=k.params.N_mod$N_states, N_execution=k.params.N_mod$N_states,
        t_timer=k.params.N_mod$t_timer,
        t_labile=k.params.N_mod$t_labile, t_nonlabile=k.params.N_mod$t_nonlabile,
        t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION, modulation=k.params.N_mod$modulation) %>%
        run(k.TRIAL_DUR) %>%
        wrap()
})
k.fix_ucm.mw.N_mod <- k.ucm.mw.N_mod %>%
    get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= k.MIN_FIXDUR & fixdur <= k.MAX_FIXDUR)

## extract the model fixations
k.fix_ucm.N_mod <- bind_rows(k.fix_ucm.on_task.N_mod, k.fix_ucm.mw.N_mod)
k.fix.N_mod <- bind_rows(k.fix_data %>% select(type, mw, fixdur),
                         k.fix_ucm.N_mod %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, k.MIN_FIXDUR, k.MAX_FIXDUR, k.BINWIDTH)))

## calculate descriptive statistics
k.fix.N_mod %>%
    arrange(type, mw) %>%
    mutate(mean=map_dbl(data, ~mean(.x$fixdur*1000)),
           sd=map_dbl(data, ~sd(.x$fixdur*1000)),
           lower=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymin),
           upper=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymax)) %>%
    select(-data, -hist) %>%
    as.data.frame

## plot model fit (histograms/means)
k.hist.N_mod <- mw_hist(k.fix.N_mod)
k.means.N_mod <- mw_means_plot(k.fix.N_mod)
(k.hist.N_mod / k.means.N_mod & coord_cartesian(xlim=c(0, 1))) +
        plot_layout(heights=c(1, .25), guides='collect') +
        plot_annotation(title='Krasich et al. (2018)')
ggsave('plots/fit_krasich_2018_N_mod.png', width=8, height=4)

## plot model fit (ecdfs/ecdf_differences)
k.ecdf.N_mod <- mw_ecdf_plot(k.fix.N_mod)
k.ecdf_diff.N_mod <- mw_ecdf_difference_plot(k.fix.N_mod, n=500)
(k.ecdf.N_mod / k.ecdf_diff.N_mod) +
    plot_annotation(title='Krasich et al. (2018)')
ggsave('plots/ecdf_krasich_2018_N_mod.png', width=8, height=8)


## Calculate log-likelihood
k.ll.N_mod <- LL_discrete(k.fix_data.on_task$fixdur, k.fix_ucm.on_task.N_mod$fixdur, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR, binwidth=k.BINWIDTH, delta=1/k.N_TRIALS) +
    LL_discrete(k.fix_data.mw$fixdur, k.fix_ucm.mw.N_mod$fixdur, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR, binwidth=k.BINWIDTH, delta=1/k.N_TRIALS)
k.ll.N_mod


## calculate cancellation rates
k.ucm.on_task.N_mod %>%
    get_cancellations() %>%
    summarize(M=mean(cancelled))

cancellation_prob(k.params.N_mod$N_states, k.params.N_mod$N_states, ## derived value
                  k.params.N_mod$t_timer, k.params.N_mod$t_labile)

k.ucm.mw.N_mod %>%
    get_cancellations() %>%
    summarize(M=mean(cancelled))

cancellation_prob(round(k.params.N_mod$N_states * k.params.N_mod$N_mod),  ## derived value
                  round(k.params.N_mod$N_states * k.params.N_mod$N_mod),
                  k.params.N_mod$t_timer, k.params.N_mod$t_labile)

##################################################################################################
##                               Optimze w/ rate modulation only
##################################################################################################
k.bopt.rate_mod <- optimize('krasich_2018_rate_mod.rds', k.LL.rate_mod, k.bounds.rate_mod)
k.params.rate_mod <- getBestPars(k.bopt.rate_mod)
k.params.rate_mod

k.fix_ucm.on_task.rate_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=k.params.rate_mod$N_states, N_labile=k.params.rate_mod$N_states,
        N_nonlabile=k.params.rate_mod$N_states, N_motor=k.params.rate_mod$N_states,
        N_execution=k.params.rate_mod$N_states, t_timer=k.params.rate_mod$t_timer,
        t_labile=k.params.rate_mod$t_labile, t_nonlabile=k.params.rate_mod$t_nonlabile,
        t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION) %>%
        run(k.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= k.MIN_FIXDUR & fixdur <= k.MAX_FIXDUR)
k.fix_ucm.mw.rate_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=k.params.rate_mod$N_states, N_labile=k.params.rate_mod$N_states,
        N_nonlabile=k.params.rate_mod$N_states, N_motor=k.params.rate_mod$N_states,
        N_execution=k.params.rate_mod$N_states, t_timer=k.params.rate_mod$t_timer,
        t_labile=k.params.rate_mod$t_labile, t_nonlabile=k.params.rate_mod$t_nonlabile,
        t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION, modulation=k.params.rate_mod$modulation) %>%
        run(k.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= k.MIN_FIXDUR & fixdur <= k.MAX_FIXDUR)

## extract the model fixations
k.fix_ucm.rate_mod <- bind_rows(k.fix_ucm.on_task.rate_mod, k.fix_ucm.mw.rate_mod)
k.fix.rate_mod <- bind_rows(k.fix_data %>% select(type, mw, fixdur),
                            k.fix_ucm.rate_mod %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, k.MIN_FIXDUR, k.MAX_FIXDUR, k.BINWIDTH)))

## calculate descriptive statistics
k.fix.rate_mod %>%
    mutate(mean=map_dbl(data, ~mean(.x$fixdur*1000)),
           sd=map_dbl(data, ~sd(.x$fixdur*1000)),
           lower=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymin),
           upper=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymax)) %>%
    select(-data, -hist) %>%
    as.data.frame

## plot model fit (histograms/means)
k.hist.rate_mod <- mw_hist(k.fix.rate_mod)
k.means.rate_mod <- mw_means_plot(k.fix.rate_mod)
(k.hist.rate_mod / k.means.rate_mod & coord_cartesian(xlim=c(0, 1))) +
        plot_layout(heights=c(1, .25), guides='collect') +
        plot_annotation(title='Krasich et al. (2018)')
ggsave('plots/fit_krasich_2018_rate_mod.png', width=8, height=4)

## plot model fit (ecdfs/ecdf_differences)
k.ecdf.rate_mod <- mw_ecdf_plot(k.fix.rate_mod)
k.ecdf_diff.rate_mod <- mw_ecdf_difference_plot(k.fix.rate_mod, n=500)
(k.ecdf.rate_mod / k.ecdf_diff.rate_mod) +
    plot_annotation(title='Krasich et al. (2018)')
ggsave('plots/ecdf_krasich_2018_rate_mod.png', width=8, height=8)



## Calculate log-likelihood
k.ll.rate_mod <- LL_discrete(k.fix_data.on_task$fixdur, k.fix_ucm.on_task.rate_mod$fixdur, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR, binwidth=k.BINWIDTH, delta=1/k.N_TRIALS) +
    LL_discrete(k.fix_data.mw$fixdur, k.fix_ucm.mw.rate_mod$fixdur, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR, binwidth=k.BINWIDTH, delta=1/k.N_TRIALS)
k.ll.rate_mod

## compare to N_mod model
k.likelihood.ratio <- -2*(k.ll.rate_mod - k.ll.N_mod)
k.likelihood.ratio
pchisq(k.likelihood.ratio, df=1, lower.tail=FALSE)



##################################################################################################
##                           Optimze on-task & mw trials separately
##################################################################################################
k.bopt.on_task <- optimize('krasich_2018_on_task.rds', k.LL.on_task, k.bounds.N.separate)
k.params.on_task <- getBestPars(k.bopt.on_task)
k.params.on_task

k.bopt.mw <- optimize('krasich_2018_mw.rds', k.LL.mw, k.bounds.separate)
k.params.mw <- getBestPars(k.bopt.mw)
k.params.mw

k.fix_ucm.on_task.separate <- mclapply(1:10000, function(i) {
    UCM(N_timer=k.params.on_task$N_states, N_labile=k.params.on_task$N_states,
        N_nonlabile=k.params.on_task$N_states, N_motor=k.params.on_task$N_states,
        N_execution=k.params.on_task$N_states, t_timer=k.params.on_task$t_timer,
        t_labile=k.params.on_task$t_labile, t_nonlabile=k.params.on_task$t_nonlabile,
        t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION) %>%
        run(k.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= k.MIN_FIXDUR & fixdur <= k.MAX_FIXDUR)
k.fix_ucm.mw.separate <- mclapply(1:10000, function(i) {
    UCM(N_timer=k.params.mw$N_states, N_labile=k.params.mw$N_states,
        N_nonlabile=k.params.mw$N_states, N_motor=k.params.mw$N_states,
        N_execution=k.params.mw$N_states, t_timer=k.params.mw$t_timer,
        t_labile=k.params.mw$t_labile, t_nonlabile=k.params.mw$t_nonlabile,
        t_motor=k.T_MOTOR, t_execution=k.T_EXECUTION) %>%
        run(k.TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= k.MIN_FIXDUR & fixdur <= k.MAX_FIXDUR)

## extract the model fixations
k.fix_ucm.separate <- bind_rows(k.fix_ucm.on_task.separate, k.fix_ucm.mw.separate)
k.fix.separate <- bind_rows(k.fix_data %>% select(type, mw, fixdur),
                            k.fix_ucm.separate %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, k.MIN_FIXDUR, k.MAX_FIXDUR, k.BINWIDTH)))

## calculate descriptive statistics
k.fix.separate %>%
    mutate(mean=map_dbl(data, ~mean(.x$fixdur*1000)),
           sd=map_dbl(data, ~sd(.x$fixdur*1000)),
           lower=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymin),
           upper=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymax)) %>%
    select(-data, -hist) %>%
    as.data.frame

## plot model fit (histograms/means)
k.hist.separate <- mw_hist(k.fix.separate)
k.means.separate <- mw_means_plot(k.fix.separate)
(k.hist.separate / k.means.separate & coord_cartesian(xlim=c(0, 1))) +
        plot_layout(heights=c(1, .25), guides='collect') +
        plot_annotation(title='Krasich et al. (2018)')
ggsave('plots/fit_krasich_2018_separate.png', width=8, height=4)

## plot model fit (ecdfs/ecdf_differences)
k.ecdf.separate <- mw_ecdf_plot(k.fix.separate)
k.ecdf_diff.separate <- mw_ecdf_difference_plot(k.fix.separate, n=500)
(k.ecdf.separate / k.ecdf_diff.separate) +
    plot_annotation(title='Krasich et al. (2018)')
ggsave('plots/ecdf_krasich_2018_separate.png', width=8, height=8)


## Calculate log-likelihood
k.ll.separate <- LL_discrete(k.fix_data.on_task$fixdur, k.fix_ucm.on_task.separate$fixdur, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR, binwidth=k.BINWIDTH, delta=1/k.N_TRIALS) +
    LL_discrete(k.fix_data.mw$fixdur, k.fix_ucm.mw.separate$fixdur, min=k.MIN_FIXDUR, max=k.MAX_FIXDUR, binwidth=k.BINWIDTH, delta=1/k.N_TRIALS)
k.ll.separate
