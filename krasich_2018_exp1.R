library(readxl)
library(httr)
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
T_MOTOR <- .03           ## duration of motor planning stage (fixed parameter)
T_EXECUTION <- .04       ## duration of saccade execution (fixed parameter)
N_TRIALS <- 2500         ## number of UCM trials to simulate per epoch to estimate LL
MIN_FIXDUR <- 0.05       ## only keep fixations over this limit
MAX_FIXDUR <- 2          ## only keep fixations under this limit
BINWIDTH <- 0.05         ## the width of histogram bins to calculate LL
TRIAL_DUR <- 15          ## trial length (in seconds)
INIT_POINTS <- 16        ## number of points to initialize optimization
KAPPA <- 5.0             ## importance of uncertainty in optimization
UTILITY_THRESH <- 0.001   ## stop optimization when utility is smaller than this value
bounds.separate <- list(N_states=c(2L, 30L),
                        t_timer=c(.1, .4),
                        t_labile=c(.1, .4),
                        t_nonlabile=c(.025, .25))
bounds.rate_mod <- c(bounds.separate, list(modulation=c(0.25, 1)))
bounds.N_mod <- c(bounds.rate_mod, list(N_mod=c(0.25, 1)))

## Load data
fix_data <- read_excel('data/MWScenesTriad_datat.xlsx') %>%
    mutate(fixdur=fixdur/1000,
           Image_Duration=Image_Duration/1000,
           ID=factor(ID),
           proberesp=factor(proberesp, levels=c('notmw', 'mw')),
           mw=ifelse(proberesp == 'mw', 1, 0),
           type='Data') %>%
    filter(proberesp != 'noprobe' & `5sBins` <= 15 &
           fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)

## count number of trials with/without MW
fix_data %>% group_by(ID, Image_Duration, proberesp) %>%
    summarize %>% group_by(proberesp) %>% count

fix_data.on_task <- fix_data[fix_data$mw == 0,]
fix_data.mw <- fix_data[fix_data$mw == 1,]

## Plot data
mw_plot('Krasich et al. (2018) Raw Data', fix_data, BINWIDTH)
ggsave('plots/hist_krasich_2018_exp1.png', width=8, height=6)


## Define log likelihood functions
LL.on_task <- LL(fix_data.on_task$fixdur, bounds.separate,
                 n_trials=N_TRIALS, trial_dur=TRIAL_DUR, min=MIN_FIXDUR, max=MAX_FIXDUR,
                 binwidth=BINWIDTH, default_params=list(t_motor=T_MOTOR, t_execution=T_EXECUTION))
LL.mw <- LL(fix_data.mw$fixdur, bounds.separate,
            n_trials=N_TRIALS, trial_dur=TRIAL_DUR, min=MIN_FIXDUR, max=MAX_FIXDUR,
            binwidth=BINWIDTH, default_params=list(t_motor=T_MOTOR, t_execution=T_EXECUTION))
LL.rate_mod <- LL.joint(fix_data.on_task$fixdur, fix_data.mw$fixdur, bounds.rate_mod,
                        n_trials=N_TRIALS, trial_dur=TRIAL_DUR, min=MIN_FIXDUR, max=MAX_FIXDUR,
                        binwidth=BINWIDTH, default_params=list(t_motor=T_MOTOR, t_execution=T_EXECUTION))
LL.N_mod <- LL.joint(fix_data.on_task$fixdur, fix_data.mw$fixdur, bounds.N_mod,
                     n_trials=N_TRIALS, trial_dur=TRIAL_DUR, min=MIN_FIXDUR, max=MAX_FIXDUR,
                     binwidth=BINWIDTH, default_params=list(t_motor=T_MOTOR, t_execution=T_EXECUTION))


##################################################################################################
##                      Optimze using rate modulation & N modulation
##################################################################################################
bopt.N_mod <- optimize('krasich_2018_N_mod.rds', LL.N_mod, bounds.N_mod)
params.N_mod <- getBestPars(bopt.N_mod)
params.N_mod

## Simulate lots of fixations from the best parameter values
ucm.on_task.N_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=params.N_mod$N_states, N_labile=params.N_mod$N_states,
        N_nonlabile=params.N_mod$N_states, N_motor=params.N_mod$N_states,
        N_execution=params.N_mod$N_states, t_timer=params.N_mod$t_timer,
        t_labile=params.N_mod$t_labile, t_nonlabile=params.N_mod$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
        run(TRIAL_DUR) %>%
        wrap()
})
ucm.mw.N_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=round(params.N_mod$N_states * params.N_mod$N_mod),
        N_labile=round(params.N_mod$N_states * params.N_mod$N_mod),
        N_nonlabile=round(params.N_mod$N_states * params.N_mod$N_mod),
        N_motor=params.N_mod$N_states, N_execution=params.N_mod$N_states,
        t_timer=params.N_mod$t_timer,
        t_labile=params.N_mod$t_labile, t_nonlabile=params.N_mod$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION, modulation=params.N_mod$modulation) %>%
        run(TRIAL_DUR) %>%
        wrap()
})

## Extract the model fixations
fix_ucm.on_task.N_mod <- get_fixations(ucm.on_task.N_mod) %>%
    left_join(get_cancellations(ucm.on_task.N_mod)) %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_ucm.mw.N_mod <- get_fixations(ucm.mw.N_mod) %>%
    left_join(get_cancellations(ucm.mw.N_mod)) %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_ucm.N_mod <- bind_rows(fix_ucm.on_task.N_mod, fix_ucm.mw.N_mod)


## plot the results
fix.N_mod <- bind_rows(fix_data %>% select(type, mw, fixdur),
                       fix_ucm.N_mod %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, MIN_FIXDUR, MAX_FIXDUR, BINWIDTH)))

fix.N_mod %>%
    mutate(mean=map_dbl(data, ~mean(.x$fixdur*1000)),
           sd=map_dbl(data, ~sd(.x$fixdur*1000)),
           lower=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymin),
           upper=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymax)) %>%
    select(-data, -hist) %>%
    as.data.frame

fit_plot('Krasich et al. (2018)', fix.N_mod)
ggsave('plots/fit_krasich_2018_N_mod.png', width=8, height=4)

difference_plot('Krasich et al. (2018)', fix.N_mod, max=MAX_FIXDUR, binwidth=BINWIDTH)
ggsave('plots/difference_krasich_2018_N_mod.png', width=8, height=4)

## Calculate log-likelihood
ll.N_mod <- LL_discrete(fix_data.on_task$fixdur, fix_ucm.on_task.N_mod$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_ucm.mw.N_mod$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)
ll.N_mod

## Estimate number of cancellations
fix_ucm.N_mod %>% group_by(type, mw) %>% summarize(M=mean(n_cancellations), SD=sd(n_cancellations))

ggplot(fix_ucm.N_mod, aes(x=n_cancellations, fill=as.factor(mw))) +
    geom_histogram(aes(y=after_stat(count / sum(count))), alpha=.5, binwidth=1, data=fix_ucm.on_task.N_mod) +
    geom_histogram(aes(y=after_stat(count / sum(count))), alpha=.5, binwidth=1, data=fix_ucm.mw.N_mod) +
    geom_vline(xintercept=mean(fix_ucm.on_task.N_mod$n_cancellations), color= pal_jco()(2)[1], size=1) +
    geom_vline(xintercept=mean(fix_ucm.mw.N_mod$n_cancellations), color= pal_jco()(2)[2], size=1) +
    scale_fill_jco(name='Probe Response', labels=c('On-task', 'Mind-wandering')) +
    xlab('Number of Cancellations') + ylab('Proportion') +
    theme_bw()
ggsave('plots/cancellations_krasich.png', width=8, height=6)

ggplot(fix_ucm.N_mod, aes(x=duration, group=n_cancellations, fill=factor(n_cancellations))) +
    scale_fill_viridis(name='Cancellations', discrete=TRUE) +
    geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=0.06, data=fix_ucm.on_task.N_mod) +
    geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=0.06, data=fix_ucm.mw.N_mod) +
    coord_cartesian(xlim=c(0, MAX_FIXDUR)) +
    facet_grid( ~ mw, labeller=labeller(mw=c('0'='On-task', '1'='Mind-wandering'))) +
    xlab('Fixation Duration (s)') + ylab('Proportion') +
    theme_bw()
ggsave('plots/cancellations_duration_krasich.png', width=8, height=6)




##################################################################################################
##                               Optimze w/ rate modulation only
##################################################################################################
bopt.rate_mod <- optimize('krasich_2018_rate_mod.rds', LL.rate_mod, bounds.rate_mod)
params.rate_mod <- getBestPars(bopt.rate_mod)
params.rate_mod

fix_ucm.on_task.rate_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=params.rate_mod$N_states, N_labile=params.rate_mod$N_states,
        N_nonlabile=params.rate_mod$N_states, N_motor=params.rate_mod$N_states,
        N_execution=params.rate_mod$N_states, t_timer=params.rate_mod$t_timer,
        t_labile=params.rate_mod$t_labile, t_nonlabile=params.rate_mod$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
        run(TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_ucm.mw.rate_mod <- mclapply(1:10000, function(i) {
    UCM(N_timer=params.rate_mod$N_states, N_labile=params.rate_mod$N_states,
        N_nonlabile=params.rate_mod$N_states, N_motor=params.rate_mod$N_states,
        N_execution=params.rate_mod$N_states, t_timer=params.rate_mod$t_timer,
        t_labile=params.rate_mod$t_labile, t_nonlabile=params.rate_mod$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION, modulation=params.rate_mod$modulation) %>%
        run(TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_ucm.rate_mod <- bind_rows(fix_ucm.on_task.rate_mod, fix_ucm.mw.rate_mod)


## plot the results
fix.rate_mod <- bind_rows(fix_data %>% select(type, mw, fixdur),
                            fix_ucm.rate_mod %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, MIN_FIXDUR, MAX_FIXDUR, BINWIDTH)))

fix.rate_mod %>%
    mutate(mean=map_dbl(data, ~mean(.x$fixdur*1000)),
           sd=map_dbl(data, ~sd(.x$fixdur*1000)),
           lower=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymin),
           upper=map_dbl(data, ~mean_cl_normal(.x$fixdur*1000)$ymax)) %>%
    select(-data, -hist) %>%
    as.data.frame


fit_plot('Krasich et al. (2018)', fix.rate_mod)
ggsave('plots/fit_krasich_2018_rate_mod.png', width=8, height=4)

difference_plot('Krasich et al. (2018)', fix.rate_mod, max=MAX_FIXDUR, binwidth=BINWIDTH)
ggsave('plots/difference_krasich_2018_rate_mod.png', width=8, height=4)

## Calculate log-likelihood
ll.rate_mod <- LL_discrete(fix_data.on_task$fixdur, fix_ucm.on_task.rate_mod$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_ucm.mw.rate_mod$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)
ll.rate_mod

## compare to N_mod model
likelihood.ratio <- -2*(ll.rate_mod - ll.N_mod)
likelihood.ratio
pchisq(likelihood.ratio, df=1, lower.tail=FALSE)



##################################################################################################
##                           Optimze on-task & mw trials separately
##################################################################################################
bopt.on_task <- optimize('krasich_2018_on_task.rds', LL.on_task, bounds.N.separate)
params.on_task <- getBestPars(bopt.on_task)
params.on_task

bopt.mw <- optimize('krasich_2018_mw.rds', LL.mw, bounds.separate)
params.mw <- getBestPars(bopt.mw)
params.mw

fix_ucm.on_task.separate <- mclapply(1:10000, function(i) {
    UCM(N_timer=params.on_task$N_states, N_labile=params.on_task$N_states,
        N_nonlabile=params.on_task$N_states, N_motor=params.on_task$N_states,
        N_execution=params.on_task$N_states, t_timer=params.on_task$t_timer,
        t_labile=params.on_task$t_labile, t_nonlabile=params.on_task$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
        run(TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_ucm.mw.separate <- mclapply(1:10000, function(i) {
    UCM(N_timer=params.mw$N_states, N_labile=params.mw$N_states,
        N_nonlabile=params.mw$N_states, N_motor=params.mw$N_states,
        N_execution=params.mw$N_states, t_timer=params.mw$t_timer,
        t_labile=params.mw$t_labile, t_nonlabile=params.mw$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
        run(TRIAL_DUR) %>%
        wrap()
}) %>% get_fixations() %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_ucm.separate <- bind_rows(fix_ucm.on_task.separate, fix_ucm.mw.separate)


## plot the results
fix.separate <- bind_rows(fix_data %>% select(type, mw, fixdur),
                          fix_ucm.separate %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, MIN_FIXDUR, MAX_FIXDUR, BINWIDTH)))

fit_plot('Krasich et al. (2018)', fix.separate, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
ggsave('plots/fit_krasich_2018_separate.png', width=8, height=4)

difference_plot('Krasich et al. (2018)', fix.separate, max=MAX_FIXDUR, binwidth=BINWIDTH)
ggsave('plots/difference_krasich_2018_separate.png', width=8, height=4)

## Calculate log-likelihood
ll.separate <- LL_discrete(fix_data.on_task$fixdur, fix_ucm.on_task.separate$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_ucm.mw.separate$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)
ll.separate

## Calculate BIC
BIC.separate <- length(bounds.separate)*log(nrow(fix_data)) - 2*ll.separate
BIC.separate
