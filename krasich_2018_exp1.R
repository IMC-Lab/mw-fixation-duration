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


## fit LATER
fix_later <- fix_data %>% group_by(mw) %>%
    LATER.fit(100000) %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_later.on_task <- fix_later %>% filter(mw==0)
fix_later.mw <- fix_later %>% filter(mw==1)


##################################################################################################
##                 Use this code to optimze on-task & mw trials separately
##################################################################################################

bopt.on_task <- optimize('krasich_2018_on_task.rds', LL.on_task, bounds.separate)
params.on_task <- getBestPars(bopt.on_task)
params.on_task

bopt.mw <- optimize('krasich_2018_mw.rds', LL.mw, bounds.separate)
params.mw <- getBestPars(bopt.mw)
params.mw


fix_ucm.on_task <- mclapply(1:10000, function(i) {
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
fix_ucm.mw <- mclapply(1:10000, function(i) {
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

fix_ucm <- bind_rows(fix_ucm.on_task, fix_ucm.mw)


## plot the results
fix <- bind_rows(fix_data %>% select(type, mw, fixdur),
                 fix_ucm %>% select(type, mw, fixdur),
                 fix_later %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, MIN_FIXDUR, MAX_FIXDUR, BINWIDTH)))

fit_plot('Krasich et al. (2018) Model Fits (Separate)', fix, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
ggsave('plots/fit_krasich_2018_separate.png', width=8, height=4)

reciprobit_plot('Krasich et al. (2018) Model Fits (Separate)', fix, min=MIN_FIXDUR)
ggsave('plots/reciprobit_krasich_2018_separate.png', width=16, height=8)


## Calculate log-likelihood of UCM and LATER
LL_discrete(fix_data.on_task$fixdur, fix_ucm.on_task$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_ucm.mw$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)

LL_discrete(fix_data.on_task$fixdur, fix_later.on_task$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_later.mw$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)


##################################################################################################
##             Use this code to optimze on-task & mw trials with rate modulation
##################################################################################################

bopt.rate_mod <- optimize('krasich_2018_rate_mod.rds', LL.rate_mod, bounds.rate_mod)
params.rate_mod <- getBestPars(bopt.rate_mod)
params.rate_mod

fix_ucm.on_task <- mclapply(1:10000, function(i) {
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
fix_ucm.mw <- mclapply(1:10000, function(i) {
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

fix_ucm <- bind_rows(fix_ucm.on_task, fix_ucm.mw)


## plot the results
fix <- bind_rows(fix_data %>% select(type, mw, fixdur),
                 fix_ucm %>% select(type, mw, fixdur),
                 fix_later %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, MIN_FIXDUR, MAX_FIXDUR, BINWIDTH)))

fit_plot('Krasich et al. (2018) Model Fits (Rate modulation)', fix, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
ggsave('plots/fit_krasich_2018_rate_mod.png', width=8, height=4)

reciprobit_plot('Krasich et al. (2018) Model Fits (Rate modulation)', fix, min=MIN_FIXDUR)
ggsave('plots/reciprobit_krasich_2018_rate_mod.png', width=16, height=8)


## Calculate log-likelihood of UCM and LATER
LL_discrete(fix_data.on_task$fixdur, fix_ucm.on_task$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_ucm.mw$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)

LL_discrete(fix_data.on_task$fixdur, fix_later.on_task$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_later.mw$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)


##################################################################################################
##        Use this code to optimze on-task & mw trials with rate modulation & N modulation
##################################################################################################

bopt.N_mod <- optimize('krasich_2018_N_mod.rds', LL.jointN_mod, bounds.N_mod)
params.N_mod <- getBestPars(bopt.N_mod)
params.N_mod

## Simulate lots of fixations from the best parameter values
ucm.on_task <- mclapply(1:10000, function(i) {
    UCM(N_timer=params.N_mod$N_states, N_labile=params.N_mod$N_states,
        N_nonlabile=params.N_mod$N_states, N_motor=params.N_mod$N_states,
        N_execution=params.N_mod$N_states, t_timer=params.N_mod$t_timer,
        t_labile=params.N_mod$t_labile, t_nonlabile=params.N_mod$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION) %>%
        run(TRIAL_DUR) %>%
        wrap()
})
ucm.mw <- mclapply(1:10000, function(i) {
    UCM(N_timer=round(params.N_mod$N_states * params.N_mod$N_mod),
        N_labile=round(params.N_mod$N_states * params.N_mod$N_mod),
        N_nonlabile=round(params.N_mod$N_states * params.N_mod$N_mod),
        N_motor=round(params.N_mod$N_states * params.N_mod$N_mod),
        N_execution=round(params.N_mod$N_states * params.N_mod$N_mod),
        t_timer=params.N_mod$t_timer,
        t_labile=params.N_mod$t_labile, t_nonlabile=params.N_mod$t_nonlabile,
        t_motor=T_MOTOR, t_execution=T_EXECUTION, modulation=params.N_mod$modulation) %>%
        run(TRIAL_DUR) %>%
        wrap()
})

fix_ucm.on_task <- get_fixations(ucm.on_task) %>%
    left_join(get_cancellations(ucm.on_task)) %>%
    mutate(mw=0, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)
fix_ucm.mw <- get_fixations(ucm.mw) %>%
    left_join(get_cancellations(ucm.mw)) %>%
    mutate(mw=1, fixdur=duration, type='UCM')  %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR)

fix_ucm <- bind_rows(fix_ucm.on_task, fix_ucm.mw)


## plot the results
fix <- bind_rows(fix_data %>% select(type, mw, fixdur),
                 fix_ucm %>% select(type, mw, fixdur),
                 fix_later %>% select(type, mw, fixdur)) %>%
    nest(data=c(fixdur)) %>%
    mutate(hist=map(data, ~ discretize(.x$fixdur, MIN_FIXDUR, MAX_FIXDUR, BINWIDTH)))

fit_plot('Krasich et al. (2018) Model Fits (N modulation)', fix, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH)
ggsave('plots/fit_krasich_2018_N_mod.png', width=8, height=4)

reciprobit_plot('Krasich et al. (2018) Model Fits (N modulation)', fix, min=MIN_FIXDUR)
ggsave('plots/reciprobit_krasich_2018_N_mod.png', width=16, height=8)


## Calculate log-likelihood of UCM and LATER
LL_discrete(fix_data.on_task$fixdur, fix_ucm.on_task$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_ucm.mw$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)

LL_discrete(fix_data.on_task$fixdur, fix_later.on_task$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS) +
    LL_discrete(fix_data.mw$fixdur, fix_later.mw$fixdur, min=MIN_FIXDUR, max=MAX_FIXDUR, binwidth=BINWIDTH, delta=1/N_TRIALS)

summary(glm(n_cancellations ~ factor(mw), fix_ucm, family='poisson'))

ggplot(fix_ucm, aes(x=n_cancellations, fill=as.factor(mw))) +
    geom_histogram(aes(y=after_stat(count / sum(count))), alpha=.5, binwidth=1, data=fix_ucm.on_task) +
    geom_histogram(aes(y=after_stat(count / sum(count))), alpha=.5, binwidth=1, data=fix_ucm.mw) +
    geom_vline(xintercept=mean(fix_ucm.on_task$n_cancellations), color= pal_jco()(2)[1], size=1) +
    geom_vline(xintercept=mean(fix_ucm.mw$n_cancellations), color= pal_jco()(2)[2], size=1) +
    scale_fill_jco(name='Probe Response', labels=c('On-task', 'Mind-wandering')) +
    xlab('Number of Cancellations') + ylab('Proportion') +
    theme_bw()
ggsave('plots/cancellations_krasich.png', width=8, height=6)

ggplot(fix_ucm, aes(x=duration, group=n_cancellations, fill=factor(n_cancellations))) +
    scale_fill_viridis(name='Cancellations', discrete=TRUE) +
    geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=0.06, data=fix_ucm.on_task) +
    geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=0.06, data=fix_ucm.mw) +
    coord_cartesian(xlim=c(0, MAX_FIXDUR)) +
    facet_grid( ~ mw, labeller=labeller(mw=c('0'='On-task', '1'='Mind-wandering'))) +
    xlab('Fixation Duration (s)') + ylab('Proportion') +
    theme_bw()

ggsave('plots/cancellations_duration_krasich.png', width=8, height=6)

discretize(fix_ucm.on_task$n_cancellations, min=-0.5, max=9.5, binwidth=1) %>%
    mutate(mw=0) %>%
    bind_rows(discretize(fix_ucm.mw$n_cancellations, min=-0.5, max=9.5, binwidth=1) %>%
              mutate(mw=1)) %>%
    ggplot(aes(x=mid, y=p, fill=factor(mw))) +
    geom_col(position='identity', alpha=.75) +
    theme_bw()

    
##pivot_wider(names_from='mw', values_from=c('count', 'p'), names_sep='.')
