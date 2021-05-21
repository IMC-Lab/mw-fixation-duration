library(httr)
library(simmer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggsci)
library(parallel)
library(ParBayesianOptimization)

## Load UCM source code
source('UCM.R')
options(mc.cores=parallel::detectCores())
set.seed(2021)

## Set optimization parameters
OPT_FILE <- 'zhang_2021.rds'
N_TRIALS <- 2500    ## number of UCM trials to simulate per epoch
MIN_FIXDUR <- 0.08  ## only keep fixations over this limit
MAX_FIXDUR <- 2    ## only keep fixations under this limit
BINWIDTH <- 0.05    ## the width of histogram bins to calculate LL
INIT_POINTS <- 16   ## number of points to initialize optimization
KAPPA <- 5.0        ## importance of uncertainty in optimization
bounds <- list(N_states=c(10L, 30L),
               t_timer=c(.1, .3),
               t_labile=c(.1, .25),
               t_nonlabile=c(.025, .1),
               MW_mod=c(0.25, 1))


## Load data
GET('https://osf.io/fpqsj/download', write_disk(tf <- tempfile(fileext = ".csv")))
fix_data <- read.csv(tf, header=TRUE) %>%
    tibble %>%
    mutate(fixdur=duration/1000,
           proberesp=factor(Attention, levels=c('On-task', 'Mind-wandering')),
           mw=ifelse(proberesp == 'Mind-wandering', 1, 0)) %>%
    ## filter as per the manuscript
    drop_na() %>%
    filter(fixdur >= MIN_FIXDUR & fixdur <= MAX_FIXDUR) %>%
    filter(x_pos > 0 & x_pos < 1024) %>%
    filter(y_pos > 0 & y_pos < 768)


## Plot data
p.1 <- ggplot(fix_data) +
    aes(x=fixdur, fill=factor(mw)) +
    geom_histogram(aes(y=after_stat(count / sum(count))), binwidth=.01, show.legend=FALSE) +
    ylab('Proportion') +
    scale_fill_jco() +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

p.2 <- ggplot(fix_data) +
    aes(x=fixdur, y=Attention, fill=factor(mw)) +
    stat_summary(fun.data=mean_cl_normal, geom='col') +
    stat_summary(fun.data=mean_cl_normal, geom='errorbar') +
    scale_fill_jco(name='Probe Response', labels=c('On-task', 'Mind-wandering')) +
    xlab('Fixation Duration (ms)') +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y=element_blank())

(p.1 / p.2 & xlim(0, 1)) +
    plot_layout(heights=c(1, .25), guides='collect') +
    plot_annotation(title='Zhang et al., (2021)')
ggsave('plots/hist_zhang_2021.png', width=8, height=6)

#' Get the log likelihood (LL) of the data given some parameterization of UCM
#'
#' @param N_states the number of state transitions for *all* stages
#' @param t_timer, t_labile, t_nonlabile the mean duration of the timer, labile, and nonlabile stages
#' @param MW_mod the factor by which rates are modulated during mind-wandering
#' @return a list containing:
#'  Score- the LL of the data
#'  LL.on_task- the LL of on-task trials
#'  LL.mw- the LL of mind-wandering trials
#'  N_TRIALS- the number of simulated UCM trials per condition
#'  min- the minimum fixation duration
#'  max- the maximum fixation duration
zhang_LL <- function(N_states, t_timer, t_labile, t_nonlabile, MW_mod) {
    LL.on_task <- 1:N_TRIALS %>%
        mclapply(function (i) {
            UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                N_motor=N_states, N_execution=N_states,
                t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile) %>%
                run(10) %>%
                wrap()
        }) %>%
        get_fixations() %>%
        pull(duration) %>%
        LL_discrete(fix_data$fixdur[fix_data$mw==0], ., min=MIN_FIXDUR, max=MAX_FIXDUR)

    LL.mw <- 1:N_TRIALS %>%
        mclapply(function (i) {
            UCM(N_timer=N_states, N_labile=N_states, N_nonlabile=N_states,
                N_motor=N_states, N_execution=N_states,
                t_timer=t_timer, t_labile=t_labile, t_nonlabile=t_nonlabile,
                modulation=MW_mod) %>%
                run(10) %>%
                wrap()
        }) %>%
        get_fixations() %>%
        pull(duration) %>%
        LL_discrete(fix_data$fixdur[fix_data$mw==1], ., min=MIN_FIXDUR, max=MAX_FIXDUR)

    list(Score=LL.on_task + LL.mw, LL.on_task=LL.on_task, LL.mw=LL.mw)
}

## Find the MLE of the UCM
bopt <- bayesOpt(zhang_LL, bounds=bounds, saveFile=OPT_FILE, initPoints=INIT_POINTS,
                 iters.n=5, kappa=KAPPA, covtype='gauss', plotProgress=TRUE)

while (bopt$scoreSummary$gpUtility[nrow(bopt$scoreSummary)] >= 0.1) {
    bopt <- addIterations(bopt, iters.n=1, plotProgress=TRUE)
}


## Run the UCM with the best-fitting parameters
ucm.MLE <- getBestPars(bopt)
ucm.on_task <- mclapply(1:5000, function(i) {
    UCM(N_timer=ucm.MLE$N_states, N_labile=ucm.MLE$N_states,
        N_nonlabile=ucm.MLE$N_states, N_motor=ucm.MLE$N_states,
        N_execution=ucm.MLE$N_states, t_timer=ucm.MLE$t_timer,
        t_labile=ucm.MLE$t_labile, t_nonlabile=ucm.MLE$t_nonlabile) %>%
        run(10) %>%
        wrap()
})
ucm.mw <- mclapply(1:5000, function(i) {
    UCM(N_timer=ucm.MLE$N_states, N_labile=ucm.MLE$N_states,
        N_nonlabile=ucm.MLE$N_states, N_motor=ucm.MLE$N_states,
        N_execution=ucm.MLE$N_states, t_timer=ucm.MLE$t_timer,
        t_labile=ucm.MLE$t_labile, t_nonlabile=ucm.MLE$t_nonlabile,
        modulation=ucm.MLE$MW_mod) %>%
        run(10) %>%
        wrap()
})

## Get simulated fixations
fix_ucm.on_task <- get_fixations(ucm.on_task)
fix_ucm.mw <- get_fixations(ucm.mw)


## Plot the simulated fixations against the raw data
rbind(data.frame(type='Data', mw=fix_data$mw, fixdur=fix_data$fixdur),
      data.frame(type='UCM', mw=0, fixdur=fix_ucm.on_task$duration),
      data.frame(type='UCM', mw=1, fixdur=fix_ucm.mw$duration)) %>%
    filter(fixdur > .05) %>%
    ggplot(aes(x=fixdur)) +
    geom_histogram(aes(y=stat(density)), binwidth=0.01) +
    facet_grid(type ~ .) +
    scale_x_continuous(limits=c(0, 1.2)) +
    xlab('Fixation Duration (s)') + ylab('Proportion') +
    theme_bw()
