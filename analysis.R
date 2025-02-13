#------------------------------
# MAIN ANALYSIS
#------------------------------

# This script contains code for the main analysis.
# Two different types of models are fitted to the pre-processed data.

#-----------HEADER-------------
# load packages
library(knitr)        # output table editing etc.
library(tidyverse)    # data management and ggplot2
library(mgcv)         # Additive regression models
library(Hmisc)
library(emmeans)
theme_set(theme_bw()) # set global theme for ggplot2 plots
source("plot_predictions.R")

# define paths
path_prepData  <- "0_data/2_preparedData/"
path_resTables <- "2_resultFiles/tables/"
path_resImages <- "2_resultFiles/images/"

dat <- read.csv2(paste0(path_prepData, "data_prepared.csv"), fileEncoding = "utf-8")
#------------------------------

#------------------------------ DATA PREPARATION
# calculate mean nLetters.WD
nLetters.WD <- round(mean(dat$nLetters.WD), 2)
# ensure that Participant and Topic are encoded as factor variables
dat <- dat %>%
  mutate(Participant = factor(Participant),
         Topic       = factor(Topic))

# MODEL 1 SPECIFICATION

model_details <- T
# response variables
response_vec <- c("TRT","FRT","RRT")
# define for which models the 'letters per word' effect should be estimated
# nonlinearly (default: all)
nonlinear_effect <- c("TRT","FRT","RRT") # set to c() to estimate all 'letters per word' effects as linear effects
no_effect <- c() # set to c("TRT") to remove the 'letters per word' effect in TRT-model
# AOI.conditions per model
AOIcond_list <- list(
  "Model1" = c(
    "CE0_1", "C1_1", "C2_1", "C3_1", "C4_1", "C5_1", "C6_1", "C7_1",
             "E1_1", "E2_1", "E3_1", "E4_1", "E5_1", "E6_1", "E7_1",
    "CE0_2", "C1_2", "C2_2", "C3_2", "C4_2", "C5_2", "C6_2", "C7_2",
             "E1_2", "E2_2", "E3_2", "E4_2", "E5_2", "E6_2", "E7_2",
             "C1_3", "C2_3", "C3_3", "C4_3", "C5_3", "C6_3", "C7_3",
             "E1_3", "E2_3", "E3_3", "E4_3", "E5_3", "E6_3", "E7_3",
    "CE0_4", "C1_4", "C2_4", "C3_4", "C4_4", "C5_4", "C6_4", "C7_4",
             "E1_4", "E2_4", "E3_4", "E4_4", "E5_4", "E6_4", "E7_4",
    "CE0_5", "C1_5", "C2_5", "C3_5", "C4_5", "C5_5", "C6_5", "C7_5",
             "E1_5", "E2_5", "E3_5", "E4_5", "E5_5", "E6_5", "E7_5",
             "C1_6", "C2_6", "C3_6", "C4_6", "C5_6", "C6_6", "C7_6",
             "E1_6", "E2_6", "E3_6", "E4_6", "E5_6", "E6_6", "E7_6",
    "CE0_7", "C1_7", "C2_7", "C3_7", "C4_7", "C5_7", "C6_7", "C7_7",
             "E1_7", "E2_7", "E3_7", "E4_7", "E5_7", "E6_7", "E7_7",
             "C1_8", "C2_8", "C3_8", "C4_8", "C5_8", "C6_8", "C7_8",
             "E1_8", "E2_8", "E3_8", "E4_8", "E5_8", "E6_8", "E7_8"
  )
)

AOIcond_vec <- sort(unique(unlist(AOIcond_list)))
dat <- dat %>% filter(AOI.condition %in% AOIcond_vec) %>%
  mutate(AOI.condition = factor(AOI.condition, levels = AOIcond_vec))

#------------------------------ DESCRIPTIVES
# Observations per AOI.condition
freq <- sapply(AOIcond_vec, function(AOIcond) {
  length(unique(dat$observation[dat$AOI.condition == AOIcond]))
})
tab <- data.frame("AOI.condition" = AOIcond_vec,
                  "Freq" = unname(freq))
# kable(tab)
#------------------------------

#-----------MODEL 1------------- 

#------------------------------ MODEL 1 ESTIMATION
model_list <- lapply(seq_len(length(AOIcond_list)), function(i) {
  m_list <- lapply(seq_len(length(response_vec)), function(j) {
    message(paste0("Estimate model ",i,"_",response_vec[j],"..."))
    dat_i <- dat %>% mutate(AOI.condition = relevel(AOI.condition, ref = AOIcond_list[[i]][1]))
    fm_covars <- "~ AOI.condition"
    if (!(response_vec[j] %in% no_effect)) { # add 'letters per word' effect
      if (response_vec[j] %in% nonlinear_effect) {
        fm_covars <- paste(fm_covars, "+ s(nLetters.WD, bs='ps', k = 4)")
        } else { # RRT and TRT
          fm_covars <- paste(fm_covars, "+ nLetters.WD")
        }
      }
    fm_covars <- paste(fm_covars, "+ s(Topic, bs='re') + s(Participant, bs='re')")
    fm <- as.formula(paste0(response_vec[j], ".WD", fm_covars))
    bam(fm, data = dat_i, method = "REML", nthreads = 10)
    })
  names(m_list) <- response_vec
  m_list
  })

names(model_list) <- paste0("Model", seq_len(length(model_list)))
#------------------------------

#------------------------------ MODEL 1 RESULTS
for (i in seq_len(length(AOIcond_list))) {
  # retrieve the y limits for plotting over all response variables
  max_fitted <- sapply(model_list[[i]], function(model) { max(model$fitted.values[model$model$AOI.condition %in% AOIcond_list[[i]]]) })
  ylim <- c(0, max(max_fitted))
  
  for (j in seq_len(length(response_vec))) {
    res_list <- DPKogHelpers::plot_predictions(model = model_list[[i]][[j]],
                                               nLetters.WD = nLetters.WD,
                                               groups_vec = AOIcond_list[[i]],
                                               model_name = paste("Model",i),
                                               ylim = ylim,
                                               xlab_size = 2)
    # 1) Plot of model predictions
    ggsave(plot = res_list$plot, filename = paste0(path_resImages, "plot_m",i,"_", response_vec[j], ".pdf"), width = 8, height = 5)
    
    # 2) Table of estimates and model predictions
    write.csv(res_list$table, file = paste0(path_resTables, "tab_m",i,"_",response_vec[j],".csv"), row.names = TRUE)
  }
}

# results of model 1 for inspection:
res_mod1_TRT <- read.csv(paste0(path_resTables, "tab_m1_TRT.csv"))
res_mod1_FRT <- read.csv(paste0(path_resTables, "tab_m1_FRT.csv"))
res_mod1_RRT <- read.csv(paste0(path_resTables, "tab_m1_RRT.csv"))

#------------------------------

#-------MODEL 2 - AOI 1-------- 

# DATA PREPARATION
# ----------------
dat$mechanism <- substr(dat$condition, 1, 1)
dat$mechanism[dat$condition == "CE0"] <- "none"
dat$info <- substr(dat$condition, 2, 2) 
dat$info[dat$condition == "CE0"] <- 0

# factors with desired level order
dat$mechanism <- factor(dat$mechanism)
dat$mechanism <- relevel(dat$mechanism, ref = "C")
dat$info <- factor(dat$info, levels = c("1", "2", "3", "4", "5", "6", "7", "0"))
dat$info_bin <- factor(ifelse(dat$info == "1", "pronominal", "lexical"),
                       levels = c("pronominal", "lexical"))

# MODEL ESTIMATION
# ----------------
mod2_fit <- function(AOI, resp) {
  dat_mod <- dat[dat$AOI == AOI,]
  dat_mod[, paste0(resp, ".L")] <- dat_mod[, resp] / dat_mod$nLetters
  # model estimation
  message(paste0("Estimate model 2_", resp, "_", AOI, "..."))
  mod <- gam(data = dat_mod,
             formula = get(paste0(resp, ".L")) ~ mechanism + info +
               s(Participant, bs = "re") + s(Topic, bs = "re"))
  return(mod)
}

mod2_TRT_1 <- mod2_fit(1, "TRT")
mod2_FRT_1 <- mod2_fit(1, "FRT")
mod2_RRT_1 <- mod2_fit(1, "RRT")

res_mod2_TRT_1 <- round(anova.gam(mod2_TRT_1)$pTerms.table, 4)
res_mod2_FRT_1 <- round(anova.gam(mod2_FRT_1)$pTerms.table, 4)
res_mod2_RRT_1 <- round(anova.gam(mod2_RRT_1)$pTerms.table, 4)

#------------------------------

#-------MODEL 2 - AOI 6-------- 

# MODEL ESTIMATION
# ----------------

mod2_TRT_6 <- mod2_fit(6, "TRT")
mod2_FRT_6 <- mod2_fit(6, "FRT")
mod2_RRT_6 <- mod2_fit(6, "RRT")

res_mod2_TRT_6 <- round(anova.gam(mod2_TRT_6)$pTerms.table, 4)
res_mod2_FRT_6 <- round(anova.gam(mod2_FRT_6)$pTerms.table, 4)
res_mod2_RRT_6 <- round(anova.gam(mod2_RRT_6)$pTerms.table, 4)

#------------------------------