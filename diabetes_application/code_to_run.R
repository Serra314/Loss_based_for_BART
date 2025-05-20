source('code/import_functions.R')
source('utils_BART_result.R')
library(fastDummies)
library(data.tree)
library(truncnorm)
library(dplyr)
library(ggplot2)
mean(rnterm(1000, 1.56))
mean(rnterm.chip(1000, 0.25, 2))
sd(rnterm(1000, 1.56))
sd(rnterm.chip(1000, 0.25, 2))

diab_data <- read.csv(file = 'diabetes+130-us+hospitals+for+years+1999-2008/diabetic_data.csv',
                      header = TRUE)
diab_data.copy <- diab_data

head(diab_data)
summary(diab_data)

######################
## DATA PREPARATION ##
######################

num_vars <- c()
fact_vars <- c()
for(i in 1:ncol(diab_data)){
  col.name <- colnames(diab_data)[i]
  if(is.numeric(diab_data[,i])){
    num_vars <- c(num_vars, col.name)
  } 
  if(is.character(diab_data[,i])){
    diab_data.copy[,i] <- factor(diab_data[,i])
    fact_vars <- c(fact_vars, col.name)
  }
}

ncol(diab_data)
length(num_vars)
length(fact_vars)
# variables to keep and discard - 50 originally
to_discard <- c('encounter_id', # identifier
             'patient_nbr', # identifier
             'weight', # categorical mostly missing (98569)
             'admission_type_id', # categorical 9 
             'discharge_disposition_id', # categorical 29 
             'admission_source_id', # categorical 21 
             'payer_code', # categorical 23 
             'medical_specialty', # categorical 84
             'diag_1', #categorical 848
             'diag_2', # categorical 923
             'diag_3', # categorical 954
             'max_glu_serum', # categorical 3 mostly missing (~99000)
             'A1Cresult', # categorical mostly missing (84748)
             'acetohexamide', # categorical mostly one category (all except 1)
             'tolbutamide', # categorical mostly one category (all except 23)
             'troglitazone', # categorical mostly one category (all except 3)
             'tolazamide', # categorical mostly one category (all except 39)
             'examide', # categorical all the same
             'citoglipton', # categorical all the same
             'glipizide.metformin', # categorical mostly one category (all except 13)
             'glimepiride.pioglitazone', # categorical mostly one category (all except 1)
             'metformin.rosiglitazone'  # categorical mostly one category (all except 1)
             )

to_keep <- c('race', # categorical 5 plus ?
             'gender', # categorical 2 plus Unknown/Invalid
             'age', # categorical 10
             'time_in_hospital', # numerical (days)
             'num_lab_procedures', # numerical
             'num_procedures', # numerical
             'num_medications', # numerical
             'number_outpatient', # numerical
             'number_emergency', #numerical
             'number_inpatient', # numerical
             'number_diagnoses', # numerical
             'metformin', # categorical 4 turn in YES NO
             'repaglinide', # categorical 4 turn in YES NO 
             'nateglinide', # categorical 4 turn in YES NO
             'chlorpropamide', # categorical 4 turn in YES NO
             'glimepiride', # categorical 4 turn in YES NO
             'glipizide', # categorical 4 turn in YES NO
             'glyburide', # categorical 4 turn in YES NO
             'pioglitazone', # categorical 4 turn in YES NO
             'rosiglitazone', # categorical 4 turn in YES NO
             'acarbose', # categorical 4 turn in YES NO
             'miglitol', # categorical 4 turn in YES NO
             'insulin', # categorical 4 turn in YES NO
             'glyburide.metformin', # categorical 4 turn in YES NO
             'change', # categorical 2
             'diabetesMed', # categorical 2
             'readmitted' # categorical 3 turn into YES NO
             )

# select columns to keep
diab_data_keep <- diab_data.copy[,to_keep]
# 27 of original 50
ncol(diab_data_keep)

## remove missing values 

idx.missing.race = diab_data_keep$race == '?'
diab_data_keep_nomis = diab_data_keep[!idx.missing.race, ]
diab_data_keep_nomis$race <- factor(as.character(diab_data_keep_nomis$race))


idx.missing.gender = diab_data_keep_nomis$gender == 'Unknown/Invalid'
diab_data_keep_nomis = diab_data_keep_nomis[!idx.missing.gender, ]
diab_data_keep_nomis$gender <- factor(as.character(diab_data_keep_nomis$gender))

## turn age into only three categories (below 30, [30, 60), over 60

new_age <- as.character(diab_data_keep_nomis$age)
summary(diab_data_keep_nomis$age)

idx.below30 <- (new_age == '[0-10)' | new_age == '[10-20)' | new_age == '[20-30)')
idx.30_60 <- (new_age == '[30-40)' | new_age == '[40-50)' | new_age == '[50-60)')
idx.over60 <- (new_age == '[60-70)' | new_age == '[70-80)' | 
                 new_age == '[80-90)' | new_age == '[90-100)')
new_age[idx.below30] <- 'below30'
new_age[idx.30_60] <- '[30-60)'
new_age[idx.over60] <- 'over60'

diab_data_keep_nomis$age <- factor(new_age)
summary(diab_data_keep_nomis)


# turn medication data in yes or no

medication_vars <- c('metformin', # categorical 4 turn in YES NO
                     'repaglinide', # categorical 4 turn in YES NO 
                     'nateglinide', # categorical 4 turn in YES NO
                     'chlorpropamide', # categorical 4 turn in YES NO
                     'glimepiride', # categorical 4 turn in YES NO
                     'glipizide', # categorical 4 turn in YES NO
                     'glyburide', # categorical 4 turn in YES NO
                     'pioglitazone', # categorical 4 turn in YES NO
                     'rosiglitazone', # categorical 4 turn in YES NO
                     'acarbose', # categorical 4 turn in YES NO
                     'miglitol', # categorical 4 turn in YES NO
                     'insulin', # categorical 4 turn in YES NO
                     'glyburide.metformin' # categorical 4 turn in YES NO
                     )
df_med_var <- diab_data_keep_nomis[,medication_vars]
# turn the variable in 1 if that medication and 0 is No
df_med_var_new <- apply(df_med_var, 2, \(x) as.numeric(x != 'No') )

diab_data_keep_nomis_refact <-diab_data_keep_nomis 
diab_data_keep_nomis_refact[,medication_vars] <- df_med_var_new

summary(diab_data_keep_nomis_refact)

# code gender as 1 for MALE 0 for FEMALE
# code change as 1 for Ch and 0 for no
# code diabetesMed as 1 for Yes na 0 for No
diab_data_keep_nomis_refact$gender = as.numeric(diab_data_keep_nomis_refact$gender == 'Male')
diab_data_keep_nomis_refact$change = as.numeric(diab_data_keep_nomis_refact$change == 'Ch')
diab_data_keep_nomis_refact$diabetesMed = as.numeric(diab_data_keep_nomis_refact$diabetesMed == 'Yes')

summary(diab_data_keep_nomis_refact)

# refactor race, age as dummies
diab_data_keep_nomis_refact <- dummy_cols(diab_data_keep_nomis_refact, 
                                          select_columns = "race", 
                                          remove_selected_columns = TRUE)

diab_data_keep_nomis_refact <- dummy_cols(diab_data_keep_nomis_refact, 
                                          select_columns = "age", 
                                          remove_selected_columns = TRUE)

#######################################
## Initialise data for model fitting ##
#######################################

# mode for males - 45917 observations
diab_final_males <- diab_data_keep_nomis_refact[diab_data_keep_nomis_refact$gender == 1,]

df_to_prep = diab_final_males

# separate Y of interest from covariates used as predictors
yy <- as.numeric(df_to_prep$readmitted != 'NO')

X_pred <- df_to_prep
X_pred$readmitted <- NULL
X_pred$gender <- NULL
X_pred <- apply(X_pred, 2, \(x) (x - min(x))/(max(x) - min(x)))
summary(X_pred)
nrow(X_pred) 


# turn X_pred into a numeric matrix for model fitting
X_pred <- as.matrix(X_pred)
class(X_pred) <- 'numeric'
rownames(X_pred) <- as.character(1:nrow(X_pred))
colnames(X_pred) <- as.character(1:ncol(X_pred)) 
head(X_pred)


# MCMC params 
# terminal node values priors - gaussian
mu.prior.mean = 0
mu.prior.var = 3


# params for computation
include.split <- FALSE
cont.unif <- TRUE
moves.prob <- c(0.4, 0.4, 0.1, 0.1)
#n.tree = 10
n.iter = 1500 # 1500 - time x 10, 20, 40 alberi 
# meeting 16 agosto 9.00
lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.56, 0.62))
chip.prior.rep <- list(fun = chipman_prior_tree, param = c(0.25, 2))

#priors check
mean(rnterm(1000, 1.56))
mean(rnterm.chip(1000, 0.25, 2))
sd(rnterm(1000, 1.56))
sd(rnterm.chip(1000, 0.25, 2))


## fit the models

loss_based_10 <- MCMC_for_BART_binary(n.tree = 10, 
                                      n.iter = n.iter, 
                                      Y = yy, X = X_pred, 
                                      include.split = include.split, 
                                      cont.unif = cont.unif, 
                                      mean.prior = 0, 
                                      var.prior = 3, 
                                      moves.prob = moves.prob, 
                                      prior.list = lb.prior.def)
save(loss_based_10, file = 'DIAB_lb_def_10.Rds')

loss_based_20 <- MCMC_for_BART_binary(n.tree = 20, 
                                      n.iter = n.iter, 
                                      Y = yy, X = X_pred, 
                                      include.split = include.split, 
                                      cont.unif = cont.unif, 
                                      mean.prior = 0, 
                                      var.prior = 3, 
                                      moves.prob = moves.prob, 
                                      prior.list = lb.prior.def)
save(loss_based_20, file = 'DIAB_lb_def_20.Rds')



loss_based_40 <- MCMC_for_BART_binary(n.tree = 40, 
                                      n.iter = n.iter, 
                                      Y = yy, X = X_pred, 
                                      include.split = include.split, 
                                      cont.unif = cont.unif, 
                                      mean.prior = 0, 
                                      var.prior = 3, 
                                      moves.prob = moves.prob, 
                                      prior.list = lb.prior.def)
save(loss_based_40, file = 'DIAB_lb_def_40.Rds')





  
chip_10 <- MCMC_for_BART_binary(n.tree = 10, 
                                      n.iter = n.iter, 
                                      Y = yy, X = X_pred, 
                                      include.split = include.split, 
                                      cont.unif = cont.unif, 
                                      mean.prior = 0, 
                                      var.prior = 3, 
                                      moves.prob = moves.prob, 
                                      prior.list = chip.prior.rep)
save(chip_10, file = 'DIAB_chip_10.Rds')

chip_20 <- MCMC_for_BART_binary(n.tree = 20, 
                                n.iter = n.iter, 
                                Y = yy, X = X_pred, 
                                include.split = include.split, 
                                cont.unif = cont.unif, 
                                mean.prior = 0, 
                                var.prior = 3, 
                                moves.prob = moves.prob, 
                                prior.list = chip.prior.rep)
save(chip_20, file = 'DIAB_chip_20.Rds')

chip_40 <- MCMC_for_BART_binary(n.tree = 40, 
                                n.iter = n.iter, 
                                Y = yy, X = X_pred, 
                                include.split = include.split, 
                                cont.unif = cont.unif, 
                                mean.prior = 0, 
                                var.prior = 3, 
                                moves.prob = moves.prob, 
                                prior.list = chip.prior.rep)
save(chip_40, file = 'DIAB_chip_40.Rds')



load('DIAB_lb_def_10.Rds')
load('DIAB_lb_def_20.Rds')
load('DIAB_lb_def_40.Rds')

load('DIAB_chip_10.Rds')
load('DIAB_chip_20.Rds')
load('DIAB_chip_40.Rds')
## times
dd_time <- data.frame(prior = rep(c('loss based', 'chipman'), each = 3),
                      ntree = rep(c(10,20,40), 2),
                      niter = 1500,
                      time_hours = c(as.numeric(loss_based_10$time, unit = "hours"),
                                     as.numeric(loss_based_20$time, unit = "hours"),
                                     as.numeric(loss_based_40$time, unit = "hours"),
                                     as.numeric(chip_10$time, unit = "hours"),
                                     as.numeric(chip_20$time, unit = "hours"),
                                     as.numeric(chip_40$time, unit = "hours")))




## NUMBER OF TERMINAL NODES




nt_lb_10.df <- nterm_BART(loss_based_10)
nt_lb_10.df$prior = 'LB ~ (m == 10)'
nt_ch_10.df <- nterm_BART(chip_10)
nt_ch_10.df$prior = 'CH ~ (m == 10)'

nt_10.df <- rbind(nt_lb_10.df, nt_ch_10.df)


nt_lb_20.df <- nterm_BART(loss_based_20)
nt_lb_20.df$prior = 'LB ~ (m == 20)'
nt_ch_20.df <- nterm_BART(chip_20)
nt_ch_20.df$prior = 'CH ~ (m == 20)'

nt_20.df <- rbind(nt_lb_20.df, nt_ch_20.df)


nt_lb_40.df <- nterm_BART(loss_based_40)
nt_lb_40.df$prior = 'LB ~ (m == 40)'
nt_ch_40.df <- nterm_BART(chip_40)
nt_ch_40.df$prior = 'CH ~ (m == 40)'

nt_40.df <- rbind(nt_lb_40.df, nt_ch_40.df)

library(ggplot2)
nt_10.df_burn <- nt_10.df[nt_10.df$idx >= 750,]
nt_20.df_burn <- nt_20.df[nt_20.df$idx >= 750,]
nt_40.df_burn <- nt_40.df[nt_40.df$idx >= 750,]


pl_nl_10 <- ggplot(nt_10.df_burn, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('') + 
  facet_wrap(facets = ~ prior, labeller = label_parsed) + 
  theme_classic() + 
  theme(panel.grid.major = element_line())

pl_nl_20 <-ggplot(nt_20.df_burn, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Number of terminal nodes') + 
  facet_wrap(facets = ~ prior, labeller = label_parsed) + 
  theme_classic() + 
  theme(panel.grid.major = element_line())

pl_nl_40 <-ggplot(nt_40.df_burn, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('') + 
  facet_wrap(facets = ~ prior, labeller = label_parsed) + 
  theme_classic() + 
  theme(panel.grid.major = element_line())

pdf('diab_10_20_40_num_term.pdf', width = 15, height = 10)
gridExtra::grid.arrange(pl_nl_10, pl_nl_20, pl_nl_40)
dev.off()
## depth 


## DEPTH


depth_lb_10.df <- depth_BART(loss_based_10)
depth_lb_10.df$prior = 'loss-based-10'
depth_ch_10.df <- depth_BART(chip_10)
depth_ch_10.df$prior = 'chipman-10'

depth_10.df <- rbind(depth_lb_10.df, depth_ch_10.df)


depth_lb_20.df <- depth_BART(loss_based_20)
depth_lb_20.df$prior = 'loss-based-20'
depth_ch_20.df <- depth_BART(chip_20)
depth_ch_20.df$prior = 'chipman-20'

depth_20.df <- rbind(depth_lb_20.df, depth_ch_20.df)


depth_lb_40.df <- depth_BART(loss_based_40)
depth_lb_40.df$prior = 'loss-based-40'
depth_ch_40.df <- depth_BART(chip_40)
depth_ch_40.df$prior = 'chipman-40'

depth_40.df <- rbind(depth_lb_40.df, depth_ch_40.df)


ggplot(depth_10.df, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Depth') + 
  facet_wrap(facets = ~ prior)

ggplot(depth_20.df, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Depth') + 
  facet_wrap(facets = ~ prior)

ggplot(depth_40.df, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Depth') + 
  facet_wrap(facets = ~ prior)





## PREDICTIONS AT OBSERVED POINTS





list_pred_lb_10 <- lapply(1:length(loss_based_10$trees), \(idx) BART_calculate_pred(loss_based_10$trees[[idx]], X_pred))
save(list_pred_lb_10, file = 'list_pred_lb_10.Rds')


list_pred_lb_20 <- lapply(1:length(loss_based_20$trees), \(idx) BART_calculate_pred(loss_based_20$trees[[idx]], X_pred))
save(list_pred_lb_20, file = 'list_pred_lb_20.Rds')


list_pred_lb_40 <- lapply(1:length(loss_based_40$trees),\(idx) BART_calculate_pred(loss_based_40$trees[[idx]], X_pred))
save(list_pred_lb_40, file = 'list_pred_lb_40.Rds')


list_pred_ch_10 <- lapply(1:length(chip_10$trees),\(idx) BART_calculate_pred(chip_10$trees[[idx]], X_pred))
save(list_pred_ch_10, file = 'list_pred_ch_10.Rds')


list_pred_ch_20 <- lapply(1:length(chip_20$trees),\(idx) BART_calculate_pred(chip_20$trees[[idx]], X_pred))
save(list_pred_ch_20, file = 'list_pred_ch_20.Rds')


list_pred_ch_40 <- lapply(1:length(chip_40$trees),\(idx) BART_calculate_pred(chip_40$trees[[idx]], X_pred))
save(list_pred_ch_40, file = 'list_pred_ch_40.Rds')


###### BINARY LOGLIKELIHOOD
load('list_pred_lb_10.Rds')
load('list_pred_lb_20.Rds')
load('list_pred_lb_40.Rds')

load('list_pred_ch_10.Rds')
load('list_pred_ch_20.Rds')
load('list_pred_ch_40.Rds')




loglik_lb_10 <- vapply(list_pred_lb_10, \(x) binary_log_lik_from_pred(x, yy), 0)
loglik_lb_20 <- vapply(list_pred_lb_20, \(x) binary_log_lik_from_pred(x, yy), 0)
loglik_lb_40 <- vapply(list_pred_lb_40, \(x) binary_log_lik_from_pred(x, yy), 0)

loglik_ch_10 <- vapply(list_pred_ch_10, \(x) binary_log_lik_from_pred(x, yy), 0)
loglik_ch_20 <- vapply(list_pred_ch_20, \(x) binary_log_lik_from_pred(x, yy), 0)
loglik_ch_40 <- vapply(list_pred_ch_40, \(x) binary_log_lik_from_pred(x, yy), 0)



df_log_lik_lb <- rbind(data.frame(iter = 1:n.iter,
                                  loglik = loglik_lb_10,
                                  ntree = 10,
                                  prior = 'loss-bases'),
                       data.frame(iter = 1:n.iter,
                                  loglik = loglik_lb_20,
                                  ntree = 20,
                                  prior = 'loss-bases'),
                       data.frame(iter = 1:n.iter,
                                  loglik = loglik_lb_40,
                                  ntree = 40,
                                  prior = 'loss-bases'))



df_log_lik_ch <- rbind(data.frame(iter = 1:n.iter,
                                  loglik = loglik_ch_10,
                                  ntree = 10,
                                  prior = 'chipman'),
                       data.frame(iter = 1:n.iter,
                                  loglik = loglik_ch_20,
                                  ntree = 20,
                                  prior = 'chipman'),
                       data.frame(iter = 1:n.iter,
                                  loglik = loglik_ch_40,
                                  ntree = 40,
                                  prior = 'chipman'))


df_plot_loglik <- rbind(df_log_lik_ch, df_log_lik_lb)#df_log_lik_lb#


ggplot(df_plot_loglik, aes(iter, loglik, color = factor(ntree))) + 
  geom_line() + 
  facet_wrap(facets = ~ prior) + 
  theme_bw()

ggplot(df_plot_loglik, aes(iter, loglik, color = factor(ntree), linetype = factor(prior))) + 
  geom_line()

# ACCURACY



acc_lb_10 <- vapply(list_pred_lb_10, \(x) accuracy_from_pred(x, yy), 0)
acc_lb_20 <- vapply(list_pred_lb_20, \(x) accuracy_from_pred(x, yy), 0)
acc_lb_40 <- vapply(list_pred_lb_40, \(x) accuracy_from_pred(x, yy), 0)

acc_ch_10 <- vapply(list_pred_ch_10, \(x) accuracy_from_pred(x, yy), 0)
acc_ch_20 <- vapply(list_pred_ch_20, \(x) accuracy_from_pred(x, yy), 0)
acc_ch_40 <- vapply(list_pred_ch_40, \(x) accuracy_from_pred(x, yy), 0)



df_acc_lb <- rbind(data.frame(iter = 1:n.iter,
                              accuracy = acc_lb_10,
                                  ntree = 10,
                                  prior = 'loss-bases'),
                       data.frame(iter = 1:n.iter,
                                  accuracy = acc_lb_20,
                                  ntree = 20,
                                  prior = 'loss-bases'),
                       data.frame(iter = 1:n.iter,
                                  accuracy = acc_lb_40,
                                  ntree = 40,
                                  prior = 'loss-bases'))



df_acc_ch <- rbind(data.frame(iter = 1:n.iter,
                              accuracy = acc_ch_10,
                                  ntree = 10,
                                  prior = 'chipman'),
                       data.frame(iter = 1:n.iter,
                                  accuracy = acc_ch_20,
                                  ntree = 20,
                                  prior = 'chipman'),
                       data.frame(iter = 1:n.iter,
                                  accuracy = acc_ch_40,
                                  ntree = 40,
                                  prior = 'chipman'))


df_plot_acc <- rbind(df_acc_ch, df_acc_lb)
df_plot_acc$loglik <- df_plot_loglik$loglik

ggplot(df_plot_acc, aes(iter, accuracy, color = factor(ntree))) + 
  geom_line() + 
  facet_wrap(facets = ~ prior)

ggplot(df_plot_acc[df_plot_acc$iter >= 750,], aes(loglik, accuracy, color = factor(ntree))) + 
  geom_point() + 
  facet_wrap(facets = ~ prior)








## rerun models with 10 trees to check if any difference
library(parallel)

runner_lb_10 <- function(idx){
  
  MCMC_for_BART_binary(n.tree = 10, 
                       n.iter = n.iter, 
                       Y = yy, X = X_pred, 
                       include.split = include.split, 
                       cont.unif = cont.unif, 
                       mean.prior = 0, 
                       var.prior = 3, 
                       moves.prob = moves.prob, 
                       prior.list = lb.prior.def)
}
n.cores <- 7
clust <- makeCluster(n.cores)
clusterExport(clust, varlist = c(ls(), 'rtruncnorm'))
list_models_lb <- parLapply(clust, 1:10, runner_lb_10)
stopCluster(clust)
save(list_models_lb, file = 'DIAB_list_models_lb_10.Rds')



runner_lb_20 <- function(idx){
  
  MCMC_for_BART_binary(n.tree = 20, 
                       n.iter = n.iter, 
                       Y = yy, X = X_pred, 
                       include.split = include.split, 
                       cont.unif = cont.unif, 
                       mean.prior = 0, 
                       var.prior = 3, 
                       moves.prob = moves.prob, 
                       prior.list = lb.prior.def)
}
n.cores <- 10
clust <- makeCluster(n.cores)
clusterExport(clust, varlist = c(ls(), 'rtruncnorm'))
list_models_lb_20_1to5 <- parLapply(clust, 1:5, runner_lb_20)
stopCluster(clust)
save(list_models_lb_20_1to5, file = 'DIAB_list_models_lb_20_chain1to5.Rds')
rm(list_models_lb_20_1to5)

clust <- makeCluster(n.cores)
clusterExport(clust, varlist = c(ls(), 'rtruncnorm'))
list_models_lb_20_6to10 <- parLapply(clust, 1:5, runner_lb_20)
stopCluster(clust)
save(list_models_lb_20_6to10, file = 'DIAB_list_models_lb_20_chain6to10.Rds')








runner_ch <- function(idx){
  MCMC_for_BART_binary(n.tree = 10, 
                       n.iter = n.iter, 
                       Y = yy, X = X_pred, 
                       include.split = include.split, 
                       cont.unif = cont.unif, 
                       mean.prior = 0, 
                       var.prior = 3, 
                       moves.prob = moves.prob, 
                       prior.list = chip.prior.rep)
}
n.cores <- 12
stt = Sys.time()
clust <- makeCluster(n.cores)
clusterExport(clust, varlist = c(ls(), 'rtruncnorm'))
list_models_ch <- parLapply(clust, 1:10, runner_ch)
stopCluster(clust)
ett = Sys.time() - stt
save(list_models_ch, file = 'DIAB_list_models_ch_10.Rds')


runner_ch_20 <- function(idx){
  MCMC_for_BART_binary(n.tree = 20, 
                       n.iter = n.iter, 
                       Y = yy, X = X_pred, 
                       include.split = include.split, 
                       cont.unif = cont.unif, 
                       mean.prior = 0, 
                       var.prior = 3, 
                       moves.prob = moves.prob, 
                       prior.list = chip.prior.rep)
}
n.cores <- 10
clust <- makeCluster(n.cores)
clusterExport(clust, varlist = c(ls(), 'rtruncnorm'))
list_models_ch_20_1to5 <- parLapply(clust, 1:5, runner_ch_20)
stopCluster(clust)
save(list_models_ch_20_1to5, file = 'DIAB_list_models_ch_20_1to5.Rds')
rm(list_models_ch_20_1to5)


clust <- makeCluster(n.cores)
clusterExport(clust, varlist = c(ls(), 'rtruncnorm'))
list_models_ch_20_6to10 <- parLapply(clust, 1:5, runner_ch_20)
stopCluster(clust)
save(list_models_ch_20_6to10, file = 'DIAB_list_models_ch_20_6to10.Rds')
rm(list_models_ch_20_6to10)


###################
# NTERM analysis ##
###################


load('DIAB_list_models_ch_20_1to5.Rds')
load('DIAB_list_models_ch_20_6to10.Rds')
load('DIAB_list_models_lb_20_chain1to5.Rds')
load('DIAB_list_models_lb_20_chain6to10.Rds')

library(dplyr)
library(ggplot2)

##LB##

nterm_lb_1to5 <- lapply(1:length(list_models_lb_20_1to5), 
                   \(idxx) nterm_BART(list_models_lb_20_1to5[[idxx]]) %>% 
                     mutate(chain = idxx) )

nterm_lb_6to10 <- lapply(1:length(list_models_lb_20_1to5), 
                        \(idxx) nterm_BART(list_models_lb_20_1to5[[idxx]]) %>% 
                          mutate(chain = 5 + idxx) )



nterm_lb.df_1to5 <- Reduce(rbind, nterm_lb_1to5)
nterm_lb.df_6to10 <- Reduce(rbind, nterm_lb_6to10)
nterm_lb.df_20 <- rbind(nterm_lb.df_1to5, nterm_lb.df_6to10)


ggplot(nterm_lb.df_20, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Number of terminal nodes') + 
  facet_wrap(facets = ~ chain, ncol = 2)


nterm_lb_av_nn_1_to_5 <- lapply(1:length(nterm_lb.df_1to5), \(idxx) data.frame(iter = 1:1500,
                                                                               av_nn = aggregate(nterm_lb_1to5[[idxx]]$nn, 
                                                                                                 list(factor(nterm_lb_1to5[[idxx]]$idx)), mean)[,2],
                                                                               chain = idxx))

nterm_lb_av_nn_6_to_10 <- lapply(1:length(nterm_lb.df_6to10), \(idxx) data.frame(iter = 1:1500,
                                                                               av_nn = aggregate(nterm_lb_6to10[[idxx]]$nn, 
                                                                                                 list(factor(nterm_lb_6to10[[idxx]]$idx)), mean)[,2],
                                                                               chain = 5 + idxx))


nterm_lb_av_nn_df_1to5 <- Reduce(rbind, nterm_lb_av_nn_1_to_5)
nterm_lb_av_nn_df_6to10 <- Reduce(rbind, nterm_lb_av_nn_6_to_10)
nterm_lb_av_nn_df <- rbind(nterm_lb_av_nn_df_1to5, nterm_lb_av_nn_df_6to10)
nterm_lb_av_nn_df$prior = 'LB'

ggplot(nterm_lb_av_nn_df, aes(iter, av_nn, color = factor(chain))) + 
  geom_line()

## CH##
nterm_ch_1to5 <- lapply(1:length(list_models_ch_20_1to5), 
                        \(idxx) nterm_BART(list_models_ch_20_1to5[[idxx]]) %>% 
                          mutate(chain = idxx) )

nterm_ch_6to10 <- lapply(1:length(list_models_ch_20_6to10), 
                         \(idxx) nterm_BART(list_models_ch_20_6to10[[idxx]]) %>% 
                           mutate(chain = 5 + idxx) )



nterm_ch.df_1to5 <- Reduce(rbind, nterm_ch_1to5)
nterm_ch.df_6to10 <- Reduce(rbind, nterm_ch_6to10)
nterm_ch.df_20 <- rbind(nterm_ch.df_1to5, nterm_ch.df_6to10)


ggplot(nterm_ch.df_20, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Number of terminal nodes') + 
  facet_wrap(facets = ~ chain, ncol = 2)


nterm_ch_av_nn_1_to_5 <- lapply(1:length(nterm_ch.df_1to5), \(idxx) data.frame(iter = 1:1500,
                                                                               av_nn = aggregate(nterm_ch_1to5[[idxx]]$nn, 
                                                                                                 list(factor(nterm_ch_1to5[[idxx]]$idx)), mean)[,2],
                                                                               chain = idxx))

nterm_ch_av_nn_6_to_10 <- lapply(1:length(nterm_ch.df_6to10), \(idxx) data.frame(iter = 1:1500,
                                                                                 av_nn = aggregate(nterm_ch_6to10[[idxx]]$nn, 
                                                                                                   list(factor(nterm_ch_6to10[[idxx]]$idx)), mean)[,2],
                                                                                 chain = 5 + idxx))


nterm_ch_av_nn_df_1to5 <- Reduce(rbind, nterm_ch_av_nn_1_to_5)
nterm_ch_av_nn_df_6to10 <- Reduce(rbind, nterm_ch_av_nn_6_to_10)
nterm_ch_av_nn_df <- rbind(nterm_ch_av_nn_df_1to5, nterm_ch_av_nn_df_6to10)
nterm_ch_av_nn_df$prior = 'CH'

ggplot(nterm_ch_av_nn_df, aes(iter, av_nn, color = factor(chain))) + 
  geom_line()

ggplot(rbind(nterm_ch_av_nn_df, nterm_lb_av_nn_df), aes(iter, av_nn, color = factor(chain))) + 
  geom_line() + 
  facet_wrap(facets = ~ prior) #+ theme(legend.position = 'none')












load('DIAB_list_models_lb_10.Rds')
load('DIAB_list_models_ch_10.Rds')

library(dplyr)
nterm_lb <- lapply(1:length(list_models_lb), 
                   \(idxx) nterm_BART(list_models_lb[[idxx]]) %>% 
                     mutate(chain = idxx) )

nterm_lb.df <- Reduce(rbind, nterm_lb)
nterm_lb.df$chain <- paste0('LB - ',nterm_lb.df$chain)
nterm_lb.df$prior = 'LB'

library(ggplot2)
ggplot(nterm_lb.df, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Number of terminal nodes') + 
  facet_wrap(facets = ~ chain, ncol = 2)

nterm_ch <- lapply(1:length(list_models_ch), 
                   \(idxx) nterm_BART(list_models_ch[[idxx]]) %>% 
                     mutate(chain = idxx) )

nterm_ch.df <- Reduce(rbind, nterm_ch)
nterm_ch.df$chain <- paste0('CH - ',nterm_ch.df$chain)
nterm_ch.df$prior = 'CH'

ggplot(nterm_ch.df, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Number of terminal nodes') + 
  facet_wrap(facets = ~ chain, ncol = 2)



nterm_df_10 <- rbind(nterm_ch.df, nterm_lb.df)
ggplot(nterm_df_10, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Number of terminal nodes') + 
  facet_wrap(facets = ~ prior, ncol = 2)


nterm_df_10$chain <- factor(nterm_df_10$chain, 
                            levels = c('CH - 1', 'LB - 1',
                                       'CH - 2', 'LB - 2',
                                       'CH - 3', 'LB - 3',
                                       'CH - 4', 'LB - 4',
                                       'CH - 5', 'LB - 5',
                                       'CH - 6', 'LB - 6',
                                       'CH - 7', 'LB - 7',
                                       'CH - 8', 'LB - 8',
                                       'CH - 9', 'LB - 9',
                                       'CH - 10', 'LB - 10'))


ggplot(nterm_df_10, aes(trees, nn)) + 
  geom_boxplot() + 
  ylab('Number of terminal nodes') + 
  facet_wrap(facets = ~ chain, ncol = 2)


load('DIAB_list_models_lb_10.Rds')

list_pred_lb_10_10chain <- lapply(list_models_lb, 
                                  \(model_list) lapply(1:length(model_list$trees), \(idx) BART_calculate_pred(model_list$trees[[idx]], X_pred)))
  
save(list_pred_lb_10_10chain, file = 'list_pred_lb_10_10chain.Rds')

rm(list_models_lb)
rm(list_pred_lb_10_10chain)





load('DIAB_list_models_ch_10.Rds')

list_pred_ch_10_10chain <- lapply(list_models_ch, 
                                  \(model_list) lapply(1:length(model_list$trees), \(idx) BART_calculate_pred(model_list$trees[[idx]], X_pred)))

save(list_pred_ch_10_10chain, file = 'list_pred_ch_10_10chain.Rds')



## pred 20 trees

## LB
load('DIAB_list_models_lb_20_chain1to5.Rds')
list_pred_lb_20_1to5 <- lapply(list_models_lb_20_1to5, 
                               \(model_list) lapply(1:length(model_list$trees), \(idx) BART_calculate_pred(model_list$trees[[idx]], X_pred)))

save(list_pred_lb_20_1to5, file = 'list_pred_lb_20_1to5.Rds')

rm(list_models_lb_20_1to5)
rm(list_pred_lb_20_1to5)


load('DIAB_list_models_lb_20_chain6to10.Rds')
list_pred_lb_20_6to10 <- lapply(list_models_lb_20_6to10, 
                                \(model_list) lapply(1:length(model_list$trees), \(idx) BART_calculate_pred(model_list$trees[[idx]], X_pred)))

save(list_pred_lb_20_6to10, file = 'list_pred_lb_20_6to10.Rds')

rm(list_models_lb_20_6to10)
rm(list_pred_lb_20_6to10)

## CH

load('DIAB_list_models_ch_20_1to5.Rds')
list_pred_ch_20_1to5 <- lapply(list_models_ch_20_1to5, 
                               \(model_list) lapply(1:length(model_list$trees), \(idx) BART_calculate_pred(model_list$trees[[idx]], X_pred)))

save(list_pred_ch_20_1to5, file = 'list_pred_ch_20_1to5.Rds')

rm(list_models_ch_20_1to5)
rm(list_pred_ch_20_1to5)


load('DIAB_list_models_ch_20_6to10.Rds')
list_pred_ch_20_6to10 <- lapply(list_models_ch_20_6to10, 
                                \(model_list) lapply(1:length(model_list$trees), \(idx) BART_calculate_pred(model_list$trees[[idx]], X_pred)))

save(list_pred_ch_20_6to10, file = 'list_pred_ch_20_6to10.Rds')

rm(list_models_ch_20_6to10)
rm(list_pred_ch_20_6to10)


################ 20 chains analysis.

## likelihood - LB
load('list_pred_lb_20_1to5.Rds')
load('list_pred_lb_20_6to10.Rds')


loglik_lb_20chain_1to5 <- lapply(list_pred_lb_20_1to5, \(model_list) 
                            vapply(model_list, \(x) binary_log_lik_from_pred(x, yy), 0))

loglik_lb_20chain_6to10 <- lapply(list_pred_lb_20_6to10, \(model_list) 
                                 vapply(model_list, \(x) binary_log_lik_from_pred(x, yy), 0))


loglik_lb_20chain_1to5.todf <- lapply(1:length(loglik_lb_20chain_1to5), \(idxx) 
                                 data.frame(iter = 1:1500,
                                            loglik = loglik_lb_20chain_1to5[[idxx]],
                                            chain = paste0('LB - ', idxx),
                                            chain_num = idxx,
                                            prior = 'LB') 
)

loglik_lb_20chain_6to10.todf <- lapply(1:length(loglik_lb_20chain_6to10), \(idxx) 
                                      data.frame(iter = 1:1500,
                                                 loglik = loglik_lb_20chain_6to10[[idxx]],
                                                 chain = paste0('LB - ', 5 + idxx),
                                                 chain_num = 5 + idxx,
                                                 prior = 'LB') 
)



loglik_lb_20chain.df <- rbind(Reduce(rbind, loglik_lb_20chain_1to5.todf),
                              Reduce(rbind, loglik_lb_20chain_6to10.todf))


ggplot(loglik_lb_20chain.df, aes(iter, loglik, color = chain)) + 
  geom_line() + 
  theme(legend.position = 'none')


## accuracy - LB

acc_lb_20chain_1to5 <-  lapply(list_pred_lb_20_1to5, \(model_list) 
                               vapply(model_list, \(x) accuracy_from_pred(x, yy), 0))



acc_lb_20chain_6to10 <-  lapply(list_pred_lb_20_6to10, \(model_list) 
                               vapply(model_list, \(x) accuracy_from_pred(x, yy), 0))


acc_lb_20 <- c(unlist(acc_lb_20chain_1to5), unlist(acc_lb_20chain_6to10))

loglik_lb_20chain.df$accuracy <- acc_lb_20

rm(list_pred_lb_20_1to5)
rm(list_pred_lb_20_6to10)

## number of terminal nodes - LB

load('DIAB_list_models_lb_20_chain1to5.Rds')
load('DIAB_list_models_lb_20_chain6to10.Rds')

nterm_lb_20_1to5 <- lapply(1:length(list_models_lb_20_1to5), 
                   \(idxx) nterm_BART(list_models_lb_20_1to5[[idxx]]) %>% 
                     mutate(chain = idxx) )

nterm_lb_20_6to10 <- lapply(1:length(list_models_lb_20_6to10), 
                           \(idxx) nterm_BART(list_models_lb_20_6to10[[idxx]]) %>% 
                             mutate(chain = 5 + idxx) )


# average number of terminal nodes per iteration
nterm_lb_20_1to5_average <- lapply(nterm_lb_20_1to5, \(df) aggregate(df$nn, list(df$idx), mean)[,2] )
nterm_lb_20_6to10_average <- lapply(nterm_lb_20_6to10, \(df) aggregate(df$nn, list(df$idx), mean)[,2] )

# max number of terminal nodes per iteration
nterm_lb_20_1to5_max <- lapply(nterm_lb_20_1to5, \(df) aggregate(df$nn, list(df$idx), max)[,2] )
nterm_lb_20_6to10_max <- lapply(nterm_lb_20_6to10, \(df) aggregate(df$nn, list(df$idx), max)[,2] )

# median number of terminal nodes per iteration
nterm_lb_20_1to5_q0.05 <- lapply(nterm_lb_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.05))[,2] )
nterm_lb_20_6to10_q0.05 <- lapply(nterm_lb_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.05))[,2] )

nterm_lb_20_1to5_q0.95 <- lapply(nterm_lb_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.95))[,2] )
nterm_lb_20_6to10_q0.95 <- lapply(nterm_lb_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.95))[,2] )

nterm_lb_20_1to5_q0.75 <- lapply(nterm_lb_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.75))[,2] )
nterm_lb_20_6to10_q0.75 <- lapply(nterm_lb_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.75))[,2] )

nterm_lb_20_1to5_q0.25 <- lapply(nterm_lb_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.25))[,2] )
nterm_lb_20_6to10_q0.25 <- lapply(nterm_lb_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.25))[,2] )

nterm_lb_20_1to5_q0.5 <- lapply(nterm_lb_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.5))[,2] )
nterm_lb_20_6to10_q0.5 <- lapply(nterm_lb_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.5))[,2] )


loglik_lb_20chain.df$average_nn <- c(unlist(nterm_lb_20_1to5_average), unlist(nterm_lb_20_6to10_average) )
loglik_lb_20chain.df$max_nn <- c(unlist(nterm_lb_20_1to5_max), unlist(nterm_lb_20_6to10_max) )
loglik_lb_20chain.df$q0.05_nn <- c(unlist(nterm_lb_20_1to5_q0.05), unlist(nterm_lb_20_6to10_q0.05) )
loglik_lb_20chain.df$q0.95_nn <- c(unlist(nterm_lb_20_1to5_q0.95), unlist(nterm_lb_20_6to10_q0.95) )
loglik_lb_20chain.df$q0.75_nn <- c(unlist(nterm_lb_20_1to5_q0.75), unlist(nterm_lb_20_6to10_q0.75) )
loglik_lb_20chain.df$q0.25_nn <- c(unlist(nterm_lb_20_1to5_q0.25), unlist(nterm_lb_20_6to10_q0.25) )
loglik_lb_20chain.df$q0.5_nn <- c(unlist(nterm_lb_20_1to5_q0.5), unlist(nterm_lb_20_6to10_q0.5) )



ggplot(loglik_lb_20chain.df) + 
  geom_line(aes(iter, average_nn, linetype = chain), color = 'black') + 
  geom_line(aes(iter, q0.05_nn, linetype = chain), color = 'red') +
  geom_line(aes(iter, q0.95_nn, linetype = chain), color = 'red')




ggplot(loglik_lb_20chain.df[order(loglik_lb_20chain.df$average_nn),], aes(loglik, 1 - accuracy, color = average_nn)) + 
  geom_point()

ggplot(loglik_lb_20chain.df, aes(iter, average_nn, color = chain)) + 
  geom_line()

save(loglik_lb_20chain.df, file = 'result_df_20trees.Rds')
rm(list_models_lb_20_1to5)
rm(list_models_lb_20_6to10)


## likelihood - CH
load('list_pred_ch_20_1to5.Rds')
load('list_pred_ch_20_6to10.Rds')


loglik_ch_20chain_1to5 <- lapply(list_pred_ch_20_1to5, \(model_list) 
                                 vapply(model_list, \(x) binary_log_lik_from_pred(x, yy), 0))

loglik_ch_20chain_6to10 <- lapply(list_pred_ch_20_6to10, \(model_list) 
                                  vapply(model_list, \(x) binary_log_lik_from_pred(x, yy), 0))


loglik_ch_20chain_1to5.todf <- lapply(1:length(loglik_ch_20chain_1to5), \(idxx) 
                                      data.frame(iter = 1:1500,
                                                 loglik = loglik_ch_20chain_1to5[[idxx]],
                                                 chain = paste0('CH - ', idxx),
                                                 chain_num = idxx,
                                                 prior = 'CH') 
)

loglik_ch_20chain_6to10.todf <- lapply(1:length(loglik_ch_20chain_6to10), \(idxx) 
                                       data.frame(iter = 1:1500,
                                                  loglik = loglik_ch_20chain_6to10[[idxx]],
                                                  chain = paste0('CH - ', 5 + idxx),
                                                  chain_num = 5 + idxx,
                                                  prior = 'CH') 
)



loglik_ch_20chain.df <- rbind(Reduce(rbind, loglik_ch_20chain_1to5.todf),
                              Reduce(rbind, loglik_ch_20chain_6to10.todf))


ggplot(loglik_ch_20chain.df, aes(iter, loglik, color = chain)) + 
  geom_line() #+ theme(legend.position = 'none')


## accuracy - LB

acc_ch_20chain_1to5 <-  lapply(list_pred_ch_20_1to5, \(model_list) 
                               vapply(model_list, \(x) accuracy_from_pred(x, yy), 0))



acc_ch_20chain_6to10 <-  lapply(list_pred_ch_20_6to10, \(model_list) 
                                vapply(model_list, \(x) accuracy_from_pred(x, yy), 0))


acc_ch_20 <- c(unlist(acc_ch_20chain_1to5), unlist(acc_ch_20chain_6to10))

loglik_ch_20chain.df$accuracy <- acc_ch_20


rm(list_pred_ch_20_1to5)
rm(list_pred_ch_20_6to10)

## number of terminal nodes - LB

load('DIAB_list_models_ch_20_1to5.Rds')
load('DIAB_list_models_ch_20_6to10.Rds')

# number of terminal nodes
nterm_ch_20_1to5 <- lapply(1:length(list_models_ch_20_1to5), 
                           \(idxx) nterm_BART(list_models_ch_20_1to5[[idxx]]) %>% 
                             mutate(chain = idxx) )

nterm_ch_20_6to10 <- lapply(1:length(list_models_ch_20_6to10), 
                            \(idxx) nterm_BART(list_models_ch_20_6to10[[idxx]]) %>% 
                              mutate(chain = 5 + idxx) )



nterm_ch_20_1to5_average <- lapply(nterm_ch_20_1to5, \(df) aggregate(df$nn, list(df$idx), mean)[,2] )
nterm_ch_20_6to10_average <- lapply(nterm_ch_20_6to10, \(df) aggregate(df$nn, list(df$idx), mean)[,2] )


# max number of terminal nodes per iteration
nterm_ch_20_1to5_max <- lapply(nterm_ch_20_1to5, \(df) aggregate(df$nn, list(df$idx), max)[,2] )
nterm_ch_20_6to10_max <- lapply(nterm_ch_20_6to10, \(df) aggregate(df$nn, list(df$idx), max)[,2] )

# median number of terminal nodes per iteration
nterm_ch_20_1to5_q0.05 <- lapply(nterm_ch_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.05))[,2] )
nterm_ch_20_6to10_q0.05 <- lapply(nterm_ch_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.05))[,2] )

nterm_ch_20_1to5_q0.95 <- lapply(nterm_ch_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.95))[,2] )
nterm_ch_20_6to10_q0.95 <- lapply(nterm_ch_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.95))[,2] )

nterm_ch_20_1to5_q0.75 <- lapply(nterm_ch_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.75))[,2] )
nterm_ch_20_6to10_q0.75 <- lapply(nterm_ch_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.75))[,2] )

nterm_ch_20_1to5_q0.25 <- lapply(nterm_ch_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.25))[,2] )
nterm_ch_20_6to10_q0.25 <- lapply(nterm_ch_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.25))[,2] )

nterm_ch_20_1to5_q0.5 <- lapply(nterm_ch_20_1to5, \(df) aggregate(df$nn, list(df$idx), \(x) quantile(x, 0.5))[,2] )
nterm_ch_20_6to10_q0.5 <- lapply(nterm_ch_20_6to10, \(df) aggregate(df$nn, list(df$idx), \(x)quantile(x, 0.5))[,2] )

loglik_ch_20chain.df$average_nn <- c(unlist(nterm_ch_20_1to5_average), unlist(nterm_ch_20_6to10_average) )
loglik_ch_20chain.df$max_nn <- c(unlist(nterm_ch_20_1to5_max), unlist(nterm_ch_20_6to10_max) )
loglik_ch_20chain.df$q0.05_nn <- c(unlist(nterm_ch_20_1to5_q0.05), unlist(nterm_ch_20_6to10_q0.05) )
loglik_ch_20chain.df$q0.95_nn <- c(unlist(nterm_ch_20_1to5_q0.95), unlist(nterm_ch_20_6to10_q0.95) )
loglik_ch_20chain.df$q0.75_nn <- c(unlist(nterm_ch_20_1to5_q0.75), unlist(nterm_ch_20_6to10_q0.75) )
loglik_ch_20chain.df$q0.25_nn <- c(unlist(nterm_ch_20_1to5_q0.25), unlist(nterm_ch_20_6to10_q0.25) )
loglik_ch_20chain.df$q0.5_nn <- c(unlist(nterm_ch_20_1to5_q0.5), unlist(nterm_ch_20_6to10_q0.5) )




ggplot(loglik_ch_20chain.df[order(loglik_ch_20chain.df$average_nn),], aes(loglik, 1 - accuracy, color = average_nn)) + 
  geom_point()

ggplot(loglik_ch_20chain.df, aes(iter, average_nn, color = chain)) + 
  geom_line()

save(loglik_ch_20chain.df, file = 'result_ch_20trees.Rds')
rm(list_models_ch_20_1to5)
rm(list_models_ch_20_6to10)

## comparison LB vs CH



ggplot(rbind(loglik_lb_20chain.df, loglik_ch_20chain.df), 
       aes(iter, average_nn, color = factor(chain_num))) + 
  geom_line() + 
  facet_wrap(facets = ~ prior) + 
  theme(legend.position = 'none')

ggplot(rbind(loglik_lb_20chain.df, loglik_ch_20chain.df), 
       aes(iter, accuracy, color = factor(chain_num))) + 
  geom_line() + 
  facet_wrap(facets = ~ prior) + 
  theme(legend.position = 'none')

ggplot(rbind(loglik_lb_20chain.df, loglik_ch_20chain.df), 
       aes(iter, loglik, color = factor(chain_num))) + 
  geom_line() + 
  facet_wrap(facets = ~ prior) + 
  theme(legend.position = 'none')

ddf <- rbind(loglik_lb_20chain.df,loglik_ch_20chain.df)
ddf <- ddf[ddf$iter > 500, ]
ddf$prior <- paste0(ddf$prior, '~(m == 20)') 
library(latex2exp)
pdf('diabetes_summary_20.pdf', width = 10, height = 5)
ggplot(ddf[order(ddf$average_nn),], 
       aes(1 - accuracy, loglik, color = average_nn)) + 
  geom_point() + 
  facet_wrap(facets = ~ prior, labeller = label_parsed) +
  theme(panel.grid = element_line()) + 
  theme_classic() + 
  ylab('Log Likelihood') + 
  xlab('Missing Rate') + 
  labs(color = TeX('$\\bar{n}_L$')) + 
  scale_color_viridis_c()
dev.off()




ggplot(rbind(loglik_lb_20chain.df, loglik_ch_20chain.df)) + 
  geom_line(aes(iter, average_nn, linetype = factor(chain_num)), color = 'black') + 
  geom_line(aes(iter, q0.05_nn, linetype = factor(chain_num)), color = 'red') + 
  geom_line(aes(iter, q0.95_nn, linetype = factor(chain_num)), color = 'red') + 
  
  facet_wrap(facets = ~ prior) + 
  theme(legend.position = 'none')


nterm_sum_ch_20 <- aggregate(loglik_ch_20chain.df$average_nn, list(iter = loglik_ch_20chain.df$iter), mean)
nterm_sum_ch_20$q_05 <- aggregate(loglik_ch_20chain.df$q0.05_nn, list(iter = loglik_ch_20chain.df$iter), mean)[,2]
nterm_sum_ch_20$q_95 <- aggregate(loglik_ch_20chain.df$q0.95_nn, list(iter = loglik_ch_20chain.df$iter), mean)[,2]
nterm_sum_ch_20$q_75 <- aggregate(loglik_ch_20chain.df$q0.75_nn, list(iter = loglik_ch_20chain.df$iter), mean)[,2]
nterm_sum_ch_20$q_25 <- aggregate(loglik_ch_20chain.df$q0.25_nn, list(iter = loglik_ch_20chain.df$iter), mean)[,2]
nterm_sum_ch_20$q_50 <- aggregate(loglik_ch_20chain.df$q0.5_nn, list(iter = loglik_ch_20chain.df$iter), mean)[,2]

nterm_sum_ch_20$loglik <- aggregate(loglik_ch_20chain.df$loglik, list(iter = loglik_ch_20chain.df$iter), mean)[,2]
nterm_sum_ch_20$accuracy <- aggregate(loglik_ch_20chain.df$accuracy, list(iter = loglik_ch_20chain.df$iter), mean)[,2]


nterm_sum_ch_20$prior <- 'CL ~ (m == 20)'

nterm_sum_lb_20 <- aggregate(loglik_lb_20chain.df$average_nn, list(iter = loglik_lb_20chain.df$iter), mean)

nterm_sum_lb_20$q_05 <- aggregate(loglik_lb_20chain.df$q0.05_nn, list(iter = loglik_lb_20chain.df$iter), mean)[,2]
nterm_sum_lb_20$q_95 <- aggregate(loglik_lb_20chain.df$q0.95_nn, list(iter = loglik_lb_20chain.df$iter), mean)[,2]
nterm_sum_lb_20$q_75 <- aggregate(loglik_lb_20chain.df$q0.75_nn, list(iter = loglik_lb_20chain.df$iter), mean)[,2]
nterm_sum_lb_20$q_25 <- aggregate(loglik_lb_20chain.df$q0.25_nn, list(iter = loglik_lb_20chain.df$iter), mean)[,2]
nterm_sum_lb_20$q_50 <- aggregate(loglik_lb_20chain.df$q0.5_nn, list(iter = loglik_lb_20chain.df$iter), mean)[,2]

nterm_sum_lb_20$loglik <- aggregate(loglik_lb_20chain.df$loglik, list(iter = loglik_lb_20chain.df$iter), mean)[,2]
nterm_sum_lb_20$accuracy <- aggregate(loglik_lb_20chain.df$accuracy, list(iter = loglik_lb_20chain.df$iter), mean)[,2]

nterm_sum_lb_20$prior <- 'LB ~ (m == 20)'


ggplot(rbind(nterm_sum_lb_20, nterm_sum_ch_20)) + 
  geom_line(aes(iter, y = x), color = 'black') + 
  #geom_line(aes(iter, y = q_05), color = 'red') +
  geom_line(aes(iter, y = q_25), color = 'red') +
  geom_line(aes(iter, y = q_50), color = 'red') +
  geom_line(aes(iter, y = q_75), color = 'red') +
  geom_line(aes(iter, y = q_95), color = 'red') +
  facet_wrap(facets = ~ prior) + 
  theme_classic()



nterm_sum_20 <- rbind(nterm_sum_lb_20, nterm_sum_ch_20)
nterm_sum_20 <- nterm_sum_20[nterm_sum_20$iter > 500,]

ggplot(nterm_sum_20) + 
  geom_point(aes(x, loglik, color = prior)) +
  theme_classic() 

ggplot(nterm_sum_20) + 
  geom_point(aes(x, 1 - accuracy, color = prior)) +
  theme_classic() 


ggplot(nterm_sum_20) + 
  geom_point(aes(x, 1 - accuracy, color = iter)) +
  theme_classic() +
  facet_wrap(facets = ~ prior, ncol = 1) + 
  scale_color_viridis_c()

ggplot(nterm_sum_20) + 
  geom_point(aes(x, loglik, color = iter)) +
  theme_classic() +
  facet_wrap(facets = ~ prior, ncol = 1) + 
  scale_color_viridis_c()


ggplot(nterm_sum_20[order(nterm_sum_20$x),]) + 
  geom_point(aes(1 - accuracy, loglik, color = x)) +
  theme_classic() + 
  facet_wrap(facets = ~ prior) + 
  scale_color_viridis_c()

ggplot(nterm_sum_20[order(nterm_sum_20$x),]) + 
  geom_point(aes(iter, loglik, color = x)) +
  theme_classic() + 
  facet_wrap(facets = ~ prior) + 
  scale_color_viridis_c()

ggplot(nterm_sum_20[order(nterm_sum_20$x),]) + 
  geom_point(aes(iter, 1 - accuracy, color = x)) +
  theme_classic() + 
  facet_wrap(facets = ~ prior) + 
  scale_color_viridis_c()

ggplot(loglik_lb_20chain.df[loglik_lb_20chain.df$iter > 0,], aes(iter, loglik, color = chain)) + 
  geom_line()

ggplot(loglik_lb_20chain.df[loglik_lb_20chain.df$iter > 0,], aes(iter, 1 - accuracy, color = chain)) + 
  geom_line()





## 10 trees analysis




load('list_pred_lb_10_10chain.Rds')
loglik_lb_10 <- lapply(list_pred_lb_10_10chain, \(model_list) 
                            vapply(model_list, \(x) binary_log_lik_from_pred(x, yy), 0))


loglik_lb_10.todf <- lapply(1:length(loglik_lb_10), \(idxx) 
                                 data.frame(iter = 1:1500,
                                            loglik = loglik_lb_10[[idxx]],
                                            chain = paste0('LB - ', idxx),
                                            chain_num = idxx,
                                            prior = 'LB') 
                                 )
loglik_lb_10.df <- Reduce(rbind, loglik_lb_10.todf)
ggplot(loglik_lb_10.df, aes(iter, loglik, linetype = chain)) + 
  geom_line() + 
  theme(legend.position = 'none')




load('list_pred_ch_10_10chain.Rds')
loglik_ch_10 <- lapply(list_pred_ch_10_10chain, \(model_list) 
                            vapply(model_list, \(x) binary_log_lik_from_pred(x, yy), 0))


loglik_ch_10.todf <- lapply(1:length(loglik_ch_10), \(idxx) 
                                 data.frame(iter = 1:1500,
                                            loglik = loglik_ch_10[[idxx]],
                                            chain = paste0('CH - ', idxx),
                                            chain_num = idxx,
                                            prior = 'CH') 
)
loglik_ch_10.df <- Reduce(rbind, loglik_ch_10.todf)
ggplot(loglik_ch_10.df, aes(iter, loglik, linetype = chain)) + 
  geom_line() + 
  theme(legend.position = 'none')


ggplot(rbind(loglik_lb_10.df, loglik_ch_10.df),
       aes(iter, loglik, linetype = factor(chain_num))) +
  geom_line() + 
  theme(legend.position = 'none') + 
  facet_wrap(facets = ~ prior)


## accuracy
acc_lb_10 <- lapply(list_pred_lb_10_10chain, \(pred_list)
                    vapply(pred_list, \(x) accuracy_from_pred(x, yy), 0)
)

loglik_lb_10.df$acc = unlist(acc_lb_10)


acc_ch_10 <- lapply(list_pred_ch_10_10chain, \(pred_list)
                    vapply(pred_list, \(x) accuracy_from_pred(x, yy), 0)
)

loglik_ch_10.df$acc = unlist(acc_ch_10)




load('DIAB_list_models_lb_10.Rds')
load('DIAB_list_models_ch_10.Rds')
library(dplyr)
# number of terminal nodes
nterm_lb_10 <- lapply(1:length(list_models_lb), 
                   \(idxx) nterm_BART(list_models_lb[[idxx]]) %>% 
                     mutate(chain = idxx) )

nterm_lb_10_average <- lapply(nterm_lb_10, \(df) aggregate(df$nn, list(df$idx), mean)[,2] )

loglik_lb_10.df$average_nn = unlist(nterm_lb_10_average)


nterm_ch_10 <- lapply(1:length(list_models_ch), 
                   \(idxx) nterm_BART(list_models_ch[[idxx]]) %>% 
                     mutate(chain = idxx) )

nterm_ch_10_average <- lapply(nterm_ch_10, \(df) aggregate(df$nn, list(df$idx), mean)[,2] )

loglik_ch_10.df$average_nn = unlist(nterm_ch_10_average)




nterm_sum_lb_10 <- aggregate(loglik_lb_10.df$average_nn, list(iter = loglik_lb_10.df$iter), mean)
nterm_sum_lb_10$loglik <- aggregate(loglik_lb_10.df$loglik, list(iter = loglik_lb_10.df$iter), mean)[,2]
nterm_sum_lb_10$accuracy <- aggregate(loglik_lb_10.df$acc, list(iter = loglik_lb_10.df$iter), mean)[,2]
nterm_sum_lb_10$prior = 'LB ~ (m == 10)'


nterm_sum_ch_10 <- aggregate(loglik_ch_10.df$average_nn, list(iter = loglik_ch_10.df$iter), mean)
nterm_sum_ch_10$loglik <- aggregate(loglik_ch_10.df$loglik, list(iter = loglik_ch_10.df$iter), mean)[,2]
nterm_sum_ch_10$accuracy <- aggregate(loglik_ch_10.df$acc, list(iter = loglik_ch_10.df$iter), mean)[,2]
nterm_sum_ch_10$prior = 'CL ~ (m == 10)'

nterm_sum_10 <- rbind(nterm_sum_lb_10, nterm_sum_ch_10)
nterm_sum_10 <- nterm_sum_10[nterm_sum_10$iter > 500,]


nterm_sum_final = bind_rows(nterm_sum_10, nterm_sum_20)
nterm_sum_final$prior = factor(nterm_sum_final$prior, 
                               levels = c('LB ~ (m == 10)', 'LB ~ (m == 20)',
                                          'CL ~ (m == 10)', 'CL ~ (m == 20)')) 



pl_acc_10 <- ggplot(nterm_sum_10) + 
  geom_point(aes(x, 1 - accuracy, color = iter)) +
  theme_classic() +
  facet_wrap(facets = ~ prior, ncol = 1, labeller = label_parsed) + 
  scale_color_viridis_c() + 
  ylim(c(0.3760, 0.38)) +
  xlab(TeX('$\\bar{n}_L$')) + 
  ylab('Missing Rate') + 
  theme(legend.position = 'none')

pl_acc_20 <- ggplot(nterm_sum_20) + 
  geom_point(aes(x, 1 - accuracy, color = iter)) +
  theme_classic() +
  facet_wrap(facets = ~ prior, ncol = 1, labeller = label_parsed) + 
  scale_color_viridis_c() + 
  ylim(c(0.3760, 0.38)) +
  ylab('Missing Rate') +
  xlab(TeX('$\\bar{n}_L$')) + 
  theme(legend.position = 'none')

pdf('diab_final_res_miss.pdf', width = 10, height = 8)
gridExtra::grid.arrange(pl_acc_10, pl_acc_20, nrow = 1)
dev.off()

pl_loglik_10 <- ggplot(nterm_sum_10) + 
  geom_point(aes(x, loglik, color = iter)) +
  theme_classic() +
  facet_wrap(facets = ~ prior, ncol = 1, labeller = label_parsed) + 
  scale_color_viridis_c() + 
  ylim(c(-29830, -29690)) +
  xlab(TeX('$\\bar{n}_L$')) + 
  ylab('Log Likelihood') + 
  theme(legend.position = 'none')

pl_loglik_20 <- ggplot(nterm_sum_20) + 
  geom_point(aes(x, loglik, color = iter)) +
  theme_classic() +
  facet_wrap(facets = ~ prior, ncol = 1, labeller = label_parsed) + 
  scale_color_viridis_c() + 
  ylim(c(-29830, -29690)) +
  ylab('Log Likelihood') +
  xlab(TeX('$\\bar{n}_L$')) + 
  theme(legend.position = 'none')
pdf('diab_final_res_loglik.pdf', width = 10, height = 8)
gridExtra::grid.arrange(pl_loglik_10, pl_loglik_20, nrow = 1)
dev.off()


loglik_lb_20chain.df$prior <- 'CL'
df_20 <- rbind(loglik_lb_20chain.df, loglik_ch_20chain.df)
df_20$prior = paste0(df_20$prior, '~(m==20)')
pl_trace_20_loglik <- ggplot(df_20[df_20$iter > 500,], aes(iter, loglik, color = factor(chain_num))) +
  geom_line() + 
  facet_wrap(facets = ~prior, labeller = label_parsed) +
  theme_classic()  +
  ylab('Log Likelihood') + 
  xlab('Iteration') + 
  theme(legend.position = 'none') 

pl_trace_20_mr <- ggplot(df_20[df_20$iter > 500,], aes(iter, 1 - accuracy, color = factor(chain_num))) +
  geom_line() + 
  facet_wrap(facets = ~prior, labeller = label_parsed) +
  theme_classic()  +
  ylab('Missing Rate') + 
  xlab('Iteration') + 
  theme(legend.position = 'none') 

pl_trace_20_avnn <- ggplot(df_20[df_20$iter > 500,], aes(iter, average_nn, color = factor(chain_num))) +
  geom_line() + 
  facet_wrap(facets = ~prior, labeller = label_parsed) +
  theme_classic()  +
  ylab(TeX('$\\bar{n}_L')) + 
  xlab('Iteration') + 
  theme(legend.position = 'none') 

pdf('diab_traceplots_20.pdf', width = 10, height = 10) 
gridExtra::grid.arrange(pl_trace_20_loglik, pl_trace_20_mr, pl_trace_20_avnn)
dev.off()




loglik_ch_10.df$prior <- 'CL'
df_10 <- rbind(loglik_lb_10.df, loglik_ch_10.df)
df_10$prior = paste0(df_10$prior, '~(m==10)')
pl_trace_10_loglik <- ggplot(df_10[df_10$iter > 500,], aes(iter, loglik, color = factor(chain_num))) +
  geom_line() + 
  facet_wrap(facets = ~prior, labeller = label_parsed) +
  theme_classic()  +
  ylab('Log Likelihood') + 
  xlab('Iteration') + 
  theme(legend.position = 'none') 

pl_trace_10_mr <- ggplot(df_10[df_10$iter > 500,], aes(iter, 1 - acc, color = factor(chain_num))) +
  geom_line() + 
  facet_wrap(facets = ~prior, labeller = label_parsed) +
  theme_classic()  +
  ylab('Missing Rate') + 
  xlab('Iteration') + 
  theme(legend.position = 'none') 

pl_trace_10_avnn <- ggplot(df_10[df_10$iter > 500,], aes(iter, average_nn, color = factor(chain_num))) +
  geom_line() + 
  facet_wrap(facets = ~prior, labeller = label_parsed) +
  theme_classic()  +
  ylab(TeX('$\\bar{n}_L')) + 
  xlab('Iteration') + 
  theme(legend.position = 'none') 

pdf('diab_traceplots_10.pdf', width = 10, height = 10) 
gridExtra::grid.arrange(pl_trace_10_loglik, pl_trace_10_mr, pl_trace_10_avnn)
dev.off()







