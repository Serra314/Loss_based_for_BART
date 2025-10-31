source('code/import_functions.R')
library(data.tree)
library(ggplot2)
library(ggpubr)

# import data
bcancer.df <- read.table(file = 'data/breast_cancer_wisconsin_original/breast-cancer-wisconsin.data',
                         sep = ',')
# description
#1. Sample code number            id number
#2. Clump Thickness               1 - 10
#3. Uniformity of Cell Size       1 - 10
#4. Uniformity of Cell Shape      1 - 10
#5. Marginal Adhesion             1 - 10
#6. Single Epithelial Cell Size   1 - 10
#7. Bare Nuclei                   1 - 10
#8. Bland Chromatin               1 - 10
#9. Normal Nucleoli               1 - 10
#10. Mitoses                       1 - 10
#11. Class:                        (2 for benign, 4 for malignant)

colnames(bcancer.df) <- c('ID', 
                          'clump.thick', 
                          'cell.size', 
                          'cell.shape',
                          'marg.adhesion',
                          'se.cell.size',
                          'bar.nuclei',
                          'bland.chromatin',
                          'normal.nucleoli',
                          'mitosis',
                          'cancer.class')

head(bcancer.df)
# create predictor df 
X_pred <- bcancer.df[,c('clump.thick', 
                        'cell.size', 
                        'cell.shape',
                        'marg.adhesion',
                        'se.cell.size',
                        'bar.nuclei',
                        'bland.chromatin',
                        'normal.nucleoli',
                        'mitosis')]
summary(X_pred)
# trasnform in numeric bar.nuclei
X_pred$bar.nuclei <- as.numeric(X_pred$bar.nuclei)
# remove NA 
idx.na <- is.na(X_pred$bar.nuclei)
X_pred <- X_pred[!idx.na,]
rownames(X_pred) <- 1:nrow(X_pred)
# create variable 1 if malignant and 0 if benign
Y_mal <- bcancer.df$cancer.class == 4
Y_mal <- Y_mal[!idx.na]

## replicate chipman priors for the example and run the analysis.
chip.prior.small <- list(fun = chipman_prior_tree, param = c(0.95, 1.5))
chip.prior.medium <- list(fun = chipman_prior_tree, param = c(0.95, 1))
chip.prior.large <- list(fun = chipman_prior_tree, param = c(0.95, 0.5))

## mean, std, pr(nl > 8) 
chip.nl.mean.small <- 2.8
chip.nl.mean.medium <- 3.7
chip.nl.mean.large <- mean(rnterm.chip(10000, 0.95, 0.5))
# these values replicates the mean
omeg.chip.small <- -log(1 - 1/chip.nl.mean.small)
omeg.chip.medium <- -log(1 - 1/chip.nl.mean.medium)
omeg.chip.large <- -log(1 - 1/chip.nl.mean.large)

# tail prob
mean(rnterm.chip(100000, 0.95, 1) >= 13)
mean(rnterm(100000, omeg.chip.medium) >= 13)



# create lists representing different priors - 3 omega values x 2 beta values 
lb.prior.small.hb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.small, 1.5))
lb.prior.small.lb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.small, 0.5))
lb.prior.medium.hb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.medium, 1.5))
lb.prior.medium.lb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.medium, 0.5))
lb.prior.large.hb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.large, 1.5))
lb.prior.large.lb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.large, 0.5))


# normalise predictors 
X_pred.norm <- as.data.frame(apply(X_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
X_pred.norm <- as.matrix(X_pred.norm)
rownames(X_pred.norm) <- 1:nrow(X_pred.norm)

###########
# TABLE 2 #
###########

# calculate correlations
cor(X_pred.norm, Y_mal)




n.chain_par <- 100
n.iter_par <- 500
incl.split_par <- FALSE
cont.unif_par <- TRUE
moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)

#############
## DEFUALT ##
#############

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.56, 0.62))
#debugonce(multichain_MCMC_binary)

chip.prior.def <- list(fun = chipman_prior_tree, param = c(0.95, 2))
mcmc_chip.def <- multichain_MCMC_binary(n.iter = n.iter_par,
                                        n.chain = n.chain_par,
                                        X = X_pred.norm,
                                        Y = Y_mal,
                                        alpha.prior = 1,
                                        beta.prior = 1,
                                        prior_list = chip.prior.def,
                                        include.split = incl.split_par,
                                        cont.unif = cont.unif_par,
                                        moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_chip.def, file = 'code/results/breastcancer.chip.def.Rds')
rm(mcmc_chip.def)

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.56, 0.62))
mcmc_lb.def <- multichain_MCMC_binary(n.iter = n.iter_par,
                                        n.chain = n.chain_par,
                                        X = X_pred.norm,
                                        Y = Y_mal,
                                        alpha.prior = 1,
                                        beta.prior = 1,
                                        prior_list = lb.prior.def,
                                      include.split = incl.split_par,
                                      cont.unif = cont.unif_par,
                                        moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_lb.def, file = 'code/results/breastcancer.lb.def.Rds')
rm(mcmc_lb.def)

###########
## SMALL ##
###########
# - CHIPMAN
mcmc_chip.small <- multichain_MCMC_binary(n.iter = n.iter_par,
                                            n.chain = n.chain_par,
                                            X = X_pred.norm,
                                            Y = Y_mal,
                                            alpha.prior = 1,
                                            beta.prior = 1,
                                            prior_list = chip.prior.small,
                                          include.split = incl.split_par,
                                          cont.unif = cont.unif_par,
                                            moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_chip.small, file = 'code/results/breastcancer.chip.small.Rds')
rm(mcmc_chip.small)
# - Lossbased - beta = 1.5
mcmc_lb.small.hb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                          n.chain = n.chain_par,
                                          X = X_pred.norm,
                                          Y = Y_mal,
                                          alpha.prior = 1,
                                          beta.prior = 1,
                                          prior_list = lb.prior.small.hb,
                                          include.split = incl.split_par,
                                          cont.unif = cont.unif_par,
                                          moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_lb.small.hb, file = 'code/results/breastcancer.lb.small.hb.Rds')
rm(mcmc_lb.small.hb)
# - Lossbased - beta = 0.5
mcmc_lb.small.lb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                           n.chain = n.chain_par,
                                           X = X_pred.norm,
                                           Y = Y_mal,
                                           alpha.prior = 1,
                                           beta.prior = 1,
                                           prior_list = lb.prior.small.lb,
                                           include.split = incl.split_par,
                                           cont.unif = cont.unif_par,
                                           moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_lb.small.lb, file = 'code/results/breastcancer.lb.small.lb.Rds')
rm(mcmc_lb.small.lb)

############
## MEDIUM ##
############
# - CHIPMAN
mcmc_chip.medium <- multichain_MCMC_binary(n.iter = n.iter_par,
                                          n.chain = n.chain_par,
                                          X = X_pred.norm,
                                          Y = Y_mal,
                                          alpha.prior = 1,
                                          beta.prior = 1,
                                          prior_list = chip.prior.medium,
                                          include.split = incl.split_par,
                                          cont.unif = cont.unif_par,
                                          moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_chip.medium, file = 'code/results/breastcancer.chip.medium.Rds')
rm(mcmc_chip.medium)
# - Lossbased - beta = 1.5
mcmc_lb.medium.hb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                           n.chain = n.chain_par,
                                           X = X_pred.norm,
                                           Y = Y_mal,
                                           alpha.prior = 1,
                                           beta.prior = 1,
                                           prior_list = lb.prior.medium.hb,
                                           include.split = incl.split_par,
                                           cont.unif = cont.unif_par,
                                           moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_lb.medium.hb, file = 'code/results/breastcancer.lb.medium.hb.Rds')
rm(mcmc_lb.medium.hb)
# - Lossbased - beta = 0.5
mcmc_lb.medium.lb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                           n.chain = n.chain_par,
                                           X = X_pred.norm,
                                           Y = Y_mal,
                                           alpha.prior = 1,
                                           beta.prior = 1,
                                           include.split = incl.split_par,
                                           cont.unif = cont.unif_par,
                                           prior_list = lb.prior.medium.lb,
                                           moves.prob = c(0.4, 0.4, 0.1, 0.1))
save(mcmc_lb.medium.lb, file = 'code/results/breastcancer.lb.medium.lb.Rds')
rm(mcmc_lb.medium.lb)

####################
## DEFAULT MODELS ##
####################
load('code/results/breastcancer.chip.def.Rds')
load('code/results/breastcancer.lb.def.Rds')


model.list.def <- list(#mcmc_chip.medium,
                   #mcmc_chip.small,
                   mcmc_chip.def,
                   mcmc_lb.def)#,
                   #mcmc_lb.medium.lb,
                   #mcmc_lb.medium.hb,
                   #mcmc_lb.small.lb,
                   #mcmc_lb.small.hb)
names(model.list.def) <- c(#'CH - a = 0.95, b = 1',
                       #'CH - a = 0.95, b = 1.5',
                       'CL - default',
                       'LB - default')#,
                       #'LB - o = 0.30, g = 0.5',
                       #'LB - o = 0.30, g = 1.5',
                       #'LB - o = 0.42, g = 0.5',
                       #'LB - o = 0.42, g = 1.5')


# extract depth, number of terminal nodes, missing rate and loglik of all the trees
depth.df <- apply_fun_models(fun_ = get_depth, 
                             mcmc.list = model.list.def,
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)
nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes, 
                             mcmc.list = model.list.def, 
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)

rapid.miss.rate <- function(tree_top, X, Y){
  pp <- get_value_tree(tree_top, X)
  pred_ <- pp >= 0.5
  sum(pred_ != Y)
}


miss.rate.df <- apply_fun_models(fun_ = \(x) rapid.miss.rate(tree_top = x, X = X_pred.norm, Y = Y_mal),
                                 mcmc.list = model.list.def,
                                 born.out.pc = 250, 
                                 n.chain = n.chain_par, 
                                 sample.pc = n.iter_par)



# miss.rate.df <- apply_fun_models_all.tog(mcmc.list = model.list.def,
#                                          fun_ = \(model.out) miss.class.rate.mcmc(tree_list = model.out$trees[1:10],
#                                                                                   X = X_pred.norm, 
#                                                                                   Y = Y_mal, 
#                                                                                   born.out.pc = 10, 
#                                                                                   n.chain = n.chain_par, 
#                                                                                   sample.pc = n.iter_par))

like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_binary(tree_top = x, Y = Y_mal, 
                                                            X = X_pred.norm), 
                            mcmc.list = model.list.def, 
                            born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)




df.sum.def <- data.frame(tree = nterm.df$x,
                     panel.name = nterm.df$panel.name,
                     loglik = like.df$y,
                     nterm = nterm.df$y,
                     miss.rate = miss.rate.df$y)


df.sum.def$nterm.cut <- cut(df.sum.def$nterm, breaks = c(3,7,10,13,25))

save(df.sum.def, file = 'code/results/breastcancer.df.summary_def.Rds')
load('code/results/breastcancer.df.summary_def.Rds')

df.loglik_complete <- apply_fun_models(\(x) cart_log_lik_binary(x, Y_mal, X_pred.norm), 
                                       model.list.def, born.out.pc = 0,
                                       n.chain = n.chain_par, sample.pc = n.iter_par)
save(df.loglik_complete, file = 'breast_cancer_loglik_complete.Rds')

nterm.df_complete <- apply_fun_models(fun_ = get_num_terminal_nodes, 
                             mcmc.list = model.list.def, 
                             born.out.pc = 0, n.chain = n.chain_par, sample.pc = n.iter_par)

df.prior_fit.lb.def <- lapply(idx_list, 
                              \(x) joint.prior.new.tree(mcmc_lb.def$trees[[x]], 1.56, 0.62))
df.prior_fit.lb.def <- data.frame(prior = do.call(c, df.prior_fit.lb.def),
                                  panel.name = 'LB - default')

df.prior_fit.chip.def <- lapply(mcmc_chip.def$trees[idx_list], \(x) chipman_prior_tree(x, 0.95, 2))
df.prior_fit.chip.def <- data.frame(prior = do.call(c, df.prior_fit.chip.def),
                                    panel.name = 'CH - default')

prior.vec <- c(df.prior_fit.chip.def$prior, df.prior_fit.lb.def$prior)
df.loglik_complete$prior = prior.vec
df.loglik_complete$posterior = df.loglik_complete$y + log(df.sum.def$prior)


####

df.loglik_complete$nterm <- nterm.df_complete$y


library(dplyr)
pl1 <- df.loglik_complete %>%
  group_by(nterm, panel.name) %>%
  summarise(
    mean_y = mean(y),
    q025 = quantile(y, 0.025),
    q975 = quantile(y, 0.975),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = nterm, y = mean_y, color = panel.name)) +
  geom_point() +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.2) +
  theme_minimal() +
  xlim(4, 16) + 
  ylim(-140, -40)

pl2 <- df.sum.def %>%
  group_by(nterm, panel.name) %>%
  summarise(
    mean_y = mean(loglik),
    q025 = quantile(loglik, 0.025),
    q975 = quantile(loglik, 0.975),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = nterm, y = mean_y, color = panel.name)) +
  geom_point() +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.2) +
  theme_minimal() +
  xlim(4, 16) + 
  ylim(-140, -40)

library(gridExtra)

grid.arrange(pl1, pl2, ncol = 2)



pl3 <- df.sum.def %>%
  group_by(nterm, panel.name) %>%
  summarise(
    mean_y = mean(loglik),
    q025 = quantile(loglik, 0.025),
    q975 = quantile(loglik, 0.975),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = nterm, y = mean_y, color = panel.name)) +
  geom_point() +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.2) +
  theme_minimal() + 
  xlim(3,17)

pl5 <- df.sum.def %>%
  group_by(nterm, panel.name) %>%
  summarise(
    mean_y = mean(loglik),
    q025 = quantile(loglik, 0.025),
    q975 = quantile(loglik, 0.975),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = nterm, y = q975 - q025, color = panel.name)) +
  geom_point() +
  theme_minimal() + 
  xlim(3,17)

pl4 <- df.sum.def %>%
  count(nterm, panel.name) %>%
  ggplot(aes(x = nterm, y = n, fill = panel.name)) +
  geom_col(position = "dodge") +
  theme_minimal()+
  xlim(3,17)
grid.arrange(pl3, pl5, pl4)

df.sum.def %>%
  group_by(nterm, panel.name) %>%
  summarise(
    n = n(),
    iqr_y = IQR(loglik),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = iqr_y, color = panel.name)) +
  geom_point() +
  theme_minimal()

# CL prior
df.sum.def.red <- df.loglik_complete[df.loglik_complete$panel.name == 'CL - default', ]
nterm.df_complete.red <- nterm.df_complete[nterm.df_complete$panel.name == 'CL - default',]

df.sum.def.red[which.max(df.sum.def.red$posterior),]
get_tree_plot(mcmc_chip.def$trees[[25461]])
nterm.df_complete.red[25461,]

df.sum.def.red[which.max(df.sum.def.red$y),]
get_tree_plot(mcmc_chip.def$trees[[10483]])
nterm.df_complete.red[10483,]


df.sum.def.red <- df.loglik_complete[df.loglik_complete$panel.name == 'LB - default', ]
nterm.df_complete.red <- nterm.df_complete[nterm.df_complete$panel.name == 'LB - default',]

df.sum.def.red[which.max(df.sum.def.red$posterior),]
get_tree_plot(mcmc_lb.def$trees[[43248]])
nterm.df_complete.red[43248,]

df.sum.def.red[which.max(df.sum.def.red$y),]
get_tree_plot(mcmc_lb.def$trees[[43362]])
nterm.df_complete.red[43362,]

df.loglik_complete[order(df.loglik_complete$y, decreasing = TRUE),]

# create data.frame for prior

nterm.chip.ecdf.def <- ecdf(rnterm.chip(10000, 0.95, 2))
df.nl.prior <- rbind(data.frame(x = 1:20, 
                                y = prior.nterm(1:20, 1.56),
                                panel.name = 'LB - default'),
                     data.frame(x = 1:20, 
                                y = diff(nterm.chip.ecdf.def(0:20)),
                                panel.name = 'CL - default'))

# create df for post quantities
df.split <- split(df.sum.def, df.sum.def$panel.name)
post.nl.mean <- vapply(df.split, \(x) mean(x$nterm), 0)
post.nl.q0.025 <- vapply(df.split, \(x) quantile(x$nterm, 0.025), 0)
post.nl.q0.975 <- vapply(df.split, \(x) quantile(x$nterm, 0.975), 0)
df.post.nl.sum <- data.frame(mean = post.nl.mean, 
                             low = post.nl.q0.025, 
                             up = post.nl.q0.975,
                             panel.name = names(df.split))




pdf(file = 'draft/figures/fig_posterior_nl_dist_and_trace_bc.pdf', width = 6, height = 3)
pl.nl.post <- ggplot(df.sum.def, aes(x = nterm, y = after_stat(density))) + 
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = seq(1,21,by = 2), color = 'grey', alpha = 0.3) + 
  geom_line(data = df.nl.prior, aes(x,y), colour = 'deepskyblue', linewidth = 0.4) +
  geom_point(data = df.nl.prior, aes(x,y), colour = 'deepskyblue', size = 0.4) + 
  geom_vline(data = df.post.nl.sum, aes(xintercept = mean)) +
  geom_vline(data = df.post.nl.sum, aes(xintercept = low), linetype = 2) +
  geom_vline(data = df.post.nl.sum, aes(xintercept = up), linetype = 2) +
  facet_wrap(facets = ~panel.name, ncol = 1) + 
  theme_classic() + 
  xlab('Number of terminal nodes') + 
  scale_x_continuous(breaks = seq(1,21,by = 2))
pl.nl.post
dev.off()

pdf(file = 'draft/figures/fig_loglik_trace_bc.pdf', width = 6, height = 3)
pl.loglik.trace <- ggplot(df.sum.def, aes(tree, loglik)) + 
  geom_line() + 
  facet_wrap(facets = ~panel.name) + 
  theme_classic() + 
  xlab('Iteration') + 
  ylab('Log likelihood')
pl.loglik.trace
dev.off()

#############
# FIGURE 9 #
############
pdf(file = 'draft/figures/fig_posterior_nl_dist_and_trace_bc.pdf', width = 6, height = 6)
ggarrange(plotlist = list(pl.nl.post, pl.loglik.trace), ncol = 1)
dev.off()



# set bornout and select iterations
remove.bornout <- function(df.sum, born.out, sample.per.chain, n.chain){
  df.split <- split(df.sum, df.sum$panel.name)
  idx.bornout.sm <- sort(unlist(outer(1:born.out, sample.per.chain*(0:(n.chain - 1)), '+')))
  df.split.red <- lapply(df.split, function(x) x[-idx.bornout.sm,] %>% 
                           mutate(new.idx = 1:(nrow(x) - length(idx.bornout.sm))))
  Reduce(rbind, df.split.red)
}


#############
# FIGURE 10 #
############
pdf(file = 'draft/figures/fig_posterior_nl_bc.pdf', width = 6, height = 3)
ggplot(df.sum.def, aes(x = nterm, y = loglik)) + 
  #geom_vline(xintercept = , color = 'grey', alpha = 0.3) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  theme_classic() + 
  ylab('log Likelihood') + 
  xlab('Number of terminal nodes')
dev.off()

#############
# FIGURE 11 #
############
pdf(file = 'draft/figures/fig_posterior_miss_bc.pdf', width = 6, height = 3)
ggplot(df.sum.def, aes(x = nterm, y = miss.rate)) + 
  #geom_vline(xintercept = , color = 'grey', alpha = 0.3) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  theme_classic() + 
  ylab('Missing rate') + 
  xlab('Number of terminal nodes')
dev.off()


#############
# FIGURE 12 #
############
pdf(file = 'draft/figures/fig_posterior_log_miss_bc.pdf', width = 6, height = 3)
ggplot(df.sum.def[order(df.sum.def$nterm),], 
       aes(x = miss.rate, y = loglik, color = nterm.cut)) + 
  #geom_vline(xintercept = , color = 'grey', alpha = 0.3) + 
  geom_point() + 
  labs(color = ~n[L]) + 
  scale_color_viridis_d() + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  theme_classic() + 
  ylab('Log Likelihood') + 
  xlab('Missing Rate')
dev.off()


ggplot(df.sum.def, aes(x = loglik, #y = after_stat(density), 
                       color = panel.name, fill = panel.name)) + 
  geom_histogram(bins = 50, alpha = 0.5) 


# best miss rate and loglik
vapply(df.split, \(x) min(x$miss.rate), 0)
vapply(df.split, \(x) max(x$loglik), 0)

#########################
## ALTERNATIVE MODELS ###
######################### 


load('code/results/breastcancer.chip.medium.Rds')
load('code/results/breastcancer.chip.small.Rds')
load('code/results/breastcancer.lb.medium.hb.Rds')
load('code/results/breastcancer.lb.medium.lb.Rds')
load('code/results/breastcancer.lb.small.hb.Rds')
load('code/results/breastcancer.lb.small.lb.Rds')

#################
## PREDICTORS ###
#################


get_pred <- function(tree_top, predictors = c()){
  if(is.null(tree_top$right) & is.null(tree_top$left)){
    return(predictors)
  } else{
    predictors <- unique(c(predictors, tree_top$cond$x.idx))
    pred.left <- get_pred(tree_top$left, predictors)
    pred.right <- get_pred(tree_top$right, predictors)
    return(unique(c(predictors, pred.left, pred.right)))
  }
}

count.pred <- function(tree.pred, poss.pred){
  vapply(poss.pred, \(x) sum(tree.pred == x), 0)
}
# frequency for each predictor
get_pred_freq <- function(tree_list, poss.pred){
  pred.list <- lapply(tree_list, get_pred)
  pred.freq.list <- lapply(pred.list, \(x) count.pred(x, poss.pred))
  pred.freq.matrix <- Reduce(rbind, pred.freq.list)
  output = colSums(pred.freq.matrix)/nrow(pred.freq.matrix)
  names(output) <- paste0('x', poss.pred)
  output
}

burnin = 250
idx_burnin <- rep(c(rep(FALSE, burnin), rep(TRUE, n.iter_par - burnin)), n.chain_par)

pred_lb_def = get_pred_freq(mcmc_lb.def$trees[idx_burnin], 1:9)
pred_chip_def = get_pred_freq(mcmc_chip.def$trees[idx_burnin], 1:9)


pred_def <- rbind(data.frame(Variable = colnames(X_pred), 
                       Value = pred_lb_def,
                       Model = 'Lb - def'), 
                  data.frame(Variable = colnames(X_pred), 
                             Value = pred_chip_def,
                             Model = 'CL- def'))
                        
                       
model_labels <- c( expression(atop("CL", alpha == 0.95 * "," ~ beta == 2)),
                   expression(atop("LB", omega == 1.56 * "," ~ gamma == 0.62)))
#############
# Figure 13 #
#############
pdf(file = 'draft/figures/pred_synth_def.pdf', width = 6, height = 3)
ggplot(pred_def, aes(x = Variable, y = Model, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 2)), color = "black", size = 4) +
  scale_fill_gradient(low = "#fef0d9", high = "#b30000") +
  scale_y_discrete(labels = model_labels) + 
  labs(title = "", x = "", y = "", fill = 'Freq') +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.6, size = 7.5),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")  # ✅ removes all outer white margins
  )
dev.off()



# predictors for alternative models

pred_chip_medium = get_pred_freq(mcmc_chip.medium$trees[idx_burnin], 1:9)
pred_lb_medium.hb = get_pred_freq(mcmc_lb.medium.hb$trees[idx_burnin], 1:9)
pred_lb_medium.lb = get_pred_freq(mcmc_lb.medium.lb$trees[idx_burnin], 1:9)

pred_chip_small = get_pred_freq(mcmc_chip.small$trees[idx_burnin], 1:9)
pred_lb_small.hb = get_pred_freq(mcmc_lb.small.hb$trees[idx_burnin], 1:9)
pred_lb_small.lb = get_pred_freq(mcmc_lb.small.lb$trees[idx_burnin], 1:9)


cbind(colnames(X_pred), pred_chip_def, pred_lb_def)
cbind(colnames(X_pred), pred_chip_medium, pred_lb_medium.hb, pred_lb_medium.lb)
cbind(colnames(X_pred), pred_chip_small, pred_lb_small.hb, pred_lb_small.lb)



model.list.alt <- list(mcmc_chip.medium,
                       mcmc_chip.small,
                       mcmc_lb.medium.lb,
                       mcmc_lb.medium.hb,
                       mcmc_lb.small.lb,
                       mcmc_lb.small.hb)
names(model.list.alt) <- c('CH - a = 0.95, b = 1',
                           'CH - a = 0.95, b = 1.5',
                           'LB - o = 0.30, g = 0.5',
                           'LB - o = 0.30, g = 1.5',
                           'LB - o = 0.42, g = 0.5',
                           'LB - o = 0.42, g = 1.5')



ggplot(pred_df, aes(x = Model, y = Variable, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 2)), color = "black", size = 4) +
  scale_fill_gradient(low = "#fef0d9", high = "#b30000") +
  scale_x_discrete(labels = model_labels) + 
  labs(title = "", x = "", y = "", fill = 'Freq') +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 7.5),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")  # ✅ removes all outer white margins
  )




# extract depth, number of terminal nodes, missing rate and loglik of all the trees
depth.df.alt <- apply_fun_models(fun_ = get_depth, 
                             mcmc.list = model.list.alt,
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)
nterm.df.alt <- apply_fun_models(fun_ = get_num_terminal_nodes, 
                             mcmc.list = model.list.alt, 
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)

rapid.miss.rate <- function(tree_top, X, Y){
  pp <- get_value_tree(tree_top, X)
  pred_ <- pp >= 0.5
  sum(pred_ != Y)
}

miss.rate.df.alt <- apply_fun_models(fun_ = \(x) rapid.miss.rate(tree_top = x, X = X_pred.norm, Y = Y_mal),
                                 mcmc.list = model.list.alt,
                                 born.out.pc = 250, 
                                 n.chain = n.chain_par, 
                                 sample.pc = n.iter_par)


like.df.alt <- apply_fun_models(fun_ = \(x) cart_log_lik_binary(tree_top = x, Y = Y_mal, 
                                                            X = X_pred.norm), 
                            mcmc.list = model.list.alt, 
                            born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)




df.sum.alt <- data.frame(tree = nterm.df.alt$x,
                         panel.name = nterm.df.alt$panel.name,
                         loglik = like.df.alt$y,
                         nterm = nterm.df.alt$y,
                         miss.rate = miss.rate.df.alt$y)


df.sum.alt$nterm.cut <- cut(df.sum.alt$nterm, breaks = c(3,7,10,13,25))

save(df.sum.alt, file = 'code/results/breastcancer.df.summary_alt.Rds')
load('code/results/breastcancer.df.summary_alt.Rds')


nterm.chip.ecdf.med <- ecdf(rnterm.chip(10000, 0.95, 1))
nterm.chip.ecdf.small <- ecdf(rnterm.chip(10000, 0.95, 1.5))


df.nl.prior.alt <- rbind(data.frame(x = 1:20, 
                                   y = diff(nterm.chip.ecdf.med(0:20)),
                                   panel.name = 'CH - a = 0.95, b = 1'),
                         data.frame(x = 1:20, 
                                   y = diff(nterm.chip.ecdf.small(0:20)),
                                   panel.name = 'CH - a = 0.95, b = 1.5'),
                         data.frame(x = 1:20, 
                                y = prior.nterm(1:20, 0.3),
                                panel.name = 'LB - o = 0.30, g = 0.5'),
                         data.frame(x = 1:20, 
                                    y = prior.nterm(1:20, 0.42),
                                    panel.name = 'LB - o = 0.42, g = 0.5'),
                            data.frame(x = 1:20, 
                                    y = prior.nterm(1:20, 0.3),
                                    panel.name = 'LB - o = 0.30, g = 1.5'),
                         data.frame(x = 1:20, 
                                    y = prior.nterm(1:20, 0.42),
                                    panel.name = 'LB - o = 0.42, g = 1.5'))

# create df for post quantities
df.split.alt <- split(df.sum.alt, df.sum.alt$panel.name)
post.nl.mean.alt <- vapply(df.split.alt, \(x) mean(x$nterm), 0)
post.nl.q0.025.alt <- vapply(df.split.alt, \(x) quantile(x$nterm, 0.025), 0)
post.nl.q0.975.alt <- vapply(df.split.alt, \(x) quantile(x$nterm, 0.975), 0)
df.post.nl.sum.alt <- data.frame(mean = post.nl.mean.alt, 
                             low = post.nl.q0.025.alt, 
                             up = post.nl.q0.975.alt,
                             panel.name = names(df.split.alt))


factor.levs <- c( 'CH - a = 0.95, b = 1',
                  'CH - a = 0.95, b = 1.5',
                  'LB - o = 0.30, g = 0.5',
                  'LB - o = 0.42, g = 0.5',
                  'LB - o = 0.30, g = 1.5',
                  'LB - o = 0.42, g = 1.5')


factor.labs <- c( 'CH - a = 0.95, b = 1' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 1$')),
                  'CH - a = 0.95, b = 1.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 1.5$')),
                  'LB - o = 0.30, g = 0.5' = parse(text=latex2exp::TeX('$LB - \\omega = 0.3, \\gamma = 0.5$')),
                  'LB - o = 0.42, g = 0.5' = parse(text=latex2exp::TeX('$LB - \\omega = 0.42, \\gamma = 0.5$')),
                  'LB - o = 0.30, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.3, \\gamma = 1.5$')),
                  'LB - o = 0.42, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.42, \\gamma = 1.5$')))
 

df.sum.alt$panel.name <- factor(df.sum.alt$panel.name,
                                levels = factor.levs,
                            labels = factor.labs,
                            ordered = TRUE)
df.nl.prior.alt$panel.name <- factor(df.nl.prior.alt$panel.name,
                                     levels = factor.levs,
                                 labels = factor.labs,
                                 ordered = TRUE)
df.post.nl.sum.alt$panel.name <- factor(df.post.nl.sum.alt$panel.name,
                                        levels = factor.levs,
                                     labels = factor.labs,
                                     ordered = TRUE)


###############################################
# Supplementary material Appendix H FIGURE 7 #
##############################################

pdf(file = 'draft/figures/fig_posterior_nl_dist_bc_nodef.pdf', width = 7, height = 8)
pl.nl.post.alt <- ggplot(df.sum.alt, aes(x = nterm, y = after_stat(density))) + 
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = seq(1,21,by = 2), color = 'grey', alpha = 0.3) + 
  geom_line(data = df.nl.prior.alt, aes(x,y), colour = 'deepskyblue', linewidth = 0.4) +
  geom_point(data = df.nl.prior.alt, aes(x,y), colour = 'deepskyblue', size = 0.4) + 
  geom_vline(data = df.post.nl.sum.alt, aes(xintercept = mean)) +
  geom_vline(data = df.post.nl.sum.alt, aes(xintercept = low), linetype = 2) +
  geom_vline(data = df.post.nl.sum.alt, aes(xintercept = up), linetype = 2) +
   facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Number of terminal nodes') + 
  scale_x_continuous(breaks = seq(1,21,by = 2))
pl.nl.post.alt
dev.off()

###############################################
# Supplementary material Appendix H FIGURE 8 #
##############################################

pdf(file = 'draft/figures/fig_posterior_nl_bc_nodef.pdf', width = 7, height = 8)
ggplot(df.sum.alt, aes(x = nterm, y = loglik)) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Number of terminal nodes') +
  ylab('Log likelihood') 
dev.off()


###############################################
# Supplementary material Appendix H FIGURE 9 #
##############################################

pdf(file = 'draft/figures/fig_posterior_miss_bc_nodef.pdf', width = 7, height = 8)
ggplot(df.sum.alt, aes(x = nterm, y = miss.rate)) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Number of terminal nodes') +
  ylab('Missing rate') 
dev.off()


###############################################
# Supplementary material Appendix H FIGURE 10 #
##############################################

pdf(file = 'draft/figures/fig_posterior_log_miss_bc_nodef.pdf', width = 7, height = 8)
ggplot(df.sum.alt[order(df.sum.alt$nterm),], aes(x = miss.rate, y = loglik, color = nterm.cut)) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Missing rate') +
  ylab('Log likelihood') + 
  scale_color_viridis_d() + 
  labs(color = ~n[L])
dev.off()









