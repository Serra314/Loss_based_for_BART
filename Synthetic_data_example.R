source('code/import_functions.R')
library(data.tree)
library(ggplot2)
# needed to combine ggplot obj
library(ggpubr)
# for plotting labels on contour lines
library(metR)
# to write in latex in ggplot
library(latex2exp)
# to create df for plotting
library(foreach)
# just for combining plots
library(inlabru)
library(manipulateWidget)
library(plotly)
library(dplyr)
library(patchwork)


#######################
# generate predictors #
#######################

set.seed(1)
data.mult <- 1
x1 <- c(runif(200*data.mult, 0.1, 0.4), runif(100*data.mult, 0.6, 0.9))
x2 <- c(runif(100*data.mult, 0.1, 0.4), runif(100*data.mult, 0.6, 0.9), runif(100*data.mult, 0.1, 0.9))
x3 <- c(runif(200*data.mult, 0.6, 0.9), runif(100*data.mult, 0.1, 0.4))
# store them in data.frame
X = data.frame(x1 = x1, x2 = x2, x3 = x3)
# normalise data
#X <- as.data.frame(apply(X, 2, \(x) (x - min(x))/(max(x) - min(x))))
X <- as.matrix(X)
rownames(X) <- 1:nrow(X)

###############
# create tree #
###############

# initialize list
tree_ex <- list()
# one split on the left branch
tree_ex$left <- list(left = list(left = NULL,
                                 right = NULL),
                     right = list(left = NULL,
                                  right = NULL))
# no split on the right branch
tree_ex$right <-  list(left = NULL,
                       right = NULL)

# assign node index
tree_ex <- assign_node_idx(tree_ex)
# plot the tree
get_tree_plot.idx(tree_ex)

# assign first splitting rule
tree_ex$cond <- list(x.idx = 1, x.val = 0.5)
# assign second splitting rule
tree_ex$left$cond <- list(x.idx = 2, x.val = 0.5)
# set mu_3
tree_ex$left$left$cond <- NULL
tree_ex$left$left$value <- 1
# set mu_4
tree_ex$left$right$cond <- NULL
tree_ex$left$right$value <- 3
# set mu_5
tree_ex$right$cond <- NULL
tree_ex$right$value <- 5
# plot the tree
tree.plot <- get_tree_plot(tree_ex)

###########################
# FIGURE 6 TOP-LEFT PANEL #
###########################
tree.plot

#########################
# generate observations #
#########################

Y.var = 0.25
# sample observations
Y_ex <- sample_CART(tree_top = tree_ex, X = X, sigma_ = Y.var)
# store in data.frame
data.df <- cbind(Y = Y_ex, X)
# plot
pl.x1 <- ggplot(data.df, aes(x1, Y)) + geom_point(size = 0.5)
pl.x2 <- ggplot(data.df) + geom_point(aes(x2, Y), size = 0.5) + 
  xlab(~x[2])
pl.x3 <- ggplot(data.df) + geom_point(aes(x3, Y), size = 0.5) + 
  xlab(~x[3])


pl.x1.g <- ggplotly(p = pl.x1 + theme_classic() + ylab(TeX('Y')) + xlab(TeX('X_1')) +
                theme(axis.text.x = element_text(angle = 0, hjust = 1))) %>% 
  config(mathjax = "cdn")
pl.x2.g <- ggplotly(p = pl.x2 + theme_classic() + ylab(TeX('Y')) + xlab(TeX('X_2')) +
                      theme(axis.text.x = element_text(angle = 0, hjust = 1))) %>% 
  config(mathjax = "cdn")
pl.x3.g <- ggplotly(p = pl.x3 + theme_classic() + ylab(TeX('Y')) + xlab(TeX('X_3')) +
                      theme(axis.text.x = element_text(angle = 0, hjust = 1))) %>% 
  config(mathjax = "cdn")

#############
# FIGURE 6 #
############
combineWidgets(tree.plot, pl.x1.g, pl.x2.g, pl.x3.g)

## expected number of terminal nodes and depth for different priors - for TABLE 1
# chip default
mean(rdepth.chip(10000, 0.95, 2))
mean(rnterm.chip(10000, 0.95, 2))

# lower alpha
mean(rdepth.chip(10000, 0.25, 2))
mean(rnterm.chip(10000, 0.25, 2))

# lower beta
mean(rdepth.chip(10000, 0.25, 0.5))
mean(rnterm.chip(10000, 0.25, 0.5))

# loss-based default
mean(rdepth(10000, 1.56, 0.62))
mean(rnterm(10000, 1.56))

# lower omega
mean(rdepth(10000, 0.5, 0.62))
mean(rnterm(10000, 0.5))

# higher gamma
mean(rdepth(10000, 0.5, 1.5))
mean(rnterm(10000, 0.5))


##################################
## MULTIPLE CHAINS parallel MCMC #
##################################

# value of the parameters max expected loss
#pp.2 <- c(1.5618883, 0.6293944)
pp.2 <- c(1.56, 0.62)

# loss-based prior with default params
lossb.prior.list <- list(fun = joint.prior.new.tree, param = pp.2)
n.chain_par <- 500
n.iter_par <- 500
incl.split_par <- FALSE
cont.unif_par <- TRUE
moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)
n.cores = 15

#set seed for reproducibility
set.seed(1)

st = Sys.time()
MCMC_lossb.prior.multi <- multichain_MCMC_known_var_windows(n.chain = n.chain_par,
                                                  n.iter = n.iter_par, 
                                                  X = X, 
                                                  Y = Y_ex, 
                                                  Y.var = Y.var, 
                                                  mu.prior.mean = 0, 
                                                  mu.prior.var = var(Y_ex), 
                                                  prior_list = lossb.prior.list, 
                                                  moves.prob = NULL, starting.tree = NULL,
                                                  include.split = incl.split_par,
                                                  cont.unif = cont.unif_par, 
                                                  n.cores = n.cores)
Sys.time() -st
#  1.23 hours
save(MCMC_lossb.prior.multi, file = 'code/results/sim_ex_lossb_default_500chain_500iter_inclFALSE.Rds')
rm(MCMC_lossb.prior.multi)
load('code/results/sim_ex_lossb_default_100chain_500iter_inclFALSE.Rds')

# chipman - default
chip.def.list <- list(fun = chipman_prior_tree, param = c(0.95, 2))

st = Sys.time()
MCMC_chip.def <- multichain_MCMC_known_var_windows(n.chain = n.chain_par,
                                           n.iter = n.iter_par, 
                                           X = X, 
                                           Y = Y_ex, 
                                           Y.var = Y.var, 
                                           mu.prior.mean = 0, 
                                           mu.prior.var = var(Y_ex), 
                                           prior_list = chip.def.list, 
                                           moves.prob = NULL, starting.tree = NULL,
                                           include.split = incl.split_par,
                                           cont.unif = cont.unif_par,
                                           n.cores = n.cores)
Sys.time() -st
# 1h and 49 mins
save(MCMC_chip.def, file = 'code/results/sim_ex_chip_def_500chain_500iter_inclFALSE.Rds')
rm(MCMC_chip.def)
load('code/results/sim_ex_chip_def_100chain_500iter_inclFALSE.Rds')


#################################
## Different parametrisations  ##
#################################

# lower omega -- so similar to chipman default
pp.low.omeg <- c(0.5, pp.2[2])
lossb.prior.list.low.omeg <- list(fun = joint.prior.new.tree, param = pp.low.omeg)
st = Sys.time()
MCMC_lossb.prior.low.omeg <- multichain_MCMC_known_var_windows(n.chain = n.chain_par,
                                                    n.iter = n.iter_par, 
                                                    X = X, 
                                                    Y = Y_ex, 
                                                    Y.var = Y.var, 
                                                    mu.prior.mean = 0, 
                                                    mu.prior.var = var(Y_ex), 
                                                    prior_list = lossb.prior.list.low.omeg, 
                                                    moves.prob = NULL, starting.tree = NULL,
                                                    include.split = incl.split_par,
                                                    cont.unif = cont.unif_par,
                                                    n.cores = n.cores)
Sys.time() -st
save(MCMC_lossb.prior.low.omeg, file = 'code/results/sim_ex_lossb_lowomeg_500chain_500iter_inclFALSE.Rds')
rm(MCMC_lossb.prior.low.omeg)

# low omega high gamma
pp.low.omeg.high.gam <- c(0.5, 1.5)
lossb.prior.list.high.gam <- list(fun = joint.prior.new.tree, param = pp.low.omeg.high.gam)
st = Sys.time()
MCMC_lossb.prior.high.gam <- multichain_MCMC_known_var_windows(n.chain = n.chain_par,
                                                       n.iter = n.iter_par, 
                                                       X = X, 
                                                       Y = Y_ex, 
                                                       Y.var = Y.var, 
                                                       mu.prior.mean = 0, 
                                                       mu.prior.var = var(Y_ex), 
                                                       prior_list = lossb.prior.list.high.gam, 
                                                       moves.prob = NULL, starting.tree = NULL,
                                                       include.split = incl.split_par,
                                                       cont.unif = cont.unif_par,
                                                       n.cores = n.cores)
Sys.time() -st
save(MCMC_lossb.prior.high.gam, file = 'code/results/sim_ex_lossb_highgam_500chain_500iter_inclFALSE.Rds')
rm(MCMC_lossb.prior.high.gam)


# chipman lower alpha to match loss-based defult

mean(rdepth.chip(10000, 0.25, 2))
mean(rdepth(10000, pp.2[1], pp.2[2]))
mean(rnterm.chip(10000, 0.25, 2))
mean(rnterm(10000, pp.2[1]))

chip.prior.low.alpha <- list(fun = chipman_prior_tree, param = c(0.25, 2))

st = Sys.time()
MCMC_chip.low.alpha <- multichain_MCMC_known_var_windows(n.chain = n.chain_par,
                                                 n.iter = n.iter_par, 
                                                 X = X, 
                                                 Y = Y_ex, 
                                                 Y.var = Y.var, 
                                                 mu.prior.mean = 0, 
                                                 mu.prior.var = var(Y_ex), 
                                                 prior_list = chip.prior.low.alpha, 
                                                 moves.prob = NULL, starting.tree = NULL,
                                                 include.split = incl.split_par,
                                                 cont.unif = cont.unif_par,
                                                 n.cores = n.cores)
Sys.time() -st
# 1h and 49 mins
save(MCMC_chip.low.alpha, file = 'code/results/sim_ex_chip_lowalpha_500chain_500iter_inclFALSE.Rds')
rm(MCMC_chip.low.alpha)
# low beta just to change
chip.prior.low.beta <- list(fun = chipman_prior_tree, param = c(0.25, 0.5))

st = Sys.time()
MCMC_chip.low.beta <- multichain_MCMC_known_var_windows(n.chain = n.chain_par,
                                                 n.iter = n.iter_par, 
                                                 X = X, 
                                                 Y = Y_ex, 
                                                 Y.var = Y.var, 
                                                 mu.prior.mean = 0, 
                                                 mu.prior.var = var(Y_ex), 
                                                 prior_list = chip.prior.low.beta, 
                                                 moves.prob = NULL, starting.tree = NULL,
                                                include.split = incl.split_par,
                                                cont.unif = cont.unif_par,
                                                n.cores = n.cores)
Sys.time() -st
# 10 mins
save(MCMC_chip.low.beta, file = 'code/results/sim_ex_chip_lowbeta_500chain_500iter_inclFALSE.Rds')
rm(MCMC_chip.low.beta)

# prob.split
# prior.split.rule



############################################
# IMPORT MODELS RESULT AND ANALYSE RESULTS #
############################################

load('code/results/sim_ex_chip_def_500chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_chip_lowbeta_500chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_chip_lowalpha_500chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_lossb_default_500chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_lossb_lowomeg_500chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_lossb_highgam_500chain_500iter_inclFALSE.Rds')


## Analyse predictors rate
## gets the predictors used in a tree - works in general
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

# boolean if the predictor is used
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

# define burn-in vector
burn_in = 250
idx.burn_in = rep(c(rep(FALSE, burn_in), rep(TRUE, n.iter_par - burn_in)), 
                  n.chain_par)


# get predictor frequency for CL prior 
pred_chip.def <- get_pred_freq(MCMC_chip.def$trees[idx.burn_in], 1:3)
pred_chip.low.alpha <- get_pred_freq(MCMC_chip.low.alpha$trees[idx.burn_in], 1:3)
pred_chip.low.beta <- get_pred_freq(MCMC_chip.low.beta$trees[idx.burn_in], 1:3)


# get predictor frequency for LB prior 
pred_lb.def <- get_pred_freq(MCMC_lossb.prior.multi$trees[idx.burn_in], 1:3)
pred_lb.high.gam <- get_pred_freq(MCMC_lossb.prior.high.gam$trees[idx.burn_in], 1:3)
pred_lb.low.omeg <- get_pred_freq(MCMC_lossb.prior.low.omeg$trees[idx.burn_in], 1:3)



pred_matrix <- cbind(pred_chip.def,
                     pred_chip.low.alpha, 
                     pred_chip.low.beta,
                     pred_lb.def, 
                     pred_lb.high.gam,
                     pred_lb.low.omeg)

variables <- c('X1', 'X2', 'X3')


colnames(pred_matrix) <- c('CH - a = 0.95, b = 2',
                    'CH - a = 0.25, b = 2',
                    'CH - a = 0.25, b = 0.5',
                    'LB - a = 1.56, g = 0.62',
                    'LB - a = 0.5, g = 1.5',
                    'LB - a = 0.5, g = 0.62')
                     
rownames(pred_matrix) <- variables
library(tidyverse)
pred_df <- as.data.frame(pred_matrix) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(-Variable, names_to = "Model", values_to = "Value")


model_labels <- c( expression(atop("CL", alpha == 0.25 * "," ~ beta == 0.5)),
                   expression(atop("CL", alpha == 0.25 * "," ~ beta == 2)),
                   expression(atop("CL", alpha == 0.95 * "," ~ beta == 2)),
                   expression(atop("LB", omega == 0.5 * "," ~ gamma == 0.62)),
                   expression(atop("LB", omega == 0.5 * "," ~ gamma == 1.5)),
                   expression(atop("LB", omega == 1.56 * "," ~ gamma == 0.62))
)


#####################
# FIGURE 8 (Bottom) #
#####################
pdf('figures/pred_synth.pdf', height = 5, width = 5*2, useDingbats = FALSE)
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

dev.off()

#### number of events/depth/loglikelihood/missing rate calculations

res.list <- list(MCMC_chip.low.alpha, 
                 MCMC_chip.low.beta,
                 MCMC_chip.def,
                 MCMC_lossb.prior.low.omeg, 
                 MCMC_lossb.prior.high.gam,
                 MCMC_lossb.prior.multi)
names(res.list) <- c('CH - a = 0.25, b = 2',
                     'CH - a = 0.25, b = 0.5',
                     'CH - a = 0.95, b = 2',
                     'LB - a = 0.5, g = 0.62',
                     'LB - a = 0.5, g = 1.5',
                     'LB - a = 1.56, g = 0.62')

df.nl <- apply_fun_models(get_num_terminal_nodes, res.list, 
                          born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)
df.dd <- apply_fun_models(get_depth, res.list, 
                          born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)



 # QUANTITIES FOR TABLE 1 #
df.dd.split <- split(df.dd, df.dd$panel.name)
cbind(mean = vapply(df.dd.split, \(x) mean(x$y), 0),
      q0.025 = vapply(df.dd.split, \(x) quantile(x$y, 0.025), 0),
      q0.975 = vapply(df.dd.split, \(x) quantile(x$y, 0.975), 0),
      pr.above10 = vapply(df.dd.split, \(x) mean(x$y > 6), 0)
)


df.split <- split(df.dd, df.dd$panel.name)
cbind(mean = vapply(df.split, \(x) mean(x$y), 0),
      pr.above10 = vapply(df.split, \(x) mean(x$y > 6), 0)
)


## Analyse difference in number of terminal nodes


# histograms
nterm.chip.ecdf.def <- ecdf(rnterm.chip(10000, 0.95, 2))
nterm.chip.ecdf.lowa <- ecdf(rnterm.chip(10000, 0.25, 2))
nterm.chip.ecdf.lowb <- ecdf(rnterm.chip(10000, 0.25, 0.5))

#pp.2 <- 
df.nl.prior <- rbind(data.frame(x = 1:20, 
                             y = prior.nterm(1:20, pp.2[1]),
                             panel.name = 'LB - a = 1.56, g = 0.62'),
                  data.frame(x = 1:20, 
                             y = prior.nterm(1:20, 0.5),
                             panel.name = 'LB - a = 0.5, g = 0.62'),
                  data.frame(x = 1:20, 
                       y = prior.nterm(1:20, 0.5),
                       panel.name = 'LB - a = 0.5, g = 1.5'),
                  data.frame(x = 1:20, 
                             y = diff(nterm.chip.ecdf.def(0:20)),
                             panel.name = 'CH - a = 0.95, b = 2'),
                  data.frame(x = 1:20, 
                             y = diff(nterm.chip.ecdf.lowa(0:20)),
                             panel.name = 'CH - a = 0.25, b = 2'),
                  data.frame(x = 1:20, 
                             y = diff(nterm.chip.ecdf.lowb(0:20)),
                             panel.name = 'CH - a = 0.25, b = 0.5'))



# summary stat 
n_samp <- 10000
df.prior.nl.samp <- rbind(data.frame(models = "CH - a = 0.25, b = 0.5 - prior",
                                     samp = rnterm.chip(n_samp, 0.25,0.5)),
                          data.frame(models = "CH - a = 0.25, b = 2 - prior",
                                     samp = rnterm.chip(n_samp, 0.25,2)),
                          data.frame(models = "CH - a = 0.95, b = 2 - prior",
                                     samp = rnterm.chip(n_samp, 0.95,2)),
                          data.frame(models = "LB - a = 0.5, g = 0.62 - prior",
                                     samp = rnterm(n_samp, 0.5)),
                          data.frame(models = "LB - a = 0.5, g = 1.5 - prior",
                                     samp = rnterm(n_samp, 0.5)),
                          data.frame(models = "LB - a = 1.56, g = 0.62 - prior",
                                     samp = rnterm(n_samp, 1.56)))


df.nl.stat.prior <- data.frame(models = aggregate(df.prior.nl.samp$samp, list(df.prior.nl.samp$models), mean)[,1],
                              mean = aggregate(df.prior.nl.samp$samp, list(df.prior.nl.samp$models), mean)[,2],
                              median = aggregate(df.prior.nl.samp$samp, list(df.prior.nl.samp$models), median)[,2],
                              pr.up8 = aggregate(df.prior.nl.samp$samp, list(df.prior.nl.samp$models), \(x){mean(x >= 8)})[,2])

df.nl.stat.post <- data.frame(models = aggregate(df.nl$y, list(df.nl$panel.name), mean)[,1],
                              mean = aggregate(df.nl$y, list(df.nl$panel.name), mean)[,2],
                              median = aggregate(df.nl$y, list(df.nl$panel.name), median)[,2],
                              pr.up8 = aggregate(df.nl$y, list(df.nl$panel.name), \(x){mean(x >= 8)})[,2])

# QUANTITIES FOR TABLE 1
rbind(df.nl.stat.prior, df.nl.stat.post)


# to plot
df.nl$panel.name <- factor(df.nl$panel.name,
                           labels = c( 'CH - a = 0.25, b = 0.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 0.5$')),
                                       'CH - a = 0.25, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 2$')),
                                       'CH - a = 0.95, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 2$')),
                                       'LB - a = 0.5, g = 0.62' = parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 0.62$')),
                                      'LB - a = 0.5, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 1.5$')),
                                      'LB - a = 1.56, g = 0.62'= parse(text=latex2exp::TeX('$LB - \\omega = 1.56, \\gamma = 0.62$'))))

df.nl.prior$panel.name <- factor(df.nl.prior$panel.name,
                           labels = c('CH - a = 0.25, b = 0.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 0.5$')), 
                                      'CH - a = 0.25, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 2$')),
                                      'CH - a = 0.95, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 2$')),
                                       'LB - a = 0.5, g = 0.62' = parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 0.62$')),
                                       'LB - a = 0.5, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 1.5$')),
                                       'LB - a = 1.56, g = 0.62'= parse(text=latex2exp::TeX('$LB - \\omega = 1.56, \\gamma = 0.62$'))))



make_single_plot <- function(df, df_prior, model_label, annot_expr, rescaling = 0.2,
                             xlab_par = expression(n[L]),
                             annot_param = list(x.loc = 16, y.loc = 0.6, size = 1)) {
  ggplot(df, aes(x = y, y = after_stat(density))) +
    geom_point(stat = "bin", binwidth = 1, colour = "black", size = 1) +
    geom_line(data = df_prior, aes(x, (y/max(y))*rescaling ), colour = "deepskyblue", linewidth = 0.4) +
    geom_point(data = df_prior, aes(x, (y/max(y))*rescaling ), colour = "deepskyblue", size = 0.4) +
    geom_vline(xintercept = get_num_terminal_nodes(tree_ex), color = "orange") +
    scale_x_continuous(breaks = seq(0, 21, by = 3), limits = c(0, 21)) +
    #scale_y_log10(limits = c(1e-3, NA)) +
    xlab(xlab_par) +
    ylab('pmf') +
    theme_classic(base_size = 13) +
    #theme(
    #  plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    #  axis.title.y = element_blank()
    #) +
    #ylim(0,0.25) + 
    annotate(
      "label",
      x = annot_param$x.loc, 
      y = annot_param$y.loc,        
      label = annot_expr,
      parse = TRUE,
      hjust = 0,
      vjust = 0,
      label.size = 0.5, size = annot_param$size,
      fill = "white"
    )
}


models <- as.character(unique(df.nl$panel.name))



annotations <- list(
  expression(atop(scriptstyle('CL - '~alpha~'=0.25,'~beta~'=0.5'),
                  atop(E(n[L]) ~ '=' ~  1.36,
                       E(n[L] * '|' * bold(y)) ~ '=' ~  7.5))),
  expression(atop(scriptstyle('CL - '~alpha~'=0.25,'~beta~'=2'),
                  atop(E(n[L]) ~ '=' ~  1.29,
                       E(n[L] * '|' * bold(y)) ~ '=' ~  7.24))),
  expression(atop(scriptstyle('CL - '~alpha~'=0.95,'~beta~'=2'),
                  atop(E(n[L]) ~ '=' ~  2.51,
                       E(n[L] * '|' * bold(y)) ~ '=' ~  8.4))),
  expression(atop(scriptstyle('LB - '~omega~'=0.5,'~gamma~'=0.62'),
                  atop(E(n[L])  ~ '=' ~   2.52,
                       E(n[L] * '|' * bold(y)) ~ '=' ~ 8.1))),
  expression(atop(scriptstyle('LB - '~omega~'=0.5,'~gamma~'=1.5'),
                  atop(E(n[L])  ~ '=' ~   2.52,
                       E(n[L] * '|' * bold(y))~ '=' ~  8))),
  expression(atop(scriptstyle('LB - '~omega~'=1.56,'~gamma~'=0.62'),
                  atop(E(n[L]) ~ '=' ~ 1.26,
                       E(n[L] * '|' * bold(y))~ '=' ~ 5.95)))
)



plots <- purrr::map2(models, annotations, function(m, a) {
  df_sub <- df.nl %>% filter(panel.name == m)
  df_prior_sub <- df.nl.prior %>% filter(panel.name == m)
  make_single_plot(df = df_sub, 
                   df_prior = df_prior_sub, 
                   model_label = m, 
                   annot_expr =  a, annot_param = list(x.loc = 9, y.loc = 0.14, size = 4))
})

############
# FIGURE 7 #
############ 
mult.dim = 1.8
pdf('figures/fig_posterior_nl_logscale_new.pdf',  width = 6*mult.dim, height = 3*mult.dim)
wrap_plots(plots, nrow = 2)
dev.off()



## Analyse difference in depth
df.depth <- df.dd#apply_fun_models(get_depth, res.list, born.out.pc = 250,)

df.prior.depth.samp <- rbind(data.frame(models = "CH - a = 0.25, b = 0.5 - prior",
                                     samp = rdepth.chip(n_samp, 0.25,0.5)),
                          data.frame(models = "CH - a = 0.25, b = 2 - prior",
                                     samp = rdepth.chip(n_samp, 0.25,2)),
                          data.frame(models = "CH - a = 0.95, b = 2 - prior",
                                     samp = rdepth.chip(n_samp, 0.95,2)),
                          data.frame(models = "LB - a = 0.5, g = 0.62 - prior",
                                     samp = rdepth(n_samp, 0.5, 0.62)),
                          data.frame(models = "LB - a = 0.5, g = 1.5 - prior",
                                     samp = rdepth(n_samp, 0.5, 1.5)),
                          data.frame(models = "LB - a = 1.56, g = 0.62 - prior",
                                     samp = rdepth(n_samp, 1.56, 0.62)))


df.depth.stat.prior <- data.frame(models = aggregate(df.prior.depth.samp$samp, list(df.prior.depth.samp$models), mean)[,1],
                               mean = aggregate(df.prior.depth.samp$samp, list(df.prior.depth.samp$models), mean)[,2],
                               median = aggregate(df.prior.depth.samp$samp, list(df.prior.depth.samp$models), median)[,2],
                               pr.up6 = aggregate(df.prior.depth.samp$samp, list(df.prior.depth.samp$models), 
                                                  \(x){mean(x > 6)})[,2])

df.depth.stat.post <- data.frame(models = aggregate(df.depth$y, list(df.depth$panel.name), mean)[,1],
                              mean = aggregate(df.depth$y, list(df.depth$panel.name), mean)[,2],
                              median = aggregate(df.depth$y, list(df.depth$panel.name), median)[,2],
                              pr.up6 = aggregate(df.depth$y, list(df.depth$panel.name), \(x){mean(x > 6)})[,2])


# QUANTITIES FOR TABLE 1
rbind(df.depth.stat.prior, df.depth.stat.post)


chip.depth.def <- ecdf(rdepth.chip(10000, 0.95, 2))
chip.depth.lowa <- ecdf(rdepth.chip(10000, 0.25, 2))
chip.depth.lowb <- ecdf(rdepth.chip(10000, 0.25, 0.5))
lb.depth.def <- ecdf(rdepth(10000, 1.56, 0.62))
lb.depth.lowa <- ecdf(rdepth(10000, 0.5, 0.62))
lb.depth.highg <- ecdf(rdepth(10000, 0.5, 1.5))

df.depth.prior <- rbind(data.frame(x = 0:10, 
                                y = diff(lb.depth.def(-1:10)),
                                panel.name = 'LB - a = 1.56, g = 0.62'),
                     data.frame(x = 0:10, 
                                y = diff(lb.depth.highg(-1:10)),
                                panel.name = 'LB - a = 0.5, g = 0.62'),
                     data.frame(x = 0:10, 
                                y = diff(lb.depth.highg(-1:10)),
                                panel.name = 'LB - a = 0.5, g = 1.5'),
                     data.frame(x = 0:10, 
                                y = diff(chip.depth.def(-1:10)),
                                panel.name = 'CH - a = 0.95, g = 2'),
                     data.frame(x = 0:10, 
                                y = diff(chip.depth.lowa(-1:10)),
                                panel.name = 'CH - a = 0.25, b = 2'),
                     data.frame(x = 0:10, 
                                y = diff(chip.depth.lowb(-1:10)),
                                panel.name = 'CH - a = 0.25, b = 0.5'))



df.depth$panel.name <- factor(df.depth$panel.name,
                           labels = c( 'CH - a = 0.25, b = 0.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 0.5$')),
                                       'CH - a = 0.25, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 2$')),
                                       'CH - a = 0.95, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 2$')),
                                       'LB - a = 0.5, g = 0.62' = parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 0.62$')),
                                       'LB - a = 0.5, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 1.5$')),
                                       'LB - a = 1.56, g = 0.62'= parse(text=latex2exp::TeX('$LB - \\omega = 1.56, \\gamma = 0.62$'))))

df.depth.prior$panel.name <- factor(df.depth.prior$panel.name,
                                 labels = c('CH - a = 0.25, b = 0.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 0.5$')), 
                                            'CH - a = 0.25, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 2$')),
                                            'CH - a = 0.95, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 2$')),
                                            'LB - a = 0.5, g = 0.62' = parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 0.62$')),
                                            'LB - a = 0.5, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 1.5$')),
                                            'LB - a = 1.56, g = 0.62'= parse(text=latex2exp::TeX('$LB - \\omega = 1.56, \\gamma = 0.62$'))))


annotations <- list(
  expression(atop(scriptstyle('CL - '~alpha~'=0.25,'~beta~'=0.5'),
                  atop(E(D) ~ '=' ~  0.35,
                       E(D * '|' * bold(y)) ~ '=' ~  4.16))),
  expression(atop(scriptstyle('CL - '~alpha~'=0.25,'~beta~'=2'),
                  atop(E(D) ~ '=' ~ 0.28,
                       E(D * '|' * bold(y)) ~ '=' ~  4.17))),
  expression(atop(scriptstyle('CL - '~alpha~'=0.95,'~beta~'=2'),
                  atop(E(D) ~ '=' ~  1.45,
                       E(D * '|' * bold(y)) ~ '=' ~  4.66))),
  expression(atop(scriptstyle('LB - '~omega~'=0.5,'~gamma~'=0.62'),
                  atop(E(D)  ~ '=' ~   1.21,
                       E(D * '|' * bold(y)) ~ '=' ~ 4.49))),
  expression(atop(scriptstyle('LB - '~omega~'=0.5,'~gamma~'=1.5'),
                  atop(E(D)  ~ '=' ~   1.2,
                       E(D * '|' * bold(y))~ '=' ~  4.33))),
  expression(atop(scriptstyle('LB - '~omega~'=1.56,'~gamma~'=0.62'),
                  atop(E(D) ~ '=' ~ 0.25,
                       E(D * '|' * bold(y))~ '=' ~ 3.40)))
)



plots <- purrr::map2(models, annotations, function(m, a) {
  df_sub <- df.depth %>% filter(panel.name == m)
  df_prior_sub <- df.depth.prior %>% filter(panel.name == m)
  make_single_plot(df = df_sub, 
                   df_prior = df_prior_sub, 
                   model_label = m, rescaling = 0.3, 
                   annot_expr =  a, xlab_par = expression(D),
                   annot_param = list(x.loc = 9, y.loc = 0.14, size = 4))
})

##############################################
# SUPPLEMENTARY MATERIAL Appendix E FIGURE 1 #
##############################################

mult.dim = 1.8
pdf('figures/fig_posterior_depth_new.pdf',  width = 6*mult.dim, height = 3*mult.dim)
wrap_plots(plots, nrow = 2)
dev.off()


# traceplots 
# number of terminal nodes
df.nl <- apply_fun_models(get_num_terminal_nodes, res.list, born.out.pc = 0)
df.nl.split <- split(df.nl, df.nl$panel.name)
df.trace <- Reduce(rbind,lapply(df.nl.split, \(df) df[1:5000,]))


df.trace$panel.name <- factor(df.trace$panel.name,
                           labels = c( 'CH - a = 0.25, b = 0.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 0.5$')),
                                       'CH - a = 0.25, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 2$')),
                                       'CH - a = 0.95, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 2$')),
                                       'LB - a = 0.5, g = 0.62' = parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 0.62$')),
                                       'LB - a = 0.5, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 1.5$')),
                                       'LB - a = 1.56, g = 0.62'= parse(text=latex2exp::TeX('$LB - \\omega = 1.56, \\gamma = 0.62$'))))
##############################################
# SUPPLEMENTARY MATERIAL Appendix F FIGURE 2 #
##############################################

pdf(file = 'figures/fig_trace_nl.pdf', width = 6, height = 6)
ggplot(df.trace, aes(x = x, y = y)) +
  geom_vline(xintercept = seq(0,5000,by = n.iter_par), color = 'grey', alpha = 0.7, linewidth = 0.3) + 
  geom_line() + 
  ylab('Number of terminal nodes') +
  xlab('Iteration') + 
  geom_hline(yintercept = get_num_terminal_nodes(tree_ex), color = 'orange', linetype = 2) + 
  facet_wrap(facets = ~panel.name, ncol = 1, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0,5000,by = n.iter_par)) + 
  theme_classic()
dev.off()


# depth

df.depth <- apply_fun_models(get_depth, res.list, born.out.pc = 0)
df.depth.split <- split(df.depth, df.depth$panel.name)
df.depth.trace <- Reduce(rbind,lapply(df.depth.split, \(df) df[1:5000,]))




df.depth.trace$panel.name <- factor(df.depth.trace$panel.name,
                           labels = c( 'CH - a = 0.25, b = 0.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 0.5$')),
                                       'CH - a = 0.25, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 2$')),
                                       'CH - a = 0.95, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 2$')),
                                       'LB - a = 0.5, g = 0.62' = parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 0.62$')),
                                       'LB - a = 0.5, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 1.5$')),
                                       'LB - a = 1.56, g = 0.62'= parse(text=latex2exp::TeX('$LB - \\omega = 1.56, \\gamma = 0.62$'))))

##############################################
# SUPPLEMENTARY MATERIAL Appendix F FIGURE 3 #
##############################################

pdf(file = 'figures/fig_trace_depth.pdf', width = 6, height = 6)
ggplot(df.depth.trace, aes(x = x, y = y)) +
  geom_vline(xintercept = seq(0,5000,by = n.iter_par), color = 'grey', alpha = 0.7, linewidth = 0.3) + 
  geom_line() + 
  ylab('Depth') +
  xlab('Iteration') + 
  geom_hline(yintercept = get_depth(tree_ex), color = 'orange', linetype = 2) + 
  facet_wrap(facets = ~panel.name, ncol = 1, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0,5000,by=n.iter_par)) + 
  theme_classic()
dev.off()



# loglikelohood 

df.loglik <- apply_fun_models(\(x) cart_log_lik(x, Y_ex, X, 2), res.list, born.out.pc = 0)
true.loglik <- cart_log_lik(tree_ex, Y_ex, X, 2)
df.loglik.split <- split(df.loglik, df.loglik$panel.name)
df.loglik.trace <- Reduce(rbind,lapply(df.loglik.split, \(df) df[1:5000,]))

df.loglik.trace$panel.name <- factor(df.loglik.trace$panel.name,
                              labels = c( 'CH - a = 0.25, b = 0.5' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 0.5$')),
                                          'CH - a = 0.25, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.25, \\beta = 2$')),
                                          'CH - a = 0.95, b = 2' = parse(text=latex2exp::TeX('$CL - \\alpha = 0.95, \\beta = 2$')),
                                          'LB - a = 0.5, g = 0.62' = parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 0.62$')),
                                          'LB - a = 0.5, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.5, \\gamma = 1.5$')),
                                          'LB - a = 1.56, g = 0.62'= parse(text=latex2exp::TeX('$LB - \\omega = 1.56, \\gamma = 0.62$'))))

save(df.loglik.trace, file = 'code/results/df_for_loglik_trace.Rds')

##############################################
# SUPPLEMENTARY MATERIAL Appendix F FIGURE 4 #
##############################################

pdf(file = 'figures/fig_trace_loglik.pdf', width = 6, height = 6)
ggplot(df.loglik.trace, aes(x = x, y = y)) +
  geom_vline(xintercept = seq(0,5000,by = 500), color = 'grey', alpha = 0.7, linewidth = 0.3) + 
  geom_line() + 
  ylab('Log Likelihood') +
  xlab('Iteration') + 
  geom_hline(yintercept = true.loglik, color = 'orange', linetype = 2) + 
  facet_wrap(facets = ~panel.name, ncol = 1, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0,5000,by=500)) + 
  theme_classic()
dev.off()


# calculate prior for each tree so to calculate posterior distribution and 
# for each parametrisation find the tree with highest posterior

df.prior_fit.lb.def <- lapply(1:length(MCMC_lossb.prior.multi$trees), 
                              \(x) joint.prior.new.tree(MCMC_lossb.prior.multi$trees[[x]], 1.56, 0.62))
df.prior_fit.lb.def <- data.frame(prior = do.call(c, df.prior_fit.lb.def),
                                  panel.name = 'LB - a = 1.56, g = 0.62')

df.prior_fit.lb.high.gam <- lapply(1:length(MCMC_lossb.prior.high.gam$trees), 
                                   \(x) joint.prior.new.tree(MCMC_lossb.prior.high.gam$trees[[x]], 1.56, 0.62))
df.prior_fit.lb.high.gam <- data.frame(prior = do.call(c, df.prior_fit.lb.high.gam),
                                  panel.name = 'LB - a = 0.5, g = 1.5')

df.prior_fit.lb.low.omeg <- lapply(1:length(MCMC_lossb.prior.low.omeg$trees), 
                                   \(x) joint.prior.new.tree(MCMC_lossb.prior.low.omeg$trees[[x]], 1.56, 0.62))
df.prior_fit.lb.low.omeg <- data.frame(prior = do.call(c, df.prior_fit.lb.low.omeg),
                                  panel.name = 'LB - a = 0.5, g = 0.62')

df.prior_fit.chip.def <- lapply(MCMC_chip.def$trees, \(x) chipman_prior_tree(x, 0.95, 2))
df.prior_fit.chip.def <- data.frame(prior = do.call(c, df.prior_fit.chip.def),
                                       panel.name = 'CH - a = 0.95, b = 2')

df.prior_fit.chip.low.alpha <- lapply(MCMC_chip.low.alpha$trees, \(x) chipman_prior_tree(x, 0.95, 2))
df.prior_fit.chip.low.alpha <- data.frame(prior = do.call(c, df.prior_fit.chip.low.alpha),
                                    panel.name = 'CH - a = 0.25, b = 2')

df.prior_fit.chip.low.beta <- lapply(MCMC_chip.low.beta$trees, \(x) chipman_prior_tree(x, 0.95, 2))
df.prior_fit.chip.low.beta <- data.frame(prior = do.call(c, df.prior_fit.chip.low.beta),
                                    panel.name = 'CH - a = 0.25, b = 0.5')


df.prior_fit <- rbind(df.prior_fit.lb.def, 
                      df.prior_fit.lb.high.gam,
                      df.prior_fit.lb.low.omeg,
                      df.prior_fit.chip.def,
                      df.prior_fit.chip.low.alpha,
                      df.prior_fit.chip.low.beta)

## find tree with highest loglikelihood
panel_names <- unique(df.loglik$panel.name)
max_lik_tree_list <- list()
for(panel in panel_names){
  df.loglik.red <- df.loglik[df.loglik$panel.name == panel,]
  df.prior_fit.red <- df.prior_fit[df.prior_fit$panel.name == panel,]
  posterior = exp(df.loglik.red$y)*df.prior_fit.red$prior
  print(df.loglik.red[which.max(posterior),])
  print(df.loglik.red[which.max(df.loglik.red$y),])
  print('-----')
}

df.loglik[df.loglik$panel.name == 'LB - a = 1.56, g = 0.62',][c(54250,108715),]
df.prior_fit.lb.def[c(54250,65339),]


# best tree LB default
plot.max.loglik.lb.def <- get_tree_plot(MCMC_lossb.prior.multi$trees[54250][[1]])
plot.mode.lb.def <- get_tree_plot(MCMC_lossb.prior.multi$trees[65339][[1]])

# best tree LB  omega = 0.5, gamma = 1.5
plot.max.loglik.lb.high.gam <-get_tree_plot(MCMC_lossb.prior.high.gam$trees[7837][[1]])
plot.mode.lb.high.gam <- get_tree_plot(MCMC_lossb.prior.high.gam$trees[183233][[1]])
# best tree LB  omega = 0.5, gamma = 0.62
plot.max.loglik.lb.low.omeg <-get_tree_plot(MCMC_lossb.prior.low.omeg$trees[233938][[1]])
plot.mode.lb.low.omeg <- get_tree_plot(MCMC_lossb.prior.low.omeg$trees[18027][[1]])

# best tree Chipman default
plot.max.loglik.chip.def <-get_tree_plot(MCMC_chip.def$trees[202364][[1]])
plot.mode.chip.def <- get_tree_plot(MCMC_chip.def$trees[116112][[1]])

# best tree Chipman alpha = 0.25, beta = 2
plot.max.loglik.chip.low.alpha <-get_tree_plot(MCMC_chip.low.alpha$trees[72375][[1]])
plot.mode.chip.low.alpha <-get_tree_plot(MCMC_chip.low.alpha$trees[114466 ][[1]])

# best tree Chipman alpha = 0.25, beta = 0.5
plot.max.loglik.chip.low.beta <-get_tree_plot(MCMC_chip.low.beta$trees[17322][[1]])
plot.mode.chip.low.beta <-get_tree_plot(MCMC_chip.low.beta$trees[92862 ][[1]])

p1 <- combineWidgets(plot.mode.lb.def, plot.mode.lb.high.gam, plot.mode.lb.low.omeg,
               plot.mode.chip.def, plot.mode.chip.low.alpha, plot.mode.chip.low.beta,
               nrow = 2, ncol = 3)

##############################################
# SUPPLEMENTARY MATERIAL Appendix G FIGURE 6 #
##############################################
saveWidget(p1, "figures/combined_plot_mode.html", selfcontained = TRUE)
webshot("figures/combined_plot_mode.html", file = "figures/combined_plot_mode.png", vwidth = 2400, vheight = 1600)

library(manipulateWidget)
library(htmlwidgets)
library(webshot2) 

p2 <- combineWidgets(plot.max.loglik.lb.def, 
                     plot.max.loglik.lb.high.gam, 
                     plot.max.loglik.lb.low.omeg,
                     plot.max.loglik.chip.def, 
                     plot.max.loglik.chip.low.alpha, 
                     plot.max.loglik.chip.low.beta, nrow = 2, ncol = 3)

##############################################
# SUPPLEMENTARY MATERIAL Appendix G FIGURE 5 #
##############################################
saveWidget(p2, "figures/combined_plot_maxloglik.html", selfcontained = TRUE)
webshot("figures/combined_plot_maxloglik.html", file = "figures/combined_plot_maxloglik.png", vwidth = 2400, vheight = 1600)



###########################
# out-of-sample forecasts #
###########################

idx.burn_in = rep(c(rep(FALSE, 250), rep(TRUE, 250)), 100)

set.seed(1)

x1.2 <- c(runif(200*data.mult, 0.1, 0.4), runif(100*data.mult, 0.6, 0.9))
x2.2 <- c(runif(100*data.mult, 0.1, 0.4), runif(100*data.mult, 0.6, 0.9), runif(100*data.mult, 0.1, 0.9))
x3.2 <- c(runif(200*data.mult, 0.6, 0.9), runif(100*data.mult, 0.1, 0.4))
# store them in data.frame
X.2 = data.frame(x1 = x1.2, x2 = x2.2, x3 = x3.2)
# normalise data
X.2 <- as.data.frame(apply(X, 2, \(x) (x - min(x))/(max(x) - min(x))))


Y.var = 0.25
# sample observations
Y_ex.2 <- sample_CART(tree_top = tree_ex, X = X, sigma_ = Y.var)
# store in data.frame
data.df.2 <- cbind(Y = Y_ex.2, X)

# sample observations
Y_ex.2 <- sample_CART(tree_top = tree_ex, X = X.2, sigma_ = Y.var)
# store in data.frame
data.df.2 <- cbind(Y = Y_ex.2, X.2)





get_value_tree_posterior <- function(tree_set, X, n.cores){
  cl <- makeCluster(n.cores)
  # Export needed variables and functions to the workers
  clusterExport(cl, varlist = c("tree_set", "X", "get_value_tree", 'g.T', 'eval_cond'), envir = environment())
  
  # Parallel version of your lapply
  Y_pred_post <- parLapply(cl, tree_set, function(x) {
    get_value_tree(x, X)
  })
  # Stop the cluster when done
  stopCluster(cl)
  return(Y_pred_post)
}


### repeated out-of-sample experiments 
Y_list = list()
X_list = list()
n.burnin = 250
idx.burn_in = rep(c(rep(FALSE, n.burnin), rep(TRUE, n.iter_par - n.burnin)), n.chain_par)

for(rep in 8:10){
  cat("------------------ REP : ", rep, " ---------------")
  set.seed(rep)
  
  x1.n <- c(runif(200*data.mult, 0.1, 0.4), runif(100*data.mult, 0.6, 0.9))
  x2.n <- c(runif(100*data.mult, 0.1, 0.4), runif(100*data.mult, 0.6, 0.9), runif(100*data.mult, 0.1, 0.9))
  x3.n <- c(runif(200*data.mult, 0.6, 0.9), runif(100*data.mult, 0.1, 0.4))
  # store them in data.frame
  X.n = data.frame(x1 = x1.n, x2 = x2.n, x3 = x3.n)
  X_list[[rep]] <- X.n
  
  # sample observations
  Y_ex.n <- sample_CART(tree_top = tree_ex, X = X.n, sigma_ = Y.var)
  
  # store in data.frame
  Y_list[[rep]] <- Y_ex.n
  
  # LOSS-BASED PRIOR # 
  
  ##############
  st = Sys.time()
  load('code/results/sim_ex_lossb_default_500chain_500iter_inclFALSE.Rds')
  tree_Set = MCMC_lossb.prior.multi$trees[idx.burn_in]
  Y_pred.lb.def = get_value_tree_posterior(tree_set = tree_Set, 
                                            X = X.n, n.cores = 10)
  print(Sys.time() - st)
  file_path = paste0('code/results/out_sample_', rep, '_lb_def.Rds')
  rm(MCMC_lossb.prior.multi)
  save(Y_pred.lb.def, file = file_path)
  
  ###############
  st = Sys.time()
  load('code/results/sim_ex_lossb_highgam_500chain_500iter_inclFALSE.Rds')
  tree_Set = MCMC_lossb.prior.high.gam$trees[idx.burn_in]
  Y_pred.lb.high.gam = get_value_tree_posterior(tree_set = tree_Set, 
                                                 X = X.n, n.cores = 10)
  print(Sys.time() - st)
  rm(MCMC_lossb.prior.high.gam)
  file_path = paste0('code/results/out_sample_', rep, '_lb_high_gam.Rds')
  save(Y_pred.lb.high.gam, file = file_path)
  
  ##########
  st = Sys.time()
  load('code/results/sim_ex_lossb_lowomeg_500chain_500iter_inclFALSE.Rds')
  tree_Set = MCMC_lossb.prior.low.omeg$trees[idx.burn_in]
  Y_pred.lb.low.omeg = get_value_tree_posterior(tree_set = tree_Set, 
                                                 X = X.n, n.cores = 10)
  print(Sys.time() - st)
  rm(MCMC_lossb.prior.low.omeg)
  file_path = paste0('code/results/out_sample_', rep, '_lb_low_omeg.Rds')
  save(Y_pred.lb.low.omeg, file = file_path)
  
  
  # CHIPMAN PRIOR
  #############
  st = Sys.time()
  load('code/results/sim_ex_chip_def_500chain_500iter_inclFALSE.Rds')
  tree_Set = MCMC_chip.def$trees[idx.burn_in]
  Y_pred.chip.def = get_value_tree_posterior(tree_set = tree_Set, 
                                              X = X.n, n.cores = 10)
  
  print(Sys.time() - st)
  rm(MCMC_chip.def)
  file_path = paste0('code/results/out_sample_', rep, '_chip_def.Rds')
  save(Y_pred.chip.def, file = file_path)
  
  ###############
  st = Sys.time()
  load('code/results/sim_ex_chip_lowalpha_500chain_500iter_inclFALSE.Rds')
  tree_Set = MCMC_chip.low.alpha$trees[idx.burn_in]
  Y_pred.chip.low.alpha = get_value_tree_posterior(tree_set = tree_Set, 
                                                    X = X.n, n.cores = 10)
  
  print(Sys.time() - st)
  rm(MCMC_chip.low.alpha)
  file_path = paste0('code/results/out_sample_', rep, '_chip_low_alpha.Rds')
  save(Y_pred.chip.low.alpha, file = file_path)
  
  ##################
  st = Sys.time()
  load('code/results/sim_ex_chip_lowbeta_500chain_500iter_inclFALSE.Rds')
  tree_Set = MCMC_chip.low.beta$trees[idx.burn_in]
  Y_pred.chip.low.beta = get_value_tree_posterior(tree_set = tree_Set, 
                                                   X = X.n, n.cores = 10)
  
  print(Sys.time() - st)
  rm(MCMC_chip.low.beta)
  file_path = paste0('code/results/out_sample_', rep, '_low_beta.Rds')
  save(Y_pred.chip.low.beta, file = file_path)
  
}

save(Y_list, file = 'code/results/Y_outofsample_synt_data.Rds')
save(X_list, file = 'code/results/X_outofsample_synt_data.Rds')

load('code/results/Y_outofsample_synt_data.Rds')
load('code/results/X_outofsample_synt_data.Rds')
mse.list = list()
for(rep in 1:10){
  load(paste0('code/results/out_sample_', rep, '_low_beta.Rds'))
  load(paste0('code/results/out_sample_', rep, '_chip_low_alpha.Rds'))
  load(paste0('code/results/out_sample_', rep, '_chip_def.Rds'))
  load(paste0('code/results/out_sample_', rep, '_lb_low_omeg.Rds'))
  load(paste0('code/results/out_sample_', rep, '_lb_high_gam.Rds'))
  load(paste0('code/results/out_sample_', rep, '_lb_def.Rds'))
  
  Y_pred.chip.def = do.call(rbind, Y_pred.chip.def)
  Y_pred.chip.low.alpha = do.call(rbind, Y_pred.chip.low.alpha)
  Y_pred.chip.low.beta = do.call(rbind, Y_pred.chip.low.beta)
  Y_pred.lb.def = do.call(rbind, Y_pred.lb.def)
  Y_pred.lb.high.gam = do.call(rbind, Y_pred.lb.high.gam)
  Y_pred.lb.low.omeg = do.call(rbind, Y_pred.lb.low.omeg)
  
  mse.chip.def = mean((colMeans(Y_pred.chip.def) - Y_list[[rep]])^2)
  mse.chip.low.alpha =mean((colMeans(Y_pred.chip.low.alpha) - Y_list[[rep]])^2)
  mse.chip.low.beta =mean((colMeans(Y_pred.chip.low.beta) - Y_list[[rep]])^2)
  mse.lb.def =mean((colMeans(Y_pred.lb.def) - Y_list[[rep]])^2)
  mse.lb.high.gam =mean((colMeans(Y_pred.lb.high.gam) - Y_list[[rep]])^2)
  mse.lb.low.omeg =mean((colMeans(Y_pred.lb.low.omeg) - Y_list[[rep]])^2)
  mse.list[[rep]] <- c(mse.chip.def, mse.chip.low.alpha, mse.chip.low.beta,
                       mse.lb.def, mse.lb.high.gam, mse.lb.low.omeg)
}


mse.list[[1]]
mse.df <- do.call(cbind, mse.list)
mse.df <- data.frame(mse.df)
mse.df$models <- c('CL def', 'CL low alpha', 'CL low beta',
                   'LB def', 'LB high gam', 'LB low omeg')

mse.df.list <- list()
for(i in 1:nrow(mse.df)){
  mse.name <- mse.df$models[i]
  mse.values <- as.numeric(mse.df[i, 1:10])
  mse.part.df <-data.frame(mse = mse.values,
                      model = rep(mse.name, 10))
  mse.df.list[[i]] <- mse.part.df
}

mse.to.plot <- do.call(rbind, mse.df.list)


model_labels <- c(
  expression(atop('CL', scriptstyle(paste(alpha, "=", 0.95, ", ", beta, '=', 2)))),
  expression(atop('CL', scriptstyle(paste(alpha, "=", 0.25, ", ", beta, "=", 2)))),
  expression(atop('CL', scriptstyle(paste(alpha, "=", 0.25, ", ", beta, "=", 0.5)))),
  expression(atop('LB', scriptstyle(paste(omega, "=", 1.56, ", ", gamma, "=", 0.62)))),
  expression(atop('LB', scriptstyle(paste(omega, "=", 0.5, ", ", gamma, "=", 1.5)))),
  expression(atop('LB', scriptstyle(paste(omega, "=", 0.5, ", ", gamma, "=", 0.62))))
)


# in sample predictions
# LB prior

st = Sys.time()
load('code/results/sim_ex_lossb_default_500chain_500iter_inclFALSE.Rds')
tree_Set = MCMC_lossb.prior.multi$trees[idx.burn_in]
Y_pred.lb.def = get_value_tree_posterior(tree_set = tree_Set, 
                                         X = X, n.cores = 10)
print(Sys.time() - st)
file_path = paste0('code/results/in_sample_lb_def.Rds')
rm(MCMC_lossb.prior.multi)
save(Y_pred.lb.def, file = file_path)

###############
st = Sys.time()
load('code/results/sim_ex_lossb_highgam_500chain_500iter_inclFALSE.Rds')
tree_Set = MCMC_lossb.prior.high.gam$trees[idx.burn_in]
Y_pred.lb.high.gam = get_value_tree_posterior(tree_set = tree_Set, 
                                              X = X, n.cores = 10)
print(Sys.time() - st)
rm(MCMC_lossb.prior.high.gam)
file_path = paste0('code/results/in_sample_lb_high_gam.Rds')
save(Y_pred.lb.high.gam, file = file_path)

##########
st = Sys.time()
load('code/results/sim_ex_lossb_lowomeg_500chain_500iter_inclFALSE.Rds')
tree_Set = MCMC_lossb.prior.low.omeg$trees[idx.burn_in]
Y_pred.lb.low.omeg = get_value_tree_posterior(tree_set = tree_Set, 
                                              X = X, n.cores = 10)
print(Sys.time() - st)
rm(MCMC_lossb.prior.low.omeg)
file_path = paste0('code/results/in_sample_lb_low_omeg.Rds')
save(Y_pred.lb.low.omeg, file = file_path)


# CHIPMAN PRIOR
#############
st = Sys.time()
load('code/results/sim_ex_chip_def_500chain_500iter_inclFALSE.Rds')
tree_Set = MCMC_chip.def$trees[idx.burn_in]
Y_pred.chip.def = get_value_tree_posterior(tree_set = tree_Set, 
                                           X = X, n.cores = 10)

print(Sys.time() - st)
rm(MCMC_chip.def)
file_path = paste0('code/results/in_sample_chip_def.Rds')
save(Y_pred.chip.def, file = file_path)

###############
st = Sys.time()
load('code/results/sim_ex_chip_lowalpha_500chain_500iter_inclFALSE.Rds')
tree_Set = MCMC_chip.low.alpha$trees[idx.burn_in]
Y_pred.chip.low.alpha = get_value_tree_posterior(tree_set = tree_Set, 
                                                 X = X, n.cores = 10)

print(Sys.time() - st)
rm(MCMC_chip.low.alpha)
file_path = paste0('code/results/in_sample_chip_low_alpha.Rds')
save(Y_pred.chip.low.alpha, file = file_path)

##################
st = Sys.time()
load('code/results/sim_ex_chip_lowbeta_500chain_500iter_inclFALSE.Rds')
tree_Set = MCMC_chip.low.beta$trees[idx.burn_in]
Y_pred.chip.low.beta = get_value_tree_posterior(tree_set = tree_Set, 
                                                X = X, n.cores = 10)

print(Sys.time() - st)
rm(MCMC_chip.low.beta)
file_path = paste0('code/results/in_sample_low_beta.Rds')
save(Y_pred.chip.low.beta, file = file_path)


Y_pred.chip.def = do.call(rbind, Y_pred.chip.def)
Y_pred.chip.low.alpha = do.call(rbind, Y_pred.chip.low.alpha)
Y_pred.chip.low.beta = do.call(rbind, Y_pred.chip.low.beta)
Y_pred.lb.def = do.call(rbind, Y_pred.lb.def)
Y_pred.lb.high.gam = do.call(rbind, Y_pred.lb.high.gam)
Y_pred.lb.low.omeg = do.call(rbind, Y_pred.lb.low.omeg)

mse.chip.def = mean((colMeans(Y_pred.chip.def) - Y_ex)^2)
mse.chip.low.alpha =mean((colMeans(Y_pred.chip.low.alpha) - Y_ex)^2)
mse.chip.low.beta =mean((colMeans(Y_pred.chip.low.beta) - Y_ex)^2)
mse.lb.def =mean((colMeans(Y_pred.lb.def) - Y_ex)^2)
mse.lb.high.gam =mean((colMeans(Y_pred.lb.high.gam) - Y_ex)^2)
mse.lb.low.omeg =mean((colMeans(Y_pred.lb.low.omeg) - Y_ex)^2)
insample.df <- data.frame( mse = c(mse.chip.def, 
                                   mse.chip.low.alpha, 
                                   mse.chip.low.beta,
                                   mse.lb.def, 
                                   mse.lb.high.gam, 
                                   mse.lb.low.omeg),
                           model = c('CL def', 'CL low alpha', 'CL low beta', 
                                     'LB def', 'LB high gam', 'LB low omeg'))

mse.to.plot$sample <- 'out-of-sample'
mse.to.plot$rep <- rep(1:10, 6)
insample.df$sample <- 'in-sample'
insample.df$rep <- 11


filled_shapes <- c(24, 21)#, 23, 24, 25, 21, 22, 23, 24, 25, 21) 

dim.mult = 1

##################
# FIGURE 8 (Top) #
##################

pdf(file = 'figures/mse_synth_outofsample_points.pdf', width = 5*2, height = 5)
ggplot(rbind(mse.to.plot, insample.df),
       aes(x = as.factor(model), y = mse, fill = factor(rep))) +
  geom_point(
    aes(shape = sample),
    #color = "black",
    size = 3,
    show.legend = TRUE,
    stroke = 1.1
  ) +
  guides(fill = "none") + 
  scale_shape_manual(values = filled_shapes) + #, guide = "none") +
  scale_x_discrete(labels = parse(text = model_labels)) +
  labs(x = NULL, y = "MSE", shape = "") +
 # guides(
  #  fill = guide_legend(
  #    override.aes = list(shape = 21, size = 4, color = "black")  # ✅ forces filled symbols in legend
  #  )
  #) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 7.5))
dev.off()




