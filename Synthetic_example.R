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
X <- as.data.frame(apply(X, 2, \(x) (x - min(x))/(max(x) - min(x))))
#X <- as.matrix(X)

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
ggarrange(pl.x1, pl.x2, pl.x3, ncol = 2, nrow = 2)

library(manipulateWidget)
library(plotly)

pl.x1.g <- ggplotly(p = pl.x1 + theme_classic() + ylab(TeX('Y')) + xlab(TeX('X_1')) +
                theme(axis.text.x = element_text(angle = 0, hjust = 1))) %>% 
  config(mathjax = "cdn")
pl.x2.g <- ggplotly(p = pl.x2 + theme_classic() + ylab(TeX('Y')) + xlab(TeX('X_2')) +
                      theme(axis.text.x = element_text(angle = 0, hjust = 1))) %>% 
  config(mathjax = "cdn")
pl.x3.g <- ggplotly(p = pl.x3 + theme_classic() + ylab(TeX('Y')) + xlab(TeX('X_3')) +
                      theme(axis.text.x = element_text(angle = 0, hjust = 1))) %>% 
  config(mathjax = "cdn")


# save image using export
combineWidgets(tree.plot, pl.x1.g, pl.x2.g, pl.x3.g)

## expected number of terminal nodes and depth for different priors
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
pp.2 <- c(1.5618883, 0.6293944)
X <- as.data.frame(apply(X, 2, \(x) (x - min(x))/(max(x) - min(x))))
X <- as.matrix(X)
rownames(X) <- 1:nrow(X)
# loss-based prior with default params
lossb.prior.list <- list(fun = joint.prior.new.tree, param = pp.2)
n.chain_par <- 100
n.iter_par <- 500
incl.split_par <- FALSE
cont.unif_par <- TRUE
moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)

st = Sys.time()
MCMC_lossb.prior.multi <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                                  n.iter = n.iter_par, 
                                                  X = X, 
                                                  Y = Y_ex, 
                                                  Y.var = Y.var, 
                                                  mu.prior.mean = 0, 
                                                  mu.prior.var = var(Y_ex), 
                                                  prior_list = lossb.prior.list, 
                                                  moves.prob = NULL, starting.tree = NULL,
                                                  include.split = incl.split_par,
                                                  cont.unif = cont.unif_par)
Sys.time() -st
#  1.23 hours
save(MCMC_lossb.prior.multi, file = 'code/results/sim_ex_lossb_default_100chain_500iter_inclFALSE.Rds')
rm(MCMC_lossb.prior.multi)
# chipman - default
chip.def.list <- list(fun = chipman_prior_tree, param = c(0.95, 2))

st = Sys.time()
MCMC_chip.def <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                           n.iter = n.iter_par, 
                                           X = X, 
                                           Y = Y_ex, 
                                           Y.var = Y.var, 
                                           mu.prior.mean = 0, 
                                           mu.prior.var = var(Y_ex), 
                                           prior_list = chip.def.list, 
                                           moves.prob = NULL, starting.tree = NULL,
                                           include.split = incl.split_par,
                                           cont.unif = cont.unif_par)
Sys.time() -st
# 1h and 49 mins
save(MCMC_chip.def, file = 'code/results/sim_ex_chip_def_100chain_500iter_inclFALSE.Rds')
rm(MCMC_chip.def)

model.list <- list(MCMC_lossb.prior.multi, MCMC_chip.def)
names(model.list) <- c('loss based def', 'chipman def')
# extract information from trees
df.nl <- apply_fun_models(get_num_terminal_nodes, model.list, born.out.pc = 0)
df.depth <- apply_fun_models(get_depth, model.list, born.out.pc = 0)

df.nl.fortrace <- df.nl[c(1:5000,50001:55000),]
df.depth.fortrace <- df.depth[c(1:5000,50001:55000),]
# traceplots
trace.nl <- ggplot(df.nl.fortrace) + 
  geom_vline(xintercept = seq(0,5000,by = 500), 
             color = 'grey', size = 0.2, alpha=0.75)+
  geom_line(aes(x, y)) + 
  geom_hline(yintercept = get_num_terminal_nodes(tree_ex), 
             color = 'darkorange') + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  xlab('Iteration') + 
  ylab(~n[L]) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,5000,by = 500)) + 
  theme(axis.text.x = element_text(angle = 30))


trace.depth <- ggplot(df.depth.fortrace) + 
  geom_vline(xintercept = seq(0,5000,by = 500),
             color = 'grey', size = 0.2, alpha=0.75)+
  geom_line(aes(x, y)) +
  geom_hline(yintercept = get_num_terminal_nodes(tree_ex), 
             color = 'darkorange',
             linetype = 2) + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  xlab('Iteration') + 
  ylab('Depth') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,5000,by = 500)) + 
  theme(axis.text.x = element_text(angle = 30))


pdf(file = 'draft/figures/fig_trace_plot_comp_def.pdf', width = 7, height = 4)
ggarrange(plotlist = list(trace.nl, trace.depth), ncol = 1)
dev.off()



# histograms - use bornout 
df.nl <- apply_fun_models(get_num_terminal_nodes, model.list, born.out.pc = 250)
df.depth <- apply_fun_models(get_depth, model.list, born.out.pc = 250)

hist.nl <- ggplot(df.nl) + 
  geom_histogram(aes(x = y, y = after_stat(density)), 
                 binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = get_num_terminal_nodes(tree_ex), 
             color = 'darkorange',
             linetype = 2) + 
  facet_wrap(facets = ~panel.name) + 
  xlab(~n[L]) + 
  ylab('PMF') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,24,by = 3))  

hist.depth <- ggplot(df.depth) + 
  geom_histogram(aes(x = y, y = after_stat(density)), 
                 binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = get_depth(tree_ex), 
             color = 'darkorange',
             linetype = 2) + 
  facet_wrap(facets = ~panel.name) + 
  xlab('Depth') + 
  ylab('PMF') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,10,by = 1))


pdf(file = 'draft/figures/fig_hist_comp_def.pdf', width = 7, height = 4)
ggarrange(plotlist = list(hist.nl, hist.depth), ncol = 1)
dev.off()





#################################
## Different parametrisations  ##
#################################

# lower omega -- so similar to chipman default
pp.low.omeg <- c(0.5, pp.2[2])
lossb.prior.list.low.omeg <- list(fun = joint.prior.new.tree, param = pp.low.omeg)
st = Sys.time()
MCMC_lossb.prior.low.omeg <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                                    n.iter = n.iter_par, 
                                                    X = X, 
                                                    Y = Y_ex, 
                                                    Y.var = Y.var, 
                                                    mu.prior.mean = 0, 
                                                    mu.prior.var = var(Y_ex), 
                                                    prior_list = lossb.prior.list.low.omeg, 
                                                    moves.prob = NULL, starting.tree = NULL,
                                                    include.split = incl.split_par,
                                                    cont.unif = cont.unif_par)
Sys.time() -st
save(MCMC_lossb.prior.low.omeg, file = 'code/results/sim_ex_lossb_lowomeg_100chain_500iter_inclFALSE.Rds')
rm(MCMC_lossb.prior.low.omeg)

# low omega high gamma
pp.low.omeg.high.gam <- c(0.5, 1.5)
lossb.prior.list.high.gam <- list(fun = joint.prior.new.tree, param = pp.low.omeg.high.gam)
st = Sys.time()
MCMC_lossb.prior.high.gam <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                                       n.iter = n.iter_par, 
                                                       X = X, 
                                                       Y = Y_ex, 
                                                       Y.var = Y.var, 
                                                       mu.prior.mean = 0, 
                                                       mu.prior.var = var(Y_ex), 
                                                       prior_list = lossb.prior.list.high.gam, 
                                                       moves.prob = NULL, starting.tree = NULL,
                                                       include.split = incl.split_par,
                                                       cont.unif = cont.unif_par)
Sys.time() -st
save(MCMC_lossb.prior.high.gam, file = 'code/results/sim_ex_lossb_highgam_100chain_500iter_inclFALSE.Rds')
rm(MCMC_lossb.prior.high.gam)


# chipman lower alpha to match loss-based defult

mean(rdepth.chip(10000, 0.25, 2))
mean(rdepth(10000, pp.2[1], pp.2[2]))
mean(rnterm.chip(10000, 0.25, 2))
mean(rnterm(10000, pp.2[1]))

chip.prior.low.alpha <- list(fun = chipman_prior_tree, param = c(0.25, 2))

st = Sys.time()
MCMC_chip.low.alpha <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                                 n.iter = n.iter_par, 
                                                 X = X, 
                                                 Y = Y_ex, 
                                                 Y.var = Y.var, 
                                                 mu.prior.mean = 0, 
                                                 mu.prior.var = var(Y_ex), 
                                                 prior_list = chip.prior.low.alpha, 
                                                 moves.prob = NULL, starting.tree = NULL,
                                                 include.split = incl.split_par,
                                                 cont.unif = cont.unif_par)
Sys.time() -st
# 1h and 49 mins
save(MCMC_chip.low.alpha, file = 'code/results/sim_ex_chip_lowalpha_100chain_500iter_inclFALSE.Rds')
rm(MCMC_chip.low.alpha)
# low beta just to change
chip.prior.low.beta <- list(fun = chipman_prior_tree, param = c(0.25, 0.5))

st = Sys.time()
MCMC_chip.low.beta <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                                 n.iter = n.iter_par, 
                                                 X = X, 
                                                 Y = Y_ex, 
                                                 Y.var = Y.var, 
                                                 mu.prior.mean = 0, 
                                                 mu.prior.var = var(Y_ex), 
                                                 prior_list = chip.prior.low.beta, 
                                                 moves.prob = NULL, starting.tree = NULL,
                                                include.split = incl.split_par,
                                                cont.unif = cont.unif_par)
Sys.time() -st
# 10 mins
save(MCMC_chip.low.beta, file = 'code/results/sim_ex_chip_lowbeta_100chain_500iter_inclFALSE.Rds')
rm(MCMC_chip.low.beta)

# prob.split
# prior.split.rule



############################################
# IMPORT MODELS RESULT AND ANALYSE RESULTS #
############################################

load('code/results/sim_ex_chip_def_100chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_chip_lowbeta_100chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_chip_lowalpha_100chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_lossb_default_100chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_lossb_lowomeg_100chain_500iter_inclFALSE.Rds')
load('code/results/sim_ex_lossb_highgam_100chain_500iter_inclFALSE.Rds')



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

ggplot(df.dd, aes(x = y, y = after_stat(density))) + 
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') +
  facet_wrap(facets = ~ panel.name, ncol = 3) + 
  geom_vline(xintercept = get_depth(tree_ex)) +
  xlab('depth')

ggplot(df.nl, aes(x = y, y = after_stat(density))) + 
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') +
  facet_wrap(facets = ~ panel.name, ncol = 3) + 
  geom_vline(xintercept = get_num_terminal_nodes(tree_ex)) + 
  xlab('nterm')


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

pdf(file = 'draft/figures/fig_posterior_nl.pdf', width = 6, height = 3)
ggplot(df.nl, aes(x = y, y = after_stat(density))) +
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') + 
  xlab('Number terminal nodes') + 
  geom_line(data = df.nl.prior, aes(x,y), colour = 'deepskyblue', linewidth = 0.4) +
  geom_point(data = df.nl.prior, aes(x,y), colour = 'deepskyblue', size = 0.4) + 
  geom_vline(xintercept = get_num_terminal_nodes(tree_ex), color = 'orange') +
  scale_x_continuous(breaks = seq(0,21,by=3)) + 
  xlim(0,21) + 
  facet_wrap(facets = ~panel.name, nrow = 2, labeller = label_parsed) + 
  theme_classic()
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
                               pr.up5 = aggregate(df.prior.depth.samp$samp, list(df.prior.depth.samp$models), 
                                                  \(x){mean(x > 6)})[,2])

df.depth.stat.post <- data.frame(models = aggregate(df.depth$y, list(df.depth$panel.name), mean)[,1],
                              mean = aggregate(df.depth$y, list(df.depth$panel.name), mean)[,2],
                              median = aggregate(df.depth$y, list(df.depth$panel.name), median)[,2],
                              pr.up8 = aggregate(df.depth$y, list(df.depth$panel.name), \(x){mean(x >= 4)})[,2])



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
pdf(file = 'draft/figures/fig_posterior_depth.pdf', width = 6, height = 3)
ggplot(df.depth, aes(x = y, y = after_stat(density))) +
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') + 
  xlab('Depth') +
  geom_line(data = df.depth.prior, aes(x,y), colour = 'deepskyblue', linewidth = 0.4) +
  geom_point(data = df.depth.prior, aes(x,y), colour = 'deepskyblue', size = 0.4) + 
  geom_vline(xintercept = get_depth(tree_ex), color = 'orange') +
  scale_x_continuous(breaks = seq(0,20,by=2)) + 
  facet_wrap(facets = ~panel.name, nrow = 2, labeller = label_parsed) + 
  theme_classic()
dev.off()

######
# use after_stat(ncount) to get the histograms normalised to have max at 1



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
pdf(file = 'draft/figures/fig_trace_nl.pdf', width = 6, height = 6)
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
pdf(file = 'draft/figures/fig_trace_depth.pdf', width = 6, height = 6)
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

pdf(file = 'draft/figures/fig_trace_loglik.pdf', width = 6, height = 6)
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

idx = 150
get_tree_plot(MCMC_chip.prior.multi$trees[[idx]])
pp <- get_pred(MCMC_chip.prior.multi$trees[[idx]])
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

# get predictor frequency for chipman
# alpha = 0.95, beta = 2
get_pred_freq(MCMC_chip.prior.multi$trees, 1:3)
# alpha = 0.95, beta = 3
get_pred_freq(MCMC_chip.high.beta.multi$trees, 1:3)

# get predictor frequency for alternative prior 
# alpha = 1, gamma = 1
get_pred_freq(MCMC_alt.prior.multi$trees, 1:3)

# alpha = 1.9, gamma = 1
get_pred_freq(MCMC_alt.high.alpha.multi$trees, 1:3)

# alpha = 2.9, gamma = 1
get_pred_freq(MCMC_alt.very.high.alpha.multi$trees, 1:3)

# alpha = 1, gamma = 1.9
get_pred_freq(MCMC_alt.high.gamma.multi$trees, 1:3)

