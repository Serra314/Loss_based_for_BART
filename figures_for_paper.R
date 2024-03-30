source('code/import_functions.R')
library(ggplot2)
library(latex2exp)

################################################
## DEPTH DISTRIBUTION UNDER CLASSIC TREE PRIOR #
################################################

depth.sim.df <- rbind(data.frame(depth = rdepth.chip(10000, 0.5, 0.5),
                                 panel.name = c('alpha = 0.5, beta = 0.5')),
                      data.frame(depth = rdepth.chip(10000, 0.5, 1),
                                 panel.name = c('alpha = 0.5, beta = 1')),
                      data.frame(depth = rdepth.chip(10000, 0.95, 1),
                                 panel.name = c('alpha = 0.95, beta = 1')),
                      data.frame(depth = rdepth.chip(10000, 0.95, 2),
                                 panel.name = c('alpha = 0.95, beta = 2')))

depth.sim.df$panel.name <- factor(depth.sim.df$panel.name,
                                  labels=c('alpha = 0.5, beta = 0.5'=parse(text=TeX('$\\alpha = 0.5, \\beta = 0.5$')),
                                           'alpha = 0.5, beta = 1'=parse(text=TeX('$\\alpha = 0.5, \\beta = 1$')),
                                           'alpha = 0.95, beta = 1'=parse(text=TeX('$\\alpha = 0.95, \\beta = 1$')),
                                           'alpha = 0.95, beta = 2'=parse(text=TeX('$\\alpha = 0.95, \\beta = 2$'))))
pdf(file = 'draft/figures/fig_depth_chipman.pdf', width = 6, height = 3)
ggplot(depth.sim.df, aes(x = depth)) + 
  geom_vline(xintercept = c(0, 4, 8, 12), color = 'grey', alpha = 0.5) +
  geom_hline(yintercept = c(0.2, 0.4), color = 'grey', alpha = 0.5) + 
  geom_histogram(aes(y = after_stat(density)), 
                 color = 'black', fill = 'white', binwidth = 1)+
  theme_classic() + 
  xlab('Depth') +
  ylab('PMF') + 
  facet_wrap(facets = ~ panel.name, labeller = label_parsed)
dev.off()
################################################
## Nterm DISTRIBUTION UNDER CLASSIC TREE PRIOR #
################################################

nterm.sim.df <- rbind(data.frame(nterm = rnterm.chip(10000, 0.5, 0.5),
                                 panel.name = c('alpha = 0.5, beta = 0.5')),
                      data.frame(nterm = rnterm.chip(10000, 0.5, 1),
                                 panel.name = c('alpha = 0.5, beta = 1')),
                      data.frame(nterm = rnterm.chip(10000, 0.95, 1),
                                 panel.name = c('alpha = 0.95, beta = 1')),
                      data.frame(nterm = rnterm.chip(10000, 0.95, 2),
                                 panel.name = c('alpha = 0.95, beta = 2')))

nterm.sim.df$panel.name <- factor(nterm.sim.df$panel.name,
                                  labels=c('alpha = 0.5, beta = 0.5'=parse(text=TeX('$\\alpha = 0.5, \\beta = 0.5$')),
                                           'alpha = 0.5, beta = 1'=parse(text=TeX('$\\alpha = 0.5, \\beta = 1$')),
                                           'alpha = 0.95, beta = 1'=parse(text=TeX('$\\alpha = 0.95, \\beta = 1$')),
                                           'alpha = 0.95, beta = 2'=parse(text=TeX('$\\alpha = 0.95, \\beta = 2$'))))
pdf(file = 'draft/figures/fig_nterm_chipman.pdf', width = 6, height = 3)
ggplot(nterm.sim.df, aes(x = nterm)) + 
  geom_vline(xintercept = c(5,10,15,20), color = 'grey', alpha = 0.5) +
  geom_hline(yintercept = c(0.2, 0.4), color = 'grey', alpha = 0.5) + 
  geom_histogram(aes(y = after_stat(density)), 
                 color = 'black', fill = 'white', binwidth = 1)+
  theme_classic() + 
  xlab('Number of terminal nodes') +
  ylab('PMF') + 
  facet_wrap(facets = ~ panel.name, labeller = label_parsed)
dev.off()

############################################################
## LOSS-BASED PRIOR DEPTH DISTRIBUTION VARYING PARAMETERS ##
############################################################

depth.alt.1.1 <- vapply(1:10000, \(x) get_depth(generate_random_binary_tree_from_prior(1, 2)$tree), 0)
depth.alt.01.01 <- vapply(1:10000, \(x) get_depth(generate_random_binary_tree_from_prior(0.1, 0.01)$tree), 0)
depth.alt.01.1 <- vapply(1:10000, \(x) get_depth(generate_random_binary_tree_from_prior(0.1, 2)$tree), 0)
depth.alt.1.01 <- vapply(1:10000, \(x) get_depth(generate_random_binary_tree_from_prior(1, 0.01)$tree), 0)

max(depth.alt.01.01)
max(depth.alt.01.1)

depth.alt.1.1 <- sample(1:10, 10000, replace = TRUE)
depth.alt.01.01 <- sample(1:10, 10000, replace = TRUE)
depth.alt.01.1 <-sample(1:10, 10000, replace = TRUE)
depth.alt.1.01 <- sample(1:10, 10000, replace = TRUE)





depth.alt.toplot <- rbind(data.frame(depth = depth.alt.01.01,
                                     panel.name = 'o = 0.1, g = 0.1'),
                          data.frame(depth = depth.alt.01.1,
                                     panel.name = 'o = 0.1, g = 1'),
                          data.frame(depth = depth.alt.1.01,
                                     panel.name = 'o = 1, g = 0.1'),
                          data.frame(depth = depth.alt.1.1,
                                     panel.name = 'o = 1, g = 1'))

depth.alt.toplot$panel.name <- factor(depth.alt.toplot$panel.name,
                                  labels=c('o = 0.1, g = 0.1'=parse(text=TeX('$\\omega = 0.1, \\gamma = 0.01$')),
                                           'o = 0.1, g = 1'=parse(text=TeX('$\\omega = 0.1, \\gamma = 2$')),
                                           'o = 1, g = 0.1'=parse(text=TeX('$\\omega = 1, \\gamma = 0.01$')),
                                           'o = 1, g = 1'=parse(text=TeX('$\\omega = 1, \\gamma = 2$'))))

pdf(file = 'draft/figures/fig_depth_pmf_alt.pdf', width = 6, height = 3)
ggplot(depth.alt.toplot, aes(x = depth)) + 
  geom_vline(xintercept = c(0,5,10,15,20,25), color = 'grey', alpha = 0.5) +
  geom_hline(yintercept = c(0.2, 0.4, 0.6), color = 'grey', alpha = 0.5) + 
  geom_histogram(aes(y = after_stat(density)), 
                 color = 'black', fill = 'white', binwidth = 1)+
  theme_classic() + 
  xlab('Depth') +
  ylab('PMF') + 
  facet_wrap(facets = ~ panel.name, labeller = label_parsed)
dev.off()


table.depth.alt.1.1 <- table(depth.alt.1.1)
table.depth.alt.01.01 <- table(depth.alt.01.01)
table.depth.alt.01.1 <- table(depth.alt.01.1)
table.depth.alt.1.01 <- table(depth.alt.1.01)


ecdf.depth.alt.toplot <- rbind(data.frame(depth = as.numeric(names(table.depth.alt.01.01)),
                                          ecdf = cumsum(table.depth.alt.01.01),
                                          panel.name = 'o = 0.1, g = 0.1'),
                               data.frame(depth = as.numeric(names(table.depth.alt.1.1)),
                                          ecdf = cumsum(table.depth.alt.1.1),
                                          panel.name = 'o = 1, g = 1'),
                               data.frame(depth = as.numeric(names(table.depth.alt.01.1)),
                                          ecdf = cumsum(table.depth.alt.01.1),
                                          panel.name = 'o = 0.1, g = 1'),
                               data.frame(depth = as.numeric(names(table.depth.alt.1.01)),
                                          ecdf = cumsum(table.depth.alt.1.01),
                                          panel.name = 'o = 1, g = 0.1'))

ecdf.depth.alt.toplot$panel.name <- factor(ecdf.depth.alt.toplot$panel.name,
                                      labels=c('o = 0.1, g = 0.1'= TeX('$\\omega = 0.1, \\gamma = 0.01$'),
                                               'o = 0.1, g = 1'=parse(text=TeX('$\\omega = 0.1, \\gamma = 2$')),
                                               'o = 1, g = 0.1'=parse(text=TeX('$\\omega = 1, \\gamma = 0.01$')),
                                               'o = 1, g = 1'=parse(text=TeX('$\\omega = 1, \\gamma = 2$'))))

pdf(file = 'draft/figures/fig_depth_ecdf_alt.pdf', width = 6, height = 3)
ggplot(ecdf.depth.alt.toplot, aes(depth, ecdf/10000, color = panel.name, linetype = panel.name)) + 
  geom_line() + 
  scale_linetype_discrete(labels=c('o = 0.1, g = 0.1'= TeX('$\\omega = 0.1, \\gamma = 0.01$'),
                                   'o = 0.1, g = 1'=TeX('$\\omega = 0.1, \\gamma = 2$'),
                                   'o = 1, g = 0.1'=TeX('$\\omega = 1, \\gamma = 0.01$'),
                                   'o = 1, g = 1'=TeX('$\\omega = 1, \\gamma = 2$'))) + 
  scale_color_discrete(
                       labels=c('o = 0.1, g = 0.1'= TeX('$\\omega = 0.1, \\gamma = 0.01$'),
                                'o = 0.1, g = 1'=TeX('$\\omega = 0.1, \\gamma = 2$'),
                                'o = 1, g = 0.1'=TeX('$\\omega = 1, \\gamma = 0.01$'),
                                'o = 1, g = 1'=TeX('$\\omega = 1, \\gamma = 2$'))) + 
  ylab('ECDF') + 
  xlab('Depth') + 
  theme_classic() +
  theme(legend.text.align = 0,
        legend.title = element_blank(),
        legend.position = c(0.75, 0.5),
        legend.background = element_rect(colour = 'grey')) 
dev.off()  

################################
#### PARAMETERS CALIBRATION ####
################################

## create expected loss function

expected.nterm <- function(omeg){
  1/(1-exp(-omeg))
}

cond.expected.delta <- function(gamm, nl){
  if(nl == 1){
    delta.v = 0
    exp.value = 0
    return(exp.value)
  } else {
    if(f.odd(nl)){
      delta.v <- seq(1,nl-2,by=2)
      exp.value <- sum(delta.v*prior.delta(delta.v, gamm, nl))
      return(exp.value)
    }
    else{
      delta.v <- seq(0,(nl-2),by=2)
      exp.value <- sum(delta.v*prior.delta(delta.v, gamm, nl))
      return(exp.value)
    }  
  }
}



expected.delta <- function(omeg, gamm, thr = 1e-5){
  nl.max <- qgeom(p = (1-thr), prob = (1 - exp(-omeg)))
  nl.v <- 1:nl.max
  prob.nl <- prior.nterm(nl.v, omeg)
  cond.delta <- vapply(nl.v, function(x) cond.expected.delta(gamm, x), 0)
  return(sum(cond.delta*prob.nl))
}




expected.loss.1 <- function(omeg, gamm){
  (omeg^2)*(expected.nterm(omeg) - 1) + (gamm)*(expected.delta(omeg,gamm) - expected.delta(omeg,4))
}

expected.loss.2 <- function(omeg, gamm){
  (omeg^2)*(expected.nterm(omeg) - 1) + (gamm*omeg)*(expected.delta(omeg,gamm) - expected.delta(omeg,4))
}


df.par.1 <- expand.grid(omega = seq(1, 2, length.out = 30),
                        gamma = seq(0.1, 4, length.out = 30))
df.par.1$exp.loss <- vapply(1:nrow(df.par.1), \(x) expected.loss.1(df.par.1$omega[x], df.par.1$gamma[x]),0) 
df.par.1$panel.name <- 'loss1'

df.par.2 <- expand.grid(omega = seq(1, 2, length.out = 30),
                        gamma = seq(0.1, 4, length.out = 30))
df.par.2$exp.loss <- vapply(1:nrow(df.par.2), \(x) expected.loss.2(df.par.2$omega[x], df.par.2$gamma[x]),0) 
df.par.2$panel.name <- 'loss2'

df.par <- rbind(df.par.1, df.par.2)
df.par$panel.name <- factor(df.par$panel.name, 
                            labels = c('loss1' = TeX('$E_L^{(1)}$'),
                                       'loss2' = TeX('$E_L^{(2)}$')))


expected.loss1.to.optim <- function(par){
  -expected.loss.1(par[1], par[2])
}

expected.loss2.to.optim <- function(par){
  -expected.loss.2(par[1], par[2])
}


pp.1 <- optim(c(1,1), expected.loss1.to.optim)$par
pp.2 <- optim(c(1,1), expected.loss2.to.optim)$par

df.pp <- data.frame(x = c(pp.1[1], pp.2[1]),
                    y = c(pp.1[2], pp.2[2]),
                    panel.name = c('loss1', 'loss2'))

df.pp$panel.name <- factor(df.pp$panel.name, 
                           labels = c('loss1' = TeX('$E_L^{(1)}$'),
                                      'loss2' = TeX('$E_L^{(2)}$')))
pdf(file = 'draft/figures/fig_exp_loss.pdf', width = 6, height = 3)
ggplot(df.par, aes(x = omega, y = gamma, fill = exp.loss, z = exp.loss)) + 
  geom_tile() + 
  geom_contour() + 
  scale_fill_viridis_c() + 
  geom_point(data = df.pp, mapping = aes(x,y,fill=NULL, z = NULL), size = 0.5) + 
  facet_wrap(facets = ~ panel.name, labeller = label_parsed) + 
  xlab(~omega) + 
  ylab(~gamma) +
  theme_classic() + 
  theme(legend.position = 'none')
dev.off()

