
source('code/import_functions.R')
library(data.tree)
library(ggplot2)
library(ggpubr)
set.seed(3)
data.multipl <- 1
x1 <- c(runif(data.multipl*200, 0.1, 0.4), runif(data.multipl*100, 0.6, 0.9))
x2 <- c(runif(data.multipl*100, 0.1, 0.4), runif(data.multipl*100, 0.6, 0.9), runif(data.multipl*100, 0.1, 0.9))
x3 <- c(runif(data.multipl*200, 0.6, 0.9), runif(data.multipl*100, 0.1, 0.4))
# store them in data.frame
X = data.frame(x1 = x1, x2 = x2, x3 = x3)


rt <- generate_random_binary_tree_depth_free(5)
rt <- assign_node_idx(rt)
get_tree_plot.idx(rt)
rt <- assign_split_rules(rt, X = X)
term.nodes.v <- c(1, 0, 0, 1, 0)
rt <- assign_term_node_values_fixed(rt, term.nodes.v)
get_tree_plot(rt)


rbind(get_terminal_nodes_idx(rt),
      vapply(get_terminal_nodes_idx(rt), \(x) nrow(get_obs_at_node(x, X, rt)), 0))

Y_b <- sample_CART_binary(rt, X) 

data.df <- cbind(Y = Y_b, X)
# plot
pl.x1 <- ggplot(data.df, aes(x1, Y_b)) + geom_point(size = 0.5)
pl.x2 <- ggplot(data.df) + geom_point(aes(x2, Y_b), size = 0.5) + 
  xlab(~x[2])
pl.x3 <- ggplot(data.df) + geom_point(aes(x3, Y_b), size = 0.5) + 
  xlab(~x[3])
ggarrange(pl.x1, pl.x2, pl.x3, ncol = 2, nrow = 2)

lossb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.56, 0.62))
mcmc.try.def <- multichain_MCMC_binary(n.iter = 1000,
                                       n.chain = 10,
                                       X = X,
                                       Y = Y_b,
                                       alpha.prior = 1,
                                       beta.prior = 1,
                                       prior_list = lossb.prior.def,
                                       moves.prob = c(0.4, 0.4, 0.1, 0.1))
  
  

lossb.prior.low <- list(fun = joint.prior.new.tree, param = c(0.1, 0.1))
mcmc.try.low <- multichain_MCMC_binary(n.iter = 1000,
                                       n.chain = 10,
                                       X = X,
                                       Y = Y_b,
                                       alpha.prior = 1,
                                       beta.prior = 1,
                                       prior_list = lossb.prior.low,
                                       moves.prob = c(0.4, 0.4, 0.1, 0.1))



model.list <- list(mcmc.try.low, mcmc.try.def)
names(model.list) <- c('low', 'def')
# extract information from trees
df.nl <- apply_fun_models(get_num_terminal_nodes, model.list, born.out.pc = 0)
df.depth <- apply_fun_models(get_depth, model.list)
df.loglik <- apply_fun_models(\(x) cart_log_lik_binary(x, Y_b, X), model.list)



ggplot(df.nl, aes(x,y)) + 
  geom_line() + 
  facet_wrap(facets = ~panel.name)

ggplot(df.nl, aes(x = y)) + 
  geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = get_num_terminal_nodes(rt)) +
  facet_wrap(facets = ~panel.name)


ggplot(df.depth, aes(x,y)) + 
  geom_line() + 
  facet_wrap(facets = ~panel.name)


ggplot(df.depth, aes(x = y)) + 
  geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = get_depth(rt)) +
  facet_wrap(facets = ~panel.name)


ggplot(df.loglik, aes(x,y)) + 
  geom_line() + 
  facet_wrap(facets = ~panel.name)


ggplot() + 
  geom_line(aes(x = 1:1000, y = dd)) + 
  geom_hline(yintercept = get_depth(rt), color = 'orange' )

nl <- vapply(mcmc_try$trees, get_num_terminal_nodes, 1)
ggplot() + 
  geom_histogram(aes(x = nl), binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = get_num_terminal_nodes(rt)) 
  

ggplot() + 
  geom_line(aes(x = 1:1000, y = nl)) + 
  geom_hline(yintercept = get_num_terminal_nodes(rt), color = 'orange' )



get_num_obs_at_term <- function(tree_top, X){
  vapply(get_terminal_nodes_idx(tree_top), 
         function(x) nrow(get_obs_at_node(node.idx = x,X = X, tree_top = tree_top)), 0)
} 

# check number of obs per terminal nodes 
tt <- vapply(mcmc_try$trees, \(x) sum(get_num_obs_at_term(x, X) == 0), 0)
mean(tt)
cbind(1:1000, tt)

# loglikelihoood
ll <- vapply(mcmc_try$trees, function(x) cart_log_lik_binary(x, Y_b, X), 0)

ggplot() + 
  geom_line(aes(x = 1:1000, y = ll)) 



X$true.prob <- get_value_tree(rt, X)
post.prob.list <- lapply(mcmc_try$trees, 
                         \(tree) get_value_tree(tree, X)) 
post.prob.matrix <- Reduce(rbind, post.prob.list)
X$mean.prob <- apply(post.prob.matrix, 2, mean)

plot(X$true.prob, X$mean.prob)
plot(Y_b, X$mean.prob)




ll <- c()
for(i in 1:100){
  print(i)
  ll <- c(ll, cart_log_lik_binary(mcmc_try$trees[[i]], Y_b, X))
}

get_tree_plot.idx(mcmc_try$trees[[25]])

get_tree_plot.idx(mcmc_try$trees[[26]])

debugonce(cart_log_lik_binary)
cart_log_lik_binary(mcmc_try$trees[[46]], Y_b, X)

r2 <- grow_move_binary(mcmc_try$trees[[25]], X, Y_b, 1, 1)
get_tree_plot.idx(r2$tree)


acceptance.prob.list_binary(move_list = mcmc_try$trees[[26]], 
                            old.tree = mcmc_try$trees[[25]], X = X, Y = Y_b, 
                            prior_input_list = lossb.prior.def)






get_num_obs_at_term(mcmc_try$trees[[827]], X)
get_num_obs_at_term(mcmc_try$trees[[828]], X)

rr <- grow_move_binary(tree_top = mcmc_try$trees[[827]], X = X, Y = Y_b, alpha.prior = 1, beta.prior = 1)

get_tree_plot.idx(rr$tree)

get_tree_plot.idx(mcmc_try$trees[[827]])
get_tree_plot(mcmc_try$trees[[828]])
get_tree_plot(mcmc_try$trees[[829]])

ggplot() + 
  geom_histogram(aes(x = ll), binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = get_num_terminal_nodes(rt)) 







get_value_tree(tree_top, term.node.obs[[2]])






