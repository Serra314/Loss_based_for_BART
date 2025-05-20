#' Evaluate a tree 
#' 
#' The function evaluate a tree for a single observation
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param x 
#' A row of a `data.frame` representing an observation on which the tree is evaluated
#' @return a scalar representing the value at the terminal node associated with `x`
#' @export
#'
#' @examples 
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(x1 = rnorm(100), x2 = sample(c('A', 'B'), 500, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 10)
#' get_tree_plot(tree_top)
#' cbind(XX[1,], tree_value = g.T(tree_top, XX[1,]))
#' cbind(XX[10,], tree_value = g.T(tree_top, XX[10,]))

g.T <- function(tree_top, x){
  if(is.null(tree_top$cond)){
    return(tree_top$value)
  }
  else(
    if(eval_cond(tree_top$cond, x)){
      g.T(tree_top$left, x)
    } else{
      g.T(tree_top$right, x)
    }
  )
}

#' Evaluate a tree on multiple observations
#' 
#' The function evaluates the tree for each row of `X` 
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param X
#' A `data.frame` representing the observations on which the tree is evaluated
#'
#' @return a numeric vector with length equal to the number of rows of `X`. 
#' Each element represents the value of the tree for the corresponding row of `X`
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(x1 = rnorm(100), x2 = sample(c('A', 'B'), 500, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 10)
#' get_tree_plot(tree_top)
#' cbind(XX[1:5,], tree_value = get_value_tree(tree_top, XX[1:5,]))
get_value_tree <- function(tree_top, X){
  vapply(1:nrow(X), \(idx) g.T(tree_top, X[idx,]), 0)
}

#' Evaluate CART log-likelihood
#'
#' The function returns the log-likelihood of a CART model assuming a Gaussian density. 
#' The log-likelihood is given by
#' \deqn{\mathcal L = \sum_i \phi(y_i, g(x_i,T), \sigma)}
#' where \eqn{\phi(y_i, g(x_i, T), \sigma)} is the logarithm of a Gaussian density 
#' evaluated at \eqn{y_i}, with mean \eqn{g(x_i, T)} and standard deviation \eqn{\sigma}
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param Y vector of observations
#' @param X
#' A `data.frame` representing the predictors 
#' @param sigma_ 
#' Standard deviation of the normal distribution
#'
#' @return a scalar representing the value of the log-likelihood
#' @export
#'
#' @examples 
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(x1 = rnorm(100), x2 = sample(c('A', 'B', 'C'), 100, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 10)
#' get_tree_plot(tree_top)
#' YY <- sample_CART(tree_top, XX, 1)
#' cart_log_lik(tree_top, YY, XX, 1)
cart_log_lik <- function(tree_top, Y, X, sigma){
  if(length(Y) != nrow(X)){
    stop('Observations (Y) and predictors (X) have not the same length.')
  }
  tree.values <- get_value_tree(tree_top, X)
  sum(dnorm(Y, mean = tree.values, sd = sigma, log = TRUE))
}


# loglik for binary observations


cart_log_lik_binary <- function(tree_top, Y, X){
  tree.at.obs <- get_value_tree(tree_top, X)
  log.prob.obs <- Y*log(tree.at.obs) + (1 - Y)*log(1 - tree.at.obs)
  return(sum(log.prob.obs))
}




# 
# cart_log_lik_binary <- function(tree_top, Y, X){
#   if(length(Y) != nrow(X)){
#     stop('Observations (Y) and predictors (X) have not the same length.')
#   }
#   term.node.idx <- get_terminal_nodes_idx(tree_top)
#   term.node.obs <- lapply(term.node.idx, function(x) get_obs_at_node(node.idx = x,
#                                                                      X = X,
#                                                                      tree_top = tree_top))
#   # calculate value terminal nodes (prob of binomial)
#   term.node.value <- vapply(term.node.obs, \(x) get_value_tree(tree_top, x[1,]), 0)
#   # calculate number of obs at terminal nodes (size of binomial)
#   nobs.at.node <- vapply(term.node.obs, nrow, 0)
#   # calculate sum of obs at terminal nodes (obs of binomial)
#   Ysum.at.node <- vapply(term.node.obs, \(x) sum(Y[as.numeric(rownames(x))]),0)
#   # calculate loglik at each terminal node
#   loglik.to.sum <- vapply(1:length(term.node.idx), \(x) dbinom(x = Ysum.at.node[x],
#                                                                size = nobs.at.node[x],
#                                                                prob = term.node.value[x],
#                                                                log = TRUE), 0)
#   sum(loglik.to.sum)
# }
