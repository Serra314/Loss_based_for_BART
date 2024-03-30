#' Loss-based marginal prior for the depth of a binary tree
#'
#' The function returns the value of the loss-based marginal prior for the depth of the tree.
#' The distribution depends on the parameter \eqn{\omega} and its given by
#' \deqn{\pi(d | \omega) = \frac{(e^{- d \omega}(e^{\omega} - 1) }{e^{\omega}}}  
#' for \eqn{\omega > 0} and \eqn{d = 0,1,2,\dots}.
#' 
#' @param d.T a vector of integers representing the depths at which the prior is evaluated
#' @param omega. a scalar representing the value of the parameter 
#'
#' @return a vector with length equal to the length of `d.T`. Each element
#' represents the value of the prior for the corresponding depth value in `d.T`
#' @export
#'
#' @examples 
#' pi.depth <- depth_loss_based(0:18, 1)
#' plot(pi.depth)
#' sum(pi.depth)
depth_loss_based <- function(d.T, omega.){
  exp(-omega.*d.T)*(exp(omega.) - 1)/exp(omega.)
}

#' Loss-based conditional prior for the number of leaves in a binary tree
#' 
#' The function calculate the value of the loss-based conditional prior for the number of leaves
#' of the tree given the depth of the tree. The distribution depends on the parameter \eqn{\epsilon}
#' and is given by
#' \deqn{\pi(l | d, \epsilon) = \frac{e^{-l \epsilon}(e^{\epsilon} - 1)}{e^{-d \epsilon } - e^{-2^{d}\epsilon }} }
#' for \eqn{\epsilon > 0}, \eqn{d = 0,1,\dots}, and \eqn{l = d+1,\dots,2^d}.
#' 
#' @param l.T a vector of integers representing the number of leaves at which the prior is evaluated
#' @param d.T a scalar representing the depth  
#' @param epsilon. a scalar representing the parameter of the distribution
#'
#' @return a vector with length equal to the length of `l.T`. Each element represents the value of the 
#' prior for the corresponding value of the `l.T`
#' @export
#'
#' @examples 
#' d = 4
#' l.v = (d+1):2^d
#' pi.l = nterm_loss_based_cond(l.v, d, 1)
#' plot(l.v, pi.l)
#' sum(pi.l)
nterm_loss_based_cond <- function(l.T, d.T, epsilon.){
  out <- rep(0, length(l.T))
  idx. <- l.T >= d.T+1 & l.T <= 2^d.T
  if(sum(idx.) != 0){
    out[idx.] <- 
      exp(-epsilon.*l.T[idx.])*(exp(epsilon.) - 1)/(exp(-epsilon.*d.T) - exp(-epsilon.*(2^d.T)))  
  }
  return(out)
}

#' Loss-based joint prior for the depth and the number of leaves of a binary tree
#' 
#' The function returns the value of the loss-based joint prior for the depth and the number of leaves of 
#' a binary tree. The distribution depends on two parameters \eqn{\omega, \epsilon > 0}, and is given by
#' \deqn{\pi(d, l | \omega, \epsilon) = \pi(d | \omega)\pi(l| d, \epsilon)} 
#' where 
#' \deqn{\pi(d | \omega) = \frac{e^{- d \omega}(e^{\omega} - 1) }{e^{\omega}}} 
#' and 
#' \deqn{\pi(l | d, \epsilon) = \frac{e^{-l \epsilon}(e^{\epsilon} - 1)}{e^{-d \epsilon } - e^{-2^{d}\epsilon }} }
#' 
#' @param d.T a vector of depth values
#' @param l.T a vector of number of leaves values
#' @param omega. parameter of the marginal distribution of the depth
#' @param epsilon. parameter of the conditional distribution of the number of leaves
#'
#' @return a vector representing the value of the joint prior distribution.
#' @export
#'
#' @examples 
#' d.l.df <- rbind(data.frame(d = 3, l = 4:9), data.frame(d = 4, l = 5:16) )
#' joint.pi <- joint_loss_based(d.l.df$d, d.l.df$l, 1, 1)
#' 
joint_loss_based <- function(d.T, l.T, omega., epsilon.){
  prior.depth <- depth_loss_based(d.T, omega.)
  prior.nterm <- vapply(1:length(l.T), \(x) nterm_loss_based_cond(l.T[x], d.T[x], epsilon.),0)
  prior.depth*prior.nterm
}

#' Loss-based joint prior for the depth of the number of leaves of a tree
#'
#' The function returns the value of the loss-based joint prior calculated on the
#' tree given as input. Given a tree with depth \eqn{d} and \eqn{l} number of 
#' leaves, the distribution depends on two parameters \eqn{\omega, \epsilon > 0} and
#' is given by 
#' \deqn{\pi(d, l | \omega, \epsilon) = \pi(d | \omega)\pi(l| d, \epsilon)} 
#' where 
#' \deqn{\pi(d | \omega) = \frac{(e^{- d \omega}(e^{\omega} - 1) }{e^{\omega}}} 
#' and 
#' \deqn{\pi(l | d, \epsilon) = \frac{e^{-l \epsilon}(e^{\epsilon} - 1)}{e^{-d \epsilon } - e^{-2^{d}\epsilon }} }

#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @param omega. scalar representing the parameter of the marginal depth distribution
#' @param epsilon. scalar representing the parameter of the conditional number of leaves distribution
#' given the depth of the tree
#'
#' @return a scalar representing the value of the joint prior for the tree provided as input
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' joint_loss_based_tree(tree_top, 1, 1)
joint_loss_based_tree <- function(tree_top, omega., epsilon.){
  depth_ <- get_depth(tree_top)
  n.term_ <- get_num_terminal_nodes(tree_top)
  joint_loss_based(d.T = depth_, l.T = n.term_, 
                   omega. = omega., epsilon. = epsilon.)
}

# priors on tree
loss_based_prior_tree <- function(tree_top, omeg, eps){
  n.term <- get_num_terminal_nodes(tree_top)
  depth <- get_depth(tree_top)
  n.tree <- enumerate_trees(depth, n.term)
  depth_loss_based(depth, omeg)*nterm_loss_based_cond(n.term, depth, eps)/n.tree
}

# priors on tree with maximum depth
loss_based_prior_tree_dm <- function(tree_top, omeg, eps, dm){
  n.term <- get_num_terminal_nodes(tree_top)
  depth <- get_depth(tree_top)
  numerator <- exp(-eps*n.term - omeg*depth)*(exp(omeg) - 1)*(exp(eps) - 1)
  denominator <- exp(-omeg*dm)*(exp(omeg*(dm + 1)) - 1)*(exp(-eps*depth) - exp(-eps*(2^depth)))
  numerator/denominator
}


# priors on depth and nterm
loss_based_prior_tree2 <- function(depth, n.term, omeg, eps){
  n.tree <- enumerate_trees(depth, n.term)
  depth_loss_based(depth, omeg)*nterm_loss_based_cond(n.term, depth, eps)/n.tree
}

# priors on depth and nterm with maximum depth
loss_based_prior_tree_dm2 <- function(depth, n.term, omeg, eps, dm){
  numerator <- exp(-eps*n.term - omeg*depth)*(exp(omeg) - 1)*(exp(eps) - 1)
  denominator <- exp(-omeg*dm)*(exp(omeg*(dm + 1)) - 1)*(exp(-eps*depth) - exp(-eps*(2^depth)))
  numerator/denominator
}

left.right.nodes <- function(tree_top, int_nodes_idx){
  off.list <- lapply(int_nodes_idx, function(x) get_offsprings(x, tree_top))
  n.left.right <- lapply(off.list, function(x) data.frame(left = get_num_terminal_nodes(x$left),
                                                          right = get_num_terminal_nodes(x$right)))
  n.left.right.df <- do.call(rbind, n.left.right)
  cbind(node.idx = int_nodes_idx, n.left.right.df)
}

prob.split <- function(valid.split, cont.unif = TRUE){
  if(length(valid.split) == 1){
    return(1)
  }
  if(is.character(valid.split)){
    1/length(valid.split)
  } else{
    if(cont.unif){
      valid.split.extremes <- range(valid.split)
      1/(valid.split.extremes[2] - valid.split.extremes[1])  
    } else{
      1/length(valid.split)
    }
    
  }
}

prior.split.rule <- function(tree_top, X, cont.unif = TRUE){
  int.nodes <- get_internal_nodes_idx(tree_top)
  if(length(int.nodes) == 0){
    return(1)
  } else{
    obs.at.nodes <- lapply(int.nodes, \(x) get_obs_at_node(node.idx = x, 
                                                           X = X, 
                                                           tree_top = tree_top, 
                                                           X.orig = X))
    lr.nodes <- left.right.nodes(tree_top, int.nodes)
    pred.idx <- vapply(int.nodes, \(x) get_node_condition(x, tree_top)$x.idx,0)
    prob.pred <- vapply(1:nrow(lr.nodes), \(x) 1/length(get_set_valid_predictor(X.sel = obs.at.nodes[[x]],
                                                                                n.term.left = lr.nodes$left[x],
                                                                                n.term.right = lr.nodes$right[x],
                                                                                obs.per.term = 1)),0)
    valid.split.list <- lapply(1:nrow(lr.nodes), \(x) get_set_valid_split(X.pred = obs.at.nodes[[x]][,pred.idx[x]],
                                                                          n.term.left = lr.nodes$left[x],
                                                                          n.term.right = lr.nodes$right[x],
                                                                          obs.per.term = 1))
    prob.split <- vapply(valid.split.list, \(x) prob.split(x, cont.unif = cont.unif), 0)
    return(prod(prob.pred*prob.split))
  }
}








