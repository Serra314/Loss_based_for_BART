#' Random samples from a CART model
#'
#' The function returns a sample from a CART model. The model is 
#' \deqn{y = g(T,x) + \epsilon}
#' where \deqn{x} is a numerical vector representing the value of the predictors, 
#' \deqn{T} is a binary tree, and \deqn{\epsilon \sim N(0, \sigma^2)}. So that,
#' \deqn{y|x,T \sim N(g(T,x), \sigma^2)}
#' 
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param X
#' A `data.frame` representing the observations on which the tree is evaluated
#' @param sigma_ 
#' Standard deviation of the normal distribution from which the samples are obtained
#' @return
#' A numerical vector with length equal to the number of rows of `X`. 
#' Each element represents a sample from \deqn{y|x,T} where \deqn{x} is the corresponding row of `X`
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(x1 = rnorm(100), x2 = sample(c('A', 'B'), 500, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 10)
#' get_tree_plot(tree_top)
#' YY <- sample_CART(tree_top, XX, 1)
#' data.df <- cbind(XX, y = YY)
#' ggplot(data.df, aes(x1, y)) + geom_point()
#' ggplot(data.df, aes(x2, y)) + geom_point()
#' ggplot(data.df[data.df$x2 == 'B',], aes(x1, y)) + geom_point()
sample_CART <- function(tree_top, X, sigma_){
  y.mean <- get_value_tree(tree_top, X)
  rnorm(length(y.mean), mean = y.mean, sd = sigma_)
}


sample_CART_binary <- function(tree_top, X){
  y.prob <- get_value_tree(tree_top, X)
  rbinom(length(y.prob), size = 1, prob = y.prob)
}


#' Assign index value to a node of a binary tree
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch.
#' @param counter a scalar used for internal calculations, the default value is 1. It does not need to be changed.
#'
#' @return a nested list as the one provided in input as `tree_top`. 
#' The difference is that the returned `list` now have a `node.idx` argument representing the value of the index at that node. 
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
assign_node_idx <- function(tree_top, counter = 1){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    
    return(list(left = tree_top$left, right = tree_top$right,
                node.idx = counter, value = tree_top$value))
  } else{
    tree_top$node.idx = counter
    counter.l = counter + 1
    tree_left = assign_node_idx(tree_top$left, counter.l)
    n.term.left = get_num_terminal_nodes(tree_top$left)
    n.node.left = 2*n.term.left 
    counter.r = counter + n.node.left
    tree_right = assign_node_idx(tree_top$right, counter.r)
    return(list(left = tree_left, right = tree_right, 
                node.idx = counter, cond = tree_top$cond))
  }
}

#' Check if a predictor is valid to be used as splitting variable
#'
#' The functions checks if a predictor is valid to be a splitting variable. 
#' A predictor is valid if it is possible to assign a given minimum number of observations
#' on the left and right branch
#' @param X.pred a vector representing the values of the predictors
#' @param n.term.left a scalar representing the number of terminal nodes on the left branch
#' @param n.term.right a scalar representing the number of terminal nodes on the right branch
#' @param obs.per.term a scalar representing the minimum number of observations per terminal node
#'
#' @return `logical` `TRUE` if it is a valid predictor and `FALSE` otherwise.
#' For categorical variables it is tree if there is one predictor value with frequency greater or equal than `obs.per.term` times
#' `n.term.left` and such that `length(X.pred)` minus the frequency is greater or equal than `n.term.right`. 
#' For numerical variables if there is an interval such that there are at least `n.term.left` times `obs.per.term` on the left and `n.term.right` times `obs.per.term` on the right.  
#'
#' @examples 
#' X.p = rep(1,10)
#' check_pred_valid(X.p, 5, 5, 1)
#' X.p = rep('A', 10)
#' check_pred_valid(X.p, 5, 5, 1)
#' X.p = rnorm(10)
#' check_pred_valid(X.p, 5, 5, 1)
#' X.p = c('A', 'A', 'A', 'B', 'B')
#' check_pred_valid(X.p, 3, 2, 1)
check_pred_valid <- function(X.pred, n.term.left, n.term.right, obs.per.term = 1){
  n.term.sub <- n.term.left + n.term.right
  #cat('number preds', length(X.pred), '\n')
  if(length(X.pred) < obs.per.term*n.term.sub){
    #stop('Not enough observations')
    # modified because breaks when length(X.pred) = 0, don't know exactly why this happens
    return(TRUE)
  }
  if(is.numeric(X.pred)){
     X.sorted <- sort(X.pred, decreasing = FALSE)
     tau.inf <- X.sorted[n.term.left*obs.per.term]
     tau.sup <- X.sorted[length(X.sorted) - (n.term.right*obs.per.term) + 1]
     return(tau.inf < tau.sup)
   } else {
     X.freq <- table(X.pred)
     X.cond <- (X.freq >= obs.per.term*n.term.left) & (length(X.pred) - X.freq >=  obs.per.term*n.term.right) 
     return(sum(X.cond) >= 1)
   }
}

#' Get the set of valid predictors
#'
#' The function returns the column indexes of the valid predictors to be used 
#' as splitting variable in the internal nodes of a binary tree. A predictor is valid if it is possible to assign a given minimum number of observations
#' on the left and right branch
#' @param X.sel `data.frame` representing the observations available at a node
#' @param n.term.left a scalar representing the number of terminal nodes on the left branch
#' @param n.term.right a scalar representing the number of terminal nodes on the right branch
#' @param obs.per.term a scalar representing the minimum number of observations per terminal node
#'
#' @return numeric vector containing the indexes of the columns of `X.sel` corresponding
#' to valid predictors to be used as splitting variables
#' @export
#'
#' @examples
#' X.s <- data.frame(x1 = rnorm(10), x2 = 'A')
#' get_set_valid_predictor(X.s, 5, 5, 1)
get_set_valid_predictor <- function(X.sel, n.term.left, n.term.right, obs.per.term = 1){
  pred.idx = 1:ncol(X.sel)
  vec.valid <- vapply(pred.idx, \(x) check_pred_valid(X.sel[,x], 
                                                      n.term.left, n.term.right,
                                                      obs.per.term), TRUE)
  valid.pred <- pred.idx[vec.valid]
  valid.pred
}

#' Get the set of valid splitting values of a splitting variable
#'
#' The function returns the set of valid splitting values of a valid spitting variable. 
#' A predictor variable is valid if it is possible to assign a given minimum number of observations
#' on the left and right branch. A splitting value is valid if 
#' 1. It appears in the predictor
#' 2. There is a minimum number of observations greater or lower than it.
#' 
#' @param X.pred vector of observed values of the predictor
#' @param n.term.left a scalar representing the number of terminal nodes on the left branch
#' @param n.term.right a scalar representing the number of terminal nodes on the right branch
#' @param obs.per.term a scalar representing the minimum number of observations per terminal node
#'
#' @return a vector of valid elements of `X.pred` to be used a splitting values for a splitting rule.
#' @export
#'
#' @examples 
#' X.p <- runif(10)
#' get_set_valid_split(X.p, 5, 5, 1)
#' X.p <- c('A', 'A', 'A', 'B', 'B')
#' get_set_valid_split(X.p, 3, 2, 1)
#' 
#' 

get_set_valid_split <- function(X.pred, n.term.left, n.term.right, obs.per.term = 1){
  if(is.numeric(X.pred)){
    X.sorted <- sort(X.pred, decreasing = FALSE)
    tau.inf <- X.sorted[n.term.left*obs.per.term]
    tau.sup <- X.sorted[length(X.sorted) - (n.term.right*obs.per.term) + 1]
    return(sort(unique(X.pred[X.pred >= tau.inf & X.pred < tau.sup]) ) )  
  } else {
    X.freq <- table(X.pred)
    X.cond <- (X.freq >= obs.per.term*n.term.left) & (length(X.pred) - X.freq >=  obs.per.term*n.term.right) 
    return(names(X.freq)[X.cond])
  }
  
}


#' Generate a valid splitting rule for an internal node of a binary tree
#'
#' The function returns a random valid splitting rule for an internal node of a binary tree. 
#' A valid splitting rule is composed by a valid splitting variable and a valid splitting value.
#' A splitting variable is valid if it is possible to assign a given minimum number of observations
#' on the left and right branch. A splitting value is valid if 
#' 1. It appears in the splitting variable
#' 2. There is a minimum number of observations greater or lower than it.
#' 
#' @param node.idx a scalar representing the index of the node at which the condition is generated. 
#' The nodes are indexed going from left to right
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node, 
#' `cond` (only for internal nodes) splitting rules. This is needed only for the nodes before the one at which we want to generate the splitting rule. 
#' @param X `data.frame` containing all the observed values of the predictor variables 
#' @param obs.per.term a scalar representing the minimum number of observations per terminal node
#' @param cont.unif `logical`, if `TRUE` the splitting values are sampled from a continuous uniform distribution on the
#'  on the range of the valid split values. If `FALSE` a uniform distribution on the obsersed valid split values is used. Default is `TRUE`.
#'
#' @return a `list` of the form `list(x.idx = ..., x.val = ...)` where 
#' `x.idx` is the index of the predictor used as splitting variable
#' `x.val` the value of the splitting variable used as splitting value
#' @export
#'
#' @examples 
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(rnorm(100), sample(c('A', 'B'), 100, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot.idx(tree_top)
#' get_tree_plot(tree_top)
#' gen_node_condition(2, tree_top, XX)
gen_node_condition <- function(node.idx, tree_top, X, 
                               obs.per.term = 1, 
                               cont.unif = TRUE,
                               for.grow = FALSE){
  obs.at.node <- get_obs_at_node(node.idx = node.idx,
                                 X = X,
                                 tree_top = tree_top,
                                 X.orig = X)
  subtree.at.node <- get_offsprings(node.idx, tree_top)
  if(for.grow){
    n.term.left <- 1
    n.term.right <- 1
  } else {
    n.term.left <- get_num_terminal_nodes(subtree.at.node$left) 
    n.term.right <- get_num_terminal_nodes(subtree.at.node$right)
  }
  valid.pred <- get_set_valid_predictor(X.sel = obs.at.node,
                                        n.term.left = n.term.left,
                                        n.term.right = n.term.right,
                                        obs.per.term = obs.per.term)
  if(length(valid.pred) == 1){
    predictor.idx = valid.pred
  } else {
    predictor.idx <- sample(valid.pred, 1)  
  }
  predictor <- obs.at.node[,predictor.idx]
  split.set <- get_set_valid_split(X.pred = predictor,
                                   n.term.left = n.term.left,
                                   n.term.right = n.term.right,
                                   obs.per.term = obs.per.term)
  if(length(split.set) == 1){
    split.value <- split.set
  } else{
    if(is.numeric(split.set)){
      if(cont.unif){
        cat('extremes split set: ', min(split.set), ',', max(split.set), '\n')
        split.value <- runif(1, min = min(split.set), max = max(split.set))
      } else {
        split.value <- sample(split.set, 1)    
      }  
    } else {
      split.value <- sample(split.set, 1)
    }
    
  }
  return(list(cond = list(x.idx = predictor.idx, x.val = split.value),
              obs.at.node = obs.at.node,
              valid.pred = valid.pred,
              valid.split = split.set,
              n.left = n.term.left,
              n.right = n.term.right))
}

#' Set the splitting rule of an internal node of a binary tree
#' 
#' The function takes in input a node index, a tree topology and a splitting rule and returns a 
#' tree with the splitting rule assigned to the node corresponding to the index.
#'
#' @param node.idx a scalar representing the index of the internal node at which the condition is referred
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node, 
#' `cond` (only for internal nodes) splitting rules. This is needed only for the nodes before the one at which we want to generate the splitting rule. 
#' @param cond `list` representing a splitting rule in the form `list(x.idx = .., x.val = ..)` 
#' where `x.idx` represents the index of the value of `x` on which the condition is evaluated, 
#' and `x.val` is the splitting value of the condition. 
#' If `x.val` is numeric the condition is `x[x.idx] <= x.val`; if `x.val` is a character or a  
#' vectore of characters the condition is `x[x.idx] %in% x.val`.
#'
#' @return a `list` equal to `tree_top` but with the new splitting rule assigned to the node specified by `node.idx`
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(rnorm(100), sample(c('A', 'B'), 100, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot.idx(tree_top)
#' cond1 <- gen_node_condition(2, tree_top, XX)$cond
#' tree_top1 <- set_node_condition(2, tree_top, cond1)
#' cond2 <- gen_node_condition(3, tree_top, XX)$cond
#' tree_top2 <- set_node_condition(3, tree_top, cond2)
#' get_tree_plot(tree_top)
#' get_tree_plot(tree_top1)
#' get_tree_plot(tree_top2)
set_node_condition <- function(node.idx, tree_top, cond){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(tree_top)
  } else {
    if(tree_top$node.idx == node.idx){
      tree_top$cond = cond
      return(tree_top)
    } else {
      tree.left <- set_node_condition(node.idx, tree_top$left, cond)
      tree.right <- set_node_condition(node.idx, tree_top$right, cond)
      return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                  cond = tree_top$cond))
    }
  }
}


#' Assign valid splitting rules to binary tree
#' 
#' The function assigns valid splitting rules based on the observed predictor variables to all internal nodes of a binary tree provided in input.
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node.
#' @param X `data.frame` of observed values of the predictor variables
#' @param obs.per.term minimum number of observations assigned to the left and right branch of each internal node.
#' @param cont.unif `logical`, if `TRUE` the splitting values are sampled from a continuous uniform distribution on the
#'  on the range of the valid split values. If `FALSE` a uniform distribution on the obsersed valid split values is used. Default is `TRUE`.
#' @return a `list` equal to `tree_top` with the additional argument `cond` for each node representing the splitting rule at that node.
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(rnorm(100), sample(c('A', 'B'), 100, replace = TRUE))
#' tree_top1 <- assign_split_rules(tree_top, XX)
#' tree_top2 <- assign_split_rules(tree_top, XX)
#' tree_top1 <- assign_term_node_values(tree_top1, 0, 2)
#' tree_top2 <- assign_term_node_values(tree_top2, 0, 2)
#' get_tree_plot(tree_top1)
#' get_tree_plot(tree_top2)
assign_split_rules <- function(tree_top, X, obs.per.term = 1, cont.unif = TRUE){
  intern.node.idx <- get_internal_nodes_idx(tree_top)
  for(node.idx in intern.node.idx){
    new.cond <- gen_node_condition(node.idx, tree_top, X, obs.per.term, cont.unif = cont.unif)$cond
    tree_top <- set_node_condition(node.idx, tree_top, new.cond)
  }
  return(tree_top)
}

#' Set value at a terminal node of a binary tree
#' 
#' The function set a value at the terminal node of a binary tree. The binary tree and the index of the terminal node are provided as input.
#' The value at the terminal node is randomly exstracted from a Gaussian distribution with mean and standard deviation provided in input.
#'
#' @param node.idx a scalar representing the index of the terminal node at which we want to set the value.
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node.
#' @param mu a scalar representing the mean of the Gaussian distribution from which the value is sampled
#' @param sigma a scalar representing the standard deviation of the Gaussian distribution from which the value is sampled
#'
#' @return a `list` representing a tree equal to `tree_top` but with a new value at the provided terminal node 
#' @export
#'
#' @examples 
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(rnorm(100), sample(c('A', 'B'), 100, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot.idx(tree_top)
#' term.node.idx <- get_terminal_nodes_idx(tree_top)
#' tree_top1 <- set_term_node_value(term.node.idx[1], tree_top, 0, 2)
#' tree_top2 <- set_term_node_value(term.node.idx[1], tree_top, 0, 2)
#' get_tree_plot(tree_top1)
#' get_tree_plot(tree_top2)
set_term_node_value <- function(node.idx, tree_top, mu, sigma){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = rnorm(1, mean = mu, sd = sigma)
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value(node.idx, tree_top$left, mu, sigma)
    tree.right <- set_term_node_value(node.idx, tree_top$right, mu, sigma)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}


set_term_node_value_binary <- function(node.idx, tree_top, alpha.prior, beta.prior){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = rbeta(1, alpha.prior, beta.prior)
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_binary(node.idx, tree_top$left, alpha.prior, beta.prior)
    tree.right <- set_term_node_value_binary(node.idx, tree_top$right, alpha.prior, beta.prior)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

set_term_node_value_fixed <- function(node.idx, tree_top, node.value){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = node.value
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_fixed(node.idx, tree_top$left, node.value)
    tree.right <- set_term_node_value_fixed(node.idx, tree_top$right, node.value)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}


#' Assign value to all terminal nodes of a binary tree
#' 
#' The function assigns a value to all the terminal nodes of a binary tree provided in input.
#' The value at the terminal node is randomly extracted from a Gaussian distribution with mean and standard deviation provided in input.
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node,
#' `cond` (only for internal nodes) `list` representing a splitting rule
#' @param mu a scalar representing the mean of the Gaussian distribution from which the value is sampled
#' @param sigma a scalar representing the standard deviation of the Gaussian distribution from which the value is sampled
#'
#' @return a `list` representing a binary tree. The tree is the same as the one provided in input as `tree_top` but with 
#' new values at the terminal nodes
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(rnorm(100), sample(c('A', 'B'), 100, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top1 <- assign_term_node_values(tree_top, 0, 2)
#' tree_top2 <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot(tree_top1)
#' get_tree_plot(tree_top2)
assign_term_node_values <- function(tree_top, mu, sigma){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    tree_top <- set_term_node_value(node.idx, tree_top, mu, sigma)
  }
  return(tree_top)
}



assign_term_node_values_binary <- function(tree_top, alpha.prior, beta.prior){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    tree_top <- set_term_node_value_binary(node.idx, tree_top, alpha.prior, beta.prior)
  }
  return(tree_top)
}


assign_term_node_values_fixed <- function(tree_top, node.values){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(idx in 1:length(term.node.idx)){
    tree_top <- set_term_node_value_fixed(term.node.idx[idx], tree_top, node.values[idx])
  }
  return(tree_top)
}


#' Generate a random binary tree topology with fixed number of leaves
#'
#' The function returns a tree topology as nested `list` representing a random binary tree with fixed number of leaves but free depth.
#' 
#' @param num_terminals a scalar representing the number of leaves of the binary tree
#'
#' @return #' A `list` representing a binary tree. The tree is a nested `list` in which each element represents a node and has elements 
#' `left` a list representing the left branch
#' `right` a list representing the right branch.
#' if `left = NULL` and `right = NULL` the node is a terminal node.
#' @export
#'
#' @examples
#' tree_top1 <- generate_random_binary_tree_depth_free(6)
#' tree_top1 <- assign_node_idx(tree_top1)
#' tree_top2 <- generate_random_binary_tree_depth_free(6)
#' tree_top2 <- assign_node_idx(tree_top2)
#' get_tree_plot.idx(tree_top1)
#' get_tree_plot.idx(tree_top2)
generate_random_binary_tree_depth_free <- function(num_terminals){
  if (num_terminals == 1) {
    return(list(left=NULL, right=NULL))
  }
  if(num_terminals == 2) {
    return(list(left = list(left = NULL, right = NULL),
                right = list(left = NULL, right = NULL)))
  }
  poss_left_nodes <- 1:(num_terminals - 1)
  if(length(poss_left_nodes) == 1){
    left_nodes = poss_left_nodes
    right_nodes = num_terminals - left_nodes
  } else{
    left_nodes = sample(poss_left_nodes, 1)
    right_nodes = num_terminals - left_nodes
  }
  left_subtree <- generate_random_binary_tree_depth_free(left_nodes)
  right_subtree <- generate_random_binary_tree_depth_free(right_nodes) 
  return(list(left=left_subtree, right=right_subtree))
}



#' Generate a random binary tree topology with fixed depth and number of leaves
#'
#' The function returns a tree topology as nested `list` representing a random binary tree with fixed depth and number of leaves. 
#'
#' @param depth a scalar representing the depth of the binary tree
#' @param num_terminals a scalar representing the number of leaves of the binary tree
#' @param display.out `logical` if `TRUE` some information on the generating process are shown. Default is set to `FALSE`.
#'
#' @return #' A `list` representing a binary tree. The tree is a nested `list` in which each element represents a node and has elements 
#' `left` a list representing the left branch
#' `right` a list representing the right branch.
#' if `left = NULL` and `right = NULL` the node is a terminal node.
#' @export
#'
#' @examples 
#' tree_top1 <- generate_random_binary_tree(4, 8)
#' tree_top1 <- assign_node_idx(tree_top1)
#' tree_top2 <- generate_random_binary_tree(4, 8)
#' tree_top2 <- assign_node_idx(tree_top2)
#' get_tree_plot.idx(tree_top1)
#' get_tree_plot.idx(tree_top2)
generate_random_binary_tree_depth_nterm <- function(depth, num_terminals,
                                                    display.out = FALSE) {
  # if (depth == 0) {
  #   return(list(left=NULL, right=NULL))
  # }
  if (num_terminals == 1) {
    return(list(left=NULL, right=NULL))
  }
  if(num_terminals == 2) {
    return(list(left = list(left = NULL, right = NULL),
                right = list(left = NULL, right = NULL)))
  }
  min_on_side = depth
  poss_left_nodes <- (max(1, num_terminals - 2^(depth - 1))):(min(2^(depth-1), num_terminals - 1))
  if(depth > (num_terminals - depth)){
    idx_rm_left_nodes <- poss_left_nodes < min_on_side & poss_left_nodes > (num_terminals - min_on_side)
    poss_left_nodes <- poss_left_nodes[!idx_rm_left_nodes]
  }
  if(display.out){
    cat('depth: ', depth, '\n')
    cat('num terminal:', num_terminals, '\n')
    cat('poss: ', poss_left_nodes, '\n')  
  }
  if(length(poss_left_nodes) == 0){
    return(list(left = NULL, right = NULL))
  }
  if(length(poss_left_nodes) == 1){
    left_nodes = poss_left_nodes
    right_nodes = num_terminals - left_nodes
  } else{
    left_nodes = sample(poss_left_nodes, 1)
    right_nodes = num_terminals - left_nodes
  }
  if(display.out){
    cat('left:', left_nodes, '\n')
    cat('----- \n')
  }
  if(left_nodes >= depth & left_nodes <= (num_terminals - depth) ){
    left_subtree <- generate_random_binary_tree_depth_nterm(depth-1, left_nodes, 
                                                            display.out = display.out)
    right_subtree <- generate_random_binary_tree_depth_nterm(depth-1, right_nodes,
                                                             display.out = display.out)  
  } else if(left_nodes > max(depth - 1, num_terminals - depth)){
    left_subtree <- generate_random_binary_tree_depth_nterm(depth-1, left_nodes,
                                                            display.out = display.out)
    right_subtree <- generate_random_binary_tree_depth_free(right_nodes)  
  } else {
    left_subtree <- generate_random_binary_tree_depth_free(left_nodes)
    right_subtree <- generate_random_binary_tree_depth_nterm(depth-1, right_nodes,
                                                             display.out = display.out)  
  }
  
  return(list(left=left_subtree, right=right_subtree))
}


# 
generate_random_binary_tree_n_delta <- function(n, delt){
  if(n == 1){
    return(list(left = NULL, right = NULL))
  }
  if(f.odd(n) != f.odd(delt)){
    stop('Please provide a valid (n,delt) couple. If n is even delt should be even. If n is odd delt should be odd.')
  } else {
    if(f.even(n) != f.even(delt)){
      stop('Please provide a valid (n,delt) couple. If n is even delt should be even. If n is odd delt should be odd.')
    }
  }
  flag = sample(c(0,1), 1)
  if(flag){
    nl = ( n + delt )/ 2
    nr = (n - delt )/2
  } else {
    nl = ( n - delt )/ 2
    nr = (n + delt )/2
  }
  tree.left <- generate_random_binary_tree_depth_free(nl)
  tree.right <- generate_random_binary_tree_depth_free(nr)
  list(left = tree.left, right = tree.right)
}

# generate random binary tree with fixed nterm and delta
generate_random_binary_tree_n_delta <- function(n, delt){
  if(n == 1){
    return(list(left = NULL, right = NULL))
  }
  if(f.odd(n) != f.odd(delt)){
    stop('Please provide a valid (n,delt) couple. If n is even delt should be even. If n is odd delt should be odd.')
  } else {
    if(f.even(n) != f.even(delt)){
      stop('Please provide a valid (n,delt) couple. If n is even delt should be even. If n is odd delt should be odd.')
    }
  }
  flag = sample(c(0,1), 1)
  if(flag){
    nl = ( n + delt )/ 2
    nr = (n - delt )/2
  } else {
    nl = ( n - delt )/ 2
    nr = (n + delt )/2
  }
  tree.left <- generate_random_binary_tree_depth_free(nl)
  tree.right <- generate_random_binary_tree_depth_free(nr)
  list(left = tree.left, right = tree.right)
}

generate_random_binary_tree_from_prior <- function(omeg, gam){
  n_ <- rnterm(1, omeg)
  if(n_ > 1){
    if(f.odd(n_)){
      delta.values <- seq(1,n_- 2,by = 2)  
    } else{
      delta.values <- seq(0,n_- 2,by = 2)
    }
    delta_ <- sample(x = delta.values, 
                     size = 1,
                     prob = prior.delta(delta.values, gam, n_))
    rt <- generate_random_binary_tree_n_delta(n_, delta_)  
  } else {
    delta_ = 0
    rt <- generate_random_binary_tree_depth_free(n_)
  }
  return(list(tree = rt, n = n_, delta = delta_))
}




# sample value at the terminal node from conjugate gaussian posterior distribution
sample.cond.mu <- function(tree_top = NULL, 
                           node.idx = NULL, 
                           mu.prior.mean, 
                           mu.prior.var, 
                           X = NULL, 
                           Y = NULL, 
                           Y.var, 
                           Y.at.node = NULL){
  if(is.null(Y.at.node)){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    Y.at.node <- Y[as.numeric(rownames(obs.at.node))]
  } 
  nobs.at.node <- length(Y.at.node)
  Y.mean <- mean(Y.at.node)
  mu.cond.mean = (Y.var/(Y.var + nobs.at.node*mu.prior.var))*mu.prior.mean + (nobs.at.node*mu.prior.var/(Y.var + nobs.at.node*mu.prior.var))*Y.mean
  mu.cond.var = 1/(1/mu.prior.var + nobs.at.node/Y.var)
  #c(mu.cond.mean, mu.cond.var)
  return(rnorm(1, mean = mu.cond.mean, sd = sqrt(mu.cond.var)))
}

# sample conditional mu assuming binary observations and beta prior
sample.cond.mu.beta <- function(tree_top = NULL, 
                           node.idx = NULL, 
                           alpha.prior, 
                           beta.prior, 
                           X = NULL, 
                           Y = NULL, 
                           Y.at.node = NULL){
  if(is.null(Y.at.node)){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    Y.at.node <- Y[as.numeric(rownames(obs.at.node))]
  } 
  nobs.at.node <- length(Y.at.node)
  succ. <- sum(Y.at.node)
  insucc. <- nobs.at.node - succ.
  return(rbeta(1, shape1 = alpha.prior + succ., shape2 = beta.prior + insucc.))
}


set_term_node_value_cond <- function(node.idx, tree_top, mu.prior.mean, mu.prior.var, 
                                     X, Y, Y.var, 
                                     Y.at.node=NULL, binary = FALSE){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = sample.cond.mu(tree_top = tree_top, 
                                      node.idx = node.idx, 
                                      mu.prior.mean = mu.prior.mean, 
                                      mu.prior.var = mu.prior.var, 
                                      X = X, 
                                      Y = Y, Y.var = Y.var, Y.at.node = Y.at.node)
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_cond(node.idx, tree_top$left, mu.prior.mean, mu.prior.var, X, Y, Y.var, 
                                          Y.at.node)
    tree.right <- set_term_node_value_cond(node.idx, tree_top$right, mu.prior.mean, mu.prior.var, X, Y, Y.var, 
                                           Y.at.node)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}


set_term_node_value_cond_binary <- function(node.idx, tree_top, alpha.prior, beta.prior, 
                                     X, Y, 
                                     Y.at.node=NULL, binary = FALSE){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = sample.cond.mu.beta(tree_top = tree_top, 
                                      node.idx = node.idx, 
                                      alpha.prior = alpha.prior,
                                      beta.prior = beta.prior,
                                      X = X, 
                                      Y = Y, Y.at.node = Y.at.node)
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_cond_binary(node.idx, tree_top$left, alpha.prior, beta.prior, X, Y, 
                                          Y.at.node)
    tree.right <- set_term_node_value_cond_binary(node.idx, tree_top$right, alpha.prior, beta.prior, X, Y, 
                                           Y.at.node)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}



assign_term_node_values_cond <- function(tree_top, mu.prior.mean, mu.prior.var, X, Y, Y.var, Y.at.node = NULL){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    Y.at.node <- Y[as.numeric(rownames(obs.at.node))]
    tree_top <- set_term_node_value_cond(node.idx = node.idx, tree_top = tree_top, 
                                         mu.prior.mean = mu.prior.mean, 
                                         mu.prior.var = mu.prior.var, 
                                         X = X, Y =Y , Y.var = Y.var, Y.at.node = Y.at.node)
  }
  return(tree_top)
}

assign_term_node_values_cond_binary <- function(tree_top, alpha.prior, beta.prior, 
                                                X, Y, Y.at.node = NULL){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    Y.at.node <- Y[as.numeric(rownames(obs.at.node))]
    tree_top <- set_term_node_value_cond_binary(node.idx = node.idx, tree_top = tree_top, 
                                               alpha.prior = alpha.prior,
                                               beta.prior, beta.prior,
                                               X = X, Y =Y, Y.at.node = Y.at.node)
  }
  return(tree_top)
}


