#' Evaluate a splitting rule
#' 
#' @param cond 
#' `list` representing a splitting rule in the form `list(x.idx = .., x.val = ..)` 
#' where `x.idx` represents the index of the value of `x` on which the condition is evaluated, 
#' and `x.val` is the splitting value of the condition. 
#' If `x.val` is numeric the condition is `x[x.idx] <= x.val`; if `x.val` is a character or a  
#' vectore of characters the condition is `x[x.idx] %in% x.val`.
#' @param x 
#' A row of a `data.frame` representing an observation on which the condition is evaluated
#'
#' @return `logical` TRUE if the condition is verified and FALSE otherwise
#'
#' @examples 
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 1, x.val = 0)
#' eval_cond(cond, x)
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 1, x.val = 2)
#' eval_cond(cond, x)
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 2, x.val = c('A', 'B')
#' eval_cond(cond, x)
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 2, x.val = 'C')
#' eval_cond(cond, x)

eval_cond <- function(cond, x){
  if(is.character(cond$x.val)){
    return(x[cond$x.idx] %in% cond$x.val)
  } else if(is.numeric(cond$x.val)){
    return(x[cond$x.idx] <= cond$x.val)
  } else{
    stop('Unknown condition value')
  }
}

#' Number of leaves for a given binary tree
#'
#' The function returns the number of leaves of the binary tree given as input 
#' @param tree a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return a scalar representing the number of leaves of `tree`
#' @export
#'
#' @examples 
#' tree_top <- generate_random_binary_tree(4,5)
#' get_num_terminal_nodes(tree_top)
get_num_terminal_nodes <- function(tree) {
  if (is.null(tree)) {
    return(0)
  }
  if (is.null(tree$left) && is.null(tree$right)) {
    return(1)
  }
  return(get_num_terminal_nodes(tree$left) + get_num_terminal_nodes(tree$right))
}

#' Depth plus 1 of a binary tree
#'
#' The function returns the depth + 1 of the binary tree given as input
#' @param tree a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return the value of the depth of the tree plus 1
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' get_depth_plus1(tree_top)
get_depth_plus1 <- function(tree) {
  if (is.null(tree)) {
    return(0)
  } else {
    left_depth <- get_depth_plus1(tree$left) 
    right_depth <- get_depth_plus1(tree$right) 
    return(max(left_depth, right_depth) + 1)
  }
}

#' Depth of a binary tree
#'
#' The function returns the depth of the binary tree given as input
#' @param tree a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return the value of the depth of the tree
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' get_depth(tree_top)
get_depth <- function(tree) {
  get_depth_plus1(tree) - 1
}

#' Observations on a node of a binary tree
#'
#' The function retrieve the observations at a given node of a binary tree provided as input
#'
#' @param node.idx a scalar representing the index of the node. The nodes are indexed going from left to right 
#' @param X a `data.frame` representing the full set of observations 
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param node.count a scalar used for internal calculations, the default value is 1. It does not need to be changed.
#'
#' @return a `data.frame` representing the observations corresponding to a node
#' @export
#'
#' @examples
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot(tree_top)
#' XX.at.node <- get_obs_at_node(node.idx = 2, X = XX, tree_top = tree_top)
#' nrow(XX.at.node)
get_obs_at_node <- function(node.idx = 2, X, tree_top, node.count = 1, X.orig){
  X.sel = X
  # this is needed because if X is matrix (which is way faster) then X[1,] is a numeric vector and nrow(X.sel) is NULL
  if(class(X.sel)[1] == 'numeric'){
    X.sel2 <- matrix(X.sel, byrow = TRUE, ncol = length(X.sel))
    colnames(X.sel2) <- names(X.sel)
    idx.obs <- which(apply(apply(X.orig, 1, \(x) x == X.sel2),2,all))
    rownames(X.sel2) <- idx.obs
    X.sel <- X.sel2
  }
  if(nrow(X.sel) == 0){
    return(X.sel)
  }
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(node.idx == node.count){
      return(X.sel)
    } else {
      return(NULL)  
    }
  } else{
    if(node.idx == node.count){
      return(X.sel)
    } else {
      # evaluate condition
      cond.sel.idx <- vapply(1:nrow(X.sel), \(x) eval_cond(tree_top$cond, 
                                                           X.sel[x,]), TRUE)
      X.sel.left <- X.sel[cond.sel.idx,]
      X.sel.right <- X.sel[!cond.sel.idx,]
      
      node.count.left = node.count + 1
      n.term.left <- get_num_terminal_nodes(tree_top$left)
      n.node.left <- 2*n.term.left 
      node.count.right = node.count + n.node.left
      
      if(node.idx < node.count.right){
        X.sel.left.new <- get_obs_at_node(node.idx = node.idx, 
                                          node.count = node.count.left,
                                          X = X.sel.left,
                                          tree_top = tree_top$left,
                                          X.orig = X.orig)
        return(X.sel.left.new)
      } else {
        X.sel.right.new <- get_obs_at_node(node.idx = node.idx,
                                           node.count = node.count.right,
                                           X = X.sel.right,
                                           tree_top = tree_top$right,
                                           X.orig = X.orig)
        return(X.sel.right.new)
      }
      
    }  
  }
}

#' Terminal node index values
#'
#' The function returns the index corresponding to the terminal nodes of 
#' the binary tree provided as input. The indexes are placed going from left
#' to right
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param counter a scalar used for internal calculations, the default value is 1. It does not need to be changed.
#'
#' @return the node indexes corresponding to the leaves of the binary tree
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
#' get_terminal_nodes_idx(tree_top)
get_terminal_nodes_idx <- function(tree_top, counter = 1){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(counter)
  } else {
    count.l <- counter + 1
    n.term.left <- get_num_terminal_nodes(tree_top$left)
    n.node.left <- 2*n.term.left 
    count.r = counter + n.node.left
    idx.left <- get_terminal_nodes_idx(tree_top$left, counter = count.l)
    idx.right <- get_terminal_nodes_idx(tree_top$right, counter = count.r)
    return(unique(c(idx.left, idx.right)))
  }
}

#' Internal node index values
#'
#' The function returns the index corresponding to the internal nodes of 
#' the binary tree provided as input. The indexes are placed going from left
#' to right
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param counter a scalar used for internal calculations, the default value is 1. It does not need to be changed.
#'
#' @return the node indexes corresponding to the internal nodes of the binary tree
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
#' get_internal_nodes_idx(tree_top)
get_internal_nodes_idx <- function(tree_top){
  n.term <- get_num_terminal_nodes(tree_top)
  n.nodes <- 2*n.term - 1
  int.idx.raw <- 1:n.nodes
  term.idx <- get_terminal_nodes_idx(tree_top)
  return(int.idx.raw[-term.idx])
}

#' Create name from splitting rule
#'
#' The function returns a string representing the condition provided as input. 
#' It is used by the `get_tree_plot()` function to fill the internal nodes.
#' @param cond `list` representing a splitting rule in the form `list(x.idx = .., x.val = ..)` 
#' where `x.idx` represents the index of the value of `x` on which the condition is evaluated, 
#' and `x.val` is the splitting value of the condition. 
#' If `x.val` is numeric the condition is `x[x.idx] <= x.val`; if `x.val` is a character or a  
#' vectore of characters the condition is `x[x.idx] %in% x.val`.
#'
#' @return a string representing the condition
#' @export
#'
#' @examples
#' cond <- list(x.idx = 1, x.val = 2)
#' node.name.from.cond(cond)
#' cond <- list(x.idx = 1, x.val = c('A', 'B'))
#' node.name.from.cond(cond)
node.name.from.cond <- function(cond){
  if(is.character(cond$x.val)){
    return(paste0('X', cond$x.idx, ' in ', paste(cond$x.val, collapse = ',')))
  } else{
    return(paste0('X', cond$x.idx, '<=', round(cond$x.val, 3)))
  }
}

#' Reshape tree for plotting
#' 
#' The function reshapes the tree topology in a format suitable to be used by 
#' the `data.tree` R-pacakge. It is used by the `get_tree_plot` function.
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return a nested `list` representing a binary tree 
#' @export
#'
#' @examples 
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' tree.to.plot(tree_top)
tree.to.plot <- function(tree_top){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(list(left = NULL, 
                right = NULL,
                node.name = round(tree_top$value,3)))
  } else {
    return(list(children = list(left = tree.to.plot(tree_top$left),
                                right = tree.to.plot(tree_top$right)
    ),
    node.name = node.name.from.cond(tree_top$cond)
    )
    )
  }
}

#' Reshape tree for plotting
#' 
#' The function reshapes the tree topology in a format suitable to be used by 
#' the `data.tree` R-pacakge. It is used by the `get_tree_plot.idx` function.
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node.
#'
#' @return a nested `list` representing a binary tree 
#' @export
#'
#' @examples 
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' tree.to.plot.idx(tree_top)
tree.to.plot.idx <- function(tree_top){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(list(node.name = tree_top$node.idx))
  } else {
    return(list(children = list(left = tree.to.plot.idx(tree_top$left),
                                right = tree.to.plot.idx(tree_top$right)),
                node.name = tree_top$node.idx)
    )
  }
}

#' Plot a binary tree with splitting rules and values at terminal nodes
#' 
#' The function plots a binary tree given as input. For the internal nodes, 
#' the splitting rules are shown, assuming that if they are true the left branch is activated.
#' For the terminal nodes, the value of the leaves is shown.
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return a plot of the tree showing at each internal node the splitting rule and at each terminal node a value.
#' @export
#'
#' @examples
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot(tree_top)
get_tree_plot <- function(tree_top){
  mod.tree.top <- tree.to.plot(tree_top)
  tree.d <- as.Node(mod.tree.top, 
                    mode = 'explicit', 
                    childrenName = 'children', 
                    nameName = 'node.name')
  plot(tree.d)
  
}

#' Plot a binary tree with nodes index value
#' 
#' The function plots a binary tree given as input. The index of each node is shown.
#' 
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node.
#' 
#' @return a plot of the binary tree showing the index value of each node
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
get_tree_plot.idx <- function(tree_top){
  mod.tree.top <- tree.to.plot.idx(tree_top)
  tree.d <- as.Node(mod.tree.top, 
                    mode = 'explicit', 
                    childrenName = 'children', 
                    nameName = 'node.name')
  plot(tree.d)
  
}

# extract offsprings of a given node
get_offsprings <- function(node.idx, tree_top){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(NULL)
  } else {
    if(tree_top$node.idx == node.idx){
      return(tree_top)
    } else {
      tree.left <- get_offsprings(node.idx, tree_top$left)
      tree.right <- get_offsprings(node.idx, tree_top$right)
      if(is.null(tree.right)){
        return(tree.left)
      } else {
        return(tree.right)
      }
    }
  }
}

# extract condition at node for switch
get_node_condition <- function(node.idx, tree_top){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(NULL)
  } else {
    if(tree_top$node.idx == node.idx){
      return(list(x.idx = tree_top$cond$x.idx, x.val = tree_top$cond$x.val))
    } else {
      cond.left <- get_node_condition(node.idx, tree_top$left)
      cond.right <- get_node_condition(node.idx, tree_top$right)
      if(is.null(cond.right)){
        return(cond.left)
      } else {
        return(cond.right)
      }
    }
  }
}


# apply a fun to a list of models
apply_fun_models <- function(fun_, mcmc.list, born.out.pc = 100, n.chain = 10, sample.pc = 1000){
  if(born.out.pc == 0){
    fun.list <- lapply(mcmc.list, \(model.res) vapply(model.res$trees, fun_, 1))
    names(fun.list) <- names(mcmc.list)
    fun.df.list <- lapply(1:length(fun.list), \(x) data.frame(x = 1:length(mcmc.list[[x]]$trees),
                                                              y = fun.list[[x]],
                                                              panel.name = names(fun.list)[x]))
    fun.df <- Reduce(rbind, fun.df.list)
  }
  else{
    idx.bornout <- sort(unlist(outer(1:born.out.pc, sample.pc*(0:(n.chain - 1)), '+')))
    ntrees.per.model <- vapply(mcmc.list, \(model.res) length(model.res$trees), 0)
    if(sum(born.out.pc > ntrees.per.model) > 0){
      stop('Born out greater than number of samples')
    }
    fun.list <- lapply(mcmc.list, \(model.res) vapply(model.res$trees[-idx.bornout], fun_, 1))
    names(fun.list) <- names(mcmc.list)
    fun.df.list <- lapply(1:length(fun.list), \(x) data.frame(x = 1:length(mcmc.list[[x]]$trees[-idx.bornout]),
                                                              y = fun.list[[x]],
                                                              panel.name = names(fun.list)[x]))
    fun.df <- Reduce(rbind, fun.df.list)
  }
  
  fun.df
}


# fun to predict
binary.pred.tree <- function(tree_top, X){
  get_value_tree(tree_top, X) >= 0.5
}

binary.pred.mcmc <- function(tree_list, X, born.out.pc = 100, sample.pc = 1000, n.chain = 10){
  if(born.out.pc == 0){
    pred.per.tree <- lapply(tree_list, \(x) binary.pred.tree(x, X))
    pred.df <- Reduce(rbind, pred.per.tree)  
  } else {
    idx.bornout <- sort(unlist(outer(1:born.out.pc, sample.pc*(0:(n.chain - 1)), '+')))
    pred.per.tree <- lapply(tree_list[-idx.bornout], \(x) binary.pred.tree(x, X))
    pred.df <- Reduce(rbind, pred.per.tree)  
  }
  return(pred.df)
}


miss.class.rate.from.matrix <- function(pred.matrix, Y){
  as.numeric(apply(pred.matrix, 1, \(x) sum(x != Y)))
}

miss.class.rate.mcmc <- function(tree_list, X, Y, born.out.pc = 100, sample.pc = 1000, n.chain = 10){
  pred.matrix <- binary.pred.mcmc(tree_list = tree_list, 
                                  X = X, 
                                  born.out.pc = born.out.pc, 
                                  sample.pc = sample.pc, 
                                  n.chain = n.chain)
  miss.class.rate.from.matrix(pred.matrix, Y)
}


apply_fun_models_all.tog <- function(mcmc.list, fun_){
  fun.list <- lapply(mcmc.list, \(x) fun_(x))
  names(fun.list) <- names(mcmc.list)
  fun.df.list <- lapply(1:length(fun.list), \(x) data.frame(x = 1:length(mcmc.list[[x]]$trees),
                                                            y = fun.list[[x]],
                                                            panel.name = names(fun.list)[x]))
  fun.df <- Reduce(rbind, fun.df.list)
  fun.df
}


