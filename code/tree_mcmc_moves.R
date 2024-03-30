# check if an idx is good for swap
check_swap <- function(node.idx, term.node.idx){
  ! (( (node.idx + 1) %in% term.node.idx ) & ((node.idx + 2) %in% term.node.idx) )
}

# function to get valid set of swap nodes
get_valid_swap_idx <- function(tree_top){
  # take terminal node idx
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  all.node.idx <- 1:max(term.node.idx)
  # take internal node idx
  intern.node.idx <- all.node.idx[-term.node.idx]
  intern.node.idx.check <- vapply(intern.node.idx, \(x) check_swap(x, term.node.idx), TRUE)
  intern.node.idx.swap <- intern.node.idx[intern.node.idx.check]
  node.idx.df <- data.frame(idx.1 = NA, idx.2 = NA)[-1,]
  for(node.idx in intern.node.idx.swap){
    subtree <- get_offsprings(node.idx, tree_top)
    cond.idx <- subtree$cond$x.idx
    if(subtree$left$node.idx %in% intern.node.idx){
      if(subtree$left$cond$x.idx != cond.idx){
        node.idx.df <- rbind(node.idx.df, 
                             data.frame(idx.1 = node.idx, 
                                        idx.2 = subtree$left$node.idx) )  
      }
    } 
    if(subtree$right$node.idx %in% intern.node.idx){
      if(subtree$right$cond$x.idx != cond.idx){
        node.idx.df <- rbind(node.idx.df, 
                             data.frame(idx.1 = node.idx, 
                                        idx.2 = subtree$right$node.idx) )  
      }
    }
  }
  return(node.idx.df)
}

# swap two nodes conditions
swap_node_condition <- function(first.node.idx, second.node.idx, tree_top){
  cond1 <- get_node_condition(first.node.idx, tree_top)
  cond2 <- get_node_condition(second.node.idx, tree_top)
  tree_top1 <- set_node_condition(first.node.idx, tree_top, cond2)
  tree_top2 <- set_node_condition(second.node.idx, tree_top1, cond1)
  return(tree_top2)
}

# perform swap move
swap_move <- function(tree_top){
  swap.node.idx <- get_valid_swap_idx(tree_top)
  if(nrow(swap.node.idx) == 0){
    return(list(tree = tree_top, move = 'swap', valid.pred = NULL, valid.split = NULL))
  } else {
    swap.idx.row <- sample(1:nrow(swap.node.idx), 1)
    swap.idx <- swap.node.idx[swap.idx.row, ]
    swap.tree_top <- swap_node_condition(swap.idx[1,1], swap.idx[1,2], tree_top)
    return(list(tree = swap.tree_top, move = 'swap', node.idx = swap.idx, valid.pred = NULL, valid.split = NULL))
  }
}

# perform change move
change_move <- function(tree_top, X, obs.per.term = 2, cont.unif = TRUE){
  internal.node.idx <- get_internal_nodes_idx(tree_top)
  change.idx <- sample(internal.node.idx, 1)
  new.cond <- gen_node_condition(change.idx, tree_top, X, obs.per.term, cont.unif)
  change.tree_top <- set_node_condition(change.idx, tree_top, new.cond$cond)
  return(list(tree = change.tree_top, move = 'change', 
              node.idx = change.idx, 
              valid.pred = new.cond$valid.pred, 
              valid.split = new.cond$valid.split,
              n.left = new.cond$n.left,
              n.right = new.cond$n.right))
}

# given an idx of terminal, create a split
grow_terminal <- function(node.idx, tree_top, new.cond, value.l = 1, value.r = 1){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      return(list(left = list(left = NULL, right = NULL, value = value.l),
                  right = list(left = NULL, right = NULL, value = value.r),
                  cond = new.cond))
    } else {
      return(tree_top)
    }
  } else {
    tree.left <- grow_terminal(node.idx, tree_top$left, new.cond, value.l, value.r)
    tree.right <- grow_terminal(node.idx, tree_top$right, new.cond, value.l, value.r)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

# this is assuming obs.per.term = 1
check_grow_valid <- function(node.idx, tree_top, X){
  obs.at.node <- get_obs_at_node(node.idx = node.idx,
                                 X = X,
                                 tree_top = tree_top, X.orig = X)
  valid.pred <- get_set_valid_predictor(X.sel = obs.at.node,
                                        n.term.left = 1,
                                        n.term.right = 1,
                                        obs.per.term = 1)
  return(length(valid.pred) > 0)
}



# perform a grow move
grow_move <- function(tree_top, X, Y, Y.var, mu.prior.mean, mu.prior.var, obs.per.term = 1,
                      cont.unif = TRUE){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  # this is because if a terminal node has only 1 obs associated cannot be grow
  nobs.at.nodes <- vapply(term.node.idx, \(x) nrow(get_obs_at_node(node.idx = x, X = X, tree_top = tree_top, X.orig = X)),0)
  term.node.idx <- term.node.idx[nobs.at.nodes > 1]
  check.valid.idx <- vapply(term.node.idx, \(x) check_grow_valid(x, tree_top, X), TRUE)
  term.node.idx <- term.node.idx[check.valid.idx]
  if(length(term.node.idx) == 0){
    return(list(tree = tree_top, node.idx = NULL, move = 'grow', 
                valid.pred = NULL, 
                valid.split = NULL))
  }  else if(length(term.node.idx) > 1){ 
    grow.idx <- sample(term.node.idx, 1, replace = FALSE)
  } else if(length(term.node.idx) == 1){
    grow.idx <- term.node.idx
  }
  new.cond.list <- gen_node_condition(grow.idx, tree_top, X, obs.per.term, cont.unif = cont.unif, for.grow = TRUE)
  obs.at.node <- new.cond.list$obs.at.node
  Y.at.parent <- Y[as.numeric(rownames(obs.at.node))]
  if(is.numeric(new.cond.list$cond$x.val)){
    idx.left <- obs.at.node[,new.cond.list$cond$x.idx] <= new.cond.list$cond$x.val
  } else {
    idx.left <- obs.at.node[, new.cond.list$cond$x.idx] %in% new.cond.list$cond$x.val
  }
  idx.right <- !idx.left
  Y.at.left <- Y.at.parent[idx.left]
  Y.at.right <- Y.at.parent[idx.right]
  term.node.value.left <- sample.cond.mu(Y.at.node = Y.at.left, Y.var = Y.var,
                                         mu.prior.mean = mu.prior.mean, 
                                         mu.prior.var = mu.prior.var)
  term.node.value.right <- sample.cond.mu(Y.at.node = Y.at.right, Y.var = Y.var,
                                          mu.prior.mean = mu.prior.mean, 
                                          mu.prior.var = mu.prior.var)
  
  tree_top_grow <- grow_terminal(node.idx = grow.idx, 
                                 tree_top = tree_top, 
                                 new.cond = new.cond.list$cond, 
                                 value.l = term.node.value.left,
                                 value.r = term.node.value.right)
  tree_top_grow <- assign_node_idx(tree_top_grow)
  return(list(tree = tree_top_grow, node.idx = grow.idx, move = 'grow', 
              valid.pred = new.cond.list$valid.pred, 
              valid.split = new.cond.list$valid.split))
}

# grow move for binary observations
grow_move_binary <- function(tree_top, X, Y, alpha.prior, beta.prior, cont.unif = TRUE, obs.per.term = 1){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  # this is because if a terminal node has only 1 obs associated cannot be grow
  nobs.at.nodes <- vapply(term.node.idx, \(x) nrow(get_obs_at_node(node.idx = x, X = X, tree_top = tree_top, X.orig = X)),0)
  term.node.idx <- term.node.idx[nobs.at.nodes > 1]
  check.valid.idx <- vapply(term.node.idx, \(x) check_grow_valid(x, tree_top, X), TRUE)
  term.node.idx <- term.node.idx[check.valid.idx]
  if(length(term.node.idx) == 0){
    return(list(tree = tree_top, node.idx = NULL, move = 'grow', 
                valid.pred = NULL, 
                valid.split = NULL))
  }  else if(length(term.node.idx) > 1){ 
    grow.idx <- sample(term.node.idx, 1, replace = FALSE)
  } else if(length(term.node.idx) == 1){
    grow.idx <- term.node.idx
  }
  new.cond.list <- gen_node_condition(grow.idx, tree_top, X, obs.per.term, cont.unif = cont.unif, for.grow = TRUE)
  obs.at.node <- new.cond.list$obs.at.node
  Y.at.parent <- Y[as.numeric(rownames(obs.at.node))]
  if(is.numeric(new.cond.list$cond$x.val)){
    idx.left <- obs.at.node[,new.cond.list$cond$x.idx] <= new.cond.list$cond$x.val
  } else {
    idx.left <- obs.at.node[, new.cond.list$cond$x.idx] %in% new.cond.list$cond$x.val
  }
  idx.right <- !idx.left
  Y.at.left <- Y.at.parent[idx.left]
  Y.at.right <- Y.at.parent[idx.right]
  term.node.value.left <- sample.cond.mu.beta(Y.at.node = Y.at.left,
                                         alpha.prior = alpha.prior,
                                         beta.prior = beta.prior)
  term.node.value.right <- sample.cond.mu.beta(Y.at.node = Y.at.right,
                                               alpha.prior = alpha.prior,
                                               beta.prior = beta.prior)
  
  tree_top_grow <- grow_terminal(node.idx = grow.idx, 
                                 tree_top = tree_top, 
                                 new.cond = new.cond.list$cond, 
                                 value.l = term.node.value.left,
                                 value.r = term.node.value.right)
  tree_top_grow <- assign_node_idx(tree_top_grow)
  return(list(tree = tree_top_grow, node.idx = grow.idx, move = 'grow', 
              valid.pred = new.cond.list$valid.pred, 
              valid.split = new.cond.list$valid.split))
}

# check if a node can be pruned 
check_prune_valid <- function(node.idx, tree_top){
  node.offspring <- get_offsprings(node.idx, tree_top)
  return(is.null(node.offspring$left$left) & 
           is.null(node.offspring$left$right) & 
           is.null(node.offspring$right$left) & 
           is.null(node.offspring$right$right))
}

# get idx of nodes that can be pruned
get_prune_idx <- function(tree_top){
  intern.node.idx <- get_internal_nodes_idx(tree_top)
  check.prune <- vapply(intern.node.idx, \(x) check_prune_valid(x, tree_top), TRUE)
  return(intern.node.idx[check.prune])
}


# prune a node given the idx
prune_terminal <- function(node.idx, tree_top, value.node){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(tree_top) 
  }
  if(tree_top$node.idx == node.idx){
    return(list(left = NULL, right = NULL, value = value.node))
  }  else {
    tree.left <- prune_terminal(node.idx, tree_top$left, value.node)
    tree.right <- prune_terminal(node.idx, tree_top$right, value.node)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

# perform prune move
prune_move <- function(tree_top, mu.prior.mean, mu.prior.var, X, Y, Y.var){
  prune.node.idx <- get_prune_idx(tree_top)
  if(length(prune.node.idx) == 1){
    prune.idx <- prune.node.idx
  } else{
    prune.idx <- sample(prune.node.idx, 1, replace = FALSE)
  }
  prune.value <- sample.cond.mu(tree_top = tree_top,
                                node.idx = prune.idx,
                                mu.prior.mean = mu.prior.mean,
                                mu.prior.var = mu.prior.var,
                                X = X,
                                Y = Y, Y.var = Y.var)
  tree_top_prune <- prune_terminal(node.idx = prune.idx, 
                                   tree_top = tree_top, 
                                   value.node = prune.value)
  tree_top_prune <- assign_node_idx(tree_top_prune)
  return(list(tree = tree_top_prune, move = 'prune', node.idx = prune.idx,
              prune.node.idx = prune.node.idx))
}

# prune move for binary observations
prune_move_binary <- function(tree_top, alpha.prior, beta.prior, X, Y){
  prune.node.idx <- get_prune_idx(tree_top)
  if(length(prune.node.idx) == 1){
    prune.idx <- prune.node.idx
  } else{
    prune.idx <- sample(prune.node.idx, 1, replace = FALSE)
  }
  prune.value <- sample.cond.mu.beta(tree_top = tree_top,
                                node.idx = prune.idx,
                                alpha.prior = alpha.prior,
                                beta.prior = beta.prior,
                                X = X,
                                Y = Y)
  tree_top_prune <- prune_terminal(node.idx = prune.idx, 
                                   tree_top = tree_top, 
                                   value.node = prune.value)
  tree_top_prune <- assign_node_idx(tree_top_prune)
  return(list(tree = tree_top_prune, move = 'prune', node.idx = prune.idx,
              prune.node.idx = prune.node.idx))
}
