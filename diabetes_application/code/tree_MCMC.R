# the transition probabilities are independent from the nodes that are actually swapped
transition.p.swap <- function(tree_top){
  # valid.swap.idx <- get_valid_swap_idx(tree_top)
  # 1/nrow(valid.swap.idx)
  1
}

# here we assume it is always continuous uniform
transition.p.change <- function(change_move_list){
  # n.int <- length(get_internal_nodes_idx(change_move_list$tree))
  # if(is.character(change_move_list$valid.split)){
  #   unique.valid.split <- unique(change_move_list$valid.split)
  #   prob.split <- 1/length(unique.valid.split)
  # } else {
  #   if(length(change_move_list$valid.split) == 1){
  #     prob.split <- 1
  #   } else {
  #     valid.set.extremes <- range(change_move_list$valid.split)
  #     prob.split <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))  
  #   }
  # }
  # (1/n.int)*(1/length(change_move_list$valid.pred))*prob.split
  1
}

transition.p.change.reverse <- function(change_move_list, X, old.tree){
  1
  # n.int <- length(get_internal_nodes_idx(change_move_list$tree))
  # old.pred <- get_node_condition(change_move_list$node.idx, old.tree)$x.idx
  # x.sel <- get_obs_at_node(change_move_list$node.idx, X, change_move_list$tree)
  # n.left <- change_move_list$n.left
  # n.right <- change_move_list$n.right
  # valid.pred <- get_set_valid_predictor(X.sel = x.sel, 
  #                                       n.term.left = n.left,
  #                                       n.term.right = n.right,
  #                                       obs.per.term = 1)
  # valid.split <- get_set_valid_split(X.pred = x.sel[,old.pred], 
  #                                    n.term.left = n.left, 
  #                                    n.term.right = n.right,
  #                                    obs.per.term = 1)
  # if(is.character(valid.split)){
  #   unique.valid.split <- unique(valid.split)
  #   prob.split <- 1/length(unique.valid.split)
  # } else {
  #   if(length(valid.split) == 1){
  #     prob.split <- 1
  #   } else {
  #     valid.set.extremes <- range(valid.split)
  #     prob.split <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))  
  #   }
  # }
  # 
  # (1/n.int)*(1/length(valid.pred))*prob.split
}

transition.p.grow <- function(grow_move_list, include.split, cont.unif = TRUE){
  n.term <- get_num_terminal_nodes(grow_move_list$tree) - 1
  if(include.split){
    n.valid.pred <- length(grow_move_list$valid.pred)
    if(is.character(grow_move_list$valid.split)){
      unique.valid.split <- unique(grow_move_list$valid.split)
      prob.split_ <- 1/length(unique.valid.split)
    } else {
      if(length(grow_move_list$valid.split) == 1){
        prob.split_ <- 1
      } else {
        if(cont.unif){
          valid.set.extremes <- range(grow_move_list$valid.split)
          prob.split_ <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))  
        } else {
          prob.split_ <- 1/length(grow_move_list$valid.split)
        }
        
      }
    }
    return((1/n.term)*(1/n.valid.pred)*prob.split_)
  } else{
    return(1/n.term)
  }
  
}

transition.p.grow.reverse <- function(grow_move_list){
  valid.prune.idx <- get_prune_idx(grow_move_list$tree)
  1/length(valid.prune.idx)
}

transition.p.prune <- function(prune_move_list){
  1/length(prune_move_list$prune.node.idx)
}

transition.p.prune.reverse <- function(prune_move_list, old.tree, X, obs.per.term,
                                       cont.unif = TRUE, include.split){
  n.term <- get_num_terminal_nodes(prune_move_list$tree)
  if(include.split){
    node.idx <- prune_move_list$node.idx
    obs.at.node <- get_obs_at_node(node.idx = node.idx, 
                                   X = X, 
                                   tree_top = prune_move_list$tree, X.orig = X)
    valid.pred <- get_set_valid_predictor(X.sel = obs.at.node, 
                                          n.term.left = 1, n.term.right = 1, 
                                          obs.per.term = 1)
    old.cond <- get_node_condition(node.idx, old.tree)
    valid.split <- get_set_valid_split(X.pred = obs.at.node[,old.cond$x.idx], 
                                       n.term.left = 1, n.term.right = 1, 
                                       obs.per.term = 1)
    if(is.character(valid.split)){
      unique.valid.split <- unique(valid.split)
      prob.split_ <- 1/length(unique.valid.split)
    } else {
      if(length(valid.split) == 1){
        prob.split_ <- 1
      } else {
        if(cont.unif){
          valid.set.extremes <- range(valid.split)
          prob.split_ <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))  
        } else{
          prob.split_ <- 1/length(valid.split)
        }
      }
    }
    return((1/n.term)*(1/length(valid.pred))*prob.split_) 
    } else{
      return(1/n.term)
  }
 
}

acceptance.prob <- function(move_list, old.tree, X, Y, omeg, eps, sigma, dm = NULL, max.depth = FALSE,
                            include.split, cont.unif = TRUE){
  if(identical(move_list$tree, old.tree)){
    return(list(prior.ratio = 1, lik.ratio = 1, 
                trans.ratio = 1, alpha = 1))
  }
  # priors
  if(max.depth){
    prior.old <- loss_based_prior_tree_dm(tree_top = old.tree, omeg = omeg, eps = eps, dm = dm)*
      prior.split.rule(old.tree, X, cont.unif = cont.unif)
    prior.new <- loss_based_prior_tree_dm(tree_top = move_list$tree, omeg = omeg, eps = eps, dm = dm)*
      prior.split.rule(move_list$tree, X, cont.unif = cont.unif)
  } else {
    prior.old <- loss_based_prior_tree(tree_top = old.tree, omeg = omeg, eps = eps)*
      prior.split.rule(old.tree, X, cont.unif = cont.unif)
    prior.new <- loss_based_prior_tree(tree_top = move_list$tree, omeg = omeg, eps = eps)*
      prior.split.rule(move_list$tree, X, cont.unif = cont.unif)
  }
  
  # likelihood 
  log.lik.old <- cart_log_lik(tree_top = old.tree, Y = Y, X = X, sigma = sigma)
  log.lik.new <- cart_log_lik(tree_top = move_list$tree, Y = Y, X = X, sigma = sigma)
  
  # transition probabilities
  if(move_list$move == 'swap'){
    prob.old.to.new <- transition.p.swap(old.tree)
    prob.new.to.old <- transition.p.swap(move_list$tree)
  } else if(move_list$move == 'change'){
    prob.old.to.new <- transition.p.change(move_list)
    prob.new.to.old <- transition.p.change.reverse(move_list, X, old.tree)
  } else if(move_list$move == 'grow'){
    prob.old.to.new <- transition.p.grow(move_list, cont.unif = cont.unif, 
                                         include.split = include.split)
    prob.new.to.old <- transition.p.grow.reverse(move_list)
  } else if(move_list$move == 'prune'){
    prob.old.to.new <- transition.p.prune(move_list)
    prob.new.to.old <- transition.p.prune.reverse(move_list, old.tree, X, cont.unif = cont.unif,
                                                  obs.per.term = 1,
                                                  include.split = include.split)
  } else {
    stop('Unknown move')
  }
  prior.ratio <- (prior.new/prior.old)
  lik.ratio <- exp(log.lik.new - log.lik.old)
  trans.ratio <- (prob.new.to.old/prob.old.to.new)
  acc.prob <- prior.ratio*lik.ratio*trans.ratio 
  return(list(prior.ratio = prior.ratio, lik.ratio = lik.ratio, 
              trans.ratio = trans.ratio, alpha = min(1, acc.prob)))
}


acceptance.prob.list <- function(move_list, old.tree, X, Y, sigma, prior_input_list, 
                                 include.split,
                                 cont.unif = TRUE){
  if(identical(move_list$tree, old.tree)){
    return(list(prior.ratio = 1, lik.ratio = 1, 
                trans.ratio = 1, alpha = 1))
  }
  
  prior.tree.old <- prior_input_list$fun(old.tree, prior_input_list$param[1], prior_input_list$param[1])
  prior.old <- prior.tree.old*prior.split.rule(old.tree, X, cont.unif = cont.unif)
  
  prior.tree.new <- prior_input_list$fun(move_list$tree, prior_input_list$param[1], prior_input_list$param[1])
  prior.new <- prior.tree.new*prior.split.rule(move_list$tree, X, cont.unif = cont.unif)
  # likelihood 
  log.lik.old <- cart_log_lik(tree_top = old.tree, Y = Y, X = X, sigma = sigma)
  log.lik.new <- cart_log_lik(tree_top = move_list$tree, Y = Y, X = X, sigma = sigma)
  
  # transition probabilities
  if(move_list$move == 'swap'){
    prob.old.to.new <- transition.p.swap(old.tree)
    prob.new.to.old <- transition.p.swap(move_list$tree)
  } else if(move_list$move == 'change'){
    prob.old.to.new <- transition.p.change(move_list)
    prob.new.to.old <- transition.p.change.reverse(move_list, X, old.tree)
  } else if(move_list$move == 'grow'){
    prob.old.to.new <- transition.p.grow(move_list, include.split = include.split)
    prob.new.to.old <- transition.p.grow.reverse(move_list)
  } else if(move_list$move == 'prune'){
    prob.old.to.new <- transition.p.prune(move_list)
    prob.new.to.old <- transition.p.prune.reverse(move_list, old.tree, X, 
                                                  include.split = include.split,
                                                  obs.per.term = 1)
  } else {
    stop('Unknown move')
  }
  prior.ratio <- (prior.new/prior.old)
  lik.ratio <- exp(log.lik.new - log.lik.old)
  trans.ratio <- (prob.new.to.old/prob.old.to.new)
  acc.prob <- prior.ratio*lik.ratio*trans.ratio 
  return(list(prior.ratio = prior.ratio, lik.ratio = lik.ratio, 
              trans.ratio = trans.ratio, alpha = min(1, acc.prob)))
}



acceptance.prob.list_binary <- function(move_list, old.tree, X, Y, prior_input_list, 
                                        include.split, cont.unif = TRUE){
  if(identical(move_list$tree, old.tree)){
    return(list(prior.ratio = 1, lik.ratio = 1, 
                trans.ratio = 1, alpha = 1))
  }
  
  prior.tree.old <- prior_input_list$fun(old.tree, prior_input_list$param[1], prior_input_list$param[1])
  prior.old <- prior.tree.old*prior.split.rule(old.tree, X, cont.unif = cont.unif)
  
  prior.tree.new <- prior_input_list$fun(move_list$tree, prior_input_list$param[1], prior_input_list$param[1])
  prior.new <- prior.tree.new*prior.split.rule(move_list$tree, X, cont.unif = cont.unif)
  # likelihood 
  log.lik.old <- cart_log_lik_binary(tree_top = old.tree, Y = Y, X = X)
  log.lik.new <- cart_log_lik_binary(tree_top = move_list$tree, Y = Y, X = X)
  
  #transition probabilities
  if(move_list$move == 'swap'){
    prob.old.to.new <- 1#transition.p.swap(old.tree)
    prob.new.to.old <- 1#transition.p.swap(move_list$tree)
  } else if(move_list$move == 'change'){
    prob.old.to.new <- 1#transition.p.change(move_list)
    prob.new.to.old <- 1#transition.p.change.reverse(move_list, X, old.tree)
  } else if(move_list$move == 'grow'){
    prob.old.to.new <- transition.p.grow(move_list, include.split = include.split)
    prob.new.to.old <- transition.p.grow.reverse(move_list)
  } else if(move_list$move == 'prune'){
    prob.old.to.new <- transition.p.prune(move_list)
    prob.new.to.old <- transition.p.prune.reverse(move_list, old.tree, X, obs.per.term = 1,
                                                  include.split = include.split)
  } else {
    stop('Unknown move')
  }
  prior.ratio <- (prior.new/prior.old)
  lik.ratio <- exp(log.lik.new - log.lik.old)
  trans.ratio <- (prob.new.to.old/prob.old.to.new)
  acc.prob <- prior.ratio*lik.ratio*trans.ratio 
  return(list(prior.ratio = prior.ratio, lik.ratio = lik.ratio, 
              trans.ratio = trans.ratio, alpha = min(1, acc.prob)))
}






tree_step <- function(move.type, old_tree, X, Y, Y.var, mu.prior.mean, mu.prior.var, 
                      prior_list, include.split,
                      obs.per.term = 1, cont.unif = TRUE){
  empty.flag = TRUE
  empty.count = 0
  if(move.type == 'swap'){
    while(empty.flag & (empty.count <= 10)){
      move.list <- swap_move(old_tree)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx[1,1],',',move.list$node.idx[1,2],'\n')
      # calculate the number of obs per term node
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}
    }
    acc.prob <- acceptance.prob.list(move_list = move.list, 
                                     old.tree = old_tree,
                                     X = X,
                                     Y = Y,
                                     sigma = Y.var, 
                                     prior_input_list = prior_list,
                                     cont.unif = cont.unif, 
                                     include.split = include.split)
    acceptance <- runif(1) <= acc.prob$alpha
  }
  if(move.type == 'change'){
    while(empty.flag& (empty.count <= 10)){
      
      move.list <- change_move(tree_top = old_tree, X = X, obs.per.term = obs.per.term, cont.unif = cont.unif)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}  
    }
    acc.prob <- acceptance.prob.list(move_list = move.list, 
                                     old.tree = old_tree,
                                     X = X,
                                     Y = Y,
                                     sigma = Y.var, 
                                     prior_input_list = prior_list,
                                     cont.unif = cont.unif, include.split = include.split)
    acceptance <- runif(1) <= acc.prob$alpha
  }
  
  if(move.type == 'grow'){
    while(empty.flag & (empty.count <= 10)){
      move.list <- grow_move(tree_top = old_tree, X = X, Y = Y, Y.var = Y.var, 
                             mu.prior.mean = mu.prior.mean, 
                             mu.prior.var = mu.prior.var,
                             obs.per.term = obs.per.term,
                             cont.unif = cont.unif)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}  
    }
    acc.prob <- acceptance.prob.list(move_list = move.list, 
                                     old.tree = old_tree,
                                     X = X,
                                     Y = Y,
                                     sigma = Y.var, 
                                     prior_input_list = prior_list,
                                     cont.unif = cont.unif, include.split = include.split)
    acceptance <- runif(1) <= acc.prob$alpha
  }
  if(move.type == 'prune'& (empty.count <= 10)){
    while(empty.flag){
      move.list <- prune_move(tree_top = old_tree, mu.prior.mean = mu.prior.mean, 
                              mu.prior.var = mu.prior.var, 
                              X = X, Y = Y, Y.var = Y.var)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}
    }
    
    acc.prob <- acceptance.prob.list(move_list = move.list, 
                                     old.tree = old_tree,
                                     X = X,
                                     Y = Y,
                                     sigma = Y.var, 
                                     prior_input_list = prior_list,
                                     cont.unif = cont.unif, include.split = include.split)
    acceptance <- runif(1) <= acc.prob$alpha
  }
  
  if(empty.count > 10){
    print('empty count exceeded')
    return('empty count exceeded')
  }
  if(acceptance){
    
    return(list(tree = move.list$tree, 
                prior.ratio = acc.prob$prior.ratio,
                lik.ratio = acc.prob$lik.ratio,
                trans.ratio = acc.prob$trans.ratio,
                acc.prob = acc.prob$alpha,
                accepted = acceptance) ) 
    
  } else{
    return(list(tree = old_tree, 
                prior.ratio = acc.prob$prior.ratio,
                lik.ratio = acc.prob$lik.ratio,
                trans.ratio = acc.prob$trans.ratio,
                acc.prob = acc.prob$alpha,
                accepted = FALSE))
  }
}


tree_step_binary <- function(move.type, old_tree, X, Y, alpha.prior, beta.prior, 
                      prior_list, cont.unif = TRUE, include.split,
                      obs.per.term = 1, empty.count.lim = 10){
  empty.flag = TRUE
  empty.count = 0
  if(move.type == 'swap'){
    while(empty.flag & (empty.count <= empty.count.lim)){
      move.list <- swap_move(old_tree)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx[1,1],',',move.list$node.idx[1,2],'\n')
      # calculate the number of obs per term node
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}
    }
    
  }
  if(move.type == 'change'){
    while(empty.flag& (empty.count <= empty.count.lim)){
      
      move.list <- change_move(tree_top = old_tree, X = X, obs.per.term = obs.per.term, 
                               cont.unif = cont.unif)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}  
    }
  }
  
  if(move.type == 'grow'){
    while(empty.flag& (empty.count <= empty.count.lim)){
      move.list <- grow_move_binary(tree_top = old_tree, X = X, Y = Y, 
                                    alpha.prior = alpha.prior,
                                    beta.prior = beta.prior,
                                    obs.per.term = obs.per.term,
                                    cont.unif = cont.unif)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}  
    }
  }
  if(move.type == 'prune'& (empty.count <= empty.count.lim)){
    while(empty.flag){
      move.list <- prune_move_binary(tree_top = old_tree, 
                                     alpha.prior = alpha.prior,
                                     beta.prior = beta.prior, 
                                     X = X, Y = Y)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}
    }
    
  }
  
  if(empty.count > empty.count.lim){
    print('empty count exceeded')
    return('empty count exceeded')
  }
  
  acc.prob <- acceptance.prob.list_binary(move_list = move.list, 
                                          old.tree = old_tree,
                                          X = X,
                                          Y = Y,
                                          prior_input_list = prior_list,
                                          include.split = include.split)
  acceptance <- runif(1) <= acc.prob$alpha
  
  if(acceptance){
    
    return(list(tree = move.list$tree, 
                prior.ratio = acc.prob$prior.ratio,
                lik.ratio = acc.prob$lik.ratio,
                trans.ratio = acc.prob$trans.ratio,
                acc.prob = acc.prob$alpha,
                accepted = acceptance) ) 
    
  } else{
    return(list(tree = old_tree, 
                prior.ratio = acc.prob$prior.ratio,
                lik.ratio = acc.prob$lik.ratio,
                trans.ratio = acc.prob$trans.ratio,
                acc.prob = acc.prob$alpha,
                accepted = FALSE))
  }
}


MCMC_known_var <- function(n.iter, X, Y, Y.var, mu.prior.mean, mu.prior.var, prior_list,
                           moves.prob = NULL, starting.tree = NULL, include.split,
                           cont.unif = TRUE, verbose = TRUE){
  if(is.null(starting.tree)){
    rt_old <- generate_random_binary_tree_depth_free(1)
    rt_old <- assign_node_idx(rt_old)
    rt_old <- assign_split_rules(rt_old, X)
    rt_old <- assign_term_node_values(rt_old, mu.prior.mean, mu.prior.var)
  } else {
    rt_old <- starting.tree
  }
  tree_list <- list()
  matrix.res <- matrix(NA, nrow = n.iter, ncol = 8)
  colnames(matrix.res) <- c('move', 'old.depth', 'old.nterm', 'prior.ratio', 
                            'lik.ratio', 'trans.ratio', 'acc.prob', 'acceptance')
  for(i in 1:n.iter){
    if(verbose){cat(i,'\n')}
    #debugonce(tree_step)
    old_depth <- get_depth(rt_old)
    if(old_depth == 0){
      move.type <- 'grow'
    } else if(old_depth == 1){
      move.type <- sample(c('change', 'grow'), size = 1, prob = c(0.8, 0.2))
    } else{
      if(is.null(moves.prob)){
        move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                            prob = c(0.4, 0.4, 0.1, 0.1))  
      } else{
        move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                            prob = moves.prob) 
      }
    }
    if(verbose){cat(move.type, '\n')}
    move.flag = TRUE
    while(move.flag){
      new.tree.list <- tree_step(move.type = move.type, 
                                 old_tree = rt_old, 
                                 X = X, Y = Y, Y.var = Y.var, 
                                 mu.prior.mean = mu.prior.mean, 
                                 mu.prior.var = mu.prior.var, 
                                 prior_list = prior_list,
                                 cont.unif = cont.unif,
                                 include.split = include.split)
      if(is.list(new.tree.list)){
        move.flag = FALSE} else{
          
          if(old_depth == 0){
            move.type <- 'grow'
          } else if(old_depth == 1){
            move.type <- sample(c('change', 'grow'), size = 1, prob = c(0.8, 0.2))
          } else{
            if(is.null(moves.prob)){
              move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                                  prob = c(0.4, 0.4, 0.1, 0.1))  
            } else{
              move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                                  prob = moves.prob) 
            }
          } 
          if(verbose){cat('new move: ', move.type, '\n')}  
        }
      }
    
    
    new.tree <- new.tree.list$tree
    new.tree <- assign_term_node_values_cond(tree_top = new.tree, 
                                             mu.prior.mean = mu.prior.mean, 
                                             mu.prior.var = mu.prior.var, 
                                             X = X, Y = Y, Y.var = Y.var)
    
    
    tree_list[[i]] <- new.tree
    matrix.res[i,] <- c(move.type, get_depth(rt_old), get_num_terminal_nodes(rt_old),
                        new.tree.list$prior.ratio, new.tree.list$lik.ratio, 
                        new.tree.list$trans.ratio, new.tree.list$acc.prob, 
                        new.tree.list$accepted)
    rt_old <- new.tree
  }
  return(list(trees = tree_list, df.res = data.frame(matrix.res)))
}



MCMC_binary <- function(n.iter, X, Y, alpha.prior, beta.prior, prior_list,
                           moves.prob = NULL, starting.tree = NULL, 
                        diag = FALSE, cont.unif = TRUE, include.split){
  if(is.null(starting.tree)){
    rt_old <- generate_random_binary_tree_depth_free(1)
    rt_old <- assign_node_idx(rt_old)
    rt_old <- assign_split_rules(rt_old, X)
    rt_old <- assign_term_node_values_binary(rt_old, alpha.prior, beta.prior)
  } else {
    rt_old <- starting.tree
  }
  tree_list <- list()
  matrix.res <- matrix(NA, nrow = n.iter, ncol = 8)
  colnames(matrix.res) <- c('move', 'old.depth', 'old.nterm', 'prior.ratio', 
                            'lik.ratio', 'trans.ratio', 'acc.prob', 'acceptance')
  for(i in 1:n.iter){
    cat(i,'\n')
    #debugonce(tree_step)
    old_depth <- get_depth(rt_old)
    if(old_depth == 0){
      move.type <- 'grow'
    } else if(old_depth == 1){
      move.type <- sample(c('change', 'grow'), size = 1, prob = c(0.8, 0.2))
    } else{
      if(is.null(moves.prob)){
        move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                            prob = c(0.4, 0.4, 0.1, 0.1))  
      } else{
        move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                            prob = moves.prob) 
      }
    }
    cat(move.type, '\n')
    move.flag = TRUE
    while(move.flag){
      #calculate obs at nodes with last tree
      if(diag == TRUE){
        nobs_diag <- vapply(get_terminal_nodes_idx(rt_old),
                            \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                      tree_top = rt_old, X.orig = X)), 0)
        cat('Before move nobs at nodes:', nobs_diag,'\n')
      }
       
      new.tree.list <- tree_step_binary(move.type = move.type, 
                                 old_tree = rt_old, 
                                 X = X, Y = Y, 
                                 alpha.prior = alpha.prior, 
                                 beta.prior = beta.prior, 
                                 prior_list = prior_list, 
                                 cont.unif = cont.unif, 
                                 include.split = include.split)
      cat('is tree char?' , is.character(new.tree.list), '\n')
      if(is.list(new.tree.list)){
        move.flag = FALSE} else{
          if(old_depth == 0){
            move.type <- 'grow'
          } else if(old_depth == 1){
            move.type <- sample(c('change', 'grow'), size = 1, prob = c(0.8, 0.2))
          } else{
            if(is.null(moves.prob)){
              move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                                  prob = c(0.4, 0.4, 0.1, 0.1))  
            } else{
              move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                                  prob = moves.prob) 
            }
          } 
          cat('new move: ', move.type, '\n')  
        }
    }
    
    
    new.tree <- new.tree.list$tree
    new.tree <- assign_term_node_values_cond_binary(tree_top = new.tree, 
                                             alpha.prior = alpha.prior,
                                             beta.prior = beta.prior,
                                             X = X, Y = Y)
    
    
    tree_list[[i]] <- new.tree
    matrix.res[i,] <- c(move.type, get_depth(rt_old), get_num_terminal_nodes(rt_old),
                        new.tree.list$prior.ratio, new.tree.list$lik.ratio, 
                        new.tree.list$trans.ratio, new.tree.list$acc.prob, 
                        new.tree.list$accepted)
    rt_old <- new.tree
  }
  return(list(trees = tree_list, df.res = data.frame(matrix.res)))
}




library(parallel)
multichain_MCMC_known_var <- function(n.iter, n.chain,
                                      X, Y, Y.var, 
                                      mu.prior.mean, mu.prior.var, 
                                      prior_list, 
                                      moves.prob = NULL, starting.tree = NULL,
                                      n.cores = 5,
                                      cont.unif = TRUE,
                                      include.split){
  chain.list <- mclapply(1:n.chain, 
                         \(x) MCMC_known_var(n.iter = n.iter, 
                                             X = X, 
                                             Y = Y, 
                                             Y.var = Y.var, 
                                             mu.prior.mean = mu.prior.mean, 
                                             mu.prior.var = mu.prior.var, 
                                             prior_list = prior_list, 
                                             moves.prob = moves.prob, 
                                             starting.tree = starting.tree,
                                             cont.unif = cont.unif,
                                             include.split = include.split),
                         mc.cores = n.cores)
  tree_list <- lapply(chain.list, function(x) x$trees)
  tree_list_comb <- Reduce(c, tree_list)
  
  df_list <- lapply(1:n.chain, 
                    function(x) chain.list[[x]]$df.res %>% mutate(chain = x))
  df.res <- Reduce(rbind, df_list)
  return(list(trees = tree_list_comb, df.res = df.res))
}


multichain_MCMC_binary <- function(n.iter, n.chain,
                                   X, Y, 
                                   alpha.prior, beta.prior, 
                                   prior_list, 
                                   moves.prob = NULL, starting.tree = NULL,
                                   n.cores = 5, cont.unif = TRUE,
                                   include.split = include.split){
  chain.list <- mclapply(1:n.chain, 
                         \(x) MCMC_binary(n.iter = n.iter, 
                                          X = X, 
                                          Y = Y, 
                                          alpha.prior= alpha.prior,
                                          beta.prior = beta.prior, 
                                          prior_list = prior_list, 
                                          moves.prob = moves.prob, 
                                          starting.tree = starting.tree,
                                          cont.unif = cont.unif,
                                          include.split = include.split),
                         mc.cores = n.cores)
  tree_list <- lapply(chain.list, function(x) x$trees)
  tree_list_comb <- Reduce(c, tree_list)
  
  df_list <- lapply(1:n.chain, 
                    function(x) chain.list[[x]]$df.res %>% mutate(chain = x))
  df.res <- Reduce(rbind, df_list)
  return(list(trees = tree_list_comb, df.res = df.res))
}


