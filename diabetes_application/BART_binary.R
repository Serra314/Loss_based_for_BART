# BART_log_lik_probit <- function(tree_top, Y, X){
#   tree.at.obs <- get_value_tree(tree_top, X)
#   prob.at.obs <- pnorm(tree.at.obs)
#   log.prob.obs <- Y*log(prob.at.obs) + (1 - Y)*log(1 - prob.at.obs)
#   return(sum(log.prob.obs))
# }

BART_log_lik_probit <- function(tree_top, Y, X){
  tree.at.obs <- get_value_tree(tree_top, X)
  log.prob.obs <- dnorm(Y, mean = tree.at.obs, sd = 1, log = TRUE)
  return(sum(log.prob.obs))
}


sample_post_beta_probit <- function(Y, mean.prior, var.prior){
  nobs <- length(Y)
  Y.mean <- mean(Y)
  Y.var <- sum((Y - Y.mean)^2)/(nobs - 1)
  mu.cond.var = 1/(nobs + 1/var.prior)
  mu.cond.mean = sum(Y)*mu.cond.var
  rnorm(1, mean = mu.cond.mean, sd = sqrt(mu.cond.var))
} 

sample_post_beta_probit_internal <- function(tree_top = NULL, 
                                             node.idx = NULL, 
                                             mean.prior, 
                                             var.prior, 
                                             X = NULL, 
                                             Y = NULL, 
                                             Y.at.node = NULL){
  if(is.null(Y.at.node)){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    Y.at.node <- Y[as.numeric(rownames(obs.at.node))]
  } 
  nobs.at.node <- length(Y.at.node)
  Y.mean <- mean(Y.at.node)
  Y.var <- sum((Y.at.node - Y.mean)^2)/(nobs.at.node - 1)
  mu.cond.var = 1/(nobs.at.node + 1/var.prior)
  mu.cond.mean = sum(Y.at.node)*mu.cond.var
  return(rnorm(1, mean = mu.cond.mean, sd = sqrt(mu.cond.var)))
}



set_term_node_value_probit_prior <- function(node.idx, tree_top, mean.prior, var.prior){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = rnorm(1, mean.prior, sqrt(var.prior))
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_probit_prior(node.idx, tree_top$left, mean.prior, var.prior)
    tree.right <- set_term_node_value_probit_prior(node.idx, tree_top$right, mean.prior, var.prior)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

set_term_node_value_cond_probit <- function(node.idx, tree_top, mean.prior, var.prior, 
                                            X, Y, 
                                            Y.at.node=NULL, binary = FALSE){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = sample_post_beta_probit_internal(tree_top = tree_top, 
                                                        node.idx = node.idx, 
                                                        mean.prior = mean.prior,
                                                        var.prior = var.prior,
                                                        X = X, 
                                                        Y = Y, Y.at.node = Y.at.node)
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_cond_probit(node.idx, tree_top$left, mean.prior, var.prior, X, Y, 
                                                 Y.at.node)
    tree.right <- set_term_node_value_cond_probit(node.idx, tree_top$right, mean.prior, var.prior, X, Y, 
                                                  Y.at.node)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}


assign_term_node_values_probit_prior <- function(tree_top, mean.prior, var.prior){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    tree_top <- set_term_node_value(node.idx, tree_top, mean.prior, var.prior)
  }
  return(tree_top)
}

assign_term_node_values_cond_probit <- function(tree_top, mean.prior, var.prior, X, Y, Y.at.node = NULL){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    Y.at.node <- Y[as.numeric(rownames(obs.at.node))]
    tree_top <- set_term_node_value_cond_probit(node.idx = node.idx, tree_top = tree_top, 
                                                mean.prior = mean.prior, var.prior = var.prior, 
                                                X = X, Y =Y , Y.at.node = Y.at.node)
  }
  return(tree_top)
}


# grow move for binary observations
grow_move_probit <- function(tree_top, X, Y, mean.prior, var.prior, cont.unif = TRUE, obs.per.term = 1){
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
  # set terminal node value
  term.node.value.left <- sample_post_beta_probit(Y = Y.at.left,
                                                  mean.prior = mean.prior,
                                                  var.prior = var.prior)
  term.node.value.right <- sample_post_beta_probit(Y = Y.at.right,
                                                  mean.prior = mean.prior,
                                                  var.prior = var.prior)

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


# prune move for binary observations
prune_move_probit <- function(tree_top, mean.prior, var.prior, X, Y){
  prune.node.idx <- get_prune_idx(tree_top)
  if(length(prune.node.idx) == 1){
    prune.idx <- prune.node.idx
  } else{
    prune.idx <- sample(prune.node.idx, 1, replace = FALSE)
  }
  prune.value <- sample_post_beta_probit_internal(tree_top = tree_top, 
                                                  node.idx = prune.idx, 
                                                  mean.prior = mean.prior,
                                                  var.prior = var.prior,
                                                  X = X, 
                                                  Y = Y)
  tree_top_prune <- prune_terminal(node.idx = prune.idx, 
                                   tree_top = tree_top, 
                                   value.node = prune.value)
  tree_top_prune <- assign_node_idx(tree_top_prune)
  return(list(tree = tree_top_prune, move = 'prune', node.idx = prune.idx,
              prune.node.idx = prune.node.idx))
}





acceptance.prob.list_probit <- function(move_list, old.tree, X, Y, prior_input_list, 
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
  log.lik.old <- BART_log_lik_probit(tree_top = old.tree, Y = Y, X = X)
  log.lik.new <- BART_log_lik_probit(tree_top = move_list$tree, Y = Y, X = X) 
  #cat('old : ', log.lik.old, '\n')
  #cat('new : ', log.lik.new, '\n')
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




tree_step_probit <- function(move.type, old_tree, X, Y, mean.prior, var.prior, 
                             prior_list, cont.unif = TRUE, include.split,
                             obs.per.term = 1, empty.count.lim = 10){
  empty.flag = TRUE
  empty.count = 0
  if(move.type == 'swap'){
    while(empty.flag & (empty.count <= empty.count.lim)){
      move.list <- swap_move(old_tree)
      #cat('move=',move.list$move, ',idx = ',move.list$node.idx[1,1],',',move.list$node.idx[1,2],'\n')
      # calculate the number of obs per term node
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        #cat('empty count: ', empty.count, '\n') }
        }
      else{empty.flag = FALSE}
    }
    
  }
  if(move.type == 'change'){
    while(empty.flag& (empty.count <= empty.count.lim)){
      
      move.list <- change_move(tree_top = old_tree, X = X, obs.per.term = obs.per.term, 
                               cont.unif = cont.unif)
      #cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        #cat('empty count: ', empty.count, '\n') }
      }
      else{empty.flag = FALSE}  
    }
  }
  
  if(move.type == 'grow'){
    while(empty.flag& (empty.count <= empty.count.lim)){
      move.list <- grow_move_probit(tree_top = old_tree, X = X, Y = Y, 
                                    mean.prior = mean.prior,
                                    var.prior = var.prior,
                                    obs.per.term = obs.per.term,
                                    cont.unif = cont.unif)
      #cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        #cat('empty count: ', empty.count, '\n') }
      }
      else{empty.flag = FALSE}  
    }
  }
  if(move.type == 'prune'& (empty.count <= empty.count.lim)){
    while(empty.flag){
      move.list <- prune_move_probit(tree_top = old_tree, 
                                     mean.prior = mean.prior,
                                     var.prior = var.prior, 
                                     X = X,
                                     Y = Y)
      #cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        #cat('empty count: ', empty.count, '\n') }
      }
      else{empty.flag = FALSE}
    }
    
  }
  
  if(empty.count > empty.count.lim){
    #print('empty count exceeded')
    return('empty count exceeded')
  }
  
  acc.prob <- acceptance.prob.list_probit(move_list = move.list, 
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


sample_Z_given_Y <- function(Y, tree_at_obs){
  nobs <- length(Y)
  zz <- rep(0, nobs)
  for(idx in 1:nobs){
    if(Y[idx] == 0){
      zz[idx] <- rtruncnorm(n = 1, a = -Inf, b = 0, mean = tree_at_obs[idx], sd = 1)
    } else {
      zz[idx] <- rtruncnorm(n = 1, a = 0, b = Inf, mean = tree_at_obs[idx], sd = 1)
    }
  }
  return(zz)
}





BART_single_step <- function(X, res_Y, mean.prior, var.prior, prior_list,
                              moves.prob = NULL, starting.tree = NULL, 
                              diag = FALSE, cont.unif = TRUE, include.split){
  if(is.null(starting.tree)){
    rt_old <- generate_random_binary_tree_depth_free(1)
    rt_old <- assign_node_idx(rt_old)
    rt_old <- assign_split_rules(rt_old, X)
    rt_old <- assign_term_node_values_probit_prior(rt_old, mean.prior, var.prior)
  } else {
    rt_old <- starting.tree
  }
  
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
  #cat(move.type, '\n')
  move.flag = TRUE
  while(move.flag){
    #calculate obs at nodes with last tree
    if(diag == TRUE){
      nobs_diag <- vapply(get_terminal_nodes_idx(rt_old),
                          \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                    tree_top = rt_old, X.orig = X)), 0)
      #cat('Before move nobs at nodes:', nobs_diag,'\n')
    }
    
    # calculate residual
    #tree_at_obs <- get_value_tree(rt_old, X)
    #zz <- sample_Z_given_Y(Y, tree_at_obs)
    #res_ <- zz - tree_at_obs
    new.tree.list <- tree_step_probit(move.type = move.type, 
                                      old_tree = rt_old, 
                                      X = X, 
                                      Y = res_Y, 
                                      mean.prior = mean.prior, 
                                      var.prior = var.prior, 
                                      prior_list = prior_list, 
                                      cont.unif = cont.unif, 
                                      include.split = include.split)
    #cat('is tree char?' , is.character(new.tree.list), '\n')
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
        #cat('new move: ', move.type, '\n')  
      }
  }
  
  
  new.tree <- new.tree.list$tree
  new.tree <- assign_term_node_values_cond_probit(tree_top = new.tree, 
                                                  mean.prior = mean.prior,
                                                  var.prior = var.prior,
                                                  X = X, Y = res_Y)
  
  
  new.tree
  
}



MCMC_probit<- function(n.iter, X, Y, mean.prior, var.prior, prior_list,
                        moves.prob = NULL, starting.tree = NULL, 
                        diag = FALSE, cont.unif = TRUE, include.split){
  if(is.null(starting.tree)){
    rt_old <- generate_random_binary_tree_depth_free(1)
    rt_old <- assign_node_idx(rt_old)
    rt_old <- assign_split_rules(rt_old, X)
    rt_old <- assign_term_node_values_probit_prior(rt_old, mean.prior, var.prior)
  } else {
    rt_old <- starting.tree
  }
  tree_list <- list()
  matrix.res <- matrix(NA, nrow = n.iter, ncol = 8)
  colnames(matrix.res) <- c('move', 'old.depth', 'old.nterm', 'prior.ratio', 
                            'lik.ratio', 'trans.ratio', 'acc.prob', 'acceptance')
  # initialise to calculate Z
  idx.y0 <- yy == 0
  idx.y1 <- yy == 1
  for(i in 1:n.iter){
    #cat(i,'\n')
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
    #cat(move.type, '\n')
    move.flag = TRUE
    while(move.flag){
      #calculate obs at nodes with last tree
      if(diag == TRUE){
        nobs_diag <- vapply(get_terminal_nodes_idx(rt_old),
                            \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                      tree_top = rt_old, X.orig = X)), 0)
        #cat('Before move nobs at nodes:', nobs_diag,'\n')
      }
      
      # calculate residual
      tree_at_obs <- get_value_tree(rt_old, X)
      zz <- sample_Z_given_Y(Y, tree_at_obs)
      res_ <- zz - tree_at_obs
      new.tree.list <- tree_step_probit(move.type = move.type, 
                                        old_tree = rt_old, 
                                        X = X, 
                                        Y = res_, 
                                        mean.prior = mean.prior, 
                                        var.prior = var.prior, 
                                        prior_list = prior_list, 
                                        cont.unif = cont.unif, 
                                        include.split = include.split)
      #cat('is tree char?' , is.character(new.tree.list), '\n')
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
          #cat('new move: ', move.type, '\n')  
        }
    }
    
    
    new.tree <- new.tree.list$tree
    new.tree <- assign_term_node_values_cond_probit(tree_top = new.tree, 
                                                    mean.prior = mean.prior,
                                                    var.prior = var.prior,
                                                    X = X, Y = res_)
    
    
    tree_list[[i]] <- new.tree
    matrix.res[i,] <- c(move.type, get_depth(rt_old), get_num_terminal_nodes(rt_old),
                        new.tree.list$prior.ratio, new.tree.list$lik.ratio, 
                        new.tree.list$trans.ratio, new.tree.list$acc.prob, 
                        new.tree.list$accepted)
    rt_old <- new.tree
  }
  return(list(trees = tree_list, df.res = data.frame(matrix.res)))
}




BART_calculate_pred <- function(tree_list, X){
  gx <- lapply(tree_list, \(x) get_value_tree(x, X))
  Reduce('+', gx)
}

BART_log_lik_probit_binary <- function(tree_list, Y, X){
  G_x <- BART_calculate_pred(tree_list, X)
  probs_at_obs <- pnorm(G_x)
  log_prob_at_obs <- Y*log(probs_at_obs) + (1 - Y)*log(1 - probs_at_obs)
  sum(log_prob_at_obs)
}

Acceptance_probit <- function(tree_list, Y, X){
  G_x <- BART_calculate_pred(tree_list, X)
  probs_at_obs <- pnorm(G_x)
  pred_at_obs <- as.numeric(probs_at_obs > 0.5)
  mean(Y == pred_at_obs)
}

probs_y01_probit <- function(tree_list, Y, X, alpha.p = 0.1){
  G_x <- BART_calculate_pred(tree_list, X)
  probs_at_obs <- pnorm(G_x)
  p0 <- probs_at_obs[Y == 0]
  p1 <- probs_at_obs[Y == 1]
  df0 <- data.frame(av.prob = mean(p0),
                    q.low = quantile(p0, alpha.p/2),
                    median = median(p0),
                    q.up = quantile(p0, 1 - alpha.p/2),
                    Y = '0')
  df1 <- data.frame(av.prob = mean(p1),
                    q.low = quantile(p1, alpha.p/2),
                    median = median(p1),
                    q.up = quantile(p1, 1 - alpha.p/2),
                    Y = '1')
  rbind(df0, df1)
}




MCMC_for_BART_binary <- function(n.tree =10,
                                 n.iter = 500, 
                                 Y, X, 
                                 include.split = FALSE,
                                 cont.unif = TRUE,
                                 mean.prior = 0,
                                 var.prior = 3,
                                 moves.prob = c(0.4, 0.4, 0.1, 0.1),
                                 prior.list = list(fun = joint.prior.new.tree, param = c(2, 0.62)),
                                 verbose = TRUE){
  idx.y0 = Y == 0
  idx.y1 = Y == 1
  #zz <- rep(0,length(Y))
  #zz[idx.y0] <- rtruncnorm(n = sum(idx.y0), a = -Inf, b = 0, mean = mean.prior, sd = 1)
  #zz[idx.y1] <- rtruncnorm(n = sum(idx.y1), a = 0, b = Inf, mean = mean.prior, sd = 1)
  
  
  
  # initialise output
  starting.tree.list <- lapply(1:n.tree, \(x) NULL)
  last.tree.list <- starting.tree.list
  idx.tree.vec <- 1:n.tree
  Y_pred_list <- lapply(1:n.tree, \(x) rep(0, length(Y)))
  nterm_tree <- lapply(1:n.tree, \(x) c() )
  depth_tree <- lapply(1:n.tree, \(x) c() )
  

  
  # params for computation
  
  trees.for.iter <- list()
  
  st_time = Sys.time()
  for(idx.iter in 1:n.iter){
    cat('Iteration: ', idx.iter, '\n')
    pred_at_obs <- Reduce('+', Y_pred_list)
    zz <- sample_Z_given_Y(Y, pred_at_obs)
    
    for(idx.tree in idx.tree.vec){
      residual_Y <- zz -  Reduce('+', Y_pred_list[idx.tree.vec != idx.tree])
      if(verbose){cat('Tree: ', idx.tree, '\n')}
      old_tree <- last.tree.list[[idx.tree]]  
      new_tree <- BART_single_step(X = X, 
                                   res_Y = residual_Y, 
                                   mean.prior  = mean.prior, 
                                   var.prior = var.prior, 
                                   prior_list = prior.list, 
                                   moves.prob = moves.prob, 
                                   starting.tree = old_tree,
                                   cont.unif = cont.unif, 
                                   include.split = include.split
      )
      last.tree.list[[idx.tree]] <- new_tree
      dt = get_depth(new_tree)
      nt = get_num_terminal_nodes(new_tree)
      if(verbose){
      cat('Depth: ', dt, '\n')
      cat('Terminal Nodes: ', nt, '\n')
      }
      Y_pred <- vapply(1:nrow(X), \(x) g.T(new_tree, X[x,]), 0)
      Y_pred_list[[idx.tree]] <- Y_pred
      
      nterm_tree[[idx.tree]] <- c(nterm_tree[[idx.tree]] , nt)
      depth_tree[[idx.tree]] <- c(depth_tree[[idx.tree]] , nt)
      if(verbose){print('----------------')}
    }
    trees.for.iter[[idx.iter]] <- last.tree.list
    #average_res <- c(average_res, mean(Y - Reduce('+', Y_pred_list)))
    if(verbose){print('##############')}
  }
  
  time_diff <- Sys.time() - st_time
  print(time_diff)
  return(list(trees = trees.for.iter, 
              time = time_diff))
}




