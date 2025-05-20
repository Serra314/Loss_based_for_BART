nterm_BART <- function(bart_tree_list){
  ll <- lapply(1:length(bart_tree_list$trees), \(idx) data.frame(idx = idx,
                                                                 nn = vapply(bart_tree_list$trees[[idx]], 
                                                                             \(x) get_num_terminal_nodes(x),
                                                                             0),
                                                                 trees = factor( 1:length(bart_tree_list$trees[[idx]]))))
  Reduce(rbind, ll)
}

depth_BART <- function(bart_tree_list){
  ll <- lapply(1:length(bart_tree_list$trees), \(idx) data.frame(idx = idx,
                                                                 nn = vapply(bart_tree_list$trees[[idx]], 
                                                                             \(x) get_depth(x),
                                                                             0),
                                                                 trees = factor( 1:length(bart_tree_list$trees[[idx]]))))
  Reduce(rbind, ll)
}

BART_calculate_pred <- function(tree_list, X){
  gx <- lapply(tree_list, \(x) get_value_tree(x, X))
  Reduce('+', gx)
}



binary_log_lik_from_pred <- function(pred, Y){
  prob_ <- pnorm(pred)
  sum(Y*log(prob_) + (1 - Y)*log(1 - prob_))
}

accuracy_from_pred <- function(pred, Y){
  prob_ <- pnorm(pred)
  point_pred <- as.numeric(prob_ > 0.5)
  mean(Y == point_pred)
}