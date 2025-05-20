
generate.tree.chip <- function(alpha., beta.){
  Gen.list <- list(1)
  depth <- 0
  n.term <- 0
  while(TRUE){
    n.nodes <- Gen.list[[depth + 1]]
    n.split <- rbinom(n = 1, 
                      size = n.nodes, 
                      prob = alpha.*((1 + depth)^(-beta.)) )
    n.term = n.term + (n.nodes - n.split)
    if(n.split == 0){
      return(list(gens = unlist(Gen.list),
                  depth = depth,
                  n.term = n.term))
    } else {
      n.nodes.next <- 2*n.split
      depth = depth + 1
      Gen.list[[depth + 1]] <- n.nodes.next
    }
  }
}

get_depth_at_node <- function(tree_top, node.idx, depth. = 0){
  if(is.null(tree_top$right) & is.null(tree_top$left)){
    if(tree_top$node.idx == node.idx){
      return(depth.)
    }
    else{
      if(tree_top$node.idx == node.idx){
        return(depth.)
      } else{
        return(NULL)  
      }
    }
  }
  else{
    if(tree_top$node.idx == node.idx){
      return(depth.)}
    depth. = depth. + 1
    r.depth = get_depth_at_node(tree_top$right, node.idx, depth.)
    l.depth = get_depth_at_node(tree_top$left, node.idx, depth.)
  }
  if(is.null(r.depth)){
    return(l.depth)
  } else{
    return(r.depth)
  }
}


chipman_prior_tree <- function(tree_top, alpha., beta.){
  # get set of internal and terminal node indexes
  term.node <- get_terminal_nodes_idx(tree_top)
  int.node <- get_internal_nodes_idx(tree_top)
  node.idx <- sort(c(term.node, int.node))
  # find depth of each node
  depth.node <- vapply(node.idx, \(x) get_depth_at_node(tree_top, x), 0)
  df.p <- data.frame(idx = node.idx, depth = depth.node)
  # logical if a node is a split
  df.p$is.int <- df.p$idx %in% int.node
  # initialize probability vector
  p.prob <- rep(NA, nrow(df.p))
  for(i in 1:nrow(df.p)){
    if(df.p$is.int[i]){
      # if node is a split
      p.prob[i] = alpha.*( (1 + df.p$depth[i])^(-beta.) )
    } else{
      # if node is not a split
      p.prob[i] = 1 - alpha.*( (1 + df.p$depth[i])^(-beta.) )
    }
  }
  # return product
  prod(p.prob)
  
}


rdepth.chip <- function(n, alpha., beta.){
  vapply(1:n, function(x) generate.tree.chip(alpha., beta.)$depth, 0) 
}

rnterm.chip <- function(n, alpha., beta.){
  vapply(1:n, function(x) generate.tree.chip(alpha., beta.)$n.term, 0) 
}

rdepth.nterm.chip <- function(n, alpha., beta.){
  out <- matrix(NA, ncol = 2, nrow = n)
  for(i in 1:n){
    gen.tree <- generate.tree.chip(alpha., beta.)
    out[i,] <- c(gen.tree$depth, gen.tree$n.term)
  }
  return(data.frame(depth = out[,1],
                    n.term = out[,2]))
}
