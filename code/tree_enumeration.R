# numbers of trees with x terminal nodes
catalan_n <- function(x){
  factorial(2*(x - 1))/(factorial(x)*factorial(x - 1))
}


# function to calculate the number of possible trees with given depth and number of terminal nodes
enumerate_trees <- function(depth, nterm){
  if(depth == 0){
    return(1)
  }
  if(depth == 1){
    return(1)
  } 
  if(nterm == 2^depth){
    return(1)
  }
  if(nterm < (depth + 1) | nterm > 2^depth){
    stop('Number of terminal nodes outside (depth+1):2^depth, please provide a valid pair.')
  }
  poss.left.nodes <- max( 1, (nterm - 2^(depth - 1)) ):min( 2^(depth-1), (nterm - 1) )
  if(depth > (nterm - depth)){
    poss.left.nodes <- poss.left.nodes[poss.left.nodes <= (nterm-depth) | 
                                         poss.left.nodes >= depth]
  }
  valid.left.right.nodes <- data.frame(nl = poss.left.nodes, nr = nterm - poss.left.nodes)
  count <- 0
  for(i in 1:nrow(valid.left.right.nodes)){
    nl <- valid.left.right.nodes$nl[i]
    nr <- valid.left.right.nodes$nr[i]
    if(nl >= depth & nr >= depth){
      count = count + enumerate_trees(depth - 1, nl)*enumerate_trees(depth - 1, nr)
    } else if(nl >= depth & nr < depth){
      count = count + enumerate_trees(depth - 1, nl)*catalan_n(nr)
    } else if(nl < depth & nr >= depth){
      count = count + catalan_n(nl)*enumerate_trees(depth - 1, nr)
    }
  }
  return(count)
}
