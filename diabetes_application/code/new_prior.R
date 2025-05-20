f.odd <- function(n){
  (1 - (-1)^n )/2
}
f.even <- function(n){
  1 - f.odd(n)
}

## if the number of terminal nodes is even
n.cost.even <- function(gam, n){
  num_ = exp(gam*n) - 1
  den_ = exp(gam*n)*(1 - exp(-2*gam))
  num_/den_
}

n.cost.odd <- function(gam, n){
  num_ <- exp(-gam)*(1 - exp(-gam*(n-1)))
  den_ <- (1 - exp(-2*gam))
  num_/den_
}

prior.nterm <- function(n, omeg){
  exp(-omeg*n)*(exp(omeg) - 1)
}

count.tree <- function(n, delt){
  nl = (n + delt)/2
  nr = (n - delt)/2
  c.nl = factorial(2*(nl - 1))/(factorial(nl)*factorial(nl-1))
  c.nr = factorial(2*(nr - 1))/(factorial(nr)*factorial(nr-1))
  if(delt == 0){
    return(c.nl*c.nr)
  } else {
    return(2*c.nl*c.nr)  
  }
}

prior.delta <- function(delt, gam, n){
  if(n == 1){
    if(delt == 1){
      return(1)
    } else{
      return(0)
    }
  }else{
    if(gam > 0){
      if(f.odd(n)){
        n.cost <- n.cost.odd(gam, n)
      } else{
        n.cost <- n.cost.even(gam, n)
      }
      return(exp(-gam*delt)/n.cost)  
    } else{
      if(f.odd(n)){
        poss.delt <- seq(1,n-2,by = 2)
      } else {
        poss.delt <- seq(0, n-2, by = 2)
      }
      return(rep(1/length(poss.delt), length(delt)))
    }
  }
}


joint.prior.new <- function(n, delt, omeg, gam){
  if(n == 1){
    n.tree = 1
    return(prior.nterm(n, omeg))
  } else {
    n.tree <- count.tree(n, delt)  
    return(prior.nterm(n, omeg)*prior.delta(delt, gam, n)/n.tree)
  }
}


joint.prior.new.tree <- function(tree_top, omeg, gam){
  nterm_ <- get_num_terminal_nodes(tree_top)
  delta_ <- abs(get_num_terminal_nodes(tree_top$left) - get_num_terminal_nodes(tree_top$right))
  joint.prior.new(n = nterm_, delt = delta_, omeg = omeg, gam = gam)
}

# CDF of number of terminal nodes
F_nterm <- function(n, omeg){
  1 - exp(-omeg*n)
}

# inverse of CDF of number of terminal nodes
F_nterm_inv <- function(p, omeg){
  -log(1 - p)/omeg
}

# sample the number of terminal nodes prior distribution
rnterm <- function(size, omeg){
  uu <- runif(size)
  ceiling(F_nterm_inv(uu, omeg))
}



# sample from the conditional distribution of delta
rdelta.cond <- function(size, gamm, nl){
  if(nl == 1){
    delta.v = 0
  } else{
    if(f.odd(nl)){
      delta.v <- seq(1,nl-2,by=2)
    } else{
      delta.v <- seq(0,nl-2,by=2)
    }  
  }
  if(length(delta.v) > 1){
    prob.v <- prior.delta(delta.v, gamm, nl)
    return(sample(delta.v, size = size, replace = TRUE, prob = prob.v))  
  } else{
    return(delta.v)
  }
}

# sample from the marginal of rdelta
rdelta <- function(size, omeg, gamm){
  output <- rep(NA, size)
  nl.v <- rnterm(size, omeg)
  for(i in seq_along(output)){
    output[i] <- rdelta.cond(1, gamm, nl.v[i])
  }
  return(output)
}


rdepth <- function(size, omeg, gamm){
  nl.v <- rnterm(size, omeg)
  delt.v <- vapply(nl.v, \(x) rdelta.cond(1, gamm, x), 0)
  depth.v <- vapply(1:size, \(x) get_depth(generate_random_binary_tree_n_delta(n = nl.v[x],
                                                                               delt = delt.v[x])), 0)
  return(depth.v)
}


