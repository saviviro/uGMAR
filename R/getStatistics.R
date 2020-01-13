
get_mu <- function(p, M, pars, parametrization) {
  if(parametrization == "mean") {
    mu <- pars[1,]
   # pars[1,] <<- vapply(1:M, function(i1) mu[i1]*(1 - sum(pars[2:(p + 1), i1])), numeric(1)) # Insert phi_0 to 'pars' in parent environment
    tmp <- pars
    tmp[1,] <- vapply(1:M, function(i1) mu[i1]*(1 - sum(pars[2:(p + 1), i1])), numeric(1)) # Insert phi_0 to 'pars' in parent environment
    assign("pars", tmp, parent.frame())
  } else {
    mu <- vapply(1:M, function(i1) pars[1, i1]/(1 - sum(pars[2:(p + 1), i1])), numeric(1))
  }
  mu
}
