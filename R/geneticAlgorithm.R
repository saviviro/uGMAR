#' @import stats
#'
#' @title Genetic algorithm for preliminary estimation of GMAR, StMAR or G-StMAR model
#'
#' @description \code{GAfit} estimates specified GMAR, StMAR or G-StMAR model using genetic algorithm. It's designed to find starting values for gradient based methods.
#'
#' @inheritParams loglikelihood_int
#' @param ngen an (optional) positive integer specifying the number of generations to be ran through in the genetic algorithm. Default is \code{min(400, max(round(0.1*length(data)), 200))}.
#' @param popsize an (optional) positive even integer specifying the population size in the genetic algorithm. Default is \code{10*d} where \code{d} is the number of parameters.
#' @param smartMu an (optional) positive integer specifying the generation after which the random mutations in the genetic algorithm are "smart".
#'  This means that mutating individuals will mostly mutate fairly close (or partially close) to the best fitting individual so far. Default is \code{min(100, round(0.5*ngen))}.
#' @param ar0scale an (optional) real valued vector of length two specifying the mean (the first element) and standard deviation (the second element) of the normal distribution
#'  from which the \eqn{\phi_{m,0}} parameters are generated in random mutations in the genetic algorithm. Default is \code{c(1.5*avg*(1-c1/c0), max(c0, 4))}, where
#'  avg is sample mean, \code{c1} is the first sample autocovariance and \code{c0} is sample variance.
#' @param sigmascale an (optional) positive real number specifying the standard deviation of the (zero mean, positive only) normal distribution
#'  from which the component variance parameters are generated in the random mutations in the genetic algorithm. Default is \code{1+sd(data)}.
#' @param initpop an (optional) list of parameter vectors from which the initial population of the genetic algorithm will be generated from.
#'   The parameter vectors should be of form...
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M}}).}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(M(p+3)+M2-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M}}).}
#'        \item{With \strong{linear constraints}:}{Replace the vectors \strong{\eqn{\phi_{m}}} with vectors \strong{\eqn{\psi_{m}}} and provide a  list of constraint
#'          matrices \strong{R} that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'          \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{M})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M})}.}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(3M+M2+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M})}.}
#'        \item{With \strong{linear constraints}:}{Replace the vector \strong{\eqn{\phi}} with vector \strong{\eqn{\psi}} and provide a constraint matrix
#'          \strong{\eqn{R}} that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{v} a degrees of
#'  freedom parameter.
#'  Note that in the case \strong{M=1} the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#'  If not specified (or \code{FALSE} as is default), the initial population will be drawn randomly.
#' @param minval a real number defining the minimum value of the log-likelihood function that will be considered.
#'   Values smaller than this will be treated as they were \code{minval} and the corresponding inidividuals will never survive.
#' @details The user should consider adjusting \code{ar0scale} and/or \code{sigmascale} accordingly to the best knowledge about the process.
#'
#'    The genetic algorithm is mostly based on the description by \emph{Dorsey R. E. ja Mayer W. J. (1995)}.
#'     It uses individually adaptive crossover and mutation rates described by \emph{Patnaik L.M. and Srinivas M. (1994)}, with slight modifications.
#'
#' @return Returns estimated parameter vector...
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M}}).}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(M(p+3)+M2-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M}}).}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above, but vectors \strong{\eqn{\phi_{m}}} replaced with vectors \strong{\eqn{\psi_{m}}}
#'          that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'          \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{M})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M})}.}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(3M+M2+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M})}.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above, but vector \strong{\eqn{\phi}} replaced with vector \strong{\eqn{\psi}}
#'          that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{\nu} a degrees of
#'  freedom parameter.
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. (2015) Gaussian Mixture Autoregressive model for univariate time series.
#'          \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Dorsey R. E. ja Mayer W. J. (1995) Genetic algorithms for estimation problems with multiple optima, nondifferentiability, and other irregular features.
#'          \emph{Journal of Business & Economic Statistics}, \strong{13}, 53-66.
#'    \item Patnaik L.M. and Srinivas M. (1994) Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item References regarding the StMAR and G-StMAR models will be updated after they are published.
#'  }

GAfit <- function(data, p, M, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, conditional=TRUE, ngen, popsize, smartMu, ar0scale, sigmascale, initpop=FALSE, epsilon, minval) {

  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p, M, GStMAR=GStMAR)
  data = checkAndCorrectData(data, p)
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  # Check constraint matrices
  if(constraints==TRUE) {
    checkConstraintMat(p, M, R, restricted=restricted)
    if(restricted==TRUE) {
      Rquals = TRUE # States wether all matrices R have equal dimensions
    } else {
      Rdims = sapply(1:M, function(i2) ncol(R[[i2]]))
      if(length(unique(Rdims))==1) {
        Rquals = TRUE
      } else {
        Rquals = FALSE
      }
    }
  } else {
    Rquals = TRUE # Always true if no constraints
  }
  d = nParams(p=p, M=M_orig, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)

  # Default settings
  if(missing(popsize)) {
    popsize = 10*d
  }
  if(missing(ngen)) {
    ngen = min(400, max(round(0.1*length(data)), 200))
  }
  if(missing(smartMu)) {
    smartMu =  min(100, round(0.5*ngen))
  }
  if(missing(ar0scale)) {
    avg = mean(data); T1 = length(data)
    c0 = (t(data-avg)%*%(data-avg))/T1
    c1 = (t((data-avg)[1:(T1-1)])%*%(data-avg)[2:T1])/(T1-1)
    ar0scale = c(1.5*avg*(1-c1/c0), max(c0, 4))
  }
  if(missing(sigmascale)) {
    sigmascale = 1+sd(data)
  }
  if(missing(epsilon)) {
    epsilon = round(log(.Machine$double.xmin)+10)
  }
  if(ngen<1) {
    stop("The number of generations can't be less than one")
  } else if(ngen%%1!=0) {
    stop("The number of generations has to be positive integer")
  } else if(popsize%%2!=0) {
    stop("The size of population has to be EVEN positive integer")
  } else if(popsize<2) {
    stop("Population size can't be smaller than two")
  } else if(smartMu%%1!=0) {
    stop("smartMu has to be (positive) integer")
  }

  # If initial population is not set, choose it randomly and form inital generation matrix G.
  if(!is.list(initpop)) {
    G = vapply(1:popsize, function(x) randomIndividual_int(p, M_orig, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R, ar0scale=ar0scale, sigmascale=sigmascale), numeric(d))
  } else {
    # If the initial individuals are set by the user: check it and create the initial generation matrix G.
    n_inds = length(initpop)
    for(i1 in 1:n_inds) {
      params = initpop[[i1]]
      if(length(params)!=d) {
        stop(paste("The parameter vector of individual", i1, "in the initial population is wrong dimension"))
      }
      if(constraints==TRUE) {
        params = reformConstrainedPars(p, M_orig, params, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
      }
      pars_tmp = reformParameters(p, M_orig, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
      params = pars_tmp$params
      pars = pars_tmp$pars
      alphas = pars_tmp$alphas
      if(StMAR==TRUE | GStMAR==TRUE) {
        dfs = pars_tmp$dfs
      }
      if(StMAR==TRUE | GStMAR==TRUE) {
        if(any(dfs<=2)) {
          stop(paste("The individual", i1, "in the initial population has degrees of freedom not larger than 2"))
        }
      }
      if(!isStationary_int(p, M, params, restricted=restricted)) {
        stop(paste("The individual", i1, "in the initial population is not stationary"))
      } else if(any(pars[p+2,]<=0)) {
        stop("Component variances has to be larger than zero")
      } else if(sum(alphas)>1+1e-6) {
        stop(paste("The mixing weights of the individual", i1," in the initial population don't sum to one"))
      } else if(constraints==FALSE) { # Sort components in the initial population so that alpha_1>...>alpha_M
        initpop[[i1]] = sortComponents(p, M_orig, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
      }
    }
    if(n_inds > popsize) {
      stop("Number of individuals in the initial population is larger than population size")
    } else if(n_inds==1) {
      G = do.call(cbind, replicate(popsize, initpop[[1]], simplify=FALSE) ) # Repeat initial 'individual' popsize times and form the generation matrix.
    } else {
      # Draw random population from the set of individuals provided
      G = sapply(1:popsize, function(x) initpop[[sample.int(n_inds, size=1)]])
    }
  }

  # Size (M*(p+3)-1)x(popsize) population matrix G that contains parameter vectors as columns.
  # Each generation will be stored in an array of generations. The log-likelihood values of each parameter matrix in each generation will be stored
  # in a separate size (ngen)x(popsize) matrix such that each row indicates a generation and each column indicates the individual in the generation.

  # Initial setup
  generations = array(dim=c(d, popsize, ngen))
  logliks = matrix(minval, nrow=ngen, ncol=popsize)
  whichRedundant_alt = numeric(M)
  wrMeasure_alt = 0

  ### Run through generations ###
  for(i1 in 1:ngen) {
    generations[,,i1] = G # Save the generation's population to the container.

    # Compute log-likelihoods
    if(i1>1) {
      # Fitness inheritance: individual has 50% changes to inherit fitness if it's a result from crossover.
      whichCO = which(1-whichNotCo==1) # Which individuals did crossover
      if(length(whichCO)>=1) {
        inherit = sample(whichCO, size=round(0.5*length(whichCO)), replace=FALSE) # Draw the ones who inherit fitness
      } else {
        inherit = numeric(0)
      }

      # survivorLiks holds the parents' loglikelihood values: for odd number they are (index, index+1) and for even (index-1, index)
      if(length(inherit)>=1 & abs(maxLiks-avgLiks)>abs(0.03*avgLiks)) { # If massive mutations: no fitness inheritance
        for(i2 in 1:length(inherit)) {
          inheritor = inherit[i2]
          if(inheritor %% 2 == 0) {
            logliks[i1, inheritor] = 0.5*(survivorLiks[inheritor-1] + survivorLiks[inheritor])
          } else {
            logliks[i1, inheritor] = 0.5*(survivorLiks[inheritor] + survivorLiks[inheritor+1])
          }
        }
        noinherit = (1:popsize)[-inherit]
      } else {
        noinherit = 1:popsize
      }
      for(i2 in noinherit) { # Fill in the rest log-likelihoods
        if(all(H[,i2]==H2[,i2])) {
          logliks[i1, i2] = survivorLiks[i2] # Individual is unchanged
        } else {
          logliks[i1, i2] = loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R,
                                              conditional=conditional, boundaries=TRUE, checks=FALSE, returnTerms=FALSE, epsilon=epsilon, minval=minval)
        }
      }
    } else { # The first round
      logliks[i1,] = vapply(1:popsize, function(i2) loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints,
                                                                      R=R, conditional=conditional, boundaries=TRUE, checks=FALSE, returnTerms=FALSE, epsilon=epsilon,
                                                                      minval=minval), numeric(1))
    }
    # If loglik is so bad it is less than "minval": setting "minval" will make the surviving probability zero in the next phase
    logliks[i1, which(logliks[i1,] < minval)] = minval

    ### Creating new generation start here ###

    # Create the reproduction pool
    if(length(unique(logliks[i1,]))==1) { # If all individuals are the same, the surviving probability is 1.
      surviveProbs = rep(1, popsize)
    } else {
      T_values = logliks[i1,] + abs(min(logliks[i1,])) # Function T as described by Dorsey R. E. ja Mayer W. J., 1995
      surviveProbs = T_values/sum(T_values) # The surviving probabilities
    }
    candidates = which(surviveProbs!=0)
    if(length(candidates)==1) {
      survivors = rep(candidates, popsize)
    } else {
      survivors = sample(x=candidates, size=popsize, replace=TRUE, prob=surviveProbs[candidates])
    }
    H = G[,survivors] # Save the survivors in reproduction pool H

    # Calculate avg, min and max of loglikelihoods of the survivors
    survivorLiks = logliks[i1, survivors]
    avgLiks = mean(survivorLiks)
    maxLiks = max(survivorLiks)

    # Individually adaptive cross-over rate as described by  M. Srinivas L.M. Patnaik (1994)
    if(maxLiks == avgLiks) {
      avgLiks = avgLiks + 0.5
    }
    indeces = seq(1,popsize-1,by=2)
    parentLiks = vapply(1:length(indeces), function(i2) max(survivorLiks[indeces[i2]], survivorLiks[indeces[i2]+1]), numeric(1)) # Max of parents
    coRate = vapply(1:(popsize/2), function(i2) max(min(1, (maxLiks - parentLiks[i2])/(maxLiks - avgLiks)), 0.4), numeric(1)) # The best genes should mix the the population too

    # For each pair of individuals draw from bernoulli wether crossover is to be made or not. If yes, do it.
    crossOvers = rbinom(popsize/2, 1, coRate) # 1 = do crossover, 0 = no cross-over
    I = round(runif(popsize/2, min=1, max=nrow(H)-1)) # Break points
    indeces = seq(1,popsize-1,by=2)
    H2 = as.vector(vapply(1:(popsize/2), function(i2) {
      i3 = indeces[i2]
      if(crossOvers[i2]==1) {
        c( c(H[1:I[i2], i3], H[(I[i2]+1):d, i3+1]), c(H[1:I[i2], i3+1], H[(I[i2]+1):d, i3]) ) # Crossover according to I
      } else {
        c(H[,i3], H[,i3+1])
      }
    }, numeric(2*d)))
    H2 = matrix(H2, ncol=popsize)

    # Mutation rate
    whichNotCo = as.vector(sapply(1:length(crossOvers), function(i2) rep(1 - crossOvers[i2], 2))) # Indicator who did NOT crossover
    if(abs(maxLiks - avgLiks)<=abs(0.03*avgLiks)) {
      muRate = rep(0.7, popsize) # Massive mutations if converging
    } else {
      muRate = 0.5*vapply(1:popsize, function(i2) min(whichNotCo[i2], (maxLiks - survivorLiks[i2])/(maxLiks - avgLiks)), numeric(1)) # Individually adaptive mutation rate
    }

    # Mutate
    mutate = which(rbinom(popsize, 1, muRate)==1)
    if(length(mutate)>=1) {
      # Smart mutations with mechanism trying to find all the relevant regimes if some regimes are wasted
      if(i1 > smartMu) {
        # The best individual so far and its mixing weights
        bestIndex = which(logliks == max(logliks), arr.ind = TRUE)[1,]
        bestind = generations[, bestIndex[2], bestIndex[1]]
        mw = mixingWeights_int(data=data, p=p, M=M_orig, params=bestind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                               constraints=constraints, R=R, checks=FALSE, epsilon=epsilon)

        # Which regimes are wasted
        whichRedundant = which((vapply(1:M, function(i2) sum(mw[,i2]>0.05)<0.01*length(data), logical(1))))

        # Which regimes are wasted with lighter criteria (implies stronger redundance)
        whichRedLight = which((vapply(1:M, function(i2) sum(mw[,i2]>0.0005)<0.01*length(data), logical(1))))

        # The redundant regime to be focused on
        if(length(whichRedundant)>0) {
          whichRedToUse = whichRedundant[1]
        } else {
          whichRedToUse = numeric(0)
        }
        wrMeasure = sum(mw[,whichRedToUse]>0.0005) # Measures the relevance of the chosen redundant regime

        # Mutate more if there are wasted regimes
        if(length(whichRedundant)>0) {
          muRate = vapply(1:popsize, function(i2) whichNotCo[i2]*max(0.1, muRate[i2]), numeric(1))
          mutate0 = which(rbinom(popsize, 1, muRate)==1)
          if(length(mutate0)>length(mutate)) {
            mutate = mutate0
          }
        }

        # Update the alternative individual if bestind fits better or as well for the job
        cond1 = wrMeasure >= wrMeasure_alt & length(whichRedundant)<=length(whichRedundant_alt)
        cond2 = length(whichRedundant)<length(whichRedundant_alt)
        if(cond1 | cond2) {
          bestind_alt = bestind
          whichRedundant_alt = whichRedundant
          whichRedLight_alt = whichRedLight
          wrMeasure_alt = wrMeasure
        }

        # Smart mutate the bestind if there are no redundant regimes or if alternative individual is not found,
        # or changes by random, or if considering G-StMAR model
        accuracy = abs(rnorm(length(mutate), mean=15, sd=10)) # Accuracy
        if(GStMAR==TRUE | length(whichRedundant)==0 | cond1 | cond2 | !length(whichRedLight_alt)<length(whichRedLight) | rbinom(1, 1, 0.5)==1) {
          if(rbinom(1, 1, 0.3)==1) {
            whichRandom = whichRedLight
          } else {
            whichRandom = whichRedundant
          }
          H2[,mutate] = vapply(1:length(mutate), function(x) smartIndividual_int(p=p, M=M_orig, params=bestind, StMAR=StMAR, GStMAR=GStMAR,
                                                                                 restricted=restricted, constraints=constraints, R=R,
                                                                                 ar0scale=ar0scale, sigmascale=sigmascale, accuracy=accuracy[x],
                                                                                 whichRandom=whichRandom), numeric(d))
        } else { # If alternative bestind is found: change to mutate close to a combination of bestindividual and alternative individual

          # Alternative ind has more somewhat relevant redundant regimes than the bestind: take the redundant regimes from best ind and
          # combine with altinds more relevant redundant regimes.
          if(Rquals==TRUE) {
            nonRedRegimes = (1:M)[-whichRedundant] # Non redundant regimes of bestind
            bestind_nonRedRegimes = sapply(nonRedRegimes, function(i2) extractRegime(p=p, M=M_orig, params=bestind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                                                                                     constraints=constraints, R=R, regime=i2)) # Regimes as columns
            all_altind_regimes = sapply(1:M, function(i2) extractRegime(p=p, M=M_orig, params=bestind_alt, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                                                                        constraints=constraints, R=R, regime=i2))

            # Calculate altind regimes' distances to bestinds non redundant regimes
            if(StMAR==TRUE) { # Remove dfs from distance comparison
              bestind_regimes0 = as.matrix(bestind_nonRedRegimes[1:(nrow(bestind_nonRedRegimes)-1),])
              all_altind_regimes0 = as.matrix(all_altind_regimes[1:(nrow(all_altind_regimes)-1),])
            } else {
              bestind_regimes0 = as.matrix(bestind_nonRedRegimes)
              all_altind_regimes0 = as.matrix(all_altind_regimes)
            }
            distanceToRegime = matrix(nrow=length(nonRedRegimes), ncol=M) # Row for each nonred-regime and column for each altind regime
            for(i2 in 1:length(nonRedRegimes)) {
              distanceToRegime[i2,] = vapply(1:M, function(i3) t(bestind_regimes0[,i2]-all_altind_regimes0[,i3])%*%(bestind_regimes0[,i2]-all_altind_regimes0[,i3]), numeric(1))
            }

            # Which regimes are closest to the non-redundant regimes (these are not used)
            closestRegimes = numeric(nrow(distanceToRegime))
            for(i2 in 1:nrow(distanceToRegime)) {
              tmp = which(distanceToRegime[i2,]==min(distanceToRegime[i2,]))
              closestRegimes[i2] = tmp # Which regime is the closest one
              distanceToRegime[,tmp] = rep(99999, nrow(distanceToRegime)) # "Exclude" already chosen regimes from the comparison
            }

            # Pick the best regime left from bestind_alt
            mw2 = mixingWeights_int(data=data, p=p, M=M_orig, params=bestind_alt, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                                    constraints=constraints, R=R, checks=FALSE, epsilon=epsilon) # Mixing weights
            wrMeasures0 = vapply(1:M, function(i2) sum(mw2[,i2]>0.0005), numeric(1)) # wrMeasures of bestind_alt
            orderAltindRegimes = order(wrMeasures0, decreasing=TRUE)
            bestRedRegime = all_altind_regimes[,orderAltindRegimes[which(!orderAltindRegimes %in% closestRegimes)][1]] # Biggest wrMeasure excluding the regimes closest to the bestind's non-redundant regime

            # Combine the regimes to a complete parameter vector
            combinedInd = changeRegime(p=p, M=M_orig, params=bestind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R,
                                       regimeParams=bestRedRegime, regime=whichRedToUse)
          } else {
            # If constraints==TRUE and matrices R differ in dimension things needs to be done differently because the regimes differ in length
            nonRedRegimes = (1:M)[-whichRedundant] # Non redundant regimes of bestind
            bestind_nonRedRegimes = lapply(nonRedRegimes, function(i2) extractRegime(p=p, M=M_orig, params=bestind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                                                                                     constraints=constraints, R=R, regime=i2)) # Regimes as columns
            all_altind_regimes = lapply(1:M, function(i2) extractRegime(p=p, M=M_orig, params=bestind_alt, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                                                                        constraints=constraints, R=R, regime=i2))
            whichRedToUse_pars = extractRegime(p=p, M=M_orig, params=bestind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R, regime=whichRedToUse)

            # Which altind regimes are the same length as the wasted regime
            bestRedRegime_index = which(vapply(1:M, function(i2) length(all_altind_regimes[[i2]]), numeric(1)) == length(whichRedToUse_pars))

            # If multiple choices
            if(length(bestRedRegime_index)>1) {
              # Find non redundant regimes that are the same dimension...
              whichAreSameLength = which(vapply(1:length(bestind_nonRedRegimes), function(i2) length(bestind_nonRedRegimes[[i2]]), numeric(1)) == length(whichRedToUse_pars))

              # ...and if any...
              if(length(whichAreSameLength)>0) {
                # ...remove dfs from distance comparison...
                if(StMAR==TRUE) {
                  bestind_regimes0 = lapply(1:length(bestind_nonRedRegimes), function(i2) bestind_nonRedRegimes[[i2]][-length(bestind_nonRedRegimes[[i2]])])
                  altind_regimes0 = lapply(1:M, function(i2) all_altind_regimes[[i2]][-length(all_altind_regimes[[i2]])])
                } else {
                  bestind_regimes0 = bestind_nonRedRegimes
                  altind_regimes0 = all_altind_regimes
                }
                #  ...and exclude the regimes closest to the same dimension non redundant regimes
                distanceToRegime = matrix(nrow=length(whichAreSameLength), ncol=length(bestRedRegime_index)) # Row for each nonred-regime and column for each altind regime
                j = 1
                for(i2 in whichAreSameLength) {
                  distanceToRegime[j,] = vapply(bestRedRegime_index, function(i3) t(bestind_regimes0[[i2]]-altind_regimes0[[i3]])%*%(bestind_regimes0[[i2]]-altind_regimes0[[i3]]), numeric(1))
                  j = j+1
                }
                # Which regimes are closest to the non-redundant regimes (these are not used)
                closestRegimes = numeric(nrow(distanceToRegime))
                for(i2 in 1:nrow(distanceToRegime)) {
                  tmp = which(distanceToRegime[i2,]==min(distanceToRegime[i2,]))
                  closestRegimes[i2] = tmp # Which regime is the closest one
                  distanceToRegime[,tmp] = rep(99999, nrow(distanceToRegime)) # "Exclude" already chosen regimes from the comparison
                }
                # Exclude the closest regimes from the comparison
                bestRedRegime_index = bestRedRegime_index[-closestRegimes]
              }

              # If still multiple choices left: choose the most relevant one
              if(length(bestRedRegime_index)>1) {
                mw2 = mixingWeights_int(data=data, p=p, M=M_orig, params=bestind_alt, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                                        constraints=constraints, R=R, checks=FALSE, epsilon=epsilon)
                wrMeasures0 = vapply(1:M, function(i2) sum(mw2[,i2]>0.0005), numeric(1)) # wrMeasures of bestind_alt
                bestRedRegime_index = bestRedRegime_index[which(wrMeasures0[bestRedRegime_index] == max(wrMeasures0[bestRedRegime_index]))[1]] # The final choice
              }
            }
            # Instert the chosen regime into the best fitting individual
            combinedInd = changeRegime(p=p, M=M_orig, params=bestind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R,
                                       regimeParams=all_altind_regimes[[bestRedRegime_index]], regime=whichRedToUse)
          }

          # If redundant regimes are left: which ones should be random
          if(rbinom(1, 1, 0.3)==1) {
            tmp = whichRedLight
          } else {
            tmp = whichRedundant
          }
          whichRandom = tmp[which(tmp != whichRedToUse)]

          # Then smart mutations with the combined individual
          H2[,mutate] = vapply(1:length(mutate), function(x) smartIndividual_int(p=p, M=M_orig, params=combinedInd, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted,
                                                                             constraints=constraints, R=R, ar0scale=ar0scale,
                                                                             sigmascale=sigmascale, accuracy=accuracy[x],
                                                                             whichRandom=whichRandom), numeric(d))
        }
      } else { # Random mutations
        H2[,mutate] = vapply(1:length(mutate), function(x) randomIndividual_int(p, M_orig, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints,
                                                                                R=R, ar0scale=ar0scale, sigmascale=sigmascale), numeric(d))
      }
    }

    # Sort components if constraints are not used
    if(constraints==FALSE) {
      H2 = vapply(1:popsize, function(i2) sortComponents(p, M_orig, H2[,i2], StMAR=StMAR, GStMAR=GStMAR, restricted=restricted), numeric(d))
    }

    # Set up for the next generation
    G = H2
    survivors_old = survivors
  }

  # Return the best individual
  bestIndex = which(logliks == max(logliks), arr.ind = TRUE)[1,]
  bestind = generations[, bestIndex[2], bestIndex[1]]
  return(bestind)
}


#' @title Create random GMAR, StMAR or G-StMAR model compatible parameter vector
#'
#' @description FOR INTERNAL USE. \code{randomIndividual_int} creates a random GMAR, StMAR or G-StMAR model compatible parameter vector.
#'
#' \code{smartIndividual_int} creates a random GMAR, StMAR or G-StMAR model compatible parameter vector close to argument \code{params}.
#'
#' @inheritParams GAfit
#' @param ar0scale a real valued vector of length two specifying the mean (the first element) and standard deviation (the second element) of the normal distribution
#'  from which the \eqn{\phi_{m,0}} parameters (for random regimes) should be generated.
#' @param sigmascale a positive real number specifying the standard deviation of the (zero mean, positive only) normal distribution
#'  from which the component variance parameters (for random regimes) should be generated.
#' @inherit GAfit return
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. (2015) Gaussian Mixture Autoregressive model for univariate time series.
#'          \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item References regarding the StMAR and G-StMAR models will be updated after they are published.
#'  }

randomIndividual_int <- function(p, M, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, ar0scale, sigmascale) {
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  # The cases of constrained models
  if(constraints==TRUE & restricted==FALSE) {
    for(j in 1:50) { # Restrict the number of attempts to find stationary individual
      ind = c()
      for(i1 in 1:M) { # Run through components
        R0 = as.matrix(R[[i1]])
        q = ncol(R0)
        scale = sum(abs(R0))
        ind = c(ind, rnorm(n=1, mean=ar0scale[1], sd=ar0scale[2]), rnorm(n=q, mean=0, sd=0.6/scale), abs(rnorm(n=1, mean=0, sd=sigmascale)))
      }
      if(M>1 & GStMAR==FALSE) { # alphas
        alphas = runif(n=M)
        alphas = alphas[order(alphas, decreasing=TRUE)]/sum(alphas)
        alphas = alphas[-M]
        ind = c(ind, alphas)
      } else if(GStMAR==TRUE) {
        alphas = runif(n=M)
        alphas = alphas/sum(alphas)
        alphasM1 = alphas[1:M1]
        alphasM2 = alphas[(M1+1):M]
        alphas = c(alphasM1[order(alphasM1, decreasing=TRUE)], alphasM2[order(alphasM2, decreasing=TRUE)])
        alphas = alphas[-M]
        ind = c(ind, alphas)
      }
      if(StMAR==TRUE) {
        ind = c(ind, 2+rgamma(M, shape=0.7, rate=0.008)) #runif(n=M, min=2, max=350)
      } else if(GStMAR==TRUE) {
        ind = c(ind, 2+rgamma(M2, shape=0.7, rate=0.008))
      }
      ind_stand = reformConstrainedPars(p=p, M=M_orig, params=ind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, R=R)
      if(isStationary_int(p=p, M=M, params=ind_stand, restricted=restricted)) {
        return(ind)
      }
    }
    stop("Failed to create random stationary parameter vector. Please make sure that all constraint matrices are sensible!")
  } else if(constraints==TRUE & restricted==TRUE) {
    for(j in 1:50) { # Restrict the number of attempts to find stationary individual
      R0 = as.matrix(R)
      q = ncol(R0)
      scale = sum(abs(R0))
      ind = c(rnorm(n=M, mean=ar0scale[1], sd=ar0scale[2]), rnorm(n=q, mean=0, sd=0.6/scale), abs(rnorm(n=M, mean=0, sd=sigmascale)))
      if(M>1 & GStMAR==FALSE) { # alphas
        alphas = runif(n=M)
        alphas = alphas[order(alphas, decreasing=TRUE)]/sum(alphas)
        alphas = alphas[-M]
        ind = c(ind, alphas)
      } else if(GStMAR==TRUE) {
        alphas = runif(n=M)
        alphas = alphas/sum(alphas)
        alphasM1 = alphas[1:M1]
        alphasM2 = alphas[(M1+1):M]
        alphas = c(alphasM1[order(alphasM1, decreasing=TRUE)], alphasM2[order(alphasM2, decreasing=TRUE)])
        alphas = alphas[-M]
        ind = c(ind, alphas)
      }
      if(StMAR==TRUE) {
        ind = c(ind, 2+rgamma(M, shape=0.7, rate=0.008)) #runif(n=M, min=2, max=350)
      } else if(GStMAR==TRUE) {
        ind = c(ind, 2+rgamma(M2, shape=0.7, rate=0.008))
      }
      ind_stand = reformConstrainedPars(p=p, M=M_orig, params=ind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, R=R)
      if(isStationary_int(p=p, M=M, params=ind_stand, restricted=restricted)) {
        return(ind)
      }
    }
    stop("Failed to create random stationary parameter vector. Please make sure that the constraint matrix is sensible!")
  }

  # The cases of unconstrained models
  if(restricted==FALSE | M==1) {
    while(TRUE) {
      ind = as.vector(vapply(1:M, function(x) c(rnorm(n=1, mean=ar0scale[1], sd=ar0scale[2]), rnorm(n=p, mean=0, sd=0.6/p),
                                                abs(rnorm(n=1, mean=0, sd=sigmascale))), numeric(p+2)))
      if(M>1 & GStMAR==FALSE) {
        alphas = runif(n=M)
        alphas = alphas[order(alphas, decreasing=TRUE)]/sum(alphas)
        alphas = alphas[-M]
        ind = c(ind, alphas)
      } else if(GStMAR==TRUE) {
        alphas = runif(n=M)
        alphas = alphas/sum(alphas)
        alphasM1 = alphas[1:M1]
        alphasM2 = alphas[(M1+1):M]
        alphas = c(alphasM1[order(alphasM1, decreasing=TRUE)], alphasM2[order(alphasM2, decreasing=TRUE)])
        alphas = alphas[-M]
        ind = c(ind, alphas)
      }
      min(2+rgamma(1000, shape=0.9, rate=0.008))
      if(isStationary_int(p=p, M=M, params=ind, restricted=restricted)) {
        if(StMAR==TRUE) {
          ind = c(ind, 2+rgamma(M, shape=0.7, rate=0.008)) #runif(n=M, min=2, max=350)
        } else if(GStMAR==TRUE) {
          ind = c(ind, 2+rgamma(M2, shape=0.7, rate=0.008))
        }
        break
      }
    }
  } else { # If restricted
    while(TRUE) {
      ind = c(rnorm(n=M, mean=ar0scale[1], sd=ar0scale[2]),  rnorm(p, mean=0, sd=0.6/p), abs(rnorm(M, mean=0, sd=sigmascale)))
      alphas = runif(n=M)
      if(GStMAR==TRUE) {
        alphas = alphas/sum(alphas)
        alphasM1 = alphas[1:M1]
        alphasM2 = alphas[(M1+1):M]
        alphas = c(alphasM1[order(alphasM1, decreasing=TRUE)], alphasM2[order(alphasM2, decreasing=TRUE)])
      } else {
        alphas = alphas[order(alphas, decreasing=TRUE)]/sum(alphas)
      }
      alphas = alphas[-M]
      ind = c(ind, alphas)
      if(isStationary_int(p=p, M=M, params=ind, restricted=restricted)) {
        if(StMAR==TRUE) {
          ind = c(ind, 2+rgamma(M, shape=0.7, rate=0.008)) #runif(n=M, min=2, max=350)
        } else if(GStMAR==TRUE) {
          ind = c(ind, 2+rgamma(M2, shape=0.7, rate=0.008))
        }
        break
      }
    }
  }
  return(ind)
}


#' @title Create somewhat random GMAR, StMAR or G-StMAR model compatible parameter vector
#'
#' @description \code{randomIndividual} creates a somewhat random GMAR, StMAR, G-StMAR model compatible parameter vector.
#'
#' \code{smartIndividual} creates a somewhat random GMAR, StMAR or G-StMAR model compatible parameter vector close to argument \code{params}.
#'   Sometimes returns exactly the given parameter vector.
#'
#' @inheritParams randomIndividual_int
#' @inherit randomIndividual_int return references
#' @examples
#' # GMAR model parameter vector
#' params22 <- randomIndividual(2, 2, ar0scale=c(0, 1), sigmascale=1)
#' smart22 <- smartIndividual(2, 2, params22, accuracy=10)
#' cbind(params22, smart22)
#'
#'
#' # Restricted GMAR parameter vector
#' params12r <- randomIndividual(1, 2, restricted=TRUE, ar0scale=c(-2, 2), sigmascale=2)
#' smart12r <- smartIndividual(1, 2, params12r, restricted=TRUE, accuracy=20)
#' cbind(params12r, smart12r)
#'
#'
#' # StMAR parameter vector: first regime is random in the "smart individual"
#' params13t <- randomIndividual(1, 3, StMAR=TRUE, ar0scale=c(3, 1), sigmascale=3)
#' smart13t <- smartIndividual(1, 3, params13t, StMAR=TRUE, accuracy=15,
#'                             ar0scale=c(3, 3), sigmascale=3, whichRandom=1)
#' cbind(params13t, smart13t)
#'
#'
#' # Restricted StMAR parameter vector
#' params22tr <- randomIndividual(2, 2, StMAR=TRUE, restricted=TRUE,
#'                                ar0scale=c(3, 2), sigmascale=0.5)
#' smart22tr <- smartIndividual(2, 2, params22tr, StMAR=TRUE, restricted=TRUE,
#'                              accuracy=30)
#' cbind(params22tr, smart22tr)
#'
#'
#' # G-StMAR parameter vector
#' params12gs <- randomIndividual(1, c(1, 1), GStMAR=TRUE, ar0scale=c(0, 1), sigmascale=1)
#' smart12gs <- smartIndividual(1, c(1, 1), params12gs, GStMAR=TRUE, accuracy=20)
#' cbind(params12gs, smart12gs)
#'
#'
#' # Restricted G-StMAR parameter vector
#' params23gsr <- randomIndividual(2, c(1, 2), GStMAR=TRUE, restricted=TRUE,
#'                                 ar0scale=c(-1, 1), sigmascale=3)
#' smart23gsr <- smartIndividual(2, c(1, 2), params23gsr, GStMAR=TRUE, restricted=TRUE,
#'                               ar0scale=c(0, 1), sigmascale=1, accuracy=20, whichRandom=2)
#' cbind(params23gsr, smart23gsr)
#'
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' R <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- randomIndividual(2, 2, constraints=TRUE, R=R,
#'                               ar0scale=c(1, 1), sigmascale=1)
#' smart22c <- smartIndividual(2, 2, params22c, constraints=TRUE, R=R, accuracy=10)
#' cbind(params22c, smart22c)
#'
#'
#' # Such constrained StMAR(3, 2) model that the second order AR coefficients
#' # are constrained to zero.
#' R0 = matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
#' R = list(R0, R0)
#' params32c <- randomIndividual(3, 2, StMAR=TRUE, constraints=TRUE, R=R,
#'                               ar0scale=c(1, 1), sigmascale=1)
#' smart32c <- smartIndividual(3, 2, params32c, StMAR=TRUE, constraints=TRUE, R=R, accuracy=10)
#' cbind(params32c, smart32c)
#'
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero. Second regime is random in the "smart individual".
#' params32trc <- randomIndividual(3, 2, StMAR=TRUE, restricted=TRUE, constraints=TRUE,
#'                                 R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2), ar0scale=c(-2, 0.5),
#'                                 sigmascale=4)
#' smart32trc <- smartIndividual(3, 2, params32trc, StMAR=TRUE, restricted=TRUE,
#'                               constraints=TRUE, R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2),
#'                               ar0scale=c(0, 0.1), sigmascale=0.1, whichRandom=2,
#'                               accuracy=20)
#' cbind(params32trc, smart32trc)
#' @export

randomIndividual <- function(p, M, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, ar0scale, sigmascale) {
  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p=p, M=M, GStMAR=GStMAR)
  if(length(ar0scale)!=2) {
    stop("ar0scale is wrong dimension")
  } else if(length(sigmascale)!=1) {
    stop("sigmascale is wrong dimension")
  } else if(ar0scale[2]<=0) {
    stop("The second element of ar0scale should be larger than zero")
  } else if(sigmascale<=0) {
    stop("sigmascale should be larger than zero")
  }
  if(constraints==TRUE) {
    checkConstraintMat(p=p, M=M, R=R, restricted=restricted)
  }
  return(randomIndividual_int(p=p, M=M, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R, ar0scale=ar0scale, sigmascale=sigmascale))
}


#' @rdname randomIndividual_int
#' @inheritParams loglikelihood
#' @param accuracy a real number larger than zero specifying how close to \code{params} the generated parameter vector should be.
#'   Standard deviation of the normal distribution from which new parameter values are drawed from will be corresponding parameter value divided by \code{accuracy}.
#' @param whichRandom an (optional) numeric vector of max length \code{M} specifying which regimes should be random instead of "smart" when
#' using \code{smartIndividual}. Does not affect on mixing weight parameters. Default in none.

smartIndividual_int <- function(p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, ar0scale, sigmascale, accuracy, whichRandom) {

  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  if(missing(whichRandom)) {
    whichRandom = numeric(0)
  }
  # The cases of constrained model
  if(constraints==TRUE & restricted==FALSE) {
    ind = c()
    j = 0
    for(i1 in 1:M) { # Run through components
      R0 = as.matrix(R[[i1]])
      q = ncol(R0)
      if(i1 %in% whichRandom) { # If regime should be random
        scale = sum(abs(R0))
        ind = c(ind, rnorm(n=1, mean=ar0scale[1], sd=ar0scale[2]), rnorm(n=q, mean=0, sd=0.6/scale), abs(rnorm(n=1, mean=0, sd=sigmascale)))
      } else {
        params_m = params[(j+1):(j+q+1)] # no sigma, alphas or dfs
        sigma_m = params[j+q+2]
        ind = c(ind, rnorm(n=q+1, mean=params_m, sd=abs(params_m/accuracy)), abs(rnorm(n=1, mean=sigma_m, sd=sigma_m/accuracy)))
      }
      j = j + q + 2
    }
    if(M>1) {
      alphas = params[(j+1):(j+M-1)]
      alphas2 = abs(rnorm(n=M-1, mean=alphas, sd=alphas/accuracy))
      ind = c(ind, alphas2)
    } else {
      alphas2 = 0
    }
    if(StMAR==TRUE | GStMAR==TRUE) {
 #     dfs = params[(j+M):(j+2*M-1)]
      d = length(params)
      if(StMAR==TRUE) {
        M2 = M # So we can just use M2 and M1 for StMAR models too
        M1 = 0
      }
      dfs = params[(d-M2+1):d] # degrees of freedom
      if(length(whichRandom)==0) {
        dfs2 = rnorm(n=M2, mean=dfs, sd=dfs/accuracy)
      } else { # If some regimes should be random
        dfs2 = numeric(0)
        for(i1 in (M1+1):M) { # Generate dfs
          if(i1 %in% whichRandom) {
            dfs2 = c(dfs2, 2+rgamma(1, shape=0.7, rate=0.008))# runif(n=1, min=2, max=350))
          } else {
            dfs2 = c(dfs2, rnorm(n=1, mean=dfs[i1-M1], sd=abs(dfs[i1-M1]/accuracy)))
          }
        }
      }
      ind = c(ind, dfs2)
      ind_stand = reformConstrainedPars(p=p, M=M_orig, params=ind, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind_stand, restricted=restricted) & all(dfs2>2)) {
        return(ind)
      } else {
        return(params) # Return the best so-far individual if smart mutation is not stationary
      }
    } else { # If StMAR==FALSE and GStMAR==FALSE
      ind_stand = reformConstrainedPars(p=p, M=M_orig, params=ind, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind_stand, restricted=restricted)) {
        return(ind)
      } else {
        return(params) # Return the best so-far individual if smart mutation is not stationary
      }
    }
  } else if(constraints==TRUE & restricted==TRUE) {
    q = ncol(as.matrix(R))
    params0 = params[1:(M+q)]
    sigmas = params[(M+q+1):(2*M+q)]
    if(length(whichRandom)==0) {
      ind = c(rnorm(n=M+q, mean=params0, sd=abs(params0/accuracy)), abs(rnorm(n=M, mean=sigmas, sd=sigmas/accuracy)))
    } else { # If some regime(s) should be random
      scale = sum(abs(R))
      ind = numeric(0)
      for(i1 in 1:M) { # Generate phi0
        if(i1 %in% whichRandom) {
          ind = c(ind, rnorm(n=1, mean=ar0scale[1], sd=ar0scale[2]))
        } else {
          ind = c(ind, rnorm(n=1, mean=params0[i1], sd=abs(params0[i1]/accuracy)))
        }
      }
      ind = c(ind, rnorm(q, mean=params0[(M+1):(M+q)], sd=abs(params0[(M+1):(M+q)]/accuracy)))
      for(i1 in 1:M) { # Generate sigmas
        if(i1 %in% whichRandom) {
          ind = c(ind, abs(rnorm(n=1, mean=0, sd=sigmascale)))
        } else {
          ind = c(ind, rnorm(n=1, mean=sigmas[i1], sd=abs(sigmas[i1]/accuracy)))
        }
      }
    }
    if(M>1) {
      alphas = params[(q+2*M+1):(3*M+q-1)]
      alphas2 = rnorm(n=M-1, mean=alphas, sd=alphas/accuracy)
      ind = c(ind, alphas2)
    }
    if(StMAR==TRUE | GStMAR==TRUE) {
      d = length(params)
      if(StMAR==TRUE) {
        M2 = M # So we can just use M2 and M1 for StMAR models too
        M1 = 0
      }
      dfs = params[(d-M2+1):d] # degrees of freedom
      if(length(whichRandom)==0) {
        dfs2 = rnorm(n=M2, mean=dfs, sd=dfs/accuracy)
      } else { # If some regimes should be random
        dfs2 = numeric(0)
        for(i1 in (M1+1):M) { # Generate dfs
          if(i1 %in% whichRandom) {
            dfs2 = c(dfs2, 2+rgamma(1, shape=0.7, rate=0.008)) #runif(n=1, min=2, max=350))
          } else {
            dfs2 = c(dfs2, rnorm(n=1, mean=dfs[i1-M1], sd=abs(dfs[i1-M1]/accuracy)))
          }
        }
      }
      ind = c(ind, dfs2)
      ind_stand = reformConstrainedPars(p=p, M=M_orig, params=ind, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind_stand, restricted=restricted) & all(dfs2>2)) {
        return(ind)
      } else {
        return(params) # Return the best so-far individual if smart mutation is not proper
      }
    } else { # If StMAR==FALSE
      ind_stand = reformConstrainedPars(p=p, M=M_orig, params=ind, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind_stand, restricted=restricted)) {
        return(ind)
      } else {
        return(params) # Return the best so-far individual if smart mutation is not proper
      }
    }
  }

  # The cases of unconstrained models
  if(restricted==FALSE | M==1) {
    pars = matrix(params[1:(M*(p+2))], ncol=M) # Component parameters by column (except alphas and dfs)
    if(length(whichRandom)==0) {
      pars2 = vapply(1:M, function(i1) rnorm(n=p+2, mean=pars[,i1], sd=abs(pars[,i1]/accuracy)), numeric(p+2))
    } else { # If some regime(s) should be random
      pars2 = numeric(0)
      for(i1 in 1:M) {
        if(i1 %in% whichRandom) {
          pars2 = cbind(pars2, c(rnorm(n=1, mean=ar0scale[1], sd=ar0scale[2]), rnorm(n=p, mean=0, sd=0.6/p),
                                 abs(rnorm(n=1, mean=0, sd=sigmascale))))
        } else {
          pars2 = cbind(pars2, rnorm(n=p+2, mean=pars[,i1], sd=abs(pars[,i1]/accuracy)))
        }
      }
    }
    pars2[p+2,] = abs(pars2[p+2,]) # Make sure variances are larger than zero
    ind = as.vector(pars2)
    if(M>1) {
      alphas = params[(M*(p+2)+1):(M*(p+3)-1)]
      alphas2 = abs(rnorm(n=M-1, mean=alphas, sd=alphas/accuracy))
      ind = c(ind, alphas2)
    } else {
      alphas2 = 0
    }
    if(StMAR==TRUE | GStMAR==TRUE) {
      d = length(params)
      if(StMAR==TRUE) {
        M2 = M # So we can just use M2 and M1 for StMAR models too
        M1 = 0
      }
      dfs = params[(d-M2+1):d] # degrees of freedom
      if(length(whichRandom)==0) {
        dfs2 = rnorm(n=M2, mean=dfs, sd=dfs/accuracy)
      } else { # If some regimes should be random
        dfs2 = numeric(0)
        for(i1 in (M1+1):M) { # Generate dfs
          if(i1 %in% whichRandom) {
            dfs2 = c(dfs2, 2+rgamma(1, shape=0.7, rate=0.008)) #runif(n=1, min=2, max=350))
          } else {
            dfs2 = c(dfs2, rnorm(n=1, mean=dfs[i1-M1], sd=abs(dfs[i1-M1]/accuracy)))
          }
        }
      }
      ind = c(ind, dfs2)
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind, restricted=restricted) & all(dfs2>2)) {
        ind = sortComponents(p=p, M=M_orig, params=ind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
        return(ind)
      } else {
        return(params) # Return params if smart mutation is not proper
      }
    } else { # If StMAR==FALSE and GStMAR==FALSE
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind, restricted=restricted)) {
        ind = sortComponents(p=p, M=M_orig, params=ind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
        return(ind)
      } else {
        return(params) # Return params if smart mutation is not proper
      }
    }
  } else { # If restricted (M=1 is already considered)
    phi0 = params[1:M]
    arcoefs = params[(M+1):(M+p)]
    sigmas = params[(M+p+1):(p+2*M)]
    alphas = params[(p+2*M+1):(3*M+p-1)]
    if(length(whichRandom)==0) {
      ind = c(rnorm(n=M, mean=phi0, sd=abs(phi0/accuracy)), rnorm(n=p, mean=arcoefs, sd=abs(arcoefs/accuracy)),
              abs(rnorm(n=M, mean=sigmas, sd=abs(sigmas/accuracy))))
    } else { # If some regime(s) should be random
      ind = numeric(0)
      for(i1 in 1:M) { # Generate phi0 parameters
        if(i1 %in% whichRandom) {
          ind = c(ind, rnorm(n=1, mean=ar0scale[1], sd=ar0scale[2]))
        } else {
          ind = c(ind, rnorm(n=1, mean=phi0[i1], sd=abs(phi0[i1]/accuracy)))
        }
      }
      ind = c(ind, rnorm(n=p, mean=arcoefs, sd=abs(arcoefs/accuracy)))
      for(i1 in 1:M) { # Generate variance parameters
        if(i1 %in% whichRandom) {
          ind = c(ind, abs(rnorm(n=1, mean=0, sd=sigmascale)))
        } else {
          ind = c(ind, rnorm(n=1, mean=sigmas[i1], sd=abs(sigmas[i1]/accuracy)))
        }
      }
    }
    alphas2 = abs(rnorm(n=M-1, mean=alphas, sd=alphas/accuracy))
    ind = c(ind, alphas2)
    if(StMAR==TRUE | GStMAR==TRUE) {
      d = length(params)
      if(StMAR==TRUE) {
        M2 = M # So we can just use M2 and M1 for StMAR models too
        M1 = 0
      }
      dfs = params[(d-M2+1):d] # degrees of freedom
      if(length(whichRandom)==0) {
        dfs2 = rnorm(n=M2, mean=dfs, sd=dfs/accuracy)
      } else { # If some regimes should be random
        dfs2 = numeric(0)
        for(i1 in (M1+1):M) { # Generate dfs
          if(i1 %in% whichRandom) {
            dfs2 = c(dfs2, 2+rgamma(1, shape=0.7, rate=0.008))  #runif(n=1, min=2, max=350))
          } else {
            dfs2 = c(dfs2, rnorm(n=1, mean=dfs[i1-M1], sd=abs(dfs[i1-M1]/accuracy)))
          }
        }
      }
      ind = c(ind, dfs2)
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind, restricted=restricted) & all(dfs2>2)) {
        ind = sortComponents(p=p, M=M_orig, params=ind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
        return(ind)
      } else {
        return(params) # Return the best so-far individual if smart mutation is not stationary
      }
    } else {
      if(sum(alphas2)<1 & isStationary_int(p=p, M=M, params=ind, restricted=restricted)) {
        ind = sortComponents(p=p, M=M_orig, params=ind, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
        return(ind)
      } else {
        return(params) # Return the best so-far individual if smart mutation is not stationary
      }
    }
  }
}


#' @rdname randomIndividual
#' @inheritParams smartIndividual_int
#' @export

smartIndividual <- function(p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, ar0scale, sigmascale, accuracy, whichRandom) {
  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p=p, M=M, GStMAR=GStMAR)
  if(constraints==TRUE) {
    checkConstraintMat(p=p, M=M, R=R, restricted=restricted)
  }
  if(accuracy<=0) {
    stop("Argument accuracy has to be larger than zero")
  } else if(length(params)!=nParams(p=p, M=M, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)) {
    stop("The parameter vector is wrong dimension")
  } else if(!missing(whichRandom)) {
    if(length(whichRandom)>sum(M)) {
      stop("Length of whichRandom should not be larger than M")
    } else if(any(!whichRandom %in% 1:sum(M))) {
      stop("The elements of whichRandom should be positive integers in the clsoed interval [1, M]")
    } else if(length(whichRandom)>0) {
      if(length(ar0scale)!=2) {
        stop("ar0scale is wrong dimension")
      } else if(length(sigmascale)!=1) {
        stop("sigmascale is wrong dimension")
      } else if(ar0scale[2]<=0) {
        stop("The second element of ar0scale should be larger than zero")
      } else if(sigmascale<=0) {
        stop("sigmascale should be larger than zero")
      }
    }
  }
  return(smartIndividual_int(p=p, M=M, params=params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R, ar0scale=ar0scale, sigmascale=sigmascale, accuracy=accuracy, whichRandom=whichRandom))
}


#' @title Extract regime from a parameter vector
#'
#' @description \code{extractRegime} extracts the specified regime from the GMAR, StMAR or G-StMAR model's parameter vector. Doesn't extract mixing weight parameter alpha.
#' @inheritParams loglikelihood
#' @param regime a positive integer in the interval [1, M] defining which regime should be extracted.
#' @return Returns a numeric vector corresponding to the regime with...
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(p+2x1)} vector \eqn{(\phi_{m,0},\phi_{m,1},...,\phi_{m,p}, \sigma_{m}^2)}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(p+3x1)} vector \eqn{(\phi_{m,0},\phi_{m,1},...,\phi_{m,p}, \sigma_{m}^2, \nu_{m})}.}
#'        \item{For \strong{G-StMAR} model:}{Same as GMAR for GMAR-components and same as StMAR for StMAR-components.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above, but vector \strong{\eqn{\phi_{m}}} replaced with
#'         vector \strong{\eqn{\psi_{m}}} that satisfies \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(2x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2)}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(3x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2, \nu_{m})}.}
#'        \item{For \strong{G-StMAR} model:}{Same as GMAR for GMAR-components and same as StMAR for StMAR-components.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above.}
#'      }
#'    }
#'  }

extractRegime <- function(p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, regime) {
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  if(restricted==FALSE) {
    j = 0 # Indicates where we at
    for(i1 in 1:M) { # Go through regimes
      if(constraints==TRUE) {
        R0 = as.matrix(R[[i1]])
        q = ncol(R0)
      } else {
        q = p
      }
      if(i1 == regime) { # Take parameters
        params0 = params[(j+1):(j+q+2)]
      }
      j = j+q+2
    }
    if(StMAR==TRUE) {
      params0 = c(params0, params[j+M-1+regime]) # dfs
    } else if(GStMAR==TRUE & regime>M1) {
      params0 = c(params0, params[j+M-1+regime-M1]) # dfs
    }
    return(params0)
  } else { # If restricted==TRUE
    if(constraints==TRUE) {
      q = ncol(as.matrix(R))
    } else {
      q = p
    }
    params0 = c(params[regime], params[M+q+regime])
    if(StMAR==TRUE) {
      params0 = c(params0, params[3*M+q-1+regime])
    } else if(GStMAR==TRUE & regime>M1) {
      params0 = c(params0, params[3*M+q-1+regime-M1])
    }
    return(params0)
  }
}


#' @title Change the specified regime of parameter vector to the given regime-parameter vector
#'
#' @description \code{changeRegime} changes the specified regime of the parameter vector to correspond the given
#'  regime-parameter vector and returns the modified parameter vector. Does not affect mixing weight parameters.
#'
#' @inheritParams loglikelihood
#' @param regimeParams a numeric vector specifying the parameter values that should be inserted to the specified regime.
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(p+2x1)} vector \eqn{(\phi_{m,0},\phi_{m,1},...,\phi_{m,p}, \sigma_{m}^2)}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(p+3x1)} vector \eqn{(\phi_{m,0},\phi_{m,1},...,\phi_{m,p}, \sigma_{m}^2, \nu_{m})}.}
#'        \item{For \strong{G-StMAR} model:}{Same as GMAR for GMAR-components and same as StMAR for StMAR-components.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above, but vector \strong{\eqn{\phi_{m}}} replaced with
#'         vector \strong{\eqn{\psi_{m}}} that satisfies \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(2x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2)}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(3x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2, \nu_{m})}.}
#'        \item{For \strong{G-StMAR} model:}{Same as GMAR for GMAR-components and same as StMAR for StMAR-components.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above.}
#'      }
#'    }
#'  }
#' @param regime a positive integer in the interval [1, M] defining which regime should be changed.
#' @return Returns modified parameter vector of the form descripted in \code{params}.

changeRegime <- function(p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, regimeParams, regime) {
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  } else {
    M1 = 0 # Just to allow cleaner code below
  }
  if(restricted==FALSE) {
    params0 = numeric(0)
    j = 0 # Indicates where we at
    for(i1 in 1:M) { # Go through regimes
      if(constraints==TRUE) {
        R0 = as.matrix(R[[i1]])
        q = ncol(R0)
      } else {
        q = p
      }
      if(i1 == regime) { # Change the parameters
        if(StMAR==TRUE | (GStMAR==TRUE & regime>M1)) {
          regimeDfs = regimeParams[length(regimeParams)]
          regimeParams = regimeParams[-length(regimeParams)] # Delete dfs
        }
        params0 = c(params0, regimeParams)
      } else { # Use same parameters
        params0 = c(params0, params[(j+1):(j+q+2)])
      }
      j = j+q+2
    }
    if(M>1) {
      params0 = c(params0, params[(j+1):(j+M-1)]) # Add alphas
    }
    if(StMAR==TRUE | GStMAR==TRUE) {
      if(StMAR==TRUE) {
        M1 = 0 # So that we can use in the both cases
      }
      dfs = params[(j+M):(j+2*M-1-M1)]
      if(StMAR==TRUE | (GStMAR==TRUE & regime>M1)) {
        dfs[regime-M1] = regimeDfs # Change the dfs
      }
      params0 = c(params0, dfs) # Add dfs
    }
    return(params0)
  } else { # Restricted == TRUE
    if(constraints==TRUE) {
      q = ncol(as.matrix(R))
    } else {
      q = p
    }
    params0 = params
    params0[regime] = regimeParams[1] # phi0
    params0[M+q+regime] = regimeParams[2] # sigma^2
    if(StMAR==TRUE) {
      params0[3*M+q-1+regime] = regimeParams[3] # dfs
    } else if(GStMAR==TRUE & regime>M1) {
      params0[3*M+q-1+regime-M1] = regimeParams[3] # dfs
    }
    return(params0)
  }
}
