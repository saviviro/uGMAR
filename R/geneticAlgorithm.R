#' @import stats
#'
#' @title Genetic algorithm for preliminary estimation of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{GAfit} estimates specified GMAR, StMAR, or G-StMAR model using a genetic algorithm.
#'  The employed genetic algorithm is designed to find starting values for gradient based methods.
#'
#' @inheritParams loglikelihood_int
#' @param ngen a positive integer specifying the number of generations to be ran through in the genetic algorithm.
#' @param popsize a positive even integer specifying the population size in the genetic algorithm.
#'  Default is \code{10*d} where \code{d} is the number of parameters.
#' @param smartMu a positive integer specifying the generation after which the random mutations in the genetic algorithm are "smart".
#'  This means that mutating individuals will mostly mutate fairly close (or partially close) to the best fitting individual so far.
#'  Default is \code{min(100, round(0.5*ngen))}.
#' @param meanscale a real valued vector of length two specifying the mean (the first element) and standard deviation (the second element)
#'  of the normal distribution from which the \eqn{\mu_{m}} mean-parameters are generated in random mutations in the genetic algorithm.
#'  Default is \code{c(mean(data), sd(data))}.
#'  Note that the genetic algorithm optimizes with mean-parametrization even when \code{parametrization=="intercept"}, but
#'  input (in \code{initpop}) and output (return value) parameter vectors may be intercept-parametrized.
#' @param sigmascale a positive real number specifying the standard deviation of the (zero mean, positive only by taking absolute value)
#'  normal distribution from which the component variance parameters are generated in the random mutations in the genetic algorithm.
#'  Default is \code{var(stats::ar(data, order.max=10)$resid, na.rm=TRUE)}.
#' @param initpop a list of parameter vectors from which the initial population of the genetic algorithm will be generated from.
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
#'          matrices \strong{C} that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{C_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
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
#'          \strong{\eqn{C}} that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{C\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight, and \eqn{v} a degrees of
#'  freedom parameter.
#'  Note that in the case \strong{M=1}, the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model,
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#'  If not specified (or \code{FALSE} as is default), the initial population will be drawn randomly.
#' @param regime_force_scale a non-negative real number specifying how much should natural selection favour individuals
#'   with less regimes that have almost all mixing weights (practically) at zero (see \code{red_criteria}), i.e., with
#'   less "redundant regimes".
#'   Set to zero for no favouring or large number for heavy favouring. Without any favouring the genetic algorithm gets more often stuck
#'   in an area of the parameter space where some regimes are wasted, but with too much favouring the best genes might never mix into the
#'   population and the algorithm might converge poorly. Default is \code{1} and it gives \eqn{2x} larger surviving probability weights for
#'   individuals with no wasted regimes compared to individuals with one wasted regime. Number \code{2} would give \eqn{3x} larger probabilities etc.
#' @param red_criteria a length 2 numeric vector specifying the criteria that is used to determine whether a regime is redundant or not.
#'   Any regime \code{m} which satisfies \code{sum(mixingWeights[,m] > red_criteria[1]) < red_criteria[2]*n_obs} will be considered "redundant".
#'   One should be careful when adjusting this argument (set \code{c(0, 0)} to fully disable the 'redundant regime' features from the algorithm).
#' @param to_return should the genetic algorithm return the best fitting individual which has the least "redundant" regimes (\code{"alt_ind"})
#'   or the individual which has the highest log-likelihood in general (\code{"best_ind"}) but might have more wasted regimes?
#' @param minval a real number defining the minimum value of the log-likelihood function that will be considered.
#'   Values smaller than this will be treated as they were \code{minval} and the corresponding individuals will never survive.
#'   The default is \code{-(10^(ceiling(log10(length(data))) + 1) - 1)}, and one should be very careful if adjusting this.
#' @param seed a single value, interpreted as an integer, or NULL, that sets seed for the random number generator in the beginning of
#'   the function call. If calling \code{GAfit} from \code{fitGSMAR}, use the argument \code{seeds} instead of passing the argument \code{seed}.
#' @details
#'    The core of the genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)}.
#'    It utilizes a slightly modified version of the individually adaptive crossover and mutation rates described
#'    by \emph{Patnaik and Srinivas (1994)} and employs (50\%) fitness inheritance discussed by \emph{Smith, Dike and Stegmann (1995)}.
#'    Large (in absolute value) but stationary AR parameter values are generated with the algorithm proposed by Monahan (1984).
#'
#'    By "redundant" or "wasted" regimes we mean regimes that have the time varying mixing weights basically at zero for all t.
#'    The model with redundant regimes would have approximately the same log-likelihood value without the redundant regimes
#'    and there is no purpose to have redundant regimes in the model.
#' @return Returns estimated parameter vector with the form described in \code{initpop}.
#' @references
#'  \itemize{
#'    \item Dorsey R. E. and Mayer W. J. 1995. Genetic algorithms for estimation problems with multiple optima,
#'          nondifferentiability, and other irregular features. \emph{Journal of Business & Economic Statistics},
#'          \strong{13}, 53-66.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'          \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'          arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item Monahan J.F. 1984. A Note on Enforcing Stationarity in Autoregressive-Moving Average Models.
#'          \emph{Biometrica} \strong{71}, 403-404.
#'    \item Patnaik L.M. and Srinivas M. 1994. Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Smith R.E., Dike B.A., Stegmann S.A. 1995. Fitness inheritance in genetic algorithms.
#'          \emph{Proceedings of the 1995 ACM Symposium on Applied Computing}, 345-350.
#'    \item There are currently no published references for G-StMAR model, but it's a straightforward generalization with
#'          theoretical properties similar to GMAR and StMAR models.
#'  }
#' @export

GAfit <- function(data, p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, parametrization=c("intercept", "mean"),
                  conditional=TRUE, ngen=200, popsize, smartMu, meanscale, sigmascale, initpop=NULL, regime_force_scale=1,
                  red_criteria=c(0.05, 0.01), to_return=c("alt_ind", "best_ind"), minval, seed=NULL) {
  set.seed(seed)
  model <- match.arg(model)
  check_model(model)
  checkPM(p=p, M=M, model=model)
  parametrization <- match.arg(parametrization)
  to_return <- match.arg(to_return)
  data <- checkAndCorrectData(data, p)
  n_obs <- length(data)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  }
  # Check constraint matrices
  if(!is.null(constraints)) {
    checkConstraintMat(p, M, restricted=restricted, constraints=constraints)
    if(restricted) {
      Cquals <- TRUE # States whether all matrices C are equal
    } else {
      Cquals <- length(unique(constraints)) == 1
    }
  } else {
    Cquals <- TRUE # Always true if no constraints
  }
  d <- nParams(p=p, M=M_orig, model=model, restricted=restricted, constraints=constraints)

  # Default settings
  if(missing(popsize)) popsize <- 10*d
  if(missing(smartMu)) smartMu <-  min(100, round(0.5*ngen))
  if(missing(meanscale)) meanscale <- c(mean(data), sd(data))
  if(missing(sigmascale)) sigmascale <- var(stats::ar(data, order.max=10)$resid, na.rm=TRUE)
  if(missing(minval)) minval <- get_minval(data)

  # Argument checks
  stopifnot(is.numeric(red_criteria) & length(red_criteria) == 2)
  if(!all_pos_ints(c(ngen, smartMu))) {
    stop("The arguments 'ngen' and 'smartMu' should be strictly positive integers")
  } else if(popsize%%2 != 0 | popsize < 2) {
    stop("The population size has to be EVEN positive integer")
  } else if(length(regime_force_scale) != 1 | regime_force_scale < 0) {
    stop("regime_force_scale should be non-negative real number")
  }

  # Initial population
  stopifnot(is.list(initpop) | is.null(initpop))
  if(is.null(initpop)) {
    nattempts <- 100
    for(i1 in 1:nattempts) {
      G <- replicate(popsize, randomIndividual_int(p, M_orig, model=model, restricted=restricted, constraints=constraints,
                                                   meanscale=meanscale, sigmascale=sigmascale), numeric(d))
      init_loks <- vapply(1:popsize, function(i2) loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], model=model, restricted=restricted,
                                                                    constraints=constraints, conditional=conditional, parametrization="mean",
                                                                    boundaries=TRUE, checks=FALSE, to_return="loglik", minval=minval), numeric(1))
      if(any(init_loks > minval)) break
      if(i1 == nattempts) stop("Failed to create initial population with good enough individuals. Consider setting up the initial population by hand using the argument 'initpop' of the function 'GAfit'.")
    }
  } else {
    # If the initial individuals are set by the user: check them and create the initial generation matrix G.
    n_inds <- length(initpop)
    for(i1 in 1:n_inds) {
      params <- initpop[[i1]]
      if(length(params) != d) {
        stop(paste("The parameter vector of individual", i1, "in the initial population is wrong dimension"))
      }
      if(parametrization == "intercept") {
        params <- change_parametrization(p=p, M=M_orig, params=params, model=model, restricted=restricted,
                                         constraints=constraints, change_to="mean")
      }
      if(!is.null(constraints)) {
        params <- reformConstrainedPars(p, M_orig, params, model=model, restricted=restricted, constraints=constraints)
      }
      pars_tmp <- reformParameters(p, M_orig, params=params, model=model, restricted=restricted)
      if(model == "StMAR" | model == "G-StMAR") {
        if(any(pars_tmp$dfs <= 2)) {
          stop(paste("The individual", i1, "in the initial population has degrees of freedom not larger than 2"))
        }
      }
      if(!isStationary_int(p, M, params=pars_tmp$params, restricted=restricted)) {
        stop(paste("The individual", i1, "in the initial population is not stationary"))
      } else if(any(pars_tmp$pars[p + 2,] <= 0)) {
        stop(paste("The individual", i1, "in the initial population has variance parameter that is not larger than zero"))
      } else if(sum(pars_tmp$alphas) > 1 + 1e-10) {
        stop(paste("The mixing weights of the individual", i1," in the initial population don't sum to one"))
      } else if(is.null(constraints)) {
        # Sort components in the initial population so that alpha_1>...>alpha_M
        initpop[[i1]] <- sortComponents(p=p, M=M_orig, params=params, model=model, restricted=restricted)
      }
    }
    G <- replicate(popsize, initpop[[sample.int(length(initpop), size=1)]])
    init_loks <- vapply(1:popsize, function(i2) loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], model=model, restricted=restricted,
                                                                  constraints=constraints, conditional=conditional, parametrization="mean",
                                                                  boundaries=TRUE, checks=FALSE, to_return="loglik", minval=minval), numeric(1))
    if(!any(init_loks > minval)) stop("The initial population does not contain good enough individuals. Check the log-likelihoods of the individuals with the function 'loglikelihood'.")
  }


  # Calculates the number of redundant regimes
  n_redundants <- function(M, mw) {
    sum(vapply(1:M, function(m) sum(mw[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
  }

  # Initial setup
  generations <- array(dim=c(d, popsize, ngen))
  logliks <- matrix(minval, nrow=ngen, ncol=popsize)
  redundants <- matrix(M, nrow=ngen, ncol=popsize) # Store the number of redundant regimes of each individual
  which_redundant_alt <- 1:M

  # Fills in the log-likelihood and number of redundant regimes
  fill_lok_and_red <- function(i1, i2, lok_and_mw) {
    if(!is.list(lok_and_mw)) {
      logliks[i1, i2] <<- minval
      redundants[i1, i2] <<- M
    } else {
      logliks[i1, i2] <<- lok_and_mw$loglik
      redundants[i1, i2] <<- n_redundants(M=M, mw=lok_and_mw$mw) # Number of redundant regimes
    }
  }

  ### Run through generations ###
  for(i1 in 1:ngen) {
    generations[, , i1] <- G # Save the generation's population to the container.

    ## Compute the log-likelihoods ##
    if(i1 > 1) {
      # Fitness inheritance: individual has 50% changes to inherit fitness if it's a result from crossover.
      I2 <- rep(I, each=2)
      which_did_co <- which(1 - which_not_co == 1) # Which individuals did crossover
      if(length(which_did_co) >= 1) {
        inherit <- sample(which_did_co, size=round(0.5*length(which_did_co)), replace=FALSE) # Draw the ones who inherit fitness
      } else {
        inherit <- numeric(0)
      }

      # survivor_liks holds the parents' log-likelihoods: for odd number they are (index, index+1) and for even (index-1, index)
      if(length(inherit) >= 1 & abs(maxLik - meanLik) > abs(0.03*meanLik)) { # If massive mutations: no fitness inheritance
        for(i2 in inherit) {
          if(i2 %% 2 == 0) {
            logliks[i1, i2] <- ((d - I2[i2])/d)*survivor_liks[i2 - 1] + (I2[i2]/d)*survivor_liks[i2]
            redundants[i1, i2] <- max(c(survivor_redundants[i2 - 1], survivor_redundants[i2]))
          } else {
            logliks[i1, i2] <- (I2[i2]/d)*survivor_liks[i2] + ((d - I2[i2])/d)*survivor_liks[i2 + 1]
            redundants[i1, i2] <- max(c(survivor_redundants[i2], survivor_redundants[i2 + 1]))
          }
        }
        noinherit <- (1:popsize)[-inherit]
      } else {
        noinherit <- 1:popsize
      }
      for(i2 in noinherit) { # Fill in the rest log-likelihoods
        if((mutate_indicator[i2] == 0 & which_not_co[i2] == 1) | all(H[,i2] == H2[,i2])) { # Individual is unchanged: no need to recalculate
          logliks[i1, i2] <- survivor_liks[i2]
          redundants[i1, i2] <- survivor_redundants[i2]
        } else {
          lok_and_mw <- loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], model=model, restricted=restricted,
                                          constraints=constraints, conditional=conditional, parametrization="mean",
                                          boundaries=TRUE, checks=FALSE, to_return="loglik_and_mw", minval=minval)
          fill_lok_and_red(i1, i2, lok_and_mw)
        }
      }
    } else { # The first round, i1 == 1
      for(i2 in 1:popsize) {
        lok_and_mw <- loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], model=model, restricted=restricted,
                                        constraints=constraints, conditional=conditional, parametrization="mean",
                                        boundaries=TRUE, checks=FALSE, to_return="loglik_and_mw", minval=minval)
        fill_lok_and_red(i1, i2, lok_and_mw)
      }
    }

    # Take care of individuals that are not good enough
    logliks[i1, which(logliks[i1,] < minval)] <- minval
    redundants[i1, which(logliks[i1,] <= minval)] <- M


    ## Natural selection and the reproduction pool ##
    if(length(unique(logliks[i1,])) == 1) { # If all individuals are the same, the surviving probability weight is 1.
      surviveProbs <- rep(1, popsize)
    } else {
      T_values <- logliks[i1,] + abs(min(logliks[i1,])) # Function T giving surviving weights, as described by Dorsey R. E. and Mayer W. J., 1995
      T_values <- T_values/(1 + regime_force_scale*redundants[i1,]) # Favour individuals with less redundant regimes
      surviveProbs <- T_values/sum(T_values) # The surviving probability weights
    }
    survivors <- sample(1:popsize, size=popsize, replace=TRUE, prob=surviveProbs)
    H <- G[,survivors] # Save the survivors in the reproduction pool H

    # Calculate mean, min, and max log-likelihoods of the survivors
    survivor_liks <- logliks[i1, survivors]
    survivor_redundants <- redundants[i1, survivors]
    meanLik <- mean(survivor_liks)
    maxLik <- max(survivor_liks)
    if(maxLik == meanLik) meanLik <- meanLik + 0.1 # We avoid dividing by zero when all the individuals are the same


    ## Cross-overs ##
    # Individually adaptive cross-over rate as described by  M. Srinivas L.M. Patnaik (1994) with the modification of
    # setting crossover rate to be at least 0.4 for all individuals (so that the best genes mix in the population too).
    indeces <- seq(1, popsize - 1, by=2)
    parentLiks <- vapply(indeces, function(i2) max(survivor_liks[i2], survivor_liks[i2 + 1]), numeric(1)) # Max of parents
    coRate <- vapply(1:length(indeces), function(i2) max(min(1, (maxLik - parentLiks[i2])/(maxLik - meanLik)), 0.4), numeric(1))

    # For each pair of individuals, draw from Bernoulli distribution whether to crossover or not. If yes, do it.
    which_co <- rbinom(popsize/2, 1, coRate) # 1 = do crossover, 0 = do not crossover
    I <- round(runif(popsize/2, min=0.5 + 1e-16, max=nrow(H) - 0.5 - 1e-16)) # Break points
    indeces <- seq(1, popsize - 1, by=2)
    H2 <- as.vector(vapply(1:(popsize/2), function(i2) {
      i3 <- indeces[i2]
      if(which_co[i2] == 1) {
        c(c(H[1:I[i2], i3], H[(I[i2] + 1):d, i3 + 1]), c(H[1:I[i2], i3 + 1], H[(I[i2] + 1):d, i3])) # Crossover according to I
      } else {
        c(H[,i3], H[,i3 + 1])
      }
    }, numeric(2*d)))
    H2 <- matrix(H2, ncol=popsize)

    # Get the best individual so far and check for reduntant regimes
    best_index0 <- which(logliks == max(logliks), arr.ind = TRUE)
    best_index <- best_index0[order(best_index0[,1], decreasing=FALSE)[1],] # The first generation when the best loglik occurred (take the first because of fitness inheritance)
    best_ind <- generations[, best_index[2], best_index[1]]
    best_mw <- mixingWeights_int(data=data, p=p, M=M_orig, params=best_ind, model=model, restricted=restricted,
                                 constraints=constraints, parametrization="mean", checks=FALSE, to_return="mw")
    which_redundant <- which((vapply(1:M, function(i2) sum(best_mw[,i2] > red_criteria[1]) < red_criteria[2]*length(data), logical(1))))

    # Keep track of "the alternative best individual" which has (weakly) less reduntant regimes than the current best one.
    if(length(which_redundant) <= length(which_redundant_alt)) {
      alt_ind <- best_ind
      which_redundant_alt <- which_redundant
    }


    ## Mutations ##
    # Mutation rates
    which_not_co <- rep(1 - which_co, each=2)
    if(abs(maxLik - meanLik) <= abs(0.03*meanLik)) {
      muRate <- rep(0.7, popsize) # Massive mutations if converging
    } else {
      # Individually adaptive mutation rates, Patnaik and Srinivas (1994); we only mutate those who did not crossover.
      muRate <- 0.5*vapply(1:popsize, function(i2) min(which_not_co[i2], (maxLik - survivor_liks[i2])/(maxLik - meanLik)), numeric(1))
    }

    # Which individuals to mutate
    mutate_indicator <- rbinom(n=popsize, size=1, prob=muRate)
    mutate <- which(mutate_indicator == 1)

    # How to generate random AR parameters
    if(!is.null(constraints) | runif(1) > 0.3) { # From normal distribution (the only option if constraints are imposed)
      forcestat <- FALSE
    } else { # Use the algorithm by Monahan (1984) to generate stationary AR parameters (slower)
      forcestat <- TRUE
    }

    # Mutate
    if(length(mutate) >= 1) {
      if(i1 > smartMu) { # Smart mutations

        # Smart mutate more if there are redundant (wasted) regimes
        if(length(which_redundant) > 0) {
          muRate <- vapply(1:popsize, function(i2) which_not_co[i2]*max(0.1, muRate[i2]), numeric(1))
          mutate0 <- which(rbinom(n=popsize, size=1, prob=muRate) == 1)
          if(length(mutate0) > length(mutate)) {
            mutate <- mutate0
          }
        }

        # Mutating accuracy
        accuracy <- abs(rnorm(length(mutate), mean=15, sd=10))

        ## 'Smart mutation': mutate close to a well fitting individual. Redundant regimes are obviously not smart
        # mutated but drawn at random ('rand_to_use' in what follows).
        if(Cquals == FALSE | length(which_redundant) <= length(which_redundant_alt) | runif(1) > 0.5) {
          # Smart mutate to alt_ind which is the best fitting individual with the least redundant regimes.
          # Note that best_ind == alt_ind when length(which_redundant) <= length(which_redundant_alt).
          ind_to_use <- alt_ind
          rand_to_use <- which_redundant_alt
        } else {
          # Alternatively, if there exists an alternative individual with strictly less redundant regimes
          # than in the best_ind, a "regime combining procedure" might take place: take a redundant regime
          # of the best_ind and replace it with a nonredundant regime taken from the alt_ind. Then do smart
          # mutation close to this new individual. For simplicity, regime combining is not considered for
          # models imposing linear constraints if the constraints are not the same for all the regimes.

          # We take a redundant regime of the best_ind ('which_to_change') and replace it with a nonredundant
          # regime taken from the alt_ind. We want to take such nonredundant regime from alt_ind that is not
          # similar to the nonredundant regimes of best_ind. In order to choose such regime, we do a trade-off:
          # we remove the degrees of freedom parameters from all StMAR type regimes and compare all the nonredundant
          # regimes of best_ind to all nonredundant regimes of alt_ind without accounting for the degree of freedom
          # parameters (and we do cross-type comparisons for G-StMAR model). After choosing a regime from alt_ind
          # based on the comparison, we add a dfs parameter if a StMAR type regime is to be changed.
          # This way we avoid the problem that similar regimes with large degrees of freedom parameters, say, 2000
          # and 5000, look like non-similar regimes in numerical distance comparison. By doing cross-type comparisons
          # in the G-StMAR model we also avoid the problem that a StMAR type regime with large dfs parameter could be
          # inserted to a parameter vector already containing a similar GMAR type regime. We loose in the trade-off
          # reliability of the comparisons between regimes that are otherwise similar but have different dfs (say, 3 and 20).

          # The first redundant regime of best_ind is to be replaced with a regime from alt_ind
          which_to_change <- which_redundant[1]

          # All nonredundant regimes of best_ind (without dfs)
          best_ind_nonRedRegimes <- sapply((1:M)[-which_redundant], function(i2) extractRegime(p=p, M=M_orig, params=best_ind, model=model,
                                                                                               restricted=restricted, constraints=constraints,
                                                                                               regime=i2, with_dfs=FALSE))
          # All nonredundant regimes of alt_ind (without dfs)
          if(length(which_redundant_alt) == 0) { # Special case for techinal reason
            nonred_altind <- 1:M
          } else {
            nonred_altind <- (1:M)[-which_redundant_alt]
          }
          alt_ind_nonRedRegimes <- sapply(nonred_altind, function(i2) extractRegime(p=p, M=M_orig, params=alt_ind, model=model,
                                                                                    restricted=restricted, constraints=constraints,
                                                                                    regime=i2, with_dfs=FALSE))

          # Calculate "distances" (parameters values are scaled to the same maginitude) between the nonredundant
          # best_ind and alt_ind regimes.
          # A row for each nonredundant best_ind regime and a column for each nonredundant alt_ind regime
          dist_to_regime <- matrix(nrow=ncol(best_ind_nonRedRegimes), ncol=ncol(alt_ind_nonRedRegimes))
          for(i2 in 1:nrow(dist_to_regime)) {
            dist_to_regime[i2,] <- vapply(1:ncol(dist_to_regime), function(i3) regime_distance(regime_pars1=best_ind_nonRedRegimes[,i2],
                                                                                               regime_pars2=alt_ind_nonRedRegimes[,i3]), numeric(1))
          }

          # Which alt_ind regime, i.e., column should be used? Choose the one with the largest 'distance' to the closest regime
          # to avoid dublicating similar regimes.
          which_reg_to_use <- which(apply(dist_to_regime, 2, min) == max(apply(dist_to_regime, 2, min)))[1]

          if(model == "StMAR" | model == "GMAR") {
            reg_to_use <- extractRegime(p=p, M=M_orig, params=alt_ind, model=model, restricted=restricted, constraints=constraints,
                                        regime=which_reg_to_use, with_dfs=TRUE)
          } else { # model == "G-StMAR"
            if(which_to_change <= M1) { # No dfs
              reg_to_use <- alt_ind_nonRedRegimes[,which_reg_to_use]
            } else { # which_to_change > M1, so we need dfs
              if(which_reg_to_use <= M1) { # We need to create new dfs
                reg_to_use <- c(alt_ind_nonRedRegimes[,which_reg_to_use], runif(n=1, min=2+1e-6, max=30))
              } else { # The chosen alt_ind regime contains dfs
                reg_to_use <- extractRegime(p=p, M=M_orig, params=alt_ind, model=model, restricted=restricted, constraints=constraints,
                                            regime=which_reg_to_use, with_dfs=TRUE)
              }
            }
          }

          # Combine the regimes to a complete parameter vector
          ind_to_use <- changeRegime(p=p, M=M_orig, params=best_ind, model=model, restricted=restricted, constraints=constraints,
                                     regimeParams=reg_to_use, regime=which_to_change)

          # If there are redundant regimes left, they should be random mutated and not smart mutated
          rand_to_use <- which_redundant[which_redundant != which_to_change]
        }
        # Do smart mutations close to ind_to_use with random regimes (if any) given by rand_to_use
        H2[,mutate] <- vapply(1:length(mutate), function(i3) smartIndividual_int(p=p, M=M_orig, params=ind_to_use, model=model,
                                                                                 restricted=restricted, constraints=constraints,
                                                                                 meanscale=meanscale, sigmascale=sigmascale,
                                                                                 accuracy=accuracy[i3], whichRandom=rand_to_use,
                                                                                 forcestat=forcestat), numeric(d))
      } else { # Random mutations
        H2[,mutate] <- vapply(1:length(mutate), function(i3) randomIndividual_int(p, M_orig, model=model, restricted=restricted,
                                                                                  constraints=constraints, meanscale=meanscale,
                                                                                  sigmascale=sigmascale, forcestat=forcestat), numeric(d))
      }
    }

    # Sort components if constraints are not used
    if(is.null(constraints)) {
      H2 <- vapply(1:popsize, function(i2) sortComponents(p, M_orig, H2[,i2], model=model, restricted=restricted), numeric(d))
    }

    # Set up for the next generation
    G <- H2
  }

  # The return value
  if(to_return == "best_ind") {
    ret <- best_ind
  } else {
    ret <- alt_ind
  }

  if(parametrization == "intercept") {
    ret <- change_parametrization(p=p, M=M_orig, params=ret, model=model, restricted=restricted, constraints=constraints, change_to="intercept")
  }
  ret
}



#' @title Calculate "distance" between two regimes
#'
#' @description \code{regime_distance} scales each regime parameter to the same magnitude
#'  and then calculates distance between scaled \code{regime_pars1} and \code{regime_pars2}.
#'
#' @param regime_pars1 a numeric vector representing a regime.
#' @param regime_pars2 a numeric vector representing another regime, same length as \code{regime_pars1}
#' @return Returns "distance" between \code{regime_pars1} and \code{regime_pars2}. Values are scaled
#'   to the same magnitude before calculating the "distance". Read the source code for details.

regime_distance <- function(regime_pars1, regime_pars2) {
  scale_fun <- function(x) {
    x <- abs(x)
    if(x < 1) {
      return(1)
    } else {
      return(10^ceiling(abs(log10(x))))
    }
  }
  scales1 <- vapply(regime_pars1, scale_fun, numeric(1))
  scales2 <- vapply(regime_pars2, scale_fun, numeric(1))
  sqrt(crossprod(regime_pars1/scales1 - regime_pars2/scales2))
}



#' @title Create random regime parameters
#'
#' @description \code{random_regime} generates random regime parameters.
#'
#' @inheritParams GAfit
#' @param forcestat use the algorithm by Monahan (1984) to force stationarity on the AR parameters (slower)?
#'   Not supported for constrained models.
#' @param m which regime? This is required for models with constraints for which a list of possibly differing
#'   constraint matrices is provided.
#' @return \describe{
#'   \item{Regular models:}{\strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m}^2)}
#'   where \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p})}.}
#'   \item{Restricted models:}{Not supported!}
#'   \item{Constrained models:}{Replace the vectors \strong{\eqn{\phi_{m}}} with vectors \strong{\eqn{\psi_{m}}}.}
#' }
#' @inherit random_arcoefs details references

random_regime <- function(p, meanscale, sigmascale, restricted=FALSE, constraints=NULL, m, forcestat=FALSE) {
  stopifnot(restricted == FALSE)
  if(!is.null(constraints)) {
    C0 <- as.matrix(constraints[[m]])
    q <- ncol(C0)
    scale <- sum(abs(C0))
    return(c(rnorm(n=1, mean=meanscale[1], sd=meanscale[2]), rnorm(n=q, mean=0, sd=0.6/scale), abs(rnorm(n=1, mean=0, sd=sigmascale))))
  } else {
    return(c(rnorm(n=1, mean=meanscale[1], sd=meanscale[2]), random_arcoefs(p=p, forcestat=forcestat, sd=0.6/p), abs(rnorm(n=1, mean=0, sd=sigmascale))))
  }
}

#' @title Create random AR coefficients
#'
#' @description \code{random_arcoefs} generates random AR coefficients.
#'
#' @inheritParams GAfit
#' @param forcestat use the algorithm by Monahan (1984) to force stationarity on the AR parameters (slower)?
#' @param sd if \code{forcestat==FALSE}, then AR parameters are drawn from zero mean normal distribution with sd given by this parameter.
#' @details If \code{forcestat==TRUE}, then the AR coefficients are relatively large, otherwise they are usually relatively small.
#' @return Returns \eqn{px1} vector containing random AR coefficients.
#' @references
#'  \itemize{
#'    \item Monahan J.F. 1984. A Note on Enforcing Stationarity in Autoregressive-Moving Average Models.
#'          \emph{Biometrica} \strong{71}, 403-404.
#'  }

random_arcoefs <- function(p, forcestat=FALSE, sd=0.6/p) {
  if(forcestat == TRUE) { # Algorithm by Mohanan 1984
    r_p <- runif(p, min=-0.90, max=0.90) # Random partial autocorrelations
    ymat <- matrix(0, nrow=p, ncol=p)
    for(k in 1:p) {
      if(k > 1) {
        for(i1 in 1:(k - 1)) {
          ymat[i1, k] <- ymat[i1, k - 1] + r_p[k]*ymat[k - i1, k - 1]
        }
      }
      ymat[k, k] <- r_p[k]
    }
    return(-ymat[,p])
  } else { # Draw from normal distribution
    return(rnorm(n=p, mean=0, sd=sd))
  }
}



#' @title Add random dfs to a vector
#'
#' @description \code{add_dfs} adds random degrees of freedom parameters to a vector.
#'
#' @param x a vector to add the dfs to
#' @param how_many how many dfs?
#' @details Read the source code for details.
#' @return Returns \code{c(x, dfs)} with \code{how_many} dfs-elements.

add_dfs <- function(x, how_many) {
  c(x, 2 + rgamma(how_many, shape=0.3, rate=0.007))
}



#' @title Create random GMAR, StMAR, or G-StMAR model compatible parameter vector
#'
#' @description \code{randomIndividual_int} creates a random GMAR, StMAR, or G-StMAR model compatible parameter vector.
#'
#' \code{smartIndividual_int} creates a random GMAR, StMAR, or G-StMAR model compatible parameter vector close to argument \code{params}.
#'
#' @inheritParams GAfit
#' @inheritParams random_regime
#' @param meanscale a real valued vector of length two specifying the mean (the first element) and standard deviation (the second element)
#'  of the normal distribution from which the \eqn{\phi_{m,0}} \strong{or} \eqn{\mu_{m}} (depending on the desired parametrization)
#'  parameters (for random regimes) should be generated.
#' @param sigmascale a positive real number specifying the standard deviation of the (zero mean, positive only by taking absolute value)
#'  normal distribution from which the component variance parameters (for random regimes) should be generated.
#' @inherit GAfit return
#' @inherit random_regime references

randomIndividual_int <- function(p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, meanscale, sigmascale, forcestat=FALSE) {
  model <- match.arg(model)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  }

  if(restricted == FALSE) {
    ind <- unlist(lapply(1:M, function(m) random_regime(p=p, meanscale=meanscale, sigmascale=sigmascale,
                                                        restricted=restricted, constraints=constraints,
                                                        m=m, forcestat=forcestat)))
  } else { # If restricted == TRUE
    if(is.null(constraints)) {
      q <- p
      scale <- p
    } else {
      C0 <- as.matrix(constraints)
      q <- ncol(C0)
      scale <- sum(abs(C0))
    }
    ind <- c(rnorm(n=M, mean=meanscale[1], sd=meanscale[2]), random_arcoefs(p=q, forcestat=forcestat, sd=0.6/scale), abs(rnorm(n=M, mean=0, sd=sigmascale)))
  }

  if(M > 1) {
    alphas <- runif(n=M)
    if(model != "G-StMAR") {
      alphas <- alphas[order(alphas, decreasing=TRUE)]/sum(alphas)
    } else { # model == "G-StMAR
      alphas <- alphas/sum(alphas)
      alphasM1 <- alphas[1:M1]
      alphasM2 <- alphas[(M1 + 1):M]
      alphas <- c(alphasM1[order(alphasM1, decreasing=TRUE)], alphasM2[order(alphasM2, decreasing=TRUE)])
    }
    alphas <- alphas[-M]
    ind <- c(ind, alphas)
  }

  if(model == "StMAR") {
    ind <- add_dfs(ind, how_many=M)
  } else if(model == "G-StMAR") {
    ind <- add_dfs(ind, how_many=M2)
  }
  ind
}


#' @title Create random GMAR, StMAR, or G-StMAR model compatible parameter vector
#'
#' @description \code{randomIndividual} creates a random GMAR, StMAR, or G-StMAR model compatible mean-parametrized parameter vector.
#'
#' \code{smartIndividual} creates a random GMAR, StMAR, or G-StMAR model compatible parameter vector close to argument \code{params}.
#'   Sometimes returns exactly the given parameter vector.
#'
#' @inheritParams randomIndividual_int
#' @inherit randomIndividual_int return references
#' @details These functions can be used, for example, to create initial populations for the genetic algorithm. Mean-parametrization
#'   (instead of intercept terms \eqn{\phi_{m,0}}) is assumed.
#' @examples
#' # GMAR model parameter vector
#' params22 <- randomIndividual(2, 2, meanscale=c(0, 1), sigmascale=1)
#' smart22 <- smartIndividual(2, 2, params22, accuracy=10)
#' cbind(params22, smart22)
#'
#'
#' # Restricted GMAR parameter vector
#' params12r <- randomIndividual(1, 2, restricted=TRUE, meanscale=c(-2, 2), sigmascale=2)
#' smart12r <- smartIndividual(1, 2, params12r, restricted=TRUE, accuracy=20)
#' cbind(params12r, smart12r)
#'
#'
#' # StMAR parameter vector: first regime is random in the "smart individual"
#' params13t <- randomIndividual(1, 3, model="StMAR", meanscale=c(3, 1), sigmascale=3)
#' smart13t <- smartIndividual(1, 3, params13t, model="StMAR", accuracy=15,
#'                             meanscale=c(3, 3), sigmascale=3, whichRandom=1)
#' cbind(params13t, smart13t)
#'
#'
#' # Restricted StMAR parameter vector
#' params22tr <- randomIndividual(2, 2, model="StMAR", restricted=TRUE,
#'                                meanscale=c(3, 2), sigmascale=0.5)
#' smart22tr <- smartIndividual(2, 2, params22tr, model="StMAR", restricted=TRUE,
#'                              accuracy=30)
#' cbind(params22tr, smart22tr)
#'
#'
#' # G-StMAR parameter vector
#' params12gs <- randomIndividual(1, c(1, 1), model="G-StMAR", meanscale=c(0, 1),
#'                                sigmascale=1)
#' smart12gs <- smartIndividual(1, c(1, 1), params12gs, model="G-StMAR", accuracy=20)
#' cbind(params12gs, smart12gs)
#'
#'
#' # Restricted G-StMAR parameter vector
#' params23gsr <- randomIndividual(2, c(1, 2), model="G-StMAR", restricted=TRUE,
#'                                 meanscale=c(-1, 1), sigmascale=3)
#' smart23gsr <- smartIndividual(2, c(1, 2), params23gsr, model="G-StMAR", restricted=TRUE,
#'                               meanscale=c(0, 1), sigmascale=1, accuracy=20, whichRandom=2)
#' cbind(params23gsr, smart23gsr)
#'
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' C <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- randomIndividual(2, 2, constraints=C, meanscale=c(1, 1),
#'                               sigmascale=1)
#' smart22c <- smartIndividual(2, 2, params22c, constraints=C, accuracy=10)
#' cbind(params22c, smart22c)
#'
#'
#' # Such constrained StMAR(3, 2) model that the second order AR coefficients
#' # are constrained to zero.
#' C0 = matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
#' C = list(C0, C0)
#' params32c <- randomIndividual(3, 2, model="StMAR", constraints=C,
#'                               meanscale=c(1, 1), sigmascale=1)
#' smart32c <- smartIndividual(3, 2, params32c, model="StMAR", constraints=C, accuracy=10)
#' cbind(params32c, smart32c)
#'
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero. Second regime is random in the "smart individual".
#' params32trc <- randomIndividual(3, 2, model="StMAR", restricted=TRUE,
#'                                 constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2),
#'                                 meanscale=c(-2, 0.5), sigmascale=4)
#' smart32trc <- smartIndividual(3, 2, params32trc, model="StMAR", restricted=TRUE,
#'                               constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2),
#'                               meanscale=c(0, 0.1), sigmascale=0.1, whichRandom=2,
#'                               accuracy=20)
#' cbind(params32trc, smart32trc)
#' @export

randomIndividual <- function(p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, meanscale, sigmascale, forcestat=FALSE) {
  model <- match.arg(model)
  check_model(model)
  checkPM(p=p, M=M, model=model)
  if(length(meanscale) != 2) {
    stop("meanscale is wrong dimension")
  } else if(length(sigmascale) != 1) {
    stop("sigmascale is wrong dimension")
  } else if(meanscale[2] <= 0) {
    stop("The second element of meanscale should be larger than zero")
  } else if(sigmascale <= 0) {
    stop("sigmascale should be larger than zero")
  }
  if(!is.null(constraints)) {
    checkConstraintMat(p=p, M=M, restricted=restricted, constraints=constraints)
  }
  for(i1 in 1:42) {
    ret <- randomIndividual_int(p=p, M=M, model=model, restricted=restricted, constraints=constraints, meanscale=meanscale,
                                sigmascale=sigmascale, forcestat=forcestat)
    ret0 <- reformConstrainedPars(p=p, M=M, params=ret, model=model, restricted=restricted, constraints=constraints)
    if(isStationary_int(p=p, M=M, params=ret0, restricted=restricted)) return(ret)
  }
  message("Failed to generate stationary parameter vector")
  ret
}


#' @rdname randomIndividual_int
#' @inheritParams loglikelihood
#' @param accuracy a real number larger than zero specifying how close to \code{params} the generated parameter vector should be.
#'   Standard deviation of the normal distribution from which new parameter values are drawn from will be corresponding parameter
#'   value divided by \code{accuracy}.
#' @param whichRandom a numeric vector of maximum length \code{M} specifying which regimes should be random instead of "smart" when
#' using \code{smartIndividual}. Does not affect mixing weight parameters. Default in none.
#' @param forcestat use the algorithm by Monahan (1984) to force stationarity on the AR parameters (slower) for random regimes?
#'   Not supported for constrained models.
#' @inherit random_regime references

smartIndividual_int <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                                meanscale, sigmascale, accuracy, whichRandom, forcestat=FALSE) {
  model <- match.arg(model)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  }
  if(missing(whichRandom)) {
    whichRandom <- numeric(0)
  }

  if(restricted == FALSE) {
    ind <- c()
    j <- 0
    for(i1 in 1:M) { # Run through components
      q <- ifelse(is.null(constraints), p, ncol(as.matrix(constraints[[i1]]))) # p for non-constrained models
      if(any(whichRandom == i1)) { # If regime should be random
        ind <- c(ind, random_regime(p=p, meanscale=meanscale, sigmascale=sigmascale,
                                    restricted=restricted, constraints=constraints,
                                    m=i1, forcestat=forcestat))
      } else {
        params_m <- params[(j + 1):(j + q + 1)] # no sigma, alphas or dfs
        sigma_m <- params[j + q + 2]
        ind <- c(ind, rnorm(n=q + 1, mean=params_m, sd=abs(params_m/accuracy)), abs(rnorm(n=1, mean=sigma_m, sd=sigma_m/accuracy)))
      }
      j <- j + q + 2
    }
  } else { # If restricted == TRUE
    q <- ifelse(is.null(constraints), p, ncol(as.matrix(constraints))) # p for non-constrained models
    params0 <- params[1:(M + q)]
    sigmas <- params[(M + q + 1):(2*M + q)]
    if(length(whichRandom) == 0) {
      ind <- c(rnorm(n=M + q, mean=params0, sd=abs(params0/accuracy)), abs(rnorm(n=M, mean=sigmas, sd=sigmas/accuracy)))
    } else { # If some regime(s) should be random
      ind <- numeric(0)
      for(i1 in 1:M) { # Generate phi0
        if(any(whichRandom == i1)) {
          ind <- c(ind, rnorm(n=1, mean=meanscale[1], sd=meanscale[2]))
        } else {
          ind <- c(ind, rnorm(n=1, mean=params0[i1], sd=abs(params0[i1]/accuracy)))
        }
      }
      ind <- c(ind, rnorm(q, mean=params0[(M + 1):(M + q)], sd=abs(params0[(M + 1):(M + q)]/accuracy)))
      for(i1 in 1:M) { # Generate sigmas
        if(any(whichRandom == i1)) {
          ind <- c(ind, abs(rnorm(n=1, mean=0, sd=sigmascale)))
        } else {
          ind <- c(ind, rnorm(n=1, mean=sigmas[i1], sd=abs(sigmas[i1]/accuracy)))
        }
      }
    }
  }

  if(M > 1) { # Need to add alphas
    if(restricted == FALSE) {
      alphas <- params[(j + 1):(j + M - 1)]
    } else {
      alphas <- params[(q + 2*M + 1):(3*M + q - 1)]
    }
    all_alphas <- c(alphas, 1 - alphas)
    alphas2 <- abs(rnorm(n=M, mean=all_alphas, sd=0.2))
    alphas2 <- alphas2/sum(alphas2)
    ind <- c(ind, alphas2[-M])
  }

  if(model == "StMAR" | model == "G-StMAR") { # Need to add dfs
    d <- length(params)
    if(model == "StMAR") {
      M2 <- M # So we can just use M2 and M1 for StMAR models too
      M1 <- 0
    }
    dfs <- params[(d - M2 + 1):d]
    if(length(whichRandom) == 0) {
      dfs2 <- rnorm(n=M2, mean=dfs, sd=dfs/accuracy)
    } else { # If some regimes should be random
      dfs2 <- numeric(0)
      for(i1 in (M1 + 1):M) { # Generate dfs
        if(any(whichRandom == i1)) {
          dfs2 <- add_dfs(dfs2, how_many=1)
        } else {
          dfs2 <- c(dfs2, rnorm(n=1, mean=dfs[i1 - M1], sd=abs(dfs[i1 - M1]/accuracy)))
        }
      }
    }
    dfs2[dfs2 <= 2] <- 2.1
    ind <- c(ind, dfs2)
  }

  ind_stand <- reformConstrainedPars(p=p, M=M_orig, params=ind, model=model, restricted=restricted, constraints=constraints)
  if(isStationary_int(p=p, M=M, params=ind_stand, restricted=restricted)) {
    return(ind)
  } else {
    return(params) # Return the given individual if smart mutation is not stationary
  }
}


#' @rdname randomIndividual
#' @inheritParams smartIndividual_int
#' @export

smartIndividual <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                            meanscale, sigmascale, accuracy, whichRandom=numeric(0), forcestat=FALSE) {
  model <- match.arg(model)
  check_model(model)
  checkPM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  if(accuracy <= 0) {
    stop("Argument accuracy has to be larger than zero")
  } else if(length(whichRandom) > 0) {
    if(any(!whichRandom %in% 1:sum(M))) {
      stop("The elements of whichRandom should be positive integers in the clsoed interval [1, M]")
    } else if(length(whichRandom) > 0) {
      if(length(meanscale) != 2) {
        stop("meanscale is wrong dimension")
      } else if(length(sigmascale) != 1) {
        stop("sigmascale is wrong dimension")
      } else if(meanscale[2] <= 0) {
        stop("The second element of meanscale should be larger than zero")
      } else if(sigmascale <= 0) {
        stop("sigmascale should be larger than zero")
      }
    }
  }
  if(!is.null(constraints)) {
    checkConstraintMat(p=p, M=M, restricted=restricted, constraints=constraints)
  }
  smartIndividual_int(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                      meanscale=meanscale, sigmascale=sigmascale, accuracy=accuracy, whichRandom=whichRandom,
                      forcestat=forcestat)
}


#' @title Extract regime from a parameter vector
#'
#' @description \code{extractRegime} extracts the specified regime from the GMAR, StMAR, or G-StMAR model parameter vector.
#'  Doesn't extract mixing weight parameter \eqn{\alpha}.
#' @inheritParams loglikelihood
#' @param regime a positive integer in the interval [1, M] defining which regime should be extracted.
#' @param with_dfs Should the degrees of freedom parameter (if any) be included?
#' @return Returns a numeric vector corresponding to the regime with...
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(p+2x1)} vector \eqn{(\phi_{m,0},\phi_{m,1},...,\phi_{m,p}, \sigma_{m}^2)}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(p+3x1)} vector \eqn{(\phi_{m,0},\phi_{m,1},...,\phi_{m,p}, \sigma_{m}^2, \nu_{m})}.}
#'        \item{For \strong{G-StMAR} model:}{Same as GMAR for GMAR-components and same as StMAR for StMAR-components.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as described above, but vector \strong{\eqn{\phi_{m}}} replaced with
#'         vector \strong{\eqn{\psi_{m}}} that satisfies \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(2x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2)}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(3x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2, \nu_{m})}.}
#'        \item{For \strong{G-StMAR} model:}{Same as GMAR for GMAR-components and same as StMAR for StMAR-components.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as described above.}
#'      }
#'    }
#'  }

extractRegime <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                          regime, with_dfs=TRUE) {
  model <- match.arg(model)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  }
  if(restricted == FALSE) {
    j <- 0 # Indicates where we at
    for(i1 in 1:M) { # Go through regimes
      if(!is.null(constraints)) {
        q <- ncol(as.matrix(constraints[[i1]]))
      } else {
        q <- p
      }
      if(i1 == regime) { # Take the parameters
        params0 <- params[(j + 1):(j + q + 2)]
      }
      j <- j + q + 2
    }
    if(model == "StMAR" & with_dfs) {
      params0 <- c(params0, params[j + M - 1 + regime]) # dfs
    } else if(model == "G-StMAR" & with_dfs) {
      if(regime > M1) {
        params0 <- c(params0, params[j + M - 1 + regime - M1]) # dfs
      }
    }
    return(params0)
  } else { # If restricted == TRUE
    if(!is.null(constraints)) {
      q <- ncol(as.matrix(constraints))
    } else {
      q <- p
    }
    params0 <- c(params[regime], params[M + q + regime])
    if(model == "StMAR" & with_dfs) {
      params0 <- c(params0, params[3*M + q - 1 + regime]) # dfs
    } else if(model == "G-StMAR" & with_dfs) {
      if(regime > M1) {
        params0 <- c(params0, params[3*M + q - 1 + regime - M1]) # dfs
      }
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
#'        \item{With \strong{linear constraints}:}{Parameter vector as described above, but vector \strong{\eqn{\phi_{m}}} replaced with
#'         vector \strong{\eqn{\psi_{m}}} that satisfies \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(2x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2)}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(3x1)} vector \eqn{(\phi_{m,0}, \sigma_{m}^2, \nu_{m})}.}
#'        \item{For \strong{G-StMAR} model:}{Same as GMAR for GMAR-components and same as StMAR for StMAR-components.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as described above.}
#'      }
#'    }
#'  }
#' @param regime a positive integer in the interval [1, M] defining which regime should be changed.
#' @return Returns modified parameter vector of the form described in \code{params}.

changeRegime <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, regimeParams, regime) {
  model <- match.arg(model)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  } else {
    M1 <- 0 # Exists for tidier code below
  }
  if(restricted == FALSE) {
    params0 <- numeric(0)
    j <- 0 # Indicates where we at
    for(i1 in 1:M) { # Go through regimes
      q <- ifelse(is.null(constraints), p, ncol(as.matrix(constraints[[i1]])))
      if(i1 == regime) { # Change the parameters
        if(model == "StMAR" | (model == "G-StMAR" & regime > M1)) {
          regimeDfs <- regimeParams[length(regimeParams)]
          regimeParams <- regimeParams[-length(regimeParams)] # Delete dfs
        }
        params0 <- c(params0, regimeParams)
      } else { # Use same parameters
        params0 <- c(params0, params[(j + 1):(j + q + 2)])
      }
      j <- j + q + 2
    }
    if(M > 1) {
      params0 <- c(params0, params[(j + 1):(j + M - 1)]) # Add alphas
    }
    if(model == "StMAR" | model == "G-StMAR") {
      dfs <- params[(j + M):(j + 2*M - 1 - M1)]
      if(model == "StMAR" | (model == "G-StMAR" & regime > M1)) {
        dfs[regime - M1] <- regimeDfs # Change the dfs
      }
      params0 <- c(params0, dfs) # Add dfs
    }
    return(params0)
  } else { # Restricted == TRUE
    q <- ifelse(is.null(constraints), p, ncol(as.matrix(constraints)))
    params0 <- params
    params0[regime] <- regimeParams[1] # phi0
    params0[M + q + regime] <- regimeParams[2] # sigma^2
    if(model == "StMAR") {
      params0[3*M + q - 1 + regime] <- regimeParams[3] # dfs
    } else if(model == "G-StMAR" & regime > M1) {
      params0[3*M + q - 1 + regime - M1] <- regimeParams[3] # dfs
    }
    return(params0)
  }
}
