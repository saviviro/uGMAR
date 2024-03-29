#' @title Calculate regime specific means \eqn{\mu_{m}}
#'
#' @description \code{get_regime_means} calculates the regime means \eqn{\mu_{m} = \phi_{m,0}/(1-\sum\phi_{i,m})}
#'   for the given GMAR, StMAR, or G-StMAR model
#'
#' @inheritParams add_data
#' @return Returns a length \code{M} vector containing the regime mean \eqn{\mu_{m}} in the m:th element.
#' @inherit is_stationary references
#' @family moment functions
#' @seealso \code{\link{cond_moments}}, \code{\link{uncond_moments}}, \code{\link{get_regime_vars}},
#'  \code{\link{get_regime_autocovs}}
#' @examples
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' gmar13 <- GSMAR(p=1, M=3, params=params13, model="GMAR")
#' get_regime_means(gmar13)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12t <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
#' get_regime_means(stmar12t)
#'
#' # G-StMAR model (similar to the StMAR model above)
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(p=1, M=c(1, 1), params=params12gs, model="G-StMAR")
#' get_regime_means(gstmar12)
#' @export

get_regime_means <- function(gsmar) {
  check_gsmar(gsmar)
  p <- gsmar$model$p
  M <- gsmar$model$M
  params <- gsmar$params
  model <- gsmar$model$model
  restricted <- gsmar$model$restricted
  constraints <- gsmar$model$constraints

  if(gsmar$model$parametrization == "intercept") {
    params <- change_parametrization(p=p, M=M, params=params, model=model, restricted=restricted,
                                     constraints=constraints, change_to="mean")
  }
  pick_phi0(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
}


#' @title Calculate regime specific autocovariances \strong{\eqn{\gamma}}\eqn{_{m,p}}
#'
#' @description \code{get_regime_autocovs} calculates the first p regime specific autocovariances \strong{\eqn{\gamma}}\eqn{_{m,p}}
#'   for the given GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams add_data
#' @return Returns a size \eqn{(pxM)} matrix containing the first p autocovariances of the components processes:
#'  i:th autocovariance in the i:th row and m:th component process in the m:th column.
#' @family moment functions
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'            \emph{Journal of Time Series Analysis}, \strong{36}(2), 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2023. A mixture autoregressive model based on Student's t-distribution.
#'          \emph{Communications in Statistics - Theory and Methods}, \strong{52}(2), 499-515.
#'    \item Virolainen S. 2022. A mixture autoregressive model based on Gaussian and Student's t-distributions.
#'          Studies in Nonlinear Dynamics & Econometrics, \strong{26}(4) 559-580.
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis. \emph{Springer}.
#'  }
#' @examples
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' gmar13 <- GSMAR(p=1, M=3, params=params13, model="GMAR")
#' get_regime_autocovs(gmar13)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12t <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
#' get_regime_autocovs(stmar12t)
#'
#' # G-StMAR model (similar to the StMAR model above)
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(p=1, M=c(1, 1), params=params12gs, model="G-StMAR")
#' get_regime_autocovs(gstmar12)
#' @export

get_regime_autocovs <- function(gsmar) {
  # Pick the relevant statistics
  check_gsmar(gsmar)
  p <- gsmar$model$p
  M <- gsmar$model$M
  pars <- pick_pars(p=p, M=M, params=gsmar$params, model=gsmar$model$model,
                    restricted=gsmar$model$restricted, constraints=gsmar$model$constraints)
  ret <- matrix(nrow=p, ncol=sum(M))

  # Go through the regimes
  for(i1 in 1:sum(M)) { # Formula from Lütkepohl (2005) eq. (2.1.39)
    phi <- pars[2:(p + 1), i1]
    if(p == 1) {
      A <- as.matrix(phi)
    } else {
      Ip1 <- diag(1, nrow=p - 1, ncol=p - 1)
      ZER <- matrix(0, nrow=p - 1, ncol=1)
      A <- rbind(phi, cbind(Ip1, ZER))
    }
    Sigma <- as.matrix(c(pars[p + 2, i1], rep(0, p^2 - 1)))
    Gamma <- matrix(solve(diag(1, nrow=p^2, ncol=p^2) - kronecker(A, A), Sigma), ncol=p, byrow=FALSE)
    gamma <- Gamma%*%as.matrix(phi)
    ret[, i1] <- gamma
  }
  ret
}


#' @title Calculate regime specific variances \eqn{\gamma_{m,0}}
#'
#' @description \code{get_regime_vars} calculates the unconditional regime specific variances \eqn{\gamma_{m,0}}
#'   for the given GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams add_data
#' @return Returns a length M vector containing the unconditional variances of the components processes:
#'   m:th element for the m:th regime.
#' @inherit get_regime_autocovs references
#' @family moment functions
#' @examples
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' gmar13 <- GSMAR(p=1, M=3, params=params13, model="GMAR")
#' get_regime_vars(gmar13)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12t <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
#' get_regime_vars(stmar12t)
#'
#' # G-StMAR model (similar to the StMAR model above)
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(p=1, M=c(1, 1), params=params12gs, model="G-StMAR")
#' get_regime_vars(gstmar12)
#' @export

get_regime_vars <- function(gsmar) {
  # Collect the relecant statistics
  reg_autocovs <- get_regime_autocovs(gsmar)
  p <- gsmar$model$p
  pars <- pick_pars(p=p, M=gsmar$model$M, params=gsmar$params, model=gsmar$model$model,
                    restricted=gsmar$model$restricted, constraints=gsmar$model$constraints)

  # Calculate the regime variances
  colSums(pars[-c(1, p + 2),]*reg_autocovs) + pars[p + 2,] # MPS 2021, eq.(4) and Theorem 1 (it is the same for GMAR and G-StMAR models also; see Virolainen 2020, Section 2.1)
}


#' @title Calculate unconditional mean, variance, and the first p autocovariances and autocorrelations
#'  of a GSMAR process.
#'
#' @description \code{uncond_moments_int} calculates the unconditional mean, variance, and the first p
#'  autocovariances and autocorrelations of the specified GSMAR process.
#'
#' @inheritParams loglikelihood_int
#' @details Differs from the function \code{uncond_moments} in arguments. This function exists for technical
#'  reasons only.
#' @return Returns a list containing the unconditional mean, variance, and the first p autocovariances and
#'  autocorrelations. Note that the lag-zero autocovariance/correlation is not included in the "first p"
#'  but is given in the \code{uncond_variance} component separately.
#' @inherit get_regime_autocovs references
#' @keywords internal

uncond_moments_int <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                              constraints=NULL, parametrization=c("intercept", "mean")) {
  # Checks and pick the relevant statistics etc
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  stopifnot(length(params) == n_params(p=p, M=M, model=model, restricted=restricted, constraints=constraints))
  gsmar <- structure(list(model=list(p=p, # "Pseudo gsmar" object for the helper functions
                                     M=M,
                                     model=model,
                                     restricted=restricted,
                                     constraints=constraints,
                                     parametrization=parametrization),
                          params=params),
                     class="gsmar")
  alphas <- pick_alphas(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params, model=gsmar$model$model,
                        restricted=gsmar$model$restricted, constraints=gsmar$model$constraints)
  reg_means <- get_regime_means(gsmar) # Regimewise unconditional means

  # Calculate the unconditional moments: KMS 2015, p.251, MPS 2021, p.6., Virolainen 2021 p.8.
  uncond_mean <- sum(alphas*reg_means)
  tmp <- sum(alphas*(reg_means - uncond_mean)^2)
  uncond_var <- sum(alphas*get_regime_vars(gsmar)) + tmp
  autocovs <- colSums(alphas*t(get_regime_autocovs(gsmar))) + tmp
  list(uncond_mean=uncond_mean,
       uncond_var=uncond_var,
       autocovs=autocovs,
       autocors=autocovs/uncond_var)
}



#' @title Calculate unconditional mean, variance, first p autocovariances and autocorrelations of the GSMAR process.
#'
#' @description \code{uncond_moments} calculates the unconditional mean, variance, and the first p autocovariances
#'  and autocorrelations of the GSMAR process.
#'
#' @inheritParams add_data
#' @inherit uncond_moments_int return references
#' @family moment functions
#' @examples
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' gmar13 <- GSMAR(p=1, M=3, params=params13, model="GMAR")
#' uncond_moments(gmar13)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12t <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
#' uncond_moments(stmar12t)
#'
#' # G-StMAR model (similar to the StMAR model above)
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(p=1, M=c(1, 1), params=params12gs, model="G-StMAR")
#' uncond_moments(gstmar12)
#' @export

uncond_moments <- function(gsmar) {
  check_gsmar(gsmar)
  uncond_moments_int(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params, model=gsmar$model$model,
                     restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                     parametrization=gsmar$model$parametrization)
}
