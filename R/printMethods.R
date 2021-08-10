#' @title Function factory for formatting values
#'
#' @description \code{format_valuef} generates functions that format
#'   values so that they print out with the desired number of digits.
#'
#' @param digits number of digits to use
#' @return Returns a function that takes an atomic vector as its argument
#'   and returns it formatted to character with \code{digits} decimals.
#' @keywords internal

format_valuef <- function(digits) {
  function(x) format(round(x, digits), nsmall=digits)
}



#' @describeIn GSMAR print method
#' @inheritParams plot.gsmar
#' @param digits number of digits to be printed (max 20)
#' @param summary_print if set to \code{TRUE} then the print will include approximate
#'  standard errors for the estimates, log-likelihood, information criteria values,
#'  modulus of the roots of the characteristic AR polynomials for each regime, and
#'  several unconditional moments.
#' @export

print.gsmar <- function(x, ..., digits=2, summary_print=FALSE) {
  gsmar <- x
  stopifnot(digits >= 0 & digits %% 1 == 0)

  # Helper functions
  format_value <- format_valuef(digits) # Function to format values for printing
  print_err <- function(val) { # Function for printing standard errors in brackets
    if(summary_print) cat(paste0(" (", format_value(val),")"))
  }
  make_string <- function(n_spaces, val) paste0(c(rep(" ", n_spaces), paste0("(", format_value(val), ")")), collapse="")
  add_string <- function(const_spaces, val1, val2) {
    if(summary_print) {
      err_string[length(err_string) + 1] <<- make_string(const_spaces + nchar(format_value(val1)) - nchar(format_value(val2)), val2)
    }
  }

  # Collect the relevant statistics etc
  p <- gsmar$model$p
  M <- gsmar$model$M
  params <- gsmar$params
  model <- gsmar$model$model
  restricted <- gsmar$model$restricted
  constraints <- gsmar$model$constraints
  all_mu <- get_regime_means(gsmar)
  npars <- length(params)
  nobs <- ifelse(is.null(gsmar$data), NA, length(gsmar$data))
  if(summary_print) all_vars <- get_regime_vars(gsmar)  # Unconditional regime variances

  if(gsmar$model$parametrization == "mean") { # Change to intercept parametrization
    params <- change_parametrization(p=p, M=M, params=params, model=model, restricted=restricted,
                                     constraints=constraints, change_to="intercept")
  }
  params <- remove_all_constraints(p=p, M=M, params=params, model=model, restricted=restricted, # Remove all constraints
                                   constraints=constraints)
  all_phi0 <- pick_phi0(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL)
  pars <- pick_pars(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M, params=params, model=model)

  # Obtain the standard errors if called from the summary method: pick them from the std error vector similarly to the parameter estimates
  if(summary_print) {
    all_ar_roots <- get_ar_roots(gsmar)
    std_errors <- remove_all_constraints(p=p, M=M, params=gsmar$std_errors, model=model, restricted=restricted,
                                         constraints=constraints) # These errors are valid ONLY IF there is NO multiplications or summations
    pars_err <- pick_pars(p=p, M=M, params=std_errors, model=model, restricted=FALSE, constraints=NULL)
    alphas_err <- pick_alphas(p=p, M=M, params=std_errors, model=model, restricted=FALSE, constraints=NULL)
    alphas_err[sum(M)] <- NA # The last mixing weight parameter is not parametrized
    dfs_err <- pick_dfs(p=p, M=M, params=std_errors, model=model)
    if(gsmar$model$parametrization == "mean") {
      mu_err <- pick_phi0(p=p, M=M, params=std_errors, model=model, restricted=FALSE, constraints=NULL)
      phi0_err <- rep(NA, sum(M))
    } else {
      mu_err <- rep(NA, sum(M))
      phi0_err <- pick_phi0(p=p, M=M, params=std_errors, model=model, restricted=FALSE, constraints=NULL)
    }
  }

  # Print the statistics related to the model
  cat("Model:\n", paste0(model, ", p = ", p, ", "))
  if(model == "G-StMAR") {
    cat(paste0("M1 = ", M[1], ", M2 = ", M[2], ", "))
  } else {
    cat(paste0("M = ", M, ", "))
  }
  cat(paste0("#parameters = " , npars, ","),
      ifelse(is.na(nobs), "\n", paste0("#observations = ", nobs, ",\n")),
      ifelse(gsmar$model$conditional, "conditional,", "exact,"),
      ifelse(gsmar$model$parametrization == "mean", "mean parametrization,", "intercept parametrization,"),
      ifelse(restricted, "AR parameters restricted,", "not restricted,"),
      ifelse(is.null(constraints), "no constraints.", "linear constraints imposed."), "\n")

  if(summary_print) {
    IC <- gsmar$IC
    form_val2 <- function(txt, val) paste(txt, format_value(val))
    cat("\n", paste(form_val2("log-likelihood:", gsmar$loglik),
                    form_val2("AIC:", IC$AIC),
                    form_val2("HQIC:", IC$HQIC),
                    form_val2("BIC:", IC$BIC),
                    sep=", "), "\n")
  }

  # Constraints are replicated for each regime so that there is no need to make special case for restricted models
  if(restricted & !is.null(constraints)) {
    constraints <- replicate(n=M, expr=constraints, simplify=FALSE)
  }

  # Print the regime statistics, parameter estimates, and standard errors.
  for(m in seq_len(sum(M))) { # Go through the regimes
    cat("\n")
    count <- 1
    err_string <- list()

    if(model == "GMAR") {
      regime_type <- "GMAR"
    } else if(model == "StMAR") {
      regime_type <- "StMAR"
      M1 <- 0
    } else {
      M1 <- M[1]
      regime_type <- ifelse(m <= M1, "GMAR", "StMAR")
    }
    cat(paste("Regime", m))
    if(model == "G-StMAR") cat(paste0(" (", regime_type, " type)"))

    if(summary_print) cat(paste("\nModuli of AR poly roots:", paste0(format_value(all_ar_roots[[m]]), collapse=", ")))

    cat(paste("\nMix weight:", format_value(alphas[m])))
    print_err(alphas_err[m])
    cat("\nReg mean:", format_value(all_mu[m]))
    print_err(mu_err[m])

    if(regime_type == "StMAR") {
      cat("\nVar param:", format_value(pars[nrow(pars), m]))
      print_err(pars_err[nrow(pars_err), m])
      cat("\n")

      cat("Df param:", format_value(dfs[m - M1]))
      print_err(dfs_err[m - M1])
    }
    if(summary_print) cat("\nReg var: ", format_value(all_vars[m])) # Unconditional regime variances
    cat("\n\n")

    cat(paste0("y = [", format_value(all_phi0[m]), "]"))
    add_string(const_spaces=4, all_phi0[m], phi0_err[m])

    for(i1 in seq_len(p)) { # Go through the lags
      cat(paste0(" + [", format_value(pars[1 + i1, m]), "]y.", i1))
      nspaces <- ifelse(i1 == 1, 3, 6)
      if(!is.null(constraints) && (any(constraints[[m]] != 1 & constraints[[m]] != 0) | any(rowSums(constraints[[m]]) > 1))) {
        # The constrained AR parameter standard errors multiplied open in 'pars_err' are valid only
        # if the constraint matrix (for the current regime) contains zeros and ones only, and
        # there is at most one one in each row (no multiplications or summations). The above if-sentece is
        # TRUE if that is the case -> we may print the standard errors
        add_string(const_spaces=nspaces, pars[1 + i1, m], NA)
        sep_AR <- TRUE
      } else {
        # Standard errors expanded from constraints are not valid in this regime.
        # -> AR parameter standard errors cannot be printed in brackets next to the AR coefficients as
        # the standard errors are related to the constraint parameters and not the AR coefficients.
        add_string(const_spaces=nspaces, pars[1 + i1, m], pars_err[1 + i1, m])
        sep_AR <- FALSE
      }
    }
    cat(" + ")

    if(regime_type == "GMAR") {
      cat(paste0("[sqrt(", format_value(pars[nrow(pars), m]), ")]eps"))
      add_string(const_spaces=11, pars[nrow(pars), m], pars_err[nrow(pars_err), m])
    } else { # StMAR or G-StMAR regime -> time-varying conditional error variance
      cat("[cond_sd]eps")
    }
    cat("\n")

    # Print the standard errors
    if(summary_print) cat(paste0(err_string, collapse=""), '\n')
    if(summary_print && sep_AR) cat(paste0("AR parameter std errors: ", paste0(format_value(pars_err[2:(p + 1), m]), collapse=", ")), "\n")
  }

  # Print unconditional moments
  if(summary_print) {
    um <- gsmar$uncond_moments
    cat("\nProcess mean:", format_value(um$uncond_mean), "\n")
    cat("Process var: ", format_value(um$uncond_var), "\n")
    cat("First p autocors:", format_value(um$autocors), "\n")
  }

  invisible(gsmar)
}


#' @title Print method from objects of class 'gsmarsum'
#'
#' @description \code{print.gsmarsum} is a print method for objects of class 'gsmarsum'
#'  created with the summary method \code{summary.gsmar}. Approximate standard errors
#'  are printed in brackets.
#'
#' @param x object of class 'gsmarsum' generated by \code{summary.gsmar}.
#' @param ... currently not in use.
#' @param digits the number of digits to be printed
#' @export

print.gsmarsum <- function(x, ..., digits) {
  if(missing(digits)) digits <- x$digits
  print.gsmar(x$gsmar, digits=digits, summary_print=TRUE)
}



#' @title Print method for class 'gsmarpred' objects
#'
#' @description \code{print.gsmarpred} is a print method for call 'gsmarpred'
#'  objects created with \code{predict.gsmar}.
#'
#' @inheritParams print.gsmarsum
#' @param x object of class \code{'gsmarpred'} generated by \code{predict.gsmar}.
#' @export

print.gsmarpred <- function(x, ..., digits=2) {
  gsmarpred <- x
  stopifnot(digits >= 0 & digits %% 1 == 0)
  format_value <- format_valuef(digits) # Function to formal values for printing

  # Different prints for different types of prediction intervals
  if(gsmarpred$pred_type == "cond_mean") {
    cat("One-step-ahead prediction by exact conditional mean, no prediction intervals.\n")
    cat("Forecast:", paste0(format_value(gsmarpred$pred), collapse=", "), "\n")

  } else if(gsmarpred$pi_type == "none") {
    cat(paste0("Prediction by ", gsmarpred$pred_type, ", no prediction intervals."), "\n")
    cat(paste0("Forecast ", gsmarpred$n_ahead, " steps ahead, based on ", gsmarpred$nsimu, " simulations.\n"))
    print(data.frame(pred=gsmarpred$pred))

  } else { # pi_type == two-sided, upper, or lower
    cat(paste0("Prediction by ", gsmarpred$pred_type, ", ", gsmarpred$pi_type,
               " prediction intervals with levels ", paste(gsmarpred$pi, collapse=", "), "."), "\n")
    cat(paste0("Forecast ", gsmarpred$n_ahead, " steps ahead, based on ", gsmarpred$nsimu, " simulations.\n"))
    cat("\n")

    # Create data frame that contains the point prediction and pi's, and then print it
    q <- gsmarpred$q
    pred_ints <- gsmarpred$pred_ints
    pred_type <- gsmarpred$pred_type
    df <- as.data.frame(lapply(1:ncol(pred_ints), function(i1) format_value(pred_ints[,i1])))
    names(df) <- q
    df[, pred_type] <- format_value(gsmarpred$pred)
    if(gsmarpred$pi_type == "two-sided") {
      new_order <- as.character(c(q[1:(length(q)/2)], pred_type, q[(length(q)/2 + 1):length(q)]))
    } else if(gsmarpred$pi_type == "upper") {
       new_order <- as.character(c(pred_type, q))
     } else {
      new_order <- names(df)
    }
     print(df[, new_order])
  }
  if(gsmarpred$pred_type != "cond_mean") {
    cat("\n Point forecasts and prediction intervals for mixing weights can be obtained with $mix_pred and $mix_pred_ints, respectively.\n")
  }

  invisible(gsmarpred)
}


#' @describeIn quantile_residual_tests Print method for class 'qrtest' objects
#' @param x object of class \code{'qrtest'} created with the function \code{quantile_residual_tests}.
#' @param ... graphical parameters passed to \code{segments} in \code{plot.qrtest}.
#'  Currently not used in \code{print.qrtest}
#' @param digits the number of digits to be print
#' @export

print.qrtest <- function(x, ..., digits=3) {
  qrtest <- x
  format_value <- format_valuef(digits) # Function to formal values for printing
  format_lag <- format_valuef(0)
  cat(paste("Normality test p-value:", format_value(qrtest$norm_res$pvalue)), "\n\n")

  # Function to print autocorrelation and cond. h.sked test results
  print_res <- function(which_res=c("ac", "ch")) {
    wres <- ifelse(which_res == "ac", 2, 3) # Index in the list 'qrtest'
    lags <- qrtest[[wres]]$lags
    for(i1 in seq_along(lags)) { # Go through tested lags
      sep <- ifelse(lags[i1] < 10, " | ", "| ")
      cat(" ", format_lag(lags[i1]), sep, format_value(qrtest[[wres]]$pvalue[i1]), "\n")
    }
  }

  # Print the statistics
  cat("Autocorrelation tests:\nlags | p-value\n")
  print_res("ac")

  cat("\nConditional hetetoskedasticity tests:\nlags | p-value\n")
  print_res("ch")
}
