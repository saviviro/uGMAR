# uGMAR 3.1.0

* Added a `NEWS.md` file to track changes to the package.
* New exported functions: get_regime_means, get_regime_autocovs, get_regime_vars, uncondMoments, get_soc, condMoments, stmarpars_to_gstmar, stmar_to_gstmar.
* Generally more functionality for conditional and unconditional moments, and convenient tools for switching to G-StMAR model from a StMAR model. 
* Implemented an algorithm by Monahan (1984) to the genetic algorithm for more thorough search of the parameter space near boundaries of the stationarity region.
* simulateGSMAR now provides better tools for forecasting. This update includes non-backward compatible changes for the return values if the argument ntimes is set to be larger than one. In additional to the samples, it now returns a list containing the mixing weights and component that was used to generate each observation.
* The arguments nCalls and nCores in fitGSMAR are now changed to ncalls and ncores for consistency.
* Fixes on minor bugs that used to cause errors in some special cases.      
* Updates on documentation

# uGMAR 3.2.0

* Added reference for the G-StMAR model.
* In the predict method arguments "ci" and "ci_type" were changed to "pi" and "pi_type" to signify "prediction interval"" as it's more correct expression than "confidence interval". Also the default prediction method is now median, and not mean.
* Changed the default number of CPU cores employed by the estimation function fitGSMAR to be at most two due to CRAN policy.
* Added the argument "seeds" to fitGSMAR allowing one to set the random number generator seed for each call to the genetic algorithm.
* Finite difference approximations for differentials regarding overly large degrees of freedom parameters now give reasonable approximation instead of numerical error.
* The maximum value for degrees of freedom parameters is now 1e5. 
* New exported function alt_gsmar that conveniently constructs a GSMAR model based on an arbitrary estimation round of fitGSMAR.
* New exported function get_foc which is the same as get_gradient but with convenient name.
* The default number of generations in the genetic algorithm is now 200 (was min(400, max(round(0.1*length(data)), 200)) before). 
* In various functions, user may now adjust the difference 'h' used in the finite difference approximations for differentials of the log-likelihood. 
* Bug fix: the summary print for gsmar objects falsely displayed standard error for the non-parametrized mixing weight
* Fixed typos etc. in documentation.

# uGMAR 3.2.1

* Fixed 'additional issue' in CRAN checks

# uGMAR 3.2.2

* New function: 'profile_logliks' for plotting profile log-likelihood functions.
* Disabled camelCase compatibility for arguments 'ncalls' and 'ncores' in fitGMAR.
* Updated the 'regime combining procedure' in the genetic algorithm to also support the G-StMAR model.
* Minor computation speed improvements.
* Tidier code for some parts.
* Improved comments and documentation.
* Bug fix: the function 'add_data' did not identify the model type correctly. 
* Bug fix: simulateGSMAR simulated some initial values from slighly wrong distribution; did not have affect on forecasts.
* Minor update on the summary print for the models

# uGMAR 3.2.3

* Updated the plot method for class 'gsmar' objects: now includes a density plot by default (can be removed).
* Updated the predict method for class 'gsmar' objects: now includes predictions for the mixing weights (can be removed from the plot).
* Fixed 'profile_logliks' to show correct headlines with mean parametrization + improved the default method for choosing the number of rows and colums in the plot-matrix.
* Now standard errors are printed correctly for models imposing all kinds of constraints. In the earlier versions, constrained AR parameter standard errors were printed incorrectly if the constraints involved multiplications or summations. 
* Removed redundant reinitialization of a PSOCK cluster in the function 'fitGSMAR'. 
* In the function quantileResidualTests the default argument for 'nsimu' is now 1 so that the tests are based on the given data only (and not on simulation).
* Added interest rate spead (10-Year minus 1-Year treasury) data.

# uGMAR 3.2.4

* Bug fix: the predict method incurred an error when plotting the results with n_ahead=1. 

# uGMAR 3.2.5

* New exported function: 'condmomentPlot' for further visualization of the model.
* Yet another bug fix: the predict method incurred an error with G-StMAR models. 
* Corrected degrees of freedom labels for G-StMAR models in the function profile_logliks.
* Updated the examples.

# uGMAR 3.2.6

* Fixed some typos in documentation
