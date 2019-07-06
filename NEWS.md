# uGMAR 3.1.0

* Added a `NEWS.md` file to track changes to the package.
* New exported functions: get_regime_means, get_regime_autocovs, get_regime_vars, uncondMoments, get_soc, condMoments, stmarpars_to_gstmar, stmar_to_gstmar.
* Generally more functionality for conditional and unconditional moments, and convenient tools for switching to G-StMAR model from a StMAR model. 
* Implemented an algorithm by Monahan (1984) to the genetic algorithm for more thorough search of the parameter space near boundaries of the stationarity region.
* simulateGSMAR now provides better tools for forecasting. This update includes non-backward compatible changes for the return values if the argument ntimes is set to be larger than one. In additional to the samples, it now returns a list containing the mixing weights and component that was used to generate each observation.
* The arguments nCalls and nCores in fitGSMAR are now changed to ncalls and ncores for consistency.
* Fixes on minor bugs that used to cause errors in some special cases.      
* Updates on documentation
* Added inflation expectation data (IE)

# uGMAR 3.1.1

* In the predict method arguments "ci" and "ci_type" were changed to "pi" and "pi_type" to signify "prediction interval"" as it's more correct expression than "confidence interval". Also the default prediction method is now median, and not mean.
* Fixed typos etc. in documentation.
* Changed the default number of CPU cores employed by the estimation function fitGSMAR to be at most two due to CRAN policy.
