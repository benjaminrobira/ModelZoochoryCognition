What do these folders contain?
This README file explains how the scripts are organised.

You have two folders: one containing files related to the Rcpp language, the other to the R language. Each contains the script in Rcpp or R.
It is important to note that these are WORKING scripts, intended to facilitate verification, which means that their cleanup is not complete. Therefore, the scripts may contain code not used in the present analyses (e.g. some spatial indices tested but not finally used). If you find an error, that would be a shame, but I would be glad to be aware. Please contact me at benjamin.robira@normalesup.org

- R folder:
  -> A: the script testing the Rcpp functions
  -> B: the script testing the spatial point distribution quantification indices (to set up benchmarks)
  -> C: the script testing the quantification of autocorrelation in fruiting date (to set up benchmarks)
  -> 0: the script listing all parameters (and their values)
  -> 1 to 5: the scripts that run the simulation for different spatio-temporal knowledge
  -> 6: the script that processes the output of simulations 1 to 5
  -> 7a, b, c, d, e: the scripts that run the simulation for other questions (change of movement rule, distance competition, ...)
  -> 8: the script that processes the output of the simulations of 7
  -> 9: the script to access the variance in fruiting dates over the course of simulations
  -> 10: the script that runs the sensitivity test on sensory range 
  -> toolbox: my private toolbox which contains various useful functions for data processing. Any reuse should include a request to me at benjamin.robira@normalesup.org (and it will surely be answered positively!). I just want to know how it is distributed. I also update it regularly, so this may become an outdated version).
  -> the script that creates the graphic legend (i.e. the black circles) on the Y-axis of some plots
  -> RiotteLambert_2017_f: Functions adapted (or simply copied) from Riotte-Lambert, Louise, Simon Benhamou, and Simon Chamaillé-Jammes. "From randomness to traplining: a framework for the study of routine movement behavior." Behavioral Ecology (2016): arw154. Their use is acknowledged in the manuscript by referring to the paper. I also contacted the first author for this.


- Rcpp folder:
 -> FunctionsRcpp.cpp : This script contains all the Rcpp code that executes the simulations (independent functions and then the wrap up).
