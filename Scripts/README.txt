What do these folders contain?
This README file explains how the scripts are organised.

You have two folders: one containing files related to the Rcpp language, the other to the R language. Each contains the script in Rcpp or R.
It is important to note that these are WORKING scripts, intended to facilitate verification, which means that their cleanup is not complete. Therefore, the scripts may contain code not used in the present analyses (e.g. some spatial indices tested but not finally used).

- R folder:
  -> A: a script testing the Rcpp functions
  -> B: a script testing the spatial point distribution quantification indices
  -> 0: the script listing all parameters (and their values)
  -> 1 to 5: the scripts that run the simulation for different spatio-temporal knowledge
  -> 6: the script that processes the output of simulations 1 to 5
  -> 7a, b, c: The scripts that run the simulation for other questions (change of movement rule, distance competition, ...)
  -> 8: the script that processes the output of the simulations of 7
  
- Rcpp folder:
 -> FunctionsRcpp.cpp : This script contains all the Rcpp code that executes the simulations (independent functions and then the wrap up).
