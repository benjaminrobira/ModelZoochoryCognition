What does these folders contain?
This README files will explain how the scripts are organised.

You have two folders: the one containing files related to Rcpp language, the other one to R language. Each contains the script in Rcpp or R.
It is important to note that those are WORKING scripts, intended to facilitate review, which means that their cleaning is not finished. Hence, the scripts might included codes unused in the present analyses (e.g. some spatial indices tested but not eventually used).

- R folder:
  -> A: A script testing the Rcpp functions
  -> B: A script testing the indices to quantify spatial point distributions
  -> 0: the script listing all parameters (and their value)
  -> 1 to 5: The scripts launching the simulation for different spatio-temporal knowledge
  -> 6: The script processing the output of the simulations

- Rcpp folder:
 -> FunctionsRcpp.cpp : This script contains all the Rcpp code running the simulations (independent functions and then the wrap up).