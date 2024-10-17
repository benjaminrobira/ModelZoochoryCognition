# How spatiotemporal cognition and movement of seed-dispersing animals influence plant distribution

This repository contains the codes and results to reproduce the paper **"How spatiotemporal cognition and movement of seed-dispersing animals influence plant distribution"** by Benjamin Robira, now published in Oikos.

It is divided into several folders:

## [:file_folder: **Raw data**](Text)

This folder contains the files used to write the article (in *Rmarkdown, Rmd*). You will find
- The "lua" folders (`author-info-blocks`, `lua-filters`, `scholarly-metadata`), including the files needed to calculate the layouts of the output pdf.
- The *Rmd* file (`Article_ModelSeedDispersal.Rmd`, `Article_ModelSeedDispersal_Revised.Rmd`) encoding the article (initial version, revised version), with a complementary file (`layout.Rmd`) used to set up the layout of the paper.
- The initial *bib* library file (`bibliography.bib`), which contains the references (bib format) used in the paper (references to *R* packages are added after, during the Rmd compilation).
- The *pdf* file (`Article_ModelSeedDispersal.pdf`), which is the article (output from the *Rmd* file).

## [:file_folder: **Scripts**](Scripts)

This folder contains files that encode the simulation (Rcpp folder) or its execution and analysis (R folder) (this is explained in the `Readme.txt` file). It is important to note that these are WORKING scripts, intended to facilitate verification, i.e. their cleanup is not complete. Therefore, the scripts may contain code not used in the present analyses (e.g. some spatial indices tested but not finally used).

* [:file_folder: **Rcpp**](Scripts/Rcpp) The *C++* code (in *R*, so based on *Rcpp*) coding for the simulation:
  `FunctionsRcpp.cpp` : This script contains all the Rcpp code that executes the simulations (independent functions and then the wrap up).

* [:file_folder: **R**](Scripts/R) The *R* to actually run the simulation, extract the results, etc.
  - A: the script testing the Rcpp functions
  - B: the script testing the spatial point distribution quantification indices (to set up benchmarks)
  - C: the script testing the quantification of autocorrelation in fruiting date (to set up benchmarks)
  - 0: the script listing all parameters (and their values)
  - 1 to 5: the scripts that run the simulation for different spatio-temporal knowledge
  - 6: the script that processes the output of simulations 1 to 5
  - 7a, b, c, d, e: the scripts that run the simulation for other questions (change of movement rule, distance competition, ...)
  - 8: the script that processes the output of the simulations of 7
  - 9: the script to access the variance in fruiting dates over the course of simulations
  - 10: the script that runs the sensitivity test on sensory range 
  - toolbox: my private toolbox which contains various useful functions for data processing. Any reuse should include a request to me at benjamin.robira@normalesup.org (and it will surely be answered positively!). I just want to know how it is distributed. I also update it regularly, so this may become an outdated version).
  - the script that creates the graphic legend (i.e. the black circles) on the Y-axis of some plots
  - RiotteLambert_2017_f: Functions adapted (or simply copied) from Riotte-Lambert, Louise, Simon Benhamou, and Simon Chamaillé-Jammes. "From randomness to traplining: a framework for the study of routine movement behavior." Behavioral Ecology (2016): arw154. Their use is acknowledged in the manuscript by referring to the paper. I also contacted the first author for this.

  
## [:file_folder: **LaunchCluster**](Scripts/LaunchCluster) Launch Cluster

This folder contains the files used to run the simulations on the Edmund Mach Foundation HPC cluster. They are the following
 - The bash commands to start the simulations on the cluster (`bashExecute.sh`). Note that this needs to be adapted for each *R* script to be run.
 - The command lines I used to create the singularity image on my own machine (`CreateSingularityImage.txt`).
 - The bash commands to update all the necessary scripts/files/folders on the cluster or on my computer (e.g. push the scripts to the cluster, import the data from the cluster, etc.). `updateAllPrimaryScripts.sh`).

To run the codes on the HPC cluster of the Edmund Mach Foundation, I used a singularity image (a kind of minicomputer that can be freely used to have exactly the same settings as when I created and used it myself). However, this image is too heavy to be stored on git. The file `R.def` shows you the settings to create this image. If you need it, just drop me an email!

:question: Please, if you encounter any problems or see any errors in the code, have any questions or would like to discuss this topic, feel free to email me at [:e-mail:](mailto:benjamin.robira@normalesup.org) benjamin.robira@normalesup.org  

I hope you find it to be of quality (and without errors)!   
Benjamin.  