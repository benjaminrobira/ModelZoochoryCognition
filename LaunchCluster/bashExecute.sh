#### ------------------------------------------ ##
#### SSH Script to run ROAD analysis on cluster ##
#### ------------------------------------------ ##

#$-cwd #Execute job in current directory

#$-M benjamin.robira@normalesup.org
#$-m abe

### say where to put output
#$-o /nfs1/$USER/
#$-e /nfs1/$USER/

# select a queue
#$-q global.q

### ask for cores
#$-pe smp 1

### set resource limits
###-l s_vmem=20G
###-l h_vmem=15G
###-l virtual_free=15G
#$-l mem_free=12G
###-l mf=20G

###Test script to see how path worked
###cd /nfs1/robirab/
rm /nfs1/robirab/ModelZoochoryCognition/workedR7a.txt
###echo $HOME $USER > /nfs1/$USER/workedR.txt
###singularity3 --version > /nfs1/robirab/workedR.txt
###cd /nfs1/robirab/
###echo $(pwd) > /nfs1/robirab/test.txt

###Launch the R script
singularity3 exec -B /nfs1:/nfs1 /nfs1/$USER/ModelZoochoryCognition/R.sif R CMD BATCH /nfs1/robirab/ModelZoochoryCognition/Scripts/R/EXECUTE_CLUSTER.R /nfs1/robirab/ModelZoochoryCognition/Rscript7a.Rout &> /nfs1/robirab/ModelZoochoryCognition/workedR7a.txt

