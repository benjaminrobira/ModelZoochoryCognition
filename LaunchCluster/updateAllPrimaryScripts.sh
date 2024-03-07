#This script modifies all necessary documents to ensure correct run on the cluster

# Connect to cluster
ssh robirab@sge-qsub.intra.ismaa.it
cd /nfs1/robirab/ModelZoochoryCognition

#Submit work
qsub /nfs1/robirab/ModelZoochoryCognition/LaunchCluster/bashExecute.sh
qstat -u robirab

#FROM PC TO CLUSTER

##Update the whole project
scp -r /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/ robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/
rsync -au /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/ robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/

##Update the R/Rcpp scripts
scp -r /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Scripts robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition

##Update the bash script
scp /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/LaunchCluster/bashExecute.sh robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/LaunchCluster

##Update the image
scp /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/R.sif robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/

# FROM CLUSTER TO PC

# Update all outputs
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition

# Update output SPACE TREE
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivitySpaceTree /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivitySpaceTree /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output

#Update output MOVING RULE
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivityMovingRule /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivityMovingRule /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output

#Update output MAIN
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/Main /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/Main /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output

#Update output SPEED
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivitySpeed /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivitySpeed /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output

#Update output LEARNING RULES
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivityLearning /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivityLearning /home/robirab/Documents//Miscellaneous/Research_Projects/ModelZoochoryCognition/Output


