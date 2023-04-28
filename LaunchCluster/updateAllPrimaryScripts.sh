#This script modifies all necessary documents to ensure correct run on the cluster

# Connect to cluster
ssh robirab@sge-qsub.intra.ismaa.it
cd nfs1/robirab/ModelZoochoryCognition

#Submit work
qsub /nfs1/robirab/ModelZoochoryCognition/LaunchCluster/bashExecute.sh
qstat -u robirab

#FROM PC TO CLUSTER

##Update the whole project
scp -r ~/Documents/Miscellaneous/ModelZoochoryCognition/ robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/
rsync -au ~/Documents/Miscellaneous/ModelZoochoryCognition/ robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/

##Update the R scripts
scp -r ~/Documents/Miscellaneous/ModelZoochoryCognition/Scripts robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition

##Update the bash script
scp ~/Documents/Miscellaneous/ModelZoochoryCognition/LaunchCluster/bashExecute.sh robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/LaunchCluster

##Update the image
scp ~/Documents/Miscellaneous/ModelZoochoryCognition/R.sif robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/

# FROM CLUSTER TO PC

# Update output SPACE TREE
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivitySpaceTree ~/Documents/Miscellaneous/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivitySpaceTree ~/Documents/Miscellaneous/ModelZoochoryCognition/Output

#Update output MOVING RULE
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivityMovingRule ~/Documents/Miscellaneous/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/SensitivityMovingRule ~/Documents/Miscellaneous/ModelZoochoryCognition/Output

#Update output MAIN
scp -r robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/Main ~/Documents/Miscellaneous/ModelZoochoryCognition/Output
rsync -au robirab@sge-qsub.intra.ismaa.it:/nfs1/robirab/ModelZoochoryCognition/Output/Main ~/Documents/Miscellaneous/ModelZoochoryCognition/Output




