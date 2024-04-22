## SUBMIT BASH FILE FOR CLUSTER JOBS
## -----------------------------------------
## Impose to read the script with bash
## -----------------------------------------

#$ -S /bin/bash

## -----------------------------------------
## Replace NPROC by the desired number of processors
## -----------------------------------------

#$ -pe openmpi 20

## -----------------------------------------
## Name of the run
## -----------------------------------------

#$ -N flame_D_tutorial

##------------------------------------------
## Define the local directory, in case it is needed.
##------------------------------------------

MYDIR=`pwd`

##------------------------------------------
## export the variable MPIR_HOME which needs to be 
## defined in the .bash_profile or .tcshrc.
## -----------------------------------------

#$ -v MPIR_HOME

## -----------------------------------------
## Change directory to this job's working directory
## -----------------------------------------

#$ -cwd

## -----------------------------------------
## Echo  Start time and present Direcory
## -----------------------------------------

echo Start time is `date`
echo Directory is `pwd`

## -----------------------------------------
## Run the executable a.out
## -----------------------------------------
echo START OF THE OUTPUT FILE
#echo $SHELL
#source /etc/profile.d/common_var.sh
#which mpirun
#ulimit -l unlimited
#more  $TMPDIR/machines
#which mpirun 
mpirun -np $NSLOTS firePimpleSMOKE -parallel >log

## -----------------------------------------
## Echo end time
## -----------------------------------------

echo End time is `date`
