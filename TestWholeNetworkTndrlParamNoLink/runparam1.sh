#!/bin/csh

#$ -M 	@ucr.edu # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  gpu 	 # Specify queue
#$ -l gpu_card=1
#$ -pe smp 4         #specifies threads??? maybe
#$ -N  TWNTP	 # Specify job name
#$ -t 1       #specify number of data input files

set data = ( data.xml )


module purge
module load gcc/6.2.0
module load gsl/2.3
module load cuda/9.1
echo -n "It is currently: ";date
echo -n "I am logged on as ";whoami
echo -n "This computer is called ";hostname
echo -n "I am currently in the directory ";pwd


./bend-model -eps=0.0001 -dt=0.001 $data[${SGE_TASK_ID}]
