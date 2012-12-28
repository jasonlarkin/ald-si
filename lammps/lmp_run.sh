#!/bin/sh
#PBS -l nodes=1:ppn=1:node
### Merge stderr with stdout
#PBS -j oe
### Queue name
#PBS -q default
###Job name
#PBS -N ankit_lammps
### Declare job-non-rerunable
#PBS -r n
#PBS -V
# This job's working directory
echo Job ID: $PBS_JOBID
echo Working directory is $PBS_O_WORKDIR cd $PBS_O_WORKDIR echo Running on host `hostname` echo Time is `date` echo Directory is `pwd` echo This job runs on the following processors:
echo `cat $PBS_NODEFILE`


RUNPATH=/home/jason/ankit/lammps/
#EXEPATH=/home/jason/lammps/lammps-30Nov10/src
EXEPATH=/home/jason/lammps/lammps-23Apr12/src

cd $RUNPATH

#bw-01
#/opt/intel/impi/3.2.2.006/bin64/mpirun -np 4 $EXEPATH/lmp_openmpi < $RUNPATH/ligand.in
#bw-02
/opt/open-mpi/tcp-gnu41/bin/mpirun -np 1 $EXEPATH/lmp_openmpi < $RUNPATH/input.dat

#python 2_replaceXYZ.py
