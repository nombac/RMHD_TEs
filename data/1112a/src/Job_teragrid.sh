#!/bin/bash
#PBS -q Queue
#PBS -m ae
#PBS -l walltime=Walltime
#PBS -M shirose@jamstec.go.jp
#PBS -l nodes=Nodes:ppn=Ppn:Resource
#PBS -A TG-MCA95C003
export F_UFMTENDIAN=3,4,22
cd ${HOME}/zeus_rad/data/Targ
mpirun -machinefile $PBS_NODEFILE -np Prcs ./a.out
