#!/bin/sh
#PBS -l memsz_job=Memory
#PBS -l cpunum_job=Prcs
#PBS -q Queue
#PBS -v MPISEPSELECT=3
F_PROGINF=DETAIL; export F_PROGINF
F_FTRACE=Ftrace; export F_FTRACE
cd ${HOME}/work/zeus_rad/data/Targ
mpirun -np Prcs /usr/lib/mpi/mpisep.sh ./a.out.mstsx8rf
