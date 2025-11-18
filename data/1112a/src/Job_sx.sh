#!/bin/sh
#PBS -l memsz_job=Memory
#PBS -q Queue
#PBS -o stdout.%s
#PBS -e stderr.%s
#PBS -v MPISEPSELECT=3
F_PROGINF=DETAIL; export F_PROGINF
F_RSVTASK=Prcs; export F_RSVTASK
F_FTRACE=Ftrace; export F_FTRACE
cd ${HOME}/zeus_rad/data/Targ
mpirun -np Prcs /usr/lib/mpi/mpisep.sh ./a.out.sx
