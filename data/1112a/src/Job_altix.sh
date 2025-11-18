#!/bin/sh
#BSUB -q Queue
#BSUB -e stderr.%J
#BSUB -o stdout.%J
#BSUB -J Targ
ulimit -c 0
cd ${HOME}/zeus_rad/data/Targ
#mpirun -np Prcs ./a.out
queue=Queue
case $queue in
P1) mpirun -np 4  /usr/bin/dplace -c40-43 -s1 ./a.out.altix ;;
P2) mpirun -np 4  /usr/bin/dplace -c40-43 -s1 ./a.out.altix ;;
P3) mpirun -np 4  /usr/bin/dplace -c44-47 -s1 ./a.out.altix ;;
P4) mpirun -np 16 /usr/bin/dplace -c48-63 -s1 ./a.out.altix ;;
esac
