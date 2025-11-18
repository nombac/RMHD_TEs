#!/bin/sh
#BSUB -q Queue
#BSUB -n Prcs
#BSUB -e stderr.%J
#BSUB -o stdout.%J
#BSUB -J Targ
#BSUB -u shirose@jamstec.go.jp
cd ${HOME}/work/zeus_rad/data/Targ
mpijob -p \"rank%g:\ \" "./a.out.mstatx00" | grep rank0
