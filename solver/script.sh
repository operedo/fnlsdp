#!/bin/bash
#$ -cwd
#$ -j y
#$ -N out_pcsdp
#$ -pe mpi 8
#$ -S /bin/bash
## Nodos (queues) donde funciona ok la libreria Intel MKL
##$ -q all.q@compute-1-1
##$ -q all.q@compute-1-4
##$ -q all.q@compute-1-5 
##$ -q all.q@compute-1-6
##$ -q all.q@compute-1-7
##$ -q all.q@compute-1-8 
##$ -q all.q@compute-1-9 
##$ -q all.q@compute-1-10
##$ -q all.q@compute-1-11
##$ -q all.q@compute-1-13


mpirun -np $NSLOTS ./fnlsdp
