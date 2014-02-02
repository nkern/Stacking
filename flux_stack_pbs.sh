#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=2:00:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/Stacking/binstack_run_table/bs_m1_run7/
#PBS -j oe
#PBS -t 0-6
#

cd /nfs/christoq_ls/nkern/Stacking
python caustic_mass_stack2D.py $PBS_ARRAYID 3 5 100 1 7






