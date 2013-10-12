#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N SELF-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=50:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/Stacking/flux_output/
#PBS -j oe
#PBS -t 0-99
#

cd /nfs/christoq_ls/nkern/Stacking
python caustic_mass_stack2D.py $PBS_ARRAYID 10 10 1 1 





