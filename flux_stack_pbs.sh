#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N SELF-STACK
#PBS -l nodes=1:ppn=1,pmem=3000mb,walltime=2:30:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/Stacking/selfstack_run_table/ss_m1_run48/
#PBS -j oe
#PBS -t 0-58
#

cd /nfs/christoq_ls/nkern/Stacking
python caustic_mass_stack2D.py $PBS_ARRAYID 36 150 50 48






