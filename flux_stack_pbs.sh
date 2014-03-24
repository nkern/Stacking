#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=2:00:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/Stacking/binstack_run_table4/bs_m0_run5/
#PBS -j oe
#

cd /nfs/christoq_ls/nkern/Stacking
python caustic_mass_stack2D.py 6 6 5 25 0 5 4 0




