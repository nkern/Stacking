#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=2:00:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/Stacking/@@data_loc@@/@@write_loc@@/
#PBS -j oe
#PBS -t 0-@@job_array@@
#

cd /nfs/christoq_ls/nkern/Stacking

python caustic_mass_stack2D.py $PBS_ARRAYID @@clus_num@@ @@gal_num@@ @@line_num@@ @@method_num@@ @@cell_num@@ @@table_num@@ @@run_los@@


