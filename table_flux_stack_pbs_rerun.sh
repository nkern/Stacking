#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=2:00:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/Stacking/mass_mix/mm_0.50_run_table@@table_num@@/@@write_loc@@/
#PBS -j oe
#

cd /nfs/christoq_ls/nkern/Stacking

python caustic_mass_stack2D.py @@run_num@@ @@clus_num@@ @@gal_num@@ @@line_num@@ @@method_num@@ @@cell_num@@ @@table_num@@ @@run_los@@


