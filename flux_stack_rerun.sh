#!/bin/bash
# This program locates all flux_stack_pbs.sh runs that failed and re-runs them.

echo "Working Directory is: $PWD"
echo -n "Is this the desired directory? (y/n): "
read accept
if [ $accept != 'y' ];
	then echo 'Quitting...'
	exit
fi 


######################
#### Begin Script ####
######################

## Initialize Configuration Arrays and Other Constants
cell_num=($(seq 1 49))				# Number of Cells
line_num=(2 5 10 15 25 50 100)			# Line of Sight Number 
gal_num=(5 10 15 25 50 100 150)			# Ngal number
clus_num=(75 30 15 10 6 3 1)			# Number of Ens Clusters done per instance
job_num=(14 14 14 14 14 14 20)			# Number of Jobs Submitted
halo_num=2100					# Number of Halos in Sample
method_num=0					# Ensemble Build Method
table_num=3					# Version of entire run table
data_loc="binstack_run_table$table_num"		# Highest Directory for Data
write_loc="bs_m0_run"				# Stem of write_loc directory

## Find failed jobs ##
ENS_DIRS=()
ENS_NUMS=()
LOS_DIRS=()
LOS_NUMS=()
# i loops over Ngal, j loops over LOS
for i in {0..6}
do
	for j in {0..6}
	do
		let "k=$i*7+($j+1)"
		dir=$write_loc$k
		echo "Working on Directory: $dir..."
		m=0
		for n in $(seq 1 $(($halo_num/${line_num[$j]}/${job_num[$j]})) $(($halo_num/${line_num[$j]})))
		do	
			n=$((n-1))
			if [ -a $dir/Ensemble_$n\_Data.pkl ]
				then
				echo -n
			else
				ENS_DIRS+=($k)
				ENS_NUMS+=($m)	
			fi
			# Check LOS systems too
			let "a=$k%7"
			if [ $a == 0 ]
				then
				echo $k
				if [ -a $dir\_los/Ensemble_$n\_Data.pkl ]
					then
					echo -n
				else
					LOS_DIRS+=($k)
					LOS_NUMS+=($m)
				fi
			fi
			m=$((m+1))
		done
	done
done
echo "ENS_DIRS=${ENS_DIRS[*]}"
echo "ENS_NUMS=${ENS_NUMS[*]}"
echo "LOS_DIRS=${LOS_DIRS[*]}"
echo "LOS_NUMS=${LOS_NUMS[*]}"
echo -n "Does this seem right? (y/n):"
read accept
if [ $accept != 'y' ];
	then echo 'Quitting...'
	exit
fi 
echo "Beginning FLUX job submission..."






