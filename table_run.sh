#!/bin/bash
# This program iteratively submits qsub flux_stack_pbs.sh to FLUX.

echo -n "Are you sure you want to run the script table_run.sh? (y/n):"
read accept
if [ $accept != 'y' ];
	then echo 'Quitting...'
	exit
fi 

##################
## Begin Script ##
##################

## Initialize Configuration Arrays and Other Constants
cell_num=($(seq 1 49))				# Number of Cells
line_num=(2 5 10 15 25 50 100)			# Line of Sight Number 
gal_num=(5 10 15 25 50 100 150)			# Ngal number
halo_num=2100					# Number of Halos in Sample
method_num=0					# Ensemble Build Method
table_num=2					# Version of entire run table
data_loc="mass_mix/mm_0.25_run_table$table_num"	# Highest Directory for Data

## Go To Stacking Directory
cd /nfs/christoq_ls/nkern/Stacking

## Double Check if Directories Exist
# Parent Directory Check
if [ -d $data_loc ]
	then 
	echo "Directory: $data_loc , exists..."
	echo -n "Do you want to continue?(y/n):"
	read accept
	if [ $accept != 'y' ]
		then echo 'Quitting...'
		exit
	fi
	else 
	mkdir $data_loc 
	echo "Created Directory: $data_loc"
	echo -n "Continue?(y/n):"
	read accept
	if [ $accept != 'y' ]
		then echo 'Quitting...'
		exit
	fi
fi

# Sub Directory Check (including LOS)
for i in ${cell_num[*]}
do
	dir="mm_m0_run"$i
	if [ -d $data_loc\/$dir ]
		then
		echo -n
		else
		mkdir $data_loc/$dir	
	fi
done


## Begin Nested Loops Through Table
for i in $(seq 0 6)
do
	for j in $(seq 0 6)
	do
		let "k=($i*7)+$j"
		echo '----------------------------------------------------------'
		echo -e "cell_num=${cell_num[$k]}\tgal_num=${gal_num[$i]}\tline_num=${line_num[$j]}"
		# Submit Job Array To FLUX by feeding "mr_flux_stack_pbs.sh" job parameters
		# $1 : clus_num			(first caustic_mass_stack2D.py positional parameter)
		# $2 : gal_num			(second positional parameter)
		# $3 : line_num			(third positional parameter)
		# $4 : method_num		(fourth positional parameter)
		# $5 : cell_num			(fifth positional parameter)
		# $6 : job_array		(number of runs per job)
		# $7 : Submission Number	(appended to data_loc, ex: binstack_run_table3)
		# $8 : write_loc		(2nd level directory to write .pkl files in)

		# Define Constants
		job_array=(13 13 13 13 13 13 20)
		let "clus_num=$halo_num/(${job_array[$j]}+1)/${line_num[$j]}"
		write_loc="mm_m0_run"${cell_num[$k]}

		# Submit FLUX JOB for Ensembles
		_clus_num="$clus_num"
		_gal_num="${gal_num[$i]}"
		_line_num="${line_num[$j]}"
		_method_num="$method_num"
		_cell_num="${cell_num[$k]}"
		_job_array="${job_array[$j]}"
		_data_loc="$data_loc"
		_write_loc="$write_loc"
		
		sed -e "s:@@write_loc@@:$_write_loc:g;s:@@data_loc@@:$_data_loc:g;s:@@job_array@@:$_job_array:g;s:@@clus_num@@:$_clus_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$_method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@table_num@@:$table_num:g;s:@@run_los@@:0:g" < table_flux_stack_pbs.sh > $_data_loc/$_write_loc/script.sh

		# Change run_los = True if line_num == 100
		let "a=${cell_num[$k]}%7"
		if [ $a == 0 ]
			then
			sed -e "s:@@write_loc@@:$_write_loc:g;s:@@data_loc@@:$_data_loc:g;s:@@job_array@@:$_job_array:g;s:@@clus_num@@:$_clus_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$_method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@table_num@@:$table_num:g;s:@@run_los@@:1:g" < table_flux_stack_pbs.sh > $_data_loc/$_write_loc/script.sh

		fi

		# Submit Script to FLUX via qsub
		qsub $_data_loc/$_write_loc/script.sh 
		echo ""

		echo '----------------------------------------------------------'
	done
done


# End of Script
