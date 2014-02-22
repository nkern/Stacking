#!/bin/bash
# This program iteratively submits qsub flux_stack_pbs.sh to FLUX.

echo -n "Are you sure you want to run the script multiple_run.sh? (y/n):"
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
submission_num=3				# Version of entire run table
data_loc="binstack_run_table$submission_num"	# Highest Directory for Data

## Go To Stacking Directory
cd /nfs/christoq_ls/nkern/Stacking

## Double Check if Directories Exist
# Parent Directory Check
if [ -d $data_loc ]
	then 
	echo "Directory: $data_loc , exists..."
	echo -n "Do you want to continue?(y/s):"
	read accept
	if [ $accept != 'y' ]
		then echo 'Quitting...'
		exit
	fi
	else 
	mkdir $data_loc 
	echo "Created Directory: $data_loc"
fi
# Sub Directory Check (including LOS)
for i in $cell_num
do
	dir="bs_m0_run"$i
	if [ -d $data_loc/$dir ]
		then
		echo -n
		else
		mkdir $data_loc/$dir	
	fi
	# LOS Check
	let "a=$i%7"
	if [ $a == 0 ] 
		then
		if [ -d $data_loc/$dir"_los" ]
			then echo -n
			else
			mkdir $data_loc/$dir"_los"
		fi
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
		let "clus_num=$halo_num/${job_array[j]}/${line_num[j]}"
		write_loc="bs_m0_run"${cell_num[$k]}
		# Submit FLUX JOB for Ensembles
		echo qsub mr_flux_stack_pbs.py $clus_num ${gal_num[i]} ${line_num[j]} $method_num ${cell_num[$k]} ${job_array[j]} $write_loc
		echo ""
		# Submit FLUX JOB for LOS if line_num == 100
		let "a=${cell_num[k]}%7"
		if [ $a == 0 ]
			then
			write_loc="bs_m0_run"${cell_num[$k]}"_los"
			qsub mr_flux_stack_pbs.py $clus_num ${gal_num[i]} ${line_num[j]} $method_num ${cell_num[$k]} ${job_array[j]} $write_loc 'True'
			echo ""
		fi
		echo '----------------------------------------------------------'
	done
done


















