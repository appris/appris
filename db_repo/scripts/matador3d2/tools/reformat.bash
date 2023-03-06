#!/bin/bash

input_dir=$1
output_dir=$2
tmp_dir=$3
num_cpus=$4
reformat_bin=$5


ls $input_dir > $tmp_dir/reformat_in_files
split -n l/$num_cpus $tmp_dir/reformat_in_files $tmp_dir/reformat_out_files

for i in `ls $tmp_dir/reformat_out_files*`; do
	for j in `cat $i`; do
		file_out=$output_dir/${j:0:-4}.a2m
		if [ ! -e "$file_out" ]; then
			$reformat_bin a3m a2m $input_dir/$j $file_out -M first -r -l 32765
		fi
	done &
done


