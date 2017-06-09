#!/bin/bash

input_dir=$1
output_dir=$2
tmp_dir=$3
num_cpus=$4

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ls $input_dir > $tmp_dir/clean_a3m_files
split -n l/$num_cpus $tmp_dir/clean_a3m_files $tmp_dir/clean_a3m_split

for i in `ls $tmp_dir/clean_a3m_split*`; do
	for j in `cat $i`; do 
		$DIR/clean_a3m.pl $input_dir/$j $output_dir/$j
	done &
done
