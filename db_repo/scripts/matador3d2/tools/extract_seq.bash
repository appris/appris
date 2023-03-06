#!/bin/bash

input_file=$1
output_file=$2

head -n 2 $input_file | awk '{if($1~/^>/){print $1"_pdb_seq"}else{print}}' > $output_file
