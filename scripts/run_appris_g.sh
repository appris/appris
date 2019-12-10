#/usr/bin/env bash

set -e


anno_dir="$1";
dir_list_file="$2";

find "$anno_dir" -type d -name "ENSG*" > "$dir_list_file"

o_dir=$(pwd);
while read gene_dir
do
  cd "$gene_dir";
  if [[ -f "appris" ]];
  then
    echo "Running singĺe-gene APPRIS in: ${gene_dir}"
	# TODO: use all APPRIS modules
    appris_bin_g -s "Homo sapiens" -m fm1scra -l info > log 2>&1;
  fi
  cd "$o_dir";
done < "$dir_list_file"

