# Get chromome sizes for trackHUBS (deleting comment line)

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -s -e "select chrom, size from hg19.chromInfo"  > hg19.chrom.sizes


