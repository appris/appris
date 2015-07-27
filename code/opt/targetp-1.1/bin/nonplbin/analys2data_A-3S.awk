# FILE analys2data_A-3S.awk
#
# IN:	files with averaged HOW network output ([mTP/SS]5to1HOW_scores.$$)
# OUT:	file with HOWLIN-input entries (only numbers...) AND
#	file with name, length and no of sequences
# NOTE:	if an input sequence is less than 100 aa long, the lacking stretch
# will be filled in with 0.000 scores.

# USAGE: nawk -f scripts/analys2data_A.awk -v w=100 -v name_file=results/name_file.$$
#		results/5to1HOW_scores.$$ > results/HOWLIN_in.$$

BEGIN {
	counter=0
	#print "# File: " name_file > name_file
}

/^ #/ {
	counter++
	name = $5 # originally $2, but name in that position truncated (why?)
	len = $3
	printf "%30s   %5.0f\n", name, len >> name_file
	i=0
	mTP_transit=0
	SS_transit=0
	no_transit=0
	while  (i<len) {
		getline
		i++
		if (i<=w) printf "%5.3f  %5.3f  ", $5, $10			
	}
	if (len<w) {
		for (j=1; j<=w-len; j++) {		  # om <100 fylls de på
			if (i+j<=w) printf "0.000  0.000  "
		}
	}	
	printf "%2.0f %2.0f %2.0f \n",mTP_transit,SS_transit,no_transit
	getline		# to get it right with the [mTP/SS]5to1HOW_scores.$$-file
	next
}

#END { 
#	printf "# Total no of entries: %7.0f\n", counter >> name_file
#}
