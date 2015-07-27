#! /usr/bin/nawk -f

BEGIN {
	nparts = 5	# No. of partitions to average over
	nf = 6		# No. of fields in each partition file

        if (format=="short") {
                printf "# SignalP-NN %-50s\n", type " predictions"
                print "# name                Cmax  pos ?  "\
                        "Ymax  pos ?  Smax  pos ?  Smean ?  D     ? " #jannick 11-08-03
        }

}	

NR=1 { 	
	# The values in $TYPE.param have been read from the command line
	# In gawk (a.o.) command line variable assignments has not taken
	# effect in the BEGIN clause
	split(Cscales, Cwt, ",")
}

/^ #/ { 
	N++
	len = substr($0,39,6)+0 
#	name = substr($0,7,30)
	name = substr($0,7,20)  #orig
	seq = ""
}

/SINGLE/ { 
	out=1 
	next
}

(out) {
	C[$1] = S[$1] = 0
	for (i=0; i<nparts; i++) {
		C[$1] += $(nf*i+5) * Cwt[i+1]
		S[$1] += $(nf*i+nf*nparts+5)
	}
	C[$1] = C[$1]/nparts
	if (C[$1]>1) C[$1]=1
	S[$1] = S[$1]/nparts
	seq = seq $2
}

out && $1 == len { 
	out=0
	datafile = datafileprefix "." N
	if (format != "short") {
		printf ">%s  length = %d\n", name, len 
		printf ">%s  length = %d\n\n", name, len > datafile
		print "# pos  aa    C       S       Y" > datafile
	}	
	Cmax = Smax = Ymax = 0
	Cmaxpos = Smaxpos = Ymaxpos = 1
	for (i=1-drange; i<=0; i++)
		S[i]=S[1]
	for (i=len+1; i<=len+drange; i++)
		S[i]=S[len]
	for (i=1; i<=len; i++) {
		diff = 0
		for (j= i-drange; j<i; j++)
			diff += S[j]
		for (j= i; j< i+drange; j++)
			diff -= S[j]
		diff = diff/drange
		Y = (diff>0) ? sqrt(C[i]*diff) : 0
		if (Y>Ymax) { Ymax=Y; Ymaxpos=i }
		if (S[i]>Smax) { Smax=S[i]; Smaxpos=i }
		if (C[i]>Cmax) { Cmax=C[i]; Cmaxpos=i }
				if (format != "short") {
			printf "%5d   %1s   %5.3f   %5.3f   %5.3f\n", 
				i, substr(seq,i,1), C[i], S[i], Y > datafile
		}
	}
	close (datafile)
	Smean = 0
	for (i=start; i<=Ymaxpos-1; i++)
		Smean += S[i]
	if (Ymaxpos>1)
		Smean = Smean/(Ymaxpos - 1)
	else
		Smean = 0
        Dmax = (Ymax/2 + Smean/2)				#Jannick 11-08-03
	 
	if (format=="short") {
		printf "%-10s  ", name
		printf "%5.3f %3d %1s  ", 
			Cmax, Cmaxpos, (Cmax>Cmaxcut) ? "Y" : "N"
		printf "%5.3f %3d %1s  ", 
			Ymax, Ymaxpos, (Ymax>Ymaxcut) ? "Y" : "N"
		printf "%5.3f %3d %1s  ", 
			Smax, Smaxpos, (Smax>Smaxcut) ? "Y" : "N"
		printf "%5.3f %1s", 
			Smean, (Smean>Smeancut) ? "Y" : "N"
		printf "%7.3f %1s\n", 				#Jannick 11-08-03
			Dmax, (Dmax>Dmaxcut) ? "Y" : "N"		#Jannick 11-08-03

		next
	}
	printf "# Measure  Position  Value  Cutoff  signal peptide?\n"
	printf "  max. C   %3d       %5.3f  %5.2f   %s\n",
		Cmaxpos, Cmax, Cmaxcut, (Cmax>Cmaxcut) ? "YES" : "NO"
	printf "  max. Y   %3d       %5.3f  %5.2f   %s\n",
		Ymaxpos, Ymax, Ymaxcut, (Ymax>Ymaxcut) ? "YES" : "NO"
	printf "  max. S   %3d       %5.3f  %5.2f   %s\n",
		Smaxpos, Smax, Smaxcut, (Smax>Smaxcut) ? "YES" : "NO"
	printf "  mean S     1-%-4d  %5.3f  %5.2f   %s\n",
		Ymaxpos-1, Smean, Smeancut, (Smean>Smeancut) ? "YES" : "NO"
	printf "       D     1-%-4d  %5.3f  %5.2f   %s\n",
		Ymaxpos-1, Dmax, Dmaxcut, (Dmax>Dmaxcut) ? "YES" : "NO"
	if (Smax>Smaxcut || Ymax>Ymaxcut || Smean>Smeancut) {
	printf "# Most likely cleavage site between pos. %d and %d: %s\n",
		Ymaxpos-1, Ymaxpos,
		substr(seq,Ymaxpos-3,3) "-" substr(seq,Ymaxpos,2)
	}
	print ""
}
