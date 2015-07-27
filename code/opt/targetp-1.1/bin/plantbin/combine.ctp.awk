#! /usr/bin/gawk -f
# FILE combine.oe.ctp.awk
# 
# 
#
# Usage: gawk -f combine.oe.ctp.awk > file-to-put-scores-in
# Outputs: (FILE file-to-put-scores-in)


#BEGIN {
	#if (format=="short") {
	#	printf "# SignalP %s predictions\n", group
	#	print "# name       Cmax  pos ?  "\
	#		"Ymax  pos ?  Smax  pos ?  Smean ?"
	#}
	#else {
	#printf "************************* CTP predictions "
	#printf "*************************\n"
	#printf "Using networks trained on %s data\n", group

	#printf "\nC = raw cleavage site score\n"
	#printf "S = raw signal peptide score\n"
	#printf "Y = combined cleavage site score:  Y = sqrt(C*(-delta S)),\n"
	#printf "    where delta S is averaged over a window " \
	#	"of 2*%d positions\n\n", drange
	#}
#}

#NR=1 { 	
	# In gawk (a.o.) command line variable assignments has not taken
	# effect in the BEGIN clause
	#Esplit(Cscales, Cwt, ",")
#}

/^ #/ { 
	N++
	#len = substr($0,30,6)+0 # KRap: in how98, names are 20 long!
	len = substr($0,39,6)+0
	name = substr($0,7,11)
	seq = ""
}

/SINGLE/ { 
	out=1 
	next
}

(out) {
        # print len # debug printout
	S[$1] = 0
	AA[$1]= $2
	if ($3== "P" || $3=="_")  TP_real[$1]=$3
	else TP_real[$1]="?"
	for (i=0; i<5; i++) {
		S[$1] += $(6*i+5)
	}
	S[$1] /= 5
	if (S[$1]>1) S[$1]=1
	#seq = seq $2
}

out && $1 == len { 
	out=0
	printf " # %-10s %5d \n",
		name, len 
	for (i=1; i<=len; i++){
	        if (S[i]>= 0.5) TP_assign = "P"
		else TP_assign = "_"
	        printf "%5d %-1s %-1s %-1s %5.3f \n", i, AA[i], TP_real[i], TP_assign, S[i] 
		# printing of TP_real lacks of course meaning when using unknown sequences
	}
	printf "\n" 
}
