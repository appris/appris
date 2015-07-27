#! /usr/bin/nawk -f
# FILE determine-CS.awk
# used by mTP_CS-1.01 to perform the choice between the output from the 3 weight matrices
# supply -v max_len=#
# NOTE: this file is CHANGED compared to version 1.0. It now presents the PRED mTP length
# and the SCORE. It is not intended for testing known sequences.

BEGIN {
Nseq=0
highest=0	# value
CS_highest=0	# position
snd_highest=0	# v
CS_snd_highest=0 #p
trd_highest=0   # v
CS_trd_highest=0 #p
maxlen=max_len	# maximum length of predicted mTP; max_len is from cmd line
}


{ 
if ($1!="#") {
Nseq++
name[Nseq]=$1		
len[Nseq]=$2
sc2[Nseq]=$3		# CS score
sc3[Nseq]=$8		#   -"-
sc10[Nseq]=$13		#   -"-
#true[Nseq]=$4		# true CS position : Not used in this version
pred2[Nseq]=$5		# predicted CS position
pred3[Nseq]=$10		#        -"-
pred10[Nseq]=$15	#        -"-
no_pred_symbol="-"
}
}

END {
for (i=1; i<=Nseq; i++ ) {
maxlen=max_len
if (sc2[i]>=sc3[i] && sc2[i]>=sc10[i]) {
	highest=sc2[i]
	CS_highest=pred2[i]	
	if (sc3[i] >= sc10[i]) {
		snd_highest=sc3[i]
		CS_snd_highest=pred3[i]
		trd_highest=sc10[i]
		CS_trd_highest=pred10[i]
	}	
	else {
		snd_highest=sc10[i]
                CS_snd_highest=pred10[i]
                trd_highest=sc3[i]
                CS_trd_highest=pred3[i]
	}
}

else {
	if (sc3[i]>=sc2[i] && sc3[i]>=sc10[i]) {
	        highest=sc3[i]
	        CS_highest=pred3[i]
		if (sc2[i] >= sc10[i]) {
			snd_highest=sc2[i]
	                CS_snd_highest=pred2[i]
	                trd_highest=sc10[i]
	                CS_trd_highest=pred10[i]
		}
		else {
	                snd_highest=sc10[i]
	                CS_snd_highest=pred10[i]
	                trd_highest=sc2[i]
	                CS_trd_highest=pred2[i]
	        }
	}
	else {
		highest=sc10[i]
	        CS_highest=pred10[i]
		if (sc2[i] >= sc3[i]) {
			snd_highest=sc2[i]
	                CS_snd_highest=pred2[i]
	                trd_highest=sc3[i]
	                CS_trd_highest=pred3[i]
		}
		else {
			snd_highest=sc3[i]
	                CS_snd_highest=pred3[i]
	                trd_highest=sc2[i]
	                CS_trd_highest=pred2[i]
		}
	}
}

# print the CS corresponding to the best score, IFF the predicted mTP length < maxlen
# also check which is the lowest of len[i] and maxlen; set maxlen till len[i] if 
if (maxlen>len[i]) maxlen=len[i]
if (CS_highest<=maxlen) print name[i]"\t"CS_highest"\t"highest
else {
	if (CS_snd_highest<=maxlen) print name[i]"\t"CS_snd_highest"\t"snd_highest
	else {
		if (CS_trd_highest<=maxlen) print name[i]"\t"CS_trd_highest"\t"trd_highest
		else print name[i]"\t-\t-"
	}
}

}
}

