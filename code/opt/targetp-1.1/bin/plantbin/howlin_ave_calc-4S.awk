#
# FILE howlin_ave_calc-4S.awk
# USAGE paste howlin.result.$$.? | nawk -f howlin_ave_calc-4S.awk -v noofseq=$antal -v Pcut=# -v Tcut=# -v Scut=# -v Ocut=#
#
# Calculate average over the 5 different NN for all 4 categories. 
# Chooses then the highest output value and makes the prediction from this, given that the node corresponding to 
# the predicted location shows a value > cutoff for that category.




BEGIN {		
		#cut_off=0.5 # cut-offs from command line in this version
                for (i=1; i<=noofseq*4; i++) score[i]=0
		no_cat=4  # number of categories
}

/^ #/ { 
	for (j=1; j<=no_cat; j++) {
                for (i=0; i<=4 ; i++) {
                      	score[($2-1)*no_cat+j]+=$(i*8+6) # collecting the 5 output values for one category of one seq.
               	}
               	score[($2-1)*no_cat+j]/=5	# averaging
		if (j<=no_cat-1) getline
	}
	
}

END {
# performing the prediction (winner-takes-all with cutoff restriction)
	for (i=1; i<=noofseq; i++) { 
		pred_score[i]=0.000
		pred[i]="?"	
		for (j=1; j<=no_cat; j++) {
			if (score[(i-1)*no_cat+j] > pred_score[i]) {
				pred_score[i]=score[(i-1)*no_cat+j]
				pred[i]=j
			}
		}
		if (pred[i]==1) {
			if (pred_score[i]>Pcut) pred[i]="C"
			else pred[i]="*"
		}
		if (pred[i]==2) {
			if (pred_score[i]>Tcut) pred[i]="M"
			else pred[i]="*"
		}
		if (pred[i]==3) {
			if (pred_score[i]>Scut) pred[i]="S"
			else pred[i]="*"
		}
		if (pred[i]==4) {
			if (pred_score[i]>Ocut) pred[i]="_"
			else pred[i]="*"
		}
		
		# Reliability class calculations
		highest=0
		snd_highest=0
		#highest_pos=0
		#snd_highest_pos=0
		for (j=1; j<=no_cat; j++) {
		        if (score[(i-1)*no_cat+j] > highest) {
        		        snd_highest=highest
        		        #snd_highest_pos=highest_pos
        		        highest=score[(i-1)*no_cat+j]
        		        #highest_pos=j
       			}
       			else {  
                		if (score[(i-1)*no_cat+j] > snd_highest) {
                		snd_highest=score[(i-1)*no_cat+j]
                		#snd_highest_pos=j
                		}   
       			}
		}
		diff= highest-snd_highest
		if (diff > 0.800) rel_class[i]=1
		else {
		        if (diff > 0.600) rel_class[i]=2
		        else {
		                if (diff > 0.400) rel_class[i]=3
		                else {
		                        if (diff > 0.200) rel_class[i]=4
		                        else rel_class[i]=5
		                }
		        }
		}					

	}

# presentation of results
	for (i=1; i<=noofseq; i++) {
		printf "%5.3f %5.3f %5.3f %5.3f %3s     %1d\n", \
		score[(i-1)*no_cat+1], score[(i-1)*no_cat+2], \
		score[(i-1)*no_cat+3], score[(i-1)*no_cat+4], pred[i], \
		rel_class[i]
	}
}
