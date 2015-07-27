# FILE in2how+fasta.awk
#!/usr/bin/nawk -f
# Translate sequences from input (fasta format or HOW format)
# to output (fasta format *and* HOW format)
# 
# These variables MUST be set:
#	informat: 	fasta | how
#	howout:		<filename>
#	fastaout:	<filename>
# These variables MAY be set:
#	maxlen:		<number>
#
# Non-letters are deleted from the sequence(s).
# Letters not in AA alphabet are converted to "X".
# Lowercase AA letters are converted to uppercase.
# Association character will be "." for all positions.

BEGIN { cut_len = 1000 		# to prevent the script getting stuck trying to deal
	len= -1 		# with too long sequences (that awk cannot handle)
	earlier_prelseq_len = 0
}

# Sequence header line (fasta format):
$0~/^>/ && informat=="fasta" {	
	if (prelseq!="")
		printseq(name,prelseq,len);

	header = $0
	
	# Use first word after ">" as name 
	if ($1==">") 
		name = $2 
	else         
		name = substr($1,2)
	prelseq="";
	prelseq2="";
	infasta=1;
	earlier_prelseq_len=0;
	next; 	# Don't include this line in sequence
}
# All other lines (fasta format):
# Add current line to preliminary sequence:
(infasta) { 
gsub("[^A-Za-z]","",$0);	# 06-03-09, KR, to handle non-aa in seq.
if (length(prelseq) < cut_len ) {
	prelseq = prelseq $0 
	len= length(prelseq)
}
else {				# this is maybe a bit opaque...
	prelseq2 = prelseq2 $0
	if (length(prelseq2) >= cut_len ) {
		earlier_prelseq_len = earlier_prelseq_len + length(prelseq2)
		prelseq2 = "";
	}
	len= length(prelseq) + length(prelseq2) + earlier_prelseq_len
}
}

# Sequence entry (how format):
#       Association lines are lost
#       stuff (such as position number) after column 80 is lost 
$0~/^[ 0-9]/ && informat=="how" {
	if (prelseq!="")
                printseq(name,prelseq,len);
	len = $1
	ctr_len = len
	if (len > cut_len) ctr_len = cut_len 
	name = $2
	prelseq="";
	header = ">" $2 " " $1 " " substr($0, 18)
	for(i=0; i<ctr_len; i+=80) {
                getline
                prelseq = prelseq substr($0,1,80)
        }	
}

END {	
	if (informat=="fasta" && prelseq!="")
		printseq(name,prelseq,len);
	if (informat=="how" && prelseq!="")
		printseq(name,prelseq,len);
}

function printseq(name,seq,len) {
	# Convert to uppercase:
	seq=toupper(seq);
	# Delete spaces, numbers, and punctuation:
	gsub(/[^A-Z]/,"",seq);
	# Convert non-AA letters to X:
	gsub(/[^ACDEFGHIKLMNPQRSTVWY]/,"X",seq);
	# Truncate:
	if (maxlen) {
		seq = substr(seq,1,maxlen)
		if (len > maxlen) len=maxlen
	}
	# Set association line to dots:
	ass=seq;
	gsub(/./,".",ass);
##	if (len > cut_len) len_to_print= len		# if true length > cut_len, only cut_len aa's
##	else len_to_print= length(seq)			# will be in the files, but in header the true
	ctr_len=length(seq);				# length will be shown
	if (ctr_len) { # If the sequence is non-empty
	
		# Print how output
		printf("%6d %-11s\n",len,substr(name,1,30)) > howout
	   	for (i=1; i<=ctr_len; i+=80)
       			printf("%-80s%7d\n",substr(seq,i,80),i+79) > howout
   		for (i=1; i<=ctr_len; i+=80)
       			printf("%-80s%7d\n",substr(ass,i,80),i+79) > howout
	
		# Print fasta output
		print header > fastaout	
   		for (i=1; i<=ctr_len; i+=80)
       			print substr(seq,i,80) > fastaout
	}
}
