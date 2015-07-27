#! /usr/sbin/perl

# USAGE: calc_of_motif_score-how-files.R-2-3-10.short-out.pl $cstype $mTPCS_max_len\
#		$MATRIX/R-10.30.-16-5.meme.matrix $TARGETPTMP/T_pred.how > \
#		$TARGETPTMP/R-10.res

$cstype=$ARGV[0];	# cstype is suppposed to be 2, 3, or 10 (corresponding to Arg position)
$mTPmaxlen=$ARGV[1];	# the longest tolerated mTP length; CS searched
			# in this area + $cs_lok figure (ie. 120+2/3/10)
$matrixfile= $ARGV[2];
$seqfile= $ARGV[3];

$Nseq=0;
$i=0;
$offset=0; # only for backward compatibility

## GET THE SCORING MATRIX
open(MATRIXFILE,$matrixfile);
while (<MATRIXFILE>) {
	chomp;
	if (/^##/) {
		($grb1,$grb2,$win_len,$grb3,$cs_lok)=split;
		print STDOUT "# $cs_lok  $matrixfile\t$seqfile\n";
		# $mTPmaxlen += $cs_lok;
	}
	if (/^ #/) {
		($grb1,$sc[1],$sc[2],$sc[3],$sc[4],$sc[5],$sc[6],$sc[7],$sc[8],$sc[9],
			$sc[10],$sc[11],$sc[12],$sc[13],$sc[14],$sc[15],$sc[16],$sc[17],
			$sc[18],$sc[19],$sc[20])=split;
		$i++;
		for ($j=1; $j<=20; $j++) {
			$MEME_matrix[($i-1)*20 + $j] = $sc[$j];
		}
	}
}
close MATRIXFILE;
$mTPmaxlen += $win_len;


## GET THE SEQUENCES
open (SEQFILE,$seqfile);
while (<SEQFILE>) {
	if ( substr($_,0,1) eq ">") {
		chomp;
		($field1)=split; 
		$name[$Nseq+1] = substr($field1,1,30);
#		print STDOUT  ">",$name[$Nseq],"\n",
#			substr($seqoneline[$Nseq],0,$mTPmaxlen),"\n"
#			if length($name[$Nseq]) > 0;
		$Nseq++;
	}
	else {
		chomp;
		$seqoneline[$Nseq] .= $_;
	}
}
# last entry...
#print STDOUT  ">",$name[$Nseq],"\n",substr($seqoneline[$Nseq],0,$mTPmaxlen),"\n";
close SEQFILE;

## CALCULATE THE SCORES, PRINT HIGHEST SCORE
for ($i=1; $i<=$Nseq; $i++) {
	$seq_len[$i]=length($seqoneline[$i]);
	$checklen[$i]=$mTPmaxlen-$win_len;
	$checklen[$i]=$seq_len[$i]-$win_len if $mTPmaxlen-$win_len > $seq_len[$i];
	$highest_score[$i]=-10000; # These 2 variables are used in the search 
	$highest_start_pos[$i]=0;  # for the highest SCORE (from the sc. matrix)
	for ($k=1; $k<=$checklen[$i]; $k++) {
		#next seqloop if $k > $mTPmaxlen; # don't check scores if we're outside allowed region
		#last if $k > $mTPmaxlen; # don't check scores if we're outside allowed region
		$score[$k]=0;
		for ($j=1; $j<=$win_len; $j++) {
		last if $k+$j-2 >= length($seqoneline[$i]);
			$score[$k]=$score[$k]+&get_score( substr($seqoneline[$i],$k+$j-2,1) ,$j);
		}
		if ($score[$k] > $highest_score[$i] ) {
			$highest_score[$i] = $score[$k];
			$highest_start_pos[$i]=$k
		}
	}
	$tp_len_assign[$i]=$highest_start_pos[$i] + $cs_lok + $offset + ($cstype - 2);
	print STDOUT "$name[$i] $seq_len[$i] $highest_score[$i] 0 $tp_len_assign[$i]\n";
}


## SUBROUTINE FOR EXTRACTING SCORES
sub get_score {
	$AA=$_[0];
	$position=$_[1];

	$aa_type=1  if ($AA eq "A"); 
        $aa_type=2  if ($AA eq "C");
        $aa_type=3  if ($AA eq "D");
        $aa_type=4  if ($AA eq "E");
        $aa_type=5  if ($AA eq "F");
        $aa_type=6  if ($AA eq "G");
        $aa_type=7  if ($AA eq "H");
        $aa_type=8  if ($AA eq "I");
        $aa_type=9  if ($AA eq "K");
        $aa_type=10 if ($AA eq "L");
        $aa_type=11 if ($AA eq "M");
        $aa_type=12 if ($AA eq "N");
        $aa_type=13 if ($AA eq "P");
        $aa_type=14 if ($AA eq "Q");
        $aa_type=15 if ($AA eq "R");
        $aa_type=16 if ($AA eq "S");
        $aa_type=17 if ($AA eq "T");
        $aa_type=18 if ($AA eq "V");
        $aa_type=19 if ($AA eq "W");
        $aa_type=20 if ($AA eq "Y");
        return -10  if ($AA eq "-"); # -10 to penalize the incorporation of nonexisting aa's
        return -10  if ($AA eq "");

	return $MEME_matrix[($position-1)*20 + $aa_type]
}
