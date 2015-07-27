#!/usr/bin/perl -w

if($#ARGV != 2) {
  printf "Usage: selex2mod.pl <protnamefile> <selexfiledir> <outfiledir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $selexfiledir = $ARGV[1]."/";
my $outfiledir = $ARGV[2]."/";

open(PROTNAMEFILE, "$protnamefile")
  || die "could not open $protnamefile\n";

while(<PROTNAMEFILE>) {
  chomp;
  my $protname = $_;

  my @seqs = ();
  my @seqnames = ();
  
  print "running $protname ... ";
  
  my $selexfile = "$selexfiledir"."$protname".".slx";
  open(SELEXFILE, "$selexfile")
  || die "could not open $selexfile\n";

  my $length = 0;
  my $seq_index = 0;
  my $first = "YES";
  while(<SELEXFILE>) {
    chomp;
    my $row = $_;
    if(length $row == 0) {
      $seq_index = 0;
      $first = "NO";
    }
    elsif(substr($row, 0, 1) eq '#' || substr($row, 0, 1) eq '%') {
      
    }
    else {
      #print "$row\n";
      my $seq = "";
      my $name = "";
      
      my $seq_pos = index($row, " ") + 1;
      #print "$seq_pos \n";
      
      while((length $seq) > 0 && (index $seq, " ") == 0) {
	my $new_seq = substr $seq, 1;
	$seq = $new_seq;
      }

      $seq = substr $row, $seq_pos;
      $seq =~ s/\s/-/g;
      $seq =~ s/_/-/g;
      $seq =~ s/\./-/g;
      #print "$seq\n";
      
      
      if($first eq 'YES') {
	$seqs[$seq_index] = $seq;
      }
      else {
	$seqs[$seq_index] .= $seq;
	#print "$seq\n";
      }
      $seq_index++;
    }
  }

  my $outfile = "$outfiledir"."$protname".".mod";
  open(OUTFILE, ">"."$outfile")
    || die "could not open $outfile\n";

  for(my $i = 0; $i <= $#seqs; $i++) {
    my $seq = $seqs[$i];
    $seq =~ s//;/g;
    $seq = substr $seq, 1;
    $seq = "<"."$seq".">\n";
    print OUTFILE $seq;
  }

  close OUTFILE;
  print "done\n";
}
