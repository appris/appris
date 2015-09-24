=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::InferMotifs

=head1 DESCRIPTION

'Funcgen::InferMotifs' 

=cut


package Bio::EnsEMBL::Funcgen::RunnableDB::InferMotifs;

use warnings;
use strict;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Motif');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (run_system_cmd);
use Data::Dumper;

sub fetch_input {   # nothing to fetch... just the parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();

  if(!$self->param('bin')){ throw "No bin defined!"; }
  $self->_bin($self->param('bin'));
  
  return 1;
}

sub run {   

  my $self = shift @_;
  
  my $afa = $self->_efgdba()->get_AnnotatedFeatureAdaptor();
  #Order features descending by score
  my @features = sort { $b->score <=> $a->score } @{$afa->fetch_all_by_FeatureSets( [ $self->_feature_set ] )};

  # Print the sequences for this bin
  my $sa = $self->_efgdba()->dnadb->get_SliceAdaptor();
  my $fasta_file = $self->_output_dir."/bin.".$self->_bin.".fasta";
  open(FO,">".$fasta_file) or throw "Couldn't open output file";
  my $bin_start = ($self->_bin - 1)*$self->_bin_size;
  for(my $i=$bin_start; $i<($bin_start+$self->_bin_size); $i++){
    my $ft = $features[$i];
    my $sr_name = $ft->seq_region_name;
    my $start = POSIX::floor($ft->summit - $self->_window_size);
    my $end = POSIX::floor($ft->summit + $self->_window_size);
    my $slice = $sa->fetch_by_region( undef, $sr_name, $start, $end);
    if(defined($slice)){ 
      print FO ">".$sr_name.",".$start.",".$end."\n".$slice->get_repeatmasked_seq()->seq."\n"; 
    } else { 
      warn $sr_name.",".$start.",".$end." could not be found\n";
    }
  }    
  close FO;

  #Create a background to try and avoid repetitive A's 
  #fasta-get-markov is part of the meme package (but is very simple to replicate if needed)
  #Use a pre-compiled background unique for the whole of human... 
  #run_system_cmd("cat ".$fasta_file." | fasta-get-markov -m 3 > ".$fasta_file.".bg");

  #Run MEME for this bin
  my $species = $self->_species;
  my $meme_file = $self->_output_dir."/bin.".$self->_bin.".meme";
  #Replace meme parameters by parameters from an analysis in the database...
  #run_system_cmd("meme.bin -nostatus -dna -text -revcomp -mod zoops -evt 0.00001 -nmotifs 10 -minsites 50 -minw 6 -maxw 20 -bfile ~/src/STAMP.v1.1/tests/${species}_masked.bg $fasta_file > ".$meme_file);
  run_system_cmd($self->_bin_folder()."/meme.bin -nostatus -dna -text -revcomp -mod zoops -evt 1e-20 -nmotifs 5 -minsites 50 -minw 6 -maxw 20 $fasta_file > ".$meme_file);

  my $basename = $self->_feature_set->name.".bin_".$self->_bin;
  if(_process_meme_to_STAMP($meme_file, $basename) == 0){
    warn "No motif found ";
  }

  #Eliminate the temp files...
  #run_system_cmd("rm -f ".$fasta_file);
  #run_system_cmd("rm -f ".$meme_file);

  return 1;
}


sub write_output {  # Nothing is written at this stage (for the moment)

  my $self = shift @_;

  return 1;

}


sub _process_meme_to_STAMP {
  my ($meme_file, $basename) = (shift, shift, shift);

  my $num_motifs = 0;
  my $log_odds = 1;
  open(FO,">".$meme_file.".tmp_TRANSFAC");
  open(FILE, $meme_file);
  while(<FILE>){
    chomp;
    if(/^MOTIF\s+.*E-value\s*=\s+(\S+)/) { $log_odds = $1; }
    if(/^BL\s+MOTIF\s+(\d+)\s+width=(\d+)\s+seqs=(\d+)\s*$/){
      my $motif_num = $1;
      my $motif_width = $2;
      my $num_seqs = $3;
      

      my @init = (0) x $motif_width;
      my %matrix;
      #Initialize with 0
      for (my $j=0; $j<$motif_width; $j++){ 
	$matrix{"A"}->[$j] = $matrix{"C"}->[$j] = $matrix{"G"}->[$j] = $matrix{"T"}->[$j] = 0;
      }

      for (my $i=1;$i<=$num_seqs;$i++){
	my $line = <FILE>;
	$line =~ /^\S+\s+\(\s*\d+\)\s+(\w+)\s+\d+\s*$/;
	my @seq = split(//,uc($1));
	if(scalar(@seq) != $motif_width){ throw "Problems parsing motif file"; } 
	for (my $j=0; $j<$motif_width; $j++){ $matrix{$seq[$j]}->[$j]++ }
      }
      my $end = <FILE>;
      if(! $end =~ /^\/\/\s*$/){ 
	throw "Error in the MEME file: expecting \/\/ after $num_seqs lines ";
      }
      print FO "DE ".$basename.".".$motif_num."\t".$log_odds."\n";
      for(my $i=0; $i<$motif_width; $i++){
	print FO $i."\t".$matrix{"A"}->[$i]."\t".$matrix{"C"}->[$i]."\t".$matrix{"G"}->[$i]."\t".$matrix{"T"}->[$i]."\n"
      }
      print FO "XX\n";

      $num_motifs++;

    }
    
  }
  close FILE;
  close FO;
  
  return $num_motifs;

}

#Private getter / setter to the bin
sub _bin {
  return $_[0]->_getter_setter('bin',$_[1]);
}

1;
