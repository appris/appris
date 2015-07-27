=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::MakeDnaseProfile

=head1 DESCRIPTION

=cut


package Bio::EnsEMBL::Funcgen::RunnableDB::MakeDnaseProfile;

use warnings;
use strict;
use Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor;
use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;

sub fetch_input {   # nothing to fetch... just the parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();
  
  my $work_dir = $self->param('work_dir');
  my $matrix = $self->param('matrix');
  my $dnase = $self->param('dnase');
  #throw "Dnase file $dnase does not exist in $work_dir" if(! -e $work_dir."/".$dnase);
  throw "Matrix file $matrix does not exist in $work_dir" if(! -e $work_dir."/matches/".$matrix);
  my $is_male = $self->param('is_male');
  throw "Need to define is_male" if(!defined($is_male));
  #throw "mapability.bedGraph file expected in $work_dir" if(! -e $work_dir."/mapability.bedGraph");

  return 1;
}

sub run {   # Create Groups and Analysis and stuff...
  my $self = shift @_;

  my $hlen = 100;
  my $work_dir = $self->param('work_dir');
  my $matrix = $self->param('matrix');
  my $dnase = $self->param('dnase');
  my $is_male = $self->param('is_male');

  my $motif_size;
  my %motifs;
  #Or remove // from matrix name
  open(FILE,$work_dir."/matches_vf/".$matrix);
  <FILE>; #title
  while(<FILE>){
    chomp;  
    #10	71532	71542	Gata1	8.55866	1	MA0035.2
    my ($chr,$start,$end,$name,$score,$strand)=split("\t");
    next if(!$is_male && ($chr =~ /Y/));
    #Ignore Mitochondria...
    next if($chr =~ /MT/);
    if(!defined($motif_size)){ $motif_size = $end - $start; }
    $motifs{$chr}{$start}{$strand}=$score;
  }
  close FILE;

  #Maybe pass this as parameters
  my $out_file = $work_dir."/output/".$matrix."_".$dnase.".counts";
  open(FO,">".$out_file);
  
  foreach my $chr (sort keys %motifs){
    my %datap = ();
    my %datan = ();
    foreach my $start (sort keys %{$motifs{$chr}}){
      for(my $i=$start-$hlen;$i<$start+$hlen+$motif_size;$i++){
	$datap{$i} = 0;
	$datan{$i} = 0;
      }
    }

    #Add a mappability score to extra filter...
    #A lot more time and didn't see much evidence of improvement...
    #my %mappability = ();
    #foreach my $start (sort keys %{$motifs{$chr}}){
    #  for(my $i=$start-$hlen;$i<$start+$hlen+$motif_size;$i++){
    #	$mappability{$i} = 0;
    #  }
    #}
    #open(FILE,$work_dir."/mapability.bedGraph");
    #while(<FILE>){
    #  chomp;
    #  my ($cur_chr,$start,$end,$score)=split();
    #  next if($cur_chr ne $chr);
    #for(my $i=$start; $i<$end; $i++){
    #	if(defined($mappability{$i})){ $mappability{$i}=$score; }
    #  }
    #}
    #close FILE;
    
    my $folder;
    if($is_male){
      $folder = "male";
    } else {
      $folder = "female";
    }
    #open(FILE,$work_dir."/dnase/".$folder."/".$dnase);
    open(FILE,"gzip -dc ${work_dir}/dnase/${folder}/${dnase}".'AlnRep*.bam.unique.tagAlign.gz |');
    while(<FILE>){
      chomp;
      my ($cur_chr,$start,$end,undef,undef,$strand)=split("\t");
      next if(($cur_chr ne $chr) || !defined($datap{$start}));
      if($strand eq '+'){
	$datap{$start}++;
      } else {
	$datan{$end-1}++;
      }
    }
    close FILE;
    
    foreach my $start (sort keys %{$motifs{$chr}}){
      foreach my $strand (keys %{$motifs{$chr}{$start}}){
	my $score = $motifs{$chr}{$start}{$strand};
	
	#filter out those that have no mappability in more than 20% of their extension...
	#my $count = 0;
	#for(my $i=$start-$hlen;$i<$start+$hlen+$motif_size;$i++){
	#  if($mappability{$i} <= 0.5){ $count++; }
	#}
	#Remove if more than 20% of region is "non_mappable" (approach from Centipede paper)
	#if(($count / (2*$hlen + $motif_size)) > 0.2){
	#  warn "Motif at chr:  $chr start: $start strand: $strand was ignored due to low mappability";
	#  next;
	#}

	print  FO $chr."\t".$start."\t".($start+$motif_size)."\t".$matrix."\t".$score."\t".$strand;
	for(my $i=$start-$hlen;$i<$start+$hlen+$motif_size;$i++){
	  if($strand eq '+'){
	    print FO "\t".$datap{$i};
	  } else {
	    print FO "\t".$datan{$i};
	  }
	}
	for(my $i=$start-$hlen;$i<$start+$hlen+$motif_size;$i++){
	  if($strand eq '+'){
	    print FO "\t".$datan{$i};
	  } else {
	    print FO "\t".$datap{$i};
	  }
	}
	print FO "\n";
      }
    }
  }
  
  close FO;

  return 1;
}


sub write_output {  # Nothing to do here
  my $self = shift @_;
  
  $self->dataflow_output_id($self->input_id, 2);

  return 1;


}

1;
