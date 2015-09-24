=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::ClusterMotifs

=head1 DESCRIPTION

'ClusterMotifs' 

=cut


package Bio::EnsEMBL::Funcgen::RunnableDB::ClusterMotifs;

use warnings;
use strict;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Motif');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (run_system_cmd);
use Data::Dumper;

sub fetch_input {   # nothing to fetch... just the parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();
  
  return 1;
}

sub run {   
  my $self = shift @_;
  my @matrices;

  my $base = $self->_output_dir."/".$self->_feature_set->name;

  my $motif_file = $base.".transfac";
  #Maybe check if files are empty or not there at all!

  #Need to clump all motif files together...
  run_system_cmd("cat ".$self->_output_dir."/*.tmp_TRANSFAC > ".$motif_file);
  
  eval {
    #RUN STAMP here to cluster all motifs... change the Jaspar directories!!
    my $bin_folder = $self->_bin_folder();
    run_system_cmd("${bin_folder}/STAMP -tf $motif_file -align SWU -cc PCC -sd ${bin_folder}/STAMP_data/Jaspar2010_1000random_PCC_SWU.scores -chp -nooverlapalign -match ${bin_folder}/STAMP_data/Jaspar_Core_PBM_PolII.Transfac -match_top 1 -ma IR -out $base");
    
    #TODO Refactor this section...
    my $is_top = 0;
    my $top_matrix; #Is the one with smallest p
    my %matrix_scores;
    my $min_p = 1;
    open(FILE,$motif_file);
    while(<FILE>){
      my $line = $_;
      if($line =~ /^DE\s+(\S+)\s+(\S+)/){ 
	my $name = $1; my $p = $2;
	$matrix_scores{$name} = $p; 
	if($p < $min_p){
	  $top_matrix = "";
	  $is_top = 1;
	  $min_p = $p;
	  #import the top matrix...
	  #Also check if it matches any known matrix...
	  my $ff = $base."_match_pairs.txt";
	  if(-e $ff){ 
	    open(FIM,$ff);
	    while(<FIM>){
	      if(/^>\s*(\S+)\s+/){
		if($1 eq $name){
		  my $fline = <FIM>; chomp($fline); my ($jaspar, $score) = split(/\s+/, $fline);
		  if($score < 0.005){ $line =~ s/\n/\t${jaspar}\n/; }
		  last;
		}
	      }
	    }
	    close FIM;
	  } else { warn $ff." was not found!"; }
	}      
      }
      if($is_top){
	$top_matrix .= $line;
	if($line =~ /^XX/){ $is_top = 0; }
      }
      
    }
    close FILE;
    if($top_matrix){
      push @matrices, _transfac_to_jaspar($top_matrix);
    } else { warn "No top matrix found!!"; }
    
    #Just check what are the most similar matrices to the clusters...
    run_system_cmd("${bin_folder}/STAMP -tf ${base}_tree_clusters.txt -align SWU -cc PCC -sd ${bin_folder}/STAMP_data/Jaspar2010_1000random_PCC_SWU.scores -nooverlapalign -match ${bin_folder}/STAMP_data/Jaspar_Core_PBM_PolII.Transfac -match_top 1 -ma IR -out ${base}_global");
        
    #Recluster based on matrix similarity...
    push @matrices, @{_recluster_motifs($self->_feature_set->name, $self->_bin_folder(),$base,\%matrix_scores)};
  };
  if ($@){ warn "No motif found: ".$@; }

  #run_system_cmd("rm -f ".$self->_output_dir."/*.tmp_TRANSFAC");

  $self->_matrix_to_store(\@matrices);

  return 1;
}


sub write_output {  # Nothing is written at this stage (for the moment)
  my $self = shift @_;

  #Store the final matrices obtained...
  open(FO,">".$self->_output_dir."/".$self->_feature_set->name.".final");
  print FO join("\n",@{$self->_matrix_to_store()});
  close FO;

  return 1;
}

#TODO Need refactoring and optimization...
sub _recluster_motifs {
  my ($fset, $bin_folder, $base, $scores) = (shift, shift, shift, shift);
  my %clusters;
  my @result_matrix = ();

  open(FILE,$base."_tree_clusters.txt");
  while(<FILE>){
    if(/^DE\s+(\S+)\s*$/){
      #try to obtain a score here!...
      my $cluster_id = $1;
      my $matrix = $_;
      while(<FILE>){ $matrix .= $_; last if(/^XX/); }
      my $cluster_members = <FILE>;
      chomp($cluster_members);
      $cluster_members =~ s/^XX\s+\Cluster_Members:\s+//;
      $clusters{$cluster_id}{"matrix"}=$matrix;
      push @{$clusters{$cluster_id}{"elements"}}, split(/\s+/,$cluster_members);
    }
  }
  close FILE;

  #recalculate scores
  foreach my $cluster (keys %clusters){
    my $score = 0;
    map { $score += $scores->{$_}; } @{$clusters{$cluster}{"elements"}};
    $score = $score / scalar(@{$clusters{$cluster}{"elements"}});
    $clusters{$cluster}{"matrix"} =~ s/\n/\t${score}\n/;
    $clusters{$cluster}{"score"} = $score;
  }

  
  my %reclust;
  open(FILE,$base."_global_match_pairs.txt");
  while(<FILE>){
    chomp;
    if(/^>\s+(\S+)\s*$/){
      my $clust_id = $1;
      my $match = <FILE>;
      chomp($match);
      my ($jaspar, $score, undef, undef) = split(/\s/,$match);
      #Make this an alpha-parameter?
      if($score < 0.005){
	#add the matrix in the first line...
	$clusters{$clust_id}{"matrix"} =~ s/\n/\t${jaspar}\n/;
	push @{$reclust{$jaspar}}, $clust_id;
      }
    }
  }
  close FILE;
  
  foreach my $jaspar (keys %reclust) { 
    if(scalar(@{$reclust{$jaspar}})>1){
      my $score = 0; 
      map { $score += $clusters{$_}{"score"}; } @{$reclust{$jaspar}};
      $score = $score / scalar(@{$reclust{$jaspar}});
      my $clust_file = $base."_cluster_".$jaspar;
      open(FO,">".$clust_file);
      map { print FO $clusters{$_}{"matrix"}; } @{$reclust{$jaspar}};
      close FO;
      run_system_cmd("${bin_folder}/STAMP -tf $clust_file -align SWU -cc PCC -sd ${bin_folder}/STAMP_data/Jaspar2010_1000random_PCC_SWU.scores -chp -nooverlapalign -match ${bin_folder}/STAMP_data/Jaspar_Core_PBM_PolII.Transfac -match_top 1 -ma IR -out $clust_file");
      my $matrix;
      open(FILE,"${clust_file}FBP.txt");
      #Re-add the averaged score here!
      while(<FILE>){ 
	if(/^DE/){
	  $matrix = "DE ".$fset."_cluster_${jaspar}\t".$score."\t".$jaspar."\n";  
	} else { $matrix .= $_;  }  
      }
      close FILE;
      push @result_matrix, _transfac_to_jaspar($matrix);
      map { delete $clusters{$_}; } @{$reclust{$jaspar}};
    }
  }
  
  map { push @result_matrix,  _transfac_to_jaspar($clusters{$_}{"matrix"}); } keys %clusters;
  
  return \@result_matrix;

}


#private function that transforms a transfac matrix to jaspar format
sub _transfac_to_jaspar{
  my ($transfac) = shift;
  my $jaspar;

  my @lines = split(/\n/,$transfac);
  pop @lines; #XX
  my $title = shift @lines;
  $title =~ s/^DE/>/;
  $jaspar = $title."\n";
  my @as; my @cs; my @gs; my @ts;
  foreach my $line (@lines){
    my (undef,$a,$c,$g,$t,undef) = split(/\s+/,$line);
    #convert it to integers if necessary to make it simpler after
    if($a<=1 && $c<=1 && $g<=1 && $t<=1){ $a = int(100*$a); $c=int(100*$c); $g=int(100*$g); $t=int(100*$t); }
    push @as, $a; push @cs, $c; push @gs, $g; push @ts, $t;
  }
  $jaspar .= "A [ ".join("\t",@as)." ]\n" ;
  $jaspar .= "C [ ".join("\t",@cs)." ]\n" ;
  $jaspar .= "G [ ".join("\t",@gs)." ]\n" ;
  $jaspar .= "T [ ".join("\t",@ts)." ]\n" ;

  return $jaspar;
}

#Private getter / setter to the matrices
sub _matrix_to_store {
  return $_[0]->_getter_setter('matrix_to_store',$_[1]);
}

1;
