=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::ImportMotifFeatures

=head1 DESCRIPTION

'Import' is the base Runnable for the Import Pipeline

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::ImportMotifFeatures;

use base ('Bio::EnsEMBL::Hive::Process');


use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (run_system_cmd generate_slices_from_names);
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Importer;
#use Data::Dumper;

#global values for the Helper... maybe pass as parameters...
$main::_debug_level = 0;
$main::_tee = 0;
$main::_no_log = 1;

sub fetch_input {   # fetch parameters...
  my $self = shift @_;

  throw "No matrix given" if ! $self->param('matrix');
  #throw "No feature type given" if ! $self->param('feature_type');
  throw "No file given" if ! $self->param('file');
  my $file_name =  $self->param('file');
    
  ### Create DB connections
  my $coredb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
						   -user => $self->param('dnadb_user'),
						   -port => $self->param('dnadb_port'),
						   -host => $self->param('dnadb_host'),
						   -dbname => $self->param('dnadb_name'),
						  );


  ### Create DB connections
  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
							-dbname  => $self->param('dbname'),
							-port    => $self->param('port'),
							-pass    => $self->param('pass'),
							-host    => $self->param('host'),
							-user    => $self->param('user'),
							-dnadb   => $coredb,
						       );


  #test db connections
  $db->dbc->db_handle;
  
  $self->param('dba', $db);
  
  if(! -e $file_name ){ throw " Could not find ".$file_name; }
  
  my $fta = $db->get_FeatureTypeAdaptor();

  #We're associating to this matrix the associated feature types
  my $bma = $db->get_BindingMatrixAdaptor();

  #Get by name...
  my @bms = @{$bma->fetch_all_by_name($self->param('matrix'))};
  if(scalar(@bms) == 0){
    throw $self->param('matrix')." could not be found";
  }
  if(scalar(@bms) > 1){
    throw "More than one matrix with the same name is not currently supported";
  }
  my $bm = $bms[0];
  $self->param('matrix', $bm);

  #check if there is already data for this matrix stored... if there is, throw an error...  
  #TODO test specifically for data in slices...
  my $count = $db->dbc->db_handle->selectrow_array("select count(1) from motif_feature where binding_matrix_id=".$bm->dbID);
  if(($count>0) && !($self->param('slices'))){ throw "Data for ".$bm->name." already exists! Remove it first"; }

  my $ft = $bm->feature_type;
  my @fts = ($ft);

  #Should we be doing this?
  push @fts, @{$fta->fetch_all_by_association($ft)};

  my @fsets;
  my $fsa = $db->get_FeatureSetAdaptor();
  map { map { push @fsets, $_->dbID; } @{$fsa->fetch_all_by_FeatureType($_)} } @fts;
  
  my $output_file_name = $self->param('output_dir')."/annotated_features_".$bm->name.".tab";

  my $query = "select distinct sr.name, af.seq_region_start, af.seq_region_end, af.annotated_feature_id from ".
    " annotated_feature af, seq_region sr where af.seq_region_id=sr.seq_region_id and ".
      " feature_set_id in (".join(",",@fsets).")";
  my $cmd = "mysql -e \"".$query."\" -quick -h".$self->param('host')." -P".$self->param('port').
    " -u".$self->param('user')." -p".$self->param('pass').' '.$self->param('dbname')." >".$output_file_name;

  print $cmd."\n";

  system($cmd) && throw "Error dumping Annotated Features for ".$bm->name." $cmd";

  $self->param('output_file_name',$output_file_name);
  
  return 1;
}

sub run {   # Check parameters and do appropriate database/file operations...

  my $self = shift @_;
  my $sa = $self->param('dba')->get_SliceAdaptor;
  my $mfa = $self->param('dba')->get_MotifFeatureAdaptor;
  
  my $bm = $self->param('matrix');

  

  my $results = $self->param('output_dir')."/overlaps_".$bm->name.".tab";
  my $cmd = "perl ".$self->param('efg_src')."/scripts/miscellaneous/cooccur.pl ".$self->param('file')." ".$self->param('output_file_name')." > ".$results;
  system($cmd) && throw "Error executing coocur: $cmd";

  my @slices;
  if($self->param('slices')){
    @slices = split(/,/,$self->param('slices'));
  }

  my %data;
  my %slice_cache;
  open(FILE,$results);
  while(<FILE>){
    chomp;
    my ($sr,$start,$end,$desc,$score,$strand,$bm_name,$sr_af,$start_af,$end_af,$af_id) = split("\t");
    # only include those completely included in the annotated feature
    next if ($start < $start_af);
    next if ($end > $end_af);
    
    #Quick hack to only import specific slices...
    my $filter = 0;
    if(scalar(@slices)>0){
      $filter=1;
      foreach my $slice (@slices){ 
	if($slice eq $sr){ $filter = 0; }
      }
    }

    next if($filter);

    if(!defined($slice_cache{$sr})){
      my $slice = $sa->fetch_by_region('toplevel',$sr);
      if($slice) { 
	$slice_cache{$sr} = $slice;
      } else { 
	warn "Slice $sr not found: silently ignoring entry"; 
      }
    }

    #maybe double check overlap?? Should be fine, though...
    #maybe also cross-check bm_name with $bm->name?
    if($bm->name ne $bm_name){ warn "Entry is for $bm_name, not for ".$bm->name." : entry ignored"; next; }
    $data{$sr}{$start}{$end}{$strand}{'score'} = $score;
    #Carefull there may be duplicates. Maybe turn into an array instead
    $data{$sr}{$start}{$end}{$strand}{'assoc_feats'}{$af_id}=1;
    
  }
  close FILE;

  #New procedure may be of loading all Filtered matches, 
  # and storing all associated annotated features...
  #Carefull the score may need to be the relative_affinity...

  #Keep the threshold;
  my $min_relative_affinity=1; 
  foreach my $sr (keys %data){
    my $slice = $slice_cache{$sr};
    foreach my $start (keys %{$data{$sr}}){
      foreach my $end (keys %{$data{$sr}{$start}}){
	foreach my $strand (keys %{$data{$sr}{$start}{$end}}){
	  
	  my $motif_slice =  $sa->fetch_by_region('toplevel',$sr, $start, $end, $strand);
	  my $relative_affinity = $bm->relative_affinity($motif_slice->seq);
	  if($relative_affinity < $min_relative_affinity){ 
	    $min_relative_affinity = $relative_affinity;  
	  }
	  my $mf = Bio::EnsEMBL::Funcgen::MotifFeature->new
	    (
	     -slice          => $slice,
	     -start          => $start,
	     -end            => $end,
	     -strand         => $strand,
	     -binding_matrix => $bm,
	     #-score          => $data{$sr}{$start}{$end}{$strand}{'score'},
	     #Only round on the MotifFeature, not the BindingMatrix
	     -score          => sprintf("%.3f", $relative_affinity),
	    );
	  if($mf){
	    $mfa->store($mf); #($mf) = store gives null? 
	    foreach my $af_id (keys %{$data{$sr}{$start}{$end}{$strand}{'assoc_feats'}}){
	      my $sql = "INSERT INTO associated_motif_feature (annotated_feature_id, motif_feature_id) VALUES ($af_id, ".$mf->dbID.")";
	      $self->param('dba')->dbc->do($sql);
	    }
	    
	    #A safer way is to fetch the af and save it through the API
	    #$mfa->store_associated_AnnotatedFeature($mf,$af); 
	  }
	}
      }
    }
  }
  
  #Update the matrix threshold (need to do it directly in SQL?)  
  my $sql = "UPDATE binding_matrix set threshold=".$min_relative_affinity." where binding_matrix_id=".$bm->dbID;
  $self->param('dba')->dbc->do($sql);

  return 1;
}


sub write_output {  
  my $self = shift @_;
  
  return 1;

}

1;
