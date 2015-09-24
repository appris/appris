=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::ImportMotifFeatures

=head1 DESCRIPTION

'Import' is the base Runnable for the Import Pipeline

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::ImportMotifFeatures;

use warnings;
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( open_file run_system_cmd );
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::MotifFeature;

use base ('Bio::EnsEMBL::Hive::Process');

#global values for the Helper... maybe pass as parameters...
$main::_debug_level = 0;
$main::_tee = 0;
$main::_no_log = 1;

sub fetch_input {   # fetch parameters...
  my $self      = shift;
  my $matrix    =  $self->param_required('matrix');
  my $file_name =  $self->param_required('file');
  throw("File does not exist:\t".$file_name) if ! -f $file_name;
    
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
  
  
  my $fta = $db->get_FeatureTypeAdaptor;
  my @bms = @{$db->get_BindingMatrixAdaptor->fetch_all_by_name($matrix)};

  #Get by name...
  if(scalar(@bms) != 1){
    throw("Failed to find unique BindingMatrix with name:\t".$matrix);
  }

  $matrix = $bms[0];
  $self->param('matrix', $matrix);

  #check if there is already data for this matrix stored... if there is, throw an error...  
  #TODO test specifically for data in slices...
  my ($count) = $db->dbc->db_handle->selectrow_array("select count(*) from motif_feature where binding_matrix_id=".$matrix->dbID);
  
  if( ($count > 0) && 
      ! $self->param('slices')){ 
    throw("Data for ".$matrix->name." already exists! Remove it first"); 
  }

  my $ft  = $matrix->feature_type;
  my @fts = ($ft);

  #Should we be doing this?
  push @fts, @{$fta->fetch_all_by_association($ft)};

  my @fsets;
  my $fsa = $db->get_FeatureSetAdaptor();
  map { map { push @fsets, $_->dbID; } @{$fsa->fetch_all_by_FeatureType($_)} } @fts;
  
  my $peaks_bed = $self->param('output_dir')."/annotated_features_".$matrix->name.".bed";


  #These are also create for pwm_filter_mappings!
  #Move this to a fan job in the hive.
  #Also need to make this dump bed and update cooccur to bedtools

  #This is not sorted!?

  my $query = "select distinct sr.name, (af.seq_region_start -1), af.seq_region_end, af.annotated_feature_id from ".
    " annotated_feature af, seq_region sr where af.seq_region_id=sr.seq_region_id and ".
      " feature_set_id in (".join(",",@fsets).")";
  my $cmd = 'mysql --skip-column-names -e "'.$query.'" -quick -h'.$self->param('host').
   ' -P'.$self->param('port').' -u'.$self->param('user').' -p'.$self->param('pass').' '.
   $self->param('dbname').' | sort -k1,1 -k2,2n -k3,3n > '.$peaks_bed;

  print $cmd."\n";

  run_system_cmd($cmd);
  $self->param('peaks_bed', $peaks_bed);
  
  return 1;
}

sub run {   # Check parameters and do appropriate database/file operations...
  my $self    = shift;
  my $sa      = $self->param('dba')->get_SliceAdaptor;
  my $bm      = $self->param('matrix');
  my $results = $self->param('output_dir')."/overlaps_".$bm->name.".tab";
 
  #To ensure full overlap we need to use -f 1 and put the annotaed features as the second file
  my $cmd = 'bedtools intersect -sorted -f 1 -wa -wb -a '.
   $self->param('file').' -b '.$self->param('peaks_bed')." > $results";
  run_system_cmd($cmd);

  my %slice_cache;

  if($self->param('slices')){
    map {$slice_cache{$_} = $sa->fetch_by_region('toplevel', $_)} 
     @{$self->param('slices')};
  }
  else{
    #Here we allow all non-ref slices
    #Assuming features are on genomic coords (i.e. inc dups for Y)
    map {$slice_cache{$_->seq_region_name} = $_} 
     @{$sa->fetch_all('toplevel', undef, 1, 1)} ;
  }

  my ($line, $slice, $rel_aff, %skipped_slices);
  #Set these here first, so they are automatically updated below.
  my $motif_features = {};
  my $af_ids         = {};
  $self->param('motif_features', $motif_features);
  $self->param('annotated_feature_ids', $af_ids);  
  my $rfile = open_file($results);
  my ($mf_sr, $mf_start, $mf_end, $score, $mf_strand, $af_id, $cache_key);

  RECORD: while(($line = $rfile->getline) && defined $line){
    chomp $line;
    ($mf_sr, $mf_start, $mf_end, undef, $score, $mf_strand, undef, undef, undef, $af_id) = split("\t", $line);
    $cache_key = join(':', ($mf_sr, $mf_start, $mf_end, $mf_strand));
    #TODO Change this to use the SliceHelper slice_cache and get_Slice methods (currently n BaseDB)
    #Quick hack to only import specific slices... 
  
    if(! exists $slice_cache{$mf_sr}){

      if(! exists $skipped_slices{$mf_sr}){
        warn "Slice not specified or found, skipping seq region:\t$mf_sr\n";
        $skipped_slices{$mf_sr} = undef;
      }
        
      next RECORD;
    }

    if(! exists $motif_features->{$cache_key}){
      #This was using the API relative afifinity method to set the threshold score
      #but this is not the same score as used by Moods which was used to filter.
      #threshold does not appear to be used by API anywhere else  so have now omitted this
      #$motif_slice   = $sa->fetch_by_region('toplevel',$mf_sr, $mf_start, $mf_end, $mf_strand);
      #$rel_aff = $bm->relative_affinity($motif_slice->seq);
      #$min_rel_aff = $rel_aff if $rel_aff < $min_rel_aff;

      $motif_features->{$cache_key} = Bio::EnsEMBL::Funcgen::MotifFeature->new
       (-slice          => $slice_cache{$mf_sr},
        -start          => ($mf_start + 1),
        -end            => $mf_end,
        -strand         => $mf_strand,
        -binding_matrix => $bm,
        -score          => sprintf("%.3f", $score));#$relative_affinity));

      $af_ids->{$cache_key} = [];
    }

    push @{$af_ids->{$cache_key}}, $af_id;     
  }

  $rfile->close;


  #throw_no_retry here is non are loaded?
  #assuming it's problem with the intersect

  #This is not the threshold calculated by the pipeline but
  #the feature with the lowest threshold above that!
  #$self->param('threshold', $min_rel_aff);
  return;
}


sub write_output {  
  my $self   = shift;
  my $mfs    = $self->param('motif_features');
  my $af_ids = $self->param('annotated_feature_ids');
  my $dbc    = $self->param('dba')->dbc;
  my $mfa    = $self->param('dba')->get_MotifFeatureAdaptor;
  my ($mf, $sql);

  foreach my $key(keys %$mfs){
    $mf = $mfa->store($mfs->{$key})->[0];

    foreach my $af_id (@{$af_ids->{$key}}){
      $sql = "INSERT INTO associated_motif_feature (annotated_feature_id, motif_feature_id) VALUES ($af_id, ".$mf->dbID.")";
      $dbc->do($sql);
    }
  }

  #Now done in pwm_file_mappings storing the actual calcuted threshold.
  #$sql = 'UPDATE binding_matrix set threshold='.$self->param('threshold').' where binding_matrix_id='.$bm->dbID;
  #$dbc->do($sql);
  return;
}

1;
