=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::Import

=head1 DESCRIPTION

'Import' is the base Runnable for the Import Pipeline

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::Import;

use base ('Bio::EnsEMBL::Hive::Process');

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (strip_param_args generate_slices_from_names strip_param_flags run_system_cmd);
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Importer;
#use Data::Dumper;

#global values for the Helper... maybe pass as parameters...
$main::_debug_level = 0;
$main::_tee = 0;
$main::_no_log = 1;

sub fetch_input {   # fetch parameters...
  my $self = shift @_;

  if(!defined($ENV{EFG_SRC})){ throw "NEED to define EFG_SRC"; }

  #Either test $ENV{EFG_DATA} or use a changed version of Importer.pm 
  
  $self->param('name',$self->param('input_set'));
  $self->param('output_dir',$self->param('output_dir')."/".$self->param('input_set'));

  $ENV{EFG_DATA} = $self->param('output_dir');

  #if($self->param('prepared')){
  #  $self->param('result_file', $self->param('output_dir')."/".$self->param('result_file'))
  #} else {
  #  
  #}

  my $Imp = Bio::EnsEMBL::Funcgen::Importer->new
    (
     -name        => $self->param('name'),
     -format      => $self->param('format'),
     -vendor      => $self->param('vendor'),
     -parser      => $self->param('parser'),
     -dbname      => $self->param('dbname'),
     -pass        => $self->param('pass'),
     -host        => $self->param('host'),
     -user        => $self->param('user'),
     -port        => $self->param('port'),
    -dnadb_host   => $self->param('dnadb_host'),
    -dnadb_port   => $self->param('dnadb_port'),
    -dnadb_user   => $self->param('dnadb_user'),
    -dnadb_name   => $self->param('dnadb_name'),
     #    -registry_pass => $self->param('registry_pass'),
     -registry_host => $self->param('registry_host'),
     -registry_user => $self->param('registry_user'),
     -registry_port => $self->param('registry_port'),
	 -release       => $self->param('registry_version'),
     #   -ssh         =>  $ssh,
     -group       => $self->param('group'),
     -location    => $self->param('location'),
     -contact     => $self->param('contact'),
     #   -array_set   => $array_set,
     -input_set_name => $self->param('input_set'),
     -input_feature_class => $self->param('input_feature_class'),
     #   -array_name  => $array_name,
     -result_set_name => $self->param('input_set'), #not implemented yet
     -feature_type_name => $self->param('feature_type'),
     -feature_analysis => $self->param('feature_analysis'),
     -cell_type_name => $self->param('cell_type'),
     #   -write_mage    => $write_mage,
     #   -update_xml => $update_xml,
     #   -no_mage => $no_mage,
     -assembly => $self->param('assembly'),
     -data_dir   => $self->param('data_dir'),
     -output_dir  => $self->param('output_dir'),
     -recover     => $self->param('recover'),
     #   -dump_fasta  => $fasta,
     #   -norm_method => ,
     -species     => $self->param('species'),
     -farm        => $self->param('farm'),
     -batch_job   => $self->param('batch_job'),
     -prepared    => $self->param('prepared'),
     -verbose     => $self->param('verbose'),
     -input_dir   => $self->param('input_dir'),
     #  -exp_date     => $exp_date,
     -result_files   => [ $self->param('result_file') ],
     -total_features => $self->param('total_features'),
     #  -old_dvd_format => $old_dvd_format,
     #  -ucsc_coords => $ucsc,
     #  -release => $release,
     _no_log => 1,
     -force  => 1,
    );
  
  if(!$Imp){ throw "Could not create importer"; }
  $self->param('importer', $Imp);

  #print "Output_dir: ".$Imp->get_dir('output')."\n";
  
  return 1;
}

sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift @_;
  
  my $Imp = $self->param('importer');
  
  if($self->param('wrap_up')){

    # run the merge script here...
    my $cmd = $ENV{EFG_SRC}."/scripts/import/merge_and_index_collections.pl ".
      " -dbhost ".$self->param('host').
	" -dbport ".$self->param('port').
	  " -dbuser ".$self->param('user').
	    " -dbpass ". $self->param('pass').
	      " -dbname ". $self->param('dbname').
		" -data_dir ". $self->param('output_dir').
		  " -result_set_name ".$self->param('input_set'),
		  		  " -dnadb_host ". $self->param('dnadb_host').
		  " -dnadb_port ". $self->param('dnadb_port').
		  " -dnadb_user ". $self->param('dnadb_user').
		  " -dnadb_name ". $self->param('dnadb_name');
    run_system_cmd($cmd);
    

    #set appropriate states...
    #dont get the result set directly from $Imp as it is not initialized...
    my $rseta = $Imp->db->get_ResultSetAdaptor();
    my ($rset) = @{$rseta->fetch_all_by_name($self->param("name"))};
    if(!$rset){ throw "Could not find ResultSet"; }
    $rseta->set_imported_states_by_Set($rset);
    return 1;
  }

  my @slices;
  my @skip_slices;
  #Allow for a list of slices as input and not only one?
  push(@slices,  $self->param('slice')) if $self->param('slice');
  #@slices = @{&generate_slices_from_names($Imp->slice_adaptor,\@slices, \@skip_slices, $toplevel, $nonref, $incdups)};
  @slices = @{&generate_slices_from_names($Imp->slice_adaptor, \@slices, \@skip_slices, 1, 0, 1)};#toplevel, nonref, incdups
  $Imp->slices(\@slices);
  
  if(!$self->param('prepared')){

    $Imp->init_experiment_import; #Define and store sets once here before setting off parallel farm jobs.

    #Preparing data...
    $Imp->read_and_import_data('prepare');

    #and now prepare all the jobs to be run...
    my @rep_out_ids;
    #Create the necessary LoadReads jobs
    #Get the slices from the preparation process itself?
    foreach my $slice (@slices){
      my $new_input_id = eval($self->input_id);
      $new_input_id->{"slice"} = $slice->seq_region_name;
      $new_input_id->{"result_file"} = $Imp->output_file;
      $new_input_id->{"total_features"} = $Imp->counts('total_features');
      push(@rep_out_ids,$new_input_id);
    }
    
    # Carefull with flow numbers... 6,7... maybe pass as parameter??
    #Wrapup job
    my ($funnel_job_id) = @{ $self->dataflow_output_id($self->input_id, 2, { -semaphore_count => scalar(@rep_out_ids) } ) };
    #All the fanned jobs...
    $self->dataflow_output_id(\@rep_out_ids, 1, { -semaphored_job_id => $funnel_job_id } );

    
  } else {
    
    #eval {
      #This eval is because a few will crash because the slice is not in the sequence!!
      #Check that this is working as it should!...
      $Imp->register_experiment();
    #};
    #if($@){ warn "Carefull with possible failure in import: $@"; }
  }

  return 1;
}


sub write_output {  # Create the relevant jobs
  my $self = shift @_;
  
  return 1;

}

1;
