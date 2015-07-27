=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::SetupMotifInference

=head1 DESCRIPTION

'SetupMotifInference' 

=cut


package Bio::EnsEMBL::Funcgen::RunnableDB::SetupMotifInference;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;
use POSIX qw(floor);


use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Motif');



sub fetch_input {   # nothing to fetch... just the parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();

  if(! -d $self->_output_dir){ 
    system('mkdir -p '. $self->_output_dir) && throw "Error creating output dir ". $self->_output_dir;
  }

  return 1;
}

sub run {   # Create Subtasks of binsize peaks each, ignoring the last set of peaks ( < binsize peaks )
  my $self = shift @_;

  my $afa = $self->_efgdba()->get_AnnotatedFeatureAdaptor();
  my @features = @{$afa->fetch_all_by_FeatureSets( [ $self->_feature_set ] )};
  my $bins = POSIX::floor(scalar(@features)/$self->param('bin_size'));
  if($bins < 1){ 
    warn "Insuficient peaks. Please select a smaller bin size."; 
  }
  warn "Number of bins is $bins";

  #Create jobs
  my @bin_input_ids;
  for (my $i=1;$i<=$bins;$i++){
    #Need to add the specific file to the input_id...
    my $new_input_id = eval($self->input_id);
    $new_input_id->{"bin"} = $i;
    push(@bin_input_ids, $new_input_id);
  }
  $self->_output_ids(\@bin_input_ids);

  return 1;
}


sub write_output {  # Nothing is written at this stage (for the moment)
  my $self = shift @_;

  if($self->_output_ids && scalar($self->_output_ids)>0){
    my ($converge_job_id) = @{ $self->dataflow_output_id($self->input_id, 3, { -semaphore_count => scalar(@{$self->_output_ids}) })  };
    $self->dataflow_output_id($self->_output_ids, 2, { -semaphored_job_id => $converge_job_id });
  }
  return 1;

}

#Private getter / setter to the output ids
sub _output_ids {
  return $_[0]->_getter_setter('output_ids',$_[1]);
}

1;
