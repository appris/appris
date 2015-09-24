
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::RunAligner

=head1 DESCRIPTION



=cut

package Bio::EnsEMBL::Funcgen::Hive::RunAligner;

use warnings;
use strict;
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( validate_package_path );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

  #We need a list of Aligner specific param requirements which are not specified
  #in $analysis->parameters
  #bwa_index_root, gender, assembly fasta_fai
  #How are we going to genericise these?
  #We can't specify ResultSet->cell_type->gender
  #and we don't want to tie the aligner to use of a ResultSet
  #Can we make these co-optional in the constructor?
  #i.e. we can pass gender 

#Some of these are depeandant on the gender being used
#so we need a method to generate the relevant 


#NOTE: There is a run_aligner.pl script which will do much the same as this module
#      without the dependancy on having access to a DB or and Analysis object
  
sub fetch_input {   # fetch parameters...
  my $self = shift;
  #Set some module defaults
  $self->param('disconnect_if_idle', 1);
  $self->check_analysis_can_run;
  
  $self->SUPER::fetch_input();

  my $query_file = $self->param_required('query_file'); #fasta or fastq

  #$self->get_param_method('output_dir', 'required'); 
  #This should have been set to a 'work' dir in PreprocessFastqs 
  #Do we even need this? The fastq chunks will already be in a work dir?  

  my $logic_name = $self->param_required('analysis');
  my $analysis   = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
  #program is no passed to Aligner so validate here
  my $aligner      = $analysis->program || 
    throw('Aligner analysis cannot have an undef program attribute:'.$analysis->logic_name);
  my $align_module = validate_package_path($analysis->module); 

  #validate program_file isn't a path?
  #would redefine bwa bin in the config for this analysis
  #this will change the bin_dir for everything else too, but we don't use bin_dir for anything else here  
  my $pfile = $analysis->program_file;
  throw('Analysis '.$analysis->logic_name.' must have a program_file defined') if ! defined $pfile; 
  my $pfile_path = ( defined $self->bin_dir ) ? $self->bin_dir.'/'.$pfile : $pfile;
  #$self->set_param_method('program_file', $pfile_path)
  
  my $ref_fasta = $self->param_silent('indexed_ref_fasta');  # This is batch flown
   
  if(! defined $ref_fasta){ 
    my $gender         = $self->param_silent('gender') || 'male';
    my $species        = $self->species; 
  
    #TODO: Check if the index file is really there? Eventually the bwa output will tell you though
    #index file suffix may change between aligners
    #best to pass just the target file, index root dir, species, gender
    #and let the Aligner construct the appropriate index file
  
    $ref_fasta = join('/', ($self->param_required('data_root_dir'),
                            $aligner.'_indexes',
                            $species,
                            $species.'_'.$gender.'_'.$self->assembly.'_unmasked.fasta'));
  }
  
  #$self->set_param_method('target_file', $ref_fasta);

  my $aligner_methods = $self->get_param_method('aligner_param_methods', 'silent');
  my %aparams;

  if($aligner_methods){
    assert_ref($aligner_methods, 'ARRAY', 'aligner_param_methods');
    
    foreach my $method(@$aligner_methods){
      
      if(! $self->can($method)){
        #This must be a param we haven't seen yet 
        $aparams{'-'.$method} = $self->param_required($method);  
      }
      else{
        #We have a method which may build this param dynamically
        $aparams{'-'.$method} = $self->$method;   
      }  

      $self->helper->debug(1, "Setting $method aligner parameter:\t".$aparams{'-'.$method});
    }
    
    #$self->set_param_method('aligner_params', $aparams);
  }
  
  
  $self->helper->debug(1, "Creating aligner:\t".$align_module); 
 
  #TODO Change this to use a run script, which can then be used outside of the pipeline
  #This is fine so long as we are not passing any objects
  #This is fine, so long as we are capturing errors properly and 
  #passing them back up the stack for reporting
  
  my $align_runnable = $align_module->new
   (-program_file      => $pfile_path,
    -parameters        => $analysis->parameters, 
    -query_file        => $query_file,
    -target_file       => $ref_fasta,
    -debug             => $self->debug,
    %aparams                                    );
    
  $self->helper->debug(1, "Setting aligner:\t".$align_runnable); 
  $self->set_param_method('aligner', $align_runnable); 
  return;
}


sub run {
  my $self = shift;
    
  if(! eval { $self->aligner->run; 1; }){
    my $err = $@;
    $self->throw_no_retry('Failed to call run on '.ref($self->aligner)."\n$err"); 
  }
  
  $self->debug(1, "Finished running ".ref($self->aligner));
  return;
}


sub write_output {  # Nothing to write
  return;
}



1;
