
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

Bio::EnsEMBL::Funcgen::Hive::Utils

=head1 DESCRIPTION



=cut

package Bio::EnsEMBL::Funcgen::Hive::Utils;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );  

use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  inject_DataflowRuleAdaptor_methods  
 ); 

sub inject_DataflowRuleAdaptor_methods{
  my $dfr_adaptor = shift;
  my $helper      = shift;
  
  assert_ref($dfr_adaptor, 'Bio::EnsEMBL::Hive::DBSQL::DataflowRuleAdaptor');
  
  if(defined $helper){
    assert_ref($helper, 'Bio::EnsEMBL::Funcgen::Utils::Helper', 'Helper');  
  }
  
  my $dfr_ref = ref($dfr_adaptor);
   
  if(! $dfr_adaptor->can('fetch_all_by_analysis_id')){      #inject method here   
  
    if($helper){
      $helper->debug(1, 'Injecting '.$dfr_ref..'::fetch_all_by_analysis_id');    
    }
    
    no strict 'refs';
  
    *{$dfr_ref.'::fetch_all_by_analysis_id'} = sub { 
      my $self    = shift;
      my $anal_id = shift;
      throw('Must provide and analysis id argument') if ! defined $anal_id;
      return $self->fetch_all("from_analysis_id=${anal_id}");
    }; 
   
    use strict;
  }
    
  if(! $dfr_adaptor->can('get_dataflow_config_by_analysis_id')){  #inject method here using fetch_all_by_analysis_id 
    

    if($helper){
      $helper->debug(1, 'Injecting '.$dfr_ref..'::get_dataflow_config_by_analysis_id');
    }
    
    no strict 'refs';
   
    *{$dfr_ref.'::get_dataflow_config_by_analysis_id'} = sub { 
      my $self    = shift;
      my $anal_id = shift;
      throw('Must provide and analysis id argument') if ! defined $anal_id;        
      my %df_config;
        
      foreach my $dfr(@{$self->fetch_all_by_analysis_id($anal_id)}){
        #$anal_id here always represents the from analysis
        my $to_analysis = $dfr->to_analysis->logic_name;
          
        #Is it valid to wire to the same analysis using two different branches
        #We need to catch this and throw          
        if(exists $df_config{$to_analysis}){
          throw('It appears that the pipeline configuration for '.$to_analysis.
            " has been wired via two separate branches:\t".$dfr->branch_code.
            ' & '.$df_config{$to_analysis}{branch}.
            "\nDynamic branch dataflow currently only supports 1 assoicated branch");  
        }
          
        $df_config{$to_analysis} ||= {branch => $dfr->branch_code, funnel=>undef };
          
        if($dfr->funnel_dataflow_rule_id){
          $df_config{$to_analysis}{funnel} = 
            $self->fetch_by_dbID($dfr->funnel_dataflow_rule_id)->to_analysis->logic_name;
        }
          
        #This will result in redundant branche entries across the to_analysis values
        #which is fine, we just need to handle this in the caller
          
        #How are we going to handle the potential of dataflowing down the same branch twice, due 
        #to the redundancy of the branches wrt to_analyis keys
        #this will warn, but not fail, and dataflow will not happen
        #so is safe-ish
          
        #This should be skipped over somehow
        #we could do this be caching refs to identify unique outputids
        #although there is nothing to stop the same output id being passed in a different array/hash
        #and hence a different ref
        #The way the runnable would and should be naturally written will prevent this
        #It will just be that a single analysis name will be used to flow to all?
        #we could add support branch_job_group by taking an array of analysis names?
        #so the args would be
        #$self->branch_job_group([[anal1, ...], [jobid1, ...], {funnel_analysis_name=>[funnel_job1, ...]}]);
        
        #branch_job_group will then check that all the analyses are on the same branch
        #and that they constitute all of the analyses i.e.
        #you don't add config without adding support and vice versa
        #this does not prevent the dynamic branch dataflow, just what is 
        #dataflown
          
        #Having the runnable directly linked to the logic names like this is 
        #normally a bad idea, but this is the sacriface made to enable
        #dynamic branch dataflow, and is better than having to proliferate
        #hardcoded dataflow code
      }
        
      return \%df_config;
    }; 
   
    use strict;
  }
  
  
  #This could actaully be injected to any adaptor as it is not dependant on
  #other DataflowRuleAdaptor functionality
  
  if(! $dfr_adaptor->can('get_semaphoring_analysis_ids_by_logic_name')){
  
    if($helper){
      $helper->debug(1, 'Injecting '.$dfr_ref..'::get_semaphoring_analysis_ids_by_logic_name');    
    }
    
    no strict 'refs';
  
    *{$dfr_ref.'::get_semaphoring_analysis_ids_by_logic_name'} = sub { 
      my $self   = shift;
      my $lname  = shift;
      throw('Must provide and analysis logic_name argument') if ! defined $lname;
      
      my $sql = 'SELECT dr.from_analysis_id from dataflow_rule dr '.
        'JOIN dataflow_rule dr1 on dr.dataflow_rule_id=dr1.funnel_dataflow_rule_id '.
        "WHERE dr.to_analysis_url='$lname'";
        
      #return array to allow testing in scalar context
      return @{$dfr_adaptor->db->dbc->sql_helper->execute_simple($sql)};
    }; 
   
    use strict;
  }
  
    
  return $dfr_adaptor;    
}