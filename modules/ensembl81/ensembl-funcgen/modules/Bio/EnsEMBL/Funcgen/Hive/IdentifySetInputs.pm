
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

Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs;

=head1 DESCRIPTION

This module simply takes a list of paramters and a 'set_type' to identify 
Sets used as inputs for various parts of the analysis pipeline.

=cut

package Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects
                                               get_set_prefix_from_Set 
                                               validate_package_path );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );



#todo -slice_import_status?


my %set_adaptor_methods = 
 (
  InputSubset => 'get_InputSubsetAdaptor',
  ResultSet   => 'get_ResultSetAdaptor',
 );


#TODO 

#1 DONE Check we have experiment/experimetal_group support for the composable queries
#  validate set type if not all SetAdaptors support this.

#2 Could do with a way of listing embargoed, when no_write is defined

#3 DONE Implement is_InputSubset_downloaded in run

#4 Genericse no_idr into something like merge_reps_only? #This is aimed at restricting to non-idr sets
# and is not yet implemented. We really need a difference flag here. As no_idr is for turning off idr
# not set selection and maintaining idr. Change this to run_idr as a boolean? This is still a bit ambiguous

#5 Review force/ignore_embargoed/ not downloaded. Should we skip these by default? Or have separate skip options?
#  Skip non-downloaded by default unless -download is specified (it might be okay to down load a few files on the farm
#  but bulk download should be done on the head node, apparently systems prefer this
#  non cpu optimised overloading of a head node rather than tying up a farm node)
#  Then we are still left with the force/ignore_embargoed options.

#6 Review throw string building once we support download

#7 Revise control handling. This needs to be restricted to experimental sub groups
#  i.e. ENCODE_Yale not ENCODE
#  DONE experiments_like => 'K562_%_ENCODE_Yale', control_experiments => 'K562_WCE_ENCODE_Yale'
#  This will also identify control sets.

#  Drop allow_inter_group controls in favour of -control_experiments?
#  Also, will specifying control_experiments allow association with any other experiment
#  i.e. check whether we are still validating the control is within the project? We probably 
#  just want to trust the -control_experiments and simply make the assocition based on cell type.
#
#  Actually this is going to cause trouble if we have already processed some experiments which would match 
#  the string. The way around this is to simply add in the control_feature_types to a normal filter query
#  

# 8 DESPERATELY NEED THIS!! 
## select e.name from experiment e left join result_set rs on rs.name like concat(e.name, '%') where e.name like "HepG2%" and rs.name is NULL;
#We also need a no_write mode which quickly lists all those which have not be aligned yet
#   So this is an inverse of the ALIGNED status query. We need this for the other entry points too.
#   How will this protect against set which have already been seeded?
#   Do we need a SEEDED status? Would this be on the Experiment or InputSubset level?
#   This is only relevant to a particular instance of the hive. So this is not really safe.
#   We need a way to just skip over/stop flow of ResultSets which arleady exist
#   This is going to be tricky as this could be seen as a failure too
#   Maybe configure_hive and DropPipeline could check/delete these states?
#   So we really need a two prong approach here. 
#   1 Throw if the ResultSet exists and is aligned. (this is currently done in DefineResultSets)
#   2 Check if the Experiment has been seeded already? This does not apply to control experiments
#     Race condition between redundant controls is already being handled by ALIGNING_CONTROL in DefineResultSets
#   This is really tying the IdentifyInputSubset analysis to DefineResultSets, but it is much more logical
#   to do the checks here, rather than having to manage a failing DefinedResultSets job and reeseding
#   We need an ignore aligned mode. This will be mutually exclusive with the rollback result_set mode 
#
# 9 Split this out into separate subs for different sets. 
#   Move a lot of this to EFGUtils or Hive Utils so we can use 
#   some of this functionality outside of the pipeline

#
# 10 DONE Fully remove InputSet support

# 11 DONE Tidy up only_replicates/handle_replicates 



my @iset_only_params = qw( allow_no_controls skip_absent
                           force_embargoed  ignore_embargoed  identify_controls  control_experiments );
# No longer iset_only as defaults are present in BaseSequenceAnalysis config
#control_feature_types 
                          
my @rset_only_params = qw( only_replicates );

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input; 
  my $set_type    = $self->get_param_method('set_type', 'silent'); 
 
  #This auto detection prevents having to do this in the SeedPipeline environment func
  if(! defined $set_type){
    my $lname = $self->analysis->logic_name;  
    
    foreach my $stype( keys %set_adaptor_methods ){
  
      if($lname =~ /$stype/){
        
        if($set_type){ #Throw if we match is not unique
          throw("Cannot auto detect unique set_type from analysis logic name:\t$lname\n",  
                'Please defined the set_type analysis_parameter or rename analysis.');
        }
        
        $set_type = $stype;
      }
    }
    
    if(! defined $set_type){
      throw("Cannot auto detect unique set_type from analysis logic name:\t$lname\n",  
            'Please defined the set_type analysis_parameter or rename analysis.'); 
    } 
    
    $self->set_type($set_type);
  }
  
  
  ### Get the set type adaptor
    
  if(! exists $set_adaptor_methods{$set_type}){
    throw("The -set_type $set_type is not supported by IdentifySetInputs.\n".
          "Valid options are ".join(' ', keys(%set_adaptor_methods)) );      
  }
  else{
    #Set the method name value to reset the value as the actual adaptor
    my $method = $set_adaptor_methods{$set_type};
    $self->set_param_method('set_adaptor', $self->out_db->$method, 'required');
  }

  ### Get/set generic optional param methods
  #todo, split this into ref types, so we can validate
  #i.e. scalar|arrayref|hashref vars
  $self->process_params([qw(set_names set_ids)], 1);#1 is optional flag
  $self->get_param_method('experiments', 'silent');
  #can't process_params for experiments as it would actually fetch the Experiments
  #and fail with a widlcard
  $self->set_param_method('constraints_hash',
                          $self->process_params([qw(feature_types cell_types states
                                                    analyses experimental_groups)],
                                                1, 1));  #optional/as array flags
  #all these are OR filters except states which is an AND filter
                      

  ### Catch mutally exclusive params  
  if(keys %{$self->constraints_hash}){ 
      
    if($self->set_names || $self->set_ids || $self->experiments){
      throw('You have specified mutually exclusive filter params for the '.
            "IdentifySetInputs analysis\nPlease specify restrict to ".
            'set_names, set_ids, experiments_like or a combination other filters '.
            '(e.g. experimental_groups experiments feature_types cell_types analyses states');
    }
  }
  elsif( (grep {defined $_ } ($self->set_names, $self->set_ids, $self->experiments)) != 1){      
    throw('You must specify some IdentifySetInputs filter params either '.
          'set_names, set_ids, experiments_like or a combination of '.
          'feature_types cell_types experiments experimental_groups states analyses');
  }

  #Fetch set type specific inputs
  my $set_input_method = '_fetch_'.$set_type.'_input';
  $self->$set_input_method;

  return;
}


#todo Should identify_controls be on by default?
#add mode to fetch expected controls form the DB, based on known control ftypes, 
#and ctypes and experimental_groups of signal experiments. This will be problematic
#until we model subgroups as experimental groups. And may not even be necessary if we associate controls
#during registration


sub _fetch_InputSubset_input{
  my $self = shift;    
  $self->_validate_param_specifity(\@rset_only_params);     
  $self->process_params(\@iset_only_params, 1);#1 is optional flag 
  
  #Do we even need this now, tracking DB is mandatory
  #and ignore/force_embargoed handle other use cases
  #$self->set_param_method('validate_InputSubset_tracking', 1); #for convenience/speed in loop below
        
  if($self->force_embargoed && $self->ignore_embargoed){
    throw('force_embargoed and ignore_embargoed are mutually exclusive parameters. '.
     'Please omit one or both');     
  }
    
  #Get feature_set analysis for control validation in run
    
  #but we should also be mindful of a peak_analysis over-ride?
  #but need to use generic params here?
  #we could replace this with control_mandatory_ftypes
  
  #Genericise feature_set_analysis_type to feature_set_analysis
  my $anal_type = $self->param_required('feature_set_analysis_type');     
  $self->set_param_method('feature_set_analysis', 
                          $self->param_silent($anal_type.'_analysis'));

  #This give the following warning
  #ParamWarning: value for param('feature_set_analysis') is used before having been initialized!
  #should we change set/get_param_method to use param_silent by default?
                               
    if(! defined $self->feature_set_analysis){
      $self->get_param_method('default_feature_set_analyses', 'silent');     
 
      if(! defined $self->default_feature_set_analyses){
        throw("Please define -${anal_type}_analysis or add to default_feature_set_analyses".
          " in the default_options config");
      }    
    }
                     
    if($self->identify_controls && $self->allow_no_controls){
      throw("Mutually exclusive parameters have been specified:\t".
        "identify_controls & allow_no_controls\nPlease omit one of these.");  
    }
    
    if($self->control_experiments && $self->allow_no_controls){
      throw("Mutually exclusive parameters have been specified:\t".
        "control_experiments & allow_no_controls\nPlease omit one of these.");  
    }
    
    if($self->experiments){
#      warn "experiments is ".$self->experiments;
      #assert_ref
      
      #if($self->experiments !~ /^[^%]+_%_[^%]+$/o){
      #  #We could be more strict here and force the % to be in the expected ftype location
      #  throw('The -experiments_like string does not have a single % mysql wildcard '.
      #    "where the expected feature type should be:\n\t".$self->experiments_like);
      #}
      
      #do we assume identify_controls here?
      #This was originally working by harvesting control subsets from other signal experiments
      #This would still work the same, as the experiments_like string would pick up the controls
      #and these will be marked as such in the DB. so we don't even need the control_feature_types
      
      #if identify controls is specified and the normal filters are set (i.e. not ids, names or experiments_like)
      #then we should probably add the control_feature_types
      
      #The problem with identify_controls at present is that allows control association withint projects (e.g. ENCODE)
      #and does not have sub project specificity (e.g. ENCODE_Yale).
      
      #Does specifying
      
      if(! $self->control_experiments){
                
      } 
    }
    
                    
    if($self->control_experiments  && $self->identify_controls){
      throw("Mutually exclusive parameters have been specified:\t".
        "control_experiments & identify_controls\nPlease omit one of these."); 
      #This prevents redundancy issues between the (unkown) identified and specified controls
      #We generally know which experiments we want to run with specified controls   
    }
    
    
    #if(! ($self->control_experiments || 
    #      $self->identify_controls ||
    #      $self->allow_no_controls) ){
    #  throw("For set_type InputSubset, one of the following must be defined:\t".
    #    'control_experiments identify_controls or allow_no_controls');
    #}
    #Currently this is valid as the controls are loaded as a separate experiment
    #However, this might not always be the case, so it has been disabled
    
    #controls *should* be identified automatically if they are pre-associated
    

    
  return;
}


sub _validate_param_specifity{
  my $self           = shift;
  my $invalid_params = shift;  
  my $err_txt        = '';
  
  foreach my $param_name(@$invalid_params){
    if(my $param = $self->param_silent($param_name)){
      $err_txt .= "\t-".$param_name." => '".$param."'\n";
    }
  }
  
  if($err_txt){    
    throw('You have specified a parameter which is not appropriate for '.$self->set_type."s:\n$err_txt");
  }  
  
  return;
}


sub _fetch_ResultSet_input{
  my $self           = shift;
  $self->_validate_param_specifity(\@iset_only_params);       
  $self->process_params(\@rset_only_params, 1);#1 is optional flag
  return;
}

#We are only using only_replicates for IdentifyReplicateResultSets
 
#Validate there is just 1 supporting signal InputSubset
#and the replicate number matches the TR suffix
#Healthy replicate paranoia here
#All InputSubsets should be replicates
#but historically we have had merged InputSubsets
    

sub _is_replicate_ResultSet{
  my $self = shift;
  my $rset = shift; 
  my $rep_set;
  
  if($rset->table_name eq 'input_subset'){
    my @issets = grep { ! $_->is_control } @{$rset->get_support};
    my $rep    = $issets[0]->replicate;
           
    if((scalar(@issets) == 1) && $rep && 
       ($rset->name =~ /_TR${rep}$/)){
      $rep_set = $rset;
    }
  }
  
  return $rep_set;
}


#todo This needs splitting out into separate sub for each set type (and then some)


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift;
  #Grab dataflow_params here as these params will not 
  #change as we iterate through sets below
  #These will only be flowed to the next job
  #where as batch_params flow across all jobs in this seeded batch
  #(for those that support/require batch params)
  my $dataflow_params = $self->dataflow_params(1);#optional flag
  my $batch_params    = $self->batch_params; 
  my $set_type        = $self->set_type;
  my $only_reps       = ($set_type eq 'ResultSet') ? $self->only_replicates : 0;
  my $set_adaptor     = $self->set_adaptor;
  my $no_write        = $self->param_silent('no_write'); 
  my $sets  = [];
  my $throw = '';
  my (@failed);
  
  ### FETCH SETS ###
 
 
  #Todo
  #1 Changed @failed to build $throw directly with more contextual info
  
  #For set_ids and set_names, catch undef return types
  if($self->set_ids || $self->set_names){
    my ($ids_or_names, $method);
  
    if($self->set_ids){
      $method       = 'fetch_by_dbID'; 
      $ids_or_names = $self->set_ids; 
    }
    else{
      $ids_or_names = $self->set_names; 
      $method = 'fetch_by_name';
    }
    
    foreach my $var(@{$ids_or_names}){
      my $set = $set_adaptor->$method($var);
      
      if(! defined $set ){
        push @failed, $var;
      }
      elsif($only_reps){#is ResultSet
        
        if(! $self->_is_replicate_ResultSet($set)){
          push @failed, $var;    
        }
      }
      elsif(($set_type eq 'InputSubet') && ! $set->replicate){
        #and this is now redundant as all InputSubsets will have a replicate?
        push @failed, $var;
      }
      else{
        push @$sets, $set;
      }
    }
  }
  else{#Must have constraints or experiments
    my $constraints = $self->constraints_hash;
  
    if($self->experiments){  #Preprocess experiment wildcards
      my @exp_names;
      my @exp_wildcards;
      
      for(@{$self->experiments}){
        
        if($_ =~ /%/){
          push @exp_wildcards, $_;
        }
        else{
          push @exp_names, $_;  
        }
      }
      
      my $sql_helper = $set_adaptor->dbc->sql_helper;
      
      foreach my $wc(@exp_wildcards){
        my $sql = 'SELECT name from experiment where name like "'.$wc.'"';
        my @wc_names = @{$sql_helper->execute_simple($sql)};
        
        if(! @wc_names){
          push @failed, $wc;  
        }
        else{
          push @exp_names, @wc_names;  
        }
      }
      
      $constraints->{experiments} =  scalars_to_objects($self->out_db, 
                                                        'Experiment',
                                                        'fetch_by_name',
                                                        \@exp_names);
    }
      
     
    #Need to account for analysis or format or type?
    #i.e. we don't want to queue up the dna meth input_subsets
    #Add string_param_exists here to validate states
    $sets = $set_adaptor->fetch_all( {constraints         => $constraints,
                                      string_param_exists => 1} ); 
    
    #Alternatively we could add a contraint for is_replicate
    #not for ResultSets at present, as they do not have the replicate field
                                        
    if($only_reps){  #Implicit that set type is ResultSet
      
      my @rep_sets;
      
      foreach my $set(@$sets){
        if(my $rep_set = $self->_is_replicate_ResultSet($set)){
          push @rep_sets, $rep_set;  
        }
      }
      
      #Redefine only_replicate sets
      @$sets = @rep_sets;
    }
  }
    
       
  if(@failed){
    $throw = 'Failed to identify some '.$self->set_type." sets. Names or IDs don't exist";
    
    if($only_reps){
      $throw .= ', and/or \'replicates\' were requested but could not find unique signal replicate'; 
    }
    $throw .= ":\n\t".join("n\t", @failed);
  } 
  
  
  ### FETCH/CACHE INPUT_SUBSET CONTROLS ###
  ### SET UP VALIDATION VARS ###
  my ($no_rel_date, $rel_date, $force, $ignore, $rel_month);
  @failed = ();
  my (@failed_ctrls, %ctrl_cache);
  my $x_grp_ctrls = $self->param_silent('allow_inter_group_controls');
  my $use_exp_id;
  
  
  #This can be moved to a _preprocess_InputSubsets methods
  #Then we can 
  
  if($set_type eq 'InputSubset'){
    
    if($self->control_experiments){
      my $exp_adaptor = $self->out_db->get_ExperimentAdaptor;
      my ($exp, $exp_group, $exp_method);
      
      #test whether we have a dbID or a name, assume first is representative of all
      if($self->control_experiments->[0] =~ /^[0-9]+$/){
        $exp_method = 'fetch_by_dbID';
      }
      else{#Assume we have a name
        $exp_method = 'fetch_by_name';
      }    
     
      foreach my $var(@{$self->control_experiments}){
        $exp = $exp_adaptor->$exp_method($var);
        
        if(! defined $exp){
          push @failed, $var;
          next; #$var
        }

        #This should really be Study! Not experimental group. So is not 100% safe
        warn "Cacheing control experiment:\t".$exp->name;
        $exp_group = $exp->get_ExperimentalGroup->name;
        
        $self->_cache_controls(
          $set_adaptor->fetch_all({constraints => {experiments => [$exp]}}),
          \%ctrl_cache, 
          \@failed_ctrls,
          undef,
          $x_grp_ctrls);
        #populates %ctrl_cache and @failed_ctrls
      }
     
      if(@failed){
        $throw .= "\nFailed to fetch some control_experiments. Names or IDs don't exist:\n\t".
                   join("n\t", @failed);   
      } 
      
      if(@failed_ctrls){
        $throw .= join("\n", @failed_ctrls); 
      }    
    }
     
    $use_exp_id = 1 if ! $self->identify_controls;
    #The subtlety here is that we may use a control for an isset which
    #which is not in the control_experiment and not directly associated
    #with that experiment. i.e. we are turning on identify_controls
    #for those experiments which aren't supported by control_experiments  
    #This is why we need to exp_name in the clabel!
    #We can simply test both keys below
    #Which will allow us to use control_experiments or experiment associated controls
    #If they are both defined, then they must match, else we don't know which ones to use!
    #Using associated controls would simply require omitting control_experiments
    #Over-riding associated controls with control_experiments would require deleting them 
    #from the DB
    #Or we could have a control over-ride flag?
    #Which would either use control_experiments preferentially or exclusively?
    
    #The other option is to simply turn on identify_controls, when we specify 
    #control_experiments?
    #This would cause problems when we have the simple case of having controls associated
    #with each experiment
        
    @failed_ctrls = ();
    $self->_cache_controls($sets,
                           \%ctrl_cache,   
                           \@failed_ctrls,
                           $use_exp_id,
                           $x_grp_ctrls);
    if(@failed_ctrls){
      $throw .= join("\n", @failed_ctrls);  
    }
    
    
    $rel_month   = $self->get_param_method('release_month',    'silent');
    $no_rel_date = (defined $rel_month) ? 0 : 1;
    #currently in american format, hence we have a release month for safety
    #rather than a full date string  
    #Don't handle days here as release day is likely to shift, 
    #and we can just manage by reseeding with force_embargoed?
  }

  

  ### BUILD THE OUTPUT ID FOR EACH SET ###
  my (%output_ids, %rep_cache, %embargoed_issets, %to_download_issets, %control_reqd);
  my $tracking_adaptor = $self->tracking_adaptor;
  
  #TODO add more STDOUT if in -no_write mode
  
  #Iterate through the sets 
 SET: foreach my $set( @$sets ){
    $self->helper->debug(1, 'Processing '.ref($set)." set:\t".$set->name.' (Experiment = '.$set->experiment->name.')');
    #"\n\tStates:\t".join(' ', @{$set->get_all_states}));
    
    foreach my $status(qw(REVOKED DISABLED)){ 
      
      if($set->has_status($status)){
        $self->helper->debug(1, "Skipping due to $status");   
        next SET; 
      }
    }
    
    if($set_type eq 'InputSubset'){
      #Horrible arg list for now
      #but we only call it from here
      #This was done to simplify this method
      $throw .= 
        $self->_cache_InputSubset_output_id($set, $tracking_adaptor, \%output_ids, \%ctrl_cache, 
                                            \%rep_cache, \%to_download_issets, \%embargoed_issets, 
                                            $self->force_embargoed, $self->ignore_embargoed, $x_grp_ctrls, 
                                            $use_exp_id, $self->allow_no_controls, $rel_month, 
                                            $no_rel_date, $batch_params, $dataflow_params, 
                                            \%control_reqd);
    }
    else{ #Not an InputSubset (currently just a ResultSet)
     #%output_ids is populated in _cache_ResultSet_output_id
      $self->_cache_ResultSet_output_id($set, \%output_ids, $batch_params, $dataflow_params);
    }
 

    #Need to submit input_subsets them based on shared controls!
    #As we need to wait for the control jobs to download/align 
    #first before we can perform any other analyses
    
    #This is only true if the control has not already been aligned!
    #as we need to semaphore downstream analysis based on success of control alignment
    #this does not however stop us from doing the signal alignments
    #we just don't want to sempahore the whole fan, as we don't want jobs blocking each other
    #only the control job
    
    #WE can't actually do this as this would require semaphore from different analysis
    
    #Can we do this with an accu output from the control?
    #No, as the control job may not finish before the signal jobs dataflow
    #Also, this would create failures of downstream analyses, when really on the 
    #control analysis has failed.
    #This is messy and also pollutes the error output, so real errors would be easily overlooked
  } 
    
  
  #Now iterate through output_ids doing the right thing
  
  my $warn_msg = '';
  
  #split this out into _branch_InputSubset|ResultSet_output_ids
  
 KEY: foreach my $key(keys %output_ids){
    #we need the key to access embargoes/to_download_issets cache 
    my $oid = $output_ids{$key};  
 
    #Also need to check the embargoed status, and skip InputSubsets which are dependant on 
    #an embargoed subset
    #no write to support listing sets before actually seeding them
    #need to be able to run this as Stand alone job, but with access
    #to config. Leo is on the case here. 
    
    if($set_type eq 'InputSubset'){
      #{%$batch_params,
      # input_subset_ids     => {},  #keys are set names or 'controls'
      # %$dataflow_params}; 
 
      my $branch = 2; #ignore or not embargoed 
      
      if(exists $oid->{input_subset_ids}->{controls}){
         #[{$cexp_isset->dbID => $cexp_isset, ...},
         #                       $cexp_isset->feature_type->name,
         #                       $exp_id,
         #                       [],  #embargoed sets
         #                       []   #to_download sets
         #        ];  
        
        if(scalar(@{$oid->{input_subset_ids}{controls}{embargoed}}) > 0){ #CTRLs are embargoed!
          #build log msg appending to throw/warn string and setting branch
          my $ctrls = delete $oid->{input_subset_ids}->{controls};
          my $msg   = "\nFound embargoed control subsets ($key):\t".
                        join(', ', @{$ctrls->{embargoed}})."\nRelated InputSubsets:\n\t".
                        join("\n\t", keys %{$oid->{input_subset_ids}})."\n";
          
          if($self->ignore_embargoed){
            $warn_msg .= "Skipping dataflow for the following InputSubets.\n$msg";   
            $branch = 0;            
            next KEY;
          }
          elsif($self->force_embargoed){
            $warn_msg .= "Forcing dataflow for embargoed controls ($key):\t".
                          join(', ', @{$ctrls->{embargoed}})."\n";
            $oid->{input_subset_ids}->{controls} = $ctrls;   #add ctrls back in!
          }
          else{
            $throw .= $msg;
            $branch = 0;
          }
        }         
      }
    
          
     #todo The more helpful messages about release_month and force/ignore_embargoed are above
     #we should print them once at the bottom
      
      if(exists $embargoed_issets{$key}){ #Some of the signal subsets are embargoed
        
        foreach my $iset_name(values %{$embargoed_issets{$key}}){
        
          if($self->ignore_embargoed){
            $warn_msg .= "\nSkipping dataflow for InputSubset with embargoed signal subsets:\t$iset_name\n";   
            delete $oid->{input_subset_ids}->{$iset_name};
          
          }
          elsif($self->force_embargoed){
            $warn_msg .= "Forcing dataflow for InputSubset with embargoed signal subsets:\t$iset_name\n";
          }
          else{
            $throw .= "Found InputSubset with embargoed signal subsets:\t$iset_name\n";            
            $branch = 0;
          }
        }
      }
        
      #todo Test %to_download_issets, and fan these to $branch to 3 if $branch is true 
      #(i.e. it's not embargoed)
      #Or do we want to download regardless of embargo status?
      #Likely yes, but leave that for now, as figuring out config
      #Also flow to branch 2 aswell after we have fanned
      #Do this by building a job group and using cache_job_group
              
      if($branch){

        if($no_write){
          my $txt   = 'Experiments identified (';
          my $ctrls = '';
          
          if(exists $oid->{input_subset_ids}->{controls}){
            $ctrls = delete $oid->{input_subset_ids}->{controls};
            
            #This currently print dbIDs
            $txt .= 'sharing controls from: '.join(' ', $ctrls->{experiment});
          }
          else{
            $txt .= 'no controls';
          } 

          #Todo add input_subset names to this output

          print STDOUT $txt."):\n\t".join("\n\t", (keys %{$oid->{input_subset_ids}}))."\n";
        }
        elsif(scalar(keys %{$oid->{input_subset_ids}}) > 1 ||
              ((scalar(keys %{$oid->{input_subset_ids}}) == 1) && 
               ! exists $oid->{input_subset_ids}->{controls} ) ){
               
          #Tidy up control cache
          #probably need to revise this when we add branching for download analyses
          if(exists $oid->{input_subset_ids}->{controls}){
            $oid->{input_subset_ids}->{controls} = [keys(%{$oid->{input_subset_ids}{controls}{input_subsets}})];      
          }
            
          $self->branch_job_group($branch, [$oid]);              
        }
        else{#We have no signal subsets to flow!
          $warn_msg .= "No input_subsets to dataflow for control group:\t$key\n";
        }
      }
    }
    else{ #ResultSet
     
      if($no_write){
        print STDOUT "ResultSets Identified\t".$key.":\n\t".
          join(', ', map {$_->{set_name}.'('.$_->{dbID}.')'} @$oid)."\n";
      }
      else{
     
        #Branches now generic wrt set type   
        #2 group of sets (optionally grouped by either controls or parent sets semaphored by 3)
        #>3 fans of those grouped sets, which require different processing
        #actually output id would change dependant on the set type, and the grouping strategy
        #grouping would have to be defined by analysis parameter
        #This will allow expansion of branches without renumbering
        #As branch 2 is always the semaphore, we should always flow this last
        
        #This does however mean that for dataflow with no semaphore we will never define
        #a funnel, and therefore only flow down branch 3 and not branch 2.
        
        my @job_group = (3, $oid);
        
        if($key ne 'ungrouped'){
          #Now define funnel job id with all set details i.e. RunIDR
          
          my $funnel_oid = {%$batch_params,
                            %$dataflow_params,
                            dbIDs     => [], 
                            set_names => [],
                            set_type  => 'ResultSet'}; 
                            
          foreach my $id(@$oid){
            push @{$funnel_oid->{dbIDs}},     $id->{dbID};
            push @{$funnel_oid->{set_names}}, $id->{set_name};
          }
          
          push @job_group, (2, [$funnel_oid]);
        }
        
        $self->branch_job_group(@job_group);  
      }
    }
  } #end foreach $oid
    
  #Only throw here rather than above so we get all the information
  #about what has failed  
  
  warn $warn_msg if $warn_msg;
  $self->throw_no_retry($throw) if $throw; #we have not actually dataflowed yet  
  return;
}


#These have been separated to enable digestion/review of caller method
#but are only ever called from run in 1 place.

sub _cache_ResultSet_output_id{
  my ($self, $rset, $oids, $batch_params, $dataflow_params) = @_; 
  my $group_key = 'ungrouped';
   
  if($self->only_replicates){
    #test is replicate here?
    #probably already tested previously     
    
    #just parent ResultSet name (which doesn't exist yet)
    #This is a little fragile as it depends on our nomenclature
    #would replace this with the get set prefix
    #then at least we are using the same code everywhere
    #($group_key = $set->name) =~ s/_TR[0-9]+$//o;
    
    $group_key = get_set_prefix_from_Set($rset).'_'.$rset->analysis->log_name;    
    #Appened analysis name here, just in case we are running
    #multiple analyses for the same set of data.
  }
  

  #if( ! exists $oids->{$group_key}){
  #  $oids->{$group_key} = {%$batch_params,
  #                         %$dataflow_params,
  #                         dbIDs => [], 
  #                         set_names => []};
  #} 

  # Is this right? We defo need some grouping here
  # For control conversion?
  # No, we always have a single id for for each job
  # So what is the group key for here?

  $oids->{$group_key} ||= [];
  push @{$oids->{$group_key}}, {%$batch_params,
                                %$dataflow_params,
                                dbID     => $rset->dbID, 
                                set_name => $rset->name,
                                set_type => 'ResultSet'}; 
  return;
}


sub _cache_InputSubset_output_id{
  my ($self, $set, $tdb_adaptor, $output_ids, $ctrl_cache, $rep_cache,
      $to_download_issets, $embargoed_issets, $force_embargoed, 
      $ignore_embargoed, $x_grp_ctrls, $use_exp_id, $allow_no_ctrls,
      $rel_month, $no_rel_date, $batch_params, $df_params, $control_reqd) = @_;
  my ($to_download, $embargoed);
  my $throw = '';#To avoid concat undef warning in caller
  #Don't return when first encountering throw, as we want the full details
  #not just the first details.

  $tdb_adaptor->fetch_tracking_info($set); #Injects tracking methods

  if( $tdb_adaptor->is_InputSubset_embargoed($set, $rel_month)){ 
 
    if(! ($force_embargoed || $ignore_embargoed)){
      my $rd_txt = '';
 
      if($no_rel_date){
        $rd_txt = "\nOr maybe you want to specify a release_month when seeding this analysis?"; 
      }
   
      $throw .= "\nFound InputSubset which is not out of embargo:\n\t".$set->name.
        "\nYou can over-ride this by specifying force_embargo or ignore_embargo.".$rd_txt;
    }
    elsif($ignore_embargoed){
      $embargoed = 1;
    }
    else{ return $throw; }
  }
      
  
  if(! defined $set->download_date){
    $to_download = 1;
    $throw .= "\nFound InputSubsets which are not downloaded:\n\t".$set->name.'('.$set->experiment->name.")\n".
      'Please specify -skip_absent, -download or if a bulk download is required, '.
      'run DOWNLOAD_SCRIPT separately?' if ! $self->skip_absent;               
  }
  elsif(! -e $set->local_url){
    $to_download = 1;
    $throw .= "\nFound InputSubsets which were downloaded on ".$set->download_date." but are not present:\n\t".
      $set->name.'('.$set->experiment->name.")\t".$set->local_url."\n".
      'Please specify -skip_absent, -download or if a bulk download is required, '.
      'run DOWNLOAD_SCRIPT separately?' if ! $self->skip_absent;       
  }
       
  #Cache embargoed/to_download ctrls
     
  #Is this too complicated?
  #The aim here is to allow submission of sets between groups
  #some of which have controls in the same group, and some which need to use control
  #from another group
      
  #The issue here is that we may identify > 1 control
  #and we want to take the one from the matching exp group first
      
  #We should probably just drop x_group_controls completely
  #and let control_experiments over-ride that?  
  #Or only use x_group_controls here should not take
      
  if($set->is_control && 
     ($embargoed || $to_download)){ #Set embargoed/to_download in the ctrl cache
    #Recreate the cache key logic here
    my $clabel = $set->cell_type->name;
  
    if(! $x_grp_ctrls){
      $clabel .= '_'.$set->experiment->experimental_group->dbID;
    }
      
    #double exists so we don't auto vivify 
    if(exists $ctrl_cache->{$clabel} &&
       exists $ctrl_cache->{$clabel}->{input_subsets}->{$set->dbID}){
            
      if($embargoed){
        $ctrl_cache->{$clabel}->{embargoed} ||= [];
        push @{$ctrl_cache->{$clabel}->{embargoed}}, $set->dbID;
      }
            
      if($to_download){
        $ctrl_cache->{$clabel}->{to_download} ||= [];
        push @{$ctrl_cache->{$clabel}->{to_download}}, $set->dbID;
      }
    }
            
                
    if($use_exp_id){
      $clabel = $set->cell_type->name.'_'.$set->experiment->dbID;
         
      if(exists $ctrl_cache->{$clabel} &&
         exists $ctrl_cache->{$clabel}->{input_subsets}->{$set->dbID}){
        #double exists so we don't auto vivify 
            
        if($embargoed){
          $ctrl_cache->{$clabel}->{embargoed} ||= [];
          push @{$ctrl_cache->{$clabel}->{embargoed}}, $set->dbID;
        }
            
        if($to_download){
          $ctrl_cache->{$clabel}->{to_download} ||= [];
          push @{$ctrl_cache->{$clabel}->{to_download}}, $set->dbID;
        }
      }
    }
  }
      
   
    
  if(! $set->is_control){    
    #We already have the ctrls in the ctrl_cache        
    #automatically use control attached to the same experiment
    #how does this interact with control_experiments?
    #This should not cache as we don't want to identify_controls between experiments
    #too much? just auto identify controls?
    my $clabel   = $set->cell_type->name;
      
    if(! $x_grp_ctrls){
      $clabel .= '_expgrp_'.$set->experiment->experimental_group->dbID;
    }
        
    $self->helper->debug(3, 'Attempting to identify control for '.$set->name.' using control key '.$clabel);
    my $non_exp_ctrl = $ctrl_cache->{$clabel} if exists $ctrl_cache->{$clabel};
    my $exp_ctrl;
        
    if($use_exp_id){
      $clabel   = $set->cell_type->name.'_exp_'.$set->experiment->dbID;
      $self->helper->debug(3, 'Attempting to identify control for '.$set->name.' using control key '.$clabel);
      $exp_ctrl = $ctrl_cache->{$clabel} if exists $ctrl_cache->{$clabel};
    }
        
                
    if( (! ($non_exp_ctrl || $exp_ctrl)) &&
       ! $allow_no_ctrls){
      $throw .= "Could not find a control for InputSubset:\n\t".$set->name.'('.
        $set->experiment->name.")\n";     
    }
    elsif($non_exp_ctrl && $exp_ctrl){
          
      if($non_exp_ctrl->{experiment} ne $exp_ctrl->{experiment}){   #Make sure these are the same!
        $throw .= "Found control clash between experiments:\t".
          $non_exp_ctrl->{experiment}.' & '.$exp_ctrl->{experiment};
      }
        
      #Can we assume these are identical?
      #They might not be control_experiments will get all the control from the experiments
      #but controls from the sets we have fetched (either identify_controls or just associated controls)
      #might not represent all the controls from an experiment
      #In this case, we can take the control_experiments
      #As this won't break any inter/intra experimental_group rules
      #Just do this below
    }
    else{  #Preferentially take possibly more complete control_experiment derived controls
      $non_exp_ctrl ||= $exp_ctrl; 
      #Now we need to group these controls with their signal sets
      #as we need to flow these in a batch  
      my $ctrl_cache_key;
        
      if(! $non_exp_ctrl){
        $ctrl_cache_key = 'no_controls';
             
        ### VALIDATE feature_set_analysis doesn't needs control (regardless of allow_no_controls )###
        #Here rather than in DefineResultSets to catch it earlier and prevent 
        #doomed job/input_id creation which may clash with future jobs submission and 
        #hence require some sort of recover/rollback flag
             
        #Do we need to be mindful of result_set_only here?
        #still don't want to create erroneous ResultSet even if we aren't running
        #the peak analysis. But wiggle track really doesn't require control
        #This is a small corner case, just bail out here for now.
        #But if we wanted to allow it would also need to test in DefineMergedOutputSet
        #Would also need to support rollback ResultSet to import with updated controls
        #this will result in collection being regenreated unecessarily :/
             
        ### GET FEATURE SET ANALYSIS ###
        my $fset_anal = $self->feature_set_analysis;
    
        if(! defined $fset_anal){
             
          if(exists $self->default_feature_set_analyses->{$set->feature_type->name}){
            $fset_anal = $self->default_feature_set_analyses->{$set->feature_type->name}; 
          }
          elsif(exists $self->default_feature_set_analyses->{$set->feature_type->class}){
            $fset_anal = $self->default_feature_set_analyses->{$set->feature_type->class};
          }
          else{
            return $throw."No default feature_set analysis available for ".$set->name.
              ' ('.$set->feature_type->name.' or '.$set->feature_type->class.").\n Please add FeatureType name or class to  default_feature_set_analyses ".
              "in BaseSequencingAnalysis::default_options or specify -feature_set_analysis\n"; 
          }
        }
              
        ### TEST WHETHER IT NEEDS CONTROLS ###          
        if(! exists $control_reqd->{$fset_anal}){
          eval { $fset_anal = scalars_to_objects($self->out_db, 
                                                 'Analysis', 
                                                 'fetch_by_logic_name', 
                                                 [$fset_anal])->[0]; };
         
          if( ! $fset_anal){
            return $throw.'Could not get analysis for '.$set->name."\n$@\n";  
          }
          
          my $peak_module;  
          eval { $peak_module = validate_package_path($fset_anal->module); };
              
          if(! $peak_module){
            return $throw.'Failed to validate analysis module for InputSubset '.$set->name."\n$@\n"; 
          }
                
          $control_reqd->{$fset_anal} = $peak_module->requires_control;
        }
            
        if($control_reqd->{$fset_anal}){
          return $throw.$set->name.' requires controls for '.$fset_anal->logic_name.
            ", but none have been identified\n";
        }
      }
      else{   #$non_exp_ctrl
        #This was the experiment dbID, but have changed this to the experiment name for 
        #readability in the logs 
        $ctrl_cache_key = $non_exp_ctrl->{experiment};
      }
      
    
       
      ### BUILD OUTPUT ID ###
      if(! exists $output_ids->{$ctrl_cache_key}){
        #Handle ignore_embargoed controls, and therefore all depedant subsets 
        $output_ids->{$ctrl_cache_key} = {%$batch_params,
                                          input_subset_ids     => {},  #keys are set names or 'controls'
                                          %$df_params}; 
        #Add ctrls first
        if($ctrl_cache_key ne 'no_controls'){
          #This should get updated with embargoed/to_download status
          #of subsequent ctrl sets which haven't yet been seen
          $output_ids->{$ctrl_cache_key}->{input_subset_ids}->{controls} = $non_exp_ctrl;  
        }       
      }
  
      #Now add sig reps based on the exp/cell_type/feature_type input_subset_ids
      #Currently exp ID is a synecdoche for this, but let's future proof
      #in case we do move to proper experiment definition i.e. a study with
      #>1 ctype/ftypes           
      #my $iset_cache_key = $set->experiment->dbID.'_'.
      #                       $set->cell_type->dbID.'_'.
      #                       $set->feature_type->dbID;         
      #Cache key is essentially the set name
      #although this is more expensive to generate than a key just based on IDs
           
      #Define set name
      #This prefix should always be the standard (i.e. not the control) prefix
      my $set_name = get_set_prefix_from_Set($set);        
    
      #Would be best to identify redundant replicates here
      #so we don't have to iterate through the cache again
      #This result in unordered/grouped throw messages
      #as InputSubsets may be processed in any order
      #IDing them here using a cache will also prevent having to
      #cache subset objects.
      #This will mean we will need to pass another cache ref.
      #Currently don't do this in _cache_controls
      #as we allow redundant replicate numbers for now
      #Although this will have to change when we change to aligning 
      #all replicate separately, then merging replicate alignments (without rmdups)
    
      #Don't care about autovivifying here
      if(defined $rep_cache->{$set_name}{$set->replicate}){
      
        if($rep_cache->{$set_name}{$set->replicate} ne $set->name){
          return $throw."Found redundant records for $set_name replicate ".$set->replicate.":\n\t".
            $set->name."\n\t".$rep_cache->{$set_name}{$set->replicate}."\n";
        }
        else{ #This should never happen
          $self->throw_no_retry("Found redundant set details for InputSubset ".$set->name); 
        }       
      }
      else{
        $rep_cache->{$set_name}{$set->replicate} = $set->name;    
      } 
           
      if($embargoed){
        $embargoed_issets->{$ctrl_cache_key}{$set_name} ||= [];
        push @{$embargoed_issets->{$ctrl_cache_key}{$set_name}}, $set->dbID;
      }
           
      if($to_download){
        #do we need to add dbids and replicate keys in here?
        $to_download_issets->{$ctrl_cache_key}{$set_name} ||= [];
        push @{$to_download_issets->{$ctrl_cache_key}{$set_name}}, $set->dbID;  
      }
                      
      $output_ids->{$ctrl_cache_key}->{input_subset_ids}->{$set_name} ||= [];
      push @{$output_ids->{$ctrl_cache_key}->{input_subset_ids}->{$set_name}}, $set->dbID;
    }        
  }
  
  
  #These are all references so are already updated in the caller     
  #return ($output_ids, $ctrl_cache, $to_download_issets, $embargoed_issets);  
  return $throw;     
}





#This currently won't pick any controls for ENCODE exps they are stored as separate experiments
#Maybe the registration code will make these assoications
#but we need to be able to pick them up here to
#so we will need to allow cross group controls
#but be explicit about the experiments and controls_experiments
#in the input_id. Otherwise we risk picking up controls from other experimental sub groups

#When using cross group controls, this checks to make sure we are only getting them form one exp
#Let's change this to use names, then we can use this as the ctrl_cache_key above

#Is the $clabel key ever used outside of this sub? Yes, it is used above to actually
#populate the embargoed and to download sets elements below

#todo Change this cache value data structure to a hash to avoid access by indexes in the rest of the code

sub _cache_controls{
  my ($self, $cexp_issets, $ctrl_cache, $failed, $use_exp_id, $x_grp_ctrls) = @_;    
  my $clabel; 
   
  foreach my $cexp_isset(@{$cexp_issets}){
       
    if($cexp_isset->is_control){   
      $clabel    = $cexp_isset->cell_type->name;
      
      #Can we change this to u
      
      my $exp_id = $cexp_isset->experiment->dbID;
    
      if($use_exp_id){
        $clabel .= '_exp_'.$exp_id;
      }
      elsif(! $x_grp_ctrls){
        $clabel .= '_expgrp_'.$cexp_isset->experiment->experimental_group->dbID;
      }
    
      if(! exists $ctrl_cache->{$clabel}){
        #Hash required here to avoid redundant controls
        #$ctrl_cache->{$clabel} = [{$cexp_isset->dbID => $cexp_isset},
        #                          $cexp_isset->feature_type->name, #Used for validation
        #                          $cexp_isset->experiment->name, #Used for validation
        #                          [],  #embargoed sets
        #                          []   #to_download sets
        #                          ];  
                                  
       $ctrl_cache->{$clabel} = { input_subsets => {$cexp_isset->dbID => $cexp_isset},
                                  feature_type  => $cexp_isset->feature_type->name, #Used for validation
                                  experiment    => $cexp_isset->experiment->name, #Used for validation
                                  embargoed     => [], 
                                  to_download   => [] };
      }
      elsif( ($cexp_isset->feature_type->name eq $ctrl_cache->{$clabel}{feature_type}) &&
             ($cexp_isset->experiment->name eq $ctrl_cache->{$clabel}{experiment}) ){
        $ctrl_cache->{$clabel}{input_subsets}{$cexp_isset->dbID} = $cexp_isset;
      }
      else{
        #May need to improve this message
        push @$failed, $clabel.' control '.$cexp_isset->name." clashes with cache:\n\t".
          $ctrl_cache->{$clabel}{feature_type}."\n\t".$ctrl_cache->{$clabel}{experiment};    
      }
    }#fi is_control
    
  }#end foreach $isset
  
  return;# $failed; 
}




#Enable a preview of what we are going to dataflow here
#This will be done using standaloneJob.pl once(if ever) this can access config from the DB
#will probably have to add a no_data_flow flag, which will print out instead of
#flow the output_ids
#Currently -no_write does this, but means we have extra IdentifyJobs in the pipe DB (both DONE and FAILED)

sub write_output {  # Create the relevant jobs
  my $self = shift;  
  $self->dataflow_job_groups;
  return;
}



1;
