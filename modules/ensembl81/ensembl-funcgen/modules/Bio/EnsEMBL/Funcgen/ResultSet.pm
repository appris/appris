#
# Ensembl module for Bio::EnsEMBL::Funcgen::ResultSet
#

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

Bio::EnsEMBL::Funcgen::ResultSet - Represents a set of experimental results (signal).


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ResultSet;

my $result_set = Bio::EnsEMBL::Funcgen::ResultSet->new
  (
   -dbid        => $dbid,
   -analysis    => $analysis,
   -table_name  => 'input_subset',
   -table_id    => $input_subset_id,
   -feature_class => 'dna_methylation',
  );



=head1 DESCRIPTION

A ResultSet object provides access to a set raw or processed signal values from an Experiment,
and can consist of any number of inputs which can either be InputSets or ExperimentalChips.
A ResultSet can represent a single or multiple merged replicates of raw or processed(normalised)
signal data.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::Set

=cut

package Bio::EnsEMBL::Funcgen::ResultSet;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use base qw( Bio::EnsEMBL::Funcgen::Set Bio::EnsEMBL::Funcgen::feature_class_Set);

#Valid enum field hashes to prevent loading NULLs
#result_set_input.table_name
my %valid_table_names =
  (
   experimental_chip => undef,
   input_subset      => undef,
   channel           => undef,
  );


=head2 new

  Arg [-ANALYSIS] :



  Example    : my $feature = Bio::EnsEMBL::Funcgen::ResultSet->new
                 (
                  -dbid               => $dbid,
                  -analysis           => $analysis,
                  -support            => \@support_objects,
                  -feature_class      => 'result',
                  -feature_type       => $ftype,
	             );
  Description: Constructor for ResultSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : Throws if feature class is not valid.
               Throws if table name not defined.
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($table_name, $table_id, $dbfile_path, $support, $rep)
    = rearrange(['TABLE_NAME', 'TABLE_ID', 'DBFILE_PATH', 'SUPPORT', 'REPLICATE'], @_);
  # TEST MANDATORY PARAMS
  # explicit type check here to avoid invalid types being imported as NULL
  # and subsequently throwing errors on retrieval
  $self->_validate_feature_class(\@_);

  # set default type until this is moved to db_file_registry.format
  # This is not possible yet as 5mC is classed as DNA not DNA Modification!!!
  $self->{replicate}      = $rep;
  $self->{table_ids}      = {};


  #Need to maintain -table_name and _add_table_id for lazy loading ability
  if(defined $table_name){
     if (! (defined $table_name &&
           exists $valid_table_names{$table_name}) ){
       throw('You need to pass a valid -table_name e.g. '.
            join(', ', keys %valid_table_names));
    }

    $self->{table_name}         = $table_name;
  }
  #else defined based on support?

  #This will currently allow a ResultSet with no support/table_ids
  if(defined $support){
    $self->add_support($support);
  }
  elsif(! defined $table_name){
    throw('You must provide either a -support or a -table_name parameter');
  }

  #Remove this when -table_id fully deprecated
  if(defined $table_id){
    deprecate('The -table_id param is now deprecated in new. Please use -support');
    $self->_add_table_id($table_id) if defined $table_id;

    if(defined $support){
      throw('Unsafe to specify -support and -table_id, please use -support');
    }
  }

  # Don't use setter here as this will call $self->adaptor->dbfile_data_root
  # which may not have been set yet.
  $self->{dbfile_path} = $dbfile_path if defined $dbfile_path;

  return $self;
}


sub _valid_feature_classes{
  return qw( result dna_methylation );
}


=head2 reset_relational_attributes

  Arg[1] : Hashref containing the following mandatory parameters:
            -analysis     => Bio::EnsEMBL::Analysis,
            -feature_type => Bio::EnsEMBL::Funcgen::FeatureType,
            -cell_type    => Bio::EnsEMBL::Funcgen::CellType,
            -support      => Arrayref of valid support objectse.g. InputSet

  Description: Resets all the relational attributes of a given ResultSet.
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;
  my ($support, $analysis, $feature_type, $cell_type, $exp) =
    rearrange(['SUPPORT', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE', 'EXPERIMENT'], %$params_hash);

  #flush table ID cache and add support
  $self->{table_ids} = undef;
  $self->{support}   = undef;
  $self->add_support($support);

  #is_stored (in corresponding db) checks will be done in store method

  assert_ref($analysis,     'Bio::EnsEMBL::Analysis');
  assert_ref($feature_type, 'Bio::EnsEMBL::Funcgen::FeatureType');
  assert_ref($cell_type,    'Bio::EnsEMBL::Funcgen::CellType');

  if( $self->experiment &&
      ! check_ref($exp, 'Bio::EnsEMBL::Funcgen::Experiment') ){
    throw("You must pass a valid Bio::EnsEMBL::Funcgen::Experiment\n".
          "Passed:\t".ref($exp));
  }

  if(defined $exp){
    #This will allow addition of an experiment
    #when it was not prevously defined
    $self->{experiment_id} = $exp->dbID;
    $self->{experiment}    = $exp;
  }


  $self->{cell_type}    = $cell_type;
  $self->{feature_type} = $feature_type;
  $self->{analysis}     = $analysis;

  #Finally undef the dbID and adaptor by default
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  return;
}


=head2 replicate

  Description: Accessor for replicate attribute
  Returntype : Int
  Exceptions : None
  Caller     : Genereal
  Status     : At risk

=cut

sub replicate{
  return shift->{replicate};
}

sub add_support{
  my($self, $support) = @_;

  if(! (defined $support &&
        (ref($support) eq 'ARRAY')) ){
    throw('You must pass an Arrayref of support objects');
  }

  #Check passed table name matches objects
  #remove this once we remove add_table_ids
  #although we still need to check all support object are the same
  my $tmp_table_name = $self->table_name;
  my $sdba = $support->[0]->adaptor;
  (my $class_name = ref($support->[0])) =~ s/.*://g;

  if(! $sdba){
    throw("No $class_name adaptor found. All ResultSet -support must be stored");
  }

  my $table_name = $sdba->_main_table->[0];

  if(defined $tmp_table_name &&
    ( $tmp_table_name ne $table_name) ){
     throw('Specified -table_name does not match -support table name. '.
          'Please omit -table_name and ensure all -support objects are of the same class');
  }
  elsif(! defined $tmp_table_name){
    $self->{table_name} = $table_name;

    #Validate table name
    if(! exists $valid_table_names{$table_name}){
      throw("The -table_name or -support objects do not refer to a valid".
            " support table:\t$table_name\nValid table names:\t".
            join(', ', keys %valid_table_names));
    }
  }

  #todo validate $ssets have same cell_type and feature_type

  #Add the table_ids and set the support
  foreach my $sset(@$support){

    if (! ( defined $sset &&
            ref($sset) &&
            $sset->isa('Bio::EnsEMBL::Funcgen::'.$class_name)) ){
         throw(ref($sset)." is not a valid support type for this $class_name ResultSet\n");
    }

    #we don't have access to is_stored_and_valid here
    #so this may throw an error

    #if the table id is already defined, it has either been set by
    #obj_from_sth or add_support
    #For both these cases we want to throw from here

    if(exists $self->{table_ids}->{$sset->dbID}){
      throw('Cannot add_support for previously added/stored ResultSet support: '.
        ref($sset).'('.$sset->dbID.')');
    }

    $self->_add_table_id($sset->dbID);
  }


  $self->{support} ||= [];
  push @{$self->{support}}, @$support;

  return $self->{support};
}


=head2 display_label

  Example    : print $set->display_label;
  Description: Getter for the display_label attribute for this Set.
  Returntype : String
  Exceptions : Warns if feature_class not explicitly handled
  Caller     : General
  Status     : At Risk

=cut

#Add project name in here
#Where is this used in the webcode? Track names?
#Not for ChIP-Seq DataSets as these use the merged class descriptions
#e.g. Open Chromatin & TFBS
#Just for dna_methylation, as these are always show on one track.

sub display_label {
  my $self = shift;

  if(! defined $self->{display_label}){

    if($self->feature_class eq 'result'){ # We have a ChIP/DNase1/FAIRE signal set
      $self->{display_label} = $self->feature_type->name.' '.$self->cell_type->name.' signal';
    }
    elsif($self->feature_class eq 'dna_methylation'){

      my $project = $self->project || '';

      if($project){
        $project = ' '.$project;
      }

      $self->{display_label} = $self->analysis->display_label.' '.$self->cell_type->name.$project;
    }
    else{
      warn 'ResultSet feature class('.$self->feature_class.') not explicitly handled in display_label';
      $self->{display_label} = $self->feature_type->name.' '.$self->cell_type->name;
    }
  }

  return $self->{display_label};
}


=head2 dbfile_path

  Arg[1]     : String (optional) - File path for flat 'DB' file for this ResultSet
                e.g. a BigWig/Bed file
  Example    : my $dbfile_path = $self->dbfile_path;
  Description: Getter/Setter for the dbfile path attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut

sub dbfile_path{
  my $self             = shift;
  $self->{dbfile_path} = shift if @_;
  return $self->{dbfile_path};
}


=head2 table_name

  Example    : if($result_set->table_name eq 'input_subset'){
                 #Do something here
               }
  Description: Getter for the table_name of a ResultSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub table_name{ return shift->{table_name}; }


=head2 _add_table_id

  Example    : $result_set->add_table_id($ec_id, $cc_id);
  Description: Caches table_id result_set_input_id to the ResultSet. In the case of an
               array ResultSet, the unique result_set_input_id is used to key into the
               result table, it also reduces redundancy and enable mapping of results to chips
               rather than just the ResultSet.  This enables result retrieval
               based on chips in the same set which  have a differing status.
               In the case of a sequencing ResultSet, this simply refers to the InputSet or InputSubset ids.
  Returntype : None
  Exceptions : Throws if no table_id defined
  Caller     : General
  Status     : At Risk

=cut

sub _add_table_id {
  my ($self, $table_id, $cc_id) = @_;

  if (! defined $table_id){
    throw('Need to pass a table_id');
  }else{

    #This allows setting of the cc_id on store
    if((exists $self->{table_ids}->{$table_id}) &&
       (defined $self->{table_ids}->{$table_id})){
      throw("You are attempting to redefine a result_set_input_id which is already defined");
    }

    $self->{table_ids}->{$table_id} = $cc_id;
  }

  return;
}


=head2 table_ids

  Example    : $result_set->feature_group_id($fg_id);
  Description: Getter and setter for the feature_group_id for this ResultSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub table_ids { return [ keys %{shift->{table_ids}} ]; }


=head2 result_set_input_ids

  Example    : my @rset_rsi_ids = @{$result_set->result_set_input_ids()};
  Description: Getter for the input ids for this ResultSet.
  Returntype : arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub result_set_input_ids { return [ values %{shift->{table_ids}} ]; }


=head2 contains

  Example    : if($result_set->contains($chip_or_channel)){...do some chip or channel erpartions here...};
  Description: Returns true if the given Channel or ExperimentalChip is part of this ResultSet
  Returntype : boolean
  Exceptions : warns if ResultSet table name is not of argument type
  Caller     : General
  Status     : At Risk

=cut

sub contains{
  my ($self, $chip_channel) = @_;
  my $contains = 0;
  my $table_name = $chip_channel->adaptor->_main_table->[0];

  if($table_name ne $self->table_name){
    warn("ResultSet(".$self->table_name().") cannot contain ${table_name}s");
  }else{
    $contains = 1 if exists $self->{table_ids}->{$chip_channel->dbID()};
  }

  return $contains;
}


=head2 get_result_set_input_id

  Arg [1]    : int - dbID (experimental_chip, channel or input_subset)
  Example    : $result_set->get_result_set_input_id($ec_id);
  Description: Retrieves a result_set_input_id from the cache given a dbID
  Returntype : int
  Exceptions : none
  Caller     : General
  Status     : At Risk

=cut

sub get_result_set_input_id{
  my ($self, $table_id) = @_;
  return (exists $self->{table_ids}->{$table_id}) ?  $self->{table_ids}->{$table_id} : undef;
}


=head2 get_support

  Example    : my @support = @{$result_set->get_support};
  Description: Retrieves a list of objects supporting the ResultSet
  Returntype : Arrayref of objects e.g. InputSet, InputSubset, ExperimentalChip or Channel
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_support {
  my $self = shift;

  if(! defined $self->{support}){
    my $adaptor_method = 'get_'.
                         join('',  (map ucfirst($_), split(/_/, $self->table_name))). #export from EFGUtils?
                         'Adaptor';

    #Adaptor may be absent if we have an unstored object and not used -support in new
    #This will disappear once -table_id and -table_name are removed from new
    my $supp_adaptor = $self->adaptor->db->$adaptor_method;

    foreach my $id(@{$self->table_ids()}){
      push @{$self->{support}}, $supp_adaptor->fetch_by_dbID($id);
    }
  }

  return $self->{support};
}




=head2 get_ExperimentalChips

  Example    : my @ecs = @{$result_set->get_ExperimentalChips()};
  Description: Retrieves a list of ExperimentalChip objects
  Returntype : Listref of ExperimentalChip object
  Exceptions : Warns if ResultSet is not associated to any ExperimentalChips
  Caller     : General
  Status     : At Risk

=cut

sub get_ExperimentalChips{
  my $self = shift;
  my $echips;

  if( $self->table_name ne 'experimental_chip' ||
      $self->table_name ne 'channel' ){
    warn 'Cannot get_ExperimentalChips for an ResultSet with table_name '.$self->table_name;
  }
  elsif($self->table_name eq 'experimental_chip'){
    $echips = $self->get_support;
  }
  else{  #table_name is channel

    if(! defined $self->{experimental_chips}){
      my %echips;
      my $chan_adaptor = $self->adaptor->db->get_ChannelAdaptor;
      my $ec_adaptor   = $self->adaptor->db->get_ExperimentalChipAdaptor;

      foreach my $chan_id( @{$self->table_ids} ){
        my $chan = $chan_adaptor->fetch_by_dbID($chan_id);
        $echips{$chan->experimental_chip_id} ||= $ec_adaptor->fetch_by_dbID($chan->experimental_chip_id);
      }

      @{$self->{experimental_chips}} = values %echips;
    }

    $echips = $self->{experimental_chips};
  }

  return $echips;
}



=head2 get_replicate_set_by_result_set_input_id

  Arg[0]     : int - chip_channel_id
  Example    : my $rep_set_name = $result_set->get_replicate_set_by_result_set_input_id($cc_id);
  Description: Retrieves the replicate set name defined by the corresponding ExperimentalChip
  Returntype : String - replicate set name
  Exceptions :
  Caller     : General
  Status     : At Risk - implement for Channels?

=cut

#Where is this used?

sub get_replicate_set_by_result_set_input_id{
  my ($self, $cc_id) = @_;

  if( ! defined $self->{_replicate_cache}){

    foreach my $ec (@{$self->get_ExperimentalChips()}){
      $self->{_replicate_cache}{$self->get_result_set_input_id($ec->dbID())} = $ec->replicate();
    }
  }

  #warn here of absent replicate info?
  return (exists $self->{_replicate_cache}{$cc_id}) ?  $self->{_replicate_cache}{$cc_id} : undef;
}

sub get_replicate_set_by_chip_channel_id{
  my ($self, $cc_id) = @_;

  deprecate('Please use get_replicate_set_by_result_set_input_id instead');
  return $self->get_replicate_set_by_result_set_input_id($cc_id);
}




=head2 log_label

  Example    : print $rset->log_label();
  Description: Get a string of the unique key fields for logging purposes
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub log_label {
  my $self = shift;

  my $label;

  if(defined $self->feature_type){
    $label = $self->feature_type->name.":";
  }else{
    $label = "Unknown FeatureType:";
  }

  if(defined $self->cell_type){
    $label .= $self->cell_type->name;
  }else{
    $label .= "Unknown CellType";
  }

  return $self->name.":".$self->analysis->logic_name.":".$label;
}


=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of ResultSet method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: name table_name feature_class get_all_states
  Args[4]    : Arrayref - Optional list of ResultSet method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: feature_type cell_type analysis get_support
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this ResultSet to another based on the defined scalar
               and storable methods.
  Returntype : Hashref of key attribute/method name keys and values which differ.
               Keys will always be the method which has been compared.
               Values can either be a error string, a hashref of diffs from a
               nested object, or an arrayref of error strings or hashrefs where
               a particular method returns more than one object.
  Exceptions : None
  Caller     : Import/migration pipeline
  Status     : At Risk

=cut

sub compare_to {
  my ($self, $obj, $shallow, $scl_methods, $obj_methods) = @_;
  $obj_methods ||= [qw(feature_type cell_type analysis get_support)];
  $scl_methods ||= [qw(name table_name feature_class get_all_states)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods, $obj_methods);
}


=head2 experiment

  Arg[1]     : Boolean - Control flag, get's control experiment instead of signal experiment
  Example    : my $exp = $result_set->experiment;
  Description: Getter for the Experiment of this ResultSet.
               Returns undef if there is more than 1 contributing Experiment
  Returntype : Bio::EnsEMBL::Funcgen::Experiment or undef
  Exceptions : Warns if cannot identify unique Experiment
  Caller     : General
  Status     : At Risk

=cut

#experiment_id is not mandatory, so we can still get this from the support
#probably best to make it mandatory in the near future, although this won't support control
#use it if it is there, else use rsi's


sub experiment{
  my $self    = shift;
  my $control = shift;
  my $attr_name = 'experiment';
  $attr_name    = 'control_'.$attr_name if $control;

  my $exp;
  if (! exists $self->{$attr_name}){ #exists as undef is valid

    if((! $control) &&
      $self->experiment_id){
      $exp = $self->adaptor->db->get_ExperimentAdaptor->fetch_by_dbID($self->experiment_id);
    }
    else{
      #These are likely InputSubsets, but could be others e.g. InputSets, ExperimentalChips etc
      my @support = @{$self->get_support};

      foreach my $support(@support){

        if($support->can('is_control')){

          if(($support->is_control && ! $control) ||
             ($control && ! $support->is_control)){
            next;
          }
        }
        elsif($control){
          throw("Cannot get control Experiment for ResultSet with support type:\t".ref($support));
        }

        $exp ||= $support->experiment;

        if($support->experiment->dbID != $exp->dbID){
          undef $exp;
          warn('Failed to get unique Experiment for ResultSet '.$self->name."\n");
          last;
        }
      }
    }

    $self->{$attr_name} = $exp;
  }

  return $self->{$attr_name};
}


### Stay of execution  ###
# These methods were to be deprecated in v81
# But due to limitations with the GRCh37 site, they must remain until the next 
# data update of that site i.e. entirely new mouse DB with bigwigs and convert
# existing human to bigwigs


# Is this even required anymore? Isn't it based on RESULT_FEATURE_SET status?


sub get_ResultFeatures_by_Slice{  # DEPRECATED in v81
  #deprecate('ResultSet not longer supports ResultFeature retrieval. Access the underlying flat file by using dbfile_path and the appropriate library e.g. Bio::DB::BigFile');
  my ($self, $slice, $max_bins, $window_size, $constraint) = @_;
  return $self->adaptor->db->get_ResultFeatureAdaptor->fetch_all_by_Slice_ResultSet($slice, $self, $max_bins, $window_size, $constraint);
}


sub get_dbfile_path_by_window_size{  # DEPRECATED IN v81
  #deprecate('Please use dbfile_path directly to get full file path');
  my ($self, $window_size, $slice) = @_;

  if($slice){

    if(! (ref($slice) && $slice->isa("Bio::EnsEMBL::Slice"))){
      throw('You must provide a valid Bio::EnsEMBL::Slice');
    }

    $window_size .= '.'.$slice->seq_region_name;
  }

  return $self->dbfile_path.'/result_features.'.$self->name.'.'.$window_size.'.col';
}


1;
