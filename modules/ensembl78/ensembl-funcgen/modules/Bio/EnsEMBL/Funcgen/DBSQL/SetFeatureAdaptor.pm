#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::SetFeatureAdaptor
#

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::DBSQL::Funcgen::SetFeatureAdaptor - A base database adaptor for SetFeature adaptors.
storing SetFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_SetFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);


=head1 DESCRIPTION

The SetFeatureAdaptor is a base adaptor for all SetFeature adaptors.
e.g. AnnotatedFeatureAdaptor, RegulatoryFeatureAdaptor etc.
It provides common methods across all feature types.



=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::SetFeature

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Bio::EnsEMBL::Funcgen::SetFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;#DBI sql_types import

use base qw( Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor );

use vars qw( @EXPORT ); #require Exporter is done in BaseAdaptor
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

use base qw( Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor );



#required for subclass adaptor store and _obj_from_sth methods

=head2 fetch_all_by_FeatureType_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [3]    : (optional) hashref - params hash
                                      associated => 1, #Also return feature which have the associated FeatureType
                                      logic_name => 'analysis.logic_name'
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice with main or associated FeatureType. 
               This is mainly used by external FeatureSets which can sometimes have more 
               than one associated FeatureType. NOTE: This is intended to work for FeatureTypes at the 
               feature level, not the more generic FeatureSet level FeatureTypes.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::SetFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

#Needs reviewing
#We are not restricting to feature_set in the associated ftypes call
#Can't this be done with one compose_query_constraint call?

sub fetch_all_by_FeatureType_FeatureSets {
  my ($self, $ftype, $fsets, $params) = @_;

  if ($self->_feature_class eq 'annotated'){
	throw("This method is not appropriate for FeatureSets with feature_class 'annotated', please use standard fetch_all_by_FeatureSets or get_all_Features on an individual FeatureSet");
	#There is feature_type_id in annotated_feature
	#Could warn and redirect
  }

  $params->{'logic_name'} ||= undef;
  $params->{'associated'} ||= undef;
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);
  my ($table_name, $table_syn) = @{$self->_main_table};

  my $constraint = $table_syn.'.feature_set_id = fs.feature_set_id AND '.
	$table_syn.'.feature_type_id='.$ftype->dbID.' AND '.$table_syn.".feature_set_id ".$self->_generate_feature_set_id_clause($fsets);

  #Should really pass the params hash on here
  $constraint = $self->_logic_name_to_constraint($constraint, $params->{logic_name});


  my @features = @{$self->generic_fetch($constraint)};


  #This is an interim solution and really needs changing!
  #Can we genericise this as a lazy loader function?
  #Isn't there already something like this in core?
  if ($params->{'associated'}){

	for my $fset(@{$fsets}){
	  #We want to bring back features from the same set
	  #We are not bringing back different class of feature i.e. AnnotatedFeatures and ExternalFeatures
	  #Let's not catch this, but let it silently fail so peaople.

	  next if $table_name ne lc($fset->feature_class).'_feature';

	  my $sql = 'SELECT table_id from associated_feature_type where table_name="'.$table_name.'" and feature_type_id='.$ftype->dbID;

  

	  my @dbids = map $_ = "@$_", @{$self->dbc->db_handle->selectall_arrayref($sql)};

	  if(@dbids){
		$constraint = " $table_syn.${table_name}_id in (".join(',',@dbids).') ' if @dbids;
		push @features, @{$self->generic_fetch($constraint, $params->{logic_name})};
	  }
	}
  }

  return \@features;
}


=head2 fetch_all_by_Feature_associated_feature_types

  Arg [1]    : Bio::EnsEMBL::Funcgen::SetFeature
  Arg [2]    : (optional) hashref - params hash, all entries optional
                                      -logic_name => 'analysis.logic_name'
                                      -include_direct_links => 1, #Also return feature which are linked by Feature->feature_type
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $sfa->fetch_all_by_Feature_associated_feature_types($feature);
  Description: Retrieves a list of all features linked via the associated FeatureTypes of 
               the given Feature in the same FeatureSet.  This is mainly used by external
               FeatureSets which can sometimes have more than one associated FeatureType.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::SetFeature objects
  Exceptions : Throws if SetFeature not stored and valid
  Caller     : General
  Status     : At Risk

=cut

#CR do we need to fetch based on main ftype vs assoc ftpes and vice versa

sub fetch_all_by_Feature_associated_feature_types {
  my ($self, $feature, $params) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::SetFeature', $feature);

  $params->{'include_direct_links'} ||= undef;
  $params->{'logic_name'}           ||= undef;

  my %dbIDs;#Use a hash to reduce dbID redundancy
  my ($table_name, $table_syn) = @{$self->_main_table};
  my $fset = $feature->feature_set;
  my ($sql, $constraint, @features);

  #Direct FeatureType
  if($params->{'include_direct_links'}){
	#Just grab dbIDs here rather than use generic fetch
	$sql = "SELECT ${table_name}_id from $table_name where feature_type_id=".$feature->feature_type->dbID.' and feature_set_id='.$fset->dbID;

	#This just sets each value to a key with an undef value
	map {$dbIDs{"@$_"} = undef} @{$self->dbc->db_handle->selectall_arrayref($sql)};
  }

 
    
  my @assoc_ftypes = @{$feature->associated_feature_types};

  if(@assoc_ftypes){

	my $ftype_ids = join(', ', (map $_->dbID, @assoc_ftypes));
  
	#Now we need to restrict the Features based on the FeatureSet of the original Feature, we could make this optional.
	$sql = "SELECT table_id from associated_feature_type aft, $table_name $table_syn where aft.table_name='".$fset->feature_class."_feature' and aft.feature_type_id in ($ftype_ids) and aft.table_id=${table_syn}.${table_name}_id and ${table_syn}.feature_set_id=".$fset->dbID;
	#This just sets each value to a key with an undef value
	
	map {$dbIDs{"@$_"} = undef} @{$self->dbc->db_handle->selectall_arrayref($sql)};
  }



  if(keys %dbIDs){
	$constraint = " $table_syn.${table_name}_id in (".join(',', keys %dbIDs).')';
	push @features, @{$self->generic_fetch($constraint, $params->{logic_name})};
  }

  return \@features;
}




=head2 fetch_all_by_Slice_FeatureType

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) string - analysis logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice, specific for a given FeatureType.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::SetFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureType {
  my ($self, $slice, $type, $logic_name) = @_;

  my $params                      = {constraints  => {feature_types => [$type]}};
  $params->{optional_constraints} = {logic_names  => [$logic_name]} if defined $logic_name;  
  my $constraint                  = $self->compose_constraint_query($params);
  my $feats                       = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
  $self->reset_true_tables; #logic_name adds analysis
  return $feats;
}

#todo fetch_all_by_Slice_CellType


=head2 fetch_all_by_FeatureSets

  Arg [1]    : Arrayref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Arg [2]    : optional - analysis.logic_name
  Example    : my $features = $set_feature_adaptor->fetch_all_by_FeatureSets(@fsets);
  Description: Retrieves a list of features specific for a given list of FeatureSets.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk - Maybe not sensible to fetch all data in one call.

=cut

sub fetch_all_by_FeatureSets {
  my ($self, $fsets, $logic_name) = @_;	
  
  my $params                      = {constraints => {feature_sets => $fsets}};
  $params->{optional_constraints} = {logic_names => [$logic_name]} if defined $logic_name;  
  my $constraint                  = $self->compose_constraint_query($params);
  my $feats                       = $self->generic_fetch($constraint);
  $self->reset_true_tables; #logic_name adds analysis
  return $feats;
}


=head2 _generate_feature_set_id_clause

  Arg [1]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Example    : my $fset_d_clause = $self->_generate_feature_set_id_clause($fsets);
  Description: Build feature_set id clause for FeatureSet methods
  Returntype : string
  Exceptions : Throws if FeatureSets are passed
               Throws if FeatureSet feature_class does not match adaptor feature_class
               Throws if FeatureSet is not valid 
  Caller     : Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor
  Status     : At Risk

=cut

#replace with _constrain_feature_sets

sub _generate_feature_set_id_clause{
  my ($self, $fsets) = @_;

  my @fs_ids;

  if(! ( (ref($fsets) eq 'ARRAY') &&
         scalar(@{$fsets}) > 0) ){
    throw('Must provide a list of Bio::EnsEMBL::Funcgen::FeatureSet objects');
  }

  my $fclass = $self->_feature_class;
  
  foreach my $fset (@{$fsets}) {
	
	if (! (ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"))){
	  throw('Not a FeatureSet object');
	}

	if($fset->feature_class ne $fclass){
	  throw('FeatureSet feature_class \''.$fset->feature_class.'\' does not match adaptor feature_class \''.$fclass.'\'');
	}

	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);

	push (@fs_ids, $fset->dbID());
  }

  return scalar(@fs_ids >1) ? 'IN ('.join(',', @fs_ids).')' : '= '.$fs_ids[0];
}



=head2 fetch_all_by_Slice_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Arg [3]    : optional - analysis.logic_name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureSets($slice, \@fsets);
  Description: Retrieves a list of features on a given slice, specific for a given list of FeatureSets.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureSets {
  my ($self, $slice, $fsets, $logic_name) = @_;

  my $params                      = {constraints  => {feature_sets => $fsets}};
  $params->{optional_constraints} = {logic_names  => [$logic_name]} if defined $logic_name;  
  my $constraint                  = $self->compose_constraint_query($params); 
  my $feats                       = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
  $self->reset_true_tables; #logic_name adds analysis
  return $feats;
}




=head2 fetch_Iterator_by_Slice_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Arg [3]    : optional - analysis.logic_name
  Arg [4]    : optional - iterator chunk length. Default is 1MB
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $iter = $ofa->fetch_Iterator_by_Slice_FeatureSets($slice, \@fsets);

           
  Description: Simple Iterator wrapper method for fetch_all_by_Slice_FeatureSets
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk

=cut


sub fetch_Iterator_by_Slice_FeatureSets{
  my ($self, $slice, $fsets, $logic_name, $chunk_length) = @_;

  return $self->fetch_Iterator_by_Slice_method
	($self->can('fetch_all_by_Slice_FeatureSets'),
	 [$slice, $fsets, $logic_name],
	 0,#Slice idx
	 $chunk_length #default is 1000000
	);
}




#This is all done at the level of the feature_set

#at risk remove this in favour of constrain method

sub _logic_name_to_constraint {
  my $self = shift;
  my $constraint = shift;
  my $logic_name = shift;

  return $constraint if (!$logic_name);

  my $aa = $self->db->get_AnalysisAdaptor();
  my $an = $aa->fetch_by_logic_name($logic_name);

  if(! $an) {
	  #warn or throw?
	  warn("No analysis associated with logic_name $logic_name");
	  return undef;
  }

  my $an_id = $an->dbID();

  $constraint .= ' AND' if($constraint);
  $constraint .= " fs.analysis_id = $an_id";
  return $constraint;
}



=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

#This is to facilitate fetches based on feature_set fields
#but we could remove this and add this via _tables and the constraint methods
#as it's not strictly a part of the data model for set features

sub _default_where_clause {
  return $_[0]->_main_table->[1].'.feature_set_id = fs.feature_set_id';
}


=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an GROUP/ORDER BY clause. 
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : At Risk

=cut

#This is required as tables constain mutliple data sets loaded at different times
#hence we may not get the expected or of features returned

sub _final_clause {
  my $self = shift;
  return ' ORDER BY '.$self->_main_table->[1].'.seq_region_id, '.$self->_main_table->[1].'.seq_region_start';
}


=head2 fetch_all_by_logic_name

  Arg [1]     : string $logic_name the logic name of the analysis of features to obtain 
  Example     : $fs = $a->fetch_all_by_logic_name('foobar'); 
  Description : Returns an arrayref of features created from the database. only 
                features with an analysis of type $logic_name will be returned. 
  Returntype  : arrayref of Bio::EnsEMBL::SetFeatures
  Exceptions  : thrown if $logic_name not defined
  Caller      : General 
  Status      : At risk

=cut 

#re-implemented for non-standard feature table i.e. points at feature_set

sub fetch_all_by_logic_name { 
  my ($self, $lname) = @_; 
  #my $logic_name = shift || throw( "Need a logic_name" ); 
  #$constraint = ' '.$self->_main_table->[1].'.feature_set_id=fs.feature_set_id and fs.analysis_id = '.$an->dbID() if($an);
  #my $an = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name); 
  #return (defined $constraint) ? $self->generic_fetch($constraint) : undef;
  my $constraint = $self->compose_constraint_query({constraints => {logic_names => [$lname]}});
  
  return $self->generic_fetch($constraint) || undef;
} 


=head2 _feature_class

  Example    : if($self->feature_class ne $fset->feature_class){ throw('some error'); }
  Description: Convenience method to retrieve the feature class for this adaptor
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub _feature_class{
  (my $fclass = $_[0]->_main_table->[0]) =~ s/_feature$//; #use the first word of the table name as the class
  return $fclass;
}


### GENERIC CONSTRAIN METHODS ###

#Need to get the relevant SetAdaptor here
#This requires a mapping between the feature type and the set type
#without access to a feature or set.
#Do we have this anywhere yet?
#This will always be feature set or result set
#and resultset adaptor won't use these constraints as it is a file adaptor
# and need separating from the BaseAdaptor

sub _constrain_cell_types {
  my $self = shift;
  my $cts  = shift;

  #Don't need to bind param this as we validate
  my $constraint = " fs.cell_type_id IN (".
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $cts, 'dbID')}).')';

  return ($constraint, {});  #{} = no further config
}



#This may be redefined in subclass adaptors if the feature table 
#implementation has a more specific feature type thatn the feature_set

#could we do this here conditional on the feature type?
#or a has_feature_type flag, defined in new?
#or grep columns for table_syn.feature_type_id here?

sub _constrain_feature_types {
  my $self = shift;
  my $fts  = shift;
 
  #Don't need to bind param this as we validate
  my $constraint = " fs.feature_type_id IN (".
        join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fts, 'dbID')}).')';  
  
  return ($constraint, {});   #{} = not futher constraint conf
}


#This currently doesn't catch array with undef value(s)
#Passing an undef $logic_name to the method would not normally result in a throw
#unless it is a mandatory method
#so we need to catch this before here really
#and make this validate as though it is mandatory

#This basically catches unsafe specification of arrays with undef values
#as these will cause errors or no results


#This should use analysis_adpator cache and return undef?
#Genericise this in BaseFeatureAdaptor by adding a call
#to _meta_info_table which would return relevant set table for SetFeatures
#else just the feature table?
#Isn't this already wrong for ResultFeature? Or is that not a SetFeature?


sub _constrain_logic_names {
  my $self        = shift;
  my $logic_names = shift;
  assert_ref($logic_names, 'ARRAY');
  
  if(! @$logic_names){
    throw('Must pass an Arraref of logic_name strings to contrain by');  
  }

  for my $lname(@$logic_names){
    if(! defined $lname){
      throw('Found undefined logic_name value');  
    }  
  }


  $self->_tables([['analysis', 'a']]);

  my $constraint = ' fs.analysis_id = a.analysis_id and a.logic_name in ("'.
    join('", "', @$logic_names).'")';
    
  return ($constraint, {});   #{} = not futher constraint conf
}

sub _constrain_feature_sets {
  my $self  = shift;
  my $fsets = shift;
 
  #Don't need to bind param this as we validate
  #match on fs.feature_set_id rather than feature_table feature_set_id
  #to avoid having to get table syn here.
  my $constraint = " fs.feature_set_id IN (".
        join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fsets, 'dbID')}).')';  
  
  #This currently does not throw if FeatureSet feature_class does not match adaptor
  
  return ($constraint, {});   #{} = not futher constraint conf
}

#See fetch_all_by_Feature_associated_feature_types above
#Not quite as straight forward
#sub _constrain_associated_feature_types {
#  
#}

1;
