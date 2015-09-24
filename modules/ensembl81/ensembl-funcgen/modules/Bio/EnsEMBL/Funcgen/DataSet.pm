#
# Ensembl module for Bio::EnsEMBL::Funcgen::DataSet
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

Bio::EnsEMBL::Funcgen::DataSet - A module to represent DataSet object.


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DataSet;

my $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(
	                                              -DBID            => $dbID,
							 					  -ADAPTOR         => $self,
                                                  -SUPPORTING_SETS => [$rset],
                                                  -FEATURE_SET     => $fset,
                                                  -DISPLAYABLE     => 1,
                                                  -NAME            => 'DATASET1',
                                                  );



=head1 DESCRIPTION

A DataSet object provides access to either or both raw results and AnnotatedFeatures
for a given experiment within a Slice, associated with set wide experimental meta data.
This was aimed primarily at easing access to data via the web API by creating
a wrapper class with convenience methods.  The focus of this class is to contain raw and
associated processed/analysed data to be displayed as a set within the browser i.e. an
experiment may have different cell lines, features or time points, these would require different DataSets.
# However a DataSet may contain mixed data types i.e. promoter & histone???? No give separate sets?
May have duplicates for raw data but only one predicted features track??
The data in this class is kept as lightweight as possible with data being loaded dynamically.


=cut

package Bio::EnsEMBL::Funcgen::DataSet;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );

use base qw(Bio::EnsEMBL::Funcgen::Storable);

#Should not be a Set as is sufficiently different
#_set_Sets_and_types also allows all Sets to be supporting
#but we should not add a DataSet as support

=head2 new


  Args [1]   : Hash of parameters.
               Mandatory:
                -NAME            => 'DATASET1',
                -SUPPORTING_SETS => [$input_set, $result_set],
                -FEATURE_SET     => $product_fset,
               Optional:
                -DISPLAYABLE     => 1,

  Example    : my $dset = Bio::EnsEMBL::Funcgen::DataSet->new(%params)
  Description: Constructor for DataSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : Throws if no -name defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($fset, $ssets, $name) =
    rearrange(['FEATURE_SET', 'SUPPORTING_SETS', 'NAME'], @_);

  #my @caller = caller();
  #if($self->dbID() && $caller[0] ne "Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor"){
  #  throw('You must use the DataSetAdaptor to generate DataSets with dbID i.e. from the DB,'.
  #        ' as this module accomodates updating which may cause incorrect data if the object'.
  #        ' is not generated from the DB');
  #}

  throw('Must defined a DataSet -name') if ! defined $name;
  $self->{name} = $name;
  $self->_set_Sets_and_types($fset, $ssets);
  return $self;
}


=head2 product_FeatureSet

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : $data_set->product_FeatureSet($fset);
  Description: Getter and setter for the main feature_set attribute for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : Throws not a valid FeatureSet or if main feature_set has already been set.
  Caller     : General
  Status     : At Risk - change to get_product_FeatureSet

=cut

sub product_FeatureSet {
  my ($self, $fset) = @_;

  if($fset){

	if (! ($fset && ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"))){
	  throw("Need to pass a valid Bio::EnsEMBL::Funcgen::FeatureSet")
	}

    if(defined $self->{'feature_set'}){
      throw("The main feature_set has already been set for this DataSet, maybe you want add_SupportingSets?");
    }
	else{
	  $self->_validate_and_set_types($fset);
	  $self->{'feature_set'} = $fset;
	}
  }

  return $self->{'feature_set'};
}




=head2 get_supporting_sets_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen:Analysis
  Arg [2]    : (optional) status - e.g 'DISPLAYABLE'
  Example    : my $anal_sets = @{$result_set->get_ResultSets_by_Analysis($analysis)};
  Description: Getter for the SupportingSet objects of a given Analysis.
  Returntype : ARRAYREF
  Exceptions : Throws if arg is not a valid stored Bio::EnsEMBL::Anaylsis
  Caller     : General
  Status     : At Risk

=cut

sub get_supporting_sets_by_Analysis {
  my ($self, $analysis, $status) = @_;


  my @rsets;


  #should we handle displayable here, and propogate to the ResultSet if update_status is set
  #is there scope to write a Funcgen::Storable, which provides convenience methods to StatusAdaptor?
  #would have to make sure Feature object also inherited from Funcgen::Storable aswell as BaseFeature


  if (! ($analysis->isa("Bio::EnsEMBL::Analysis") && $analysis->dbID())){
	  throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::ResultSet");
  }

  #will have to generate new array of object here if we want to filter displayable
  #This may result in returning a ref to the stored ResultSets for no status
  #And a ref to the abstracted/filtered i.e. non-stored ResultSets if we have a status
  #This could cause problems if people want to edit the real ResultSets via the refs
  #If we edit the ResultSets like this, we would still store via their adaptor
  #so would need to refresh DataSet anyway.

  #should ResultSet/Adaptor contain all the fetch_methods, and leave DataSet as a kind of organisational class as a single point of access.
  #DataSetAdaptor to perform the ordering according to feature/celltype
  #This will still not resolve the complex data sets which can be accomodated by the DB.
  #Maybe we can keep the data sets as simple as there are and confer the association by tracking back to the experiment?
  #Would there only ever be one experiment for a complex data_set?


  #Can have more than one experiment for a compound feature set, would we ever want to display raw data?
  #This is actually an easier problem unless we are displaying two feature types(i.e. complex and compound)

  #could we have >1 rset with the same analysis?

  foreach my $anal_rset(@{$self->{'supporting_sets'}->{$analysis->dbID()}}){

	  if(! defined $status){
		  push @rsets, $anal_rset;
	  }
	  elsif($anal_rset->has_status($status)){
		  push @rsets, $anal_rset;
	  }
  }

  return \@rsets;
}


=head2 get_supporting_sets

  Arg [1]    : String - Set type e.g result, input or feature. Optional.
  Arg [2]    : String - Status e.g DISPLAYABLE. Optional
  Example    : my @status_sets = @{$data_set->get_supporting_sets('result', $status)};
  Description: Getter for the ResultSets for this DataSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_supporting_sets{
  my ($self, $set_type, $status)  = @_;
  #swap the args here
  #Add analysis here and make above method wrapper

  my %valid = (
    input  =>  1,
    feature => 1,
    result =>  1,
    );
  
  if(defined $set_type){
    if( !exists $valid{$set_type}){
     throw("You have specified an invalid supporting set type:\t$set_type");
    }
   }

  my @ssets;

  foreach my $anal_id(keys %{$self->{'supporting_sets'}}){
    foreach my $sset(@{$self->{'supporting_sets'}->{$anal_id}}){

	  if(defined $status &&
		 (! $sset->has_status($status))){
		next;
	  }

	  if(defined $set_type && ($sset->set_type ne $set_type)){
		next;
	  }

	  push @ssets, $sset;
    }
  }

#use Data::Dumper; print Dumper(\@ssets); die;
  return \@ssets;
}


=head2 get_displayable_supporting_sets

  Example    : my @displayable_rsets = @{$result_set->get_displayable_supporting_sets()};
  Description: Convenience method for web display
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_displayable_supporting_sets{
  my ($self, $set_type) = @_;

  return $self->get_supporting_sets( $set_type, 'DISPLAYABLE' );
}


=head2 get_displayable_product_FeatureSet

  Example    : my $fset = $data_set->get_displayable_product_FeatureSet();
  Description: Convenience method for web display
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_displayable_product_FeatureSet{
  my $self = shift;

  return  $self->product_FeatureSet->has_status('DISPLAYABLE') ?  $self->product_FeatureSet() : undef;
}


=head2 name

  Example    : print "This is my data set:\t".$dset->name."\n";
  Description: Getter for the name of this DataSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return $_[0]->{name}; }


#The following attributes are generated dynamically from the
#consituent Result/FeatureSets

=head2 cell_type

  Example    : my $dset_ctype_name = $dset->cell_type->name();
  Description: Getter for the cell_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub cell_type { return $_[0]->{cell_type}; }


=head2 feature_type

  Example    : my $dset_ftype_name = $dset->feature_type->name();
  Description: Getter for the feature_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_type { return $_[0]->{feature_type}; }


=head2 display_label

  Example    : print $dset->display_label();
  Description: Getter for the display_label attribute for this DataSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_label {
  my $self = $_[0];

  if(! defined $self->{display_label}){

	if($self->product_FeatureSet->feature_class ne 'annotated'){
	  $self->{'display_label'} = $self->name;
	}
	else{
	  $self->{'display_label'}  = $self->feature_type->name()." - ";
	  $self->{'display_label'} .= ($self->cell_type->display_label ||
								   $self->cell_type->name)." Enriched Sites";
	}
  }

  return $self->{display_label};
}


=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of DataSet method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: name table_name feature_class get_all_states
  Args[4]    : Arrayref - Optional list of DataSet method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: feature_type cell_type analysis get_support
  Example    : my %shallow_diffs = %{$dset->compare_to($other_rset, 1)};
  Description: Compare this DataSet to another based on the defined scalar
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

#shallow should really omit states here
#this woudl imply states should be objects which is overkill
#added skip_states flag for now

sub compare_to {
  my ($self, $obj, $shallow, $scl_methods, $obj_methods, $skip_states) = @_;

  $obj_methods ||= [qw(product_FeatureSet get_supporting_sets)];
  $scl_methods ||= [qw(name get_all_states)];

  if($skip_states){
    pop @$scl_methods;
  }

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods,
                                  $obj_methods);
}

=head2 reset_relational_attributes

  Arg[1] : Hashref containing the following mandatory parameters:
            -analysis     => Bio::EnsEMBL::Analysis,
            -feature_type => Bio::EnsEMBL::Funcgen::FeatureType,
            -cell_type    => Bio::EnsEMBL::Funcgen::CellType,
            -support      => Arrayref of valid support objectse.g. InputSet

  Description: Resets all the relational attributes of a given DataSet.
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;
  my ($ssets, $feature_set) = rearrange(['SUPPORTING_SETS', 'FEATURE_SET'],
                                %$params_hash);

  #This also sets cell_type and feature_type
  $self->_set_Sets_and_types($feature_set, $ssets);

  #also flush dynamically set display label
  #just in case there is a difference (there shouldn't be!)
  undef $self->{display_label};

  #Finally undef the dbID and adaptor by default
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  return;
}


#Currently does not support DataSets with only mixed type support
#i.e. not FeatureSet and CellTypes or FeatureTypes difference between
#supporting sets


sub _set_Sets_and_types{
  my ($self, $fset, $ssets) = @_;

  assert_ref($ssets, 'ARRAY', 'Supporting sets');
  
  if(scalar(@$ssets) < 1){
    throw('Must pass an arrayref of supporting Sets e.g. InputSubsets, ResultSets or FeatureSets');  
  }
  
  
  my ($ftype, $ctype, $fclass, $ftype_name, $ctype_name);

  if(defined $fset){
    assert_ref($fset, 'Bio::EnsEMBL::Funcgen::FeatureSet');
    $ftype      = $fset->feature_type;
    $ctype      = $fset->cell_type;
    $fclass     = $fset->feature_class;
    $ftype_name = $ftype->name;
    $ctype_name = $ctype->name;
  }


  ### Validate supporting sets
  #Need to validate cell_type if the fset feature_class is not regulatory
  #Need to validate feature_type if the fset feature class is not regulatory or segmentation
  $self->{'supporting_sets'} = {};

  foreach my $set(@$ssets){

    if(! (defined $set &&
          ref($set) &&
          $set->isa('Bio::EnsEMBL::Funcgen::Set') )){
      throw('All -supporting_sets for a DataSet must be a '.
            "Bio::EnsEMBL::Funcgen::Set\n\te.g.InputSet, ResultSet or FeatureSet '".ref($set)."'");
    }
    
    if($fclass){
      
      if($fclass ne 'regulatory'){ #Validate cell type
        
        if($set->cell_type->name ne $ctype_name){
          throw('Cannot add '.$set->cell_type->name.
                " support to a $ctype_name"." $fclass DataSet");
        }
        
        if($fclass ne 'segmentation'){# and ne regulatory validate ftype
          #We don't have segmentation data sets just yet
          
          if($set->feature_type->name ne $ftype_name){
            throw('Cannot add '.$set->feature_type->name.
                  " support to a $ftype_name"." $fclass DataSet");
          }
        }
      }
    }
    else{
      $ftype_name ||= $set->feature_type->name;
      $ctype_name ||= $set->cell_type->name;
  
      if($ftype_name ne $set->feature_type->name){
        throw('Unable to set distinct FeatureType for a mixed support '.
              "DataSet without a product FeatureSet:\t".$self->name);
      }
      
      if($ctype_name ne $set->cell_type->name){
        throw('Unable to set distinct CellType for a mixed support '.
              "DataSet without a product FeatureSet:\t".$self->name);
      }
    }
    
    $self->{'supporting_sets'}->{$set->analysis->dbID()} ||= [];
    push @{$self->{'supporting_sets'}->{$set->analysis->dbID()}}, $set;
  }

  if(! defined $ftype){#and $ctype by proxy
    $ftype = $ssets->[0]->feature_type;
    $ctype = $ssets->[0]->cell_type;
  }

  #Reset dynamically defined attrs
  $self->{feature_type}    = $ftype;
  $self->{cell_type}       = $ctype;
  $self->{feature_set}     = $fset;
  return;
}


### DEPRECATED



sub add_supporting_sets { #deprecated in v72
  my ($self, $sets) = @_;

	throw('add_supporting_set is deprecated, please use -supporting_sets in DataSet::new');

  #should we handle displayable here, and propogate to the ResultSet if update_status is set
  #is there scope to write a Funcgen::Storable, which provides convenience methods to StatusAdaptor?
  #would have to make sure Feature object also inherited from Funcgen::Storable aswell as BaseFeature

  throw("Supporting sets need to be a reference to an ARRAY:\t".$sets) if ref($sets) ne 'ARRAY';

  foreach my $set(@$sets){

	if(!(ref($set) &&  $set->isa('Bio::EnsEMBL::Funcgen::Set') && $set->set_type ne 'data' && $set->dbID)){
	  throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::Set which is not a DataSet:\t$set");
	}
	#set type cannot be data at present as it does not inherit from Set.pm



	#Only validate if we are dealing with result type data
	#As we can have various cell/feature_types for compound analyses e.g. RegulatoryFeatures

	$self->_validate_and_set_types($set) if $set->set_type() ne 'feature';

	#should ResultSet/Adaptor contain all the fetch_methods, and leave DataSet as a kind of organisational class as a single point of access.
	#DataSetAdaptor to perform the ordering according to feature/celltype
	#This will still not resolve the complex data sets which can be accomodated by the DB.
	#Maybe we can keep the data sets as simple as there are and confer the association by tracking back to the experiment?
	#Would there only ever be one experiment for a complex data_set?


	#Can have more than one experiment for a compound feature set, would we ever want to display raw data?
	#This is actually an easier problem unless we are displaying two feature types(i.e. complex and compound)


	$self->{'supporting_sets'}->{$set->analysis->dbID()} ||= ();
	push @{$self->{'supporting_sets'}->{$set->analysis->dbID()}}, $set;
  }

  return;
}


sub _validate_and_set_types{ #in v72, remove immediately as is *private*?
  my ($self, $set) = @_;

  #slightly dodgy bypassing methods, but extendable

  #This currently restricts all set types to one cell and feature type
  #this is incorrect for feature_set types as we want to munge several feature and possibly cell types
  #into one combined data set.
  #this should set it to the FeatureSet type if is feature_set data_set
  #this only works as we only validate supporting_sets if type is not feature

  for my $type('feature_type', 'cell_type'){

    if(defined $self->{$type}){

#Need to test isa here?  Why is this passing the defined test if not set?
      if($set->{$type}->name() ne $self->{$type}->name()){

        throw(ref($set)." $type(".$set->{$type}->name().
            ") does not match DataSet $type(".$self->{$type}->name().")");
      }
    }
    else{
      $self->{$type} = $set->{$type};
    }
  }

  return;
}


1;

