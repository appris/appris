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

=cut

package Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Analysis;

use base qw( Bio::EnsEMBL::Funcgen::Parsers::BaseImporter );



# Base functionality for external_feature parsers

#Make this inherit from Helper?
#Then change all the prints to logs

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  #validate and set type, analysis and feature_set here
  my ($type, $db, $archive, $import_fsets) = rearrange(['TYPE', 'DB', 'ARCHIVE', 'IMPORT_SETS'], @_);

  #What is ExternalParser specific here?
  #archive?
  #type? is this even used?


  #throw('You must define a type of external_feature to import') if(! defined $type);

  if (! ($db && ref($db) &&
		 $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	throw('You must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }

  throw('You can only specify either -clobber|rollback or -archive, but not both') if($self->rollback && $archive);

  $self->{display_name_cache} = {};
  $self->{'db'} = $db;
  #$self->{type} = $type;
  $self->{archive} = $archive if defined $archive;


  #This is not fully implemented yet and need to be validated against the config feature_set
  #pass something like set1,set2 and split and validate each.
  #Or do this in the calling script?

  throw('-import_sets not fully implemented yet') if defined $import_fsets;
  $self->{'import_sets'} = (defined $import_fsets) ? @{$import_fsets} : undef;

  $self->log("Parsing and loading $type ExternalFeatures");

  return $self;

}


=head2 import_sets

  Args       : None
  Example    : foreach my $import_set_name(@{$self->import_sets}){ ... do the import ... }
  Description: Getter for the list of import feature set names, defaults to all in parser config.
  Returntype : Arrayref of import feature_set names
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub import_sets{
  my $self = shift;

  return $self->{'import_sets'} || [keys %{$self->{static_config}{feature_sets}}];
}


=head2 set_feature_sets

  Args       : None
  Example    : $self->set_feature_sets;
  Description: Imports feature sets defined by import_sets.
  Returntype : None
  Exceptions : Throws if feature set already present and rollback or archive not set
  Caller     : General
  Status     : Medium Risk

=cut

#This is done after validate and store feature_types
#Updating this will require making all external parsers use 'static_config'

sub set_feature_sets{
  my $self = shift;

  throw('Must provide a set feature_set config hash') if ! defined $self->{static_config}{feature_sets};


  my $fset_adaptor = $self->db->get_FeatureSetAdaptor;
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor;

  foreach my $fset_name(@{$self->import_sets}){

	$self->log("Defining FeatureSet:\t$fset_name");
	my $fset = $fset_adaptor->fetch_by_name($fset_name);

	#we don't need data sets for external_feature sets!
	#Compare against config after we have merged with defined_anld_validate etc

	if(defined $fset){
	  $self->log("Found previous FeatureSet $fset_name");

	  if($self->rollback){

		$self->rollback_FeatureSet($fset);#Need to pass \@slices here?
	  }
	  elsif($self->{'archive'}){
		my $archive_fset =  $fset_adaptor->fetch_by_name($fset_name."_v".$self->{'archive'});

		if(defined $archive_fset){
		  throw("You are trying to create an archive external feature_set which already exists:\t${fset_name}_v".$self->{archive});
		}

		my $sql = "UPDATE feature_set set name='$fset_name}_v".$self->{archive}."' where name='$fset_name'";
		$self->db->dbc->do($sql);
		undef $fset;
	  }else{
		throw("You are trying to create an external feature_set which already exists:\t$fset_name\nMaybe to want to rollback or archive?");
	  }
	}

	#Assume using static config for now
	#Will need to resolve this when it become generic
	#Maybe we set outside of config!
	#simply as analyses, feature_sets and feature_types?
	my $fset_config = 	$self->{static_config}{feature_sets}{$fset_name}{feature_set};

	if(! defined $fset){
	  my ($name, $analysis, $ftype, $display_label, $desc);

	  my $fset_analysis_key = (exists ${$fset_config}{-analysis})      ? '-analysis'      : '-ANALYSIS';
	  my $fset_name_key     = (exists ${$fset_config}{-name})          ? '-name'          : '-NAME';
	  my $fset_ftype_key    = (exists ${$fset_config}{-feature_type})  ? '-feature_type'  : '-FEATURE_TYPE';
	  my $fset_dlabel_key   = (exists ${$fset_config}{-display_label}) ? '-display_label' : '-DISPLAY_LABEL';
	  my $fset_desc_key     = (exists ${$fset_config}{-description})   ? '-description'   : '-DESCRIPTION';

	  my $fset_fclass_key   = (exists ${$fset_config}{-feature_class}) ? $fset_config->{-feature_class} : 'external';
	  my $display_name      = (exists ${$fset_config}{$fset_dlabel_key}) ? $fset_config->{$fset_dlabel_key} : $fset_name;
	  #fset config name be different from key name
	  my $fs_name           = (exists ${$fset_config}{$fset_name_key}) ? $fset_config->{$fset_name_key} : $fset_name;
	  #warn if they are different?


	  #Can't just deref config hash here as we need to deref the nested feature_type and analysis attrs

	  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(-name          => $fs_name,
													 -feature_class => $fset_fclass_key,
													 -analysis      => ${$fset_config->{$fset_analysis_key}},
													 -feature_type  => ${$fset_config->{$fset_ftype_key}},
													 -display_label => $display_name,
													 -description   => $fset_config->{$fset_desc_key}		 );
	  ($fset) = @{$self->db->get_FeatureSetAdaptor->store($fset)}; 
	}

	#Now replace config hash with object
	#Will this reset in hash or just locally?
	#$fset_config = $fset;
	$self->{static_config}{feature_sets}{$fset_name}{feature_set} = $fset;
  }

  return;
}

#Can't use this anymore as we have to use static_config for all parsers which use set_feature_sets

#
#=head2 validate_and_store_feature_types
#
#  Args       : None
#  Example    : $self->validate_and_store_feature_types;
#  Description: Imports feature types defined by import_sets.
#  Returntype : None
#  Exceptions : None
#  Caller     : General
#  Status     : High Risk - Now using BaseImporter::validate_and_store_config for vista
#
#=cut
#
##Change all external parsers to use BaseImporter::validate_and_store_config
#
#sub validate_and_store_feature_types{
#  my $self = shift;
#
#  #This currently only stores ftype associated with the feature_sets
#  #Havent't we done this already in the InputSet parser
#  #Need to write BaseImporter and inherit from there.
#
#  #InputSet does all loading, but depends on 'user_config'
#  #Where as we are using hardcoded config here
#  #Which are import_sets currently defaults to feature_sets keys
#
#  #we could simply call this static_config and let user_config over-write static config with warnings?
#  #on an key by key basis? (top level only?)
#
#
#  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
#
#  foreach my $import_set(@{$self->import_sets}){
#
#	my $ftype_config = ${$self->{static_config}{feature_sets}{$import_set}{feature_type}};
#	my $ftype = $ftype_adaptor->fetch_by_name($ftype_config->{'name'});
#
#	$self->log("Validating $import_set FeatureType:\t".$ftype_config->{'name'});
#
#	if(! defined $ftype){
#	  $self->log("FeatureType '".$ftype_config->{'name'}."' for external feature_set ".$self->{'type'}." not present");
#	  $self->log("Storing using type hash definitions");
#
#	  $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(
#													   -name => $ftype_config->{'name'},
#													   -class => $ftype_config->{'class'},
#													   -description => $ftype_config->{'description'},
#													  );
#	  ($ftype) = @{$ftype_adaptor->store($ftype)};
#	}
#
#	#Replace hash config with object
#	$self->{static_config}{feature_types}{$ftype_config->{'name'}} = $ftype;
#  }
#
#  return;
#}
#







1;
