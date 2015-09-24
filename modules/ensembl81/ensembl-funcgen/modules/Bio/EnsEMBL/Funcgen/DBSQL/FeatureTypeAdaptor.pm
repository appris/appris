#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor - Adaptor for fetching and
storing FeatureType objects.

=head1 SYNOPSIS

my $ft_adaptor = $db->get_FeatureTypeAdaptor;

my $feature_type = $ft_adaptor->fetch_by_name("H3K4me3");


=head1 DESCRIPTION

The FeatureTypeAdaptor is a database adaptor for storing and retrieving
Funcgen FeatureType objects.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::FeatureType
Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


#Regulatory evidence information

my %regulatory_evidence_info =
  (
   core =>
   {
    name      => 'Open chromatin & TFBS',
    long_name => 'Open chromatin & Transcription factor binding sites',
    label     => 'DNase1 & TFBS',
    classes   => ['Open Chromatin',
                  'Transcription Factor',
                  'Transcription Factor Complex'],
   },

   non_core =>
   {
    name      => 'Histones & polymerases',
    long_name => 'Histone modifications & RNA polymerases',
    label     => 'Hists & Pols',
    classes   => ['Polymerase',  'Histone'],
   }
  );


#Create class => evidence_type hash
my %regulatory_evidence_classes;

foreach my $evidence_type(keys %regulatory_evidence_info){

  foreach my $class(@{$regulatory_evidence_info{$evidence_type}{classes}}){
    $regulatory_evidence_classes{$class} = $evidence_type;
  }
}


=head2 fetch_by_name

  Arg [1]    : String - name of FeatureType
  Arg [2]    : String (optional) - class of FeatureType
  Arg [3]    : Bio::EnsEMBL::Analysis (optional) - Analysis used to generate FeatureType
  Example    : my $ft = $ft_adaptor->fetch_by_name('H3K4me2');
  Description: Fetches a FeatureType with the given name. Optionally returns an array if
               it exists across classes and called in a list context.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType object (or ARRAY of objects if called 
               in a list context)
  Exceptions : Throws if more than one FeatureType for a given name found and not called 
               in list context.
               Throws if Analysis is defined but not valid.
  Caller     : General
  Status     : At risk

=cut

#Remove support for analysis
#Possibly remove in favour of fetch_all_by_name, due to potential name redundancy between classes.

sub fetch_by_name{
  my ($self, $name, $class, $analysis) = @_;

  throw("Must specify a FeatureType name") if ! defined $name;

  my $constraint = ' name = ? ';
  $constraint   .= ' AND class = ? ' if $class;

  if($analysis){
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis);
    $constraint .= ' AND analysis_id = ? ';
  }

  $self->bind_param_generic_fetch($name,           SQL_VARCHAR);
  $self->bind_param_generic_fetch($class,          SQL_VARCHAR) if defined $class;
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER) if defined $analysis;
  my @fts = @{$self->generic_fetch($constraint)};
  #Can get > 1 if name isredundant between classes

  if( (! wantarray) && (scalar @fts >1) ){
    throw("Found more than one FeatureType:$name\n".
          "Please specify a class and/or analysis argument to disambiguate");
  }

  return wantarray ? @fts : $fts[0];
}


=head2 fetch_all_by_name

  Arg [1]    : String - name of FeatureType
  Arg [2]    : String (optional) - Class of FeatureType
  Arg [3]    : Bio::EnsEMBL::Analysis (optional ) - Analysis used to generate FeatureType
  Example    : my $ft = $ft_adaptor->fetch_by_name('H3K4me2');
  Description: Fetches all FeatureType objects with the given name.
  Returntype : Arrayref Bio::EnsEMBL::Funcgen::FeatureType object 
  Exceptions : Throws if name not defined.
               Throws if Analysis is defined but not valid.
  Caller     : General
  Status     : At risk 

=cut

#Remove support for analysis

sub fetch_all_by_name{
  my ($self, $name, $class, $analysis) = @_;

  throw("Must specify a FeatureType name") if ! defined $name;

  my $constraint = ' name = ? ';
  $constraint   .= ' AND class = ? ' if $class;

  if($analysis){
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis);
    $constraint .= ' AND analysis_id = ? ';
  }

  $self->bind_param_generic_fetch($name,           SQL_VARCHAR);
  $self->bind_param_generic_fetch($class,          SQL_VARCHAR) if defined $class;
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER) if defined $analysis;
  
  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Analysis
  Example    : my @fts = @{$ft_adaptor->fetch_all_by_Analysis($analysis);
  Description: Fetches all FeatureTypes for a given Analysis.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if Analysis not valid
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_Analysis{
  my ($self, $analysis) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis);
  my $constraint = ' analysis_id = ? ';

  #Use bind param method to avoid injection
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER);

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_class

  Arg [1]    : String - class of FeatureType
  Example    : my @fts = @{$ft_adaptor->fetch_all_by_class('Histone')};
  Description: Fetches all FeatureTypes of a given class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if class not defined
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_class{
  my ($self, $class) = @_;

  throw('Must specify a FeatureType class') if(! defined $class);
  my $constraint = ' class = ? ';

  #Use bind param method to avoid injection
  $self->bind_param_generic_fetch($class, SQL_VARCHAR);

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_association

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : my $assoc_ftypes = $ft_adaptor->fetch_all_by_association($ext_feature);
  Description: Fetches all associated FeatureTypes for a given Storable.
               Note: Where appropriate, the main FeatureType for a Storable is
               accessible via the feature_type method.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if arg is not valid or stored
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_association{
  my ($self, $storable) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Storable', $storable);

  $self->_tables([['associated_feature_type', 'aft']]);
  my $table_name = $storable->adaptor->_main_table->[0];

  my $constraint = 'aft.feature_type_id=ft.feature_type_id AND aft.table_name="'.$table_name.
	'" AND aft.table_id='.$storable->dbID;

  my $feature_types =  $self->generic_fetch($constraint);
  $self->reset_true_tables;

  return $feature_types;
}


=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  return (['feature_type', 'ft']);
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  return qw( ft.feature_type_id ft.name        ft.class
             ft.analysis_id     ft.description ft.so_accession
             ft.so_name
           );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Channel objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@result, $ft_id, $name, $class, $anal_id, $desc, $so_acc, $so_name, %analysis_hash);
	my $anal_a = $self->db->get_AnalysisAdaptor;

	$sth->bind_columns(\$ft_id, \$name, \$class, \$anal_id, \$desc, \$so_acc, \$so_name);

	$analysis_hash{0} = undef;

	  while ( $sth->fetch() ) {

	  $anal_id ||= 0;

	  if (! exists $analysis_hash{$anal_id}){
		$analysis_hash{$anal_id} = $anal_a->fetch_by_dbID($anal_id);
	  }

	  my $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
		(
		 -dbID         => $ft_id,
		 -NAME         => $name,
		 -CLASS        => $class,
		 -ANALYSIS     => $analysis_hash{$anal_id},
		 -DESCRIPTION  => $desc,
		 -SO_ACCESSION => $so_acc,
		 -SO_NAME      => $so_name,
		 -ADAPTOR      => $self
		);

	  push @result, $ftype;
	}

	return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::FeatureType objects
  Example    : $chan_a->store($c1, $c2, $c3);
  Description: Stores given Channel objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if FeatureType not valid
               Throws if Analysis defined but not valid
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my ($self, @args) = @_;

  #Prepare once for all ftypes
  my $sth = $self->prepare('INSERT INTO feature_type'.
                           '(name, class, analysis_id, description, so_accession, so_name)'.
                           'VALUES (?, ?, ?, ?, ?, ?)');

  #Process each ftype
  foreach my $ft (@args) {

    #Validate ftype
    if ( ! (ref($ft) && $ft->isa('Bio::EnsEMBL::Funcgen::FeatureType') )) {
      throw('Can only store FeatureType objects, skipping $ft');
    }

    #Validate analysis
    my $anal_id;
    my $analysis = $ft->analysis;

    if($analysis){
      $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis);
      $anal_id = $analysis->dbID;
    }

    #Is already stored?
    if ( $ft->dbID && ( $ft->adaptor == $self) ) {
      warn "Skipping previous stores FeatureType:\t".$ft->name.'('.$ft->dbID.")\n";
    } else {
      #Check for previously stored FeatureType
      my $s_ft = $self->fetch_by_name($ft->name, $ft->class, $ft->analysis);

      if (! $s_ft) {
        $sth->bind_param(1, $ft->name,          SQL_VARCHAR);
        $sth->bind_param(2, $ft->class,         SQL_VARCHAR);
        $sth->bind_param(3, $anal_id,           SQL_INTEGER);
        $sth->bind_param(4, $ft->description,   SQL_VARCHAR);
        $sth->bind_param(5, $ft->so_accession,  SQL_VARCHAR);
        $sth->bind_param(6, $ft->so_name,       SQL_VARCHAR);


        $sth->execute();
        $ft->dbID($self->last_insert_id);
        $ft->adaptor($self);
      }
      else {
        $ft = $s_ft;

        #Check other fields match
        my @failed_methods;

        for my $method (qw(description so_accession so_name)) {
          #Allow nulls/undefs to match empty strings
          my $ft_val   = $ft->$method   || '';
          my $s_ft_val = $s_ft->$method || '';

          if ($ft_val ne $s_ft_val) {
            push @failed_methods, "$method does not match between existing(${s_ft_val}) and new(${ft_val}) FeatureTypes";
          }
        }

        if (@failed_methods) {
          #Could throw, but maybe easier to patch after import?
          warn("Used existing FeatureType with disparities:\n\t".join("\n\t", @failed_methods));
        }
      }
    }
  }

  return \@args;
}

=head2 get_regulatory_evidence_classes

  Args       : String (optional) - Evidence type e.g. 'core' or 'non_core'
  Example    :
  Description:
  Returntype : Arrayref of Strings
  Exceptions : None
  Caller     : web code
  Status     : At risk - remove in favour of get_regulatory_evidence_info

=cut

#Actually returns list if type defined or otherwise array from map

#Remove this?

sub get_regulatory_evidence_classes{
  my ($self, $type) = @_;

  if(defined $type &&
     ! exists $regulatory_evidence_info{$type}){
    throw("The evidence type passed($type) must be one of:\t".
          join(' ',keys(%regulatory_evidence_info)));
  }

  return ($type) ? $regulatory_evidence_info{$type}{classes} :
    (map { @{$_->{classes}} } values %regulatory_evidence_info);
}


=head2 get_regulatory_evidence_label_by_class

  Args       : String (optional) - FeatureType class e.g. Histone
  Example    : my $ftype_label = $ft_adaptor->get_regulatory_label_by_class('Polymerase')};
  Description: Returns the shorts labels used for the grouped tracks i.e.
                   'Hists & Pols'
                   'Dnase1 & TFBS'
  Returntype : String
  Exceptions : None
  Caller     : web code
  Status     : At risk

=cut


#Often quicker to grab regulatory_evidence_info hash first
#and access directly in caller e.g.
#my $ft_class_label =  $ftype_info->
#      {$ft_a->get_regulatory_evidence_type($ft_class)}{label};

sub get_regulatory_evidence_label_by_class{
  my ($self, $fclass) = @_;

  return $self->get_regulatory_evidence_info
    ($self->get_regulatory_evidence_type($fclass))->{label};
}


=head2 get_regulatory_evidence_info

  Args       : String (optional) - Regulatory evidence type i.e. core or non_core
  Example    : my %info = %{$ft_adaptor->get_regulatory_evidence_info('core')};
  Description: Returns all regulatory evidence info keyed on the regulatory
               evidence type. If the type arg is omited, returns entire info hash
               keyed on evidence types
  Returntype : HASHREF
  Exceptions : None
  Caller     : web code
  Status     : At risk

=cut

sub get_regulatory_evidence_info{
  my ($self, $class) = @_;
  my $hash_ref;

  if(defined $class){

    if(! exists $regulatory_evidence_info{$class}){
      warn "FeatureType class $class does not have any regulatory evidence info\n";
    }
    else{
      $hash_ref = $regulatory_evidence_info{$class};
    }
  }
  else{
    $hash_ref = \%regulatory_evidence_info;
  }

  return $hash_ref;
}



=head2 get_regulatory_evidence_type

  Args       : String - FeatureType class
  Example    : my $evidence_Type = $ft_adaptor->get_regulatory_evidence_info($ftype->class);
  Description: Returns the regulatory evidence type i.e, core or non_core
  Returntype : String
  Exceptions : Throws if no class arguments defined
               Warns if the class is not considered as evidence for the regulatory build
  Caller     : General
  Status     : At risk

=cut

sub get_regulatory_evidence_type{
  my ($self, $class) = @_;

  my $evidence_type;

  if(defined $class){

    if (! exists $regulatory_evidence_classes{$class}){
      warn "FeatureType class $class does not have any regulatory evidence info\n";
    }
    else{
      $evidence_type = $regulatory_evidence_classes{$class};
    }
  }
  else{
    throw('You must pass a FeatureType class to get_regulatory_evidence_type');
  }

  return $evidence_type;
}



=head2 fetch_all_by_evidence_type

  Arg [1]    : String - Regulatory build evidence type i.e. core or non_core
  Example    : my @core_fts = @{$ft_adaptor->fetch_all_by_evidence_type('core')};
  Description: Fetches all FeatureTypes which can be considered for the regulatory
               build based on their evidence type. Core FeatureTypes are used to
               construct the core regions of RegulatoryFeatures (e.g. TFs, DNase1
               etc.), and non_core are used to construct the bound regions.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if evidence_type is not defined or valid.
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_evidence_type{
  my ($self, $etype) = @_;

  if(! (defined $etype &&
        exists $regulatory_evidence_info{$etype}) ){
    throw("$etype must be one of the valid evidence types:\t".
          join(' ', keys %regulatory_evidence_info));
  }

  #Don't need to bind_param_generic_fetch here as the constraint is built
  #from valid internal values

  my $constraint = ' class IN ("'.join('", "', @{$regulatory_evidence_info{$etype}{classes}}).'")';
  return $self->generic_fetch($constraint);
}


1;

