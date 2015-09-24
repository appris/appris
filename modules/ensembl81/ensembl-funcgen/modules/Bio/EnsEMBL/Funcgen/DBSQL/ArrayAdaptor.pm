#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor - A database adaptor for fetching and
storing Funcgen Array objects.

=head1 SYNOPSIS

my $oaa = $db->get_ArrayAdaptor();

my $array = $oaa->fetch_by_name_vendor('HG-U133A', 'AFFY');
my @arrays = @{$oaa->fetch_all()};

=head1 DESCRIPTION

The ArrayAdaptor is a database adaptor for storing and retrieving
Funcgen Array objects.


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::Array
Bio::EnsEMBL::Funcgen::ArrayChip

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#sql_types bareword import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_by_array_chip_dbID

  Arg [1]    : Int - dbID of an ArrayChip
  Example    : my $array = $array_adaptor->fetch_by_array_chip_dbID($ac_dbid);
  Description: Retrieves Array object based on one of it's constituent ArrayChip dbIDs
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

#Changed to use simple query extension
#Removed 1 query
#3.7 % or 1.04 times faster

sub fetch_by_array_chip_dbID {
  my ($self, $ac_dbid) = @_;

  throw('Must provide an ArrayChip dbID') if ! $ac_dbid;

  #Extend query tables
  $self->_tables([['array_chip', 'ac']]);

  #Extend query and group
  my $array = $self->generic_fetch('ac.array_chip_id='.$ac_dbid.' and ac.array_id=a.array_id GROUP by a.array_id')->[0];
  $self->reset_true_tables;

  return $array;
}


=head2 fetch_by_name_vendor

  Arg [1]    : String - Name of an array e.g. MoGene-1_0-st-v1
  Arg [2]    : String - (optional) Name of vendor e.g. AFFY
  Example    : #You can list the possible array name and vendor values as follows
               map print $_->vendor.' '.$_->name, @{$array_adaptor->fetch_all};

	           #Use some of the above values with this method
               my $array = $array_adaptor->fetch_by_name_vendor('MoGene-1_0-st-v1', 'AFFY');
  Description: Retrieves a named Array object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : Throws is name argument not defined
               Throws if vendor argument not defined and more than one arrays is found
  Caller     : General
  Status     : Stable

=cut

sub fetch_by_name_vendor {
    my ($self, $name, $vendor) = @_;

    throw("Must provide and name") if (! $name);
	my @arrays;

	if(! $vendor){
	  @arrays = @{$self->generic_fetch("a.name = '$name'")};

	  if(scalar(@arrays) > 1){
		throw("There is more than one array with this name please specify a vendor argument as one of:\t".join(' ', (map $_->vendor, @arrays)));
	  }
	}
	else{
	  #name vendor is unique key so will only ever return 1
	  @arrays = @{$self->generic_fetch("a.name = '$name' and a.vendor='".uc($vendor)."'")};
	}

    return $arrays[0];
}


=head2 fetch_by_name_class

  Arg [1]    : String - Name of an array
  Arg [2]    : String - Class of array e.g. AFFY_UTR
  Example    : #You can list the possible array name and class values as follows
               map print $_->class.' '.$_->name, @{$array_adaptor->fetch_all};

               #Use some of the above values with this method
               my $array = $array_adaptor->fetch_by_name_class('HuGene_1_0_st_v1', 'AFFY_ST');
  Description: Retrieves Array object from the database based on name and class.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : Throws is name and class not passed
  Caller     : General
  Status     : Stable

=cut

sub fetch_by_name_class {
  my ($self, $name, $class) = @_;
  throw("Must provide and array and class e.g.'HuGene_1_0_st_v1', 'AFFY_ST'") if (! ($name && $class));

  #name class is unique key so will only ever return 1
  return $self->generic_fetch("a.name = '$name' and a.class='".uc($class)."'")->[0];
}


=head2 fetch_all_by_class

  Arg [1]    : String - Class e.g. ILLUMINA_WG
  Example    : #You can list the possible array class values as follows
               map print $_->class, @{$array_adaptor->fetch_all};

               #Use some one the above values with this method
               my $array = $array_adaptor->fetch_all_by_class('AFFY_ST');
  Description: Retrieves Array object from the database based class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : Throws if nor class passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_class {
    my ($self, $class) = @_;

    throw("Must provide and array class e.g.'AFFY_ST'") if (! defined $class);
    return $self->generic_fetch("a.class='".uc($class)."'");
}


=head2 fetch_all_by_type

  Arg [1]    : List of strings - type(s) (e.g. OLIGO, PCR)
  Example    : my @arrays = @{$array_adaptor->fetch_all_by_type('OLIGO')};
  Description: Fetch all arrays of a particular type.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : Throws if no type is provided
  Caller     : General
  Status     : at risk - needs deprecating and changing to fetch_all_by_types

=cut

#Is this used in the API anywhere?

sub fetch_all_by_type {
  my ($self, @types) = @_;

  throw('Need type as parameter') if ! @types;

  my $constraint;
  if (scalar @types == 1) {
    $constraint = qq( a.type = '$types[0]' );
  } else {
    $constraint = join q(','), @types;
    $constraint = qq( a.type IN ('$constraint') );
  }

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Experiment

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my @arrays = @{$array_adaptor->fetch_all_by_Experiment($exp)};
  Description: Fetch all arrays associated with a given Experiment
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : none
  Caller     : General
  Status     : Stable

=cut

#Changed to use simple query extension
#Removed 2 queries
#96.4% or 26.7 times faster!!!

sub fetch_all_by_Experiment{
  my ($self, $exp) = @_;

 $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exp);

  #Extend query tables
  $self->_tables([['array_chip', 'ac'], ['experimental_chip', 'ec']]);

  #Extend query and group
  my $arrays = $self->generic_fetch($exp->dbID.'=ec.experiment_id and ec.array_chip_id=ac.array_chip_id and ac.array_id=a.array_id GROUP by a.array_id');

  $self->reset_true_tables;

  return $arrays;
}


=head2 fetch_all_by_ProbeSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::ProbeSet
  Example    : my @arrays = @{$aa->fetch_all_by_ProbeSet($probeset)};
  Description: Fetch all arrays containing a given ProbeSet
               This is a convenience method to hide the 2 adaptors required
               for this call.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : none
  Caller     : General
  Status     : at risk

=cut


#Changed to use simple query extension
#Removed 1 query and hash loop
#This is only 1.04 times faster or ~ 4%

sub fetch_all_by_ProbeSet{
  my ($self, $pset) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ProbeSet', $pset);

  #Extend query tables
  $self->_tables([['array_chip', 'ac'], ['probe', 'p']]);

  #Extend query and group
  my $arrays =  $self->generic_fetch('p.probe_set_id='.$pset->dbID.' and p.array_chip_id=ac.array_chip_id and ac.array_id=a.array_id GROUP BY a.array_id');
  # ORDER BY NULL');#Surpresses default order by group columns. Actually slower? Result set too small?

  $self->reset_true_tables;

  return $arrays;
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
  return (['array', 'a']);
}


=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Stable

=cut

sub _columns {
  return qw( a.array_id a.name a.format a.vendor a.description a.type a.class );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (@result, $array_id, $name, $format, $vendor, $description, $type, $class);

  $sth->bind_columns(\$array_id, \$name, \$format, \$vendor, \$description, \$type, \$class);

  while ( $sth->fetch() ) {

    my $array = Bio::EnsEMBL::Funcgen::Array->new
	  (
	   -dbID        => $array_id,
	   -adaptor     => $self,
	   -name        => $name,
	   -format      => $format,
	   -vendor      => $vendor,
	   -description => $description,
	   -type        => $type,
	   -class       => $class,
	  );

    push @result, $array;


  }
  return \@result;
}




=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::Array objects
  Example    : $oaa->store($array1, $array2, $array3);
  Description: Stores given Array objects in the database. This
               method checks for arrays previously stored and updates
               and new array_chips accordingly.
  Returntype : Listref of stored Array objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


#This works slightly differently as arary_chip are not stored as objects,
#yet we need to retrieve a dbID for the array before we know about all the array_chips

sub store {
  my $self = shift;
  my @args = @_;

  my ($sarray);

  my $sth = $self->prepare("
			INSERT INTO array
			(name, format, vendor, description, type, class)
			VALUES (?, ?, ?, ?, ?, ?)");


  foreach my $array (@args) {
    if ( !$array->isa('Bio::EnsEMBL::Funcgen::Array') ) {
      warning('Can only store Array objects, skipping $array');
      next;
    }

    if (!( $array->dbID() && $array->adaptor() == $self )){
      #try and fetch array here and set to array if valid
      $sarray = $self->fetch_by_name_vendor($array->name(), $array->vendor());#this should be name_vendor?

      if( ! $sarray){
		#sanity check here
		throw("Array name must not be longer than 30 characters") if (length($array->name) > 40);
		$sth->bind_param(1, $array->name(),         SQL_VARCHAR);
		$sth->bind_param(2, $array->format(),       SQL_VARCHAR);
		$sth->bind_param(3, $array->vendor(),       SQL_VARCHAR);
		$sth->bind_param(4, $array->description(),  SQL_VARCHAR);
		$sth->bind_param(5, $array->type(),         SQL_VARCHAR);
		$sth->bind_param(6, $array->class(),        SQL_VARCHAR);


		$sth->execute();
		$array->dbID($self->last_insert_id);
		$array->adaptor($self);
      }
      else{
		#warn("Array already stored, using previously stored array\n");# validating array_chips\n");
		$array = $sarray;
      }
    }
  }

  return \@args;
}




=head2 fetch_probe_count_by_Array

  Args       : None
  Example    : my $probe_count = @{$array_adaptor->fetch_probe_count_by_Array($array)};
  Description: Counts probes on given array
  Returntype : ints
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_probe_count_by_Array{
  my ($self, $array) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Array', $array);

  my ($count) = @{$self->dbc->db_handle->selectrow_arrayref('select count(distinct probe_id) from array_chip ac, probe p where ac.array_id='.$array->dbID.' and ac.array_chip_id=p.array_chip_id')};

  return $count;
}


=head2 fetch_Probe_dbIDs_by_Array

Arg [1]    : Bio::EnsEMBL::Funcgen::Array
Example    : my @dbids = @{$array_adaptor->fetch_Probe_dbIDs_by_Array($array)}
Description: Fetches a arrayref of Probe dbIDs for a given Array
Returntype : arrayref of Probe dbIDs
Exceptions : None
Caller     : General
Status     : at risk

=cut

sub fetch_Probe_dbIDs_by_Array{
  my ($self, $array) = @_;
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Array', $array);

  my $sql = sprintf(qq/
    SELECT
      distinct p.probe_id
    FROM
      probe p
    WHERE
      array_chip_id in( %s ) /, join( ',', @{$array->get_array_chip_ids} ) );

	my $sth = $self->prepare( $sql );
  $sth->execute || die ($sth->errstr);
  my @dbids =  map{$_->[0]} @{$sth->fetchall_arrayref};

  return \@dbids;
}

=head2 fetch_Probe_name2dbID_by_Array

Arg [1]    : Bio::EnsEMBL::Funcgen::Array
Example    : my %name2dbid = %{$array_adaptor->fetch_Probe_name2dbID_by_Array($array)}
Description: Fetches a hashref of Probe dbIDs keyed by probe name
             for all Probes of a given Array
Returntype : Hashref for Probe dbIDs keyed by probe name
Exceptions : None
Caller     : General
Status     : at risk

=cut

sub fetch_Probe_name2dbID_by_Array{
  my ($self, $array) = @_;
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Array', $array);

  my $sql = sprintf(qq/
SELECT p.name, p.probe_id
FROM   probe p
WHERE  array_chip_id in( %s ) /, join( ',', @{$array->get_array_chip_ids} ) );

  my $sth = $self->prepare( $sql );
  $sth->execute || die ($sth->errstr);
  my %mapping;
  map{$mapping{$_->[0]}=$_->[1]} @{$sth->fetchall_arrayref};
  return \%mapping;
}



sub check_status_by_class{
  my ($self, $status, $class) = @_;

  foreach my $array(@{$self->fetch_all_by_class($class)}){

	foreach my $ac(@{$array->get_ArrayChips}){

	  if(! $ac->has_status($status)){
		throw('Found '.$class.' ArrayChip '.$ac->name." without $status status");
	  }
	}
  }

  return;
}



1;

