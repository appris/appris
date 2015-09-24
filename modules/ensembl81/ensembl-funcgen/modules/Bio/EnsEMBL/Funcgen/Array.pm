#
# Ensembl module for Bio::EnsEMBL::Funcgen::Array
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

Bio::EnsEMBL::Funcgen::Array - A module to represent a nucleotide microarray.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Array;

my $array = Bio::EnsEMBL::Funcgen::Array->new
 (
	-NAME        => 'Array-1',
  -FORMAT      => 'Tiled',
  -SIZE        => '1',
  -VENDOR      => 'Nimblegen',
  -DESCRIPTION => $desc,
  -TYPE        => 'OLIGO',
  -CLASS       => 'VENDOR_FORMAT'
 );

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $array_adaptor = $db_adaptor->get_ArrayAdaptor();
my $array = $array_adaptor->fetch_by_name($array_name)

=head1 DESCRIPTION

An Array object represents a nucleotide (OLIGO, PCR etc.) microarray. The data
(currently the name, format, size, species, vendor and description) are stored
in the array table.

=cut


package Bio::EnsEMBL::Funcgen::Array;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );

use base qw(Bio::EnsEMBL::Funcgen::Storable);


# Possible types for OligoArray objects
#This should match the vendor enum values?
#%VALID_TYPE = (
#	'AFFY'  => 1,
#	'OLIGO' => 1,
#);


=head2 new

  Arg [-NAME]        : String - the name of this array
  Arg [-VENDOR]      : String - the vendor of this array (AFFY, NIMBLEGEN etc)
  Arg [-TYPE]        : String - type of array e.g. OLIGO, PCR
  Arg [-FORMAT]      : String - the format of this array (TILED, TARGETTED, GENE etc)
  Arg [-DESCRIPTION] : String - description of the array 

  Example    : my $array = Bio::EnsEMBL::Funcgen::Array->new
                (
								 -NAME        => 'Array-1',
								 -FORMAT      => 'Tiled',
								 -VENDOR      => 'Nimblegen',
                 -TYPE        => 'OLIGO',
								 -DESCRIPTION => $desc,
                 -CLASS       => 'VENDOR_FORMAT',#e.g. AFFY_UTR, ILLUMINA_WG
							  );
  Description: Creates a new Bio::EnsEMBL::Funcgen::Array object.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : Throws if mandatory params not set/valid
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  
  my ($name, $format, $vendor, $type, $desc, $aclass)
    = rearrange( ['NAME', 'FORMAT', 'VENDOR', 'TYPE', 'DESCRIPTION', 'CLASS'], @_ );
  
  #mandatory params?
  #name, format, vendor
  #enum on format?

  my @stack = caller();

  if($self->dbID() && $stack[0] ne "Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor"){
    throw("You must use the ArrayAdaptor($stack[0]) to generate Arrays with a dbID i.e. from the DB, as this module accomodates updating which may cause incorrect data if the object is not generated form the DB");
  } 


  throw("Must provide a vendor parameter") if ! $vendor;
  throw("Must provide a name parameter")   if ! $name;
   
  $self->name($name);
  $self->format($format)    if defined $format;

  if(defined $format && $format eq 'EXPRESSION' && ! defined $class){
	throw('You must defined a class if you are importing and array with an EXPRESSION format');
  }

  $self->class(uc($aclass))     if defined $aclass;
  $self->vendor($vendor);
  $self->description($desc) if defined $desc;
  $self->type($type)        if defined $type;
  
  return $self;
}

=head2 get_all_Probes

  Args       : None
  Example    : my $probes = $array->get_all_Probes();
  Description: Returns all probes on an array. Needs a database connection.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Probe objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_all_Probes {
	my $self = shift;

	if ( $self->dbID() && $self->adaptor() ) {
		my $opa = $self->adaptor()->db()->get_ProbeAdaptor();
		my $probes = $opa->fetch_all_by_Array($self);
		return $probes;
	} else {
		warning('Need database connection to retrieve Probes');
		return [];
	}
}

=head2 get_all_Probe_dbIDs

  Args       : None
  Example    : my @dbids = @{$array->get_all_Probe_dbIDs};
  Description: Returns an array ref of all the Probe database IDs for this array
  Returntype : arrayref of ints
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_all_Probe_dbIDs {
  my $self = shift;

  if(!  $self->{probe_dbids}){
	#check for adaptor here?
	
	if(! $self->adaptor){
	  throw('Must have set an adaptor to get_all_Probe_dbIDs');
	}
	
	$self->{probe_dbids} = $self->adaptor->fetch_Probe_dbIDs_by_Array($self);
  }

  return  $self->{probe_dbids};
}


=head2 get_all_ProbeSets

  Args       : None
  Example    : my $probesets = $array->get_all_ProbeSets();
  Description: Returns all probesets on an array. Needs a database connection.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeSets objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_ProbeSets {
	my $self = shift;

	if ( $self->dbID() && $self->adaptor() ) {
		my $opsa = $self->adaptor()->db()->get_ProbeSetAdaptor();
		my $probesets = $opsa->fetch_all_by_Array($self);
		return $probesets;
	} else {
		warning('Need database connection to retrieve ProbeSets');
		return [];
	}
}


#All the array_chip methods will be migrated to ArrayChip.pm

=head2 get_array_chip_ids

  Example    : my @ac_ids = @{$array->get_array_chip_ids()};
  Description: Returns all array_chip_ids for this array.
  Returntype : Listref of array_chip ids
  Exceptions : Throws if none retrieved
  Caller     : General
  Status     : At Risk

=cut

sub get_array_chip_ids {
  my $self = shift;

  my @ac_ids;


  $self->get_ArrayChips();

  #should we get_ArrayChips is we have none cached?
  #this may cause problem


  foreach my $achip(values %{$self->{'array_chips'}}){
    push @ac_ids, $achip->dbID();
  }

  if(! @ac_ids){
    throw("No array_chip_ids available"); # should this be warn?
  }
  
  return \@ac_ids;
}

=head2 get_design_ids

  Example    : my @design_ids = @{$array->get_design_ids()};
  Description: Returns a the design_ids for each array_chip contained within this array
  Returntype : list
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut



sub get_design_ids{
  my $self = shift;
  return [keys %{$self->{'array_chips'}}];
}
    


=head2 name

  Arg [1]    : (optional) string - the name of this array
  Example    : my $name = $array->name();
  Description: Getter, setter of the name attribute for Array
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub name{
  my $self = shift;
  
  $self->{'name'} = shift if @_;
  
  #do we need this?
  #if ( !exists $self->{'name'} && $self->dbID() && $self->adaptor() ) {
  #  $self->adaptor->fetch_attributes($self);
  #}
  
  return $self->{'name'};
}


=head2 type

  Arg [1]    : (optional) string - the type of this array
  Example    : $array->type('OLIGO');
  Description: Getter, setter of the type attribute for Array
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub type{
  my $self = shift;
  
  $self->{'type'} = shift if @_;
    
  return $self->{'type'};
}


=head2 format

  Arg [1]    : (optional) string - the format of the array
  Example    : my $format = $array->format();
  Description: Getter, setter of format attribute for
               Array objects e.g. Tiled, Targetted etc...
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub format {
  my $self = shift;
  
  $self->{'format'} = shift if @_;
  
  #do we need this?
  #if ( !exists $self->{'format'} && $self->dbID() && $self->adaptor() ) {
  #  $self->adaptor->fetch_attributes($self);
  #}
  
  return $self->{'format'};
}

=head2 class

  Arg [1]    : (optional) string - the class of the array
  Example    : my $class = $array->class('AFFY_UTR');
  Description: Getter, setter of class attribute for
               Array objects e.g. AFFY_UTR, AFFY_ST
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub class {
  my $self = shift;
  
  $self->{'class'} = shift if @_;
  
  return $self->{'class'};
}


=head2 vendor

  Arg [1]    : (optional) string - the name of the array vendor
  Example    : my $vendor = $array->vendor();
  Description: Getter, setter of vendor attribute for
               Array objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub vendor {
  my $self = shift;
  $self->{'vendor'} = shift if @_;
  
  #do we need this?
  #if ( !exists $self->{'vendor'} && $self->dbID() && $self->adaptor() ) {
  #  $self->adaptor->fetch_attributes($self);
  #}

  return $self->{'vendor'};
}

=head2 description

  Arg [1]    : (optional) string - the description of the array
  Example    : my $desc = $array->description();
  Description: Getter, setter of description attribute for
               Array objects. 
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub description {
  my $self = shift;
  $self->{'description'} = shift if @_;
  
  #do we need this?
  #if ( !exists $self->{'description'} && $self->dbID() && $self->adaptor() ) {
  #  $self->adaptor->fetch_attributes($self);
  #}

  return $self->{'description'};
}

=head2 probe_count

  Example    : my $num_probes = $array->probe_count();
  Description: Return number of probes on array
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub probe_count {
  my ($self)  = @_;
  #Do we want a distinct flag here?

  if(! defined $self->{'probe_count'}){
	$self->{'probe_count'} = $self->adaptor->fetch_probe_count_by_Array($self);
  }
  
  return $self->{'probe_count'};
}



=head2 get_ArrayChips

  Example    : my @achips = @{$array->get_ArrayChips()};
  Description: Getter, setter and lazy loader of array_chip hashes
  Returntype : Arrays of ArrayChip objects
  Exceptions : Throws exception if none found for array_id
  Caller     : General
  Status     : High Risk - migrate to ArrayChip.pm

=cut

sub get_ArrayChips {
  my $self = shift;
 
  #lazy loaded as we won't want this for light DB
  #should do meta check and want here

  if ( ! exists $self->{'array_chips'}){

    if( $self->dbID() && $self->adaptor() ) {
      #$self->adaptor->fetch_attributes($self);
      #need to do this differently as we're accessing a different table
      $self->{'array_chips'} = {};

      foreach my $achip(@{$self->adaptor->db->get_ArrayChipAdaptor->fetch_all_by_array_id($self->dbID())}){
	$self->{'array_chips'}{$achip->design_id} = $achip;
	#%{$self->{'array_chips'}} = %{$self->adaptor->db->get_ArrayAdaptor->_fetch_array_chips_by_array_dbID($self->dbID())};
      }
    }
    else{
      throw("Need array dbID and DB connection to retrieve array_chips");
    }
  }

  return [ values %{$self->{'array_chips'}} ];
}

=head2 get_ArrayChip_by_design_id

  Arg [1]    : (mandatory) int - design_id
  Example    : my %ac = %{$array->get_ArrayChip_by_design_id('1234')};
  Description: Getter for array_chip hashes
  Returntype : Hashref
  Exceptions : Throws exception if no design_id defined, warns if not part of array
  Caller     : General
  Status     : At risk

=cut

sub get_ArrayChip_by_design_id{
  my ($self, $design_id) = @_;


  #warn "This needs to get the array chip if not defined?? but we're using it to test whether is has been stored same problem as probe_design?";

  my ($achip);
  throw("Must supply a valid array chip design_id") if (! defined $design_id);

  if(defined $self->{'array_chips'}{$design_id}){
    $achip = $self->{'array_chips'}{$design_id};
  }else{
    #No we use this to check whether it has been stored with the array
    #warn("should this throw? Array does not contain ArrayChip:$design_id\n"); 
  }

  return $achip;
}


=head2 add_ArrayChip

  Arg [1]    : mandatory - Bio::EnsEMBL::Funcgen::ArrayChip
  Example    : $array->add_ArrayChip($array_chip);
  Description: Setter for array chips
  Returntype : None
  Exceptions : Throws if arg not a Bio::EnsEMBL::Funcgen::ArrayChip, or Array not stored
  Caller     : General
  Status     : Ar risk

=cut

sub add_ArrayChip{
  my ($self, $array_chip) = @_;

  if ($self->adaptor){
    $self->adaptor->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ArrayChip', $array_chip, 'array_chip');
    $self->get_ArrayChips if ! $self->{array_chips};

    if(! exists $self->{array_chips}{$array_chip->design_id}){
      $self->{'array_chips'}{$array_chip->design_id} = $array_chip;
    }
  }
  else{
    throw('Array needs an adaptor before adding an array_chip');
  }

  return;
}


### DEPRECATED METHODS ###



1;

