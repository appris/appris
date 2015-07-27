#
# Ensembl module for Bio::EnsEMBL::Funcgen::ArrayChip
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

Bio::EnsEMBL::Funcgen::ArrayChip - A simple module to represent the concept/template of 
a chip/slide within an array,  of which the physical manifestation is an ExperimentalChip.

=head1 SYNOPSIS

 use Bio::EnsEMBL::Funcgen::ArrayChip;

 my $ec = Bio::EnsEMBL::Funcgen::ArrayChip->new(
                    							 -ARRAY_ID  => $array->dbID(),
					                    		 -NAME      => $desc,
                                                 -DESIGN_ID => $design_id,
							                   );

#add more methods here?


=head1 DESCRIPTION

An ArrayChip object represent the concept of an array chip/slide withing a given array/chipset.
The data for ArrayChips is stored in the array_chip table.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::ArrayChipAdaptor
Bio::EnsEMBL::Funcgen::Array

=cut

package Bio::EnsEMBL::Funcgen::ArrayChip;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-ARRAY_ID]  : Int - the dbID of the parent array
  Arg [-ARRAY]     : Bio::EnsEMBL::Funcgen::Array
  Arg [-DESIGN_ID] : String - the unqiue deisng ID defined by the array vendor
  Arg [-NAME]      : String - the name of the array chip


  Example    : my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new
                 (
							    -ARRAY_ID  => $array->dbID(),
							    -NAME      => $desc,
                  -DESIGN_ID => $design_id,
							   );	
  Description: Creates a new Bio::EnsEMBL::Funcgen::ArrayChip object.
  Returntype : Bio::EnsEMBL::Funcgen::ArrayChip
  Exceptions : Throws if mandatory parameters are not set.
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = $class->SUPER::new(@_);

  my ($array_id, $name,  $design_id, $array)
    = rearrange( ['ARRAY_ID', 'NAME', 'DESIGN_ID', 'ARRAY'], @_ );

   
  #Remove array_id so we can remove checking below?

  throw("Must define a name($name) and design_id($design_id)") if(! $name || ! $design_id);


  #Make these mutually exclusive to avoid checking
  if($array_id && $array){
    throw('Must provide either -array or -array_id but not both');
  }

  if(defined $array){
	
    if(!(ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array'))){
      throw('Array paramter must be a valid Bio::EnsEMBL::Funcgen::Array');
    }
	
    $self->{array} = $array;
  }

  $self->{array_id}  = $array_id;
  $self->{name}      = $name;
  $self->{design_id} = $design_id;

  return $self;
}


=head2 array_id

  Example    : my $array_id = $array_chip->array_id();
  Description: Getter array_id attribute.
  Returntype : Int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub array_id {  return shift->{array_id}; }

=head2 name

  Example    : my $ac_name = $array_chip->name();
  Description: Getter for the name attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return shift->{name}; }

=head2 design_id

  Example    : my $design_id = $array_chip->design_id();
  Description: Getter for the design_id attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub design_id {  return shift->{design_id}; }


=head2 get_Array

  Example    : my $array = $array_chip->get_array();
  Description: Getter for the array attribute
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Array {
  my $self = shift;

  if(! defined $self->{'array'}){
    $self->{'array'} = $self->adaptor->db->get_ArrayAdaptor()->fetch_by_dbID($self->array_id());
  }

  return $self->{'array'};
}



1;

