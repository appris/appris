#
# Ensembl module for Bio::EnsEMBL::Funcgen::CellType
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

Bio::EnsEMBL::Funcgen::CellType - A module to represent a CellType.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::CellType;

#Fetch from adaptor
my $ctype = $cell_type_adaptor->fetch_by_name($ctype_name);

#Create from new
my $ctype = Bio::EnsEMBL::Funcgen::CellType->new
  (
   -name          => 'H1-ESC',
   -display_label => 'H1-ESC',
   -description   => 'Human Embryonic Stem Cell',
   -efo_id        => 'efo:EFO_0003042',
   -tissue        => 'embryonic stem cell',
  );

print $ctype->name.' is a '.$ctype->description."\n";

#H1-ESC is a Human Embryonic Stem Cell


=head1 DESCRIPTION

This is a simple class to represent information about a cell type.  This may represent
harvested cells, a cell line or a more generic tissue type.

=cut

package Bio::EnsEMBL::Funcgen::CellType;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use parent qw(Bio::EnsEMBL::Funcgen::Storable);

my %valid_genders = (male   => 1,
                     female => 1,
                     hermaphrodite => 1,
                     mixed => 1 
                    );

=head2 new

  Arg [1]    : String - name of CellType
  Arg [2]    : String (optional) - display label of CellType. Defaults to name
  Arg [3]    : String (optional) - description of CellType
  Arg [4]    : String (optional) - gender e.g. male, female or hermaphrodite
  Arg [5]    : String (optional) - Experimental Factor Ontology ID e.g. efo:EFO_0002869

  Example              : my $ct = Bio::EnsEMBL::Funcgen::CellType->new
                                    (
                                     -name          => "U2OS",
                                     -display_label => "U20S",
                                     -description   => "Human Bone Osteosarcoma Epithelial Cells",
                                     -gender        => 'female',
                                     -efo_id        => 'efo:EFO_0002869',
                                    );

  Description: Constructor method for CellType class
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : Throws if name not defined or gender is invalid
  Caller     : General
  Status     : Stable

=cut

#-type/class => "TISSUE", enum? Mandatory.
#remove display label?

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($name, $dlabel, $desc, $gender, $efo_id, $tissue) = rearrange
    (['NAME', 'DISPLAY_LABEL', 'DESCRIPTION','GENDER', 'EFO_ID', 'TISSUE'], @_);

  throw("Must supply a CellType name") if ! defined $name;

  if (defined $gender) {

    if ( ! exists $valid_genders{lc($gender)} ) { #enum will not force this so validate here
      throw("Gender not valid, must be one of:\t".join(' ', keys %valid_genders));
    }

    $self->{gender} = $gender;
  }

  #Set explicitly to enable faster getter only methods
  $self->{name}          = $name;
  $self->{display_label} = $dlabel || $name;
  $self->{description}   = $desc   if defined $desc;
  $self->{efo_id}        = $efo_id if defined $efo_id;
  $self->{tissue}        = $tissue if defined $tissue;

  return $self;
}

=head2 name

  Example    : my $name = $ct->name;
  Description: Getter of name attribute for CellType objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return $_[0]->{name}; }


=head2 gender

  Example    : my $gender = $ct->gender();
  Description: Getter for the gender attribute for CellType objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub gender {  return $_[0]->{gender}; }


=head2 description

  Example    : my $desc = $ct->description();
  Description: Getter of description attribute for CellType objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub description {  return $_[0]->{description}; }


=head2 display_label

  Example    : my $display_label = $ct->display_label();
  Description: Getter of display_label attribute for CellType objects.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {  return $_[0]->{display_label}; }


=head2 efo_id

  Example    : my $efo_id = $ft->efo_id;
  Description: Getter of the Experimental Factor Ontology ID
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub efo_id{  return $_[0]->{efo_id}; }


=head2 tissue

  Example    : my $tissue = $ft->tissue;
  Description: Getter of the tissue attribute for a given cell type
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub tissue{   return $_[0]->{tissue}; }

=head2 reset_relational_attributes

  Arg[1]     : Hashref containing linked objects - KEPT for compliance
               with similar methods in other modules
  Arg[2]     : Flag to avoid reseting of adaptor and dbID

  Description: Undefs Adaptor and dbID
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : None
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = shift;

  # undef the dbID and adaptor
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  return;
}

=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of CellType method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: name table_name feature_class get_all_states
  Args[4]    : Arrayref - Optional list of CellType method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: feature_type cell_type analysis get_support
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this CellType to another based on the defined scalar
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

  $scl_methods ||= [qw(name display_label description gender efo_id tissue)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods,
                                  $obj_methods);
}


1;

