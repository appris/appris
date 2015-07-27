#
# Ensembl module for Bio::EnsEMBL::Funcgen::ExperimentalGroup
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

  Please email comments or questions to the public Ensembl  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::ExperimentalGroup - A module to represent
an ExperimentalGroup. i.e. the authors of an experiment.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExperimentalGroup;



=head1 DESCRIPTION

This is a simple class to represent information about an ExperimentalGroup,
containing a name and a more detailed description
This module is part of the Ensembl project: http://www.ensembl.org/

=cut


package Bio::EnsEMBL::Funcgen::ExperimentalGroup;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-name]: string - name of ExperimentalGroup
  Arg [-location]: (optional) string - location of ExperimentalGroup
  Arg [-contact]: (optional) string - contact of ExperimentalGroup
  Arg [-url]: (optional) string - url containing information for the ExperimentalGroup
  Arg [-description]: (optional) string - descriptiom of ExperimentalGroup
  Arg [-project]: (optional) boolean - True if this is part of a large project (eg. ENCODE)
  Example    : my $group = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
                                                               -name  => "EBI",
                                                               -location  => "Hinxton",
                                                               -contact  => "http://lists.ensembl.org/mailman/listinfo/dev",
                                                               -url  => "http://www.ebi.ac.uk/",
                                                               -description => "European Bioinformatics Institute",
                                                               -is_project => 0,
                                                             );
  Description: Constructor method for ExperimentalGroup class
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalGroup
  Exceptions : Throws if name not defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;

  my $obj_class = ref($caller) || $caller;
  my $self = $obj_class->SUPER::new(@_);

  my (
      $name,
      $location,
      $contact,
      $url,
      $desc,
      $is_project
     ) = rearrange([
       'NAME', 'LOCATION', 'CONTACT', 'URL', 'DESCRIPTION', 'IS_PROJECT'
          ], @_);


  if($name){
    $self->name($name);
  }else{
    throw("Must supply a Group name\n");
  }
  $self->location($location) if $location;
  $self->contact($contact) if $contact;
  $self->url($url) if $url;
  $self->description($desc) if $desc;
  $self->is_project($is_project) if $is_project;

  return $self;
}



=head2 name

  Arg [1]    : string - name
  Example    : my $name = $ft->name();
  Description: Getter and setter of name attribute for ExperimentalGroup objects
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub name {
    my $self = shift;
    $self->{'name'} = shift if @_;
    return $self->{'name'};
}

=head2 description

  Arg [1]    : (optional) string - description
  Example    : my $desc = $group->description();
  Description: Getter and setter of description attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if @_;
    return $self->{'description'};
}

=head2 location

  Arg [1]    : (optional) string - location
  Example    : my $location = $group->location();
  Description: Getter and setter of location attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub location {
    my $self = shift;
    $self->{'location'} = shift if @_;
    return $self->{'location'};
}


=head2 contact

  Arg [1]    : (optional) string - contact
  Example    : my $contact = $group->contact();
  Description: Getter and setter of contact attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub contact {
    my $self = shift;
    $self->{'contact'} = shift if @_;
    return $self->{'contact'};
}


=head2 url

  Arg [1]    : (optional) string - url
  Example    : my $url = $group->url();
  Description: Getter and setter of url attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub url {
    my $self = shift;
    $self->{'url'} = shift if @_;
    return $self->{'url'};
}

=head2 is_project

  Arg [1]    : (optional) Boolean - is_project
  Example    : $group->is_project();
  Description: Getter and setter of is_project attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : High Risk

=cut

sub is_project {
    my $self = shift;
    $self->{'is_project'} = shift if @_;
    return $self->{'is_project'};
}


=head2 reset_relational_attributes

  Arg[1]     : UNDEF - placeholder for compliance with similar methods in
               other modules
  Arg[2]     : Flag to avoid reseting of adaptor and dbID

  Description: Undefs Adaptor and dbID
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;

  # undef the dbID and adaptor
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  return;
}


=head2

Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
Args[2]    : Boolean - Optional 'shallow' - no object methods compared
Args[3]    : Arrayref - Optional list of InputSubset method names each
             returning a Scalar or an Array or Arrayref of Scalars.
             Defaults to: name location contact description url is_project
Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
Description: Compare this ExperimentalGroup to another based on the defined scalar
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

  $scl_methods ||= [qw(name location contact description url is_project)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods, $obj_methods);
}

1;

