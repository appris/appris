#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::GFF
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

#Could this be based on a Generic Flat file parser?

=head1 NAME

Bio::EnsEMBL::Funcgen::Parsers::GFF

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::GFF";
  push @INC, $parser_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally set by the Importer as the parent class.  GFF contains meta 
data and methods specific to data in bed format, to aid 
parsing and importing of experimental data.

=cut

package Bio::EnsEMBL::Funcgen::Parsers::GFF;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

use base qw( Bio::EnsEMBL::Funcgen::Parsers::InputSet );

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for GFF class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::GFF
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;
  
  my $class = ref($caller) || $caller;

  #define default fields here and pass
  #We also need to be able to take custom attrs mappings

  #keys are array index of field, key are Feature paramter names
  #reverse this?
  #Unless we have a compound field which we name accordingly
  #And then call e.g. parse_attrs
  #Which will return a hash with the relevant Feature attributes

  #Is splitting this up simply going to make the parse slower due to acessor methods?

  #Pass or just set directly here?
  #<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]


  #Some of these may be highly redundant due to the nature of the data.
  #We can hash things to lessen the impact but we're still going to be checking if exists for each one
  #No way around this?  Unless it is marginally faster to set a permanent type and then only check a boolean.
  #Yes there is, this is the exhaustive GFF definition, we can just redefine or delete some entries dynamically to
  #avoid ever considering a particular field index.


  #Don't need any of this? Can we simply define process fields?
  #This will remove the ability to define custom formats
  #But then again we can only have custom format if it has ensembl compliant data
  #i.e. no preprocessing has to be done before populating the feature_params hash

  #my %fields = (
#				0 => 'fetch_slice',
#				1 => 'get_source',
#				2 => 'get_feature_type',
#				3 => '-start',
	#			4 => '-end',
#				5 => '-strand',#Will most likely be , need to convert to -.+ > -1 0 1
				#6 => 'frame',#will most likely be .
#				7 => 'get_attributes',
#			   );

  #We want to be able to define mappings between attributes and fields
  #we're basically just dealing with display_label for annotated_feature
  #e.g -display_label_format => ID+ACC
  #Or maybe format of several fields and attrs + text?
  #We need a separator which will not be used in the GFF attr names
  #we also need to be able to differentiate
  #First check standard GFF field, then check attrs
  ##No no no, just have method, generate display label
  #forget this for now and just use one field

  my $display_label_field = 'ID';#default

  #We still need to define the field name here as a global hash to allow this display_label_field look up.


  my $self  = $class->SUPER::new(@_);#, -fields => \%fields);
  
  ($display_label_field) = rearrange(['DISPLAY_LABEL_FIELD'], @_);

  #We need to define meta header method, starting with '##'
  #Also need to skip comments '#' at begining or end of line
  #Do we also need to skip field header? No methinks not.

  #Define result method
 # $self->{'file_ext'} => 'gff';#Could use vendor here?
  
  #define this if we want to override the generic method in Simple
  #$self->{'config'}{'results_data'} => ["and_import_gff"];  

  $self->display_label_field($display_label_field);


  return $self;
}


=head2 set_config

  Example    : my $self->set_config;
  Description: Sets attribute dependent config
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub set_config{
  my $self = shift;

  $self->SUPER::set_config;

  #GFF specific stuff here.

  return;
}

#Need to implement this!
sub parse_line{
  my ($self, $line) = @_;

  #return if $line ~=

 #my %fields = (
#				0 => 'fetch_slice',
#				1 => 'get_source',
#				2 => 'get_feature_type',
#				3 => '-start',
	#			4 => '-end',
#				5 => '-strand',#Will most likely be , need to convert to -.+ > -1 0 1
				#6 => 'frame',#will most likely be .
#				7 => 'get_attributes',
#			   );



  my ($chr, $start, $end, $pid, $score) = split/\t/o, $line;			
  
  #we need to return feature_params and seq if defined?

}
 


1;
