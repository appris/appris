#
# Ensembl module for Bio::EnsEMBL::Funcgen::Utils::Encode
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

Bio::EnsEMBL::Funcgen::Utils::Encode - provides some handy methods to deal 
with Encode regions.

=head1 SYNOPSIS

    use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);

    my $encode_regions = &get_encode_regions($dnadb, $assembly_version);

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Utils::Encode;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use POSIX;
use Data::Dumper;

use base qw( Exporter );

use vars qw(@EXPORT_OK);
@EXPORT_OK = qw(get_encode_regions);

=head2 get_encode_regions

  Arg [1]     : Bio::EnsEMBL::DBSQL::DBAdaptor
  Description : fetch encode region coordinates from core db
  Returns     : reference to list of Bio::EnsEMBL::MiscFeatures Encode region slices
  Exceptions  : Throws if no Bio::EnsEMBL::DBSQL::DBAdaptor defined
  Example     : get_encode_regions($dnadb);

=cut

sub get_encode_regions {

    my $db = shift;
    
    throw("Need to pass a valid Bio::EnsEMBL::DBSQL::DBAdaptor") 
        if (! ($db  && $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")));

    my $sa = $db->get_SliceAdaptor();
    my $tls = $sa->fetch_all('toplevel');
    #map { print Dumper $_->name } @$tls;

    my $mfa = $db->get_MiscFeatureAdaptor();

    my @encode_regions;
    map { 
        push @encode_regions, @{$mfa->fetch_all_by_Slice_and_set_code($_, 'encode')};
    } @$tls;

    return \@encode_regions;


    #my %encode_regions;
    #map { $encode_regions{$_->display_id} = sprintf
    #          ("%s:%s:%s:%d:%d:%d",
    #           $_->slice->coord_system_name, 
    #           $_->slice->coord_system()->version(),
    #           $_->slice->seq_region_name,
    #           $_->start, $_->end, $_->strand);
    #  } @encode_regions;
    #
    #return \%encode_regions;

}

1;
