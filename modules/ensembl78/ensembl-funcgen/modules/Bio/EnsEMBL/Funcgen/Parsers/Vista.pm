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

=cut

package Bio::EnsEMBL::Funcgen::Parsers::Vista;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw( Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser );

# Parse data from LBL enhancers, see http://enhancer.lbl.gov/cgi-bin/imagedb.pl?show=1;search.result=yes;form=search;search.form=no;action=search;search.sequence=1
# e.g.
#
# >chr16:84987588-84988227 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]
# AACTGAAGGGACCCCGTTAGCATAtaaacaaaaggtggggggtagccccgagcctcttct
# ctgacagccagtggcggcagtgatgaatttgtgaagttatctaattttccactgttttaa
# ttagagacttgggctctgaggcctcgcagctggcttctttgtgctgtattctgttgcctg
# acagag

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_, type => 'Vista');

  #Set default feature_type and feature_set config
  $self->{static_config}{feature_types} = {(
                                            'VISTA Target'   => {
                                                                 -name        => 'VISTA Target',
                                                                 -class       => 'Search Region',
                                                                 -description => 'VISTA target region',
                                                                },
                                            'VISTA Enhancer' => {
                                                                 -name         => 'VISTA Enhancer',
                                                                 -class        => 'Enhancer',
                                                                 -description  => 'Enhancer identified by positive VISTA assay',
                                                                 -so_name      => 'enhancer',
                                                                 -so_accession => 'SO:0000165',
                                                                },
                                            'VISTA Target - Negative' => {
                                                                          -name        => 'VISTA Target - Negative',
                                                                          -class => 'Search Region',
                                                                          -description => 'Enhancer negative region identified by VISTA assay',
                                                                         },
                                           )};

  $self->{static_config}{analyses} = {
                                      VISTA =>  {
                                                 -logic_name => 'VISTA',
                                                 -description   => 'VISTA Enhancer Assay (http://enhancer.lbl.gov/)',
                                                 -display_label => 'VISTA',
                                                 -displayable   => 1,
                                                },
                                     };

  #This is used as the entry point to store/validate
  #So all of the above needs to be referenced in here
  $self->{static_config}{feature_sets} = {
                                          'VISTA enhancer set' =>
                                          {
                                           #Stored in this order

                                           #Entries here are flexible
                                           #Can be omited if defined in feature_set
                                           #top level analyses/feature_types definition required if no DB defaults available
                                           #These can be a ref to the whole or subset of the top level analyses/feature_types hash
                                           #A key with an empty hash or undef(with or without a matching key in the top level analyses/feature_types hash

                                           #analyses      => $self->{static_config}{analyses},
                                           feature_types => $self->{static_config}{feature_types},

                                           #feature_type and analysis values must be string key to top level hash
                                           #This wont work for feature_types as they are not unique by name!!!!!!
                                           #This is why we have top level hash where we can define a unique compound key name

                                           feature_set   =>
                                           {
                                            -feature_type      => 'VISTA Target', #feature_types config key name not object
                                            -display_label     => 'VISTA Enhancers',
                                            -description       => 'Experimentally validated enhancers',
                                            -analysis          => 'VISTA', #analyses config key name not object
                                           },
                                          }
                                         };

  #$self->validate_and_store_feature_types;
  $self->validate_and_store_config([keys %{$self->{static_config}{feature_sets}}]);
  $self->set_feature_sets;

  return $self;
}



# Parse file and return hashref containing:
#
# - arrayref of features
# - arrayref of factors



sub parse_and_load{
  my ($self, $files, $old_assembly, $new_assembly) = @_;

  if (scalar(@$files) != 1) {
    throw('You must provide a unique file path to load VISTA features from:\t'.join(' ', @$files));;
  }

  my $file = $files->[0];
  $self->log_header("Parsing and loading LBNL VISTA enhancer data from:\t$file");

  my $extfeat_adaptor = $self->db->get_ExternalFeatureAdaptor;
  my $fset_config      = $self->{static_config}{feature_sets}{'VISTA enhancer set'};
  my $feature_positive = $fset_config->{'feature_types'}{'VISTA Enhancer'};
  my $feature_negative = $fset_config->{'feature_types'}{'VISTA Target - Negative'};
  my $set              = $fset_config->{feature_set};

  use Bio::EnsEMBL::Registry;
  my %id_prefixes = (
                     homo_sapiens => 'hs',
                     mus_musculus => 'mm',
                    );

  my $species = Bio::EnsEMBL::Registry->get_alias($self->db->species);

  if ( (! defined $species) ||
       (! exists $id_prefixes{$species}) ) {
    throw("Failed to get a VISTA ID prefix for species alias:\t$species");
  }

  $species = $id_prefixes{$species};

  ### Read file
  open (FILE, "<$file") || die "Can't open $file";
  my $cnt = 0;
  my $skipped = 0;


  while (<FILE>) {

    next if ($_ !~ /^>/o);      # only read headers

    # OLD >chr16:84987588-84988227 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]
    # from v66 >Mouse|chr12:112380949-112381824 | element 3 | positive  | neural tube[4/4] | hindbrain (rhombencephalon)[4/4] | forebrain[4/4]

    #ID naming scheme change from LBNL-1 to hs1 or mm1
    #But the flat file and url use two different naming schemes!
    #VISTA URL is: where experiment id is element number and species_id 1 = human and 2 = mouse
    #http://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=presentation&show=1&experiment_id=1&organism_id=1

    #Add links to cell_type for tissues in @expression_pattern?
    #This would be vista specific cell_type_annotation? Or we could just have associated_cell_type? (without annotation)
    #Just link for now

    my (undef, $coords, $element, $posneg, @expression_patterns) = split /\s*\|\s*/o; #was \s+
    # parse co-ordinates & id
    my ($chr, $start, $end) = $coords =~ /chr([^:]+):(\d+)-(\d+)/o;
    my ($element_number) = $element =~ /\s*element\s*(\d+)/o;

    # seq_region ID and co-ordinates
    my $chr_slice;

    if ($old_assembly) {
      $chr_slice = $self->slice_adaptor->fetch_by_region('chromosome', $chr, undef, undef, undef, $old_assembly);
    } else {
      $chr_slice = $self->slice_adaptor->fetch_by_region('chromosome', $chr);
    }

    if (!$chr_slice) {
      warn "Can't get slice for chromosme $chr\n";
      next;
    }

    my $seq_region_id = $chr_slice->get_seq_region_id;
    throw("Can't get seq_region_id for chromosome $chr") if (!$seq_region_id);

    # Assume these are all on the positive strand? Is this correct?
    my $strand = 1;


    my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
      (
       -start         => $start, #is this in UCSC coords?
       -end           => $end,  #is this in UCSC coords?
       -strand        => $strand,
       -feature_type  => $posneg eq 'positive' ? $feature_positive : $feature_negative,
       -slice         => $self->slice_adaptor->fetch_by_region('chromosome', $chr, undef, undef, $strand, $old_assembly),
       -display_label => $species.$element_number, #"LBNL-$element_number",
       -feature_set   => $set,
      );


    # project if necessary
    if ($new_assembly) {

      $feature = $self->project_feature($feature, $new_assembly);

      if (! defined $feature) {
        $skipped ++;
        next;
      }
    }

    $cnt ++;
    $extfeat_adaptor->store($feature);
  }

  close FILE;

  $self->log('Parsed '.($cnt+$skipped).' features');
  $self->log("Loaded $cnt features");
  $self->log("Skipped $skipped features");

  #Now set states
  foreach my $status (qw(DISPLAYABLE MART_DISPLAYABLE)) {
    $set->adaptor->store_status($status, $set);
  }


  return;

}


1;
