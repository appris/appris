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

package Bio::EnsEMBL::Funcgen::Parsers::biotiffin;

use strict;
use warnings;
use File::Basename;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Funcgen::ExternalFeature;

use base qw( Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser );


# To get files for bioTIFFIN, download the following GFF file (e.g. via wget):
#
#  http://td-blade.gurdon.cam.ac.uk/tad26/fly-tiffinScan-tiffin12.dm3.gff.gz

# Thomas Down <thomas.down@gurdon.cam.ac.uk>

# 
# 3R      MotifScanner    TIFDMEM0000001  936391  936401  0.0     +       0
# 3R      MotifScanner    TIFDMEM0000001  13455911        13455921        0.0     -       0
# 3R      MotifScanner    TIFDMEM0000001  17062830        17062840        0.0     +       0



sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_, type => 'BioTiffin');

  #Set default feature_type and feature_set config

  #We need to capture version/release/data of external feature sets.
  #This can be nested in the description?  Need to add description to feature_set?

  $self->{static_config}{feature_types} = 
	{
	 'BioTIFFIN Motif'   => {
							 name        => 'BioTIFFIN Motif',
							 class       => 'Regulatory Motif',
							 description => 'BioTIFFIN motif',
							}
	};
  
  $self->{static_config}{analyses} = 
	{
	 'BioTIFFIN Motif' => { 
						   -logic_name    => 'BioTIFFIN Motif',
						   -description   => 'BioTIFFIN regulatory motif database',
						   -display_label => 'BioTIFFIN motifs',
						   -displayable   => 1,
						  },
	};
  
  $self->{static_config}{feature_sets} = 
	{
	 'BioTIFFIN Motif' => 
	 {
	  feature_set => {
					  -feature_type => 'BioTIFFIN Motif',
					  -analysis     => 'BioTIFFIN Motif',
					 },
	  xrefs => 0,
	 }
	};
 

  #Move xref flag here?
  $self->{config} =  {
					  'BioTIFFIN Motif' => {
										file  => $ENV{'EFG_DATA'}.'/input/BioTIFFIN/fly-tiffinScan-tiffin12.dm3.gff',
										gff_attrs => {
													  'ID' => 1,
													 },
									   },					  
					 };
  
  $self->validate_and_store_config([keys %{$self->{static_config}{feature_sets}}]);
  $self->set_feature_sets;

  return $self;
}



# Parse file and return hashref containing:
#
# - arrayref of features
# - arrayref of factors




sub parse_and_load {
  my ($self, $files, $old_assembly, $new_assembly) = @_;

  if(scalar(@$files) != 1){
	throw('You must provide a unique file path to load VISTA features from:\t'.join(' ', @$files));
  }


  my %slice_cache;
  my $extf_adaptor  = $self->db->get_ExternalFeatureAdaptor;
  my $dbentry_adaptor = $self->db->get_DBEntryAdaptor;
  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'BioTIFFINProjection');#do we need this?
  my $species = $self->db->species;
  if(! $species){
	throw('Must define a species to define the external_db');
  }
  #Just to make sure we hav homo_sapiens and not Homo Sapiens
  ($species = lc($species)) =~ s/ /_/;

  
  if(scalar @{$self->import_sets} != 1){
	throw('biotiffin parser currently only supports one import FeatureSet');
  }

  my ($import_set) = @{$self->import_sets};


  #foreach my $import_set(@{$self->import_sets}){
  $self->log_header("Parsing $import_set data");
  
  my %motif_cache; # name -> factor_id
  my $config = $self->{'config'}{$import_set};
  my $fset =  $self->{static_config}{feature_sets}{$import_set}{feature_set};
  my %gff_attrs =  %{$config->{'gff_attrs'}};
  
  
  # Parse motifs.txt file
  #my $file =  $config->{'file'};
  my $file = $files->[0];
  my $skipped = 0;
  my $motif_cnt = 0;
  my $factor_xref_cnt = 0;
  my $feature_cnt = 0;
  my $feature_target_cnt = 0;
	
  open (FILE, "<$file") || die("Can't open $file\n$!\n");

 LINE: while (my $line = <FILE>) {
	  chomp $line;
	  
	  #GFF3
		#3R      MotifScanner    TIFDMEM0000001  936391  936401  0.0     +       0
		#3R      MotifScanner    TIFDMEM0000001  13455911        13455921        0.0     -       0
    #3R      MotifScanner    TIFDMEM0000001  17062830        17062840        0.0     +       0
    #3R      MotifScanner    TIFDMEM0000001  17973965        17973975        0.0     +       0

	  #seq_name, source, feature, start, end, score, strand, frame, [attrs]
	  my ($chromosome, $program, $feature, $start, $end, $score, $strand, undef) = split /\t/o, $line;

	  if(! exists $slice_cache{$chromosome}){
				
				if($old_assembly){
						$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', 
																																							$chromosome, 
																																							undef, 
																																							undef, 
																																							undef, 
																																							$old_assembly);
				} else {
						$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', $chromosome);
				}
	  }
		
	  if(! defined  $slice_cache{$chromosome}){
				warn "Can't get slice $chromosome for motif $feature;\n";
				$skipped++;
				next;
	  }
				
		if(! exists $motif_cache{$feature}){
				
				$motif_cache{$feature} = $ftype_adaptor->fetch_by_name($feature);
				
				if(! defined $motif_cache{$feature}){
						
						($motif_cache{$feature}) = @{$ftype_adaptor->store(Bio::EnsEMBL::Funcgen::FeatureType->new
																															 (
																																-name  => $feature,
																																-class => $fset->feature_type->class,
																																-description => $fset->feature_type->description,
																															 ))};						
						
						$motif_cnt ++;
				}
	  }
		
		my $feature_type = $motif_cache{$feature};
		
	  #Now build actual feature

	  $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
				(
				 -display_label => $feature,
				 -start         => $start,
				 -end           => $end,
				 -strand        => (($strand eq '+') ? 1 : -1),
				 -feature_type  => $feature_type,
				 -feature_set   => $fset,
				 -slice         => $slice_cache{$chromosome},
				);
		
		
	  # project if necessary
	  if ($new_assembly) {
				$feature = $self->project_feature($feature, $new_assembly);
		
				if(! defined $feature){
						$skipped ++;
						next;
				}
	  }
	  
	  ($feature) = @{$extf_adaptor->store($feature)};
	  $feature_cnt++;


	}
	
	close FILE;
	
	$self->log("Loaded ".$fset->name);
	$self->log("$motif_cnt feature types");
	$self->log("$feature_cnt features");
	$self->log("Skipped $skipped features");

#}

  return;
}

1;
