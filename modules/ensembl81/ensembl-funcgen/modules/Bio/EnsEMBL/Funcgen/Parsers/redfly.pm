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

=cut

package Bio::EnsEMBL::Funcgen::Parsers::redfly;

use strict;
use warnings;
use File::Basename;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Funcgen::ExternalFeature;

use base qw( Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser );


# To get files for REDfly, download the following 2 GFF3 files (e.g. via wget):
#
# http://redfly.ccr.buffalo.edu/datadumps/tbfs_dump.gff
# http://redfly.ccr.buffalo.edu/datadumps/crm_dump.gff

#contact?

#TFBS
#2L 	REDfly	regulatory_region	2456365	2456372	.	.	.	ID="Unspecified_dpp:REDFLY:TF000068"; Dbxref="Flybase:FBgn0000490", "PMID:8543160", "REDfly:644, "FlyBase:"; Evidence="footprint/binding assay"; Factor="Unspecified"; Target="dpp";
#2L 	REDfly	regulatory_region	2456352	2456369	.	.	.	ID="dl_dpp:REDFLY:TF000069"; Dbxref="Flybase:FBgn0000490", "PMID:8458580", "REDfly:645, "FlyBase:FBgn0000463"; Evidence="footprint/binding assay"; Factor="dl"; Target="dpp";

#CRMs
#2L 	REDfly	regulatory_region	2455781	2457764	.	.	.	ID="dpp_intron2"; Dbxref="Flybase:FBgn0000490", "PMID:8167377", "REDfly:247; Evidence="reporter construct (in vivo)"; Ontology_term="FBbt:00005304";
#2L 	REDfly	regulatory_region	2445769	2446581	.	.	.	ID="dpp_dpp813"; Dbxref="Flybase:FBgn0000490", "PMID:7821226", "REDfly:246; Evidence="reporter construct (in vivo)"; Ontology_term="FBbt:00005653","FBbt:00001051";




sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_, type => 'REDfly');

  #Set default feature_type and feature_set config

  #We need to capture version/release/data of external feature sets.
  #This can be nested in the description?  Need to add description to feature_set?

  $self->{static_config}{feature_types} = 
	{
	 'REDfly TFBS'   => {
						 -name        => 'REDfly TFBS',
						 -class       => 'Transcription Factor',
						 -description => 'REDfly transciption factor binding site',
						},
	 'REDfly CRM' => {
					  -name        => 'REDfly CRM',
					  -class       => 'Regulatory Motif',
					  -description => 'REDfly cis regulatory motif',
					 },
	};
  
  
  $self->{static_config}{analyses} = 
	{
	 'REDfly TFBS' => { 
					   -logic_name    => 'REDfly TFBS',
					   -description   => 'REDfly transcription factor binding sites (http://redfly.ccr.buffalo.edu/)',
					   -display_label => 'REDfly TFBS',
					   -displayable   => 1,
					  },
	 
	  'REDfly CRM' => { 
					   -logic_name    => 'REDfly CRM',
					   -description   => 'REDfly cis regulatory motif (http://redfly.ccr.buffalo.edu/)',
					   -display_label => 'REDfly CRM',
					   -displayable   => 1,
					  },
	};
  
  $self->{static_config}{feature_sets} = 
	{
	 'REDfly TFBSs' => {
						feature_set =>
						{
						 -feature_type  => 'REDfly TFBS',
						 -analysis      => 'REDfly TFBS',
						},
						xrefs => 1,
					   },
	 
	 'REDfly CRMs' => {
					   feature_set   =>
					   {
						-feature_type  => 'REDfly CRM',
						-analysis      => 'REDfly CRM',
					   },
					   
					   xrefs => 1,
					  },
	};
 

  #Move xref flag here?
  $self->{config} =  {
					  'REDfly CRMs' => {
										#file  => $ENV{'EFG_DATA'}.'/input/REDFLY/redfly_crm.gff',
										gff_attrs => {
													  'ID' => 1,
													 },
									   },
					  
					  'REDfly TFBSs' => {
										 #file  => $ENV{'EFG_DATA'}.'/input/REDFLY/redfly_tfbs.gff',
										 gff_attrs => {
													   'ID' => 1,
													   'Factor' => 1,
													   'Target' => 1,
													  },
										 desc_suffix => ' binding site',
										}
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

  if(scalar(@$files) != 2){
	throw('You must currently define a crm and tfbs file to load redfly features from:\t'.join(' ', @$files));
  }

  #More validation of files here?
  $self->{config}{'REDfly CRMs'}{file}  = grep(/crm/,  @$files);
  $self->{config}{'REDfly TFBSs'}{file} = grep(/tfbs/, @$files);

  my %slice_cache;
  my $extf_adaptor  = $self->db->get_ExternalFeatureAdaptor;
  my $dbentry_adaptor = $self->db->get_DBEntryAdaptor;
  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'REDflyProjection');#do we need this?
  my $species = $self->db->species;

  if(! $species){
	throw('Must define a species to define the external_db');
  }

  #Just to make sure we hav homo_sapiens and not Homo Sapiens
  ($species = lc($species)) =~ s/ /_/;


  foreach my $import_set(@{$self->import_sets}){
	$self->log_header("Parsing $import_set data");

	my %factor_cache; # name -> factor_id
	my %target_cache;
	my $config = $self->{'config'}{$import_set};
	my $fset =  $self->{static_config}{feature_sets}{$import_set}{feature_set};
	my %gff_attrs =  %{$config->{'gff_attrs'}};
	
	
	# Parse motifs.txt file
	my $file =  $config->{'file'};
	my $skipped = 0;
	my $factor_cnt = 0;
	my $factor_xref_cnt = 0;
	my $feature_cnt = 0;
	my $feature_target_cnt = 0;
	
	open (FILE, "<$file") || die("Can't open $file\n$!\n");
	<FILE>; # skip header

	LINE: while (my $line = <FILE>) {
	  next if ($line =~ /^\s*\#/o || $line =~ /^\s*$/o);
	  chomp $line;
	  my %attr_cache;#Can we move this outside the loop and rely on it being reset each time?


	  #GFF3
	  #Is this format valid, missing " after REDfly xref
	  #2L 	REDfly	regulatory_region	2456365	2456372	.	.	.	ID="Unspecified_dpp:REDFLY:TF000068"; Dbxref="Flybase:FBgn0000490", "PMID:8543160", "REDfly:644, "FlyBase:"; Evidence="footprint/binding assay"; Factor="Unspecified"; Target="dpp";
	  #seq_name, source, feature, start, end, score, strand, frame, [attrs]
	  my ($chromosome, undef, $feature, $start, $end, undef, undef, undef, $attrs) = split /\t/o, $line;
	  my @attrs = split/\;\s+/o, $attrs;


	  #UCSC coords
	  $start ++;
	  $end ++;



	  foreach my $gff_attr(keys %gff_attrs){

		if(($attr_cache{$gff_attr}) = grep {/^${gff_attr}\=/} @attrs){
		  $attr_cache{$gff_attr} =~ s/(^${gff_attr}\=\")(.*)(\")/$2/;
		  
		  #warn "attr cache is  $attr_cache{$gff_attr} ";
	
		}
		else{
		  warn "Skipping import, unable to find mandatory $gff_attr attribute in:\t$line";
		  next LINE;
		}
	  }
	
   
	  #For TFBS
	  #Factor = coding gene name display_label
	  #Target = Target gene?
	  #Ignore other xrefs for name, just put ID in feature as display_label

	  #These are mixed up! and where not getting any coding xrefs!


	  #For CRM
	  #Can we split the ID and have Reguatory XREF?
	  #e.g.  ID="dpp_dpp813"; => dpp




	  #This can be moved to the BaseExternalParser

	  if(! exists $slice_cache{$chromosome}){
	
		if($old_assembly){
		  $slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', 
																			$chromosome, 
																			undef, 
																			undef, 
																			undef, 
																			$old_assembly);
		}else{
		  $slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', $chromosome);
		}
	  }

	  if(! defined  $slice_cache{$chromosome}){
		warn "Can't get slice $chromosome for motif $attr_cache{'ID'};\n";
		$skipped++;
		next;
	  }
	

	#get feature_type first

	#we are not maintaining this link in the DB!
	#Do we need another xref for this or a different table?

	  my $feature_type;
	  
	  #TFBSs
	  if(exists $attr_cache{'Factor'}){

		if(! exists $factor_cache{$attr_cache{'Factor'}}){
		  
		  $factor_cache{$attr_cache{'Factor'}} = $ftype_adaptor->fetch_by_name($attr_cache{'Factor'});
		  
		  if(! defined $factor_cache{$attr_cache{'Factor'}}){
		  
			#Would need to add CODING DBEntry here!
			#Will this work on a scalar ref to a hash?
			my $desc = (exists $config->{'desc_suffix'}) ? $attr_cache{'Factor'}.$config->{'desc_suffix'} : undef;
			
			($factor_cache{$attr_cache{'Factor'}}) = @{$ftype_adaptor->store(Bio::EnsEMBL::Funcgen::FeatureType->new
																			 (
																			  -name  => $attr_cache{'Factor'},
																			  -class => $fset->feature_type->class,
																			  -description => $desc,
																			 ))};
			
			$feature_type = $factor_cache{$attr_cache{'Factor'}};
			$factor_cnt ++;
			my $stable_id = $self->get_core_stable_id_by_display_name($self->db->dnadb, $attr_cache{'Factor'});
			
			#Handle release/version in xref version as stable_id version?
			
			if(! defined $stable_id){
			  warn "Could not generate CODING xref for feature_type:\t". $attr_cache{'Factor'}."\n";
			}else{
			  #warn "got $stable_id for ".$attr_cache{'Factor'};
			  my $dbentry = Bio::EnsEMBL::DBEntry->new(
													   -dbname                 => $species.'_core_Gene',
													   #-release                => $self->db->dnadb->dbc->dbname,
													   -status                 => 'KNOWNXREF',#This is for the external DB
													   #-display_label_linkable => 1,
													   -#db_display_name        => $self->db->dnadb->dbc->dbname,
													   -db_display_name        => 'EnsemblGene',
													   -type                   => 'MISC',#Is for the external_db 
													   -primary_id             => $stable_id,
													   -display_id             => $attr_cache{'Factor'},
													   -info_type              => 'MISC',
													   -into_text              => 'GENE',
													   -linkage_annotation     => 'REDfly Coding'
													   -analysis               => $fset->analysis,

													   #-description            => 'cisRED motif gene xref',#This is now generic and no longer resitricted to REDfly
													   #could have version here if we use the correct dnadb to build the cache
													  );
			  
			  $dbentry_adaptor->store($dbentry, $factor_cache{$attr_cache{'Factor'}}->dbID, 'FeatureType', 1);#1 is ignore release flag
			  $factor_xref_cnt ++;
			}
		  }
		}
	  }
	  else{
		#CRMs
		$feature_type = $fset->feature_type;
	  }


	  #Now build actual feature
	  $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
	  (
	   -display_label => $attr_cache{'ID'},
	   -start         => $start,
	   -end           => $end,
	   -strand        => 0,
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

	  my $target = (exists $attr_cache{'Target'}) ?  $attr_cache{'Target'} : (split/_/, $attr_cache{'ID'})[0];
	  my $stable_id;

	  if($target ne 'Unspecified'){
		$stable_id = $self->get_core_stable_id_by_display_name($self->db->dnadb, $target);
	  }


	  if(! defined $stable_id){
		warn "Could not generate TARGET xref for feature:\t". $attr_cache{'ID'}."\n" if $target ne 'Unspecified';
	  }
	  else{
		#Handle release/version in xref version as stable_id version?		
		my $dbentry = Bio::EnsEMBL::DBEntry->new(
												 -dbname                 => $species.'_core_Gene',
												 #-release                => $self->db->dnadb->dbc->dbname,
												 -status                 => 'KNOWNXREF',
												 #-display_label_linkable => 1,
												 -#db_display_name        => $self->db->dnadb->dbc->dbname,
												 -db_display_name        => 'EnsemblGene',
												 -type                   => 'MISC',#
												 -primary_id             => $stable_id,
												 -display_id             => $target,
												 -info_type              => 'MISC',
												 -info_text              => 'GENE',
												 -linkage_annotation     => $fset->feature_type->name.' Target',
												 -analysis               => $fset->analysis,

												 #could have version here if we use the correct dnadb to build the cache
											);

		$dbentry_adaptor->store($dbentry, $feature->dbID, 'ExternalFeature', 1);#1 is ignore release flag
		
		$feature_target_cnt ++;
	  }
	}
	
	close FILE;

	$self->log("Loaded ".$fset->name);
	$self->log("$factor_cnt feature types");
	$self->log("$factor_xref_cnt feature type coding xrefs");
	$self->log("$feature_cnt features");
	$self->log("$feature_target_cnt feature target xrefs");
	$self->log("Skipped $skipped features");

  }

  return;
}




1;
