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

Bio::EnsEMBL::Funcgen::MAGE
  
=head1 SYNOPSIS

my $imp = Bio::EnsEMBL::Funcgen::Importer->new(%params);
$imp->register_experiment();


=head1 DESCRIPTION

B<This program> is a base main class for all MAGE type array importers(e.g. Nimblegen).

=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Parsers::MAGE;

use strict;
use warnings;
use File::Path;
use Bio::MAGE::XMLUtils;
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( get_date open_file run_system_cmd );
use Bio::EnsEMBL::Funcgen::Utils::Helper;

use base qw( Bio::EnsEMBL::Funcgen::Utils::Helper );



################################################################################

=head2 new

 Description : Constructor method

 Arg  [1]    : hash containing optional attributes:

 ReturnType  : Bio::EnsEMBL::Funcgen::MAGE
 Example     : my $Exp = Bio::EnsEMBL::Nimblegen->new(%params);
 Exceptions  : throws if mandatory params are not set or DB connect fails
 Caller      : General
 Status      : Medium - potential for %params names to change, remove %attrdata?

=cut

################################################################################

sub new{
  my ($caller) = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  #This needs to be an Importer!
  throw("This is base class for the experiment Bio::EnsEMBL::Funcgen::Parsers, needs to inherit from Bio::EnsEMBL::Funcgen::Importer") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));


  #Are we not passing any Helper params?
  #Log file etc is set in the run script

  my ($write_mage, $no_mage, $vendor)
	= rearrange(['WRITE_MAGE', 'NO_MAGE', 'VENDOR'], @_);


  #$self->{'update_xml'} = $update_xml || 0;
  $self->{'write_mage'} = $write_mage || 0;
  $self->{'no_mage'} = $no_mage || 0;
  $self->{'vendor'} = $vendor;


  if ($self->vendor ne 'NIMBLEGEN'){
	$self->{'no_mage'} = 1;
	warn "Hardcoding no_mage for non-NIMBLEGEN imports";
  }

 
  if($self->{'no_mage'} && $self->{'write_mage'}){
	throw('-no_mage and -write_mage options are mutually exclusive, please select just one');
  }
    
  return $self;
}

=head2 process_experiment_config

  Example    : $self->init_experiment_import();
  Description: Initialises import by creating working directories 
               and by storing the Experiemnt
  Returntype : none
  Exceptions : warns and throws depending on recover and Experiment status 
  Caller     : general
  Status     : at risk - merge with register exeriment

=cut

#This is actually processing the tab2mage file & writing XML

sub process_experiment_config{
  my $self = shift;

   #Here, this is where we need to call the a Parser from Importer to do this for only MAGE experiments
  #validate_import?

  #This is only used for the first test below.
  my $exp_adaptor = $self->db->get_ExperimentAdaptor();  
  my $xml = $exp_adaptor->fetch_mage_xml_by_experiment_name($self->name);# if $self->{'write_xml'};

  #DO NOT CHANGE THIS LOGIC!
  #write mage if we specify or we don't have a the final xml or the template
  #recovery is turned on to stop exiting when previously stored chips are found from the 'write_mage' run.
  #This does mean that if you import without running the write_mage step
  #you could potentially be overwriting someone elses experiment info
  #No way of getting around this, need to make warning obvious, add to end of log!!!
  #We always want to write and update xml and ResultSets if we're running the 2nd stage of the import
  #Why would we ever want to skip the validate process?
  #Leave for now as this is working as we want it
  #But propose to remove skip functionality
  
  if( ! $self->{'no_mage'}){

	if($self->{'write_mage'} || !( -f $self->get_config('tab2mage_file') || $xml)){
	  $self->{'write_mage'} = 1;
	  $self->backup_file($self->get_config('tab2mage_file'));
	}
	#elsif($xml && (! $self->{'update_xml'})){#Changed this so we always update
	#elsif(! $self->{'update_xml'}){



	  #Here, we need to always update_xml
	  #If we are doing the 2nd stage
	  #Currently this is skipping as we haven't explicitly set it
	  #To remove this...
	  #what we need to do is check that we don't test for update_xml, 
	  # i.e. assuming that we're running the second stage of the import.
	  # Therefore we need a boolean to set whether it is the first stage..else update_xml implicit
	  # write mage is explicit flag
	  # Or if we have not tab2mage file?
	  # we can then override this explicitly with update_xml?
	  # WE're never likely edit the xml directly, so we always want to validate and update
	  # so write mage flag become update_experiment?  No this is no obvious behaviour
	  # We need to warn about removing the write_mage flag after we have updated it
	  # Otherwise we will never get to 2nd stage


	#No mage is still valid as we may want to jus import and experiment
	#Before receiving correct meta data
	#When we can then rerun the import with -write_mage to update the resultsets

	#  $self->{'recover'} = 1;
	#  $self->{'skip_validate'} = 1;
	#}
	elsif( -f $self->get_config('tab2mage_file')){#Run Tab2Mage

	  $self->backup_file($self->get_config('mage_xml_file'));
	  my $cmd = 'tab2mage.pl -e '.$self->get_config('tab2mage_file').' -k -t '.$self->get_dir('output').
		' -c -d '.$self->get_dir('results');
	  
	  $self->log('Reading tab2mage file');
	  my $t2m_exit_code = run_system_cmd($cmd, 1);#no exit flag due to non-zero exit codes	  
	  warn "tab2mage exit code is  $t2m_exit_code"; 
	  
	  if(! ($t2m_exit_code > -1) && ($t2m_exit_code <255)){
		throw("tab2mage failed.  Please check and correct:\t".$self->get_config('tab2mage_file')."\n...and try again");
	  }
	  
	  $self->{'recover'} = 1;
	}
  }

  return;
}

=heead init_tab2mage_export

  Example    : $self->init_tab2mage_export;
  Description: Writes the standard experiment section of the tab2mage file
  Returntype : FileHandle
  Exceptions : ???
  Caller     : general
  Status     : at risk

=cut

sub init_tab2mage_export{
  my $self = shift;

  $self->backup_file($self->get_config('tab2mage_file')) if(-f $self->get_config('tab2mage_file'));

  my $t2m_file = open_file($self->get_config('tab2mage_file'), '>');

  #reformat this
  my $exp_section = "experiment section\ndomain\t".(split/@/, $self->contact())[1]."\naccession\t\n".
	"quality_control\tbiological_replicate\nexperiment_design_type\tbinding_site_identification\n".
	  "name\t".$self->name()."\nrelease_date\t\nsubmission_date\t\nsubmitter\t???\n".
		"submitter_email\t???\ninvestigator\t???\ninvestigator_email\t???\norganization\t???\naddress\t".
		  "???\npublication_title\t\nauthors\t\njournal\t\nvolume\t\nissue\t\npages\t\nyear\t\npubmed_id\t\n";

  my $protocol_section = "Protocol section\naccession\tname\ttext\tparameters\n";

  foreach my $protocol(sort (keys %{$self->get_config('protocols')})){
	$protocol_section .= $self->get_config('protocols')->{$protocol}->{'accession'}.
	  "\t".$self->get_config('protocols')->{$protocol}->{'name'}.
		"\t".$self->get_config('protocols')->{$protocol}->{'text'}."\t";

	$protocol_section .= (defined $self->get_config('protocols')->{$protocol}->{'parameters'}) ?
	  $self->get_config('protocols')->{$protocol}->{'parameters'}."\t\n" : "\t\n";
  }

  #File[raw]	Array[accession]	Array[serial]	Protocol[grow]	Protocol[treatment]	Protocol[extraction]	Protocol[labeling]	Protocol[hybridization]	Protocol[scanning]	BioSource	Sample	Extract	LabeledExtract	Immunoprecipitate	Hybridization	BioSourceMaterial	SampleMaterial	ExtractMaterial	LabeledExtractMaterial	Dye	BioMaterialCharacteristics[Organism]	BioMaterialCharacteristics[BioSourceType]	BioMaterialCharacteristics[StrainOrLine]	BioMaterialCharacteristics[CellType]	BioMaterialCharacteristics[Sex]	FactorValue[StrainOrLine]	FactorValue[Immunoprecipitate]


  #Need to do this bit better?
  #have array of fields.  We can then populate a hash in the read method based on field names, then use the array to print in order

  my $hyb_header = "\nHybridization section\n".join("\t", @{$self->hybridisation_fields()});

  print $t2m_file $exp_section."\n".$protocol_section."\n".$hyb_header."\n";

  return $t2m_file;
}


#Move to MAGE package?

sub hybridisation_fields{
  my $self = shift;

  return ['File[raw]', 'Array[accession]', 'Array[serial]', 
		  (map 'Protocol['.$_.']', (sort (keys %{$self->get_config('protocols')}))),
		  'BioSource', 'Sample', 'Extract', 'LabeledExtract', 'Immunoprecipitate', 'Hybridization', 
		  'BioSourceMaterial', 'SampleMaterial', 'ExtractMaterial', 'LabeledExtractMaterial',
		  'Dye', 'BioMaterialCharacteristics[Organism]', 'BioMaterialCharacteristics[BioSourceType]',	
		  'BioMaterialCharacteristics[StrainOrLine]', 'BioMaterialCharacteristics[CellType]', 
		  'BioMaterialCharacteristics[Sex]', 'FactorValue[StrainOrLine]', 'FactorValue[Immunoprecipitate]'];
}



#=head2 register_experiment
#  
#  Example    : $imp->register_experiment()
#  Description: General control method, performs all data import and normalisations
#  Arg [1]    : optional - dnadb DBAdaptor
#  Returntype : none
#  Exceptions : throws if arg is not Bio::EnsEMBL::DBSQL::DBAdaptor
#  Caller     : general
#  Status     : Medium
#
#=cut

#write/validate_mage

sub write_validate_experiment_config{
  my $self = shift;
    

  if($self->{'write_mage'} || $self->{'no_mage'}){
	$self->read_data("array");

	if(! $self->{'no_mage'}){
	  $self->log("PLEASE CHECK AND EDIT AUTOGENERATED TAB2MAGE FILE:\t".$self->get_config('tab2mage_file'));
	  #we could make this print only if it was set by the user, not by the Importer
	  $self->log('REMEMBER TO REMOVE -write_mage FLAG BEFORE UPDATING');
	  exit;
	}
  }
  elsif(! $self->{'no_mage'}){#This should be a no_channel flag, set dependent on import mode(gff_chip, gff_chan)
	#Need to accomodate chip level imports in validate?
	
	if (! $self->{'skip_validate'}){

	  $self->log("Validating mage file:\t".$self->get_config('mage_xml_file'));


	  #Updating ResultSets:
	  #Given that we might want to add a chip to an experiment we will also need to update the tab2MAGE
	  #mage_xml and ResultSets accordingly.
	  
	  #This should happen if we specify update_xml
	  #Should recovery also always force update?
	  #Considering the two run modes, write tab2mage & validate and import
	  #There is a subtle difference between recovery and update mage_xml
	  #Do we always run in recovery mode for the validate&import step?
	  #Yes we do, so can't guarantee the this means we want to update.
	  #So we need to change update_xml to update to reflect the changed functionality on ResultSets
	  
	  #If we run an update without on then chips will be loaded but xml and ResultSets will not be altered :(
	  #If we're running the 2nd stage we should always be updating the xml anyway!!!!
	  #As there is no reason to rerun the validate & import step without it.(unless we're debugging maybe)
	  #So why should we ever run without it?
	  
	  #To update ResultSets we validate as normal and then update where appropriate
	  #What has precedence? Replicate name?
	  #Update echip types as appropriate
	  #What if this invalidates original rsets?
	  #Then list sets not covered for removal by script?
	  
	  
	  
	  my (%echips, @log);
	  my $rset_adaptor = $self->db->get_ResultSetAdaptor;
	  my $chan_anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name('RawValue');
	  
	  #need to change this to default analysis
	  #There we issues with setting VSN_GLOG as an env var
	  #as this is tested for and the norm was skipped for some reason?
	  my $chip_anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($self->norm_method());
	  
	  #Try import sets like this first, so we know ther is new data
	  my $chan_rset = $self->get_import_ResultSet($chan_anal, 'channel');
	  my $rset =  $self->get_import_ResultSet($chip_anal, 'experimental_chip');
	  
	  
	  #Else get them anyway and log
	  if(! $rset){
		
		if($chan_rset){
		  $self->log('Identified partial Channel only import, updating MAGE-XML');
		}
		else{
		  ($chan_rset) = @{$rset_adaptor->fetch_all_by_name_Analysis($self->experiment->name.'_IMPORT', $chan_anal)};
		  #Don't need to test for >1 here as this has already been done in get_import_ResultSet
		  $self->log('All ExperimentalChips imported, updating MAGE-XML only');
		}
		
		($rset) = @{$rset_adaptor->fetch_all_by_name_Analysis($self->experiment->name.'_IMPORT', $chip_anal)};
	  }


  #This will never happen now due to the change tab2mage rules in init_experiment
  #Remove?
  if(! $rset){
	throw('Cannot find ResultSet, are you trying to import a new experiment which already has a tab2mage file present?  Try removing the file, or specifying the -write_mage flag to parse_and_import.pl');
  }

  if(! -l $self->get_dir('output').'/MAGE-ML.dtd'){
	system('ln -s '.$ENV{'EFG_DATA'}.'/MAGE-ML.dtd '.$self->get_dir('output').'/MAGE-ML.dtd') == 0 || 
	  throw('Failed to link MAGE-ML.dtd');
  }
 

  $self->log('VALIDATING MAGE XML');
  my $reader = Bio::MAGE::XML::Reader->new();
	  my $mage_xml ||= $self->get_config('mage_xml_file');
  $self->{'mage'} = $reader->read($mage_xml);
  
  #this should only ever return 1 for an import
  foreach my $mage_exp(@{$self->{'mage'}->getExperiment_package->getExperiment_list()}){	

	if($mage_exp->getName() ne $self->name()){
	  $self->log('MAGE experiment name ('.$mage_exp->getName().') does not match import name ('.$self->name().')');
	}
	
	#add more experiment level validation here?

	foreach my $assay (@{$mage_exp->getBioAssays()}){

	  if($assay->isa('Bio::MAGE::BioAssay::PhysicalBioAssay')){#channel
		$self->log('Validating PhysicalBioAssay "'.$assay->getName()."'\n");#hyb name(this is the file name for measured assays

		my $bioassc = $assay->getBioAssayCreation();#This is a Hybridisation
		my $array = $bioassc->getArray();#this is an ArrayChip
		my $design_id = $array->getArrayDesign->getIdentifier();
		my $chip_uid = $array->getArrayIdentifier();
	
	
		foreach my $echip(@{$rset->get_ExperimentalChips()}){
		

		  if($echip->unique_id() eq $chip_uid){
			$self->log("Found ExperimentalChip:\t".$chip_uid);

			if(! exists $echips{$chip_uid}){
			  $echips{$chip_uid} = {(
									 total_biorep            => undef,
									 total_biotechrep        => undef,
									 experimental_biorep     => undef,
									 experimental_biotechrep => undef,
									 total_dye               => undef,
									 experimental_dye        => undef,
									 cell_type               => undef,
									 feature_type            => undef,
									)};
			}

			#Validate ArrayChip
			my ($achip) = @{$self->db->get_ArrayChipAdaptor->fetch_all_by_ExperimentalChips([$echip])};

			if($achip->design_id() ne $design_id){
			  push @log, "ArrayDesign Identifier (${design_id}) does not match ArrayChip design ID (".
				$achip->design_id().")\n\tSkipping channel and replicate validation";
			  #skip the channel/replicate validation here?	  
			} 
			else {			#validate channels and replicate names
					
			  foreach my $src_biomat (@{$bioassc->getSourceBioMaterialMeasurements()}) { #Channel materials(X1)?
				my $biomat = $src_biomat->getBioMaterial();	#LabelledExtract (IP/Control)
				#we could sub this passing $echip and biomat?
				#messy to pass regexs and populate correct echip hash attrs
				#also messy to populate log
				#keeping nested loop also prevents further obfuscation
				#do we need to do all the defined checks, or maybe just the first one?
				#Then we can skip all following warning?

				foreach my $treat (@{$biomat->getTreatments()}) {
				  #As there is effectively one more level of material extraction for the IP channel
				  #this loop will returns materials an iteration out of sync for each channel

				  foreach my $ssrc_biomat (@{$treat->getSourceBioMaterialMeasurements()}) {	#Channel measurement(x1)
					my $sbiomat = $ssrc_biomat->getBioMaterial();
					#This will either be techrep name for control of IP name for experimental channel
					#SOM0035_BR1_TR2 IP  #Immunoprecicpitate
					#SOM0035_BR1_TR2     #Extract
				
					if ($sbiomat->getName() =~ /BR[0-9]+_TR[0-9]+$/o) { #Total

					  if (! defined $echips{$chip_uid}{'total_biotechrep'}) {
						$echips{$chip_uid}{'total_biotechrep'} = $sbiomat->getName();
					  }
					  else{
						push @log, "Found two TOTAL Channels on same chip with biotechreps:\t".$sbiomat->getName().
						  " and ".$echips{$chip_uid}{'total_biotechrep'};
					  }
					}else{#Experimental

					  #get feature type from assay
					  my $fv_ref = $assay->getBioAssayFactorValues();
					  if(! defined $fv_ref){
						throw('No FactorValues found, you must populate the "Immunoprecipitate" field. Maybe you forgot to specify -feature_type?'); 
					  }

					  my ($feature_type);
					  
					  foreach my $fvalue(@{$fv_ref}){
						
						if($fvalue->getValue()->getCategory() eq 'Immunoprecipitate'){
						  $feature_type = $fvalue->getName();
						  $feature_type =~ s/anti\s*-\s*//;
						  $feature_type =~ s/\s*antibody\s*//;
						}
					  }
					  $echips{$chip_uid}{'feature_type'} = $feature_type;
					}

					foreach my $ttreat (@{$sbiomat->getTreatments()}) {
									  
					  foreach my $tsrc_biomat (@{$ttreat->getSourceBioMaterialMeasurements()}) {
						my $tbiomat = $tsrc_biomat->getBioMaterial();
						#SOM0035_BR1_TR2     #Extract (exp)
						#SOM0035_BR1              #Sample (total)
									   	
						if ($tbiomat->getName() =~ /BR[0-9]+_TR[0-9]+$/o) { #experimental
						  
						  if (! defined $echips{$chip_uid}{'experimental_biotechrep'}) {
							$echips{$chip_uid}{'experimental_biotechrep'} = $tbiomat->getName();
						  }
						  else{
							push @log, "Found two EXPERIMENTAL Channels on same chip with biotechreps:\t".$tbiomat->getName().
							  " and ".$echips{$chip_uid}{'experimental_biotechrep'};
						  }
						
						  my $dye = $biomat->getLabels()->[0]->getName();
							
						  foreach my $chan (@{$echip->get_Channels()}) {
							  
							if ($chan->type() eq 'EXPERIMENTAL') {
								
							  if (uc($dye) ne uc($chan->dye())) {
								push @log, "EXPERIMENTAL channel dye mismatch:\tMAGE = ".uc($dye).' vs DB '.uc($chan->dye);
							  } else {
								$echips{$chip_uid}{'experimental_dye'} = uc($dye);
							  }
							}
						  }
						} 
						else { #control
											  
						  if (! defined $echips{$chip_uid}{'total_biorep'}) {
							$echips{$chip_uid}{'total_biorep'} = $tbiomat->getName();
						  }
						  else{
							push @log, "Found two TOTAL Channels on same chip with biotechreps:\t".$tbiomat->getName().
							  " and ".$echips{$chip_uid}{'total_biorep'};
						  }
						  
						  my $dye = $biomat->getLabels()->[0]->getName();
						
						  foreach my $chan (@{$echip->get_Channels()}) {
						  
							if ($chan->type() eq 'TOTAL') {
							
							  if (uc($dye) ne uc($chan->dye())) {
								push @log, "TOTAL channel dye mismatch:\tMAGE = ".uc($dye).' vs DB '.uc($chan->dye);
							  } 
							  else {
								$echips{$chip_uid}{'total_dye'} = uc($dye);
							  }
							}
						  }
						}
						#could do one more iteration and get Source and FeatureType?
						#we should really extend this, and then update the EC cell_type and feature_types
						#these features might not be biotmats tho...need to check


						foreach my $ftreat (@{$tbiomat->getTreatments()}) {
									  
						  foreach my $fsrc_biomat (@{$ftreat->getSourceBioMaterialMeasurements()}) {
							my $fbiomat = $fsrc_biomat->getBioMaterial();
							#EXPERIMENTAL - biorep 
							#TOTAL        - source/cell type
							my $cell_type;

							if($fbiomat->getName() =~ /BR[0-9]+$/o){#EXPERIMETNAL
							
							  if(! defined $echips{$chip_uid}{'experimental_biorep'}){
								$echips{$chip_uid}{'experimental_biorep'} = $fbiomat->getName();
							  }
							  else{
								push @log, "Found two Experimental Channels on same chip with bioreps:\t".$fbiomat->getName().
								  " and ".$echips{$chip_uid}{'experimental_biorep'};
							  }


							  #last treatment/measurement/biomat level should go here
							  #as TOTAL channel does not have another level and will fail
							  foreach my $xtreat (@{$fbiomat->getTreatments()}) {
								
								foreach my $xsrc_biomat (@{$xtreat->getSourceBioMaterialMeasurements()}) {
								  my $xbiomat = $xsrc_biomat->getBioMaterial();
								  
								  foreach my $char(@{$xbiomat->getCharacteristics()}){
									$cell_type = $char->getValue() if($char->getCategory() eq 'CellType');
								  }
								}
							  }

							}else{#this should be BioSource
							  #which should have CellType as characteristic
							  #we could change tab2mage and have this as a factor value, 
							  #but don't want to start messing with "standard" format
						
							  foreach my $char(@{$fbiomat->getCharacteristics()}){
								$cell_type = $char->getValue() if($char->getCategory() eq 'CellType');
							  }
							}
						
							#can have cell_type validation here
							if(! defined $echips{$chip_uid}{'cell_type'}){
							  $echips{$chip_uid}{'cell_type'} = $cell_type;
							}
							elsif( $echips{$chip_uid}{'cell_type'} ne $cell_type){
							  push @log, "Found Channels on same chip (${chip_uid}) with different cell types:\t".
								$cell_type." and ".$echips{$chip_uid}{'cell_type'};
							}
						  }
						}
					  }
					}
				  }
				}
			  }
			}
		  }						#end of echip
		}						#end of foreach echip
	  }							#end of physbioassay	
	}							#end of foreach assay
  }								#end of foreach exp



  #we should fail here with log before we update the result sets

   #we need to build rep names
  #we're currently using sample labels, in the tab2mage file
  #altho' previous sets have been using exp name
  #these have been manually patched afterwards

  #More desirable to have exp name as rset name, but no way of doing BR validation
  #based on sample label, if we don't have it in the tab2mage
  #if we change it in the DB then we need to update the tab2mage

  #no way to do this when generating tab2mage as the user hasn't yet defined the reps
  #we could just make reps based on sample labels
  #then we just assume that alterations made by the user are correct
  #as we can no longer validate using sample labels
  #can still validate using cell/feature type

  #no longer need vendor specific validation as this will be done in tab2mage generation


  #We need to validate reps here
  #then update ec records as appropriate and then create rsets

  my (%bio_reps, %tech_reps);
  my $ct_adaptor = $self->db->get_CellTypeAdaptor();
  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
 
#select rs.*, ec.*, c.* from result_set rs, chip_channel cc, channel c, experimental_chip ec where rs.result_set_id=cc.result_set_id and cc.table_name='experimental_chip' and cc.table_id=ec.experimental_chip_id and cc.table_id=c.experimental_chip_id order by name;

  foreach my $echip (@{$rset->get_ExperimentalChips()}) {

	my ($biorep, $biotechrep);

	if (! exists $echips{$echip->unique_id()}) {
	  push @log, "No MAGE entry found for ExperimentalChip:\t".$echip->unique_id();
	} 
	else {

	  foreach my $chan_type('total', 'experimental'){
		
		$biorep = $echips{$echip->unique_id()}{$chan_type.'_biorep'};
		$biotechrep = $echips{$echip->unique_id()}{$chan_type.'_biotechrep'};

		if (! defined $biotechrep) {
		  push @log, 'ExperimentalChip('.$echip->unique_id().') Extract field do not meet naming convention(SAMPLE_BRN_TRN)';
		}							#! defined biorep? will never occur at present
		elsif ($biotechrep !~ /\Q$biorep\E/) {
		  push @log, "Found Extract(techrep) vs Sample(biorep) naming mismatch\t${biotechrep}\tvs\t$biorep";
		} 
		
		if ( ! $echips{$echip->unique_id()}{$chan_type.'_dye'}) {
		  push @log, "No ".uc($chan_type)." channel found for ExperimentalChip:\t".$echip->unique_id();
		}

	  }

	  #Is this is really implicit in the test above
	  if($echips{$echip->unique_id()}{'experimental_biorep'} ne $echips{$echip->unique_id()}{'total_biorep'}){
		push @log, "Found biorep mismatch between channels of ExperimentalChip ".$echip->unique_id().":\n".
		  "\tEXPERIMENTAL\t".$echips{$echip->unique_id()}{'experimental_biorep'}."\tTOTAL\t".
			$echips{$echip->unique_id()}{'total_biorep'};
	  }

	  #Is this is really implicit in the test above
	  if($echips{$echip->unique_id()}{'experimental_biotechrep'} ne $echips{$echip->unique_id()}{'total_biotechrep'}){
		push @log, "Found biotechrep mismatch between channels of ExperimentalChip ".$echip->unique_id().":\n".
		  "\tEXPERIMENTAL\t".$echips{$echip->unique_id()}{'experimental_biotechrep'}."\tTOTAL\t".
			$echips{$echip->unique_id()}{'total_biotechrep'};
	  }

		   
	}


	#Now we need to validate ec has same feature/cell type as other ecs in this br
	#this does not handle import sets which ARE allowed to have same name but different types

	#warn "Processing ".$echip->unique_id()." $biorep $biotechrep";

	
	if(exists $bio_reps{$biorep}){


	  if(! defined $bio_reps{$biorep}{'cell_type'}){
		push @log, "Found undefined CellType for biorep $biorep";
	  }
	  elsif($bio_reps{$biorep}{'cell_type'}->name() ne  $echips{$echip->unique_id()}{'cell_type'}){
		push @log, "Found CellType mismatch between $biorep and ExperimentalChip ".$echip->unique_id();
	  }
	  
	  
	  if(! defined $bio_reps{$biorep}{'feature_type'}){
		push @log, "Found undefined FeatureType for biorep $biorep";
	  }
	  elsif($bio_reps{$biorep}{'feature_type'}->name() ne  $echips{$echip->unique_id()}{'feature_type'}){
		push @log, "Found FeatureType mismatch between $biorep and ExperimentalChip ".$echip->unique_id();
	  }
	  
	  #warn "$biorep exists with\t".$bio_reps{$biorep}{'cell_type'}->name().' '.$bio_reps{$biorep}{'feature_type'}->name();
	  
	  #We need to set the tech rep here too!
	  #Do we need to validate this also, as above.
	  #This would be overkill due to the inherant nature of the TR to BR relationship

	  if(! exists $tech_reps{$biotechrep}){
		 $tech_reps{$biotechrep}{'cell_type'} = $bio_reps{$biorep}{'cell_type'};
		 $tech_reps{$biotechrep}{'feature_type'} = $bio_reps{$biorep}{'feature_type'};
	  }


	}else{

	  #warn "Creating new BR $biorep and TR $biotechrep";

	  if(defined $echips{$echip->unique_id()}{'cell_type'}){

		my $cell_type = $ct_adaptor->fetch_by_name($echips{$echip->unique_id()}{'cell_type'});

		if(! defined $cell_type){
		  push @log, 'CellType '.$echips{$echip->unique_id()}{'cell_type'}.' does not exist in the database, please use the import_type.pl script';
		}else{
		  $bio_reps{$biorep}{'cell_type'} = $cell_type;
		  $tech_reps{$biotechrep}{'cell_type'} = $cell_type;
		#  warn "Setting ".$echip->unique_id()." $biorep $biotechrep ".$cell_type->name;
		}
	  }else{
		warn "No CellType specified for ExperimentalChip:\t".$echip->unique_id()."\n";
	  }


	  if(defined $echips{$echip->unique_id()}{'feature_type'}){
		my $feature_type = $ft_adaptor->fetch_by_name($echips{$echip->unique_id()}{'feature_type'});

		if(! defined $feature_type){
		  push @log, 'FeatureType '.$echips{$echip->unique_id()}{'feature_type'}.' does not exist in the database, please use the import_type.pl script';
		}
		else{
		  $bio_reps{$biorep}{'feature_type'} = $feature_type;
		  $tech_reps{$biotechrep}{'feature_type'} = $feature_type;

		  #warn "Setting ".$echip->unique_id()." $biorep $biotechrep ".$feature_type->name;
		}
	  }else{
		warn "No FeatureType specified for ExperimentalChip:\t".$echip->unique_id()."\n";
	  }
	}

	push @{$tech_reps{$biotechrep}{'echips'}}, $echip->unique_id();
	push @{$bio_reps{$biorep}{'echips'}}, $echip->unique_id();	
  }
  



  if (@log) {
	$self->log("MAGE VALIDATION REPORT\n::\t".join("\n::\t", @log));
	throw("MAGE VALIDATION FAILED\nPlease correct tab2mage file and try again:\t".$self->get_config('tab2mage_file'));
  } else {
	$self->log('MAGE VALDIATION SUCCEEDED');
  }


  #we also need to build the tech rep results sets(not displayable)
  #do we need to have result sets for each biorep too?
  #update ExperimentalChip replicate info
  my (%rsets);
  my %types = (
			   feature => {},
			   cell    => {},
			  );


 
  #This needs to update and split the import/top level sets so they are of same types
  #update ec type here as we have ec context
  #careful not to update multiple times, just once for each ec
 
  my $eca = $self->db->get_ExperimentalChipAdaptor();
  

  foreach my $echip (@{$rset->get_ExperimentalChips()}) {
	my ($cell_type, $feature_type);

	#Set biorep info and rset
	foreach my $biorep (keys %bio_reps){

	  foreach my $chip_uid(@{$bio_reps{$biorep}{'echips'}}){

		if($chip_uid eq $echip->unique_id()){
		  $echip->biological_replicate($biorep);
		  $cell_type = $bio_reps{$biorep}{'cell_type'};
		  $feature_type = $bio_reps{$biorep}{'feature_type'};

		  if(! defined $rsets{$biorep}){
			
			$rsets{$biorep} = Bio::EnsEMBL::Funcgen::ResultSet->new
			  (
			   -NAME         => $biorep,#this may not be unique, prepend with exp name? Force method to use Experiment_and_name?
			   -ANALYSIS     => $rset->analysis(),
			   -TABLE_NAME   => 'experimental_chip',
			   -FEATURE_TYPE => $feature_type,
			   -CELL_TYPE    => $cell_type,
			  );

			#record cell and feature types
			$types{'feature'}{$feature_type->name()} = $feature_type;
			$types{'cell'}{$cell_type->name()} = $cell_type;
			$self->log("Created BioRep ResultSet:\t".$rsets{$biorep}->log_label);
		  }
		  
		  $rsets{$biorep}->add_table_id($echip->dbID(), $rset->get_chip_channel_id($echip->dbID()));
		}
	  }
	}

	#reset echip types
	$echip->feature_type($feature_type);
	$echip->cell_type($cell_type);

	
	#set tech rep info and rset
	foreach my $techrep(keys %tech_reps){
	  
	  foreach my $chip_uid(@{$tech_reps{$techrep}{'echips'}}){
		
		if($chip_uid eq $echip->unique_id()){
		  $echip->technical_replicate($techrep);

		  if(! defined $rsets{$techrep}){
			$rsets{$techrep} = Bio::EnsEMBL::Funcgen::ResultSet->new
			  (
			   -NAME       => $techrep,#this may not be unique, prepend with exp name? Force method to use Experiment_and_name?
			   -ANALYSIS   => $rset->analysis(),
			   -TABLE_NAME => 'experimental_chip',
			   -FEATURE_TYPE => $tech_reps{$techrep}{'feature_type'},
			   -CELL_TYPE    => $tech_reps{$techrep}{'cell_type'},
			  );

			$self->log("Created TechRep ResultSet:\t".$rsets{$techrep}->log_label);
		  }
		  $rsets{$techrep}->add_table_id($echip->dbID(), $rset->get_chip_channel_id($echip->dbID()));
		}
	  }
	}

	$echip->adaptor->update_replicate_types($echip);#store rep info
  }


  ### Reset/Update/Clean import sets type fields
  my $sql;

  if(scalar keys %{$types{'feature'}} >1){
	$self->log('Resetting IMPORT FeatureType to NULL for multi-FeatureType Experiment');
	$sql = "UPDATE result_set set feature_type_id='NULL' where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';

  }else{
	my ($ftype) = values %{$types{'feature'}};

	if(! defined $rset->feature_type()){
	  $self->log('Updating IMPORT FeatureType to '.$ftype->name());
	  $sql = "UPDATE result_set set feature_type_id=".$ftype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
	}
	elsif($rset->feature_type->dbID ne $ftype->dbID()){
	  $self->log('WARNING: FeatureType mismatch.  Updating IMPORT FeatureType('.$rset->feature_type->name().') to match meta('.$ftype->name.')');
	  $sql = "UPDATE result_set set feature_type_id=".$ftype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';

	}
  }

  $self->db->dbc->do($sql) if $sql;

  undef $sql;

  if(scalar keys %{$types{'cell'}} >1){
	$self->log('Resetting IMPORT CellType to NULL for multi-CellType Experiment');
	my $sql = "UPDATE result_set set cell_type_id='NULL' where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
  }else{
	my ($ctype) = values %{$types{'cell'}};

	if(! defined $rset->cell_type()){
	  $self->log('Updating IMPORT CellType to '.$ctype->name());
	  $sql = "UPDATE result_set set cell_type_id=".$ctype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
	}
	elsif($rset->cell_type->dbID ne $ctype->dbID()){
	  $self->log('WARNING: CellType mismatch.  Updating IMPORT CellType('.$rset->cell_type->name().') to match meta('.$ctype->name.')');
	  $sql = "UPDATE result_set set cell_type_id=".$ctype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
	}
  }

  $self->db->dbc->do($sql) if $sql;

  ### Generate new top level sets here based on br type combos
  #we risk duplicating sets here if import set is set to one cell/featuretype
  #duplicate anyway, as import is really just for easy tracking of all chips during import

  my %toplevel_sets;
  my $toplevel_cnt = 1;
  #could tidy up toplevel_sets implmentation

  foreach my $new_rset(values %rsets){
	
	my $ftype_name = (defined $new_rset->{'feature_type'}) ? $new_rset->{'feature_type'}->name() : undef;
	my $ctype_name = (defined $new_rset->{'cell_type'}) ? $new_rset->{'cell_type'}->name() : undef;

	if(! exists $toplevel_sets{$ftype_name}){
	  $toplevel_sets{$ftype_name} = {};
	  $toplevel_sets{$ftype_name}{'feature_type'} = $new_rset->{'feature_type'};
	}



	if(! exists $toplevel_sets{$ftype_name}{$ctype_name}){
	  $toplevel_sets{$ftype_name}{$ctype_name}{'cell_type'} = $new_rset->{'cell_type'};
	  $toplevel_sets{$ftype_name}{$ctype_name}{'rsets'} = [$new_rset];
	}else{
	  push @{$toplevel_sets{$ftype_name}{$ctype_name}{'rsets'}}, $new_rset;
	}
  }



  #build toplevel sets for each feature/cell type combo using constituent rsets
  foreach my $ftype_name(keys %toplevel_sets){
	
	foreach my $ctype_name(keys %{$toplevel_sets{$ftype_name}}){
	  
	  next if $ctype_name eq 'feature_type';#skip feature type

	  #we need to give these a different key so we're not overwriting in the rset hash
	  $rsets{$self->experiment->name().'_'.$toplevel_cnt} = Bio::EnsEMBL::Funcgen::ResultSet->new
		(
		 -NAME       => $self->experiment->name(),
		 -ANALYSIS   => $rset->analysis(),
		 -TABLE_NAME => 'experimental_chip',
		 -FEATURE_TYPE => $toplevel_sets{$ftype_name}{'feature_type'},
		 -CELL_TYPE    => $toplevel_sets{$ftype_name}{$ctype_name}{'cell_type'},
		);

	  $self->log("Created toplevel ResultSet for:\t". $rsets{$self->experiment->name().'_'.$toplevel_cnt}->log_label);

	  #add consituent table ids
	  foreach my $new_rset(@{$toplevel_sets{$ftype_name}{$ctype_name}{'rsets'}}){
		
		foreach my $ec_id(@{$new_rset->table_ids()}){

		  #Only add it if it has not already been added
		  if(!  $rsets{$self->experiment->name().'_'.$toplevel_cnt}->get_chip_channel_id($ec_id)){
			$rsets{$self->experiment->name().'_'.$toplevel_cnt}->add_table_id($ec_id, $new_rset->get_chip_channel_id($ec_id));
		  }
		}
	  }
	  $toplevel_cnt++;
	}
  }

  #ResultSet update strategy
  #To avoid messyness in resolving result_set differences
  #Simply delete all that are not used as supporting sets
  #and load new ones, log old supporting rsets for manual
  #reassignment and rollback.
  #If we have clash between an old set and a new set, rename old
  #set and log
  #We might not always have the previous data files.
  #But we might want to maintain all the previous rsets and just add a new one
  #At present this would require acquiring the previous Tab2Mage file
  #and adding the new data to it.
  #We could do with a way to merge data already in the DB with new meta data to form a new Tab2Mage file
  #and validate that
  

  my @previous_rep_sets;
  my @supporting_rset_dsets;


  #Get non-import Sets
  map {push @previous_rep_sets, $_ if $_->name !~ /_IMPORT$/} 
	@{$rset_adaptor->fetch_all_by_Experiment_Analysis($self->experiment, $chip_anal)};
 
	  
  #rollback_ResultSet if possible?
	  #This is just checking if they are supporting, not actually rolling them back
  if(@previous_rep_sets){
	$self->log('Found previously stored ResultSets');

	foreach my $prev_rset(@previous_rep_sets){
	  #This should not rollback anything, just return skipped sets
	  #i.e. sets which have a product feature set
	  #It also used to delete the supporting set records which maybe important for redefining the DataSet below
	  my $rset_dset = $self->rollback_ResultSet($prev_rset);
	  push @supporting_rset_dsets, $rset_dset if @$rset_dset;
	}
  }

  #Note: If we remove chips from an experiment, they are only removed from the non-import sets
  #To fully remove them, you need to use the rollback_experiment.pl script with -chip_ids
  #can we log this in get_import_ResultSet?

  $self->log('Storing ResultSets');

  #Store new tech, biol and toplevel type rsets
  foreach my $new_rset(values %rsets){
	my $replace_txt;

	#Rename old set if we have a name/anal/type clash
	foreach my $prs(@supporting_rset_dsets){
	
	  my ($pset, $dset) = @$prs;

	  if($pset->log_label eq $new_rset->log_label){
		my $new_name = "OLD_".$rset->log_label;
		$self->log("Found update supporting ResultSet clash, renaming to:\t${new_name}");
		$self->unlink_ResultSet_DataSet($rset, $dset, $new_name);

		#This pset dbID has already been removed
		#Will get updated with new rset dbID when updating DataSet
		$replace_txt = 'Proposed ResultSet(dbID) replacement for DataSet('.$dset->name."):\t".$pset->dbID.' > ';
	  }
	}


	$new_rset->add_status('DAS_DISPLAYABLE');
	my ($new_rset) = @{$rset_adaptor->store($new_rset)};

	if(defined $replace_txt){
	  $self->log($replace_txt.$new_rset->dbID);
	}
  }

  my $xml_file = open_file($self->get_config('mage_xml_file'));

  #slurp in changing separator to null so we get it all in one string.
  $self->experiment->mage_xml(do{ local ($/); <$xml_file>});
  close($xml_file);

  $self->experiment($self->db->get_ExperimentAdaptor->update_mage_xml_by_Experiment($self->experiment()));
  	}
  }
  
  return;
}


1;
