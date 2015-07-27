#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::Nimblegen
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

Bio::EnsEMBL::Funcgen::Parsers::Nimblegen

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::Nimblegen";
  push @INC, $parser_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a parser class which should not be instatiated directly, it 
normally set by the Importer as the parent class.  Nimblegen contains meta 
data and methods specific to NimbleGen arrays to aid parsing and importing of 
experimental data.

=cut

package Bio::EnsEMBL::Funcgen::Parsers::Nimblegen;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( species_chr_num open_file );
use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Channel;

use base qw( Bio::EnsEMBL::Funcgen::Parsers::MAGE );

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Nimblegen class
  Returntype : Bio::EnsEMBL::Funcgen::Parser::Nimblegen
  Exceptions : throws if Experiment name not defined or if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new(@_);

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
	
  #should we provide args override for all of these?

  $self->{'config'} =  
    {(
      #order of these data arrays is important!
	  #Remove these method arrays, snd just run them serially?
      array_data   => ['experiment'],#Rename this!!
      probe_data   => ["probe"],
      results_data => ["and_import_results"],
      sample_key_fields => ['DESIGN_ID', 'CHIP_ID', 'DYE', 'PROMOT_SAMPLE_TYPE'],# 'SAMPLE_LABEL'],label now optional
      # 'SAMPLE_DESCRIPTION removed due to naming disparities
      ndf_fields      => ['CONTAINER', 'PROBE_SEQUENCE', 'MISMATCH','FEATURE_ID', 'PROBE_ID'],#MISMATCH is always 0!
      pos_fields      => ['CHROMOSOME', 'PROBE_ID', 'POSITION', 'COUNT'],
      result_fields   => ['PROBE_ID', 'PM', 'X', 'Y'],
      notes_fields   => ['DESIGN_ID', 'DESIGN_NAME', 'DESCRIPTION'],
      norm_method => 'VSN_GLOG',
      dye_freqs => {(
		     Cy5 => 635,
		     Cy3 => 532,
		    )},


	  #Need to make these definable?
	  #have protocolfile arg and just parse tab2mage protocol section format
	  protocols => {(
					 grow          => {(
										accession => 'GROW_NIMB',
										name      => 'GROW NIMBLEGEN CULTURE CONDITIONS',
										text      => 'Nimblegen culture conditions description here. Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
										paramters => undef,
									   )},
					 treatment     => {(
										accession => 'CROSSLINK_NIMB',
										name      => 'NIMBLEGEN CHROMATIN PREPARATION',
										text      => 'Nimblegen X-linking and DNA extraction protocol.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
										paramters => undef,
									   )},
					 extraction    => {(
										accession => 'CHROMATIN_IP_NIMB',
										name      => 'NIMBLEGEN CHROMATIN IMMUNOPRECIPITATION and DNA RECOVERY',
										text      => 'Nimblegen chromatin immunoprecipitation and DNA extraction protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
										paramters => undef,
									   )},
					 labeling      => {(
										accession => 'LABELLING_NIMB',
										name      => 'NIMBLEGEN LABELLING',
										text      => 'Nimblegen labelling protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
										paramteres => undef,
									   )},
					 hybridization => {(
										accession => 'HYBRIDISATION_NIMB',
										name      => 'NIMBLEGEN HYBRIDISATION',
										text      => 'Nimblegen chip hybridisation protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
										parameters => undef,
									   )},
					 scanning      => {(
										accession => 'SCANNING_NIMB',
										name      => 'NIMBLESCAN',
										text      => 'Nimblegen Nimblescan protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
										paramters => undef,
									   )},
					)},

     )};
	
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

  #This should be general config for all types of import
  #dirs are not set in config to enable generic get_dir method access
  #This is really just setting paths rather than config rename?

  #This is generic for all imports



  if($self->{'old_dvd_format'}){
	$self->{'design_dir'} = $self->get_dir('input').'/DesignFiles';
  }else{
	$self->{'design_dir'} = $self->get_dir('input').'/Design_information';
  }

  
  if($self->{'old_dvd_format'}){
	$self->{'config'}{'notes_file'} = $self->get_dir('input').'/DesignNotes.txt';
  }else{
	$self->{'config'}{'notes_file'} = $self->get_dir('design').'/DesignNotes.txt';
  }
  
  $self->{'config'}{'chip_file'} = $self->get_dir('input').'/SampleKey.txt';



  #Experiment(output) specific
  #This should already be set in the run script
  #As we could get log write errors before we have created the output dir otherwise
  $self->{'output_dir'} ||= $self->get_dir("data").'/output/'.$self->{'param_species'}.'/'.$self->vendor().'/'.$self->name();

  $self->{'config'}{'tab2mage_file'} = $self->get_dir('output').'/E-TABM-'.$self->name().'.txt';

  $self->{'config'}{'mage_xml_file'} = $self->get_dir('output').'/{UNASSIGNED}.xml';

  if($self->{'old_dvd_format'}){
	$self->{'results_dir'} = $self->get_dir('input').'/PairData';
  }else{
	$self->{'results_dir'} = $self->get_dir('input').'/Raw_data_files';
  }
	
  return;
}




=head2 read_array_data

  Example    : $imp->read_array_data();
  Description: Parses NimbleGen DesignNotes.txt files to create and store new Arrays
  Returntype : none
  Exceptions : None
  Caller     : general
  Status     : At risk - Can this be generic? Can we force the creation of a DesignNotes file on other formats?

=cut


sub read_array_data{
  my ($self, $notes_file) = @_;


  $notes_file ||= $self->get_config('notes_file');
  my ($line, $array, $array_chip, @data, %hpos);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
  my $ac_adaptor = $self->db->get_ArrayChipAdaptor();

  #Slurp file to string, sets local delimtter to null and subs new lines
  my $fh = open_file($notes_file);
  #($design_desc = do { local ($/); <$fh>;}) =~ s/\r*\n$//;
  #close($fh);
  
  #Would be better if we only import the design info for the chips listed in the SampleKey.txt file
  #Some cyclical dependency going on here :|

  while ($line = <$fh>){
    
    $line =~ s/\r*\n//;#chump
    @data =  split/\t/o, $line;

    #We need to have a DESIGN vendor type?
    #also need to be able to set file path independently of config
    
    if($. == 1){
	  %hpos = %{$self->set_header_hash(\@data, $self->get_config('notes_fields'))};
      next;
    }

    ### CREATE AND STORE Array and ArrayChips  
    if(! defined $array ){
      #This is treating each array chip as a separate array, unless arrayset is defined
      #AT present we have no way of differentiating between different array_chips on same array???!!!
      #Need to add functionality afterwards to collate array_chips into single array
           
      #This will use a stored array if present

      $array = Bio::EnsEMBL::Funcgen::Array->new
	(
	 -NAME        => $self->array_name() || $data[$hpos{'DESIGN_NAME'}],
	 -FORMAT      => uc($self->format()),
	 -VENDOR      => uc($self->vendor()),
	 -TYPE        => 'OLIGO',
	 -DESCRIPTION => $data[$hpos{'DESCRIPTION'}],#need to trim the array chip specific description here
	);

      ($array) = @{$oa_adaptor->store($array)};

      $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							  -ARRAY_ID  => $array->dbID(),
							  -NAME      => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							  #add description?
							 );

      #This will use a stored array_chip if present
      ($array_chip) = @{$ac_adaptor->store($array_chip)};
      $array->add_ArrayChip($array_chip);
        
    }
    elsif((! $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])) && ($self->array_set())){
      
      $self->log("Generating new ArrayChip(".$data[$hpos{'DESIGN_NAME'}].") for same Array:\t".$array->name()."\n");
      
      $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							  -ARRAY_ID  => $array->dbID(),
							  -NAME        => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							 );
      
      ($array_chip) = @{$ac_adaptor->store($array_chip)};
      $array->add_ArrayChip($array_chip);
      
    }
    elsif(! $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])){
      throw("Found experiment with more than one design without -array_set");
    }
  }
  

  $self->add_Array($array);

  close($fh);
  
  return;

}


=head2 read_experiment_data

  Example    : $imp->read_array_chip_data();
  Description: Parses and imports array & experimental chip meta data/objects
  Returntype : none
  Exceptions : throws if more than one array/design found and not an "array set"
  Caller     : Importer
  Status     : At risk 

=cut

sub read_experiment_data{
  my $self = shift;

  $self->read_array_data();

  my $t2m_file = $self->init_tab2mage_export() if $self->{'write_mage'};


  my ($design_desc, $line, $tmp_uid, $channel, $echip, $sample_label);
  my ($sample_desc, %hpos, @data, %uid_reps, %did_reps, %sample_reps);
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $chan_adaptor = $self->db->get_ChannelAdaptor();
  my $br_cnt = 1;
  my $tr_cnt = 1;
  
  #Currently 1 design = 1 chip = 1 array /DVD
  #Different designs are not currently collated into a chip_set/array in any ordered manner
  #Register each design as an array and an array_chip
  #May want to group array_chips into array/chip sets by association though the API
  
 
  warn("Harcoded for one array(can have multiple chips from the same array) per experiment\n");
  my $fh = open_file($self->get_config("chip_file"));
  $self->log("Reading chip data");


  #warn "Do we need to validate each line here against the header array?";

  while ($line = <$fh>){
    next if $line =~ /^\s+\r*\n/;
    $line =~ s/\r*\n//;#chump
    @data =  split/\t/o, $line;
    #we could validate line against scalar of header array
     #ORD_ID  CHIP_ID DYE     DESIGN_NAME     DESIGN_ID       SAMPLE_LABEL    SAMPLE_SPECIES  SAMPLE_DESCRIPTION      TISSUE_TREATMENT        PROMOT_SAMPLE_TYPE
    if ($. == 1){
      %hpos = %{$self->set_header_hash(\@data, $self->get_config('sample_key_fields'))};

      #we need to set the sample description field name, as it can vary :(((
      @data = grep(/SAMPLE_DESCRIPTION/, keys %hpos);
      $sample_desc = $data[0];
 
      throw("More than one sample description(@data) in ".$self->get_config("chip_file")."\n") if(scalar @data >1);
      next;
    }
    
    #Need to handle array class here i.e. two channel arrays will have two lines
    #validate species here
    #look up alias from registry and match to self->species
    #registry may not be loaded for local installation


    ### CREATE AND STORE ExperimentalChips
    if ((! $tmp_uid) || ($data[$hpos{'CHIP_ID'}] ne $tmp_uid)){

      #Test both channels are available, i.e. the SampleKey has two TOTAL channels
      if($echip){

		for my $type('TOTAL', 'EXPERIMENTAL'){

		  my $test_chan =  $chan_adaptor->fetch_by_type_experimental_chip_id($type, $echip->dbID());
		  throw("ExperimentalChip(".$echip->unique_id().
				") does not have a $type channel, please check the SampleKey.txt file") if ! $test_chan;
		  
		}
	  }
      
      $tmp_uid = $data[$hpos{'CHIP_ID'}];
      $echip = $ec_adaptor->fetch_by_unique_id_vendor($data[$hpos{'CHIP_ID'}], 'NIMBLEGEN');
	  
      if($echip){

		if(! $self->recovery()){
		  throw("ExperimentalChip(".$echip->unqiue_id().
				" already exists in the database\nMaybe you want to recover?");
		}
	  }else{

		$echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
		  (
		   -EXPERIMENT_ID  => $self->experiment->dbID(),
		   -DESCRIPTION    => $data[$hpos{$sample_desc}],
           -FEATURE_TYPE   => $self->feature_type,
           -CELL_TYPE      => $self->cell_type,
		   -ARRAY_CHIP_ID  => $self->arrays->[0]->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])->dbID(),
		   -UNIQUE_ID      => $data[$hpos{'CHIP_ID'}],
           #-BIOLOGICAL_REPLICATE => ,
           #-TECHNICAL_REPLICATE => ,
		  );
		
		($echip) = @{$ec_adaptor->store($echip)};
		$self->experiment->add_ExperimentalChip($echip);
      }
    }


    ### CREATE AND STORE Channels
	my $type = uc($data[$hpos{'PROMOT_SAMPLE_TYPE'}]);
	my $sample_label = (! exists $hpos{'SAMPLE_LABEL'}) ? '' :  $data[$hpos{'SAMPLE_LABEL'}];


	$type = 'TOTAL' if ($type ne 'EXPERIMENTAL');
    $channel = $chan_adaptor->fetch_by_type_experimental_chip_id($type, $echip->dbID());

    if($channel){
      if(! $self->recovery()){
		throw("Channel(".$echip->unqiue_id().":".uc($data[$hpos{'PROMOT_SAMPLE_TYPE'}]).
			  " already exists in the database\nMaybe you want to recover?");
	  }else{
		#push @{$self->{'_rollback_ids'}}, $channel->dbID();
		#No point in doing this as all Channels mey be pre-registered in recovery mode
		#Hence all will be rolled back
	  }
    }else{
      #Handles single/mutli
   	  
      $channel =  Bio::EnsEMBL::Funcgen::Channel->new
		(
		 -EXPERIMENTAL_CHIP_ID => $echip->dbID(),
		 -DYE                  => $data[$hpos{'DYE'}],
		 -SAMPLE_ID            => $sample_label,
		 -TYPE                 => $type,
		);
    
      #-SPECIES              => $self->species(),#on channel/sample to enable multi-species chip/experiment
      #would never happen on one chip?  May happen between chips in one experiment
	  
      ($channel) = @{$chan_adaptor->store($channel)};

	}

	#we need to build the channel level tab2mage line here

	#For each BR there will be two sample_labels, one for each channel
	#These will be used across multiple chips.
	#If two chips the same design ID and the same sample labels, then we have a technical replicate
	#else if they have different sample labels then we have another biological replicate
	#We have a problem of associating channels to the same BR with differing sample labels
	#This is solved by checking whether the chip ID has already been registered in a BR

	#This fails if more than one sample label is used for any given BR
	#This will result in the BR being split into the number of unique sample label pairs(Experimental/Control  channel)
	#This also fails if the same sample label has been used for two different BRs

	if($self->{'write_mage'}){
	  #my $sample_name = ($sample_label eq '') ? '???' : substr($sample_label, 0, (length($sample_label)-1));
	  my $ctype_name = (defined $self->cell_type()) ? $self->cell_type->name() : '???';
	  my $ftype_name = (defined $self->feature_type()) ? $self->feature_type->name() : '???';
	  my $ctype_desc = (defined $self->cell_type()) ? $self->cell_type->description() : '???';

	  #define reps


	  #we need one to get the biorep based on the sample label and making sure the unique ID are the same
	  #we need to define the tech rep by matching the sample label and the making sure the design_id isn't already used

	  #Is this doing the BR assignment properly?

	  if(exists $sample_reps{$sample_label}){#Found chip in a previously seen BR
		#Register the BR of this chip ID
		$uid_reps{$data[$hpos{'CHIP_ID'}]}{'br'} = $sample_reps{$sample_label};

	  }
	  elsif(exists $uid_reps{$data[$hpos{'CHIP_ID'}]}){#Found the other channel
		$sample_reps{$sample_label} = $uid_reps{$data[$hpos{'CHIP_ID'}]}{'br'};
	  }
	  else{#assign new br
		$sample_reps{$sample_label} = $br_cnt;                  #Assign BR to sample label
		$uid_reps{$data[$hpos{'CHIP_ID'}]}{'br'} = $br_cnt;     #Assign BR to chip id
		$br_cnt++;
	  }



	  #Something is going awry here. The TR is not being reset for some new BRs
	  
	  if(! exists $uid_reps{$data[$hpos{'CHIP_ID'}]}{'tr'}){
		#we only assign a new tr here if this design has not been seen in any of the reps
		#i.e. we need to get the first tr which does not contain this design_id

		my $create_rep = 1;
		my $tr;
		my @chip_ids;
		my $br = 	$uid_reps{$data[$hpos{'CHIP_ID'}]}{'br'};

		foreach my $chip_id(keys %uid_reps){
		  
		  push @chip_ids, $chip_id if($uid_reps{$chip_id}{'br'} == $br);
		}

		#This is looping through all the TRs for all the design IDs
	   
		foreach my $rep(sort keys %did_reps){
		  #Doesn't exist for the given BR?
		  #So we need to get all the chip_ids for a given br
		  #Check wether it exists and wether it exists in did_reps and check wether the chip_id value is part of the BR set
		  #else we add it
																				 
		  if(! exists $did_reps{$rep}{$data[$hpos{'DESIGN_ID'}]}){
			#Not seen in a TR of this $rep yet
			$create_rep = 0;
		  }elsif(! grep(/$did_reps{$rep}{$data[$hpos{'DESIGN_ID'}]}/, @chip_ids)){
			#Not seen in this BR with this TR $rep
			$create_rep = 0;
		  }

		  if(! $create_rep){
			#Design ID not seen so add to this TR
			$did_reps{$rep}{$data[$hpos{'DESIGN_ID'}]} = $data[$hpos{'CHIP_ID'}]; #don't really need to assign this
			$tr = $rep;
			last;#do not remove this or we get wierd TR incrementing
		  }
		}

		
		if($create_rep){
		  #Get the next TR value for this given BR
		  my @trs;

		  foreach my $rep(keys %did_reps){
			
			foreach my $chip_id(values %{$did_reps{$rep}}){
			  
			  #Push TR if chip_id is present in this BR
			  push @trs, $rep if(grep(/$chip_id/, @chip_ids));
			}
		  }
		  ($tr) = sort {$b<=>$a} @trs;

		  $tr ||=0;
		  $tr++;
		  
		  #register design ID to chip ID mapping for this TR
		  $did_reps{$tr}{$data[$hpos{'DESIGN_ID'}]} = $data[$hpos{'CHIP_ID'}];
		}

		#register TR for this chip ID
		$uid_reps{$data[$hpos{'CHIP_ID'}]}{'tr'} = $tr;
	  }

	  my $br = $self->experiment->name().'_BR'.	$uid_reps{$data[$hpos{'CHIP_ID'}]}{'br'};
	  my $tr = $br.'_TR'.$uid_reps{$data[$hpos{'CHIP_ID'}]}{'tr'};


	  #File[raw]
	  my $tsm_line = $echip->unique_id().'_'.$self->get_config('dye_freqs')->{$data[$hpos{'DYE'}]}.'_pair.txt';
	  #Array[accession] # Should this be left blank for AE accession?
	  $tsm_line .= "\t".$data[$hpos{'DESIGN_ID'}];
	  #Array[serial]
	  $tsm_line .= "\t".$echip->unique_id();

	  #Protocol(s)[grow][treatment][extraction][labelling][hybridisation][scanning]
	  foreach my $protocol(sort (keys %{$self->get_config('protocols')})){
		$tsm_line .= "\t".$self->get_config('protocols')->{$protocol}->{'accession'};
	  }

	  
	  #BioSource
	  $tsm_line .= "\t$ctype_name";
	  #Sample
	  $tsm_line .= "\t$br";
	  #Extract 
	  $tsm_line .= "\t$tr";

	  #LabeledExtract & Immunoprecipitate
	  if($type eq 'EXPERIMENTAL'){
		$tsm_line .= "\t$sample_label - IP of $tr with anti $ftype_name (Ab vendor, Ab ID)";
		$tsm_line .= "\t$tr IP";
	  }else{
		$tsm_line .= "\t$sample_label - Input control DNA of $tr\t";
	  }
		
	  #Hybridization	
	  #U2OS BR1_TR1 ChIP H3KAc 46092 hyb
	  $tsm_line .= "\t$ctype_name $tr ChIP $ftype_name ".$echip->unique_id().' hyb';
	  
	  #BioSourceMaterial    SampleMaterial	ExtractMaterial	LabeledExtractMaterial
	  $tsm_line .= "\tcell\tgenomic_DNA\tgenomic_DNA\tsynthetic_DNA";
	  
	  #Dye	
	  $tsm_line .= "\t".$data[$hpos{'DYE'}];

	  #BioMaterialCharacteristics[Organism]
	  $tsm_line .= "\t".$self->species();
	  
	  #BioMaterialCharacteristics[BioSourceType]
	  $tsm_line .= "\tfrozen_sample";

	  #BioMaterialCharacteristics[StrainOrLine]	
	  $tsm_line .= "\t$ctype_name";

	  #BioMaterialCharacteristics[CellType]
	  $tsm_line .= "\t$ctype_name";

	  #BioMaterialCharacteristics[Sex]
	  $tsm_line .= "\t???";
	  #FactorValue[StrainOrLine]	
	  $tsm_line .= "\t$ctype_name";
	  #FactorValue[Immunoprecipitate]
	  $tsm_line .=  ($type eq 'EXPERIMENTAL') ? "\tanti-${ftype_name} antibody\n" : "\t\n";

	  print $t2m_file $tsm_line;

	}
   
  }
  
  close($t2m_file) if $self->{'write_mage'};
  close($fh);

  return;
}




=head2 read_probe_data

  Example    : $imp->read_probe_data();
  Description: Parses and imports probes, probe sets and features of a given array
               No duplicate handling or probe caching is performed due to memory 
               issues, this is done in resolve_probe_data.
  Returntype : none
  Exceptions : none
  Caller     : Importer
  Status     : Medium

=cut


#Assumes one chip_design per experimental set.
sub read_probe_data{
  my ($self) = shift;
  
  my ($fh, $line, @data, %hpos, %probe_pos);#, %duplicate_probes);
  $self->log("Parsing and importing ".$self->vendor()." probe data (".localtime().")", 1);

  ### Read in
  #ndf file: probe_set, probe and probe_feature(.err contains multiple mappings)
  #pos file: probe chromosome locations
  
  #Need to change how probe_names are generated for nimblegen?
  #native probe_ids may not be unique, but should be when combined with the seq_id which is currently being used as the xref_id
  #Handle with API!!
  
  #READ REGION POSITIONS
  #We need to handle different coord systems and possibly different assmemblies
  my $slice_a = $self->db->get_SliceAdaptor();
  my $cs = $self->db->get_FGCoordSystemAdaptor()->fetch_by_name('chromosome');


  #TIED FILE CACHE!!!
  #We need to rebuild the cache from the DB before we start adding new probe info
  #We only need to rebuild cache if we find a chip that hasn't been imported?
  #No, we just need to import without cache, then re-do the resolve step
  #Are we still going to get disconnects when we dump the cache?

  
  #warn "Read probe data should only read in the array chips which are specified by the ExperimentalChip?  Not just what is present in the DesignNotes file?";

  
  foreach my $array(@{$self->arrays()}){
    
    foreach my $achip(@{$array->get_ArrayChips()}){

      my (@log, %probe_pos, $fasta_file, $f_out);
	  #do we need to fetch probe by seq and array?
	  #this would also id non-unique seqs in design

      #warn "We need to account for different cs feature amppings here

	  if($achip->has_status('IMPORTED')){
		$self->log("Skipping fully imported ArrayChip:\t".$achip->design_id());
		next;
      }elsif($self->recovery()){
		$self->log("Rolling back ArrayChip:\t".$achip->design_id());
		$self->rollback_ArrayChips([$achip]);
      }
      
      $self->log("Importing ArrayChip:".$achip->design_id());
       
      #Always use pos file, ndf file cannot be guranteed to contain all location info
      #pos file also gives a key to which probes should be considered 'EXPERIMENTAL'
      
      #CACHE PROBE POSITIONS
      $fh = open_file($self->get_dir("design")."/".$achip->name().".pos");
      
      #don't % = map ! Takes a lot longer than a while ;)
      while($line = <$fh>){
		$line =~ s/\r*\n//o;#Not using last element
		@data =  split/\t/o, $line;
	
		#SEQ_ID	CHROMOSOME	PROBE_ID	POSITION	COUNT
		if ($. == 1){
		  %hpos = %{$self->set_header_hash(\@data, $self->get_config('pos_fields'))};
		  next;
		}
		
		#Skip probe if there is a duplicate
		if(exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){
		  
		  if($data[$hpos{'CHROMOSOME'}] eq $probe_pos{$data[$hpos{'PROBE_ID'}]}->{chr} &&
			 ($data[$hpos{'POSITION'}]+1) eq $probe_pos{$data[$hpos{'PROBE_ID'}]}->{start}){
			#log or warn here?
			
			#Not matching probe length here
			
			next;
			#Do we need to skip this in the ndf file too?

		  }
		  else{
			throw("Found duplicate mapping for ".$data[$hpos{'PROBE_ID'}].
				  " need implement duplicate logging/cleaning");
		  }
		  #need to build duplicate hash to clean elements from hash
		  # $duplicate_probes{$data[$hpos{'PROBE_ID'}]} = 1;
		  #next;
		}
		
		
		my $random = 0;

		if(! $self->cache_slice($data[$hpos{'CHROMOSOME'}])){
		  push @log, "Skipping feature import for probe ".$data[$hpos{'PROBE_ID'}]." with non-standard region ".$data[$hpos{'CHROMOSOME'}];

		  #should we try and resolve the random chrs here?
		  #at least store probe/set/result and skip feature
		  #can we just import as chr with no start postition?
		  
		  #if ($data[$hpos{'CHROMOSOME'}] =~ /_random/){
		  
		  #	if(! $self->cache_slice($data[$hpos{'CHROMOSOME'}])){
		   #push @log, "Skipping probe ".$data[$hpos{'PROBE_ID'}]." with non-standard region ".$data[$hpos{'CHROMOSOME'}];
			#}else{
		  #we should really log this in a seprate file to avoid overloading the lgo file
		  # push @log, "Importing random probe ".$data[$hpos{'PROBE_ID'}]." on ".$data[$hpos{'CHROMOSOME'}]." omitting position";
		  
		#}

		#}
		  
		}
		
		#This is not handling probes with random chrs

		$probe_pos{$data[$hpos{'PROBE_ID'}]} = {(
												 chr => $data[$hpos{'CHROMOSOME'}],
												 start => ($data[$hpos{'POSITION'}] +1),#default UCSC->Ensembl coord conversion
												)};
		
      }


	  
      
      #Remove duplicate probes 
	  
	  $self->log("Built position cache from : ".$achip->name().".pos", 1);
      close($fh);
      
	  $self->log("Importing design probes from : ".$achip->name().".ndf");
      #OPEN PROBE IN/OUT FILES
      $fh = open_file($self->get_dir("design")."/".$achip->name().".ndf");
      #Need to set these paths in each  achip hash, file names could be tablename.chip_id.txt
 
	  #Need to add dbname/port/species/vendor to this path?
	  #.efg DropDatabase should also clean the fasta dumps and caches for a given DB
	  
	  if($self->dump_fasta()){
		$fasta_file = $self->get_dir('fastas').'/'.$achip->name().".fasta";
		$self->backup_file($fasta_file);
		$f_out = open_file($fasta_file, '>');
	  }


      my ($length, $ops, $op, $of, %pfs);

      #should define mapping_method arg to allows this to be set to LiftOver/EnsemblMap
      my $anal = $self->db->get_AnalysisAdaptor()->fetch_by_logic_name("VendorMap");
			
		
      my $strand = 0;	#default for nimblegen, should be config hash?

      #my $cig_line = "50M";	#default for nimblegen, should be config hash?
	  #probe length can change within design, should be built from length


      my $fasta = "";
      
      #$self->Timer()->mark("Starting probe loop");
      

	  #This is leaking about 30-60MB for each normal density chip?
	  #need Devel::Monitor here?


      while($line = <$fh>){
		$line =~ s/\r*\n//;
		@data =  split/\t/o, $line;
		my $loc = "";
		my $class = "EXPERIMENTAL";
		
		#PROBE_DESIGN_ID	CONTAINER	DESIGN_NOTE	SELECTION_CRITERIA	SEQ_ID	PROBE_SEQUENCE	MISMATCH	MATCH_INDEX	FEATURE_ID	ROW_NUM	COL_NUM	PROBE_CLASS	PROBE_ID	POSITION	DESIGN_ID	X	Y
		#2067_0025_0001  BLOCK1          0       chrX    TTAGTTTAAAATAAACAAAAAGATACTCTCTGGTTATTAAATCAATTTCT      0       52822449        52822449        1       25      experimental    chrXP10404896   10404896        2067    25      1
		
		if ($. == 1){	
		  %hpos = %{$self->set_header_hash(\@data, $self->get_config('ndf_fields'))};
		  next;
		}
		
		
		if (! exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){
		  push @log, "Skipping non-experimental probe:\t".$data[$hpos{'PROBE_ID'}];
		  next;
		}
		
		#Which non-experimental probes might we want to store?
		#if($data[$hpos{'CONTAINER'}] =~ /control/io){
		#	$class = "CONTROL";
		#}
		#elsif($data[$hpos{'CONTAINER'}] =~ /random/io){
		#	$class = "RANDOM";
		#}
		#elsif($data[$hpos{'PROBE_CLASS'}] !~ /experimental/io){
		#	$class = "OTHER";
		#}
		#elsif(! exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){	#HACKY HACKY HACK HACK!! Needed for valid region retrival
		#  $class = "OTHER";
		#}
		#SPIKE INS?
		
		
		
		#This assumes all probes in feature/probeset are next to each other!!!!!!!!!!
		
		
		if($data[$hpos{'FEATURE_ID'}] != $data[$hpos{'MATCH_INDEX'}]){#Probe set data
		  #print "Generating new probeset:\tFeature id:\t".$data[$hpos{'FEATURE_ID'}]."\tmatchindex:\t".$data[$hpos{'MATCH_INDEX'}]."\n";
		  
		  if($ops && ($data[$hpos{'FEATURE_ID'}] ne $ops->name())){
			#THis is where we chose to update/validate
			#Do we need to pass probes if they're already stored..may aswell to reduce mysql load?
			#No point as we have to query anyway
			$self->store_set_probes_features($achip->dbID(), \%pfs, $ops);
			throw("ops still defined in caller") if defined $ops;
		  }
		  
		  $ops = Bio::EnsEMBL::Funcgen::ProbeSet->new(
													  -NAME => $data[$hpos{'FEATURE_ID'}],
													  -SIZE => undef,
													  -FAMILY  => $data[$hpos{'CONTAINER'}],
													  #xref_id => $data[$hpos{'SEQ_ID'}],#Need to populate xref table
													 );
		  
		  #should we store straight away or build a probeset/probe/feature set, and then store and validate in turn?
		  #Store directly have separate method to validate and update?
		  #would need to check if one exists before storing anyway, else we could potentially duplicate the same probe/probeset from a different array
		  #remember for affy we need duplicate probe records with identical probe ids, probeset records unique across all arrays
		  
		  undef %pfs
		}
		elsif($. > 2){#may have previous ops set, but next has no ops, or maybe just no ops's at all
		  $self->store_set_probes_features($achip->dbID(), \%pfs, $ops);
		  throw("ops still defined in caller") if defined $ops;
		}
	
		
		###PROBES
		#should we cat $xref_id to $probe_id here to generate unique id?
		#would be messy to handle in the code, but would have to somewhere(in the retrieval code)
		
		$length = length($data[$hpos{'PROBE_SEQUENCE'}]);	
		#$probe_string .= "\t${psid}\t".$data[$hpos{'PROBE_ID'}]."\t${length}\t$ac_id\t${class}\n";
		
		$op = Bio::EnsEMBL::Funcgen::Probe->new(
												-NAME          => $data[$hpos{'PROBE_ID'}],
												-LENGTH        => $length,
												-ARRAY         => $array,
												-ARRAY_CHIP_ID => $achip->dbID(),
												-CLASS         => $class,
											   );
		
		
		%{$pfs{$data[$hpos{'PROBE_ID'}]}} = (
											 probe => $op,
											 features => [],
											);
		
				
		###PROBE FEATURES
		#How can we be certain that we have the same mapping in the DB?
		#Put checks in here for build?
		#Need to handle controls/randoms here
		#won't have features but will have results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#The format of the pos file looks like it should have all the data required, but 
		# chromsome is missing, first undef :(
		#Need to use $pos here instead of .pos file
		#However have problems defining probe class, as not populated in test set
		#derive from container! :(
		#Ignore controls/random as they won't have a region
		#Also need to handle multiple mappings? 
		#As results reference probes, no features, can have multiple features on different builds
		
		#if($class eq "EXPERIMENTAL"){
		
		#if(exists $regions{$data[$hpos{'SEQ_ID'}]}){
		#  $fid++;
		#  $pf_string .= "\t".$regions{$data[$hpos{'SEQ_ID'}]}{'seq_region_id'}."\t".$data[$hpos{'POSITION'}]."\t".
		#    ($data[$hpos{'POSITION'}] + $length)."\t${strand}\t".$regions{$data[$hpos{'SEQ_ID'}]}{'coord_system_id'}.
		#	"\t${pid}\t${anal_id}\t".$data[$hpos{'MISMATCH'}]."\t${cig_line}\n";
		#$loc .= $regions{$data[$hpos{'SEQ_ID'}]}{'seq_region_id'}.":".$data[$hpos{'POSITION'}].
		#  "-".($data[$hpos{'POSITION'}] + $length).";" if ($self->{'_dump_fasta'});
		# }
		# else{ 
		#   die("No regions defined for ".$data[$hpos{'SEQ_ID'}]." ".$data[$hpos{'PROBE_ID'}].
		#" with family ".$data[$hpos{'CONTAINER'}]);
		#  }
		
	
				
	
		if ($self->dump_fasta()){
		  #(my $chr = $probe_pos{$data[$hpos{'PROBE_ID'}]}->{'chr'}) =~ s/chr//;

		  #$loc .= $chr.":".$probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'}."-".
		#	($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'}+ $length).";";
	
		  #$fasta .= ">".$data[$hpos{'PROBE_ID'}]."\t".$data[$hpos{'CHROMOSOME'}].
		  #	"\t$loc\n".$data[$hpos{'PROBE_SEQUENCE'}]."\n";



		  #filter controls/randoms?  Or would it be sensible to see where they map
		  #wrap seq here?
		  #$fasta .= ">".$data[$hpos{'PROBE_ID'}]."\n".$data[$hpos{'PROBE_SEQUENCE'}]."\n";


		  #To use this for mapping, we really need the dbID nr fasta
		  #This can be generated after the import, or maybe during resolve?
		  #This is also currently done on a chip level, where as the cache is resolved at the array level
		  #We could simply cat the files before resolving the fasta file
		  #Need to do this otherwise we risk overwriting the fasta file with incomplete data.
		  #Can we validate sequence across probes with same name in this step?
		  #Just use probe name for now.
		  
		  #We could cat and sort the fastas to make sure we have the same sequences
		  #Need to dump the design_id in the fasta header
		  #This would also reduce IO on the DB as identical probe will be consecutive, hence just one query to get the id.

		  #Changed th format an content of this to facilitate dbID nr fasta file generation and sequence validation
		  
		  $fasta .= ">".$data[$hpos{'PROBE_ID'}]."\t".$achip->design_id."\n".$data[$hpos{'PROBE_SEQUENCE'}]."\n";

		  #Print fasta every 10000 lines
		  if(! ($. % 10000)){
			print $f_out $fasta;
			$fasta = '';
		  }
		}


		#Hack!!!!!! Still importing probe (and result?)
		next if(!  $self->cache_slice($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'chr'}));
		#warn("Skipping non standard probe (".$data[$hpos{'PROBE_ID'}].") with location:\t$loc\n");
		  				
		
		$of = Bio::EnsEMBL::Funcgen::ProbeFeature->new
		  (
		   -START         => $probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'},
		   -END           =>($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'}  + $length),
		   -STRAND        => $strand,
		   -SLICE         => $self->cache_slice($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'chr'}),
		   -ANALYSIS      => $anal,
		   -MISMATCHCOUNT => $data[$hpos{'MISMATCH'}],#Is this always 0 for import? remove from header hash?
		   -PROBE         => undef,	#Need to update this in the store method
		  );
				
	
		push @{$pfs{$data[$hpos{'PROBE_ID'}]}{'features'}}, $of;
	
      }
      
      #need to store last data here
      $self->store_set_probes_features($achip->dbID(), \%pfs, $ops);
      $self->log(join("\n", @log));
      $achip->adaptor->store_status("IMPORTED", $achip);
	  $self->log("ArrayChip:\t".$achip->design_id()." has been IMPORTED");
      
      if ($self->dump_fasta()){
		print $f_out $fasta;
		close($f_out);
      }
    
	  $self->log("Imported design from:\t".$achip->name().".ndf", 1);
		


	  #$self->{'_probe_cache'} = undef;#As we can't get Y and Y info from the DB, this is only possible as the results files contain X and Y info
	}
	
    #Should we build hash of probe_names:probe_feature_ids here for results import
    #Should we dump this as a lookup file for easier recoverability
    #This is the biggest step in the import
    #Building the hash would be fastest but least recoverable
    #Would have to write recover statments for each step i.e. build the most recent data structure required for the next import step
    #Individual queries for each result would take ages
    #This is all assuming there are no random records in the table i.e. ID series is linear with no gaps.
    
    
    #Could throw a random validation check in every X entries?
    #This would only work for non-parallel imports
    #periodic import when hit new probe_set, but no new data printed
		
  }
  
  
  $self->log("Finished parsing probe data");
  #Total probe_sets:\t$psid\n".
  #	       "Total probes:\t$pid\nTotal probe_features:\t$fid");
 

  $self->resolve_probe_data();
 
  return;
}




=head2 read_and_import_results_data

  Example    : $imp->read_results_data();
  Description: Parses and dumps raw results to file
  Returntype : none
  Exceptions : none
  Caller     : Importer
  Status     : at risk 

=cut

 
sub read_and_import_results_data{
  my $self = shift;
  
  $self->log("Parsing ".$self->vendor()." results");
  my (@header, @data, @design_ids, @lines);
  my ($fh, $pid, $line, $file);
  my $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("RawValue");
  my $result_set = $self->get_import_ResultSet($anal, 'channel');


  if ($result_set) {			#we have some new data to import
    
    foreach my $echip (@{$self->experiment->get_ExperimentalChips()}) {
	  
      #if( ! $echip->has_status('IMPORTED')){
		
	  foreach my $chan (@{$echip->get_Channels()}) {
		
		if ( ! $chan->has_status('IMPORTED')) {

		  #we need to set write_mage here

		  my $array = $echip->get_ArrayChip->get_Array();
		  
		  $self->get_probe_cache_by_Array($array) || throw('Failed to get the probe cache handle for results import, resolve cache here?');			
		  my ($probe_elem, $score_elem, %hpos);
		  my $cnt = 0;
		  my $r_string = "";
		  my $chan_name = $echip->unique_id()."_".$self->get_config('dye_freqs')->{$chan->dye()};
		  my $cc_id = $result_set->get_chip_channel_id($chan->dbID());
		  
		  
		  #if ($self->recovery()) {
		#	$self->log("Rolling back results for channel:\t${chan_name}");
		#	$self->db->rollback_results($cc_id);
		#  }
		  
		  #open/backup output
		  my $out_file = $self->get_dir("raw")."/result.".$chan_name.".txt";	
		  $self->backup_file($out_file);
		  my $r_out = open_file($out_file, '>');
		  
		  (my $alt_chan_name = $chan_name) =~ s/\_/\_1h\_/;
		  my $found = 0;
		  
		FILE: foreach my $name($chan_name, $alt_chan_name){
			
			foreach my $suffix ("_pair.txt", ".pair", ".txt") {
			  
			  $file = $self->get_dir("results")."/".$name.$suffix;
			  
			  if (-f $file) {
				$found = 1;
				last FILE;
			  }
			}
		  }
		  
		  throw("Could not find result file for Channel(${chan_name}) in ".$self->get_dir('results')) if ! $found;
		  
		  #open/slurp input
		  $self->log("Reading result for channel $chan_name:\t$file", 1);
		  $fh = open_file($file);	
		  @lines = <$fh>;
		  close($fh);
		  
					
		  ###PROCESS HEADER
		  
		  foreach my $i (0..$#lines) {
			
			if ($lines[$i] =~ /PROBE_ID/o) {
			  $lines[$i] =~ s/\r*\n//o;
			  @data = split/\t/o, $lines[$i];
			  
			  %hpos = %{$self->set_header_hash(\@data, $self->get_config('result_fields'))};
			  
			  #remove header
			  splice @lines, $i, 1;
			  
			  last;				#finished processing header
			}
		  }
		  
		  #we need to sort the result files based on the unique key(name at present, should replace with seq at some point)
		  @lines = sort {(split/\t/o, $a)[$hpos{'PROBE_ID'}] cmp (split/\t/o, $b)[$hpos{'PROBE_ID'}]} @lines;
		  
		  $self->log('Parsing results', 1);
		  
		  foreach $line(@lines) {
			
			#can we preprocess effectively?
			next if $line =~ /^#/;
			next if $line =~ /NGS_CONTROLS/;
			next if $line =~ /V_CODE/;
			next if $line =~ /H_CODE/;
			next if $line =~ /RANDOM/;
			
			$line =~ s/\r*\n//o;
			@data = split/\t/o, $line;
			
			###PROCESS HEADER
			#if ($line =~ /PROBE_ID/o){
			#  
			#  %hpos = %{$self->set_header_hash(\@data, $self->get_config('result_fields'))};
			#  next;#finished processing header
			#}
			
			###PROCESS DATA
			#Is this string concat causing the slow down, would it befaster to use an array and print a join?
			
			if ($pid = $self->get_probe_id_by_name_Array($data[$hpos{'PROBE_ID'}], $array)) {
			  $cnt ++;
			  $r_string .= '\N'."\t${pid}\t".$data[$hpos{'PM'}]."\t${cc_id}\t".$data[$hpos{'X'}]."\t".$data[$hpos{'Y'}]."\n";
			} else {
			  warn "Found unfiltered non-experimental probe in input $data[$hpos{'PROBE_ID'}]";
			}
			
			###PRINT SOME RESULTS
			if ($cnt > 10000) {
			  $cnt = 0;
			  print $r_out $r_string;
			  $r_string ="";
			  #could we fork here and import in the background?
			}
			
		  } 
		  #PRINT/CLOSE Channel file
		  print $r_out $r_string;
		  close($r_out);
		  $self->log("Finished parsing $chan_name result", 1);
		  
		  #Import directly here to avoid having to reparse all results if we crash!!!!
		  $self->log("Importing:\t$out_file");
		  $self->db->load_table_data("result",  $out_file);
		  $self->log("Finished importing:\t$out_file", 1);
		  $chan->adaptor->store_status('IMPORTED', $chan);
		  
		  
		}
	  }
	}
  } 
  else {
	$self->log("Skipping results parse and import");
  }
  
  $self->log("Finished parsing and importing results");
  
  return;
}
  


1;
