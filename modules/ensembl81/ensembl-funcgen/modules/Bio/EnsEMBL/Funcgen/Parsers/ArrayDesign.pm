#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::ArrayDesign
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

Bio::EnsEMBL::Funcgen::Parsers::ArrayDesign

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::ArrayDesign";
  push @INC, $parser_type;
  my $imp = $class->SUPER::new(@_);  my $imp = Bio::EnsEMBL::Funcgen::Importer->new(%params);

  $imp->set_config();


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally inherited from the Importer.  ArrayDesign contains meta data and methods 
specific to handling array designs only (i.e. no experimental data), which have 
been produced from the eFG array design software.

=cut

package Bio::EnsEMBL::Funcgen::Parsers::ArrayDesign;

use strict;
use warnings;
use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Channel;
use Bio::EnsEMBL::Utils::Exception         qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( species_chr_num open_file );

use base qw( Bio::EnsEMBL::Funcgen::Utils::Helper );



=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for ArrayDesign class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::ArrayDesign
  Exceptions : throws if Experiment name not defined or if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new();

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));


  $self->{'config'} =  
    {(     
      probe_data   => ["probe"],
      prb_fields   => ['SEQ_ID',  'POSITION', 'LENGTH', 'PROBE_SEQUENCE', 'PROBE_ID', 'UNIQUENESS_SCORE', 'TM', 'MAS_CYCLES'],
      notes_fields => ['DESIGN_ID', 'DESIGN_NAME', 'DESCRIPTION'],
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
  my ($self) = @_;
  
  #placeholder method
  #set paths

  return;
}


=head2 read_array_data

  Example    : $imp->read_array_data();
  Description: Parses NimbleGen style DesignNotes.txt format files to create and store new Arrays
  Returntype : none
  Exceptions : None
  Caller     : general
  Status     : At risk - Can this be generic? Can we force the creation of a DesignNotes file on other formats?

=cut

#this is currently OLIGO specific.

sub read_array_data{
  my ($self, $design_notes) = @_;

  $self->log("Reading and importing array data");

  throw('You need to pass the path to a DesignNotes.txt file') if ! defined $design_notes;
  $self->{'design_notes'} = $design_notes;

  my ($line, $array, $array_chip, @data, %hpos);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
  my $ac_adaptor = $self->db->get_ArrayChipAdaptor();
  my $fh = open_file("<", $self->{'design_notes'});
  
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
      
      $self->log("Generating new ArrayChip(".$data[$hpos{'DESIGN_NAME'}].". for same Array ".$array->name()."\n");
      
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





=head2 read_probe_data

  Example    : $imp->read_probe_data();
  Description: Parses and imports probes, probe sets and features for a given array design
  Returntype : none
  Exceptions : throws is not tiling format
  Caller     : Importer
  Status     : at risk

=cut


#Assumes one chip_design per experimental set.
sub read_probe_data{
  my ($self, $array_file) = @_;
  
  $self->log("Reading and importing probe data");


  my ($fh, $line, @data, @log, %hpos, %probe_pos);#, %duplicate_probes);
  my $aa = $self->db->get_AnalysisAdaptor();
  my $manal = $aa->fetch_by_logic_name('MASCycles');
  my $uanal = $aa->fetch_by_logic_name('UScore');
  my $tmanal= $aa->fetch_by_logic_name('NimblegenTM');

  $array_file ||= $self->array_file();
  $self->log("Parsing ".$self->vendor()." probe data (".localtime().")");
  throw("ArrayDesign only accomodates a tiling design with no feature/probesets") if ($self->format() ne 'TILED');

  ### Read in
  # eFG prb file, not chiip info yet so only one ArrayChip per design
  # potential to have pos file here for probes built on generic slices of genome

  #We need to handle different coord systems and possibly different assmemblies
  my $slice_a = $self->db->get_SliceAdaptor();
  my $cs = $self->db->get_FGCoordSystemAdaptor()->fetch_by_name_schema_build_version(
										     'chromosome', 
										     $self->db->_get_schema_build($self->db->dnadb())
										    );


  #sanity check we're only dealing with one array/chip
  my @arrays = @{$self->arrays()};
  if(scalar(@arrays) != 1){
    throw("Array DESIGN imports only accomodate one Array per import, please check ".$self->{'design_notes'});
  }

  my @achips = @{$arrays[0]->get_ArrayChips()};
  if(scalar(@achips) != 1){
    throw("Array DESIGN imports only accomodates one ArrayChip per import, please check ".$self->{'design_notes'});
  }

  my $achip = $achips[0];

  #foreach my $array(@{$self->arrays()}){

    
  #  foreach my $achip(@{$array->get_ArrayChips()}){

  $self->log("Importing array design(".$achip->name().") from ".$array_file);


  if($achip->has_status('IMPORTED')){
    $self->log("Skipping fully imported ArrayChip:\t".$achip->design_id());
    return;
  }elsif($self->recovery()){
    $self->log("Rolling back partially imported ArrayChip:\t".$achip->design_id());
    $self->db->rollback_ArrayChip([$achip]);
  }
      
  $self->log("Importing ArrayChip:".$achip->design_id());
            
  
  #OPEN PROBE IN/OUT FILES
  $fh = open_file("<", $array_file);
  my $f_out = open_file(">", $self->get_dir("output")."/probe.".$achip->name()."fasta")	if($self->{'_dump_fasta'});
  my ($op, $of, %pfs);

  #should define mapping_method arg to allows this to be set to LiftOver/EnsemblMap
  my $anal = $self->db->get_AnalysisAdaptor()->fetch_by_logic_name("TileMap");##???
  my $strand = 0;	#default for TileMap, should be config hash?
  my $fasta = "";
       
  while($line = <$fh>){
    $line =~ s/\r*\n//;
    @data =  split/\t/o, $line;
    my $loc = "";
      
    #SEQ_ID  WINDOW_START    WINDOW_END      POSITION        LENGTH  PROBE_SEQUENCE  TM      UNIQUENESS_SCORE        MAS_CYCLES
    #X       3000001 3000100 3000041 56      TGACATCTTCAGTTCTTTACATAGTTTTCATATTAGTCCTCTATCAGATGTGGAGT        73.09   132     15

    
    if ($. == 1){	
      %hpos = %{$self->set_header_hash(\@data, $self->get_config('prb_fields'))};
      next;
    }
		  
	
    #This assumes tiling format with no feature/probe sets

		  
    if(%pfs){
      $self->store_set_probes_features($achip->dbID(), \%pfs);
      undef %pfs;
    }
	  
				
    #PROBE
    $op = Bio::EnsEMBL::Funcgen::Probe->new(
					    -NAME          => $data[$hpos{'PROBE_ID'}],
					    -LENGTH        => $data[$hpos{'LENGTH'}],
					    -ARRAY         => $arrays[0],
					    -ARRAY_CHIP_ID => $achip->dbID(),
					    -CLASS         => 'DESIGN',
					   );

    $op->add_Analysis_score($manal, $data[$hpos{'MAS_CYCLES'}]);
    $op->add_Analysis_score($tmanal, $data[$hpos{'TM'}]);
    $op->add_Analysis_CoordSystem_score($uanal, $cs, $data[$hpos{'UNIQUENESS_SCORE'}]);	
    
    #would need to pass cs to store USCORE, CYCLES and TM are seq/anal dependent not cs
    #options:
    #associate cs dependent scores with features
    #we are duplicating cs in probe_design table
    
    #do we need another object? ProbeDesign
    #-cs would be empty for all but uscore
    #-mas_cycles analysis_id
    #-uscore analysis_id cs_id
    #-tm analysis_id  (anal most likely wont change so highly redundant)
    #this would produce 3 records for each probe, with cs being empty for two and anal being redundant for the other
    #would however provide for extensible design attributes
    #would only need one tm and mas_cycles for each probe irrespective of cs
    #could have separate table probe_design_feature?
    #mmm doesn't have location, just cs
    #just have empty cs fields, calls by cs would have to be in ('cs_id', 'NULL')
    #or just make uscore method dependent on cs.
    #or split table?
    #ProbeDesign::add_analysis_score
    #ProbeDesign::add_coord_sys_analysis_score
    
    #can we add this directly to the probe?
    #separate retrieval in ProbeAdaptor so we're not joining everytime
    #this would require generic get_analysis/analysis_coord_system_attribute method
    


    %{$pfs{$data[$hpos{'PROBE_ID'}]}} = (
					 probe => $op,
					 features => [],
					);
	
				
    #PROBE FEATURE
    if(!  $self->cache_slice($data[$hpos{'SEQ_ID'}])){
      warn("Skipping non-standard probe chromosome");
      undef %pfs;
      next;
    }

    my $end = ($data[$hpos{'POSITION'}] + $data[$hpos{'LENGTH'}]);

    if ($self->{'_dump_fasta'}){
      $loc .= $data[$hpos{'SEQ_ID'}].":".$data[$hpos{'POSITION'}]."-${end};";
    }	
				
	
    $of = Bio::EnsEMBL::Funcgen::ProbeFeature->new
      (
       -START         => $data[$hpos{'POSITION'}],
       -END           => $end,
       -STRAND        => $strand,
       -SLICE         => $self->cache_slice($data[$hpos{'SEQ_ID'}]),
       -ANALYSIS      => $anal,
       -MISMATCHCOUNT => 0,
       -CIGAR_LINE    => $data[$hpos{'LENGTH'}].'M',
       -PROBE         => undef,#Need to update this in the store method
      );
				
	
    push @{$pfs{$data[$hpos{'PROBE_ID'}]}{'features'}}, $of;
	
    if($self->{'_dump_fasta'}){			
      #filter controls/randoms?  Or would it be sensible to see where they map
      #wrap seq here?
      $fasta .= ">".$data[$hpos{'PROBE_ID'}]."\t".$data[$hpos{'CHROMOSOME'}].
	"\t$loc\n".$data[$hpos{'PROBE_SEQUENCE'}]."\n";
    }
  }
      
  #need to store last data here
  $self->store_set_probes_features($achip->dbID(), \%pfs);
  $self->log(join("\n", @log));
  $achip->adaptor->set_status("IMPORTED", $achip);
  $self->log("ArrayChip:\t".$achip->design_id()." has been IMPORTED");
  
  if ($self->{'_dump_fasta'}){
    print $f_out $fasta if($self->{'_dump_fasta'});
    close($f_out);
  }
   
  $self->log("Finished parsing probe data");
  #Total probe_sets:\t$psid\n".
  #	       "Total probes:\t$pid\nTotal probe_features:\t$fid");
  
  return;
}


1;
