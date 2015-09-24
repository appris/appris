#
# Ensembl module for Bio::EnsEMBL::Funcgen::ExperimentalChip
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

Bio::EnsEMBL::Funcgen::ExperimentalChip - A module to represent a physical unique experimental chip.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExperimentalChip;

my $ec = Bio::EnsEMBL::Funcgen::ExperimentalChip->new(
							 -dbID           => $ec_id,
							 -unique_id      => $c_uid,
						         -experiment_id  => $exp_id,
	                                                 -array_chip_id  => $ac_id,
							 -feature_type   => $ftpye,
                                                         -cell_type      => $ctype,
                                                         -chip_set_id    => 1,
							 );


=head1 DESCRIPTION

An ExperimentalChip object represent a physical array chip/slide used in an experiment. The data
(currently the unique_id, experiment_id, array_chip_id, and description) are stored
in the experimental_chip table.

=cut


package Bio::EnsEMBL::Funcgen::ExperimentalChip;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use base qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-unique_id]     : int - the unique id of this individual experimental chip
  Arg [-experiment_id] : int - the experiment dbID
  Arg [-array_chip_id] : int - the dbID or the array_chip
  Arg [-feature_type ] : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [-cell_type ]    : Bio::EnsEMBL::Funcgen::CellType
  Arg [-biological_replicate ]    : string - the name to define the biological replicate set
  Arg [-technical_replicate ]    : string - the name to define the technical replicate set


  Example    : my $array = Bio::EnsEMBL::Funcgen::ExperimentalChip->new(
									-dbID          => $ec_id,
									-unique_id     => $c_uid,
									-experiment_id => $exp_id,
									-array_chip_id => $ac_id,
									-feature_type  => $ftype,
              	                    -cell_type     => $ftype,
	                                -biological_replicate => 'BIOREP1',
                                    -technical_replicate  => 'techrep_1',
								       );
  Description: Creates a new Bio::EnsEMBL::Funcgen::ExperimentalChip object.
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : None ? should throw if mandaotry params not set
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  
  #can we lc these?
  my ($c_uid, $exp_dbid, $ac_id, $ftype, $ctype, $brep, $trep)
    = rearrange( ['UNIQUE_ID', 'EXPERIMENT_ID', 'ARRAY_CHIP_ID', 'FEATURE_TYPE', 'CELL_TYPE', 'BIOLOGICAL_REPLICATE', 'TECHNICAL_REPLICATE'], @_ );
  
  
  $self->unique_id($c_uid)          if defined $c_uid;
  $self->experiment_id($exp_dbid)          if defined $exp_dbid;
  $self->array_chip_id($ac_id)    if defined $ac_id;
  $self->feature_type($ftype)   if defined $ftype;
  $self->cell_type($ctype)   if defined $ctype;
  $self->biological_replicate($brep) if defined $brep;
  $self->technical_replicate($trep) if defined $trep;
  return $self;
}

=head2 get_Experiment

  Args       : None
  Example    : my $exp = $exp_chip->get_Experiment();
  Description: Returns the Experiment which this ExperimentalChip belongs to.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Experiment {
  my $self = shift;
  
  if (! $self->{'experiment'}){
	
    if ($self->dbID() && $self->adaptor() ) {
	  $self->{'experiment'} = $self->adaptor->db->get_ExperimentAdaptor->fetch_by_dbID($self->experiment_id);
    } else {
      warning('Need database connection to retrieve Experiment');
    }
  }


  return $self->{'experiment'};
}



=head2 get_Channels

  Args       : None
  Example    : my $channels = $exp_chip->get_Channels();
  Description: Returns all channels on a ExperimentalChip. Needs a database connection.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Channels {
  my $self = shift;

  if (! $self->{'channels'}){

    if ($self->dbID() && $self->adaptor() ) {
      foreach my $channel (@{$self->adaptor->db->get_ChannelAdaptor->fetch_all_by_ExperimentalChip($self)}){
	$self->add_Channel($channel);
      }	
    } else {
      warning('Need database connection to retrieve Channels');
    }
  }
  
  return [values %{$self->{'channels'}}];
}


=head2 add_Channel

  Args       : Bio::EnsEMBL::Funcgen::Channel
  Example    : $exp_chip->add_channel($chan);
  Description: Sets ad channel object for the ExperimentalChip
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel objects
  Exceptions : warns if Channel already set
  Caller     : General
  Status     : At Risk

=cut

sub add_Channel{
  my ($self, $chan) = @_;


  if(! ($chan  && $chan->isa("Bio::EnsEMBL::Funcgen::Channel") && $chan->dbID())){
    throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::Channel object");
  }
  

  $self->{'channels'} ||= {};

  if (exists $self->{'channels'}->{$chan->dbID()}){
    #should this throw?
    #This currently prevents haveing to check whether a channel has already been added
    #If we were duplicating then we probably would have a different dbID
    warn("You cannot add the same Channel to an ExperimentalChip more than once");
  }else{


    ##change this to key on freq?

    $self->{'channels'}{$chan->dbID()} = $chan;
  }
  
  return;
}



=head2 get_channel_ids

  Args       : None
  Example    : my @channel_ids = @{$array->get_channel_ids()};
  Description: Returns all channel ids for an ExperimentalChip. Needs a database connection.
  Returntype : List of ints
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_channel_ids{
  my $self = shift;

  $self->get_Channels();

  return [keys %{$self->{'channels'}}];
}

=head2 get_Channel_by_dye

  Args       : string - dye used in channel
  Example    : my $chan = $echip->get_Channel_by_dye("CY5");
  Description: Returnsthe channel corresponding to the frequency specified
  Returntype : Bio::EnsEMBL::Funcgen::Channel
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Channel_by_dye{
  my ($self, $dye) = @_;

  my @chans;

  foreach my $chan(@{$self->get_Channels()}){
    push @chans, $chan if uc($chan->dye()) eq uc($dye);
  }

  throw("Found more than one Channels with the same dye") if(scalar(@chans) > 1);


  return (@chans) ? $chans[0] : undef;
}

=head2 contains_Channel

  Args [1]   : Bio::EnsEMBL::Funcgen::Channel
  Example    : if(! $echip->contains_Channel($chan){..add channel ..};
  Description: Checks whether this Channel has already been added to the ExperimentalChip
  Returntype : Boolean
  Exceptions : Throws if arg not a valid stored Bio::EnseMBL::Funcgen::Channel
  Caller     : General
  Status     : At Risk

=cut

sub contains_Channel{
  my ($self, $chan) = @_;
  
  if(! ($chan  && $chan->isa("Bio::EnsEMBL::Funcgen::Channel") && $chan->dbID())){
    throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::Channel object");
  }
  
  $self->get_Channels();

  my $contains = 0;

  $contains = 1 if(exists $self->{'channels'}->{$chan->dbID()});

  return $contains;
}



=head2 unique_id

  Arg [1]    : (optional) int - the unique chip id for this ExperimentalChip
  Example    : my $c_uid = $array->unique_id();
  Description: Getter, setter unique_id attribute.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : at Risk

=cut

sub unique_id {
  my $self = shift;
  $self->{'unique_id'} = shift if @_;
  return $self->{'unique_id'};
}

=head2 feature_type

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::FeatureType
  Example    : $ec->feature_type($ftype);
  Description: Getter/Setter thefeature_type attribute.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : Throws if arg is not a Bio::EnsEMBL::FeatureType
  Caller     : General
  Status     : At Risk

=cut

sub feature_type {
  my $self = shift;

  if(@_){
    throw("Must pass a valid Bio::EnsEMBL::Funcgen::FeatureType object") if (! $_[0]->isa("Bio::EnsEMBL::Funcgen::FeatureType"));
    $self->{'feature_type'} = shift;
  }

  return $self->{'feature_type'};
}


=head2 cell_type

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::CellType
  Example    : $ec->cell_type($ctype);
  Description: Getter/Setter the cell_type attribute.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : Throws if arg is not a Bio::EnsEMBL::CellType
  Caller     : General
  Status     : At Risk

=cut

sub cell_type {
  my $self = shift;

  if(@_){
    throw("Must pass a valid Bio::EnsEMBL::Funcgen::CellType object") if (! $_[0]->isa("Bio::EnsEMBL::Funcgen::CellType"));
    $self->{'cell_type'} = shift;
  }

  return $self->{'cell_type'};
}



=head2 biological_replicate

  Arg [1]    : (optional) string - the name or number of the chip biological replicate set
  Example    : $ec->biological_replicate('SAMPLENAME_BR1');
  Description: Getter, setter for the biological_replicate attribute.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub biological_replicate {
  my $self = shift;
  $self->{'biological_replicate'} = shift if @_;
  return $self->{'biological_replicate'};
}

=head2 technical_replicate

  Arg [1]    : (optional) string - the name or number of the chip technical replicate set
  Example    : $ec->technical_replicate('SAMPLENAME_BR1_TR1');
  Description: Getter, setter for the technical_replicate attribute.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub technical_replicate {
  my $self = shift;
  $self->{'technical_replicate'} = shift if @_;
  return $self->{'technical_replicate'};
}



=head2 experiment_id

  Arg [1]    : (optional) int - the experiment dbID
  Example    : my $exp_id = $array->experiment_id();
  Description: Getter, setter experiment_id attribute
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub experiment_id {
  my $self = shift;
  $self->{'experiment_id'} = shift if @_;
  return $self->{'experiment_id'};
}

=head2 array_chip_id

  Arg [1]    : (optional) int - the array_chip dbID
  Example    : my $ac_id = $ec->array_chip_id();
  Description: Getter, setter array_chip_id attribute
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub array_chip_id {
  my $self = shift;
  $self->{'array_chip_id'} = shift if @_; 
  return $self->{'array_chip_id'};
}

=head2 get_ArrayChip

  Example    : my $array_chip = $exp_chip->get_ArrayChip();
  Description: Getter for the array_chip attribute
  Returntype : Bio::EnsEMBL::Funcgen::ArrayChip
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_ArrayChip {
  my $self = shift;

  if(! defined $self->{'array_chip'}){
    $self->{'array_chip'} = $self->adaptor->db->get_ArrayChipAdaptor()->fetch_by_dbID($self->array_chip_id());
  }

  return $self->{'array_chip'};
}



1;

