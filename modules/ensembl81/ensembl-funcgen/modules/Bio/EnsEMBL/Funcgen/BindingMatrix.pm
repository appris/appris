# Ensembl module for Bio::EnsEMBL::Funcgen::BindingMatrix

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

Bio::EnsEMBL::Funcgen::BindingMatrix - A module to represent a BindingMatrix. 
In EFG this represents the binding affinities of a Transcription Factor to DNA.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::BindingMatrix;

# A C G T frequency matrix
my $freqs = "1126 6975 6741 2506 7171 0 11 13 812 867 899 1332
4583 0 99 1117 0 12 0 0 5637 1681 875 4568
801 181 268 3282 0 0 7160 7158 38 2765 4655 391
661 15 63 266 0 7159 0 0 684 1858 742 880";

my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(-name         => 'MA0095.2',
                                                       -description  => "MEF2C Jaspar Matrix",
                                                       -frequencies  => $freqs
                                                       -analysis     => $jaspar_analysis,
                                                       -feature_type => $MEF2C_ftype ;

print $matrix->relative_affinity("TGGCCACCA")."\n";

print $matrix->threshold."\n";

=head1 DESCRIPTION

This class represents information about a BindingMatrix, containing the name 
(e.g. the Jaspar ID, or an internal name), and description. A BindingMatrix 
is always associated to an Analysis (indicating the origin of the matrix e.g. 
Jaspar) and a FeatureType (the binding factor).   

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor
Bio::EnsEMBL::Funcgen::MotifFeature

=cut

package Bio::EnsEMBL::Funcgen::BindingMatrix;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw( parse_matrix_line 
                                                      reverse_complement_matrix );
  
use base qw( Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  Arg [-name]       : Scalar - Name of matrix
  Arg [-analysis]   : Bio::EnsEMBL::Analysis - analysis describing how the matrix was obtained
  Arg [-frequencies]: String or Arrayref (Mandatory) - A string or a 2d array of frequencies
                      representing the ACGT bases (in that order) across the sequence.
  Arg [-threshold]  : Scalar (optional) - Numeric minimum relative affinity for binding sites of this matrix
  Arg [-description]: Scalar (optional) - Descriptiom of matrix
  Example    : my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
                                                               -name  => "MA0122.1",
                                                               -analysis => $analysis,
                                                               -description => "Jaspar Matrix",
                                                                );
  Description: Constructor method for BindingMatrix class
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : Throws if name or/and type not defined
  Caller     : General
  Status     : Medium risk

=cut


# What we want to do is allow a 2d matrix array to be passed or the frequencies string
# Each would populate the other dynamically, no need to calc weights until they are required
# as this will slow down construction i.e. track display. 

sub new {
  my $caller    = shift;
  my $obj_class = ref($caller) || $caller;
  my $self      = $obj_class->SUPER::new(@_);
  
  my ( $name, $analysis, $freqs, $desc, $ftype, $thresh ) = rearrange
   (['NAME', 'ANALYSIS', 'FREQUENCIES', 'DESCRIPTION', 'FEATURE_TYPE', 'THRESHOLD'], @_);
  
  throw('Must supply a -name parameter')        if ! defined $name;
  throw('Must supply a -frequencies parameter') if ! defined $freqs;
  assert_ref($analysis, 'Bio::EnsEMBL::Analysis', 'Analysis');
  assert_ref($ftype, 'Bio::EnsEMBL::Funcgen::FeatureType', 'FeatureType');

  $self->{name}         = $name;
  $self->{analysis}     = $analysis;
  $self->{feature_type} = $ftype;

  if(check_ref($freqs, 'ARRAY')){
    $self->{freq_matrix} = $freqs;
    map { check_ref($_, 'ARRAY'); } @$freqs; # 2d array check
    $self->_validate_matrix;  # Also sets length
  }
  else{  # Assume string
    # defer validation to speed up track display
    $self->{tmp_frequencies}  = $freqs;    
  }

  $self->{description} = $desc if defined $desc;
  $self->{threshold}   = $thresh if defined $thresh;
  return $self;
}


=head2 feature_type

  Example    : my $ft_name = $matrix->feature_type->name;
  Description: Getter for the feature_type attribute for this matrix.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_type { return shift->{feature_type}; }


=head2 name

  Example    : my $name = $matrix->name();
  Description: Getter for the name attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return shift->{name}; }


=head2 description

  Example    : my $desc = $matrix->description;
  Description: Getter for the description attribute 
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub description { return shift->{description}; }


=head2 threshold

  Arg [1]    : Scalar (optional) - Numeric threshold
  Example    : if($score >= $matrix->threshold) { # Do something here }
  Description: Getter/setter for threshold attribute 
  Returntype : Scalar - numeric
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

# Do we ever really want to set this after construction? 
# thresholds are updated via direct sql, so unlikely unless
# someone if loading these via another path

sub threshold {
  my $self = shift;
  $self->{threshold} = shift if @_;
  return $self->{threshold};
}


=head2 analysis
  Example    : $matrix->analysis()->logic_name();
  Description: Getter for the feature_type attribute for this matrix.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub analysis { return shift->{analysis}; }


=head2 frequencies

  Arg[1]     : Boolean - Reverse complement flag
  Example    : print '>'.$matrix->name."\n".$matrix->frequencies;
  Description: Getter for the frequencies attribute. Returns 4 lines 
               of ACGT base frequencies for each position in the matrix.
  Returntype : Scalar - string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub frequencies {
  my $self      = shift;
  my $revcomp   = shift;
  my $attr_name = $revcomp ? 'rc_frequencies' : 'frequencies';

  if(! defined $self->{$attr_name}){
    $self->_build_matrix if defined $self->{tmp_frequencies};
    $self->{$attr_name} = join("\n", map { join(' ', @{$_} )} @{$self->matrix($revcomp)});
  }
 
  return $self->{$attr_name}; 
}


sub _build_matrix{
  my $self = shift;
  my @tmp  = split("\n", $self->{tmp_frequencies});
  shift @tmp if $tmp[0] =~ /^>/;
  @{$self->{freq_matrix}} = map { [ parse_matrix_line($_) ] } @tmp;
  $self->_validate_matrix;
  delete $self->{tmp_frequencies};
  return;
}


sub _validate_matrix{
  my $self = shift;

  if(scalar(@{$self->{freq_matrix}}) != 4){
    throw('Matrix '.$self->name." has an invalid number of lines:\n\t".
      join("\n\t", @{$self->{freq_matrix}}));
  }

  my $length = scalar @{$self->{freq_matrix}->[0]};
  $self->{length} = $length;

  for my $i(1..3){
    if($length != scalar(@{$self->{freq_matrix}->[$i]})){
      throw("Matrix has lines with differing sizes:\n\t".
        join("\n\t", @{$self->{freq_matrix}}));
    }
  }

  return;
}



sub frequency_matrix{ return shift->matrix(undef, shift); }

sub weight_matrix{ return shift->matrix(1, shift); }

=head2 matrix

  Arg[1]     : Boolean - Returns weights, default is frequencies
  Arg[2]     : Boolean - Returns reverse complement version of matrix.
  Example    : 
  Description: Getter for the frequency or weight matrix. 
               Array of ACGT base rows in that e.g.
                 [[a1, a2, a3, ...], # A values
                  [c1, c2, c3, ...], # C values
                  [g1, g2, g3, ...], # G values
                  [t1, t2, t3, ...]] # T values
  Returntype : Arrayref - 2d array
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

# Change this to also get weight matrix?

sub matrix{
  my $self      = shift;
  my $weights   = shift;
  my $revcomp   = shift;
  my $attr_name = $weights ? 'weight_matrix' : 'freq_matrix';
  
  if(! defined $self->{$attr_name}){

    if(! $weights){
      if(defined $self->{tmp_frequencies}){
        $self->_build_matrix;  
      }
    }
    elsif(! defined $self->{weights}){
      $self->_process_frequency_matrix;
    }
  }

  if($revcomp){
    $self->{$attr_name.'_rc'} = reverse_complement_matrix($self->{$attr_name});
    $attr_name .= '_rc' ;
  }

  return $self->{$attr_name};
}


=head2 frequencies_revcomp
  Example    : $matrix->frequencies_revcomp();
  Description: Getter for the reverse complement frequencies attribute.
               The attribute represents the reverse complement of frequencies
  Returntype : Scalar - string
  Caller     : General
  Status     : At Risk
=cut
# CThis was not revcomp, just rev!!!
# Remove this in favour of passing a revcomp arg to frequencies
#sub frequencies_revcomp { return shift->{frequencies_revcomp}; }


=head2 relative_affinity

  Arg [1]    : String - Binding site sequence (strand matched to the matrix alignment)
  Arg [2]    : Boolean (optional) - Linear scale results (default is log scale)
  Example    : $matrix->relative_affinity($sequence);
  Description: Calculates the binding affinity of a given sequence relative to
               the optimal site for the matrix. The site is taken as if it were 
               in the proper orientation. Considers a purely random background e.g.
                   p(A)=p(C)=p(G)=p(T)
               Returns value between 0 and 1    
  Returntype : Scalar - numeric or undef
  Exceptions : Throws if the sequence length does not have the matrix length
               or if the sequence has unclear bases (N is not accepted)
  Caller     : General
  Status     : At Risk

=cut

sub relative_affinity {
  my $self    = shift;
  my $seq     = shift;
  my $linear  = shift;
  throw('Must provide a sequence argument') if ! defined $seq;
  ($seq = uc($seq)) =~ s/\s+//g; #be forgiving of case and spaces
  
  if($seq =~ /[^ACGT]/){
    warn($self->name." sequence contains in-valid [^acgtACGT] characters:\t$seq");
    return; # Not undef which can be true in list context
  }
  
  if(length($seq) != $self->length){
    throw('Specified sequence does not match matrix length('.
      $self->length."):\n".length($seq)." $seq");
  }
  
  my $rel_aff;
  my $log_odds = 0;
  my @bases         = split(//, $seq);
  my $weight_matrix = $self->weights;

  for my $i(0..$#bases){
    $log_odds += $weight_matrix->{$bases[$i]}->[$i];  
  }
  
  #This log scale may be quite unrealistic, but useful just for comparison.
  if(! $linear){
    $rel_aff = ($log_odds - $self->min_affinity) / 
                 ($self->max_affinity - $self->min_affinity);
  } else {
    $rel_aff = (exp($log_odds) - exp($self->min_affinity)) / 
                 (exp($self->max_affinity) - exp($self->min_affinity));
  }

  return $rel_aff;
}


=head2 is_position_informative

  Arg [1]    : Scalar - 1 based integer position within the matrix
  Arg [2]    : Scalar (optional) - Threshold [0-2] for information content [default is 1.5]
  Example    : $matrix->is_position_informative($pos);
  Description: Returns true if position information content is over threshold
  Returntype : Boolean
  Exceptions : Throws if position or threshold out of bounds
  Caller     : General
  Status     : At High Risk

=cut

sub is_position_informative {
  my ($self, $position, $revcomp, $threshold) = @_;
  throw('Must pass a position argument') if ! defined $position;

  if($revcomp){ 
    $position = $self->length - $position + 1; 
  }

  if(($position < 1) || ($position > $self->length)){
    throw("Position($position) is out of bounds: 1 - ".$self->length);
  }

  if(! defined $threshold){
    $threshold = 1.5;
  }
  elsif(($threshold < 0) || ($threshold > 2)){
    throw("Threshold($threshold) must be within 0 - 2") 
  }

  return ($self->info_content->[$position-1] >= $threshold);
}


=head2 length

  Example    : my $length = $bm->length;
  Description: Returns the length of the the matrix
  Returntype : Scalar - numeric (int)
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub length { 
  my $self = shift;
  # Length is set by _validate_matrix if matrix is passed to constructor
  # So if it is not defined need to _build_matrix from frequencies
  $self->_build_matrix if ! defined $self->{length};
  return $self->{length}; 
}


=head2 weights

  Example    : 
  Description: Private getter for the weights matrix based on frequencies
  Returntype : HASHREF with the weights of this binding matrix 
  Exceptions : None
  Caller     : General + relative_affinity
  Status     : At Risk

=cut

sub weights {
  my $self    = shift;
  $self->_process_frequency_matrix if ! defined $self->{weights};
  return $self->{weights};  
}


=head2 _process_frequency_matrix

  Arg[1]     : Scalar (optional) - Pseudo count/lapace estimator for base weight 
                 computation. This allows computation for bases with frequencies 
                 of 0. Default is 0.1.
  Example    : $self->_process_frequency_matrix;
  Description: Private function to calculate the matrix information content, 
                 base weights and  min and max relative affinity.
  ReturnType : Arrayref
  Caller     : weights and information_content
  Status     : At Risk

=cut

# We can allow distinct background per nucleotide, instead of 0.25 for all... pass as parameter?
# This should probably be a meta entry, loaded by the adaptor, or defaults to 0.25 for all.
# But if the matrix was obtained using in-vivo data, it shouldn't matter the organism nucleotide bias.
# We're using 0.1 as pseudo-count... the matrix cannot have very few elements... (e.g. <30 not good)

# Is 0.1 optimal here?
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/pdf/gkn1019.pdf

sub _process_frequency_matrix{
  my $self   = shift;
  my $pseudo = shift;
  $pseudo ||= 0.1;

  my ($As, $Cs, $Gs, $Ts) = @{$self->matrix}; 
  my ($tmp_min, $tmp_max, $total_info);
  $self->{weights} = {A => [], C => [], G => [], T => []};  # This gets populated in _compute_base_weight
  $self->{info_content} = [];
  $self->{min_affinity} = 0; 
  $self->{max_affinity} = 0;

  for my $i(0..$#{$As}){ 
    $total_info = $As->[$i] + $Cs->[$i] + $Gs->[$i] + $Ts->[$i] + (4 * $pseudo);
    ($tmp_min, undef, undef, $tmp_max) = 
      sort {$a <=> $b} ($self->_compute_base_weight('A', $As->[$i], $total_info, $pseudo),
                        $self->_compute_base_weight('C', $Cs->[$i], $total_info, $pseudo),
                        $self->_compute_base_weight('G', $Gs->[$i], $total_info, $pseudo),
                        $self->_compute_base_weight('T', $Ts->[$i], $total_info, $pseudo));
    $self->{min_affinity} += $tmp_min;
    $self->{max_affinity} += $tmp_max; 

    # Slight redundancy here with calculating $base_freq + $pseudo

    my $fas = ($As->[$i] + $pseudo) / $total_info;
    my $fcs = ($Cs->[$i] + $pseudo) / $total_info;
    my $fgs = ($Gs->[$i] + $pseudo) / $total_info;
    my $fts = ($Ts->[$i] + $pseudo) / $total_info;    
    my $ic_i = 2 + ($fas * log($fas) / log(2)) + ($fcs * log($fcs) / log(2)) + 
      ($fgs * log($fgs) / log(2)) + ($fts * log($fts)/log(2));
    push @{$self->{info_content}}, $ic_i;
  }

  $self->{weight_matrix} = [$self->{weights}{A},
                            $self->{weights}{C},
                            $self->{weights}{G},
                            $self->{weights}{T}];
  return;
}


sub _compute_base_weight{
  my ($self, $base, $base_freq, $total_info, $pseudo) = @_;
  my $weight = log(( ($base_freq + $pseudo) / ($total_info + (4 * $pseudo) )) / 0.25); 
  push @{$self->{weights}{$base}}, $weight; 
  return $weight;
}


=head2 info_content

  Example    : info_content($as,$cs,$gs,$ts,$pseudo);
  Description: Private function to calculate the matrix information content per position
  ReturnType : Arrayref
  Caller     : self
  Status     : At Risk

=cut

sub info_content {
  my $self    = shift;
  my $revcomp = shift;
  $self->_process_frequency_matrix if ! defined $self->{info_content};
  return $revcomp ? [ reverse(@{$self->{info_content}}) ] : $self->{info_content};
}


=head2 max_affinity

  Example    : my $max = $matrix->maximim_affinity;
  Description: Getter for maximum binding affinity attribute
  Returntype : Scalar - numeric 
  Exceptions : None
  Caller     : relative_affinity
  Status     : At Risk

=cut

sub max_affinity {  
  my $self = shift;
  $self->_process_frequency_matrix if ! defined $self->{max_affinity};
  return $self->{max_affinity}; 
}


=head2 _min_affinity

  Example    : my $min = $matrix->min_affinity;
  Description: Getter of minimum binding affinity attribute
  Returntype : Scalar - numeric
  Exceptions : None
  Caller     : relative_affinity
  Status     : At Risk

=cut

sub min_affinity { 
  my $self = shift;
  $self->_process_frequency_matrix if ! defined $self->{min_affinity};
  return $self->{min_affinity}; 
}

1;
