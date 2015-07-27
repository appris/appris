#
# Ensembl module for Bio::EnsEMBL::Funcgen::BindingMatrix
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

Bio::EnsEMBL::Funcgen::BindingMatrix - A module to represent a BindingMatrix. 
In EFG this represents the binding affinities of a Transcription Factor to DNA.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::BindingMatrix;

my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
                                                        -name  => "MA0122.1",
                                                        -type => "Jaspar",
                                                        -description => "Nkx3-2 Jaspar Matrix",
                                                       );
$matrix->frequencies("A  [ 4  1 13 24  0  0  6  4  9 ]
                      C  [ 7  4  1  0  0  0  0  6  7 ]
                      G  [ 4  5  7  0 24  0 18 12  5 ]
                      T  [ 9 14  3  0  0 24  0  2  3 ]");

print $matrix->relative_affinity("TGGCCACCA")."\n";

print $matrix->threshold."\n";

=head1 DESCRIPTION

This class represents information about a BindingMatrix, containing the name 
(e.g. the Jaspar ID, or an internal name), and description. A BindingMatrix 
is always associated to an Analysis (indicating the origin of the matrix e.g. 
Jaspar) and a FeatureType (the binding factor).   

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor

=cut

package Bio::EnsEMBL::Funcgen::BindingMatrix;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw );
  
use base qw(Bio::EnsEMBL::Funcgen::Storable);

=head2 new

  Arg [-name]: string - name of Matrix
  Arg [-analysis]: Bio::EnsEMBL::Analysis - analysis describing how the matrix was obtained
  Arg [-frequencies]: (optional) string - frequencies representing the binding affinities of a Matrix
  Arg [-threshold]: (optional) float - minimum relative affinity for binding sites of this matrix
  Arg [-description]: (optional) string - descriptiom of Matrix
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

sub new {
  my $caller = shift;

  my $obj_class = ref($caller) || $caller;
  my $self = $obj_class->SUPER::new(@_);
  
  my ( $name, $analysis, $freq, $desc, $ftype, $thresh ) = rearrange
	( [
	   'NAME', 'ANALYSIS', 'FREQUENCIES', 'DESCRIPTION', 'FEATURE_TYPE', 'THRESHOLD'
	  ], @_);
  
  
  if(! defined $name){
    throw("Must supply a name\n");
  }

  if(! ((ref $analysis) && $analysis->isa('Bio::EnsEMBL::Analysis') )){
    throw("You must define a valid Bio::EnsEMBL::Analysis");
    #leave is stored test to adaptor
  }
  
  if(! (ref($ftype) && $ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType'))){
	throw("You must define a valid Bio::EnsEMBL::Funcgen::FeatureType");
	#leave is stored test to adaptor
  }

  $self->name($name);
  $self->{analysis} = $analysis;
  $self->{feature_type} = $ftype;
  $self->frequencies($freq) if $freq;
  $self->description($desc) if $desc;
  $self->threshold($thresh) if $thresh;

  return $self;
}


=head2 feature_type

  Example    : my $ft_name = $matrix->feature_type()->name();
  Description: Getter for the feature_type attribute for this matrix.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub feature_type {
  my $self = shift;
  
  return $self->{'feature_type'};
}


=head2 name

  Arg [1]    : (optional) string - name
  Example    : my $name = $matrix->name();
  Description: Getter and setter of name attribute
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub name {
    my $self = shift;
    $self->{'name'} = shift if @_;
    return $self->{'name'};
}

=head2 description

  Arg [1]    : (optional) string - description
  Example    : my $desc = $matrix->description();
  Description: Getter and setter of description attribute 
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if @_;
    return $self->{'description'};
}

=head2 threshold

  Arg [1]    : (optional) float - threshold
  Example    : my $thresh = $matrix->threshold();
  Description: Getter and setter of threshold attribute 
  Returntype : float
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub threshold {
    my $self = shift;
    $self->{'threshold'} = shift if @_;
    return $self->{'threshold'};
}


=head2 analysis
  Example    : $matrix->analysis()->logic_name();
  Description: Getter for the feature_type attribute for this matrix.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub analysis {
  my $self = shift;

  return $self->{'analysis'};
}


=head2 frequencies

  Arg [1]    : (optional) string - frequencies
  Example    : $matrix->frequencies($frequencies_string);
  Description: Getter and setter of frequencies attribute

  The attribute is a string representing the matrix binding
  affinities in the Jaspar format. E.g. 
  	">
  	[ ]
  	"
  
  Returntype : string
  Exceptions : Throws if the string attribute is not a properly
  formed matrix in the Jaspar format
  Caller     : General
  Status     : At Risk

=cut

sub frequencies {
  my $self = shift;
 
  my $frequencies = shift if @_; 
  if($frequencies){
  	$self->_weights($frequencies);
  	$self->{'frequencies'} = $frequencies;  	
  }
  return $self->{'frequencies'};
}

=head2 frequencies_revcomp

  Example    : $matrix->frequencies_revcomp();
  Description: Getter for the reverse complement frequencies attribute

  The attribute represents the reverse complement of frequencies
  
  Returntype : string
  Caller     : General
  Status     : At Risk

=cut

sub frequencies_revcomp {
  my $self = shift;
 
  return $self->{'frequencies_revcomp'};
}


=head2 relative_affinity

  Arg [1]    : string - Binding Site Sequence
  Arg [2]    : (optional) boolean - 1 if results are to be in linear scale (default is log scale)
  Example    : $matrix->relative_affinity($sequence);
  Description: Calculates the binding affinity of a given sequence
	relative to the optimal site for the matrix
	The site is taken as if it were in the proper orientation
        Considers a purely random background p(A)=p(C)=p(G)=p(T)
  Returntype : double
  Exceptions : Throws if the sequence length does not have the matrix length
  or if the sequence has unclear bases (N is not accepted)
  Caller     : General
  Status     : At Risk

=cut

sub relative_affinity {
  my ($self, $sequence, $linear) = (shift, shift, shift);
  $sequence =~ s/^\s+//;
  $sequence =~ s/\s+$//;
  
  throw "No sequence given" if !$sequence;
  $sequence = uc($sequence);
  if($sequence =~ /[^ACGT]/){
    throw "Sequence $sequence contains invalid characters: Only Aa Cc Gg Tt accepted";	
  }
  
  my $weight_matrix = $self->_weights;
  my $matrix_length = scalar(@{$weight_matrix->{'A'}});
  if(length($sequence) != $matrix_length){
    throw "Sequence $sequence does not have length $matrix_length";
  }
  
  my $log_odds = 0;
  my @bases = split(//,$sequence);
  for(my $i=0;$i<$matrix_length;$i++){
    $log_odds += $weight_matrix->{$bases[$i]}->[$i];	
  }
  
  #This log scale may be quite unrealistic... but usefull just for comparisons...
  if(!$linear){
    return ($log_odds - $self->_min_bind) / ($self->_max_bind - $self->_min_bind);
  } else {
    return (exp($log_odds) - exp($self->_min_bind)) / (exp($self->_max_bind) - exp($self->_min_bind));
  }

}

=head2 is_position_informative

  Arg [1]    : Int - 1bp position within the matrix
  Arg [2]    : (optional) double - threshold [0-2] for information content [default is 1.5]
  Example    : $matrix->is_position_informative($pos);
  Description: Returns true if position information content is over threshold
  Returntype : boolean
  Exceptions : Throws if position or threshold out of bounds
  Caller     : General
  Status     : At High Risk

=cut

sub is_position_informative {
  my ($self, $position, $threshold) = @_;
  
  throw "Need a position" if(!defined($position));
  throw "Position out of bounds" if(($position<1) || ($position > $self->length));
  if(!defined($threshold)){ $threshold = 1.5; }
  throw "Threshold out of bounds" if(($threshold<0) || ($threshold>2));
  return ($self->{'ic'}->[$position-1] >= $threshold);
}



=head2 length

  Example    : $bm->length();
  Description: Returns the length of the the matrix (e.g. 19bp long)
  Returntype : int with the length of this binding matrix 
  Exceptions : none
  Caller     : General
  Status     : At Risk

=cut

sub length {
  my $self = shift;

  my $weight_matrix = $self->_weights;

  return scalar(@{$weight_matrix->{'A'}});
}

=head2 _weights

  Arg [1]    : (optional) string - frequencies
  Example    : _weights($frequencies);
  Description: Private Getter Setter for the weight matrix based on frequencies
  Returntype : HASHREF with the weights of this binding matrix 
  Exceptions : Throws if the frequencies attribute string does not correspond 
  to 4 rows of an equal number of integer numbers.
  Caller     : Self
  Status     : At Risk

=cut

sub _weights {
	my $self = shift;
	
 	#for the moment use equiprobability and constant pseudo-count
	my $pseudo = 0.1;

 	#TODO allow for it to be passed as parameters?
  	my $frequencies = shift if @_; 
  	if($frequencies){
  		$frequencies =~ s/^(>.*?\n)//;
		my $header = $1;

  		my ($a,$c,$g,$t) = split(/\n/,$frequencies);
  		my @As = split(/\s+/,_parse_matrix_line('[A\[\]]',$a));
  		my @Cs = split(/\s+/,_parse_matrix_line('[C\[\]]',$c));
  		my @Gs = split(/\s+/,_parse_matrix_line('[G\[\]]',$g));
  		my @Ts = split(/\s+/,_parse_matrix_line('[T\[\]]',$t));
		if((scalar(@As)!=scalar(@Cs)) || (scalar(@As)!=scalar(@Gs)) || (scalar(@As)!=scalar(@Ts)) ){
			throw "Frequencies provided are not a valid frequency matrix"		
		}
		$self->_calc_ic(\@As,\@Cs,\@Gs,\@Ts,$pseudo);

		#Create the reverse complement
		my @revT = reverse(@As);
		my @revA = reverse(@Ts);
		my @revC = reverse(@Gs);
		my @revG = reverse(@Cs);
		my $revcomp = $header;
		$revcomp.= "A [ ".join("\t",@revA)." ]\n";
		$revcomp.= "C [ ".join("\t",@revC)." ]\n";
		$revcomp.= "G [ ".join("\t",@revG)." ]\n";
		$revcomp.= "T [ ".join("\t",@revT)." ]\n";
		$self->{'frequencies_revcomp'} = $revcomp;

  		my @totals;
  		for(my $i=0;$i<scalar(@As);$i++){ 
  			$totals[$i]=$As[$i]+$Cs[$i]+$Gs[$i]+$Ts[$i];
  		}
  		
		my %weights;			
		#We can allow distinct background per nucleotide, instead of 0.25 for all... pass as parameter
		#But if the matrix was obtained using in-vivo data, it shouldn't matter the organism nucleotide bias..
		#We're using 0.1 as pseudo-count... the matrix cannot have very few elements... (e.g. <30 not good)
		my @was; for(my $i=0;$i<scalar(@As);$i++){ $was[$i] = log((($As[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };
		$weights{'A'} = \@was;
		my @wcs; for(my $i=0;$i<scalar(@Cs);$i++){ $wcs[$i] = log((($Cs[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };
		$weights{'C'} = \@wcs;
		my @wgs; for(my $i=0;$i<scalar(@Gs);$i++){ $wgs[$i] = log((($Gs[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };
		$weights{'G'} = \@wgs;
		my @wts; for(my $i=0;$i<scalar(@Ts);$i++){ $wts[$i] = log((($Ts[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };	     
		$weights{'T'} = \@wts;
	
		$self->{'weights'} = \%weights;

		my $max = 0; my $min = 0;
		for(my $i=0;$i<scalar(@As);$i++){
			my $col = [ $was[$i], $wcs[$i], $wgs[$i], $wts[$i] ];
			$min += _min($col);
			$max += _max($col);
		}

		#Log scale
		$self->_max_bind($max);
		$self->_min_bind($min);
	}
	
	return $self->{'weights'};	

}

=head2 _calc_ic

  Example    : _calc_ic($as,$cs,$gs,$ts,$pseudo);
  Description: Private function to calculate the matrix information content per position
  Caller     : self
  Status     : At Risk

=cut

sub _calc_ic {
  my ($self,$as, $cs, $gs, $ts,$pseudo) = (shift,shift, shift, shift, shift, shift);
  my @ic = ();
  for (my $i=0;$i<scalar(@$as);$i++){
    my $total_i = $as->[$i] + $cs->[$i] + $gs->[$i] + $ts->[$i] + (4*$pseudo);
    my $fas = ($as->[$i] + $pseudo) / $total_i;
    my $fcs = ($cs->[$i] + $pseudo) / $total_i;
    my $fgs = ($gs->[$i] + $pseudo) / $total_i;
    my $fts = ($ts->[$i] + $pseudo) / $total_i;    
    my $ic_i = 2 + ($fas * log($fas)/log(2)) + ($fcs * log($fcs)/log(2)) + ($fgs * log($fgs)/log(2)) + ($fts * log($fts)/log(2));
    push @ic, $ic_i;
  }
  $self->{'ic'} = \@ic;
}

sub _parse_matrix_line {
	 my ($pat,$line) = (shift,shift);
	 $line=~s/$pat//g; $line=~s/^\s+//; $line=~s/\s+$//;	
	 return $line;
}

sub _max { return _min_max(shift, 0); }

sub _min { return _min_max(shift, 1); }

sub _min_max {
	my ($list,$min) = (shift, shift);
	my $min_max = $list->[0];
	map { if($min ? $_ < $min_max : $_ > $min_max){ $min_max = $_; } } @$list;
	return $min_max;
}


=head2 _max_bind

  Arg [1]    : (optional) double - maximum binding affinity
  Example    : $matrix->_max_bind(10.2);
  Description: Private Getter and setter of max_bind attribute (not to be called directly)
  Returntype : float with the maximum binding affinity of the matrix 
  Exceptions : None
  Caller     : Self
  Status     : At Risk

=cut

sub _max_bind {
  my $self = shift;
  
  $self->{'max_bind'} = shift if @_;

  return $self->{'max_bind'};
}

=head2 _min_bind

  Arg [1]    : (optional) double - minimum binding affinity
  Example    : $matrix->_min_bind(-10.2);
  Description: Private Getter and setter of min_bind attribute (not to be called directly)
  Returntype : float with the minimum binding affinity of the matrix 
  Exceptions : None
  Caller     : Self
  Status     : At Risk

=cut

sub _min_bind {
  my $self = shift;
  
  $self->{'min_bind'} = shift if @_;

  return $self->{'min_bind'};
}

1;
