=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::APPRIS - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::APPRIS->new(
    -result                          => <Anlysis result of the gene/transcript>,
    
    -functional_residues_score        => <Num. functional residues>,
    -homologous_structure_score       => <Score of the homologous structure>,
    -vertebrate_conservation_score    => <Score of the vertebrate conservation>,
    -domain_score                     => <the number of domains:
                                            'Num. whole domains' + 
                                            'Num. possible damaged domains' + 
                                            'Num. damaged domains' + 
                                            'Num. wrong domains'>,
    -transmembrane_helices_score      => <Number of transmembre helices:
                                             'Num. transmembrane helices' +
                                             'Num. damaged helices'>,
    -peptide_score                    => <Peptide signals sequence>,
    -mitochondrial_score              => <Mitochondrial signals sequence>,
    -unusual_evolution_score          => <number of unusual exons:
                                            'Num. unusual exons (Consensus)' =
                                            'Num. unusual exons (MAF)' + 
                                            'Num. unusual exons (Prank)' + 
                                            'Num. unusual exons (Kaling)'>,
    -peptide_evidence_score          => <Num. peptides>,
    -principal_isoform_score         => <Score analysed>,    
    
    -functional_residues_signal      => <Annotation analysed>,
    -homologous_structure_signal     => <Annotation analysed>,
    -vertebrate_conservation_signal  => <Annotation analysed>,
    -domain_signal                   => <Annotation analysed>,
    -transmembrane_helices_signal    => <Annotation analysed>,
    -peptide_signal                  => <Annotation analysed>,
    -mitochondrial_signal            => <Annotation analysed>,    
    -unusual_evolution_signal        => <Annotation analysed>,
    -peptide_evidence_signal         => <Annotation analysed>,
    -principal_isoform_signal        => <Annotation analysed>,
    
    -reliability                      => <Reliability of Annotation analysed>,
        
  );

=head1 DESCRIPTION

A representation of analysis of APPRIS within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::APPRIS;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new
  
  Arg [-result]  : 
       string - the anlysis result of the gene/transcript
  Arg [-functional_residues_score]:
        string - num. functional residues
  Arg [-homologous_structure_score]:
        string - score of the homologous structure
  Arg [-vertebrate_conservation_score]:
        string - the number of domains
  Arg [-score_con_vertebrate]:
        string - score of the vertebrate conservation
  Arg [-transmembrane_helices_score]:
        string - number of transmembre helices
  Arg [-peptide_score]:
        string - score of peptide signal sequence        
  Arg [-mitochondrial_score]:
        string - score of mitochondrial signal sequence
  Arg [-unusual_evolution_score]:
        string - number of unusual exons
  Arg [-peptide_evidence_score]:
        string - num. peptides 
  Arg [-principal_isoform_score]:
        string - the appris score for this analysis
  Arg [-functional_residues_signal]:
        string - the appris annotation for this analysis
  Arg [-homologous_structure_signal]:
        string - the appris annotation for this analysis
  Arg [-vertebrate_conservation_signal]:
        string - the appris annotation for this analysis
  Arg [-domain_signal]:
        string - the appris annotation for this analysis
  Arg [-transmembrane_helices_signal]:
        string - the appris annotation for this analysis
  Arg [-peptide_signal]:
        string - the appris annotation for this analysis
  Arg [-mitochondrial_signal]:
        string - the appris annotation for this analysis
  Arg [-unusual_evolution_signal]:
        string - the appris annotation for this analysis
  Arg [-peptide_evidence_signal]:
        string - the appris annotation for this analysis
  Arg [-principal_isoform_signal]:
        string - the appris annotation for this analysis
  Arg [-reliability]: (optional)
        string - the reliability of the final analysis
        
  Example    : $analysis = APPRIS::Analysis::APPRIS->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::APPRIS
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
	my ($caller) = shift;

	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;
	
	my (
		$result,							
		
		$functional_residues_score,			$homologous_structure_score,
		$vertebrate_conservation_score,		$domain_score,
		$transmembrane_helices_score,		$peptide_score,
		$mitochondrial_score,				$unusual_evolution_score,
		$peptide_evidence_score,			$principal_isoform_score,
		
		$functional_residues_signal,		$homologous_structure_signal,
		$vertebrate_conservation_signal,	$domain_signal,
		$transmembrane_helices_signal,		$peptide_signal,
		$mitochondrial_signal,				$unusual_evolution_signal,
		$peptide_evidence_signal,			$principal_isoform_signal,				
		$reliability,

	)
	= rearrange( [
		'result',							
		
		'functional_residues_score',		'homologous_structure_score',		
		'vertebrate_conservation_score',	'domain_score',		
		'transmembrane_helices_score',		'peptide_score',
		'mitochondrial_score',				'unusual_evolution_score',
		'peptide_evidence_score',			'principal_isoform_score',
		
		'functional_residues_signal',		'homologous_structure_signal',		
		'vertebrate_conservation_signal',	'domain_signal',		
		'transmembrane_helices_signal',		'peptide_signal',
		'mitochondrial_signal',				'unusual_evolution_signal',			
		'peptide_evidence_signal',			'principal_isoform_signal',
		'reliability',
	],
	@_
	);
	
	$self->result($result);
	
 	$self->functional_residues_score($functional_residues_score);
 	$self->homologous_structure_score($homologous_structure_score);
 	$self->vertebrate_conservation_score($vertebrate_conservation_score);
 	$self->domain_score($domain_score);
 	$self->transmembrane_helices_score($transmembrane_helices_score);
 	$self->peptide_score($peptide_score);
 	$self->mitochondrial_score($mitochondrial_score);
 	$self->unusual_evolution_score($unusual_evolution_score);
 	$self->peptide_evidence_score($peptide_evidence_score);
 	$self->principal_isoform_score($principal_isoform_score);
 	
 	$self->functional_residues_signal($functional_residues_signal);
 	$self->homologous_structure_signal($homologous_structure_signal);
 	$self->vertebrate_conservation_signal($vertebrate_conservation_signal);
 	$self->domain_signal($domain_signal);
 	$self->transmembrane_helices_signal($transmembrane_helices_signal);
 	$self->peptide_signal($peptide_signal);
 	$self->mitochondrial_signal($mitochondrial_signal);
 	$self->unusual_evolution_signal($unusual_evolution_signal);
 	$self->peptide_evidence_signal($peptide_evidence_signal);
 	$self->principal_isoform_signal($principal_isoform_signal);
 	$self->reliability($reliability); 	 	
		
	return $self;
}


=head2 result

  Arg [1]    : (optional) String - the result of the analysis
  Example    : $analysis->result('result');
  Description: Getter/setter for the result of analysis 
               that for this analysis
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub result {
	my ($self) = shift;
	$self->{'result'} = shift if(@_);
	return $self->{'result'};
}

=head2 functional_residues_score

  Arg [1]    : (optional) String - the num. functional residues
  Example    : $analysis->functional_residues_score($value);
  Description: Getter/setter for the num. functional residues
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub functional_residues_score {
	my ($self) = shift;
	$self->{'functional_residues_score'} = shift if(@_);
	return $self->{'functional_residues_score'};
}

=head2 homologous_structure_score

  Arg [1]    : (optional) String - the score of the homologous structure
  Example    : $analysis->homologous_structure_score($value);
  Description: Getter/setter for the score of the homologous structure
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub homologous_structure_score {
	my ($self) = shift;
	$self->{'homologous_structure_score'} = shift if(@_);
	return $self->{'homologous_structure_score'};
}

=head2 vertebrate_conservation_score

  Arg [1]    : (optional) String - the score of the vertebrate conservation
  Example    : $analysis->vertebrate_conservation_score($value);
  Description: Getter/setter for the number of domains
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub vertebrate_conservation_score {
	my ($self) = shift;
	$self->{'vertebrate_conservation_score'} = shift if(@_);
	return $self->{'vertebrate_conservation_score'};
}

=head2 domain_score

  Arg [1]    : (optional) String - the number of domains
  Example    : $analysis->domain_score($value);
  Description: Getter/setter for the number of domains
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub domain_score {
	my ($self) = shift;
	$self->{'domain_score'} = shift if(@_);
	return $self->{'domain_score'};
}

=head2 transmembrane_helices_score

  Arg [1]    : (optional) String - the number of transmembre helices
  Example    : $analysis->transmembrane_helices_score($value);
  Description: Getter/setter for the number of transmembre helices
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transmembrane_helices_score {
	my ($self) = shift;
	$self->{'transmembrane_helices_score'} = shift if(@_);
	return $self->{'transmembrane_helices_score'};
}

=head2 peptide_score

  Arg [1]    : (optional) String - the score of peptide signal sequence
  Example    : $analysis->peptide_score($value);
  Description: Getter/setter for the score of peptide signal sequence
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptide_score {
	my ($self) = shift;
	$self->{'peptide_score'} = shift if(@_);
	return $self->{'peptide_score'};
}

=head2 mitochondrial_score

  Arg [1]    : (optional) String - the score of mitochondrial signal sequence
  Example    : $analysis->mitochondrial_score($value);
  Description: Getter/setter for the score of mitochondrial signal sequence
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub mitochondrial_score {
	my ($self) = shift;
	$self->{'mitochondrial_score'} = shift if(@_);
	return $self->{'mitochondrial_score'};
}

=head2 unusual_evolution_score

  Arg [1]    : (optional) String - the number of unusual exons
  Example    : $analysis->unusual_evolution_score($value);
  Description: Getter/setter for the number of unusual exons
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub unusual_evolution_score {
	my ($self) = shift;
	$self->{'unusual_evolution_score'} = shift if(@_);
	return $self->{'unusual_evolution_score'};
}

=head2 peptide_evidence_score

  Arg [1]    : (optional) String - the score of mitochondrial signal sequence
  Example    : $analysis->peptide_evidence_score($value);
  Description: Getter/setter for the score of mitochondrial signal sequence
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptide_evidence_score {
	my ($self) = shift;
	$self->{'peptide_evidence_score'} = shift if(@_);
	return $self->{'peptide_evidence_score'};
}

=head2 principal_isoform_score

  Arg [1]    : (optional) String - the analysed principal_isoform_score to set
  Example    : $analysis->principal_isoform_score($method);
  Description: Getter/setter for the analysed principal_isoform_score
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub principal_isoform_score {
	my ($self) = shift;
	$self->{'principal_isoform_score'} = shift if(@_);
	return $self->{'principal_isoform_score'};
}

=head2 functional_residues_signal

  Arg [1]    : (optional) String - the analysed functional_residues_signal to set
  Example    : $analysis->functional_residues_signal($method);
  Description: Getter/setter for the analysed functional_residues_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub functional_residues_signal {
	my ($self) = shift;
	$self->{'functional_residues_signal'} = shift if(@_);
	return $self->{'functional_residues_signal'};
}

=head2 homologous_structure_signal

  Arg [1]    : (optional) String - the analysed homologous_structure_signal to set
  Example    : $analysis->homologous_structure_signal($method);
  Description: Getter/setter for the analysed homologous_structure_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub homologous_structure_signal {
	my ($self) = shift;
	$self->{'homologous_structure_signal'} = shift if(@_);
	return $self->{'homologous_structure_signal'};
}

=head2 vertebrate_conservation_signal

  Arg [1]    : (optional) String - the analysed vertebrate_conservation_signal to set
  Example    : $analysis->vertebrate_conservation_signal($method);
  Description: Getter/setter for the analysed vertebrate_conservation_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub vertebrate_conservation_signal {
	my ($self) = shift;
	$self->{'vertebrate_conservation_signal'} = shift if(@_);
	return $self->{'vertebrate_conservation_signal'};
}

=head2 domain_signal

  Arg [1]    : (optional) String - the analysed domain_signal to set
  Example    : $analysis->domain_signal($method);
  Description: Getter/setter for the analysed domain_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub domain_signal {
	my ($self) = shift;
	$self->{'domain_signal'} = shift if(@_);
	return $self->{'domain_signal'};
}

=head2 transmembrane_helices_signal

  Arg [1]    : (optional) String - the analysed transmembrane_helices_signal to set
  Example    : $analysis->transmembrane_helices_signal($method);
  Description: Getter/setter for the analysed transmembrane_helices_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transmembrane_helices_signal {
	my ($self) = shift;
	$self->{'transmembrane_helices_signal'} = shift if(@_);
	return $self->{'transmembrane_helices_signal'};
}

=head2 peptide_signal

  Arg [1]    : (optional) String - the analysed peptide_signal to set
  Example    : $analysis->peptide_signal($method);
  Description: Getter/setter for the analysed peptide_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptide_signal {
	my ($self) = shift;
	$self->{'peptide_signal'} = shift if(@_);
	return $self->{'peptide_signal'};
}

=head2 mitochondrial_signal

  Arg [1]    : (optional) String - the analysed mitochondrial_signal to set
  Example    : $analysis->mitochondrial_signal($method);
  Description: Getter/setter for the analysed mitochondrial_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub mitochondrial_signal {
	my ($self) = shift;
	$self->{'mitochondrial_signal'} = shift if(@_);
	return $self->{'mitochondrial_signal'};
}

=head2 unusual_evolution_signal

  Arg [1]    : (optional) String - the analysed unusual_evolution_signal to set
  Example    : $analysis->unusual_evolution_signal($method);
  Description: Getter/setter for the analysed unusual_evolution_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub unusual_evolution_signal {
	my ($self) = shift;
	$self->{'unusual_evolution_signal'} = shift if(@_);
	return $self->{'unusual_evolution_signal'};
}

=head2 peptide_evidence_signal

  Arg [1]    : (optional) String - the analysed peptide_evidence_signal to set
  Example    : $analysis->peptide_evidence_signal($method);
  Description: Getter/setter for the analysed peptide_evidence_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptide_evidence_signal {
	my ($self) = shift;
	$self->{'peptide_evidence_signal'} = shift if(@_);
	return $self->{'peptide_evidence_signal'};
}

=head2 principal_isoform_signal

  Arg [1]    : (optional) String - the analysed principal_isoform_signal to set
  Example    : $analysis->principal_isoform_signal($method);
  Description: Getter/setter for the analysed principal_isoform_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub principal_isoform_signal {
	my ($self) = shift;
	$self->{'principal_isoform_signal'} = shift if(@_);
	return $self->{'principal_isoform_signal'};
}

=head2 reliability

  Arg [1]    : (optional) String - the reliability of the analysis
  Example    : $analysis->reliability('100');
  Description: Getter/setter for the reliability of analysis 
               that for this analysis
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub reliability {
	my ($self) = shift;
	$self->{'reliability'} = shift if(@_);
	return $self->{'reliability'};
}

sub DESTROY {}

1;
