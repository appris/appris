=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Exporter

=head1 SYNOPSIS

  $registry = APPRIS::Exporter->new(
    -dbhost  => 'localhost',
    -dbname  => 'homo_sapiens_encode_3c',
    -dbuser  => 'jmrodriguez'
    );

  $gene = $registry->fetch_by_stable_id($stable_id);

  @genes = @{ $registry->fetch_by_chr_start_end('X', 1, 10000) };

=head1 DESCRIPTION

All Adaptors are stored/registered using this module.
This module should then be used to get the adaptors needed.

The registry can be loaded from a configuration file or from database info.

=head1 METHODS

=cut

package APPRIS::Exporter;

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

use APPRIS::Export::GTF;
use APPRIS::Export::GFF3;
use APPRIS::Export::SEQ;
use APPRIS::Export::CDS;
use APPRIS::Export::BED;
use APPRIS::Export::BED12;
use APPRIS::Export::JSON;
use APPRIS::Export::SVG;
use APPRIS::Export::TXT;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Utils::Constant qw(
        $API_VERSION
        $VERSION
        $DATE
);

my ($API_VERSION) = $APPRIS::Utils::Constant::API_VERSION;
my ($VERSION) = undef; #$APPRIS::Utils::Constant::VERSION;
my ($DATE) = $APPRIS::Utils::Constant::DATE;

{
    # Encapsulated class data
    #___________________________________________________________
    my %_attr_data = # DEFAULT
		(
		);
    #_____________________________________________________________

	# Classwide default value for a specified object attribute
	sub _default_for {
		my ($self, $attr) = @_;
		$_attr_data{$attr};
	}

	# List of names of all specified object attributes
	sub _standard_keys {
		keys %_attr_data;
	}	
}
=head2 new

  Example :

    APPRIS::Exporter->new()

  Description: Will load the correct versions of the appris
               databases for the software release it can find on a
               database instance into the registry.

  Exceptions : None.
  Status     : Stable

=cut

sub new {
	my ($caller, %args) = @_;
	
	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;

	foreach my $attrname ($self->_standard_keys) {
		my ($attr) = "-".$attrname;
		if (exists $args{$attr} && defined $args{$attr}) {
			$self->{$attrname} = $args{$attr};
		} else {
			$self->{$attrname} = $self->_default_for($attrname);
		}
	}

	return $self;
}

=head2 software_version
  
  get the software version.
  
  Args       : none
  ReturnType : int
  Status     : At Risk
  
=cut
  
sub software_version {
	my ($self) = @_;
	return $API_VERSION;
}

=head2 date
  
  get the date of exported data.
  
  Args       : none
  ReturnType : string
  Status     : At Risk
  
=cut
  
sub date {
	my ($self) = @_;
	return $DATE;
}

=head2 version
  
  get the version of exported data.
  
  Args       : none
  ReturnType : string
  Status     : At Risk
  
=cut
  
sub version {
	my ($self) = @_;
	return $VERSION;
}

=head2 get_gtf_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : List of sources
  Example    : $annot = $exporter->get_gtf_annotations($feature,$sources);
  Description: Retrieves text as GFF format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_gtf_annotations {
    my ($self, $feature, $sources) = @_;
    my ($output) = '';

	if ($feature and (ref($feature) ne 'ARRAY')) {
    	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				$output .= APPRIS::Export::GTF::get_trans_annotations($transcript, $sources);
			}
    	}
    	elsif ($feature->isa("APPRIS::Transcript")) {
    		$output .= APPRIS::Export::GTF::get_trans_annotations($feature, $sources);
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
    }
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					$output .= APPRIS::Export::GTF::get_trans_annotations($transcript, $sources);
				}
	    	}
	    	elsif ($feat->isa("APPRIS::Transcript")) {
	    		$output .= APPRIS::Export::GTF::get_trans_annotations($feat, $sources);
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}    		
    	}
	}
    else {
		warning('Argument must be define');
   	}
	return $output;
}

=head2 get_gff3_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : List of sources
  Example    : $annot = $exporter->get_gff3_annotations($feature,$sources);
  Description: Retrieves text as GFF format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_gff3_annotations {
    my ($self, $feature, $sources) = @_;
    my ($output) = '';

	if ($feature and (ref($feature) ne 'ARRAY')) {
    	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				$output .= APPRIS::Export::GFF3::get_trans_annotations($transcript, $sources);
			}
    	}
    	elsif ($feature->isa("APPRIS::Transcript")) {
    		$output .= APPRIS::Export::GFF3::get_trans_annotations($feature, $sources);
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
    }
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					$output .= APPRIS::Export::GFF3::get_trans_annotations($transcript, $sources);
				}
	    	}
	    	elsif ($feat->isa("APPRIS::Transcript")) {
	    		$output .= APPRIS::Export::GFF3::get_trans_annotations($feat, $sources);
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}    		
    	}
	}
    else {
		warning('Argument must be define');
   	}
	return $output;
}

=head2 get_bed_annotations

  Arg [1]    : Listref of APPRIS::Gene or undef
  Arg [2]    : String - $sources
               List of sources ('all', ... )
  Arg [3]    : String - $position (optional)
               genome position (chr22:20116979-20137016)
  Example    : $annot = $exporter->get_bed_annotations($feature,'chr22:20116979-20137016','no','appris');
  Description: Retrieves text as BED format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_bed_annotations {
    my ($self, $feature, $sources, $position) = @_;
    my ($output) = '';

	# IMPORTANT:
	# If we have a list of objects, BED annotations only for the first transcript.
	#if ($feature and (ref($feature) eq 'ARRAY') ) { 
	#	$feature = $feature->[0];
	#}
	
	# Get position if its not defined
	unless (defined $position) {
		my ($chr);
		my ($start);
		my ($end);
		if ($feature and (ref($feature) ne 'ARRAY')) {			
	    	if ($feature->isa("APPRIS::Gene") or $feature->isa("APPRIS::Transcript")) {
				$chr = $feature->chromosome;
				$start = $feature->start;
				$end = $feature->end;	    		
	    	}
		}
		elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    		foreach my $feat (@{$feature}) {
	    		if ($feat->isa("APPRIS::Gene") or $feat->isa("APPRIS::Transcript")) {	    			
					$chr = $feat->chromosome unless (defined $chr);
					$start = $feat->start unless (defined $start);
					$end = $feat->end unless (defined $end);
	    			
	    			$start = $feat->start if ($feat->start < $start);
	    			$end = $feat->end if ($feat->end > $end);
	    		}
    		}
		}
	    if ( defined $chr and defined $start and defined $end ) {
	    	$position = "$chr:$start\-$end";
	    }		
	}
	if (defined $position and ($position ne '')) {
		$output .= APPRIS::Export::BED::get_annotations($feature, $position, $sources);		
	}
    else {
		warning('Argument must be define');
    }
	return $output;
}

=head2 get_bed12_annotations

  Arg [1]    : Listref of APPRIS::Gene or undef
  Arg [2]    : String - $sources
               List of sources ('all', ... )
  Arg [3]    : String - $position (optional)
               genome position (chr22:20116979-20137016)
  Example    : $annot = $exporter->get_bed_annotations($feature,'chr22:20116979-20137016','no','appris');
  Description: Retrieves text as BED format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_bed12_annotations {
    my ($self, $feature, $sources, $position) = @_;
    my ($output) = '';
    
	# IMPORTANT:
	# If we have a list of objects, BED annotations only for the first transcript.
	#if ($feature and (ref($feature) eq 'ARRAY') ) { 
	#	$feature = $feature->[0];
	#}
	
	# Get position if its not defined
	unless (defined $position) {
		my ($chr);
		my ($start);
		my ($end);
		if ($feature and (ref($feature) ne 'ARRAY')) {			
	    	if ($feature->isa("APPRIS::Gene") or $feature->isa("APPRIS::Transcript")) {
				$chr = $feature->chromosome;
				$start = $feature->start;
				$end = $feature->end;	    		
	    	}
		}
		elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    		foreach my $feat (@{$feature}) {
	    		if ($feat->isa("APPRIS::Gene") or $feat->isa("APPRIS::Transcript")) {	    			
					$chr = $feat->chromosome unless (defined $chr);
					$start = $feat->start unless (defined $start);
					$end = $feat->end unless (defined $end);
	    			
	    			$start = $feat->start if ($feat->start < $start);
	    			$end = $feat->end if ($feat->end > $end);
	    		}
    		}
		}
	    if ( defined $chr and defined $start and defined $end ) {
	    	$position = "$chr:$start\-$end";
	    }		
	}
	if (defined $position and ($position ne '')) {
		$output .= APPRIS::Export::BED12::get_annotations($feature, $position, $sources);		
	}
    else {
		warning('Argument must be define');
    }
	return $output;
}

=head2 get_json_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : List of sources
  Example    : $annot = $exporter->get_json_annotations($feature,$sources);
  Description: Retrieves text as JSON format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_json_annotations {
    my ($self, $feature, $sources) = @_;
    my ($output) = '';

	my ($gtf_output) = $self->get_gtf_annotations($feature, $sources);
	$output .= APPRIS::Export::JSON::get_annotations($gtf_output);
	
	return $output;
}

=head2 get_tsv_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : List of sources
  Arg [3]    : List of residues positions
  Example    : $annot = $exporter->get_tsv_annotations($feature,$params);
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

#sub get_tsv_annotations {
#    my ($self,$feature,$sources) = @_;
#    my ($output) = '';
#
#    if ( defined $feature ) {
#    	$output .= APPRIS::Export::TXT::get_tsv_annotations($feature,$sources);
#    }
#    else {
#		warning('Argument must be define');
#   	}
#	return $output;
#}
sub get_tsv_annotations {
    my ($self, $feature, $sources, $res) = @_;
    my ($output) = '';

	if ($feature and (ref($feature) ne 'ARRAY')) {
    	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				$output .= APPRIS::Export::TXT::get_trans_annotations($transcript, $sources, $res);
			}
    	}
    	elsif ($feature->isa("APPRIS::Transcript")) {
    		$output .= APPRIS::Export::TXT::get_trans_annotations($feature, $sources, $res);
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
    }
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					$output .= APPRIS::Export::TXT::get_trans_annotations($transcript, $sources, $res);
				}
	    	}
	    	elsif ($feat->isa("APPRIS::Transcript")) {
	    		$output .= APPRIS::Export::TXT::get_trans_annotations($feat, $sources, $res);
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}    		
    	}
	}
    else {
		warning('Argument must be define');
   	}
	return $output;
}

=head2 get_raw_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : List of sources
  Example    : $annot = $exporter->get_raw_annotations($feature,$params);
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_raw_annotations {
    my ($self,$feature,$sources) = @_;
    my ($output) = '';

    if ( defined $feature ) {
    	$output .= APPRIS::Export::TXT::get_raw_annotations($feature,$sources);
    }
    else {
		warning('Argument must be define');
   	}
	return $output;
}

=head2 get_json_results

  Arg [1]    : Hash of Raw result or undef
  Arg [2]    : List of sources  
  Example    : $annot = $exporter->get_json_results($feature);
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_json_results {
    my ($self,$report,$sources) = @_;
    my ($output) = '';

    if (defined $report) {
    	$output .= APPRIS::Export::JSON::create_array($report,$sources);
    }
    else {
		warning('Argument must be define');
   	}
	return $output;
}

=head2 get_raw_results

  Arg [1]    : Hash of Raw result or undef
  Arg [2]    : List of sources  
  Example    : $annot = $exporter->get_raw_results($feature);
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_raw_results {
    my ($self,$report,$sources) = @_;
    my ($output) = '';

    if (defined $report) {
    	$output .= APPRIS::Export::TXT::get_raw_results($report,$sources);
    }
    else {
		warning('Argument must be define');
   	}
	return $output;
}

=head2 get_tsv_results

  Arg [1]    : Hash of Raw result or undef
  Arg [2]    : List of sources  
  Example    : $annot = $exporter->get_tsv_results($feature);
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_tsv_results {
    my ($self,$report,$sources) = @_;
    my ($output) = '';

    if (defined $report) {
    	$output .= APPRIS::Export::TXT::get_tsv_results($report,$sources);
    }
    else {
		warning('Argument must be define');
   	}
	return $output;
}

=head2 get_seq_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : String $type
               type of sequence ('na' or 'aa')
  Arg [3]    : (optional) String $format
               format of output (by default 'fasta')
  Example    : $annot = $exporter->get_seq_annotations($feature,'aa');
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_seq_annotations {
	my ($self, $feature, $type, $format) = @_;
    my ($string) = '';
    my ($output) = '';

	if ($feature and (ref($feature) ne 'ARRAY')) {
    	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				$string .= APPRIS::Export::SEQ::get_trans_annotations($transcript, $type);
			}
    	}
    	elsif ($feature->isa("APPRIS::Transcript")) {
    		$string .= APPRIS::Export::SEQ::get_trans_annotations($feature, $type);
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
    }
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					$string .= APPRIS::Export::SEQ::get_trans_annotations($transcript, $type);
				}
	    	}
	    	elsif ($feat->isa("APPRIS::Transcript")) {
	    		$string .= APPRIS::Export::SEQ::get_trans_annotations($feat, $type);
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}    		
    	}
	}
    else {
		warning('Argument must be define');
   	}
   	
	if ( defined $format and ($format eq 'json') and ($string ne '') ) {
		my ($report);
		my ($in) = Bio::SeqIO-> new(
								-fh     => IO::String->new($string),
								-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() ) {
			my ($seq_id) = $seq->id; $seq_id =~ s/^([^\|]*)\|[^\$]*$/$1/g;				
			my ($seq_seq) = $seq->seq;
			push(@{$report->{'seqs'}}, {
				'id'		=> $seq_id,
				'seq'		=> $seq->seq,
				'length'	=> length($seq->seq),
			});
		}
		$report->{'num'} = scalar(@{$report->{'seqs'}});
		$output .= APPRIS::Export::JSON::create($report);
	}
	elsif ( $string ne '' ) {		
		$output = $string;
	} 	
	return $output;
	
} # end get_seq_annotations

=head2 get_res_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : String - $soure List of sources  
  Arg [3]    : Int    - $res Residue position  
  Example    : $annot = $exporter->get_res_annotations($feature,'aa');
  Description: Retrieves residues of methods.
  Returntype : String or undef

=cut

sub get_res_annotations {
	my ($self, $feature, $sources, $inres) = @_;
	my ($report);
    my ($output) = '';

	if ($feature and (ref($feature) ne 'ARRAY')) {
    	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				my ($rep) = APPRIS::Export::RES::get_trans_annotations($transcript, $sources, $inres);				
				push(@{$report}, $rep) if ( defined $rep );
			}
    	}
    	elsif ($feature->isa("APPRIS::Transcript")) {
    		my ($rep) = APPRIS::Export::RES::get_trans_annotations($feature, $sources, $inres);
				push(@{$report}, $rep) if ( defined $rep );
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
    }
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					my ($rep) = APPRIS::Export::RES::get_trans_annotations($transcript, $sources, $inres);
				push(@{$report}, $rep) if ( defined $rep );
				}
	    	}
	    	elsif ($feat->isa("APPRIS::Transcript")) {
	    		my ($rep) = APPRIS::Export::RES::get_trans_annotations($feat, $sources, $inres);
				push(@{$report}, $rep) if ( defined $rep );
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}    		
    	}
	}
    else {
		warning('Argument must be define');
   	}
	#if ( defined $report ) {
		$output .= APPRIS::Export::JSON::create($report);
	#}
   	
	return $output;
	
} # end get_res_annotations

=head2 get_cds_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : String - $soure List of sources  
  Arg [3]    : Int    - $res Residue position  
  Example    : $annot = $exporter->get_cds_annotations($feature,'aa');
  Description: Retrieves residues of methods.
  Returntype : String or undef

=cut

sub get_cds_annotations {
	my ($self, $feature, $sources, $inres) = @_;
	my ($report);
    my ($output) = '';

	if ($feature and (ref($feature) ne 'ARRAY')) {
    	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				my ($rep) = APPRIS::Export::CDS::get_trans_annotations($transcript, $sources, $inres);				
				push(@{$report}, $rep) if ( defined $rep );
			}
    	}
    	elsif ($feature->isa("APPRIS::Transcript")) {
    		my ($rep) = APPRIS::Export::CDS::get_trans_annotations($feature, $sources, $inres);
				push(@{$report}, $rep) if ( defined $rep );
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
    }
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					my ($rep) = APPRIS::Export::CDS::get_trans_annotations($transcript, $sources, $inres);
				push(@{$report}, $rep) if ( defined $rep );
				}
	    	}
	    	elsif ($feat->isa("APPRIS::Transcript")) {
	    		my ($rep) = APPRIS::Export::CDS::get_trans_annotations($feat, $sources, $inres);
				push(@{$report}, $rep) if ( defined $rep );
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}    		
    	}
	}
    else {
		warning('Argument must be define');
   	}
	#if ( defined $report ) {
		$output .= APPRIS::Export::JSON::create($report);
	#}
   	
	return $output;
	
} # end get_cds_annotations

=head2 get_img_seq_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Example    : $annot = $exporter->get_img_seq_annotations($feature);
  Description: Retrieves image of sequences.
  Returntype : String or undef

=cut

sub get_img_seq_annotations {
	my ($self, $feature) = @_;
    my ($output) = '';

	if ($feature and (ref($feature) ne 'ARRAY')) {
    	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				$output .= APPRIS::Export::SVG::get_trans_annotations($transcript);
			}
    	}
    	elsif ($feature->isa("APPRIS::Transcript")) {
    		$output .= APPRIS::Export::SVG::get_trans_annotations($feature);
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
    }
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
    	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					$output .= APPRIS::Export::SVG::get_trans_annotations($transcript);
				}
	    	}
	    	elsif ($feat->isa("APPRIS::Transcript")) {
	    		$output .= APPRIS::Export::SVG::get_trans_annotations($feat);
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}    		
    	}
	}
    else {
		warning('Argument must be define');
   	}
	return $output;
}

sub DESTROY {}

1;
