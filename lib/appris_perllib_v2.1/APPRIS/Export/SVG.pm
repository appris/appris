=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::SEQ - Utility functions for info exporting

=head1 SYNOPSIS

  use APPRIS::Export::SEQ
    qw(
       get_trans_annotations
     );

  or to get all methods just

  use APPRIS::Export::SEQ;

  eval { get_trans_annotations($feature,$params) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

Retrieves sequences of transcript as fasta format.

=head1 METHODS

=cut

package APPRIS::Export::SVG;

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use SVG;

use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Export::RES qw( get_trans_residues );
use APPRIS::Export::SEQ qw(
	$NUM_RESIDUES
	get_trans_annotations 
);

my ($TYPE) = 'aa';
my ($FORMAT) = 'fasta';
my ($NUM_RESIDUES) = $APPRIS::Export::SEQ::NUM_RESIDUES;

my ($SVG_STYLE) = qq|
	text {
		font-family:Courier New;
	}
	rect.functional_residue {
		fill: #d29123;
		fill-opacity: 1;
	}
	rect.functional_domain {
		fill: #769c02;
		fill-opacity: 0.5;
	}
	rect.transmembrane_signal {
		fill: #ff00ff;
		fill-opacity: 0.4;
	}
	rect.peptide_mitochondrial_signal, rect.signal_peptide, rect.mitochondrial_signal {
		fill: #004142;
		fill-opacity: 0.4;
	}
	.homologous_structure {
		/*fill: #c60806;*/
		fill: #800000;
	}
|;

my ($SVG_STYLE_INDEX) = qq|
	text.index {
		font-size:11px;		
	}
|;
my ($SVG_STYLE_SEQ) = qq|
	text.seq {
		font-size: 14px;
		/*font-weight: bold;*/
	}
|;
  
my ($PANEL_SEP) = 20; # Pixels that separe the line panel (aa + method)

my ($TEXT_X_INIT_I) = 10; # Init pixel of x-coord for index num
my ($TEXT_X_INIT_S) = 45; # Init pixel of x-coord for seq
my ($TEXT_Y_INIT) = 15; # Init pixel of y-coord

my ($M_LINE_X_INIT) = $TEXT_X_INIT_S; # Init pixel of x-coord for method lines
my ($M_LINE_Y_INIT) = $TEXT_Y_INIT + 2; # Init pixel of y-coord for method lines

my ($CHAR_WIDTH) = 9;
my ($CHAR_HEIGHT) = 10;
my ($M_LINE_SEP) = 3;

my ($METHOD_LINE) = {
	'transmembrane_signal'			=> 1,
	'peptide_mitochondrial_signal'	=> 2,
	'signal_peptide'				=> 2,
	'mitochondrial_signal'			=> 2,	
};
my ($METHOD_SUFFIX) = {
	'functional_residue'			=> 'f',
	'homologous_structure'			=> 'm',
	'functional_domain'				=> 's',
	'transmembrane_signal'			=> 't',
	'peptide_mitochondrial_signal'	=> 'c',
	'signal_peptide'				=> 'sp',
	'mitochondrial_signal'			=> 'tp',
#	'vertebrate_conservation'		=> 'o',
#	'neutral_evolution'				=> 'i',
};
my ($w_proportion) = $CHAR_WIDTH+1;
my ($h_proportion) = $PANEL_SEP+2;
my ($SVG_WIDTH) = $w_proportion*$NUM_RESIDUES;
my ($SVG_HEIGHT) = '100%';

my ($CONTROL_UNIQUE_RES);

=head2 get_trans_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : String $type
               type of sequence ('na' or 'aa')
  Arg [3]    : (optional) String $format
               format of output (by default 'fasta')
  Example    : $annot = $exporter->get_trans_annotations($feature,'aa');
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_trans_annotations {
    my ($feature) = @_;
    my ($output) = '';

	
	if (defined $feature and $feature->stable_id) {
		my ($id) = $feature->stable_id;
		
		# Help to control possible errors
		my ($length);
		if ($feature->translate and $feature->translate->sequence) {
			$length = length($feature->translate->sequence);
			my ($total_lines) = int($length / $NUM_RESIDUES) +1;
			$SVG_HEIGHT = $h_proportion*$total_lines;			
		}

		my ($fasta_seq) = APPRIS::Export::SEQ::get_trans_annotations($feature,$TYPE,$FORMAT);
		if ( defined $fasta_seq and ($fasta_seq ne '') ) {
			my ($num_line) = 0;
			my ($svg) = SVG->new(
							width 				=> $SVG_WIDTH,
							height				=> $SVG_HEIGHT,
							preserveAspectRatio => "xMidYMid meet",
							viewBox 			=> "0 0 $SVG_WIDTH $SVG_HEIGHT"
			);
   			
			my ($style) = $svg->style(type => 'text/css');
			$style->CDATA($SVG_STYLE.$SVG_STYLE_INDEX.$SVG_STYLE_SEQ);
	        		
			my ($method_positions) = APPRIS::Export::RES::get_trans_annotations($feature);
			
			my ($method_attrs) = get_method_figures($id, $method_positions, $length);
			while ( my ($type, $attrs) = each(%{$method_attrs}) ) {
				if ( $type eq 'line' ) {
					foreach my $attr (@{$attrs}) {
						my ($l) = $svg->line(%{$attr});
					}	
				}
				elsif ( $type eq 'rect' ) {
					foreach my $attr (@{$attrs}) {
						my ($l) = $svg->rect(%{$attr});
					}					
				}
			}
			
						
			# Add aminoacids
			my ($pos) = 1;			
			while ( $fasta_seq =~ /([^\n]*)\n/g ) {
				my ($line) = $1;
				unless ( $line =~ /^\>/ ) {
					my ($text_id) = 'text' . ($num_line+1);
					my ($index_id) = 'index' . ($num_line+1);

					my ($x_i) = $TEXT_X_INIT_I;
					my ($y_i) = $TEXT_Y_INIT + ($num_line*$PANEL_SEP);
					my ($y_s) = $y_i;						
					
					my ($index_cont) = (($num_line*$NUM_RESIDUES)+1).":";
					my ($index) = $svg->text(id => $index_id, class => 'index', x => $x_i, y => $y_i);
					$index->text(-type => 'span')->cdata($index_cont);
					
					my ($text) = $svg->text(id => $text_id, class => 'seq', y => $y_s);					
					my (@aas) = split(//,$line);
					for ( my $a = 0; $a < scalar(@aas); $a++ ) {
						my ($c) = $aas[$a];
						my ($x_s) = $TEXT_X_INIT_S + ($CHAR_WIDTH*$a);
						my (%attr) = (
							id => $pos,
							x => $x_s,
							-type => 'span'
						);
						
						if ( exists $method_attrs->{'class'} and exists $method_attrs->{'class'}->{$pos} and defined $method_attrs->{'class'}->{$pos} ) {
							$attr{class} = $method_attrs->{'class'}->{$pos}; 
						}
						$text->text(%attr)->cdata($c);
						
						$pos++;
					}					
					$num_line++;				
				}
			}
			#print STDERR "\nSVG:\n".$svg->xmlify;
			$output .= $svg->xmlify;
		}
	    else {
			throw('Argument has to be sequence');
	   	}
	}
    else {
		throw('Argument has to be sequence');
   	}
	return $output;
}

sub get_method_figures {
	my ($id,$method_positions,$length) = @_;	
	my ($attrs);
	
	if ( exists $method_positions->{$id} and defined $method_positions->{$id} ) {
		my ($method_pos) = $method_positions->{$id};
		while (my ($met,$pos) = each(%{$method_pos}) ) {
			my ($m_line) = $METHOD_LINE->{$met};			
			foreach my $p (@{$pos}) {
				my ($p_s) = $p->{'start'};
				my ($p_e) = $p->{'end'};
				
				# Control the length of residues (bad for Matador3d: ENST00000300482)
				if ( $p_e > $length ) {
					$p_e = $length;
				}
				
				for ( my $p = $p_s; $p <= $p_e; $p++ ) {
					my ($id) = $p;
					
					my ($num_pos_s) = $p % $NUM_RESIDUES;
					my ($num_line_s) = int( $p / $NUM_RESIDUES );
					if ($num_pos_s == 0) { $num_line_s -= 1; $num_pos_s = $NUM_RESIDUES; };
										
					$attrs = get_method_attrs($id, $met, $num_pos_s, $num_line_s, $attrs);
				}
			}
		}
	}
		
	return $attrs;
}

sub get_method_attrs {
	my ($pos, $met, $n_pos, $line, $attrs) = @_;
	
	if ( $met eq 'functional_residue' ) {
		my ($i) = $METHOD_SUFFIX->{$met}.'-'.$pos;
		unless ( exists $CONTROL_UNIQUE_RES->{$i} ) {
			$CONTROL_UNIQUE_RES->{$i} = 1;
			my ($w) = $CHAR_WIDTH;
			my ($h) = $CHAR_HEIGHT;
			my ($x) = $TEXT_X_INIT_S + ($CHAR_WIDTH*($n_pos-1));
			my ($y) = $TEXT_Y_INIT + ($line*$PANEL_SEP) - $h + 0.5;
			my ($attr) = {
				'id'		=> $i,
				'class'		=> $met,
				'x'			=> $x,
				'y'			=> $y,
				'width'		=> $w,
				'height'	=> $h
			};
			push(@{$attrs->{'rect'}},$attr);
		}	
	}
	elsif ( $met eq 'functional_domain' ) {
		my ($i) = $METHOD_SUFFIX->{$met}.'-'.$pos;
		unless ( exists $CONTROL_UNIQUE_RES->{$i} ) {
			$CONTROL_UNIQUE_RES->{$i} = 1;
			my ($w) = $CHAR_WIDTH;
			my ($h) = $CHAR_HEIGHT;
			my ($x) = $TEXT_X_INIT_S + ($CHAR_WIDTH*($n_pos-1));
			my ($y) = $TEXT_Y_INIT + ($line*$PANEL_SEP) - $h + 0.5;
			my ($attr) = {
				'id'		=> $i,
				'class'		=> $met,
				'x'			=> $x,
				'y'			=> $y,
				'width'		=> $w,
				'height'	=> $h
			};
			push(@{$attrs->{'rect'}},$attr);			
		}		
	}
	elsif (	$met eq 'transmembrane_signal' ) {
		my ($i) = $METHOD_SUFFIX->{$met}.'-'.$pos;
		unless ( exists $CONTROL_UNIQUE_RES->{$i} ) {
			$CONTROL_UNIQUE_RES->{$i} = 1;
			my ($w) = $CHAR_WIDTH;
			my ($h) = $CHAR_HEIGHT;
			my ($x) = $TEXT_X_INIT_S + ($CHAR_WIDTH*($n_pos-1));
			my ($y) = $TEXT_Y_INIT + ($line*$PANEL_SEP) - $h + 0.5;
			my ($attr) = {
				'id'		=> $i,
				'class'		=> $met,
				'x'			=> $x,
				'y'			=> $y,
				'width'		=> $w,
				'height'	=> $h
			};
			push(@{$attrs->{'rect'}},$attr);
		}
	}	
	elsif (	($met eq 'peptide_mitochondrial_signal') or ($met eq 'signal_peptide') or ($met eq 'mitochondrial_signal') ) {
		$met = 'peptide_mitochondrial_signal'; 
		my ($i) = $METHOD_SUFFIX->{$met}.'-'.$pos;
		unless ( exists $CONTROL_UNIQUE_RES->{$i} ) {
			$CONTROL_UNIQUE_RES->{$i} = 1;	
			my ($w) = $CHAR_WIDTH;
			my ($h) = $CHAR_HEIGHT;
			my ($x) = $TEXT_X_INIT_S + ($CHAR_WIDTH*($n_pos-1));
			my ($y) = $TEXT_Y_INIT + ($line*$PANEL_SEP) - $h + 0.5;
			my ($attr) = {
				'id'		=> $i,
				'class'		=> $met,
				'x'			=> $x,
				'y'			=> $y,
				'width'		=> $w,
				'height'	=> $h
			};
			push(@{$attrs->{'rect'}},$attr);
		}
	}	
	elsif ( $met eq 'homologous_structure' ) {
		$attrs->{'class'}->{$pos} = $met;									
	}

	return $attrs;	
}


1;
