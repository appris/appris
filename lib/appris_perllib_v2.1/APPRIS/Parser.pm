=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

Parser - Module to handle results of methods of APPRIS pipeline.

=head1 SYNOPSIS

  use APPRIS::Parser
    qw(
       parse_firestar
     );

  or to get all methods just

  use APPRIS::Parser;

  eval { parse_firestar($result) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

Module to handle results of methods of APPRIS pipeline.

=head1 METHODS

=cut

package APPRIS::Parser;

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use IO::String;
use Scalar::Util qw(reftype);

use APPRIS::Gene;
use APPRIS::Transcript;
use APPRIS::Translation;
use APPRIS::Exon;
use APPRIS::CDS;
use APPRIS::Codon;
use APPRIS::XrefEntry;
use APPRIS::Analysis;
use APPRIS::Analysis::Region;
use APPRIS::Analysis::Firestar;
use APPRIS::Analysis::Matador3D;
use APPRIS::Analysis::SPADE;
use APPRIS::Analysis::INERTIA;
use APPRIS::Analysis::CRASH;
use APPRIS::Analysis::THUMP;
use APPRIS::Analysis::CORSAIR;
use APPRIS::Analysis::PROTEO;
use APPRIS::Analysis::APPRIS;
use APPRIS::Utils::File qw(parse_file);
use APPRIS::Utils::ProCDS qw(sort_cds get_coords_from_residue);
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Utils::Constant qw(
	$OK_LABEL
	$NO_LABEL
	$UNKNOWN_LABEL

	$FIRESTAR_ACCEPT_LABEL
	$FIRESTAR_REJECT_LABEL
);

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	parse_infiles
	parse_gencode
	parse_transl_data
	parse_firestar_rst
	parse_firestar
	parse_matador3d_rst
	parse_matador3d
	parse_spade_rst
	parse_spade
	parse_inertia
	parse_corsair_rst
	parse_corsair
	parse_crash_rst
	parse_crash
	parse_thump_rst
	parse_thump
	parse_proteo_rst
	parse_proteo
	parse_appris_rst
	parse_appris
	parse_appris_methods
	create_appris_entity
);

sub parse_infiles($;$;$);
sub _parse_dataline($);
sub _parse_indata($);
sub _parse_refseq_gbk($);
sub _parse_inseq_transc($);
sub _parse_inseq_transl($);
sub parse_gencode($;$;$);
sub parse_transl_data($);
sub _get_id_version($);
sub _parse_gencode_data($);
sub _parse_gencode_seq_transc($);
sub _parse_gencode_seq_transl($);
sub _parse_seq_data($);
sub _fetch_transc_objects($$;$;$);
sub _fetch_transl_objects($$;$);
sub parse_firestar_rst($);
sub parse_firestar($$);
sub parse_matador3d_rst($);
sub parse_matador3d($$);
sub parse_spade_rst($);
sub parse_spade($$);
sub parse_inertia($$);
sub _parse_inertia_file($$\$);
sub _parse_omega_file($$\$);
sub parse_corsair_rst($);
sub parse_corsair($$);
sub parse_crash_rst($);
sub parse_crash($$);
sub parse_thump_rst($);
sub parse_thump($$);
sub parse_proteo_rst($);
sub parse_proteo($$);
sub parse_appris_rst($$);
sub parse_appris($$$);
sub parse_appris_methods($$$$$$$$;$;$;$);
sub create_appris_entity($$$$$$$$$$$$$$);


=head2 parse_infiles

  Arg [1]    : string $type
               Type of infiles
  Arg [2]    : string $data_file
               File of data as GTF format
  Arg [3]    : string $transc_file (optional)
               File of transcript sequence as Fasta format
  Arg [4]    : string $transl_file (optional)
               File of translation sequence as Fasta format
  Example    : use APPRIS::Parser qw(parse_infiles);
               parse_infiles($data_file, $transc_file, $transl_file);
  Description: Parse input files of gencode/ensembl.
  Returntype : Listref of APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_infiles($;$;$)
{
	my ($data_file, $transc_file, $transl_file) = @_;
	my ($entity_list, $raw_data_list);
	my ($data_cont);
	my ($transc_cont);
	my ($transl_cont);
	my ($index_genes);
	my ($index) = 0;	
		
	# Parse data/sequences
	if ( _is_refseq($data_file) ) {
		$data_cont = _parse_indata_refseq($data_file);
		$transc_cont = _parse_inseq_transc($transc_file) if (defined $transc_file);
		$transl_cont = _parse_inseq_transl($transl_file) if (defined $transl_file);
	}
	else {
		$data_cont = _parse_indata($data_file);
		$transc_cont = _parse_inseq_transc($transc_file) if (defined $transc_file);
		$transl_cont = _parse_inseq_transl($transl_file) if (defined $transl_file);
	}
	
	# Scan genes
	while ( my ($gene_id, $gene_features) = each(%{$data_cont}) )
	{
		my ($xref_identities);
		my ($transcripts, $index_transcripts) = _fetch_transc_objects($gene_id, $gene_features->{'transcripts'}, $transc_cont, $transl_cont);
		
		# Create gene object
		my ($gene) = APPRIS::Gene->new
		(
			-stable_id	=> $gene_id,
			-chr		=> $gene_features->{'chr'},
			-start		=> $gene_features->{'start'},
			-end		=> $gene_features->{'end'},
			-strand		=> $gene_features->{'strand'},
			-biotype	=> $gene_features->{'biotype'},
			-status		=> $gene_features->{'status'},
			-source		=> $gene_features->{'source'},
			-level		=> $gene_features->{'level'},
			-version	=> $gene_features->{'version'}
		);

		# Xref identifiers
		if ( exists $gene_features->{'external_id'} and defined $gene_features->{'external_id'} ) {
			$gene->external_name($gene_features->{'external_id'});
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $gene_features->{'external_id'},
						-dbname			=> 'External_Id'
					)
			);
		}
		$gene->xref_identify($xref_identities) if (defined $xref_identities);
		$gene->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
		push(@{$entity_list},$gene) if (defined $gene);
		
		# Get raw data
		if ( exists $gene_features->{'raw'} and ($gene_features->{'raw'} ne '') ) {
			$raw_data_list->{$gene_id} =  $gene_features->{'raw'};					
		}
		
		# Get index genes
		$index_genes->{$gene_id} = $index; $index++; # Index the list of transcripts
	}
		
	return ($entity_list, $index_genes, $raw_data_list);
}

=head2 parse_refseq

  Arg [1]    : string $data_file
               File of refseq data as GenBank format
  Arg [2]    : string $transc_file (optional)
               File of refseq transcript sequence as Fasta format
  Arg [3]    : string $transl_file (optional)
               File of refseq translation sequence as Fasta format
  Example    : use APPRIS::Parser qw(parse_refseq);
               parse_refseq($data_file, $transc_file, $transl_file);
  Description: Parse data of refseq.
  Returntype : Listref of APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_refseq($)
{
	my ($data_file, $transc_file, $transl_file) = @_;
	my ($entity_list, $raw_data_list);
	my ($data_cont);
	my ($transc_cont);
	my ($transl_cont);
	my ($index_genes);
	my ($index) = 0;	
		
	# Parse data/sequences
	$data_cont = _parse_refseq_gbk($data_file);
	$transc_cont = _parse_inseq_transc($transc_file) if (defined $transc_file);
	$transl_cont = _parse_inseq_transl($transl_file) if (defined $transl_file);

	# Scan genes
	while ( my ($gene_id, $gene_features) = each(%{$data_cont}) )
	{
		my ($xref_identities);
		my ($transcripts, $index_transcripts) = _fetch_transc_objects($gene_id, $gene_features->{'transcripts'}, $transc_cont, $transl_cont);
		
		# Create gene object
		my ($gene) = APPRIS::Gene->new
		(
			-stable_id	=> $gene_id,
			-chr		=> $gene_features->{'chr'},
			-start		=> $gene_features->{'start'},
			-end		=> $gene_features->{'end'},
			-strand		=> $gene_features->{'strand'},
			-biotype	=> $gene_features->{'biotype'},
			-status		=> $gene_features->{'status'},
			-source		=> $gene_features->{'source'},
			-level		=> $gene_features->{'level'},
			-version	=> $gene_features->{'version'}
		);

		# Xref identifiers
		if ( exists $gene_features->{'external_id'} and defined $gene_features->{'external_id'} ) {
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $gene_features->{'external_id'},
						-dbname			=> 'External_Id'
					)
			);
		}
		$gene->xref_identify($xref_identities) if (defined $xref_identities);
		$gene->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
		push(@{$entity_list},$gene) if (defined $gene);
		
		# Get raw data
		if ( exists $gene_features->{'raw'} and ($gene_features->{'raw'} ne '') ) {
			$raw_data_list->{$gene_id} =  $gene_features->{'raw'};					
		}
		
		# Get index genes
		$index_genes->{$gene_id} = $index; $index++; # Index the list of transcripts
	}
		
	return ($entity_list, $index_genes, $raw_data_list);
}

=head2 parse_gencode

  Arg [1]    : string $data_file
               File of gencode data as GTF format
  Arg [2]    : string $transc_file (optional)
               File of gencode transcript sequence as Fasta format
  Arg [3]    : string $transl_file (optional)
               File of gencode translation sequence as Fasta format
  Example    : use APPRIS::Parser qw(parse_gencode);
               parse_gencode($data_file, $transc_file, $transl_file);
  Description: Parse data of gencode.
  Returntype : Listref of APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_gencode($;$;$)
{
	my ($data_file, $transc_file, $transl_file) = @_;
	my ($entity_list, $raw_data_list);
	my ($data_cont);
	my ($transc_cont);
	my ($transl_cont);
	my ($index_genes);
	my ($index) = 0;	
		
	# Parse data/sequences
	$data_cont = _parse_gencode_data($data_file);
	$transc_cont = _parse_gencode_seq_transc($transc_file) if (defined $transc_file);
	$transl_cont = _parse_gencode_seq_transl($transl_file) if (defined $transl_file);

	# Scan genes
	while ( my ($gene_id, $gene_features) = each(%{$data_cont}) )
	{
		my ($xref_identities);
		my ($transcripts, $index_transcripts) = _fetch_transc_objects($gene_id, $gene_features->{'transcripts'}, $transc_cont, $transl_cont);
		
		# Create gene object
		my ($gene) = APPRIS::Gene->new
		(
			-stable_id	=> $gene_id,
			-chr		=> $gene_features->{'chr'},
			-start		=> $gene_features->{'start'},
			-end		=> $gene_features->{'end'},
			-strand		=> $gene_features->{'strand'},
			-biotype	=> $gene_features->{'biotype'},
			-status		=> $gene_features->{'status'},
			-source		=> $gene_features->{'source'},
			-level		=> $gene_features->{'level'},
			-version	=> $gene_features->{'version'}
		);

		# Xref identifiers
		if ( exists $gene_features->{'external_id'} and defined $gene_features->{'external_id'} ) {
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $gene_features->{'external_id'},
						-dbname			=> 'External_Id'
					)
			);
		}
		$gene->xref_identify($xref_identities) if (defined $xref_identities);
		$gene->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
		push(@{$entity_list},$gene) if (defined $gene);
		
		# Get raw data
		if ( exists $gene_features->{'raw'} and ($gene_features->{'raw'} ne '') ) {
			$raw_data_list->{$gene_id} =  $gene_features->{'raw'};					
		}
		
		# Get index genes
		$index_genes->{$gene_id} = $index; $index++; # Index the list of transcripts
	}
		
	return ($entity_list, $index_genes, $raw_data_list);
}

=head2 parse_transl_data

  Arg [1]    : string $transl_file
               File of translation sequence as Fasta format
  Example    : use APPRIS::Parser qw(parse_transl_data);
               parse_transl_data($transl_file);
  Description: Parse data from sequence file.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_transl_data($)
{
	my ($transl_file) = @_;
	my ($entity_list);
		
	# Parse data/sequences
	my ($data_cont) = _parse_seq_data($transl_file) if (defined $transl_file);
	
	# Scan genes
	while ( my ($gene_id, $gene_features) = each(%{$data_cont}) )
	{
		my ($transcripts);
		my ($index_transcripts);
		my ($t_index) = 0;	
		
		# Scan transcripts
		if ( exists $gene_features->{'transcripts'} ) {
			while (my ($transcript_id, $transcript_features) = each(%{$gene_features->{'transcripts'}}) )
			{
				my ($xref_identities);				
				
				# Create transcript object
				my ($transcript) = APPRIS::Transcript->new
				(
					-stable_id	=> $transcript_id,
					-source		=> $gene_features->{'source'}
				);
				
				# add gene id
				if ( defined $gene_id ) {
					push(@{$xref_identities},
							APPRIS::XrefEntry->new
							(
								-id				=> $gene_id,
								-dbname			=> 'Gene_Id'
							)
					);
				}
				
				# add Xref
				if ( exists $transcript_features->{'name'} and defined $transcript_features->{'name'} ) {
					push(@{$xref_identities},
							APPRIS::XrefEntry->new
							(
								-id				=> $transcript_features->{'name'},
								-dbname			=> 'External_Id'
							)
					);
				}
				$transcript->xref_identify($xref_identities) if (defined $xref_identities);				
				
				# add Xref
				if ( exists $transcript_features->{'ccdsid'} and defined $transcript_features->{'ccdsid'} ) {
					push(@{$xref_identities},
							APPRIS::XrefEntry->new
							(
								-id				=> $transcript_features->{'ccdsid'},
								-dbname			=> 'CCDS'
							)
					);
				}
				$transcript->xref_identify($xref_identities) if (defined $xref_identities);

				# add Xref
				if ( exists $transcript_features->{'transc_ids'} and defined $transcript_features->{'transc_ids'} ) {
					foreach my $tids ( split(/\+/, $transcript_features->{'transc_ids'}) ) {
						if ( $tids =~ /^([^\:]+)\:([^\$]+)$/ ) {
							my ($source) = $1;
							my ($ids) = $2;
							foreach my $id ( split(/\,/, $ids) ) {
								if ( $id ne '' and $id ne '-' and $id ne '?' ) {
									push(@{$xref_identities},
											APPRIS::XrefEntry->new
											(
												-id				=> $id,
												-dbname			=> ucfirst($source).'_Transcript_Id'
											)
									);
								}
							}
						}
					}
				}
				$transcript->xref_identify($xref_identities) if (defined $xref_identities);				

				# Add translation
				my ($translate);
				if ( exists $transcript_features->{'seq'} and defined $transcript_features->{'seq'} ) {
					my ($translation_seq) = $transcript_features->{'seq'};
					
					# Create translation object
					$translate = APPRIS::Translation->new
					(
						-stable_id	=> $transcript_id,
					);			
					$translate->sequence($translation_seq);
				}
				
				$transcript->translate($translate) if (defined $translate);
					
				push(@{$transcripts}, $transcript) if (defined $transcript);
				$index_transcripts->{$transcript_id} = $t_index; $t_index++; # Index the list of transcripts		
			}
		}	
		
		# Create gene object
		my ($gene) = APPRIS::Gene->new
		(
			-stable_id	=> $gene_id,
			-source		=> $gene_features->{'source'}
		);		
		$gene->external_name($gene_features->{'name'}) if ( exists $gene_features->{'name'});

		# add Xref
		my ($xref_identities);
		if ( exists $gene_features->{'name'} and defined $gene_features->{'name'} ) {
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $gene_features->{'name'},
						-dbname			=> 'External_Id'
					)
			);
		}
		$gene->xref_identify($xref_identities) if (defined $xref_identities);

		# add Xref
		if ( exists $gene_features->{'gene_ids'} and defined $gene_features->{'gene_ids'} ) {
			foreach my $gids ( split(/\+/, $gene_features->{'gene_ids'}) ) {
				if ( $gids =~ /^([^\:]+)\:([^\$]+)$/ ) {
					my ($source) = $1;
					my ($ids) = $2;
					foreach my $id ( split(/\,/, $ids) ) {
						if ( $id ne '' and $id ne '-' and $id ne '?' ) {
							push(@{$xref_identities},
									APPRIS::XrefEntry->new
									(
										-id				=> $id,
										-dbname			=> ucfirst($source).'_Gene_Id'
									)
							);
						}
					}
				}
			}
		}
		$gene->xref_identify($xref_identities) if (defined $xref_identities);

		
		$gene->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
		push(@{$entity_list},$gene) if (defined $gene);		
	}

	return $entity_list;
}

=head2 parse_firestar_rst

  Arg [1]    : string $result
               Parse firestar result
  Example    : use APPRIS::Parser qw(parse_firestar);
               parse_firestar($result);
  Description: Parse output of firestar.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_firestar_rst($)
{
	my ($result) = @_;
	my ($transcript_result) = '';
	my ($init_trans_result) = 0;
	my ($cutoffs);	
		
	my (@results) = split( '\n', $result);
	foreach my $line (@results)
	{
		if ( $line eq '######') {
			# Init trans result
			$transcript_result = '';
			$init_trans_result = 1;
		}
		if ( $init_trans_result ) {
			$transcript_result .= $line."\n";		
		}

		#>>>OTTHUMT00000171822      29      46,47,48,49,67,68,69,70,98,99,100,101,166,170,171,394,411,412,413,448,452,480,481,492,495,497,498,499,504
		if ( $line=~/^>>>([^\t]+)\t+([^\t]+)\t+([^\$]+)$/ )
		{
			my ($id) = $1;
			my (@residue_list) = split(',', $3);

			# Get the peptide position and scores
			my ($residue_list_report);
			foreach my $residue_position (@residue_list)
			{
				if ( defined $residue_position and $residue_position ne '' )
				{
					#349     GLLCGGSAGSTVA   PLP[0.71,5.7,99.5]|Cat_Site_Atl[1.00,4,XXX]
					if ( $transcript_result =~ /$residue_position\t+([^\t]*)\t+([^\n]*)[^\>]*>>>$id/ )
					{
						my ($domain) = $1;
						my ($ligands) = $2; $ligands =~ s/^\s*//; $ligands =~ s/\s*$//; $ligands =~ s/\s+/\|/g; #$ligands = join('|',split(/\s+/,$ligands));
						push(@{$residue_list_report},{
								'residue'	=> $residue_position,
								'domain'	=> $domain,
								'ligands'	=> $ligands,
						});
					}
				}			
			}

			if(defined $residue_list_report and scalar(@{$residue_list_report})>0)
			{
				$cutoffs->{$id}->{'residues'} = $residue_list_report;
			}
			
			# Save result for each transcript
			$cutoffs->{$id}->{'result'} = $transcript_result;
			
			# Init trans result
			$transcript_result = ''; 
		}
		#C>>     OTTHUMT00000171822      3      46,47,48
		if ( $line=~/^C>>\t+([^\t]+)\t+([^\t]+)\t+([^\$]+)$/ )
		{
			my ($id) = $1;
			my (@residue_list) = split(',', $3);

			# Get the peptide position and scores
			foreach my $residue_position (@residue_list)
			{
				if ( defined $residue_position and $residue_position ne '' )
				{
					#349     GLLCGGSAGSTVA   PLP[0.71,6,99.5]|Cat_Site_Atl[1.00,4,XXX]					
					if ( $transcript_result =~ /$residue_position\t+([^\t]*)\t+([^\n]*)[^\>]*>>\t+$id/ )
					{
						my ($domain) = $1;
						my ($ligands) = $2; $ligands =~ s/^\s*//; $ligands =~ s/\s*$//; $ligands =~ s/\s+/\|/g; #$ligands = join('|',split(/\s+/,$ligands));
						push(@{$cutoffs->{$id}->{'residues'}},{
								'residue'	=> $residue_position,
								'domain'	=> $domain,
								'ligands'	=> $ligands,				
						});		
					}
				}
			}
			
			# Sort residues
			if ( exists $cutoffs->{$id}->{'residues'} and defined $cutoffs->{$id}->{'residues'} and
				 (scalar(@{$cutoffs->{$id}->{'residues'}}) > 0) )
			{
				my (@res_list) = @{$cutoffs->{$id}->{'residues'}};
				my (@sort_res_list) = sort { $a->{'residue'} <=> $b->{'residue'} } @res_list;
				$cutoffs->{$id}->{'residues'} = \@sort_res_list;				
			}

			# Save result for each transcript
			my ($trans_result) = '';
			$trans_result = $cutoffs->{$id}->{'result'} if ( exists $cutoffs->{$id}->{'result'} and defined $cutoffs->{$id}->{'result'} );
			$cutoffs->{$id}->{'result'} = $trans_result . $transcript_result;
			
			# Init result variable for the next
			$transcript_result = '';
		}
		
		#ACCEPT: ID\tTOTAL_SCORE\tTOTAL_MOTIFS\n
		# BEGIN: DEPRECATED
		#if ( $line =~ /^ACCEPT:\s*([^\t]+)\t([^\t]+)\t([^\n]+)\n*/ )
		if ( $line =~ /^F>>\t+([^\t]+)\t([^\t]+)\t([^\n]+)\n*/ )
		# END: DEPRECATED
		{
			my ($id) = $1;
			my ($total_score) = $2;
			my ($total_residues) = $3;

			if ( defined $id and ($id ne '') )
			{
				unless ( defined $total_residues and $total_residues ne '' )
					{ $total_residues = 0; }					
				$cutoffs->{$id}->{'num_residues'} = $total_residues;

				# Save result for each transcript
				my ($trans_result) = '';
				if ( exists $cutoffs->{$id}->{'result'} and defined $cutoffs->{$id}->{'result'} ) {
					$trans_result .= $cutoffs->{$id}->{'result'};
					$trans_result .= "----------------------------------------------------------------------\n";					
				}				
				$cutoffs->{$id}->{'result'} = $trans_result . $line;
			}
		}
		#REJECT: ID\tTOTAL_SCORE\tTOTAL_MOTIFS\n
		#if ( $line =~ /^REJECT:\s*([^\t]+)\t([^\t]+)\t([^\n]+)\n*/ )
		#{
		#	my ($id) = $1;
		#	my ($total_score) = $2;
		#	my ($total_residues) = $3;
		#
		#	if ( defined $id and ($id ne '') )
		#	{
		#		unless ( defined $total_residues and $total_residues ne '' )
		#			{ $total_residues = 0; }
		#			
		#		$cutoffs->{$id}->{'num_residues'} = $total_residues;
		#		$cutoffs->{$id}->{'functional_residue'} = $APPRIS::Utils::Constant::FIRESTAR_REJECT_LABEL;
		#		
		#		# Save result for each transcript
		#		my ($trans_result) = '';
		#		if ( exists $cutoffs->{$id}->{'result'} and defined $cutoffs->{$id}->{'result'} ) {
		#			$trans_result .= $cutoffs->{$id}->{'result'};
		#			$trans_result .= "----------------------------------------------------------------------\n";					
		#		}				
		#		$cutoffs->{$id}->{'result'} = $trans_result . $line;				
		#	}
		#}
	}
	$cutoffs->{'result'} = $result;
	
	return $cutoffs;
}

=head2 parse_firestar

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse firestar result
  Example    : use APPRIS::Parser qw(parse_firestar);
               parse_firestar($result);
  Description: Parse output of firestar.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_firestar($$)
{
	my ($gene, $result) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;	
	
	# Create hash object from result
	my ($cutoffs) = parse_firestar_rst($result);

	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
				
		# create method object
		#if ( exists $cutoffs->{$transcript_id} and exists $cutoffs->{$transcript_id}->{'functional_residue'} ) {
		if ( exists $cutoffs->{$transcript_id} ) {
			my ($report) = $cutoffs->{$transcript_id};
			my ($regions);			
			if ( $transcript->translate ) {
				if ( exists $report->{'residues'} ) {
					foreach my $re (@{$report->{'residues'}}) {
						my ($region);
						if ( $transcript->translate->cds ) {  # with CDS coords from GTF files
							my ($pro_coord) = get_coords_from_residue($transcript, $re->{'residue'});
							$region = APPRIS::Analysis::FirestarRegion->new (
											-residue	=> $re->{'residue'},
											-domain		=> $re->{'domain'},
											-ligands	=> $re->{'ligands'},
											-start		=> $pro_coord->start,
											-end		=> $pro_coord->end,
											-strand		=> $pro_coord->strand,
							);
						}
						else {
							$region = APPRIS::Analysis::FirestarRegion->new (
											-residue	=> $re->{'residue'},
											-domain		=> $re->{'domain'},
											-ligands	=> $re->{'ligands'},
							);							
						}
						push(@{$regions}, $region) if ( defined $region );
					}
				}
			}
						
			# create Analysis object (for trans)			
			my ($method) = APPRIS::Analysis::Firestar->new (
							-result					=> $report->{'result'},
							-num_residues			=> $report->{'num_residues'},
			);
			$method->residues($regions) if (defined $regions); 
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->firestar($method);
				$analysis->number($analysis->number+1);
			}			
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::Firestar->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->firestar($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

#=head2 parse_matador3d_rst
#
#  Arg [1]    : string $result
#               Parse firestar result
#  Example    : use APPRIS::Parser qw(parse_firestar);
#               parse_firestar($result);
#  Description: Parse output of firestar.
#  Returntype : APPRIS::Gene or undef
#  Exceptions : return undef
#  Caller     : generally on error
#
#=cut
#
#sub parse_matador3d_rst($)
#{
#	my ($result) = @_;
#	my ($cutoffs);
#		
#	my (@results) = split('>',$result);	
#	foreach my $transcript_result (@results)
#	{
#		#>ENST00000308249        3.65
#		#- 196:262[3]    1.33
#        #	196:242[190:240] 0.33[0.33*1*1]	 1Q33_A[35.9]
#        if ( $transcript_result =~ /^([^\t]+)\t+([^\n]+)\n+/ )
#		{
#			my ($id) = $1;
#			my ($structure_score) = $2;
#			$cutoffs->{$id}->{'score'} = $structure_score;
#
#			my ($alignment_list_report);			
#			my (@trans_alignments) = split('- ', $transcript_result);
#
#			for (my $i = 1; $i < scalar(@trans_alignments); $i++) { # jump the first line (ENST00000308249        3.65)
#				if ( $trans_alignments[$i] =~ /^(\d+)\:(\d+)\[(\d+)\]\t+([^\n]*)\n*([^\$]*)$/ )
#				{				
#					my ($cds_start) = $1;
#					my ($cds_end) = $2;
#					my ($cds_order) = $3;
#					my ($cds_score) = $4;
#					my ($trans_mini_alignments) = $5;
#
#					if(defined $cds_start and defined $cds_end and defined $cds_order and defined $cds_score and defined $trans_mini_alignments)
#					{
#						my ($alignment_report) = {
#							'cds_id'	=> $cds_order,
#							'start'		=> $cds_start,
#							'end'		=> $cds_end,
#							'score'		=> $cds_score,
#							'type'		=> 'exon',
#						};
#						push(@{$alignment_list_report}, $alignment_report);
#						
#						# get mini-alignments
#						while ( $trans_mini_alignments =~ /^\s+([^\n]*)\n*/mg ) # per line
#						{
#							if ( $1 =~ /^(\d+)\:(\d+)\[(\d+)\:(\d+)\]\t+([^\[]+)\[([^\]]+)\]\t+([^\$]*)$/ ) { #	196:242[190:240] 0.33[0.33*1*1]	 1Q33_A[35.9]
#								my ($mini_cds_start) = $1;
#								my ($mini_cds_end) = $2;
#								my ($align_start) = $3;
#								my ($align_end) = $4;								
#								my ($mini_cds_score) = $5;
#								my ($mini_cds_info) = $6;
#								my ($mini_pdb_list) = $7;
#								my ($mini_alignment_report) = {
#									'alignment_start'		=> $align_start,
#									'alignment_end'		=> $align_end,									
#									'cds_id'	=> $cds_order,
#									'start'		=> $mini_cds_start,
#									'end'		=> $mini_cds_end,
#									'score'		=> $mini_cds_score,
#									'info'		=> $mini_cds_info,
#									'type'		=> 'mini-exon',								
#								};
#								my ($mini_pdb_ident);
#								my ($mini_pdb);
#								my ($mini_ident);
#								my ($mini_ext_id);
#
#								if ( $mini_pdb_list =~ /^([^\t]+)\t+([^\$]+)$/ ) {
#									$mini_pdb_ident = $1;
#									$mini_ext_id = $2;
#								}
#								else {
#									$mini_pdb_ident = $mini_pdb_list;
#								}
#								if ( defined $mini_pdb_ident and ($mini_pdb_ident =~ /^([^\[]+)\[([^\]]+)\]$/) ) {
#									$mini_pdb = $1;
#									$mini_ident = $2;	
#								}									
#								$mini_alignment_report->{'pdb_id'} = $mini_pdb if (defined $mini_pdb);
#								$mini_alignment_report->{'identity'} = $mini_ident if (defined $mini_ident);
#								$mini_alignment_report->{'external_id'} = $mini_ext_id if (defined $mini_ext_id);								
#								push(@{$alignment_list_report}, $mini_alignment_report);							
#							}
#						}						
#					}
#				}
#			}
#			if ( defined $alignment_list_report and (scalar(@{$alignment_list_report}) > 0) )
#			{
#				$cutoffs->{$id}->{'alignments'} = $alignment_list_report;
#			}
#			$transcript_result =~ s/\n*#[^#]+#\n+#[^#]+#\n+#[^#]+#//mg;
#			$cutoffs->{$id}->{'result'}='>'.$transcript_result;			
#		}
#	}
#	$cutoffs->{'result'} = $result;
#	
#	return $cutoffs;
#}
#
#=head2 parse_matador3d
#
#  Arg [1]    : APPRIS::Gene $gene 
#               APPRIS::Gene object
#  Arg [2]    : string $result
#               Parse matador3d result
#  Example    : use APPRIS::Parser qw(parse_matador3d);
#               parse_matador3d($result);
#  Description: Parse output of matador3d.
#  Returntype : APPRIS::Gene or undef
#  Exceptions : return undef
#  Caller     : generally on error
#
#=cut
#
#sub parse_matador3d($$)
#{
#	my ($gene, $result) = @_;
#
#	my ($stable_id) = $gene->stable_id;
#	my ($transcripts);
#	my ($index_transcripts);
#	my ($index) = 0;
#	
#	# Create hash object from result
#	my ($cutoffs) = parse_matador3d_rst($result);
#	
#	# Create APPRIS object
#	foreach my $transcript (@{$gene->transcripts}) {			
#		my ($transcript_id) = $transcript->stable_id;
#		my ($transcript_ver);
#		my ($transc_eid) = $transcript_id;
#		if ( $transcript->version ) {
#			$transcript_ver = $transcript->version;			
#			$transc_eid = $transcript_id.'.'.$transcript_ver;
#		}
#		my ($analysis);
#		
#		# create method object
#		if ( exists $cutoffs->{$transcript_id} ) {
#			my ($report) = $cutoffs->{$transcript_id};
#			my ($regions);			
#			if ( $transcript->translate ) {
#				my ($translate) = $transcript->translate;
#				
#				if ( exists $report->{'alignments'} ) {
#					my ($strand) = $transcript->strand;
#					foreach my $residue (@{$report->{'alignments'}}) {
#						my ($region);
#						if ( $transcript->translate->cds ) { # with CDS coords from GTF files
#							my ($pro_coord_start) = get_coords_from_residue($transcript, $residue->{'start'});
#							my ($pro_coord_end) = get_coords_from_residue($transcript, $residue->{'end'});
#							$residue->{'trans_strand'} = $strand;
#							if ( $strand eq '-' ) {
#								$residue->{'trans_end'} = $pro_coord_start->{'start'};                                                
#								$residue->{'trans_start'} = $pro_coord_end->{'end'};                                              
#							}
#							else {
#								$residue->{'trans_start'} = $pro_coord_start->{'start'};                                                
#								$residue->{'trans_end'} = $pro_coord_end->{'end'};                                              
#							}
#							$region = APPRIS::Analysis::Matador3DRegion->new (
#											-cds_id		=> $residue->{'cds_id'},
#											-pstart		=> $residue->{'start'},
#											-pend		=> $residue->{'end'},
#											-score		=> $residue->{'score'},
#											-start		=> $residue->{'trans_start'},
#											-end		=> $residue->{'trans_end'},
#											-strand		=> $residue->{'trans_strand'},										
#							);
#							$region->type($residue->{'type'}) if (exists $residue->{'type'} and defined $residue->{'type'});
#							$region->alignment_start($residue->{'alignment_start'}) if (exists $residue->{'alignment_start'} and defined $residue->{'alignment_start'});
#							$region->alignment_end($residue->{'alignment_end'}) if (exists $residue->{'alignment_end'} and defined $residue->{'alignment_end'});
#							$region->pdb_id($residue->{'pdb_id'}) if (exists $residue->{'pdb_id'} and defined $residue->{'pdb_id'});
#							$region->identity($residue->{'identity'}) if (exists $residue->{'identity'} and defined $residue->{'identity'});
#							$region->external_id($residue->{'external_id'}) if (exists $residue->{'external_id'} and defined $residue->{'external_id'});
#						}
#						else {
#							$region = APPRIS::Analysis::Matador3DRegion->new (
#											-cds_id		=> $residue->{'cds_id'},
#											-pstart		=> $residue->{'start'},
#											-pend		=> $residue->{'end'},
#											-score		=> $residue->{'score'},
#							);
#							$region->type($residue->{'type'}) if (exists $residue->{'type'} and defined $residue->{'type'});
#							$region->alignment_start($residue->{'alignment_start'}) if (exists $residue->{'alignment_start'} and defined $residue->{'alignment_start'});
#							$region->alignment_end($residue->{'alignment_end'}) if (exists $residue->{'alignment_end'} and defined $residue->{'alignment_end'});
#							$region->pdb_id($residue->{'pdb_id'}) if (exists $residue->{'pdb_id'} and defined $residue->{'pdb_id'});
#							$region->identity($residue->{'identity'}) if (exists $residue->{'identity'} and defined $residue->{'identity'});
#							$region->external_id($residue->{'external_id'}) if (exists $residue->{'external_id'} and defined $residue->{'external_id'});
#						}
#						push(@{$regions}, $region) if ( defined $region );
#					}
#				}
#			}
#			
#			# create Analysis object (for trans)			
#			my ($method) = APPRIS::Analysis::Matador3D->new (
#							-result					=> $report->{'result'},
#							-score					=> $report->{'score'}
#			);
#			if (defined $regions and (scalar(@{$regions}) > 0) ) {
#				$method->alignments($regions);
#				$method->num_alignments(scalar(@{$regions}));
#			}			
#			$analysis = APPRIS::Analysis->new();
#			if (defined $method) {
#				$analysis->matador3d($method);
#				$analysis->number($analysis->number+1);
#			}			
#		}
#				
#		# create Transcript object
#		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
#		$transcript->version($transcript_ver) if (defined $transcript_ver);
#		$transcript->analysis($analysis) if (defined $analysis);
#		push(@{$transcripts}, $transcript);
#		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
#	}
#
#	# create Analysis object (for gene)
#	my ($method2) = APPRIS::Analysis::Matador3D->new( -result => $cutoffs->{'result'} );	
#	my ($analysis2) = APPRIS::Analysis->new();
#	if (defined $method2) {
#		$analysis2->matador3d($method2);
#		$analysis2->number($analysis2->number+1);
#	}
#	
#	# create Gene object
#	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
#	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
#	$entity->analysis($analysis2) if (defined $analysis2);	
#
#	return $entity;
#}


=head2 parse_matador3d_rst

  Arg [1]    : string $result
               Parse firestar result
  Example    : use APPRIS::Parser qw(parse_firestar);
               parse_firestar($result);
  Description: Parse output of firestar.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_matador3d_rst($)
{
	my ($result) = @_;
	my ($cutoffs);
		
	my (@results) = split("\n",$result);	
	foreach my $transcript_result (@results)
	{
		#ENSG00000007216	ENST00000459818	0
		#ENSG00000007216	ENST00000314669	258.8	4f35_D;124.3;13.9;405-552;148	4f35_D;102.8;10.0;45-162;118	4r1i_A;31.7;6.6;220-345;126
        if ( $transcript_result =~ /^([^\s]*)\s+([^\s]+)\s+([^\$]+)$/ )
		{
			my ($gid) = $1; # can be undefined
			my ($id) = $2;
			my (@trans_score_alignments) = split("\t", $3);
			my ($structure_score) = $trans_score_alignments[0];
			my (@trans_alignments) = ( defined $trans_score_alignments[1] ) ?  @trans_score_alignments[1..$#trans_score_alignments] : undef;
			$cutoffs->{$id}->{'score'} = $structure_score;

			my ($alignment_list_report);
			for (my $i = 0; $i < scalar(@trans_alignments); $i++) { # 4f35_D;124.3;13.9;405-552;148
				if ( defined $trans_alignments[$i] and $trans_alignments[$i] =~ /^([^\;]+)\;+([^\;]+)\;+([^\;]+)\;+([^\-]*)\-([^\;]+)\;+([^\$]*)$/ )
				{
					my ($trans_align_pdb_id) = $1;
					my ($trans_align_score) = $2;
					my ($trans_align_bias) = $3;
					my ($trans_align_start) = $4;
					my ($trans_align_end) = $5;
					my ($trans_align_abs_pos) = $6;

					if( defined $trans_align_pdb_id and defined $trans_align_score and defined $trans_align_bias and 
						defined $trans_align_start and defined $trans_align_end and defined $trans_align_abs_pos
					){
						my ($alignment_report) = {
							'start'		=> $trans_align_start,
							'end'		=> $trans_align_end,
							'score'		=> $trans_align_score,
							'bias'		=> $trans_align_bias,
							'pdb_id'	=> $trans_align_pdb_id,
						};
						push(@{$alignment_list_report}, $alignment_report);					
					}
				}
			}			
			if ( defined $alignment_list_report and (scalar(@{$alignment_list_report}) > 0) )
			{
				$cutoffs->{$id}->{'alignments'} = $alignment_list_report;
			}
			$transcript_result =~ s/\n*#[^#]+#\n+#[^#]+#\n+#[^#]+#//mg;
			$cutoffs->{$id}->{'result'}=$transcript_result;
		}
	}
	$cutoffs->{'result'} = $result;
	
	return $cutoffs;
}

=head2 parse_matador3d

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse matador3d result
  Example    : use APPRIS::Parser qw(parse_matador3d);
               parse_matador3d($result);
  Description: Parse output of matador3d.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_matador3d($$)
{
	my ($gene, $result) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;
	
	# Create hash object from result
	my ($cutoffs) = parse_matador3d_rst($result);
	
	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
		
		# create method object
		if ( exists $cutoffs->{$transcript_id} ) {
			my ($report) = $cutoffs->{$transcript_id};
			my ($regions);			
			if ( $transcript->translate ) {
				my ($translate) = $transcript->translate;
				
				if ( exists $report->{'alignments'} ) {
					my ($strand) = $transcript->strand;
					foreach my $residue (@{$report->{'alignments'}}) {
						my ($region);
						if ( $transcript->translate->cds ) { # with CDS coords from GTF files
							my ($pro_coord_start) = get_coords_from_residue($transcript, $residue->{'start'});
							my ($pro_coord_end) = get_coords_from_residue($transcript, $residue->{'end'});
							$residue->{'trans_strand'} = $strand;
							if ( $strand eq '-' ) {
								$residue->{'trans_end'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_start'} = $pro_coord_end->{'end'};                                              
							}
							else {
								$residue->{'trans_start'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_end'} = $pro_coord_end->{'end'};                                              
							}
							$region = APPRIS::Analysis::Matador3DRegion->new (
											-pstart		=> $residue->{'start'},
											-pend		=> $residue->{'end'},
											-score		=> $residue->{'score'},
											-bias		=> $residue->{'bias'},
											-pdb_id		=> $residue->{'pdb_id'},
											-start		=> $residue->{'trans_start'},
											-end		=> $residue->{'trans_end'},
											-strand		=> $residue->{'trans_strand'}								
							);
						}
						else {
							$region = APPRIS::Analysis::Matador3DRegion->new (
											-pstart		=> $residue->{'start'},
											-pend		=> $residue->{'end'},
											-score		=> $residue->{'score'},
											-bias		=> $residue->{'bias'},
											-pdb_id		=> $residue->{'pdb_id'}											
							);
						}
						push(@{$regions}, $region) if ( defined $region );
					}
				}
			}
			
			# create Analysis object (for trans)			
			my ($method) = APPRIS::Analysis::Matador3D->new (
							-result					=> $report->{'result'},
							-score					=> $report->{'score'}
			);
			if (defined $regions and (scalar(@{$regions}) > 0) ) {
				$method->alignments($regions);
				$method->num_alignments(scalar(@{$regions}));
			}			
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->matador3d($method);
				$analysis->number($analysis->number+1);
			}			
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::Matador3D->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->matador3d($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}


=head2 parse_spade_rst

  Arg [1]    : string $result
               Parse spade result
  Example    : use APPRIS::Parser qw(parse_spade);
               parse_firestar($result);
  Description: Parse output of firestar.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_spade_rst($)
{
	my ($result) = @_;
	my ($cutoffs);
	
	my (@results) = split('>',$result);	
	foreach my $transcript_result (@results)
	{
		#>ENST00000356093        571.9   4       2       0       0
		#domain  5       81      5       81      PF09379.3       FERM_N  Domain  1       80      80      74.3    4.7e-21 1       CL0072  [ext:ENST00000373800]
		#domain_possibly_damaged 83      192     83      192     PF00373.11      FERM_M  Domain  1       117     117     74.7    4.9e-21 1       No_clan [ext:ENST00000373800]
		#domain_possibly_damaged 198     281     196     285     PF09380.3       FERM_C  Domain  3       85      90      81.9    2.2e-23 1       CL0266  [ext:ENST00000373800]
		#domain  289     333     289     335     PF08736.4       FA      Family  1       45      47      66.0    1.5e-18 1       No_clan [ext:ENST00000373800]
		#domain  425     473     425     473     PF04382.6       SAB     Domain  1       48      48      93.1    4.6e-27 1       No_clan [discarded]
		#domain  506     619     505     619     PF05902.6       4_1_CTD Domain  2       114     114     181.9   2.3e-54 1       No_clan
        if ( $transcript_result=~/^([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\n]+)\n+/ )
		{
			my ($id) = $1;
			my ($bitscore) = $2;
			my ($num_domains) = $3;
			my ($num_possibly_damaged_domains) = $4;
			my ($num_damaged_domains) = $5;
			my ($num_wrong_domains) = $6;
			$cutoffs->{$id}->{'bitscore'} = $bitscore;
			$cutoffs->{$id}->{'num_domains'} = $num_domains;
			$cutoffs->{$id}->{'num_possibly_damaged_domains'} = $num_possibly_damaged_domains;
			$cutoffs->{$id}->{'num_damaged_domains'} = $num_damaged_domains;
			$cutoffs->{$id}->{'num_wrong_domains'} = $num_wrong_domains;
	
			# <type_domain>
			# <alignment start> <alignment end> <envelope start> <envelope end>
			# <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value>
			# <significance> <clan> <predicted_active_site_residues> -optional values-
			# <external id> || <discarded> -optional values-
			my ($alignment_list);
			while ( $transcript_result =~ /([domain|domain_possibly_damaged|domain_damaged|domain_wrong][^\n]*)\n*/mg )
			{
				my ($value_line) = $1;
				my (@value_list) = split(/\t+/,$value_line);
				if ( scalar(@value_list) > 12 )
				{
					my ($alignment_report);
					
					my ($type_domain) = $value_list[0];
					my ($alignment_start) = $value_list[1];
					my ($alignment_end) = $value_list[2];
					my ($envelope_start) = $value_list[3];
					my ($envelope_end) = $value_list[4];
					my ($hmm_acc) = $value_list[5];
					my ($hmm_name) = $value_list[6];
					my ($hmm_type) = $value_list[7];
					my ($hmm_start) = $value_list[8];
					my ($hmm_end) = $value_list[9];
					my ($hmm_length) = $value_list[10];
					my ($bit_score) = $value_list[11];
					my ($e_value) = $value_list[12];
					
					# required values
					$alignment_report->{'type_domain'} = $type_domain;
					$alignment_report->{'alignment_start'} = $alignment_start;
					$alignment_report->{'alignment_end'} = $alignment_end;
					$alignment_report->{'envelope_start'} = $envelope_start;
					$alignment_report->{'envelope_end'} = $envelope_end;
					$alignment_report->{'hmm_acc'} = $hmm_acc;
					$alignment_report->{'hmm_name'} = $hmm_name;
					$alignment_report->{'hmm_type'} = $hmm_type;
					$alignment_report->{'hmm_start'} = $hmm_start;
					$alignment_report->{'hmm_end'} = $hmm_end;
					$alignment_report->{'hmm_length'} = $hmm_length;
					$alignment_report->{'bit_score'} = $bit_score;
					$alignment_report->{'e_value'} = $e_value;
					
					# optional values taking into account the value of external_id
					if(defined $value_list[14] and !($value_list[14] =~ /\[[^\]]*\]/)) {
						$alignment_report->{'significance'} = $value_list[14];
					}
					elsif(defined $value_list[14] and ($value_list[14] =~ /\[([^\]]*)\]/)) {
						my ($val) = $1;
						if ( $val =~ /discarded/ ) { $alignment_report->{'discarded'} = 1 }
						elsif ( $val =~ /ext\:([^\$]*)$/ ) { $alignment_report->{'external_id'} = $1 }
					}
					
					if(defined $value_list[15] and !($value_list[15] =~ /\[[^\]]*\]/)) {
						$alignment_report->{'clan'} = $value_list[15];
					}
					elsif(defined $value_list[15] and ($value_list[15] =~ /\[([^\]]*)\]/)) {
						my ($val) = $1;
						if ( $val =~ /discarded/ ) { $alignment_report->{'discarded'} = 1 }
						elsif ( $val =~ /ext\:([^\$]*)$/ ) { $alignment_report->{'external_id'} = $1 }
					}
									
					if(defined $value_list[16] and !($value_list[16] =~ /\[[^\]]*\]/)) {
						$alignment_report->{'predicted_active_site_residues'} = $value_list[16];
					}
					elsif(defined $value_list[16] and ($value_list[16] =~ /\[([^\]]*)\]/)) {
						my ($val) = $1;
						if ( $val =~ /discarded/ ) { $alignment_report->{'discarded'} = 1 }
						elsif ( $val =~ /ext\:([^\$]*)$/ ) { $alignment_report->{'external_id'} = $1 }
					}

					if(defined $value_list[17] and ($value_list[17] =~ /\[([^\]]*)\]/)) {
						my ($val) = $1;
						if ( $val =~ /discarded/ ) { $alignment_report->{'discarded'} = 1 }
						elsif ( $val =~ /ext\:([^\$]*)$/ ) { $alignment_report->{'external_id'} = $1 }
					}

					push(@{$alignment_list}, $alignment_report);					
				}
			}
			if ( defined $alignment_list and (scalar(@{$alignment_list}) > 0) )
			{
				$cutoffs->{$id}->{'domains'} = $alignment_list;
			}
			$transcript_result =~ s/\n*#[^#]+#\n+#[^#]+#\n+#[^#]+#//mg;
			$cutoffs->{$id}->{'result'} = '>'.$transcript_result;			
		}
        elsif ( $transcript_result =~ /^([^\s]+)(\s+[^\n]*\n+#HMM[^\n]*\n+#MATCH[^\n]*\n+#PP[^\n]*\n+#SEQ[^\n]*\n+)/ )
		{
			#>ENST00000270190     10    144     10    145 PF04118.7   Dopey_N           Family     1   136   309    199.8   3.3e-59   1 No_clan  
			##HMM       kdskqkkyasevekaLksFetlqEWADyisfLskLlkalqkkqeklsyvpskllvskrLaqcLnpsLPsGVHqkaLevYelIfekigketLskdlalylsGlfpllsyasisvkplllellekyllpLekalrpll
			##MATCH     +d++++ y+s +ekaL++Fe+++EWAD+is+L+kL+kalq ++ ++s +p++ll+skrLaqcL+p+LPsGVH kaLe+Ye+If+++g+++L+kdl+ly+ Glfpll++a++sv+p+ll+l+eky+lpL+k l+p l
			##PP        5899************************************.*****************************************************************************************999977
			##SEQ       NDYRYRSYSSVIEKALRNFESSSEWADLISSLGKLNKALQ-SNLRYSLLPRRLLISKRLAQCLHPALPSGVHLKALETYEIIFKIVGTKWLAKDLFLYSCGLFPLLAHAAVSVRPVLLTLYEKYFLPLQKLLLPSL

			my ($id) = $1;
			$transcript_result =~ s/\n*#[^#]+#\n+#[^#]+#\n+#[^#]+#//mg;						
			$cutoffs->{$id}->{'result'} .= '>'.$transcript_result;
		}	
	}
	$cutoffs->{'result'} = $result;
	
	return $cutoffs; 
}

=head2 parse_spade

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse spade result
  Example    : use APPRIS::Parser qw(parse_spade);
               parse_spade($result);
  Description: Parse output of spade.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_spade($$)
{
	my ($gene, $result) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;

	# Create hash object from result
	my ($cutoffs) = parse_spade_rst($result);

	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
		
		# create method object
		if ( exists $cutoffs->{$transcript_id} ) {
			my ($report) = $cutoffs->{$transcript_id};
			my ($regions);
			if ( $transcript->translate ) {				
				if ( exists $report->{'domains'} ) {
					my ($strand) = $transcript->strand;
					foreach my $residue (@{$report->{'domains'}}) {
						my ($region);
						if ( $transcript->translate->cds ) { # with CDS coords from GTF file
							my ($pro_coord_start) = get_coords_from_residue($transcript, $residue->{'alignment_start'});
							my ($pro_coord_end) = get_coords_from_residue($transcript, $residue->{'alignment_end'});
							$residue->{'trans_strand'} = $strand;
							if ( $strand eq '-' ) {
								$residue->{'trans_end'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_start'} = $pro_coord_end->{'end'};                                              
							}
							else {
								$residue->{'trans_start'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_end'} = $pro_coord_end->{'end'};                                              
							}
							$region = APPRIS::Analysis::SPADERegion->new (
											-start								=> $residue->{'trans_start'},
											-end								=> $residue->{'trans_end'},
											-strand								=> $residue->{'trans_strand'},						
											-type_domain						=> $residue->{'type_domain'},						
											-alignment_start					=> $residue->{'alignment_start'},					
											-alignment_end						=> $residue->{'alignment_end'},
											-envelope_start						=> $residue->{'envelope_start'},					
											-envelope_end						=> $residue->{'envelope_end'},
											-hmm_start							=> $residue->{'hmm_start'},					
											-hmm_end							=> $residue->{'hmm_end'},
											-hmm_length							=> $residue->{'hmm_length'},					
											-hmm_acc							=> $residue->{'hmm_acc'},
											-hmm_name							=> $residue->{'hmm_name'},					
											-hmm_type							=> $residue->{'hmm_type'},
											-bit_score							=> $residue->{'bit_score'},
											-evalue								=> $residue->{'e_value'}
							);
							$region->significance($residue->{'significance'}) if (exists $residue->{'significance'} and defined $residue->{'significance'});
							$region->clan($residue->{'clan'}) if (exists $residue->{'clan'} and defined $residue->{'clan'});
							$region->predicted_active_site_residues($residue->{'predicted_active_site_residues'}) if (exists $residue->{'predicted_active_site_residues'} and defined $residue->{'predicted_active_site_residues'});
							$region->external_id($residue->{'external_id'}) if (exists $residue->{'external_id'} and defined $residue->{'external_id'});
							$region->discarded($residue->{'discarded'}) if (exists $residue->{'discarded'} and defined $residue->{'discarded'});
						}
						else {
							$region = APPRIS::Analysis::SPADERegion->new (
											-type_domain						=> $residue->{'type_domain'},						
											-alignment_start					=> $residue->{'alignment_start'},					
											-alignment_end						=> $residue->{'alignment_end'},
											-envelope_start						=> $residue->{'envelope_start'},					
											-envelope_end						=> $residue->{'envelope_end'},
											-hmm_start							=> $residue->{'hmm_start'},					
											-hmm_end							=> $residue->{'hmm_end'},
											-hmm_length							=> $residue->{'hmm_length'},					
											-hmm_acc							=> $residue->{'hmm_acc'},
											-hmm_name							=> $residue->{'hmm_name'},					
											-hmm_type							=> $residue->{'hmm_type'},
											-bit_score							=> $residue->{'bit_score'},
											-evalue								=> $residue->{'e_value'}
							);
							$region->significance($residue->{'significance'}) if (exists $residue->{'significance'} and defined $residue->{'significance'});
							$region->clan($residue->{'clan'}) if (exists $residue->{'clan'} and defined $residue->{'clan'});
							$region->predicted_active_site_residues($residue->{'predicted_active_site_residues'}) if (exists $residue->{'predicted_active_site_residues'} and defined $residue->{'predicted_active_site_residues'});
							$region->external_id($residue->{'external_id'}) if (exists $residue->{'external_id'} and defined $residue->{'external_id'});
							$region->discarded($residue->{'discarded'}) if (exists $residue->{'discarded'} and defined $residue->{'discarded'});
						}
						push(@{$regions}, $region);
					}
				}
			}
			
			# create Analysis object (for trans)			
			my ($method) = APPRIS::Analysis::SPADE->new (
							-result							=> $report->{'result'},
							-num_domains					=> $report->{'num_domains'},
							-num_possibly_damaged_domains	=> $report->{'num_possibly_damaged_domains'},
							-num_damaged_domains			=> $report->{'num_damaged_domains'},
							-num_wrong_domains				=> $report->{'num_wrong_domains'},
							-bitscore						=> $report->{'bitscore'}
			);
			$method->regions($regions) if (defined $regions and (scalar(@{$regions}) > 0) );
			
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->spade($method);
				$analysis->number($analysis->number+1);
			}			
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);		
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::SPADE->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->spade($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

=head2 parse_inertia

  Arg [1]    : String $id 
               The stable ID of the gene to retrieve
  Arg [2]    : string $inertia
               INERTIA result
  Arg [3]    : string $mafft
               MAFFT Omega result
  Arg [4]    : string $prank
               Prank Omega result
  Arg [5]    : string $kalign
               Kalign Omega result
  Example    : use APPRIS::Parser qw(parse_inertia);
               parse_inertia($id, $mafft, $prank, $kalign);
  Description: Parse the output of inertia.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

#sub parse_inertia($$;$;$;$;$)
#{
#	my ($gene, $inertia_i, $mafft_i, $prank_i, $kalign_i, $compara_i) = @_;
sub parse_inertia($$)
{
	my ($gene, $inertia_i) = @_;
	
	my ($stable_id) = $gene->stable_id;	
	my ($transcripts);
	my ($cutoffs);
	my ($index_transcripts);
	my ($index) = 0;	
	
	_parse_inertia_file('inertia', $inertia_i, $cutoffs) if (defined $inertia_i);
#	_parse_omega_file('mafft', $mafft_i, $cutoffs) if (defined $mafft_i);
#	_parse_omega_file('prank', $prank_i, $cutoffs) if (defined $prank_i);
#	_parse_omega_file('kalign', $kalign_i, $cutoffs) if (defined $kalign_i);
#	_parse_omega_file('compara', $compara_i, $cutoffs) if (defined $compara_i);

	# Create APPRIS objects
	foreach my $transcript (@{$gene->transcripts}) {
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);

		# create method object
		if ( exists $cutoffs->{$transcript_id} and exists $cutoffs->{$transcript_id}->{'inertia'} ) {
			my ($report) = $cutoffs->{$transcript_id};
			my ($report2) = $report->{'inertia'};
			my ($method);
			
			# create inertia object
			my ($regions);			
			foreach my $residue (@{$report2->{'residues'}})
			{
				push(@{$regions},
					APPRIS::Analysis::INERTIARegion->new
					(
						-start				=> $residue->{'start'},
						-end				=> $residue->{'end'},
						-strand				=> $residue->{'strand'},
						-unusual_evolution	=> $residue->{'unusual_evolution'}
					)
				);
			}
			if ( exists $report2->{'unusual_evolution'} and defined $report2->{'unusual_evolution'} ) {
				$method = APPRIS::Analysis::INERTIA->new
				(
					-result					=> $report2->{'result'},
					-unusual_evolution		=> $report2->{'unusual_evolution'}
				);
				$method->regions($regions) if (defined $regions);				
			}
			# create omega objects
			while ( my ($type, $report2) = each(%{$report}) )
			{
				if ( ($type eq 'mafft') or ($type eq 'prank') or ($type eq 'kalign') or ($type eq 'compara') )
				{					
					my ($omega);
					my ($regions);
					
					foreach my $residue (@{$report2->{'residues'}})
					{
						push(@{$regions},
							APPRIS::Analysis::OmegaRegion->new
							(
								-start				=> $residue->{'start'},
								-end				=> $residue->{'end'},
								-omega_mean			=> $residue->{'omega_mean'},
								-st_deviation		=> $residue->{'st_deviation'},
								-p_value			=> $residue->{'p_value'},
								-difference_value	=> $residue->{'difference_value'},
								-unusual_evolution	=> $residue->{'unusual_evolution'}
							)
						);			
					}
					$omega = APPRIS::Analysis::Omega->new
					(
						-average			=> $report2->{'omega_average'},
						-st_desviation		=> $report2->{'omega_st_desviation'},
						-result				=> $report2->{'result'},
						-unusual_evolution	=> $report2->{'unusual_evolution'}
					);
					$omega->regions($regions) if (defined $regions);
					
					$method->mafft_alignment($omega) if ($method and $omega and ($type eq 'mafft') );
					$method->prank_alignment($omega) if ($method and $omega and ($type eq 'prank') );
					$method->kalign_alignment($omega) if ($method and $omega and ($type eq 'kalign') );
					$method->compara_alignment($omega) if ($method and $omega and ($type eq 'compara') );
				}			
			}
			
			# create Analysis object (for trans)
			$analysis = APPRIS::Analysis->new();			
			if (defined $method) {				
				$analysis->inertia($method);
				$analysis->number($analysis->number+1);
			}
		}
					
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);		
		$transcript->analysis($analysis) if (defined $analysis);
		if ( defined $transcript ) {
			push(@{$transcripts}, $transcript);
			$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
		}
	}
	
	# create Analysis object (for gene)
	my ($analysis2) = APPRIS::Analysis->new();
	if ( defined $inertia_i and ($inertia_i ne '') ) {
		my ($method2) = APPRIS::Analysis::INERTIA->new( -result => $inertia_i );	
		if (defined $method2) {
#			if ( defined $mafft_i and ($mafft_i ne '') ) {	
#				my ($method3) = APPRIS::Analysis::Omega->new( -result => $mafft_i );	
#				if (defined $method3) {
#					$method2->mafft_alignment($method3);
#				}
#			}		
#			if ( defined $prank_i and ($prank_i ne '') ) {	
#				my ($method4) = APPRIS::Analysis::Omega->new( -result => $prank_i );	
#				if (defined $method4) {
#					$method2->prank_alignment($method4);
#				}
#			}		
#			if ( defined $kalign_i and ($kalign_i ne '') ) {	
#				my ($method5) = APPRIS::Analysis::Omega->new( -result => $kalign_i );	
#				if (defined $method5) {
#					$method2->kalign_alignment($method5);
#				}
#			}
#			if ( defined $compara_i and ($compara_i ne '') ) {	
#				my ($method6) = APPRIS::Analysis::Omega->new( -result => $compara_i );	
#				if (defined $method6) {
#					$method2->compara_alignment($method6);
#				}
#			}
			$analysis2->inertia($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

=head2 parse_corsair_rst

  Arg [1]    : string $result
               Parse corsair result
  Example    : use APPRIS::Parser qw(parse_corsair);
               parse_firestar($result);
  Description: Parse output of firestar.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_corsair_rst($)
{
	my ($result) = @_;
	my ($cutoffs);
	
	my (@results) = split('>',$result);	
	foreach my $transcript_result (@results)
	{
		#>ENST00000518498.1      0.5
		#Homo sapiens    100.00  0.5
		#        -43735403:43735526[1:42]        0.5
		#                -Homo sapiens   100.00  0.5
		#        -43733595:43733741[43:91]       0.5 {1-ENST00000291525.10}
		#                -Homo sapiens   100.00  0.5
		#        -43732369:43732379[92:94]       0.5 {1-ENST00000291525.10}
		#                -Homo sapiens   100.00  0.5
        if ( $transcript_result =~ /^([^\t]+)\t+([^\n]+)\n+/ )
		{
			my ($id) = $1;
			my ($score) = $2;
			$cutoffs->{$id}->{'score'} = $score;
									
			my ($alignment_list_report);			
			my (@trans_alignments) = split('- ', $transcript_result);
			
			my ($cds_order) = 1;
			for (my $i = 1; $i < scalar(@trans_alignments); $i++) { # jump the first line #ENST00000518498.1      0.5 #Homo sapiens    100.00  0.5
		        #- 47588222:47588303[1:28]       1.5
		        #        Pan troglodytes 96.43   1
		        #        Homo sapiens    100.00  0.5
		        #- 47581830:47581981[28:78]      1.5 {6.1-ENST00000291672.5}
		        #        Pan troglodytes 98.04   1
		        #        Homo sapiens    100.00  0.5
                if ( $trans_alignments[$i] =~ /^(\d+)\:(\d+)\[(\d+)\:(\d+)\]\t+([^\n]*)\n*([^\$]*)$/ )
				{
					my ($trans_start) = $1;
					my ($trans_end) = $2;
					my ($trans_strand) = '.';
					my ($cds_start) = $3;
					my ($cds_end) = $4;
					my ($score_list) = $5;
					my ($sp_report) = $6;
					$sp_report =~ s/^\s*//g; $sp_report =~ s/\s*$//g; $sp_report =~ s/\n*#[^\$]*$//g;
					my ($cds_score);
					my ($cds_maxscore);
					if ( $score_list =~ /^([\d|\.]+)$/ ) {
						$cds_score = $1;
					}
					elsif ( $score_list =~ /^([^\s]+)\s+\{([^\-]*)\-/ ) {
						$cds_score = $1;
						$cds_maxscore = $2;
					}
					if(defined $trans_start and defined $trans_end and defined $trans_strand and defined $cds_start and defined $cds_end and defined $cds_order and defined $cds_score )
					{
						my ($alignment_report) = {
							'cds_id'		=> $cds_order,
							'trans_start'	=> $trans_start,
							'trans_end'		=> $trans_end,
							'trans_strand'	=> $trans_strand,
							'start'			=> $cds_start,
							'end'			=> $cds_end,
							'score'			=> $cds_score,
							'maxscore'		=> $cds_maxscore,
							'sp_report'		=> $sp_report,
							'type'			=> 'exon',
						};
						push(@{$alignment_list_report}, $alignment_report);
						$cds_order ++;							
					}
				}
			}
			if ( defined $alignment_list_report and (scalar(@{$alignment_list_report}) > 0) )
			{
				$cutoffs->{$id}->{'alignments'} = $alignment_list_report;
			}
			$transcript_result =~ s/\n*#[^#]+#\n+#[^#]+#\n+#[^#]+#//mg;
			$cutoffs->{$id}->{'result'}='>'.$transcript_result;	

		}
	}
	$cutoffs->{'result'} = $result;
	
	return $cutoffs;
}

=head2 parse_corsair

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse corsair result
  Example    : use APPRIS::Parser qw(parse_corsair);
               parse_corsair($result);
  Description: Parse output of corsair.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_corsair($$)
{
	my ($gene, $result) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;

	# Create hash object from result
	my ($cutoffs) = parse_corsair_rst($result);

	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
		
		# create method object
		if ( exists $cutoffs->{$transcript_id} ) {
			my ($report) = $cutoffs->{$transcript_id};			
			my ($regions);			
			if ( $transcript->translate ) {
				my ($translate) = $transcript->translate;				
				if ( exists $report->{'alignments'} ) {
					my ($strand) = $transcript->strand;
					foreach my $residue (@{$report->{'alignments'}}) {
						my ($region);
						if ( $transcript->translate->cds ) { # with CDS coords from GTF files
							$region = APPRIS::Analysis::CORSAIRRegion->new (
											-cds_id		=> $residue->{'cds_id'},
											-pstart		=> $residue->{'start'},
											-pend		=> $residue->{'end'},
											-score		=> $residue->{'score'},
											-start		=> $residue->{'trans_start'},
											-end		=> $residue->{'trans_end'},
											-strand		=> $residue->{'trans_strand'},										
							);
							$region->type($residue->{'type'}) if (exists $residue->{'type'} and defined $residue->{'type'});
							$region->maxscore($residue->{'maxscore'}) if (exists $residue->{'maxscore'} and defined $residue->{'maxscore'});
							$region->sp_report($residue->{'sp_report'}) if (exists $residue->{'sp_report'} and defined $residue->{'sp_report'});
						}
						else {
							$region = APPRIS::Analysis::CORSAIRRegion->new (
											-cds_id		=> $residue->{'cds_id'},
											-pstart		=> $residue->{'start'},
											-pend		=> $residue->{'end'},
											-score		=> $residue->{'score'},
							);
							$region->type($residue->{'type'}) if (exists $residue->{'type'} and defined $residue->{'type'});
							$region->maxscore($residue->{'maxscore'}) if (exists $residue->{'maxscore'} and defined $residue->{'maxscore'});
							$region->sp_report($residue->{'sp_report'}) if (exists $residue->{'sp_report'} and defined $residue->{'sp_report'});							
						}
						push(@{$regions}, $region);						
					}
				}
			}

			# create Analysis object (for trans)			
			my ($method) = APPRIS::Analysis::CORSAIR->new (
							-result							=> $report->{'result'},
							-score							=> $report->{'score'}							
			);
			if (defined $regions and (scalar(@{$regions}) > 0) ) {
				$method->alignments($regions);
				$method->num_alignments(scalar(@{$regions}));
			}			
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->corsair($method);
				$analysis->number($analysis->number+1);
			}			
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new
		(
			-stable_id	=> $transcript_id,
			-version	=> $transcript_ver,
		);
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::CORSAIR->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->corsair($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

=head2 parse_crash_rst

  Arg [1]    : string $result
               Parse crash result
  Example    : use APPRIS::Parser qw(parse_crash);
               parse_firestar($result);
  Description: Parse output of firestar.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_crash_rst($)
{
	my ($result) = @_;
	my ($cutoffs);
	
	my (@results) = split('>',$result);	
	foreach my $transcript_result (@results)
	{
        #id      start   end     s_mean  d_score c_max   s_prob  sp_score        peptide_signal  localization    reliability     tp_score        mitochondrial_signal
		#>ENST00000479548        1       20      0.902   0.832   0.963   0.709   2       YES     S       2       -3      NO
        if ( $transcript_result =~ /^([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+($APPRIS::Utils::Constant::OK_LABEL|$APPRIS::Utils::Constant::NO_LABEL|$APPRIS::Utils::Constant::UNKNOWN_LABEL)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+($APPRIS::Utils::Constant::OK_LABEL|$APPRIS::Utils::Constant::NO_LABEL|$APPRIS::Utils::Constant::UNKNOWN_LABEL)\n+/ )        
		{
			my ($id) = $1;
			my ($pstart) = $2;
			my ($pend) = $3;
			my ($s_mean) = $4;
			my ($d_score) = $5;
			my ($c_max) = $6;
			my ($s_prob) = $7;
			my ($sp_score) = $8;
			my ($peptide_signal) = $9;
			my ($localization) = $10;
			my ($reliability) = $11;
			my ($tp_score) = $12;
			my ($mitochondrial_signal) = $13;

			$cutoffs->{$id}->{'result'} = '>'.$transcript_result;
			$cutoffs->{$id}->{'sp_score'} = $sp_score;
			$cutoffs->{$id}->{'tp_score'} = $tp_score;
			$cutoffs->{$id}->{'peptide_signal'} = $peptide_signal;
			$cutoffs->{$id}->{'mitochondrial_signal'} = $mitochondrial_signal;
			$cutoffs->{$id}->{'pstart'} = $pstart;
			$cutoffs->{$id}->{'pend'} = $pend;
			$cutoffs->{$id}->{'s_mean'} = $s_mean;
			$cutoffs->{$id}->{'s_prob'} = $s_prob;
			$cutoffs->{$id}->{'d_score'} = $d_score;
			$cutoffs->{$id}->{'c_max'} = $c_max;
			$cutoffs->{$id}->{'reliability'} = $reliability;
			$cutoffs->{$id}->{'localization'} = $localization;
		}
	}
	$cutoffs->{'result'} = $result;
	
	return $cutoffs;
}

=head2 parse_crash

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse crash result
  Example    : use APPRIS::Parser qw(parse_crash);
               parse_crash($result);
  Description: Parse output of crash.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_crash($$)
{
	my ($gene, $result) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;

	# Create hash object from result
	my ($cutoffs) = parse_crash_rst($result);

	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
		
		# create method object
		if ( exists $cutoffs->{$transcript_id} and exists $cutoffs->{$transcript_id}->{'peptide_signal'} and exists $cutoffs->{$transcript_id}->{'mitochondrial_signal'} ) {
			my ($report) = $cutoffs->{$transcript_id};
			my ($regions);
			if ( $transcript->translate ) {						
				if ( exists $report->{'pstart'} and exists $report->{'pend'} and ($report->{'pstart'} ne '-') and ($report->{'pend'} ne '-') ) {
					my ($region);
					if ( $transcript->translate->cds ) { # with CDS coords from GTF file
						my ($strand) = $transcript->strand;
						my ($pro_coord_start) = get_coords_from_residue($transcript, $report->{'pstart'});
						my ($pro_coord_end) = get_coords_from_residue($transcript, $report->{'pend'});
						$report->{'trans_strand'} = $strand;
						if ( $strand eq '-' ) {
							$report->{'trans_end'} = $pro_coord_start->{'start'};                                                
							$report->{'trans_start'} = $pro_coord_end->{'end'};                                              
						}
						else {
							$report->{'trans_start'} = $pro_coord_start->{'start'};                                                
							$report->{'trans_end'} = $pro_coord_end->{'end'};                                              
						}						
						$region = APPRIS::Analysis::CRASHRegion->new (
										-start						=> $report->{'trans_start'},
										-end						=> $report->{'trans_end'},
										-strand						=> $report->{'trans_strand'},
										-pstart						=> $report->{'pstart'},						
										-pend						=> $report->{'pend'},						
										-s_mean						=> $report->{'s_mean'},					
										-s_prob						=> $report->{'s_prob'},
										-d_score					=> $report->{'d_score'},					
										-c_max						=> $report->{'c_max'},
										-reliability				=> $report->{'reliability'},					
										-localization				=> $report->{'localization'},
						);						
					}
					else {
						$region = APPRIS::Analysis::CRASHRegion->new (
										-pstart						=> $report->{'pstart'},						
										-pend						=> $report->{'pend'},						
										-s_mean						=> $report->{'s_mean'},					
										-s_prob						=> $report->{'s_prob'},
										-d_score					=> $report->{'d_score'},					
										-c_max						=> $report->{'c_max'},
										-reliability				=> $report->{'reliability'},					
										-localization				=> $report->{'localization'},
						);						
					}
					push(@{$regions}, $region) if ( defined $region );
				}
			}
			
			# create Analysis object (for trans)			
			my ($method) = APPRIS::Analysis::CRASH->new (
							-result						=> $report->{'result'},
							-sp_score					=> $report->{'sp_score'},
							-tp_score					=> $report->{'tp_score'},
							-peptide_signal				=> $report->{'peptide_signal'},
							-mitochondrial_signal		=> $report->{'mitochondrial_signal'}
			);
			$method->regions($regions) if (defined $regions and (scalar(@{$regions}) > 0) );
			
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->crash($method);
				$analysis->number($analysis->number+1);
			}			
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);		
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::CRASH->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->crash($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

=head2 parse_thump_rst

  Arg [1]    : string $result
               Parse thump result
  Example    : use APPRIS::Parser qw(parse_thump);
               parse_firestar($result);
  Description: Parse output of firestar.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_thump_rst($)
{
	my ($result) = @_;
	my ($cutoffs);
	
	my (@results) = split('>',$result);	
	foreach my $transcript_result (@results)
	{
		#>ENST00000300482    length 1503 a.a.
		#helix number 1 start: 798       end: 818
		#helix number 2 start: 829       end: 841        damaged
		#helix number 3 start: 868       end: 886
		#helix number 4 start: 895       end: 908        damaged
		#helix number 5 start: 935       end: 955
		#helix number 6 start: 1024      end: 1041		
        if ( $transcript_result =~ /^([^\t]+)\t+length\s+([^\s]+)\s+a\.a\.\n+/ )
		{
			# get the helix coordinates
			my ($id) = $1;
			my ($transmembrane_list);
			$cutoffs->{$id}->{'num_tmh'} = 0;
			$cutoffs->{$id}->{'num_damaged_tmh'} = 0;
			while ($transcript_result =~ /^helix number \d+ start:\s+(\d+)\s+end:\s+(\d+)(\s*damaged|)$/mg ) {
				my ($start) = $1;
				my ($end) = $2;
				my ($damaged) = $3;
				if ( defined $start and defined $end ) {
					my ($helix_report) = {
						'start'	=> $start,
						'end'	=> $end
					};
					if ( $damaged =~ /damaged/ ) {
						$helix_report->{'damaged'} = 1;
						$cutoffs->{$id}->{'num_damaged_tmh'}++;
					}
					else {
						$cutoffs->{$id}->{'num_tmh'}++;						
					}
					push(@{$transmembrane_list}, $helix_report);					
				}
			}
			$cutoffs->{$id}->{'tmhs'} = $transmembrane_list if (defined $transmembrane_list); 

			# save result for each transcript
			$transcript_result =~ s/\n*#[^#]+#\n+#[^#]+#\n+#[^#]+#//mg;
			$cutoffs->{$id}->{'result'}='>'.$transcript_result;				
		}
	}
	$cutoffs->{'result'} = $result;
	
	return $cutoffs;
}

=head2 parse_thump

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse thump result
  Example    : use APPRIS::Parser qw(parse_thump);
               parse_thump($result);
  Description: Parse output of thump.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_thump($$)
{
	my ($gene, $result) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;
	
	# Create hash object from result
	my ($cutoffs) = parse_thump_rst($result);
	
	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
		
		# create method object
		if ( exists $cutoffs->{$transcript_id} ) {
			my ($report) = $cutoffs->{$transcript_id};
			my ($regions);
			if ( $transcript->translate ) {				
				if ( exists $report->{'tmhs'} ) {
					my ($strand) = $transcript->strand;
					foreach my $residue (@{$report->{'tmhs'}}) {
						my ($region);
						if ( $transcript->translate->cds ) { # with CDS coords from GTF file
							my ($pro_coord_start) = get_coords_from_residue($transcript, $residue->{'start'});
							my ($pro_coord_end) = get_coords_from_residue($transcript, $residue->{'end'});
							$residue->{'trans_strand'} = $strand;
							if ( $strand eq '-' ) {
								$residue->{'trans_end'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_start'} = $pro_coord_end->{'end'};                                              
							}
							else {
								$residue->{'trans_start'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_end'} = $pro_coord_end->{'end'};                                              
							}			
							$region = APPRIS::Analysis::THUMPRegion->new (
											-start								=> $residue->{'trans_start'},
											-end								=> $residue->{'trans_end'},
											-strand								=> $residue->{'trans_strand'},						
											-pstart								=> $residue->{'start'},					
											-pend								=> $residue->{'end'},
											-damaged							=> $residue->{'damaged'}										
							);
						}
						else {
							$region = APPRIS::Analysis::THUMPRegion->new (
											-pstart								=> $residue->{'start'},					
											-pend								=> $residue->{'end'},
											-damaged							=> $residue->{'damaged'}										
							);
						}	
						push(@{$regions}, $region);
					}
				}
			}
			# create Analysis object (for trans) Note: we only create an analysis object when trans has got translation 			
			my ($method) = APPRIS::Analysis::THUMP->new (
							-result							=> $report->{'result'},
							-num_tmh						=> $report->{'num_tmh'},
							-num_damaged_tmh				=> $report->{'num_damaged_tmh'}
			);
			$method->regions($regions) if (defined $regions and (scalar(@{$regions}) > 0) );
			
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->thump($method);
				$analysis->number($analysis->number+1);
			}
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);		
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::THUMP->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->thump($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

=head2 parse_proteo_rst

  Arg [1]    : string $result
               Parse proteo result
  Example    : use APPRIS::Parser qw(parse_proteo);
               parse_proteo($result);
  Description: Parse output of proteo.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_proteo_rst($)
{
	my ($result) = @_;
	my ($cutoffs);
		
	my (@results) = split('\n',$result);	
	#$peptide, $gene, $transcripts with, $transcripts without, $no_expts
	#You should delete all peptides with the "$transcripts from other genes" and "$transcripts without" column.
	#VSLENYFSLLNEK,ENSG00000000003.11,ENST00000373020.5,ENST00000612152.1;ENST00000614008.1,6
	#DWTDTNYYSEK,ENSG00000000003.11,ENST00000612152.1;ENST00000373020.5;ENST00000614008.1,,4
	foreach my $transcript_result (@results)
	{
		next unless (defined $transcript_result);
		
		$transcript_result =~ s/\n*$//;
		my (@split_line) = split(",", $transcript_result); # version with commas
		#my (@split_line) = split(/\t/, $transcript_result); # version with tabs
		next if ( scalar(@split_line) <= 0 );
		
		my ($mapped_peptide_sequence) = $split_line[0];
		my ($mapped_gene) = $split_line[1]; # Gene
		my ($mapped_transcript_list) = $split_line[2]; # GeneTranscripts
		my ($mapped_transc_without_list) = $split_line[3]; # GeneTranscriptsWithOut
		my ($num_experiments) = $split_line[4]; # TotalNumExperiments
		$mapped_peptide_sequence =~ s/^\"//;			$mapped_peptide_sequence =~ s/\"$//;
		$mapped_peptide_sequence=~s/Z$//;		
		$mapped_transcript_list =~ s/^\"//;				$mapped_transcript_list =~ s/\"$//;		
		 		
		unless(
			defined $mapped_gene and $mapped_gene ne '' and
			defined $mapped_peptide_sequence and $mapped_peptide_sequence ne '' and
			defined $mapped_transcript_list and $mapped_transcript_list ne '' and
			defined $num_experiments and $num_experiments ne ''			
		){
			next;
		}
		
		# Discard peptides with num of experiments less than 2
		next if ( $num_experiments < 2 );

		foreach my $mapped_transc_id (split(";",$mapped_transcript_list))
		{
			if ( defined $mapped_transc_id and $mapped_transc_id ne '' ) {
				if ( $mapped_transc_id =~ /^ENS/ ) { $mapped_transc_id =~ s/\.\d*$// } # delete suffix in Ensembl ids
				my ($peptide) = {						
								'sequence'			=> $mapped_peptide_sequence,
								'num_experiments'	=> $num_experiments,
				};
				$cutoffs->{$mapped_transc_id}->{'gene_id'} = $mapped_gene;
				push(@{$cutoffs->{$mapped_transc_id}->{'peptides'}}, $peptide);
				
				$cutoffs->{$mapped_transc_id}->{'peptide_evidence'} = $APPRIS::Utils::Constant::OK_LABEL;
				
				# save result for each transcript
				unless ( exists $cutoffs->{$mapped_transc_id}->{'result'} ) {
					$cutoffs->{$mapped_transc_id}->{'result'} = $transcript_result."\n";
					$cutoffs->{$mapped_transc_id}->{'num_peptides'} = 1;
					$cutoffs->{$mapped_transc_id}->{'num_experiments'} = $num_experiments;
				}
				else {
					$cutoffs->{$mapped_transc_id}->{'result'} .= $transcript_result."\n";
					$cutoffs->{$mapped_transc_id}->{'num_peptides'}++;
					$cutoffs->{$mapped_transc_id}->{'num_experiments'} += $num_experiments;
				}			
			}
		}		
	}
	$cutoffs->{'result'} = $result;
	
	return $cutoffs;
}

=head2 parse_proteo

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse proteo result
  Example    : use APPRIS::Parser qw(parse_proteo);
               parse_proteo($result);
  Description: Parse output of proteo.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_proteo($$)
{
	my ($gene, $result) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;
	
	# Create hash object from result
	my ($cutoffs) = parse_proteo_rst($result);

	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;			
			$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
		
		# create method object
		if ( exists $cutoffs->{$transcript_id} and exists $cutoffs->{$transcript_id}->{'peptide_evidence'} ) {
			my ($report) = $cutoffs->{$transcript_id};			
			my ($regions);
			if ( $transcript->translate and $transcript->translate->sequence ) {				
				if ( exists $report->{'peptides'} ) {
					my ($strand) = $transcript->strand;
					foreach my $residue (@{$report->{'peptides'}}) {
						my ($region);
						if ( $transcript->translate->cds ) { # with CDS coords from GTF file
							my ($protein_sequence) = $transcript->translate->sequence;
							my ($mapped_pep_sequence) = $residue->{'sequence'};
							my ($index_peptide_position) = index($protein_sequence, $mapped_pep_sequence);
							next if($index_peptide_position == -1); # we have not found the sequence
							my ($start_peptide_position) = $index_peptide_position + 1;
							my ($stop_peptide_position) = $start_peptide_position + length($mapped_pep_sequence)-1;
							
							my ($pro_coord_start) = get_coords_from_residue($transcript, $start_peptide_position);
							my ($pro_coord_end) = get_coords_from_residue($transcript, $stop_peptide_position);
							$residue->{'trans_strand'} = $strand;
							if ( $strand eq '-' ) {
								$residue->{'trans_end'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_start'} = $pro_coord_end->{'end'};                                              
							}
							else {
								$residue->{'trans_start'} = $pro_coord_start->{'start'};                                                
								$residue->{'trans_end'} = $pro_coord_end->{'end'};                                              
							}
							$region = APPRIS::Analysis::PROTEORegion->new (
											-start								=> $residue->{'trans_start'},
											-end								=> $residue->{'trans_end'},
											-strand								=> $residue->{'trans_strand'},						
											-pstart								=> $start_peptide_position,					
											-pend								=> $stop_peptide_position,
											-sequence							=> $residue->{'sequence'},
											-num_experiments					=> $residue->{'num_experiments'},
											-experiments						=> $residue->{'experiments'},
							);							
						}
						else {
							my ($protein_sequence) = $transcript->translate->sequence;
							my ($mapped_pep_sequence) = $residue->{'sequence'};
							my ($index_peptide_position) = index($protein_sequence, $mapped_pep_sequence);
							next if($index_peptide_position == -1); # we have not found the sequence
							my ($start_peptide_position) = $index_peptide_position + 1;
							my ($stop_peptide_position) = $start_peptide_position + length($mapped_pep_sequence)-1;							
							$region = APPRIS::Analysis::PROTEORegion->new (
											-pstart								=> $start_peptide_position,					
											-pend								=> $stop_peptide_position,
											-sequence							=> $residue->{'sequence'},
											-num_experiments					=> $residue->{'num_experiments'},
											-experiments						=> $residue->{'experiments'},
							);							
						}
						push(@{$regions}, $region);
					}
				}
			}
			
			# create Analysis object (for trans) Note: we only create an analysis object when trans has got translation 			
			my ($method) = APPRIS::Analysis::PROTEO->new (
							-result							=> $report->{'result'},
							-peptide_evidence				=> $report->{'peptide_evidence'},
							-num_peptides					=> $report->{'num_peptides'},
							-num_experiments				=> $report->{'num_experiments'}
			);
			$method->peptides($regions) if (defined $regions and (scalar(@{$regions}) > 0) );
			
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->proteo($method);
				$analysis->number($analysis->number+1);
			}
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);		
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::PROTEO->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->proteo($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

=head2 parse_appris_rst

  Arg [1]    : string $result
               Parse thump result
  Example    : use APPRIS::Parser qw(parse_appris_rst);
               parse_appris_rst($result);
  Description: Parse output of appris.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_appris_rst($$)
{
	my ($result, $result_label) = @_;
	my ($cutoffs);
	
	# Parse result with scores
	my (@resultscore) = split("\n", $result);
	foreach my $transcript_result2 (@resultscore)
	{
		next if ($transcript_result2 =~ /^#/); # skip the first line

		my (@rst2) = split("\t", $transcript_result2);

        # gene_id	gene_name	transcript_id	translation	status	biotype	no_codons ccds_id tsl (7)
        # fun_res
        # con_struct
        # vert_signal
        # dom_signal
        # u_evol
        # exon_signal
        # tmh_signal
        # pep_signal
        # mit_signal
        # prin_isoform
        
        if ( scalar(@rst2) == 19 ) {
			my ($gene_id) = $rst2[0];
			my ($gene_name) = $rst2[1];
			my ($transc_id) = $rst2[2];
			my ($translation) = $rst2[3];
			my ($biotype) = $rst2[4];
			my ($no_codons) = $rst2[5];
			my ($ccds_id) = $rst2[6];
			my ($tsl) = $rst2[7];
			my ($aa_length) = $rst2[8];

			my ($fun_res_annot) = $rst2[9];
			my ($con_struct_annot) = $rst2[10];
			my ($vert_con_annot) = $rst2[11];
			my ($dom_annot) = $rst2[12];
			my ($tmh_annot) = $rst2[13];
			my ($spep_mit_annot) = $rst2[14];
			my (@aux_crash_annot) = split(',', $spep_mit_annot);
			my ($spep_annot) = $aux_crash_annot[0];
			my ($mit_annot) = $aux_crash_annot[1];
			my ($u_evol_annot) = $rst2[15];
			my ($n_pep_annot) = $rst2[16];
			my ($p_isof_annot) = $rst2[17];
			my ($relia_annot) = $rst2[18];
			
			$cutoffs->{$transc_id}->{'functional_residues_score'} = $fun_res_annot;
			$cutoffs->{$transc_id}->{'homologous_structure_score'} = $con_struct_annot;
			$cutoffs->{$transc_id}->{'vertebrate_conservation_score'} = $vert_con_annot;
			$cutoffs->{$transc_id}->{'domain_score'} = $dom_annot;			
			$cutoffs->{$transc_id}->{'transmembrane_helices_score'} = $tmh_annot;
			$cutoffs->{$transc_id}->{'peptide_score'} = $spep_annot;
			$cutoffs->{$transc_id}->{'mitochondrial_score'} = $mit_annot;
			$cutoffs->{$transc_id}->{'unusual_evolution_score'} = $u_evol_annot;
			$cutoffs->{$transc_id}->{'peptide_evidence_score'} = $n_pep_annot;
			$cutoffs->{$transc_id}->{'principal_isoform_score'} = $p_isof_annot;
			$cutoffs->{$transc_id}->{'reliability'} = $relia_annot;
						
			# join transcript results
			if ( defined $transcript_result2 and ($transcript_result2 ne '') ) {
					$cutoffs->{$transc_id}->{'result'} = $transcript_result2;
			}
        }
	}
	$result =~ s/^#[^\n]*\n+//g; # delete the comment line to be homogenous
	# join results
	if ( defined $result and $result ne '') {
		$cutoffs->{'result'} = $result;
	}	

	# Parse result with labels
	my (@resultlabel) = split("\n", $result_label);
	foreach my $transcript_result2 (@resultlabel)
	{
		next if ($transcript_result2 =~ /^#/); # skip the first line

		my (@rst2) = split("\t", $transcript_result2);

        # gene_id	gene_name	transcript_id	translation	status	biotype	no_codons ccds_id (7)
        # fun_res
        # con_struct
        # vert_signal
        # dom_signal
        # u_evol
        # exon_signal
        # tmh_signal
        # pep_signal
        # mit_signal
        # prin_isoform
        
        if ( scalar(@rst2) == 19 ) {        	
			my ($gene_id) = $rst2[0];
			my ($gene_name) = $rst2[1];
			my ($transc_id) = $rst2[2];
			my ($translation) = $rst2[3];
			my ($biotype) = $rst2[4];
			my ($no_codons) = $rst2[5];
			my ($ccds_id) = $rst2[6];
			my ($tsl) = $rst2[7];
			my ($aa_length) = $rst2[8];

			my ($fun_res_annot) = $rst2[9];
			my ($con_struct_annot) = $rst2[10];
			my ($vert_con_annot) = $rst2[11];
			my ($dom_annot) = $rst2[12];
			my ($tmh_annot) = $rst2[13];
			my ($spep_mit_annot) = $rst2[14];
			my (@aux_crash_annot) = split(',', $spep_mit_annot);
			my ($spep_annot) = $aux_crash_annot[0];
			my ($mit_annot) = $aux_crash_annot[1];			
			my ($u_evol_annot) = $rst2[15];
			my ($n_pep_annot) = $rst2[16];
			my ($prin_isoform_annot) = $rst2[17];
			my ($relia_annot) = $rst2[18];
			
			$cutoffs->{$transc_id}->{'functional_residues_signal'} = $fun_res_annot;
			$cutoffs->{$transc_id}->{'homologous_structure_signal'} = $con_struct_annot;
			$cutoffs->{$transc_id}->{'vertebrate_conservation_signal'} = $vert_con_annot;
			$cutoffs->{$transc_id}->{'domain_signal'} = $dom_annot;			
			$cutoffs->{$transc_id}->{'transmembrane_helices_signal'} = $tmh_annot;
			$cutoffs->{$transc_id}->{'peptide_signal'} = $spep_annot;
			$cutoffs->{$transc_id}->{'mitochondrial_signal'} = $mit_annot;
			$cutoffs->{$transc_id}->{'unusual_evolution_signal'} = $u_evol_annot;
			$cutoffs->{$transc_id}->{'peptide_evidence_signal'} = $n_pep_annot;
			$cutoffs->{$transc_id}->{'principal_isoform_signal'} = $prin_isoform_annot;
			$cutoffs->{$transc_id}->{'reliability'} = $relia_annot;
			
			# join transcript results
			if ( defined $transcript_result2 and ($transcript_result2 ne '') ) {
					$cutoffs->{$transc_id}->{'result'} .= "\n#------------------\n";
					$cutoffs->{$transc_id}->{'result'} .= $transcript_result2;
			}
        }
	}
	$result_label =~ s/^#[^\n]*\n+//g; # delete the comment line to be homogenous
	# join results
	if ( exists $cutoffs->{'result'} and defined $result_label and $result_label ne '') {
		$cutoffs->{'result'} .= "\n#------------------\n";
		$cutoffs->{'result'} .= $result_label;
	}
	
	return $cutoffs;	
}

=head2 parse_appris

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse appris result
  Arg [3]    : string $result_label
               Parse appris label result
  Example    : use APPRIS::Parser qw(parse_appris);
               parse_appris($result,$result2);
  Description: Parse output of appris.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_appris($$$)
{
	my ($gene, $result, $result_label) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;
	
	# Create hash object from result
	my ($cutoffs) = parse_appris_rst($result, $result_label);
	
	# Create APPRIS object
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($transcript_ver);
		my ($transc_eid) = $transcript_id;
		if ( $transcript->version ) {
			$transcript_ver = $transcript->version;
			# NOTE: For this case, we don't use ensembl version
			#$transc_eid = $transcript_id.'.'.$transcript_ver;
		}
		my ($analysis);
		
		# create method object
		if ( exists $cutoffs->{$transcript_id} ) {
			my ($report) = $cutoffs->{$transcript_id};
			
			# create Analysis object (for trans)			
			my ($method) = APPRIS::Analysis::APPRIS->new (
							-functional_residues_score			=> $report->{'functional_residues_score'},
							-homologous_structure_score			=> $report->{'homologous_structure_score'},
							-vertebrate_conservation_score		=> $report->{'vertebrate_conservation_score'},
							-domain_score						=> $report->{'domain_score'},
							-transmembrane_helices_score		=> $report->{'transmembrane_helices_score'},
							-peptide_score						=> $report->{'peptide_score'},
							-mitochondrial_score				=> $report->{'mitochondrial_score'},							
							-unusual_evolution_score			=> $report->{'unusual_evolution_score'},
							-peptide_evidence_score				=> $report->{'peptide_evidence_score'},
							-principal_isoform_score			=> $report->{'principal_isoform_score'},
							-functional_residues_signal			=> $report->{'functional_residues_signal'},
							-homologous_structure_signal		=> $report->{'homologous_structure_signal'},
							-vertebrate_conservation_signal		=> $report->{'vertebrate_conservation_signal'},
							-domain_signal						=> $report->{'domain_signal'},
							-transmembrane_helices_signal		=> $report->{'transmembrane_helices_signal'},
							-peptide_signal						=> $report->{'peptide_signal'},
							-mitochondrial_signal				=> $report->{'mitochondrial_signal'},
							-unusual_evolution_signal			=> $report->{'unusual_evolution_signal'},
							-peptide_evidence_signal			=> $report->{'peptide_evidence_signal'},
							-principal_isoform_signal			=> $report->{'principal_isoform_signal'},
							-reliability						=> $report->{'reliability'},
							-result								=> $report->{'result'},							
			);
			$analysis = APPRIS::Analysis->new();
			if (defined $method) {
				$analysis->appris($method);
				$analysis->number($analysis->number+1);
			}			
		}
				
		# create Transcript object
		my ($transcript) = APPRIS::Transcript->new( -stable_id	=> $transcript_id );
		$transcript->version($transcript_ver) if (defined $transcript_ver);		
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
	}

	# create Analysis object (for gene)
	my ($method2) = APPRIS::Analysis::APPRIS->new( -result => $cutoffs->{'result'} );	
	my ($analysis2) = APPRIS::Analysis->new();
	if (defined $method2) {
		$analysis2->appris($method2);
		$analysis2->number($analysis2->number+1);
	}
	
	# create Gene object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	

	return $entity;
}

=head2 parse_appris_methods

  Arg [1]    : APPRIS::Gene $gene 
               APPRIS::Gene object
  Arg [2]    : string $result
               Parse firestar result
  Arg [3]    : string $result
               Parse matador3d result
  Arg [4]    : string $result
               Parse spade result
  Arg [5]    : string $result
               Parse corsair result
  Arg [6]    : string $result
               Parse crash result
  Arg [7]    : string $result
               Parse thump result
  Arg [8]    : string $result
               Parse inertia result
  Arg [9]    : string $result (optional)
               Parse inertia_maf result
  Arg [10]    : string $result (optional)
               Parse inertia_prank result
  Arg [11]    : string $result (optional)
               Parse inertia_kalign result               
  Arg [12]    : string $result (optional)
               Parse inertia_compara result               
  Arg [13]    : string $result (optional)
               Parse proteo result
  Arg [14]    : string $result (optional)
               Parse appris result
  Arg [15]    : string $result (optional)
               Parse appris_label result
  Example    : use APPRIS::Parser qw(parse_appris_methods);
               parse_appris_methods([$results]);
  Description: Parse output of methods of appris.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub parse_appris_methods($$$$$$$$;$;$;$)
{
	my (
		$gene,
		$firestar_result,
		$matador3d_result,
		$spade_result,
		$corsair_result,
		$crash_result,
		$thump_result,
		$inertia_result,
		$proteo_result,
		$appris_result,
		$appris_lb_result
	) = @_;
	
	my ($stable_id) = $gene->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;
	
	my ($firestar);
	my ($matador3d);
	my ($corsair);
	my ($spade);
	my ($thump);
	my ($crash);
	my ($inertia);
	my ($proteo);
	my ($appris);
	
	# get the reports for every method
	if ( defined $firestar_result ) {
		$firestar = parse_firestar($gene, $firestar_result);
	}
	if ( defined $matador3d_result ) {
		$matador3d = parse_matador3d($gene, $matador3d_result);
	}
	if ( defined $corsair_result ) {
		$corsair = parse_corsair($gene, $corsair_result);
	}	
	if ( defined $spade_result ) {
		$spade = parse_spade($gene, $spade_result);
	}
	if ( defined $thump_result ) {
		$thump = parse_thump($gene, $thump_result);
	}	
	if ( defined $crash_result ) {
		$crash = parse_crash($gene, $crash_result);
	}
	if ( defined $inertia_result ) {
		#$inertia = parse_inertia($gene, $inertia_result, $inertia_maf_result, $inertia_prank_result, $inertia_kalign_result, $inertia_compara_result);
		$inertia = parse_inertia($gene, $inertia_result);
	}
	if ( defined $proteo_result ) {
		$proteo = parse_proteo($gene, $proteo_result);
	}
	if ( defined $appris_result ) {
		$appris = parse_appris($gene, $appris_result, $appris_lb_result);
	}

	# get the results for each transcript
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($index) = $gene->{'_index_transcripts'}->{$transcript_id};
		my ($analysis) = APPRIS::Analysis->new();		

		# get firestar
		if ( $firestar and $firestar->transcripts->[$index] and $firestar->transcripts->[$index]->analysis ) {
			my ($result) = $firestar->transcripts->[$index];
			if ( $result->analysis->firestar ) {
				my ($method) = $result->analysis->firestar;
				if (defined $method) {
					$analysis->firestar($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}			
		# get matador3d
		if ( $matador3d and $matador3d->transcripts->[$index] and $matador3d->transcripts->[$index]->analysis ) {
			my ($result) = $matador3d->transcripts->[$index];
			if ( $result->analysis->matador3d ) {
				my ($method) = $result->analysis->matador3d;
				if (defined $method) {
					$analysis->matador3d($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}
		# get corsair
		if ( $corsair and $corsair->transcripts->[$index] and $corsair->transcripts->[$index]->analysis ) {
			my ($result) = $corsair->transcripts->[$index];
			if ( $result->analysis->corsair ) {
				my ($method) = $result->analysis->corsair;
				if (defined $method) {
					$analysis->corsair($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}		
		# get spade
		if ( $spade and $spade->transcripts->[$index] and $spade->transcripts->[$index]->analysis ) {
			my ($result) = $spade->transcripts->[$index];
			if ( $result->analysis->spade ) {
				my ($method) = $result->analysis->spade;
				if (defined $method) {
					$analysis->spade($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}
		# get thump
		if ( $thump and $thump->transcripts->[$index] and $thump->transcripts->[$index]->analysis ) {
			my ($result) = $thump->transcripts->[$index];
			if ( $result->analysis->thump ) {
				my ($method) = $result->analysis->thump;
				if (defined $method) {
					$analysis->thump($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}		
		# get crash
		if ( $crash and $crash->transcripts->[$index] and $crash->transcripts->[$index]->analysis ) {
			my ($result) = $crash->transcripts->[$index];
			if ( $result->analysis->crash ) {
				my ($method) = $result->analysis->crash;
				if (defined $method) {
					$analysis->crash($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}
		# get inertia
		if ( $inertia and $inertia->transcripts->[$index] and $inertia->transcripts->[$index]->analysis ) {
			my ($result) = $inertia->transcripts->[$index];
			if ( $result->analysis->inertia ) {
				my ($method) = $result->analysis->inertia;
				if (defined $method) {
					$analysis->inertia($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}
		# get proteo
		if ( $proteo and $proteo->transcripts->[$index] and $proteo->transcripts->[$index]->analysis ) {
			my ($result) = $proteo->transcripts->[$index];
			if ( $result->analysis->proteo ) {
				my ($method) = $result->analysis->proteo;
				if (defined $method) {
					$analysis->proteo($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}		
		# get appris
		if ( $appris and $appris->transcripts->[$index] and $appris->transcripts->[$index]->analysis ) {
			my ($result) = $appris->transcripts->[$index];
			if ( $result->analysis->appris ) {
				my ($method) = $result->analysis->appris;
				if (defined $method) {
					$analysis->appris($method);
					$analysis->number($analysis->number+1);
				}
			}			
		}
				
		# create object
		my ($transcript) = APPRIS::Transcript->new
		(
			-stable_id	=> $transcript_id,
		);
		$transcript->analysis($analysis) if (defined $analysis);
		push(@{$transcripts}, $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts		
	}
	
	# get the results for gene
	my ($analysis2) = APPRIS::Analysis->new();
	
	# get firestar
	if ( $firestar and $firestar->analysis and $firestar->analysis->firestar ) {
		my ($method2) = $firestar->analysis->firestar;
		if (defined $method2) {
			$analysis2->firestar($method2);
			$analysis2->number($analysis2->number+1);
		}			
	}
	# get matador3d
	if ( $matador3d and $matador3d->analysis and $matador3d->analysis->matador3d ) {
		my ($method2) = $matador3d->analysis->matador3d;
		if (defined $method2) {
			$analysis2->matador3d($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get corsair
	if ( $corsair and $corsair->analysis and $corsair->analysis->corsair ) {
		my ($method2) = $corsair->analysis->corsair;
		if (defined $method2) {
			$analysis2->corsair($method2);
			$analysis2->number($analysis2->number+1);
		}
	}	
	# get spade
	if ( $spade and $spade->analysis and $spade->analysis->spade ) {
		my ($method2) = $spade->analysis->spade;
		if (defined $method2) {
			$analysis2->spade($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get thump
	if ( $thump and $thump->analysis and $thump->analysis->thump ) {
		my ($method2) = $thump->analysis->thump;
		if (defined $method2) {
			$analysis2->thump($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get crash
	if ( $crash and $crash->analysis and $crash->analysis->crash ) {
		my ($method2) = $crash->analysis->crash;
		if (defined $method2) {
			$analysis2->crash($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get inertia
	if ( $inertia and $inertia->analysis and $inertia->analysis->inertia ) {
		my ($method2) = $inertia->analysis->inertia;
		if (defined $method2) {
			$analysis2->inertia($method2);
			$analysis2->number($analysis2->number+1);
		}
	}	
	# get proteo
	if ( $proteo and $proteo->analysis and $proteo->analysis->proteo ) {
		my ($method2) = $proteo->analysis->proteo;
		if (defined $method2) {
			$analysis2->proteo($method2);
			$analysis2->number($analysis2->number+1);
		}
	}	
	
	# create object
	my ($entity) = APPRIS::Gene->new( -stable_id => $stable_id );
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	
	
	return $entity;	
}

=head2 create_appris_entity

  Arg [1]    : string $data_file
               File of gencode data as GTF format
  Arg [2]    : string $transc_file (optional)
               File of gencode transcript sequence as Fasta format
  Arg [3]    : string $transl_file (optional)
               File of gencode translation sequence as Fasta format
  Arg [4]    : string $result
               Parse firestar result
  Arg [5]    : string $result
               Parse matador3d result
  Arg [6]    : string $result
               Parse spade result
  Arg [7]    : string $result
               Parse proteo result
  Arg [8]    : string $result
               Parse crash result
  Arg [9]    : string $result
               Parse thump result
  Arg [10]   : string $result
               Parse inertia result
  Arg [11]   : string $result
               Parse inertia result
  Arg [12]   : string $result (optional)
               Parse appris result
  Arg [13]   : string $result (optional)
               Parse appris_label result
  Example    : use APPRIS::Parser qw(create_appris_entity);
               create_appris_entity([$results]);
  Description: Parse output of methods of appris.
  Returntype : APPRIS::Gene or undef
  Exceptions : return undef
  Caller     : generally on error

=cut

sub create_appris_entity($$$$$$$$$$$$$$)
{
	my (
		$data_file,
		$transc_file,
		$transl_file,
		$firestar_result,
		$matador3d_result,
		$spade_result,
		$corsair_result,
		$crash_result,
		$thump_result,
		$inertia_result,
		$proteo_result,
		$appris_result,
		$appris_lb_result,
		$ids
	) = @_;
	
	# get entity from input files
	my ($entity);
	my ($entities);
	if ( defined $data_file and defined $transc_file and defined $transl_file ) {
		($entities) = parse_infiles($data_file, $transc_file, $transl_file);
	}
	elsif ( defined $transl_file ) {
		($entities) = parse_transl_data($transl_file);
	}
	if ( defined $entities and UNIVERSAL::isa($entities, 'ARRAY') and (scalar(@{$entities}) > 0) ) {
		$entity = $entities->[0];
	} else {
		return undef;
	}
	
	my ($stable_id) = $entity->stable_id;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;
	
	my ($firestar);
	my ($matador3d);
	my ($corsair);
	my ($spade);
	my ($cexonic);
	my ($thump);
	my ($crash);
	my ($inertia);
	my ($proteo);
	my ($appris);
	
	# get the reports for every method
	if ( defined $firestar_result ) {
		$firestar = parse_firestar($entity, $firestar_result);
	}
	if ( defined $matador3d_result ) {
		$matador3d = parse_matador3d($entity, $matador3d_result);
	}
	if ( defined $corsair_result ) {
		$corsair = parse_corsair($entity, $corsair_result);
	}	
	if ( defined $spade_result ) {
		$spade = parse_spade($entity, $spade_result);
	}
	if ( defined $thump_result ) {
		$thump = parse_thump($entity, $thump_result);
	}	
	if ( defined $crash_result ) {
		$crash = parse_crash($entity, $crash_result);
	}
	if ( defined $inertia_result ) {
		$inertia = parse_inertia($entity, $inertia_result);
	}
	if ( defined $proteo_result ) {
		$proteo = parse_proteo($entity, $proteo_result);
	}
	if ( defined $appris_result ) {
		$appris = parse_appris($entity, $appris_result, $appris_lb_result);
	}

	# get the results for each transcript
	foreach my $transcript (@{$entity->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($index);
		if ( !defined $ids ) { # by default all
			$index = $entity->{'_index_transcripts'}->{$transcript_id};
		}
		elsif ( defined $ids and $ids =~ /$transcript_id/ ) {
			$index = $entity->{'_index_transcripts'}->{$transcript_id};
		}
		if ( defined $index ) {
			my ($analysis) = APPRIS::Analysis->new();
			# get firestar
			if ( $firestar and $firestar->transcripts->[$index] and $firestar->transcripts->[$index]->analysis ) {
				my ($result) = $firestar->transcripts->[$index];
				if ( $result->analysis->firestar ) {
					my ($method) = $result->analysis->firestar;
					if (defined $method) {
						$analysis->firestar($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}			
			# get matador3d
			if ( $matador3d and $matador3d->transcripts->[$index] and $matador3d->transcripts->[$index]->analysis ) {
				my ($result) = $matador3d->transcripts->[$index];
				if ( $result->analysis->matador3d ) {
					my ($method) = $result->analysis->matador3d;
					if (defined $method) {
						$analysis->matador3d($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}
			# get corsair
			if ( $corsair and $corsair->transcripts->[$index] and $corsair->transcripts->[$index]->analysis ) {
				my ($result) = $corsair->transcripts->[$index];
				if ( $result->analysis->corsair ) {
					my ($method) = $result->analysis->corsair;
					if (defined $method) {
						$analysis->corsair($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}		
			# get spade
			if ( $spade and $spade->transcripts->[$index] and $spade->transcripts->[$index]->analysis ) {
				my ($result) = $spade->transcripts->[$index];
				if ( $result->analysis->spade ) {
					my ($method) = $result->analysis->spade;
					if (defined $method) {
						$analysis->spade($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}
			# get thump
			if ( $thump and $thump->transcripts->[$index] and $thump->transcripts->[$index]->analysis ) {
				my ($result) = $thump->transcripts->[$index];
				if ( $result->analysis->thump ) {
					my ($method) = $result->analysis->thump;
					if (defined $method) {
						$analysis->thump($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}		
			# get crash
			if ( $crash and $crash->transcripts->[$index] and $crash->transcripts->[$index]->analysis ) {
				my ($result) = $crash->transcripts->[$index];
				if ( $result->analysis->crash ) {
					my ($method) = $result->analysis->crash;
					if (defined $method) {
						$analysis->crash($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}
			# get inertia
			if ( $inertia and $inertia->transcripts->[$index] and $inertia->transcripts->[$index]->analysis ) {
				my ($result) = $inertia->transcripts->[$index];
				if ( $result->analysis->inertia ) {
					my ($method) = $result->analysis->inertia;
					if (defined $method) {
						$analysis->inertia($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}
			# get proteo
			if ( $proteo and $proteo->transcripts->[$index] and $proteo->transcripts->[$index]->analysis ) {
				my ($result) = $proteo->transcripts->[$index];
				if ( $result->analysis->proteo ) {
					my ($method) = $result->analysis->proteo;
					if (defined $method) {
						$analysis->proteo($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}		
			# get appris
			if ( $appris and $appris->transcripts->[$index] and $appris->transcripts->[$index]->analysis ) {
				my ($result) = $appris->transcripts->[$index];
				if ( $result->analysis->appris ) {
					my ($method) = $result->analysis->appris;
					if (defined $method) {
						$analysis->appris($method);
						$analysis->number($analysis->number+1);
					}
				}			
			}
					
			# add analysis into transcript
			$transcript->analysis($analysis) if (defined $analysis);
			push(@{$transcripts}, $transcript);
			$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
		}
	}
	
	# get the results for gene
	my ($analysis2) = APPRIS::Analysis->new();
	
	# get firestar
	if ( $firestar and $firestar->analysis and $firestar->analysis->firestar ) {
		my ($method2) = $firestar->analysis->firestar;
		if (defined $method2) {
			$analysis2->firestar($method2);
			$analysis2->number($analysis2->number+1);
		}			
	}
	# get matador3d
	if ( $matador3d and $matador3d->analysis and $matador3d->analysis->matador3d ) {
		my ($method2) = $matador3d->analysis->matador3d;
		if (defined $method2) {
			$analysis2->matador3d($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get corsair
	if ( $corsair and $corsair->analysis and $corsair->analysis->corsair ) {
		my ($method2) = $corsair->analysis->corsair;
		if (defined $method2) {
			$analysis2->corsair($method2);
			$analysis2->number($analysis2->number+1);
		}
	}	
	# get spade
	if ( $spade and $spade->analysis and $spade->analysis->spade ) {
		my ($method2) = $spade->analysis->spade;
		if (defined $method2) {
			$analysis2->spade($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get thump
	if ( $thump and $thump->analysis and $thump->analysis->thump ) {
		my ($method2) = $thump->analysis->thump;
		if (defined $method2) {
			$analysis2->thump($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get crash
	if ( $crash and $crash->analysis and $crash->analysis->crash ) {
		my ($method2) = $crash->analysis->crash;
		if (defined $method2) {
			$analysis2->crash($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get inertia
	if ( $inertia and $inertia->analysis and $inertia->analysis->inertia ) {
		my ($method2) = $inertia->analysis->inertia;
		if (defined $method2) {
			$analysis2->inertia($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	# get proteo
	if ( $proteo and $proteo->analysis and $proteo->analysis->proteo ) {
		my ($method2) = $proteo->analysis->proteo;
		if (defined $method2) {
			$analysis2->proteo($method2);
			$analysis2->number($analysis2->number+1);
		}
	}
	
	
	# add analysis into gene
	$entity->transcripts($transcripts, $index_transcripts) if (defined $transcripts and defined $index_transcripts);
	$entity->analysis($analysis2) if (defined $analysis2);	
	
	return $entity;	
}

# *********************** #
# *** PRIVATE METHODS *** #
# *********************** #

# Get the id and the version from Ensembl identifiers
sub _get_id_version($)
{
	my ($i_id) = @_;
	my ($id, $version) = (undef,undef);
	
	if ( $i_id =~ /^(ENS[^\.]*)\.(\d*)$/ ) {
		($id, $version) = ($1, $2);
	}
	else {
		$id = $i_id;
	}
	return ($id, $version);
		
} # End _get_id_version

# Is RefSeq gene data
sub _is_refseq($)
{
	my ($file) = @_;
	my ($is_refseq) = 0;	
	

	open (IN_FILE, $file) or throw('Can not open file');
	while ( my $line = <IN_FILE> )
	{
		if ( $line =~ /^#/ ) {
			if ( $line =~ /processor NCBI annotwriter/ ) { $is_refseq = 1 }
		}
		else {
			my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$attributes) = split("\t", $line);
			if ( (lc($source) =~ /refseq/) or (lc($source) =~ /bestrefseq/) or (lc($source) =~ /curated genomic/) or (lc($source) =~ /gnomon/) ) { if ( $attributes =~ /ID\=/ ) { $is_refseq = 1 } }
			else { last; }
		}
	}
	close(IN_FILE);
	
	return $is_refseq;
	
} # End _is_refseq

sub _parse_dataline($)
{
	my ($line) = @_;
	my ($fields);
	
	if ( defined $line and ($line ne '') and ($line !~ /^#/) ) {
		my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$attributes) = split("\t", $line);
		next unless(defined $chr and 
					defined $source and
					defined $type and
					defined $start and
					defined $end and
					defined $score and
					defined $strand and
					defined $phase and 
					defined $attributes);
		if (defined $chr and $chr=~/chr(\w*)/) { $chr = $1 if(defined $1) }		
		
		#store ids and additional information in second hash
		my ($attribs);
		my (@add_attributes) = split(";", $attributes);				
		for ( my $i=0; $i<scalar @add_attributes; $i++ )
		{
			my ($c_type);
			my ($c_value);
			if ( $add_attributes[$i] =~ /^(.+)\s+\"([^\"]*)\"$/ ) {
				$c_type  = $1;
				$c_value = $2;			
			}
			elsif ( $add_attributes[$i] =~ /^(.+)=(.+)$/ ) {
				$c_type  = $1;
				$c_value = $2;
			}
			if(	defined $c_type and !($c_type=~/^\s*$/) and
				defined $c_value and !($c_value=~/^\s*$/))
			{
				$c_type =~ s/^\s*//g;
				$c_value =~ s/"//g;
				if(!exists($attribs->{$c_type})) {
					$attribs->{$c_type} = $c_value;						
				}
				else {
					if ( $c_type eq 'tag' ) {
						$attribs->{$c_type} .= ','.$c_value;
					}
				}				
			}
		}
		
		#store nine columns in hash
		$fields = {
				chr        => $chr,
				source     => $source,
				type       => $type,
				start      => $start,
				end        => $end,
				score      => $score,
				strand     => $strand,
				phase      => $phase,
				attrs      => $attribs,
		};
	}
			
	return $fields;
	
} # End _parse_dataline

# Parse GFT file
sub _parse_indata($)
{
	my ($file) = @_;
	my ($data);	
	return $data unless (-e $file and (-s $file > 0) );
	
	open (IN_FILE, $file) or throw('Can not open file');
	while ( my $line = <IN_FILE> )
	{
		my ($fields) = _parse_dataline($line);
		next unless ( defined $fields );
		my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$attribs) = (
			$fields->{'chr'},
			$fields->{'source'},
			$fields->{'type'},
			$fields->{'start'},
			$fields->{'end'},
			$fields->{'score'},
			$fields->{'strand'},
			$fields->{'phase'},
			$fields->{'attrs'}
		);
		
		# Always we have Gene Id
		if(	exists $attribs->{'gene_id'} and defined $attribs->{'gene_id'} )
		{
			my ($gene_id, $gene_version) = _get_id_version($attribs->{'gene_id'});
			my ($transcript_id, $trans_version) = (undef, undef);
			if ( exists $attribs->{'transcript_id'} and defined $attribs->{'transcript_id'} ) {
				($transcript_id, $trans_version) = _get_id_version($attribs->{'transcript_id'});				
			}

			if (defined $gene_id and ($type eq 'gene') ) # Gene Information
			{
				$data->{$gene_id}->{'chr'} = $chr if(defined $chr);			
				$data->{$gene_id}->{'start'} = $start if(defined $start);
				$data->{$gene_id}->{'end'} = $end if(defined $end);
				$data->{$gene_id}->{'strand'} = $strand if(defined $strand);

				if (defined $source)
				{
					$data->{$gene_id}->{'source'} = $source; 					
				}
				if(exists $attribs->{'gene_status'} and defined $attribs->{'gene_status'})
				{
					$data->{$gene_id}->{'status'} = $attribs->{'gene_status'};	
				}			
				if(exists $attribs->{'gene_type'} and defined $attribs->{'gene_type'})
				{
					$data->{$gene_id}->{'biotype'} = $attribs->{'gene_type'};	
				}
				elsif(exists $attribs->{'gene_biotype'} and defined $attribs->{'gene_biotype'})
				{
					$data->{$gene_id}->{'biotype'} = $attribs->{'gene_biotype'};	
				}
				if(exists $attribs->{'gene_name'} and defined $attribs->{'gene_name'})
				{
					$data->{$gene_id}->{'external_id'} = $attribs->{'gene_name'};	
				}
				if(exists $attribs->{'havana_gene'} and defined $attribs->{'havana_gene'})
				{
					$data->{$gene_id}->{'havana_gene'} = $attribs->{'havana_gene'};	
				}
				if(exists $attribs->{'level'} and defined $attribs->{'level'})
				{
					$data->{$gene_id}->{'level'} = $attribs->{'level'};	
				}
				if (defined $gene_version)
				{
					$data->{$gene_id}->{'version'} = $gene_version;
				}
				elsif(exists $attribs->{'gene_version'} and defined $attribs->{'gene_version'})
				{
					$data->{$gene_id}->{'version'} = $attribs->{'gene_version'};
				}				
			}
			elsif (defined $gene_id and defined $transcript_id and ($type eq 'transcript') ) # Transcript Information
			{
				my ($transcript);
				$transcript->{'chr'} = $chr if(defined $chr);
				$transcript->{'start'} = $start if(defined $start);
				$transcript->{'end'} = $end if(defined $end);
				$transcript->{'strand'} = $strand if(defined $strand);					

				if (defined $source)
				{
					$transcript->{'source'} = $source;
				}
				if(exists $attribs->{'transcript_status'} and defined $attribs->{'transcript_status'})
				{
					$transcript->{'status'} = $attribs->{'transcript_status'};	
				}			
				if(exists $attribs->{'transcript_type'} and defined $attribs->{'transcript_type'})
				{
					$transcript->{'biotype'} = $attribs->{'transcript_type'};	
				}
				elsif(exists $attribs->{'transcript_biotype'} and defined $attribs->{'transcript_biotype'})
				{
					$transcript->{'biotype'} = $attribs->{'transcript_biotype'};	
				}				
				if(exists $attribs->{'transcript_name'} and defined $attribs->{'transcript_name'})
				{
					$transcript->{'external_id'} = $attribs->{'transcript_name'};	
				}
				if(exists $attribs->{'havana_transcript'} and defined $attribs->{'havana_transcript'})
				{
					$transcript->{'havana_transcript'} = $attribs->{'havana_transcript'};	
				}
				if(exists $attribs->{'protein_id'} and defined $attribs->{'protein_id'})
				{
					my ($x) = $attribs->{'protein_id'};
					if ( $x =~ /^ENS/ ) { $x =~ s/\.\d*$// } # delete suffix in Ensembl ids
					$transcript->{'protein_id'} = $x;
				}
				if(exists $attribs->{'level'} and defined $attribs->{'level'})
				{
					$transcript->{'level'} = $attribs->{'level'};	
				}
				if (defined $trans_version)
				{
					$transcript->{'version'} = $trans_version;
				}
				elsif(exists $attribs->{'transcript_version'} and defined $attribs->{'transcript_version'})
				{
					$transcript->{'version'} = $attribs->{'transcript_version'};
				}								
				if(exists $attribs->{'ccdsid'} and defined $attribs->{'ccdsid'})
				{
					$transcript->{'ccdsid'} = $attribs->{'ccdsid'};	
				}
				elsif(exists $attribs->{'ccds_id'} and defined $attribs->{'ccds_id'})
				{
					$transcript->{'ccdsid'} = $attribs->{'ccds_id'};	
				}
				if(exists $attribs->{'tag'} and defined $attribs->{'tag'})
				{
					$transcript->{'tag'} = $attribs->{'tag'};
				}
				if(exists $attribs->{'transcript_support_level'} and defined $attribs->{'transcript_support_level'})
				{
					my ($t) = $attribs->{'transcript_support_level'};
					if ( $t =~ /^([0-9]*)\s*/ ) { $transcript->{'tsl'} = $1 }					
				}
				elsif(exists $attribs->{'tsl'} and defined $attribs->{'tsl'})
				{
					$transcript->{'tsl'} = $attribs->{'tsl'};
				}
					
				$data->{$gene_id}->{'transcripts'}->{$transcript_id} = $transcript if(defined $transcript);
			}
			elsif (defined $gene_id and defined $transcript_id and ($type eq 'exon') ) # Exon Information
			{
				my ($exon);
				$exon->{'start'} = $start if(defined $start);
				$exon->{'end'} = $end if(defined $end);
				$exon->{'strand'} = $strand if(defined $strand);
				
				if(exists $attribs->{'exon_id'} and defined $attribs->{'exon_id'})
				{
					my ($x) = $attribs->{'exon_id'};
					if ( $x =~ /^ENS/ ) { $x =~ s/\.\d*$// } # delete suffix in Ensembl ids
					$exon->{'exon_id'} = $x;
				}				
						
				push(@{$data->{$gene_id}->{'transcripts'}->{$transcript_id}->{'exons'}},$exon);
			}			
			elsif (defined $gene_id and defined $transcript_id and ($type eq 'CDS') ) # CDS Information
			{
				my ($cds);
				$cds->{'start'} = $start if(defined $start);
				$cds->{'end'} = $end if(defined $end);
				$cds->{'strand'} = $strand if(defined $strand);
				$cds->{'phase'} = $phase if(defined $phase);
				
				if( (exists $attribs->{'exon_id'} and defined $attribs->{'exon_id'}) or (exists $attribs->{'cds_id'} and defined $attribs->{'cds_id'}) )
				{
					my ($x) = (exists $attribs->{'exon_id'} and defined $attribs->{'exon_id'}) ? $attribs->{'exon_id'} : $attribs->{'cds_id'};
					if ( $x =~ /^ENS/ ) { $x =~ s/\.\d*$// } # delete suffix in Ensembl ids
					$cds->{'exon_id'} = $x;
				}				
				
				push(@{$data->{$gene_id}->{'transcripts'}->{$transcript_id}->{'cds'}},$cds);
				
				if(exists $attribs->{'protein_id'} and defined $attribs->{'protein_id'} and !exists $data->{$gene_id}->{'transcripts'}->{$transcript_id}->{'protein_id'} ) # ensembl exception
				{
					my ($x) = $attribs->{'protein_id'};
					if ( $x =~ /^ENS/ ) { $x =~ s/\.\d*$// } # delete suffix in Ensembl ids
					$data->{$gene_id}->{'transcripts'}->{$transcript_id}->{'protein_id'} = $x;
				}				
			}
			elsif (defined $gene_id and defined $transcript_id and ($type eq 'start_codon') ) # Codon Information
			{
				my ($codon);
				$codon->{'type'}='start';
				$codon->{'start'} = $start if(defined $start);
				$codon->{'end'} = $end if(defined $end);
				$codon->{'strand'} = $strand if(defined $strand);
				$codon->{'phase'} = $phase if(defined $phase);
				
				push(@{$data->{$gene_id}->{'transcripts'}->{$transcript_id}->{'codons'}},$codon) if(defined $codon);
			}
			elsif (defined $gene_id and defined $transcript_id and ($type eq 'stop_codon') ) # Codon Information
			{
				my ($codon);
				$codon->{'type'}='stop';
				$codon->{'start'} = $start if(defined $start);
				$codon->{'end'} = $end if(defined $end);
				$codon->{'strand'} = $strand if(defined $strand);
				$codon->{'phase'} = $phase if(defined $phase);
				
				push(@{$data->{$gene_id}->{'transcripts'}->{$transcript_id}->{'codons'}},$codon) if(defined $codon);
			}
			$data->{$gene_id}->{'raw'} .= $line; # Save Raw Data
		}
		else
		{
			throw('Wrong entity');
		}
	}
	close(IN_FILE);
	
	return $data;
	
} # End _parse_indata

# Parse GFF3 file from RefSeq
sub _parse_indata_refseq($)
{
	my ($file) = @_;
	my ($data);
	my ($cache_transcId);
	return $data unless (-e $file and (-s $file > 0) );
	
	#my ($molecule) = ',transcript,primary_transcript,mRNA,ncRNA,rRNA,tRNA,V_gene_segment,cDNA_match,C_gene_segment,D_gene_segment,D_loop,J_gene_segment,long_terminal_repeat,match,';
	
	my $_extract_dbXref = sub {
		my ($data,$patt) = @_;
		my ($match);
		if ( defined $data and $data =~ /$patt\:([^\,]*)/ ) {
			$match = $1;
		}
		return $match;
	};
	
	open (IN_FILE, $file) or throw('Can not open file');
	while ( my $line = <IN_FILE> )
	{
		my ($fields) = _parse_dataline($line);
		next unless ( defined $fields );
		my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$attribs) = (
			$fields->{'chr'},
			$fields->{'source'},
			$fields->{'type'},
			$fields->{'start'},
			$fields->{'end'},
			$fields->{'score'},
			$fields->{'strand'},
			$fields->{'phase'},
			$fields->{'attrs'}
		);

		# Only BestRefSeq features
		#next unless ( ($source =~ /BestRefSeq/) or ($source =~ /Curated Genomic/) or ($source =~ /Gnomon/) );
		
		# Always we have Gene Id
		if(	exists $attribs->{'ID'} and defined $attribs->{'ID'} )
		{			
			my ($ID) = $attribs->{'ID'};
			my ($Parent);
			if ( exists $attribs->{'Parent'} and defined $attribs->{'Parent'} ) {
				($Parent) = $attribs->{'Parent'};				
			}
			my ($gene_id) = $_extract_dbXref->($attribs->{'Dbxref'}, 'GeneID');
			my ($accesion_id) = $_extract_dbXref->($attribs->{'Dbxref'}, 'Genbank');
			my ($ccds_id) = $_extract_dbXref->($attribs->{'Dbxref'}, 'CCDS');
			
			if ( defined $gene_id ) {
				if ( $type eq 'gene' ) # Gene Information
				{
					$data->{$gene_id}->{'chr'} = $chr if(defined $chr);			
					$data->{$gene_id}->{'start'} = $start if(defined $start);
					$data->{$gene_id}->{'end'} = $end if(defined $end);
					$data->{$gene_id}->{'strand'} = $strand if(defined $strand);
					
					if (defined $source)
					{
						$data->{$gene_id}->{'source'} = $source; 					
						if ( ($source eq 'BestRefSeq') or ($source eq 'Curated Genomic') ) { $data->{$gene_id}->{'level'} = '1' }
						elsif ( $source eq 'Gnomon' ) { $data->{$gene_id}->{'level'} = '3' }
					}
					if(exists $attribs->{'Name'} and defined $attribs->{'Name'})
					{
						$data->{$gene_id}->{'external_id'} = $attribs->{'Name'};	
					}
				}
				elsif ( defined $accesion_id and ($type eq 'exon') ) # Exon Information
				{
					my ($exon);
					$exon->{'start'} = $start if(defined $start);
					$exon->{'end'} = $end if(defined $end);
					$exon->{'strand'} = $strand if(defined $strand);
					
					$exon->{'exon_id'} = $ID;			
									
					push(@{$data->{$gene_id}->{'transcripts'}->{$accesion_id}->{'exons'}},$exon);
				}			
				elsif ( defined $accesion_id and ($type eq 'CDS') ) # CDS Information
				{
					my ($cds);
					$cds->{'start'} = $start if(defined $start);
					$cds->{'end'} = $end if(defined $end);
					$cds->{'strand'} = $strand if(defined $strand);
					$cds->{'phase'} = $phase if(defined $phase);
									
					$cds->{'cds_id'} = $ID;
					
					if ( exists $cache_transcId->{$Parent} and defined $cache_transcId->{$Parent} ) {
						my ($transc_id) = $cache_transcId->{$Parent};
						$data->{$gene_id}->{'transcripts'}->{$transc_id}->{'protein_id'} = $accesion_id;
						$data->{$gene_id}->{'transcripts'}->{$transc_id}->{'ccdsid'} = $ccds_id if ( defined $ccds_id );
						
						push(@{$data->{$gene_id}->{'transcripts'}->{$transc_id}->{'cds'}},$cds);				
					}				
				}
				elsif ( defined $accesion_id ) # Any type of molecule: mRNA, transcript, ncRNA, etc.
				{				
					my ($transcript);
					$transcript->{'chr'} = $chr if(defined $chr);
					$transcript->{'start'} = $start if(defined $start);
					$transcript->{'end'} = $end if(defined $end);
					$transcript->{'strand'} = $strand if(defined $strand);					
	
					if (defined $source)
					{
						$transcript->{'source'} = $source;
						if ( ($source eq 'BestRefSeq') or ($source eq 'Curated Genomic') ) { $transcript->{'level'} = '1' }
						elsif ( $source eq 'Gnomon' ) { $transcript->{'level'} = '3' }
					}
					if(exists $attribs->{'Name'} and defined $attribs->{'Name'})
					{
						$transcript->{'external_id'} = $attribs->{'Name'};	
					}
					if(exists $attribs->{'tag'} and defined $attribs->{'tag'})
					{
						$transcript->{'tag'} = $attribs->{'tag'};
					}
					if(exists $attribs->{'tsl'} and defined $attribs->{'tsl'})
					{
						$transcript->{'tsl'} = $attribs->{'tsl'};
					}				
					
					# cache transc Ids
					$cache_transcId->{$ID} = $accesion_id;
									
					# HARD-CORE attrs!!
					#if ( ($transcript_id =~ /^NM\_/) or ($transcript_id =~ /^NR\_/) or ($transcript_id =~ /^NP\_/) or ($transcript_id =~ /^YP\_/) ) {
					if ( ($accesion_id =~ /^NM\_/) or ($accesion_id =~ /^XM\_/) ) {
						$transcript->{'status'} = 'KNOWN';
					}
					else {
						$transcript->{'status'} = 'UNKNOWN';
					}
					$transcript->{'biotype'} = $type;
									
					# NOTE: HARD-CORE!!! We have decided the all mRNA from RefSeq have start/stop codons
					for my $type ('start','stop') {
						my ($codon);
						$codon->{'type'}=$type;
						#$codon->{'start'} = $start if(defined $start);
						#$codon->{'end'} = $end if(defined $end);
						#$codon->{'strand'} = $strand if(defined $strand);
						#$codon->{'phase'} = $phase if(defined $phase);
						$codon->{'start'} = $transcript->{'start'};
						$codon->{'end'} = $transcript->{'end'};
						$codon->{'strand'} = $transcript->{'strand'};
						$codon->{'phase'} = 0;
						push(@{$transcript->{'codons'}},$codon) if(defined $codon);
					}
					
					$data->{$gene_id}->{'transcripts'}->{$accesion_id} = $transcript if(defined $transcript);				
				}
						
				$data->{$gene_id}->{'raw'} .= $line; # Save Raw Data
			}
		}
		else
		{
			throw('Wrong entity');
		}
	}
	close(IN_FILE);
	
	return $data;
	
} # End _parse_indata_refseq

sub _parse_refseq_gbk($)
{
	
} # End _parse_refseq_gbk

sub _parse_inseq_transc($)
{
	my ($file) = @_;
	my ($data);

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() )
		{
			my ($sequence_id);
			if ( $seq->id =~ /^[sp|tr]\|([^|]*)\|([^|]*)/ ) { # UniProt sequences
				my ($id1) = $1;
				my ($id2) = $2;
				$sequence_id = $id1;
			}
			elsif ( $seq->id =~ /^gi\|[^|]*\|[^|]*\|([^|]*)/ ) { # RefSeq sequences
				my ($id1) = $1;
				$sequence_id = $id1;
			}
			elsif ( $seq->id =~ /^([^|]*)\|([^|]*)/ ) { # GENCODE sequences
				my ($id1) = $1;
				my ($id2) = $2;
				$sequence_id = $id1;
			}
			elsif ( $seq->desc =~ / transcript:([^\s]+)\s*/ ) { # Ensembl sequences
				$sequence_id = $1;
			}
			elsif ( $seq->id =~ /^(.*)$/ ) { # General case
				my ($id1) = $1;
				$sequence_id = $id1;
			}
			else {
				$sequence_id  = $seq->id if ( defined $seq->id );
				$sequence_id .= $seq->desc if ( defined $seq->desc );				
				$sequence_id  =~ s///mg;
			}
						
			if ( defined $sequence_id ) {
				$sequence_id =~ s/\s*//;
				if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// } # delete suffix in Ensembl ids
				if(exists $data->{$sequence_id}) {
					throw("Duplicated sequence: $sequence_id");
				}
				else {
					my ($sequence) = $seq->seq; # control short sequences
					my ($seq_len) = length($sequence);
					if ( $seq_len > 2 ) {
						$data->{$sequence_id} = $seq->seq;						
					}
					else {
						warning("Short sequence: $sequence_id");
					}
				}				
			}
		}		
	}
	return $data;
	
} # End _parse_inseq_transc

sub _parse_inseq_transl($)
{
	my ($inseqs) = @_;
	my ($data);

	$inseqs =~ s/\r//g; # delete jump line (REST services)

	my ($in);
	if ( defined $inseqs and (ref(\$inseqs) eq "SCALAR") and -e $inseqs and (-s $inseqs > 0) ) { # ensembl case
		$in = Bio::SeqIO->new(
							-file => $inseqs,
							-format => 'Fasta'
		);
	}
	elsif ( defined $inseqs and (ref(\$inseqs) eq "SCALAR") ) { # seq case
		my ($stringfh) = IO::String->new($inseqs);
		$in = Bio::SeqIO-> new(
								-fh     => $stringfh,
								-format => 'Fasta'
		);		
	}
	if ( defined $in ) {
		while ( my $seq = $in->next_seq() )
		{
			my ($sequence_id);
			if ( $seq->id =~ /^[sp|tr]\|([^|]*)\|([^|]*)/ ) { # UniProt sequences
				my ($id1) = $1;
				my ($id2) = $2;
				$sequence_id = $id1;
			}
			elsif ( $seq->id =~ /^gi\|[^|]*\|[^|]*\|([^|]*)/ ) { # RefSeq sequences
				my ($id1) = $1;
				$sequence_id = $id1;
			}
			elsif ( $seq->id =~ /^([^|]*)\|([^|]*)/ ) { # GENCODE sequences
				my ($id1) = $1;
				my ($id2) = $2;
				$sequence_id = $id1;
				if ( $id2 =~ /^ENSTR?[\d*]/ or $id2 =~ /^ENSMUST[\d*]/ ) { # From GENCODE (human and mouse), the second id is Ensembl Transcript id
					$sequence_id = $id2;
				}
			}
			elsif ( $seq->desc =~ / transcript:([^\s]+)\s*/ ) { # Ensembl sequences
				$sequence_id = $1;
			}
			elsif ( $seq->id =~ /^(.*)$/ ) { # General case
				my ($id1) = $1;
				$sequence_id = $id1;
			}			
			else {
				$sequence_id  = $seq->id if ( defined $seq->id );
				$sequence_id .= $seq->desc if ( defined $seq->desc );				
				$sequence_id  =~ s///mg;
			}
			
			if ( defined $sequence_id ) {
				$sequence_id =~ s/\s*//;
				if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// } # delete suffix in Ensembl ids
				if(exists $data->{$sequence_id}) {
					warning("Duplicated sequence: $sequence_id");
					#return undef;
				}
				else {
					my ($sequence) = $seq->seq; # control short sequences
					my ($seq_len) = length($sequence);
					if ( $seq_len > 2 ) {
						$data->{$sequence_id} = $seq->seq;						
					}
					else {
						warning("Short sequence: $sequence_id");
					}
				}				
			}
		}		
	}
	return $data;
	
} # End _parse_inseq_transl

sub _parse_seq_data($)
{
	my ($file) = @_;
	my ($data);
	my ($source) = 'wserver';

	if (-e $file and (-s $file > 0) ) {
		
		# get main id from file
		my ($dirname,$basename) = parse_file($file, 'transl.fa');
		if ( defined $basename ) {
			my (@dirs) = split('/',$dirname);
			my ($main_id) = $dirs[scalar(@dirs)-1];
			
			# get data from sequence file
			my ($in) = Bio::SeqIO->new(
								-file => $file,
								-format => 'Fasta'
			);
			while ( my $seq = $in->next_seq() )
			{
				if ( $seq->id =~ /([^|]*)/ )
				{
					my ($sequence_id) = $1;
					$source = 'sequence';
					
					my ($transc_id, $transl_id, $gene_id, $gene_name, $ccds_id) = (undef,undef,undef,undef,undef);
					my (@ids) = split('\|', $seq->id);
					if ( scalar(@ids) > 4 ) {
						$transc_id = $ids[0];
						$transl_id = $ids[1];
						$gene_id = $ids[2];
						$gene_name = $ids[3];
						$ccds_id = $ids[4];
						
						$main_id = $gene_id;
					}
					
					my ($gene_ids, $transc_ids) = (undef,undef);
					if ( $seq->desc =~ /gene_ids\>([^\s]+)/ ) { $gene_ids = $1 }
					if ( $seq->desc =~ /transc_ids\>([^\s]+)/ ) { $transc_ids = $1 }
					
					unless ( exists $data->{$main_id} ) {
						$data->{$main_id}->{'source'} = $source;
						$data->{$main_id}->{'transcripts'} = undef;
						if ( defined $gene_name ) { $data->{$main_id}->{'name'} = $gene_name }
						if ( defined $gene_ids )  { $data->{$main_id}->{'gene_ids'} = $gene_ids }
					}
					
					if(exists $data->{$main_id}->{'transcripts'}->{$sequence_id}) {
						throw("Duplicated sequence: $sequence_id");
					}
					else {
						my ($sequence) = $seq->seq;
						my ($seq_len) = length($sequence);

						if ( $seq_len > 2 ) { # control short sequences
							$data->{$main_id}->{'transcripts'}->{$sequence_id}->{'seq'} = $sequence;
							if ( defined $ccds_id and $ccds_id ne '' and $ccds_id ne '-' ) {
								$data->{$main_id}->{'transcripts'}->{$sequence_id}->{'ccdsid'} = $ccds_id
							}
							if ( defined $gene_name and $gene_name ne '' and $gene_name ne '-' ) {
								$data->{$main_id}->{'transcripts'}->{$sequence_id}->{'name'} = $gene_name;
							}
							if ( defined $transc_ids and $transc_ids ne '' and $transc_ids ne '-' ) {
								$data->{$main_id}->{'transcripts'}->{$sequence_id}->{'transc_ids'} = $transc_ids;
							}
						}
						else {
							warning("Short sequence: $sequence_id");
						}
					}						
				}
			}
		}
	}
	
	return $data;
	
} # End _parse_seq_data

# Create APPRIS::Transcript object from gencode data (GTF)
sub _fetch_transc_objects($$;$;$)
{
	my ($gene_id, $gene_features, $transc_seq, $transl_seq) = @_;
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;	
	
	# Scan transcripts
	while (my ($transcript_id, $transcript_features) = each(%{$gene_features}) )
	{
		my ($xref_identities);
		my ($sequence);
		my ($exons);

		# Create transcript object
		my ($transcript) = APPRIS::Transcript->new
		(
			-stable_id	=> $transcript_id,
			-chr		=> $transcript_features->{'chr'},
			-start		=> $transcript_features->{'start'},
			-end		=> $transcript_features->{'end'},
			-strand		=> $transcript_features->{'strand'},
			-biotype	=> $transcript_features->{'biotype'},
			-status		=> $transcript_features->{'status'},
			-source		=> $transcript_features->{'source'},
			-level		=> $transcript_features->{'level'},
			-version	=> $transcript_features->{'version'},
			-tsl		=> $transcript_features->{'tsl'},
			-tag		=> $transcript_features->{'tag'}
		);
			
		# Xref identifiers
		if ( defined $gene_id ) {
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $gene_id,
						-dbname			=> 'Gene_Id'
					)
			);
		}
		if ( exists $transcript_features->{'external_id'} and defined $transcript_features->{'external_id'} ) {
			$transcript->external_name($transcript_features->{'external_id'}); 
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $transcript_features->{'external_id'},
						-dbname			=> 'External_Id'
					)
			);
		}
		if ( exists $transcript_features->{'protein_id'} and defined $transcript_features->{'protein_id'} ) {
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $transcript_features->{'protein_id'},
						-dbname			=> 'Protein_Id'
					)
			);
		}
		if ( exists $transcript_features->{'ccdsid'} and defined $transcript_features->{'ccdsid'} ) {
			push(@{$xref_identities},
					APPRIS::XrefEntry->new
					(
						-id				=> $transcript_features->{'ccdsid'},
						-dbname			=> 'CCDS'
					)
			);
		}
		
		# Get transcript sequence
		if ( defined $transc_seq and
			 exists $transc_seq->{$transcript_id} and defined $transc_seq->{$transcript_id} ) {
				$sequence = $transc_seq->{$transcript_id};
		}
			
		# Get exon ids
		if ( exists $transcript_features->{'exons'} and scalar(@{$transcript_features->{'exons'}} > 0) )
		{
			my ($aux_exons);
			foreach my $exon (@{$transcript_features->{'exons'}})
			{
				my ($exon_id) = $transcript_id;
				if (exists $exon->{'exon_id'} and defined $exon->{'exon_id'}) {
					$exon_id = $exon->{'exon_id'};
				}
				push(@{$aux_exons},
					APPRIS::Exon->new
					(
						-stable_id	=> $exon_id,
						-start		=> $exon->{'start'},
						-end		=> $exon->{'end'},
						-strand		=> $exon->{'strand'},
					)
				);
			}
			$exons = sort_cds($aux_exons, $transcript_features->{'strand'}); # sort exons
		}

		# Add translation
		my ($translate) = _fetch_transl_objects($transcript_id, $transcript_features, $transl_seq);
		
		
		$transcript->xref_identify($xref_identities) if (defined $xref_identities);
		$transcript->sequence($sequence) if (defined $sequence);
		$transcript->exons($exons) if (defined $exons);
		$transcript->translate($translate) if (defined $translate);
			
		push(@{$transcripts}, $transcript) if (defined $transcript);
		$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts		
	}
	return ($transcripts,$index_transcripts);
	
} # End _fetch_transc_objects

# Create APPRIS::Translation object from gencode data
sub _fetch_transl_objects($$;$)
{
	my ($transcript_id, $transcript_features, $transl_seq) = @_;
	my ($translate);
	my ($protein_id);
	my ($sequence);
	my ($cds);
	my ($codons);	

	# Get protein id
	if ( exists $transcript_features->{'protein_id'} and defined $transcript_features->{'protein_id'} ) {
			$protein_id = $transcript_features->{'protein_id'};
	}

	# Get translate sequence
	if ( defined $transl_seq ) {
		 if ( exists $transl_seq->{$transcript_id} and defined $transl_seq->{$transcript_id} ) {
			$sequence = $transl_seq->{$transcript_id};
		}
		elsif ( defined $protein_id and exists $transl_seq->{$protein_id} and defined $transl_seq->{$protein_id} ) {
			$sequence = $transl_seq->{$protein_id};
		}
	}
		
	# Get cds
	if ( exists $transcript_features->{'cds'} and scalar(@{$transcript_features->{'cds'}} > 0) )
	{
		my ($aux_cds);
		foreach my $cds (@{$transcript_features->{'cds'}})
		{
			my ($exon_id);
			if (exists $cds->{'exon_id'} and defined $cds->{'exon_id'}) {
				$exon_id = $cds->{'exon_id'};
			}			
			push(@{$aux_cds},
				APPRIS::CDS->new
				(
					-start		=> $cds->{'start'},
					-end		=> $cds->{'end'},
					-strand		=> $cds->{'strand'},
					-phase		=> $cds->{'phase'},
					-stable_id	=> $exon_id,
				)
			);
		}
		$cds = sort_cds($aux_cds, $transcript_features->{'strand'}); # sort exons
	}	

	# Get codons
	if ( exists $transcript_features->{'codons'} and scalar(@{$transcript_features->{'codons'}} > 0) )
	{
		foreach my $codon (@{$transcript_features->{'codons'}})
		{
			push(@{$codons},
				APPRIS::Codon->new
				(
					-type		=> $codon->{'type'},
					-start		=> $codon->{'start'},
					-end		=> $codon->{'end'},
					-strand		=> $codon->{'strand'},
					-phase		=> $codon->{'phase'},
				)
			);
		}
	}
	
	if ( defined $sequence ) {
		
		# Create object
		$translate = APPRIS::Translation->new
		(
			-stable_id	=> $transcript_id,
		);			
		$translate->sequence($sequence);
		$translate->protein_id($protein_id) if (defined $protein_id);
		$translate->cds($cds) if (defined $cds);
		$translate->codons($codons) if (defined $codons);
		$translate->cds_sequence($translate);	
	}
	
	return $translate;
	
} # End _fetch_transl_objects

# Parser result file of INERTIA method
sub _parse_inertia_file($$\$)
{
	my ($type, $result, $ref_cutoffs) = @_;

	my ($transcript_id);	
	my (@results) = split( '\n', $result);
	
	foreach my $line (@results)
	{
		next if( $line =~ /^#/ ); # Skip comment line
		$line.="\n"; # Due we are spliting by '\n'
		
		if ( $line =~ /^>([^\t]+)\t+([^\n]+)\n+$/ )
		{
			$transcript_id = $1;
			my ($unusual_evolution) = $2;

			${$ref_cutoffs}->{$transcript_id}->{$type}->{'unusual_evolution'} = $unusual_evolution;			
			unless ( exists ${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} ) {
				${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} = $line;
			} else {
				${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} .= $line;
			}			
		}		
		elsif ( defined $transcript_id and ($line =~ /^\t+([^\:]+)\:([^\_]+)\:([^\t]+)\t([^\n]+)\n+$/) )
		{
			my ($start) = $1;
			my ($end) = $2;
			my ($strand) = $3;
			my ($exon_annotation) = $4;
		
			my ($exon_report) = {
							'start'					=> $start,
							'end'					=> $end,
							'strand'				=> $strand,
							'unusual_evolution'		=> $exon_annotation						
			};
			push( @{${$ref_cutoffs}->{$transcript_id}->{$type}->{'residues'}}, $exon_report );
			${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} .= $line;
		}		
	}	
} # End _parse_inertia_file

# Parser result file of Omega-INERTIA
sub _parse_omega_file($$\$)
{
	my ($type, $result, $ref_cutoffs) = @_;

	my (@results) = split( '\n', $result);
	
	foreach my $line (@results)
	{
		next if( $line =~/^#/ ); # Skip comment line
		$line.="\n"; # Due we are spliting by '\n'
				
		# omega_average omega_exon_id   start_exon      end_exon        strand_exon     difference_value        p_value st_desviation   exon_annotation transcript_list
		if ( $line =~ /^([^\t]+)\t+([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\t]+)\t+([^\n]+)\n+$/ )
		{
			my ($omega_mean) = $1;
			my ($omega_exon_id) = $2;
			my ($start) = $3;
			my ($end) = $4;
			my ($strand) = $5;
			my ($d_value) = $6;
			my ($p_value) = $7;
			my ($st_desviation) = $8;
			my ($exon_annotation) = $9;
			my ($exon_transcrits_list) = $10;
		
			# Get the trasncipt with omega exons
			my (@exon_transcrits);
			if ( $exon_transcrits_list ne 'NULL' ) {
				@exon_transcrits = split(';',$exon_transcrits_list);			
			}
			
			foreach my $transcript_id (@exon_transcrits)
			{
				my ($omega_exon_report) = {
							'omega_exon_id'			=> $omega_exon_id,
							'start'					=> $start,
							'end'					=> $end,
							'strand'				=> $strand,
							'omega_mean'			=> $omega_mean,
							'st_deviation'			=> $st_desviation,
							'difference_value'		=> $d_value,
							'p_value'				=> $p_value,
							'unusual_evolution'		=> $exon_annotation						
				};
				push( @{${$ref_cutoffs}->{$transcript_id}->{$type}->{'residues'}}, $omega_exon_report );
				unless ( exists ${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} ) {
					${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} = $line;
				} else {
					${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} .= $line;
				}
			}
		}

		# # omega_average omega_exon_id   start_exon      end_exon        difference_value        p_value st_desviation   exon_annotation transcript_list
		if ( $line =~ /^>([^\t]+)\t+([^\n]+)\n+$/ )
		{
			my ($transcript_id) = $1;
			my ($unusual_evolution) = $2;

			${$ref_cutoffs}->{$transcript_id}->{$type}->{'unusual_evolution'} = $unusual_evolution;			
			${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} .= "----------------------------------------------------------------------\n";
			${$ref_cutoffs}->{$transcript_id}->{$type}->{'result'} .= $line;

			${$ref_cutoffs}->{$transcript_id}->{$type}->{'omega_average'} = 0; # DEPRECATED
			${$ref_cutoffs}->{$transcript_id}->{$type}->{'omega_st_desviation'} = 0; # DEPRECATED			
		}		
	}
} # End _parse_omega_file

1;