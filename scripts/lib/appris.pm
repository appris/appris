=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

appris - common process

=head1 SYNOPSIS

=head1 DESCRIPTION


=head1 METHODS

=cut

package appris;

use strict;
use warnings;
use Bio::SeqIO;
use File::Temp;
use Config::IniFiles;
use MIME::Lite;

use APPRIS::Parser qw( parse_gencode parse_infiles );
use APPRIS::Utils::File qw( getTotalStringFromFile printStringIntoFile );
use APPRIS::Utils::Argument qw( rearrange );
use APPRIS::Utils::Exception qw( info throw warning deprecate );

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	create_gene_list
	get_gene_list
	get_gene_list_from_seq
	retrieve_gene_list
	create_appris_input
	create_ensembl_input
	create_appris_seqinput
	reduce_input
	create_gencode_data
	create_indata
	create_seqdata
	send_email
);


=head2 create_gene_list

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_gene_list
{
	my ( $itype, $id,  $list,  $pos,  $gdata, $gtransc, $gtransl, $econf, $ever, $spe) = rearrange( 
		['type', 'id', 'list', 'pos', 'gdata','gtransc','gtransl','econf','ever', 'spe' ],
	@_ );
	my ($report, $data) = (undef, undef);
	
	if ( $itype =~ /datafile/ and defined $gdata )
	{
		my ($data_fh);
		if ( defined $list ) {
			$report = get_gene_list($list);
			$data_fh = reduce_input($gdata, undef, $report);
			if ( UNIVERSAL::isa($data_fh,'File::Temp') ) {
				$gdata = $data_fh->filename;
			}
		}		
		elsif ( defined $pos ) {
			$data_fh = reduce_input($gdata, $pos);
			if ( UNIVERSAL::isa($data_fh,'File::Temp') ) {
				$gdata = $data_fh->filename;
			}
		}
		$data = create_indata($gdata, $gtransc, $gtransl);
						
		# delete tmp file
		if ( defined $data_fh and UNIVERSAL::isa($data_fh,'File::Temp') ) {
			$data_fh->unlink_on_destroy(1);
		}
	}	
	elsif ( $itype =~ /ensembl/ and defined $econf and defined $ever and defined $spe )
	{
		if ( defined $id ) {
			$report->{$id} = 1;
		}
		elsif ( defined $list ) {
			$report = get_gene_list($list);
		}
		elsif ( defined $pos ) {
			$report = create_ensembl_input($pos, $econf, $ever, $spe);			
		}
		else {
			$report = create_ensembl_input(undef, $econf, $ever, $spe);			
		}
	}
	elsif ( $itype =~ /sequence/ and defined $gtransl )
	{
		$data = create_seqdata($gtransl);
	}
	
	File::Temp::cleanup();
	
	return ($report, $data);
}

=head2 get_gene_list

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub get_gene_list($)
{
	my ($file) = @_;
	my ($list);
	
	my ($genes) = getTotalStringFromFile($file);
	foreach my $gene_id (@{$genes}) {
		$gene_id =~ s/\s*//mg;
		$list->{$gene_id} = 1 if ( $gene_id ne '');
	}
	return $list;
} # end get_gene_list

=head2 get_gene_list_from_seq

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub get_gene_list_from_seq($)
{
	my ($file) = @_;
	my ($list);
	my ($in) = Bio::SeqIO->new(
						-file => $file,
						-format => 'Fasta'
	);
	if ( defined $in ) {
		while ( my $seq = $in->next_seq() )
		{
			my ($id);
			# At the moment, only for UniProt/neXtProt cases
			if ( $seq->id =~ /^(sp|tr)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences
				$id = $2;
			}
			elsif ( $seq->id =~ /^nxp:([^\s]*)/ ) { # neXtProt sequences
				$id = $1;
			}
			if ( defined $id ) {
				my ($gene_id) = ( $id =~ /([^\-]*)/ ) ? $1 : $id;
				$list->{$gene_id} = 1;				
			}
		}		
	}
	return $list;
} # end get_gene_list

=head2 retrieve_gene_list

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub retrieve_gene_list($;$)
{
	my ($data_file, $position) = @_;
	my ($list);
		
	# customize data for a position
	if ( defined $position ) {
		foreach my $pos ( split(',', $position) ) {
			if ( defined $pos ) {
				eval {
					my ($cmd) =	"awk '{if( (\$1==\"$pos\" || \$1==\"chr$pos\") && \$3==\"gene\") {print \$10}}' $data_file | sed 's/[\"|;]//g' | sed 's/\.[0-9]*\$//'";
					info("** script: $cmd\n");
					my (@aux_list) = `$cmd`;
					foreach my $gene_id (@aux_list) {
						if ( defined $gene_id and ($gene_id ne '') ) {
							$gene_id =~ s/\s*//mg;
							$list->{$gene_id} = 1 if ( $gene_id ne '');
						}
					}
				};
				throw("creating input for $pos") if($@);
			}		
		}		
	}
	return $list;
} # end retrieve_gene_list

=head2 create_appris_input

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_appris_input($$)
{
	my ($gene,$in_files) = @_;
	my ($create) = undef;
	my ($data_cont) = '';
	my ($pdata_cont) = '';
	my ($transc_cont) = '';
	my ($transl_cont) = '';
	my ($cdsseq_cont) = '';
	
	# gene vars
	my ($chr) = $gene->chromosome;
	my ($gene_id) = $gene->stable_id;
	my ($g_start) = $gene->start;
	my ($g_end) = $gene->end;
	my ($g_strand) = $gene->strand;
	my ($g_phase) = '.';
	my ($gene_name) = $gene->external_name ? $gene->external_name : $gene_id;
	my ($g_source) = '';
	if ( defined $gene->source ) {
		if ( (lc($gene->source) =~ /ensembl/) or (lc($gene->source) =~ /havana/) or (lc($gene->source) =~ /gencode/) ) {
			$g_source = 'GENCODE';
		}
		elsif ( (lc($gene->source) =~ /refseq/) or (lc($gene->source) =~ /bestrefseq/) or (lc($gene->source) =~ /curated genomic/) or (lc($gene->source) =~ /gnomon/) ) {
			$g_source = 'REFSEQ';
		}			
	}
	
	# get gene annots
	my ($data_gene_attrs) = "gene_id \"$gene_id\"; gene_name \"$gene_name\"";
	if ( $gene->biotype ) { $data_gene_attrs .= '; gene_type "'.$gene->biotype.'"' }
	$data_cont .=	$chr."\t".
					$g_source."\t".
					'gene'."\t".
					$g_start."\t".
					$g_end."\t".
					'.'."\t".
					$g_strand."\t".
					$g_phase."\t".
					$data_gene_attrs."\n";

	# scan transcript/translation/cds seq/cds info
	if ( $gene->transcripts ) {
		foreach my $transcript (@{$gene->transcripts}) {		
			my ($transcript_id) = $transcript->stable_id;
			my ($transcript_eid) = $transcript_id;
			#if ( $transcript->version ) {
			#	$transcript_eid = $transcript_id.'.'.$transcript->version;
			#}
			my ($t_start) = $transcript->start;
			my ($t_end) = $transcript->end;
			my ($t_strand) = $transcript->strand;
			my ($t_phase) = '.';
			my ($transcript_name) = $transcript->external_name ? $transcript->external_name : $transcript_id;
			my ($data_transc_attrs) = "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; transcript_name \"$transcript_name\"";
			my ($ccds_id) = '-';				
			if ( $transcript->xref_identify ) {
				foreach my $xref_identify (@{$transcript->xref_identify}) {								
					if ($xref_identify->dbname eq 'CCDS') { $ccds_id = $xref_identify->id; last; }
				}
				if ( $ccds_id ne '-' ) { $data_transc_attrs .= '; ccds_id "'.$ccds_id.'"'; }
			}
			if ( $transcript->biotype ) { $data_transc_attrs .= '; transcript_type "'.$transcript->biotype.'"' }
			if ( $transcript->tsl ) { $data_transc_attrs .= '; tsl "'.$transcript->tsl.'"' }
			if ( $transcript->tag ) { $data_transc_attrs .= '; tag "'.$transcript->tag.'"' }
			$data_cont .=	$chr."\t".
							$g_source."\t".
							'transcript'."\t".
							$t_start."\t".
							$t_end."\t".
							'.'."\t".
							$t_strand."\t".
							$t_phase."\t".
							$data_transc_attrs."\n";
							
			# get transcript seq/exon
			if ( $transcript->sequence ) {
				my ($seq) = $transcript->sequence;
				my ($len) = length($transcript->sequence);
				$transc_cont .= ">$transcript_eid|$gene_id|$gene_name|$len\n";
				$transc_cont .= $seq."\n";
				if ( $transcript->exons ) {
					foreach my $exon (@{$transcript->exons}) {
						my ($exon_start) = $exon->start;
						my ($exon_end) = $exon->end;
						my ($exon_strand) = $exon->strand;
						my ($exon_phase) = '.';
						my ($exon_id) = defined $exon->stable_id ? $exon->stable_id : '-';						
						$data_cont .=	$chr."\t".
										$g_source."\t".
										'exon'."\t".
										$exon_start."\t".
										$exon_end."\t".
										'.'."\t".
										$exon_strand."\t".
										$exon_phase."\t".
										"gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; exon_id \"$exon_id\";\n";			
					}
				}
			}
			# get translation seq/cds/cds_seq/codons
			if ( $transcript->translate ) {
				my ($translate) = $transcript->translate;
				if ( $translate->sequence ) {
					my ($seq) = $translate->sequence;
					my ($len) = length($translate->sequence);
					my ($translate_id) = $transcript_eid;
					$translate_id = $translate->protein_id if ( defined $translate->protein_id );
					# mask short sequences
					#if ( $len <= 2 ) { $seq .= 'X'; }
					$transl_cont .= ">$transcript_eid|$translate_id|$gene_id|$gene_name|$ccds_id|$len\n";
					$transl_cont .= $seq."\n";					
				}
				if ( $translate->cds and $translate->cds_sequence ) {
					my ($exons) = $transcript->exons;					
					for (my $icds = 0; $icds < scalar(@{$translate->cds}); $icds++) {
						my ($cds) = $translate->cds->[$icds];
						my ($exon_id) = (defined $exons->[$icds] and $exons->[$icds]->stable_id) ? $exons->[$icds]->stable_id : '-';
						my ($pro_cds) = $translate->cds_sequence->[$icds];
	
						my ($cds_start) = $cds->start;
						my ($cds_end) = $cds->end;
						my ($cds_strand) = $cds->strand;
						my ($cds_phase) = $cds->phase;
						
						my ($pro_cds_start) = $pro_cds->start;
						my ($pro_cds_end) = $pro_cds->end;
						my ($pro_cds_end_phase) = $pro_cds->end_phase;
						my ($pro_cds_seq) = $pro_cds->sequence;
							
						# delete the residue that is shared by two CDS						
						#if (defined $pro_cds_end_phase and $pro_cds_end_phase != 0) { $pro_cds_seq =~ s/\w{1}$//; }
							
						if (defined $pro_cds_seq and $pro_cds_seq ne '') {
							my ($len) = length($pro_cds_seq);
							$cdsseq_cont .= ">$transcript_eid|$gene_id|$gene_name|$len|$exon_id|$chr|$cds_start|$cds_end|$cds_strand|$cds_phase\n";
							$cdsseq_cont .= $pro_cds_seq."\n";							
						}													
						if (defined $cds_start and defined $cds_end) {					
							$data_cont .=	$chr."\t".
											$g_source."\t".
											'CDS'."\t".
											$cds_start."\t".
											$cds_end."\t".
											'.'."\t".
											$cds_strand."\t".
											$cds_phase."\t".
											"gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; cds_id \"$exon_id\";\n";
						}						
						if (defined $pro_cds_start and defined $pro_cds_end) {
							$pdata_cont .=	'SEQ'."\t".
											$g_source."\t".
											'Protein'."\t".
											$pro_cds_start."\t".
											$pro_cds_end."\t".
											'.'."\t".
											'.'."\t".
											$pro_cds_end_phase."\t".
											"ID=$exon_id;Parent=$transcript_eid;Gene=$gene_id;Note=cds_coord>$cds_start-$cds_end:$cds_strand\n";					
						}
					}
				}
				if ( $translate->codons ) {
					foreach my $codon (@{$translate->codons}) {
						my ($codon_type) = $codon->type.'_codon';
						my ($codon_start) = $codon->start;
						my ($codon_end) = $codon->end;
						my ($codon_strand) = $codon->strand;
						my ($codon_phase) = $codon->phase;
						$data_cont .=	$chr."\t".
										$g_source."\t".
										$codon_type."\t".
										$codon_start."\t".
										$codon_end."\t".
										'.'."\t".
										$codon_strand."\t".
										$codon_phase."\t".
										"gene_id \"$gene_id\"; transcript_id \"$transcript_id\";\n";			
					}
				}
			}
		}
	}

	# create files
	if ( $data_cont ne '' ) {
		my ($output_file) = $in_files->{'data'};		
		my ($printing_file_log) = printStringIntoFile($data_cont, $output_file);
		throw("creating $output_file file") unless ( defined $printing_file_log );
	}
	if ( $pdata_cont ne '' ) {
		my ($output_file) = $in_files->{'pdata'};
		my ($printing_file_log) = printStringIntoFile($pdata_cont, $output_file);
		throw("creating $output_file file") unless ( defined $printing_file_log );
	}
	if ( $transc_cont ne '' ) {
		my ($output_file) = $in_files->{'transc'};
		my ($printing_file_log) = printStringIntoFile($transc_cont, $output_file);
		throw("creating $output_file file") unless ( defined $printing_file_log );
	}
	if ( $transl_cont ne '' ) {
		my ($output_file) = $in_files->{'transl'};
		my ($printing_file_log) = printStringIntoFile($transl_cont, $output_file);
		throw("creating $output_file file") unless ( defined $printing_file_log );
	}
	if ( $cdsseq_cont ne '' ) {
		my ($output_file) = $in_files->{'cdsseq'};
		my ($printing_file_log) = printStringIntoFile($cdsseq_cont, $output_file);
		throw("creating $output_file file") unless ( defined $printing_file_log );
	}
	
	# determine if appris has to run
	if ( ($data_cont ne '') and ($pdata_cont ne '') and ($transl_cont ne '') and ($cdsseq_cont ne '') ) {
		$create = 1;
	}
	
	return $create;
	
} # end create_appris_input

=head2 create_ensembl_input

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_ensembl_input($$$$)
{
	my ($position, $conf_file, $e_version, $species) = @_;
	my ($gene_list);
	
	# create ensembl registry from default config file
	my ($cfg) = new Config::IniFiles( -file => $conf_file );
	my ($e_core_param) = {
			'-host'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host'),
			'-user'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user'),
			'-pass'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'pass'),
			'-verbose'    => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose'),
			'-db_version' => $e_version,
			'-species'    => $species,
	};
	info(
		"\t-host        => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host')."\n".
		"\t-user        => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user')."\n".
		"\t-pass        => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'pass')."\n".
		"\t-verbose     => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose')."\n".
		"\t-db_version  => ".$e_version."\n".										
		"\t-species     => ".$species."\n"
	);
	require Bio::EnsEMBL::Registry;
	my ($registry) = 'Bio::EnsEMBL::Registry';
	eval {
		$registry->load_registry_from_db(%{$e_core_param});
	};
	throw("can not load ensembl registry: $!\n") if $@;
	
	# create slice adaptor
	my ($slice_adaptor) = $registry->get_adaptor( $species, 'Core', 'Slice' );	
	
	# create ensembl input for each given position
	my (@slices);
	if ( defined $position ) {
		@slices = split(',', $position);
	}
	else {
 		foreach my $slice ( @{$slice_adaptor->fetch_all('chromosome')} ) { 			
 			push(@slices, $slice->seq_region_name);
 		}		
	}
	
	# get gene list depending on given position
	foreach my $pos ( @slices ) {
		if ( defined $pos ) {
			$pos =~ s/^chr//;		
			my ($slice) = $slice_adaptor->fetch_by_region( 'chromosome', $pos );
			if ( defined $slice ) {
				my ($genes) = $slice->get_all_Genes();
				if ( defined $genes ) {
					while ( my $gene = shift @{$genes} ) {
						my ($gene_id) = $gene->stable_id;
						my ($gene_eid) = $gene_id;
						#if ( $gene->version ) {
						#	$gene_eid = $gene_id.'.'.$gene->version;
						#}
						$gene_list->{$gene_eid} = $gene->biotype();
					}					
				}
			}
		}		
	}
	
	return $gene_list;
	
} # end create_ensembl_input

=head2 create_appris_seqinput

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_appris_seqinput($$)
{
	my ($gene,$in_files) = @_;
	my ($create) = undef;
	my ($transl_cont) = '';
	
	# gene vars
	my ($gene_id) = $gene->{'id'};
	my ($gene_name) = ( exists $gene->{'name'} ) ? $gene->{'name'} : $gene_id;
	
	# scan translation info
	if ( exists $gene->{'varsplic'} ) {
		foreach my $isof_id (sort( keys(%{$gene->{'varsplic'}}) ) ) {
			my ($isof) = $gene->{'varsplic'}->{$isof_id};
			if ( exists $isof->{'seq'} ) {
				my ($ccds_id) = ( exists $isof->{'ccds'} ) ? $isof->{'ccds'} : '-';
				my ($seq) = $isof->{'seq'};
				my ($len) = length($seq);
				$transl_cont .= ">$isof_id|$isof_id|$gene_id|$gene_name|$ccds_id|$len\n";
				$transl_cont .= $seq."\n";				
			}
		}
	}

	# create files
	if ( $transl_cont ne '' ) {
		my ($output_file) = $in_files->{'transl'};
		my ($printing_file_log) = printStringIntoFile($transl_cont, $output_file);
		throw("creating $output_file file") unless ( defined $printing_file_log );
	}
	
	# determine if appris has to run
	if ( $transl_cont ne '' ) {
		$create = 1;
	}
	
	return $create;
	
} # end create_appris_seqinput

=head2 reduce_input

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub reduce_input($;$;$)
{
	my ($data_file, $position, $gene_list) = @_;
	my ($data_tmpfile) = File::Temp->new( UNLINK => 0, SUFFIX => '.appris.dat' );
	my ($data_tmpfilename) = $data_tmpfile->filename;
		
	# customize data for a position
	if ( defined $position ) {
		foreach my $pos ( split(',', $position) ) {
			if ( defined $pos ) {
				eval {
					my ($cmd) =	"awk '{if(\$1==\"$pos\" || \$1==\"chr$pos\") {print \$0}}' $data_file >> $data_tmpfilename";
					info("** script: $cmd\n");
					my (@cmd_out) = `$cmd`;
				};
				throw("creating input for $pos") if($@);
			}		
		}		
	}
	# get customized data for a gene list
	elsif ( defined $gene_list ) {
		my ($g_cond) = '';	
		foreach my $g ( keys(%{$gene_list}) ) {
			$g_cond .= ' $10 ~ /'.$g.'/ ||'; # for rel7 version		
		}
		if ( $g_cond ne '' ) {
	    	$g_cond =~ s/\|\|$//;
	    	$g_cond = '('.$g_cond.')';			
			eval {
				my ($cmd) = "awk '{if(\$9==\"gene_id\" && $g_cond) {print \$0}}' $data_file >> $data_tmpfilename";
				info("** script: $cmd\n");
				my (@cmd_out) = `$cmd`;
			};
			throw("creating gencode data for genelist") if($@);			
		}		
	}	
	
	return $data_tmpfile;
	
} # end reduce_input

=head2 create_gencode_data

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_gencode_data($;$;$)
{
	my ($data_file, $transcripts_file, $translations_file) = @_;
	
	my ($data) = parse_gencode($data_file, $transcripts_file, $translations_file);
	unless ( defined $data ) {
		throw("can not create gencode object\n");
	}
	
	return $data;
	
} # end create_gencode_data

=head2 create_indata

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_indata($;$;$)
{
	my ($data_file, $transcripts_file, $translations_file) = @_;
	
	my ($data) = parse_infiles($data_file, $transcripts_file, $translations_file);
	unless ( defined $data ) {
		throw("can not create gencode object\n");
	}
	
	return $data;
	
} # end create_indata

=head2 create_seqdata

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_seqdata($)
{
	my ($file) = @_;
	my ($data);

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() ) {
			my ($s_id) = $seq->id;
			my ($s_desc) = $seq->desc;
			my ($s_seq) = $seq->seq;
			my ($isof_id);
			my ($transl_id);
			my ($gene_id);
			my ($gene_name);
			my ($ccds_id);
			my ($seq_length);
			my ($a_genes_id);
			my ($a_transl_id);
			my ($a_genes_name);
			my ($a_ccds_id);
			
			# At the moment, only for UniProt/neXtProt cases
			if ( $s_id =~ /^(sp|tr)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences
				$isof_id = $2;
				my (@desc) = split('GN=', $s_desc);
				if ( scalar(@desc) >= 2 ) {
					if ( $desc[1] =~ /([^\s]*)/ ) { $gene_name = $1 }					
				}				
			}
			elsif ( $s_id =~ /^(sp_a|tr_a)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences with extra values
				$isof_id = $2;
				my ($name) = $3;
				$gene_id = $4;
				$gene_name = $5;
				$ccds_id = $6;
				$seq_length = $7;
			}
			elsif ( $s_id =~ /^ge_a\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^\$]*)$/ ) { # GENCODE sequences with extra values
				my ($pep_id) = $1;
				$isof_id = $2;
				$gene_id = $3;
				$gene_name = $4;
				$ccds_id = $5;
				$seq_length = $6;
			}
			elsif ( $s_id =~ /^en_a\|([^\s]*)/ ) { # ENSEMBL sequences with extra values
				my ($pep_id) = $1;
				if ( $s_desc =~ /gene:([^\s]+).*transcript:([^\s]+).*gene_symbol:([^\s]+).*ccds:([^\s]+)/ ) {
					$gene_id = $1;
					$isof_id = $2;
					$gene_name = $3;
					$ccds_id = $4;
				}
			}
			elsif ( $s_id =~ /^gi_a\|[^|]*\|[^|]*\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^\s]*)\s*([^\$]*)$/ ) { # RefSeq sequences with extra values
				my ($pep_id) = $1;
				$isof_id = $2;
				$gene_id = $3;
				$gene_name = $4;
				$ccds_id = $5;
				$seq_length = $6;
			}
			elsif ( $s_id =~ /^nxp:([^\s]*)/ ) { # neXtProt sequences
				$isof_id = $1;
				my (@desc) = split('Gname=', $s_desc);
				if ( scalar(@desc) >= 2 ) {
					if ( $desc[1] =~ /([^\s]*)/ ) { $gene_name = $1 }					
				}				
			}
			elsif ( $s_id =~ /^appris\|([^\s]*)/ ) { # APPRIS sequences FASTA file
				$isof_id = $1;
				if ( $s_desc =~ /xref_genes:([^\:]+)\:([^\/]+).*transc:([^\s]+).*genes:([^\s]+).*gene_names:([^\s]+).*ccds:([^\s]+)/ ) {
					$gene_id = $1;
					$gene_name = $2;
					$a_transl_id = $3;
					$a_genes_id = $4;
					$a_genes_name = $5;
					$a_ccds_id = $6;	
				}
			}
			if ( defined $isof_id ) {
				$gene_id =~ s/\.[0-9]+$//g; $isof_id =~ s/\.[0-9]+$//g; # delete version suffix
				unless ( defined $gene_id ) {
					$gene_id = ( $isof_id =~ /([^\-]*)/ ) ? $1 : $isof_id;					
				}
				unless ( exists $data->{$gene_id} ) {
					$data->{$gene_id} = {
						'id'		=> $gene_id,
						'varsplic'	=> {}
					};													
					if ( defined $gene_name and $gene_name ne '' and $gene_name ne '-' ) {
						$data->{$gene_id}->{'name'} = $gene_name;
					}
				}
				$data->{$gene_id}->{'varsplic'}->{$isof_id} = {
					'desc'		=> $s_desc,
					'seq' 		=> $s_seq
				};
				if ( defined $ccds_id and $ccds_id ne '' and $ccds_id ne '-' ) {
					$data->{$gene_id}->{'varsplic'}->{$isof_id}->{'ccds'} = $ccds_id;
				}
				if ( defined $a_transl_id and $a_transl_id ne '' and $a_transl_id ne '-' ) {
					$data->{$gene_id}->{'varsplic'}->{$isof_id}->{'a_transl'} = $a_transl_id;
				}
				if ( defined $a_genes_id and $a_genes_id ne '' and $a_genes_id ne '-' ) {
					$data->{$gene_id}->{'varsplic'}->{$isof_id}->{'a_genes'} = $a_genes_id;
				}
				if ( defined $a_ccds_id and $a_ccds_id ne '' and $a_ccds_id ne '-' ) {
					$data->{$gene_id}->{'varsplic'}->{$isof_id}->{'a_ccds'} = $a_ccds_id;
				}
			}
		}		
	}
	return $data;
	
} # End create_seqdata

=head2 send_email

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub send_email($$$$$$)
{
	my ($econf, $email, $wserver, $species, $methods, $outpath) = @_;
	
	# Attach results if webserver
	my ($gzip_file);
	unless ( defined $wserver ) { # local execution
		my (@inpath_n) = split('/', $outpath);		
		if ( scalar(@inpath_n) > 0 ) {
			my ($num) = scalar(@inpath_n);
			$wserver = $inpath_n[$num-1];
		}
	}
	else
	{ # webserver
		# create tmp dir for email sending
		my ($email_tmp_dir) = $outpath."/email";
		eval {
			my ($cmd) =	"mkdir $email_tmp_dir";
			my (@cmd_out) = `$cmd`;
		};
		throw("creating tmp dir $email_tmp_dir") if($@);
		
		# create tmp files from methods
		foreach my $method ( split(',', $methods) ) {
			eval {
				my ($cmd) =	"find $outpath -type f -name '*.$method' ! -size 0 -exec cp -v {} $email_tmp_dir/$wserver.$method \\;";
				my (@cmd_out) = `$cmd`;
			};
			throw("creating tmp files $method in $email_tmp_dir") if($@);
		}
		
		# create gzip file
		my ($gzip_file) = "$email_tmp_dir/$wserver.tar.gz";
		eval {
			my ($cmd) =	"cd $email_tmp_dir && tar -cf - * | gzip -9 > $gzip_file";
			my (@cmd_out) = `$cmd`;
		};
		throw("compressing file $gzip_file in $email_tmp_dir") if($@);
	}
	
	my ($cfg) = new Config::IniFiles( -file => $econf );
	my ($cfg_from) = $cfg->val('APPRIS_EMAIL', 'from');
	my ($cfg_subject) = $cfg->val('APPRIS_EMAIL', 'subject');
	
	# send email
	my ($subject,$content);
	if ( defined $ENV{APPRIS_WSERVER_REPORT_URL} ) {
		$subject = $cfg_subject." Your jobid $wserver has finished";
		$content = "<p>The results for your query $wserver are now available at ";	
		my ($rst_url) = $ENV{APPRIS_WSERVER_REPORT_URL} . '/' . $wserver;
		$content	 .= "<a href='$rst_url'>$rst_url</a></p>";
		$content	 .= "</p>";
		#$content	 .= "<p>Please remember that your results will be stored in our server only for 3 weeks, than will be deleted !!</p>";		
	}
	else {
		$subject = $cfg_subject." Your jobid $species/$wserver has finished";
		$content = "<p>The results for your query $species/$wserver are now available.";	
	}
	_semail($cfg_from, $email, $subject, $content, $gzip_file);
	
	# rm email dir
	#if ( -d $email_tmp_dir ) {
	#	$logger->info("-- delete email dir\n");
	#	eval {
	#		my ($cmd) =	"rm -rf $email_tmp_dir";
	#		$logger->debug("** script: $cmd\n");
	#		my (@cmd_out) = `$cmd`;
	#	};
	#	$logger->error("deleting file\n") if($@);
	#}
					
} # end send_email

sub _semail($$$$;$)
{
	my ($from, $to, $subject, $content, $files) = @_;
	my ($msg) = MIME::Lite->new(
        From     	=> $from,
        To 			=> $to,
        Subject 	=> $subject,
		Type		=> 'multipart/mixed',
	);
	$msg->attach(
		Type    	=> 'text/html',
		Encoding 	=> 'quoted-printable',
		Data    	=> $content,
	);
	if ( defined $files ) {
		foreach my $file (split(';', $files)) {
			$msg->attach(
				Type		 => 'x-gzip',
				Path    	 => $file,
				Disposition	=> 'attachment'
			);		
		}		
	}
 	eval { $msg->send }; die "MIME::Lite->send failed: $@\n" if $@;
 	
} # end _semail


1;
