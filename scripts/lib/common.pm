package common;

use strict;

use FindBin;
use File::Temp;
use Bio::SeqIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( getTotalStringFromFile );
use APPRIS::Utils::Exception qw( info throw warning deprecate );

use lib "$ENV{APPRIS_PROGRAMS_SRC_DIR}/appris";
use appris;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	get_seq_report
	get_longest_seq
	get_main_report
	get_label_report
	get_prin_report
	src_appris_decision
);

###################
# Global variable #
###################
use vars qw(
	$CONFIG_FILE
);

$CONFIG_FILE = $ENV{APPRIS_CODE_CONF_DIR}."/pipeline.ini";

=head2 get_seq_report

  Arg[1]      : String $seq_file  - File path of translation sequences
  Description : Retrieves report of sequences.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub get_seq_report($)
{
	my ($file) = @_;
	my ($report);

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() )
		{
			if ( $seq->id=~/^([^|]*)\|([^|]*)\|([^|]*)/ )
			{
				my ($prot_id) = $1;
				my ($trans_id) = $2;
				my ($gene_id) = $3;
				# delete suffix in Ensembl ids
				if ( $prot_id =~ /^ENS/ ) { $prot_id =~ s/\.\d*$// }
				if ( $trans_id =~ /^ENS/ ) { $trans_id =~ s/\.\d*$// }
				if ( $gene_id =~ /^ENS/ ) { $gene_id =~ s/\.\d*$// }
				if(exists $report->{$gene_id}->{'transcripts'}->{$trans_id}) {
					if ( $report->{$gene_id}->{'transcripts'}->{$trans_id} eq $seq->seq ) {
						warning("Duplicated sequence: $trans_id");	
					}
					else {
						throw("Duplicated id but different sequence!!!: $trans_id");
					}
				}
				else {
					$report->{$gene_id}->{'transcripts'}->{$trans_id} = $seq->seq; 
				}
			}
		}		
	}
	return $report;
	
} # end get_seq_report

=head2 get_longest_seq

  Arg[1]      : String $seq_file  - File path of translation sequences
  Description : Retrieves report of sequences.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub get_longest_seq($)
{
	my ($seq_report) = @_;
	my ($report);
	
	while ( my ($gene_id, $g_report) = each (%{$seq_report}) ) {
		if ( exists $g_report->{'transcripts'} ) {
			my ($aux_report);
			while ( my ($transcript_id, $translation_seq) = each (%{$g_report->{'transcripts'}}) ) {
				$report->{$gene_id}->{'longest_len'}->{$transcript_id} = '0';
				my ($translation_length) = length($translation_seq);	
				if (defined $aux_report) {
					if ( $translation_length > $aux_report->{'translation_length'} ) { # maximum length
						# discard the last assignments
						foreach my $t (@{$aux_report->{'transcript_id'}}) {
							$report->{$gene_id}->{'longest_len'}->{$t}			= '0';								
						}
						$report->{$gene_id}->{'longest_len'}->{$transcript_id}	= '1'; #assign the new value
						$aux_report->{'translation_length'}						= $translation_length;
						$aux_report->{'transcript_id'}							= undef;
						push(@{$aux_report->{'transcript_id'}}, $transcript_id);
					}
					elsif ( $translation_length == $aux_report->{'translation_length'} ) { # equal length
						$report->{$gene_id}->{'longest_len'}->{$transcript_id}	= '1';
						$aux_report->{'translation_length'}						= $translation_length;
						push(@{$aux_report->{'transcript_id'}}, $transcript_id);
					}					
				}
				else { # init
					$report->{$gene_id}->{'longest_len'}->{$transcript_id}	= '1';
					$aux_report->{'translation_length'}						= $translation_length;
					push(@{$aux_report->{'transcript_id'}}, $transcript_id);
				}
			}
		}
    }		
	return $report;
		
} # end get_longest_seq

=head2 get_main_report

  Arg[1]      : String $main_file - File path of main annotations
  Arg[2]      : String $seq_file  - File path of translation sequences
  Arg[2]      : String $data_file - File path of data (Optional)
  Description : Retrieves report of main annotations.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub get_main_report($$;$)
{
	my ($main_file, $seq_file, $data_file) = @_;
	my ($report);
	
	# get content of appris annots
	my ($content) = getTotalStringFromFile($main_file);
	
	# get report from gencode
	my ($data_cont);
	if ( defined $data_file ) {
		$data_cont = APPRIS::Parser::_parse_indata($data_file);
	}
	
	# create sequence report
	my ($seq_report) = get_seq_report($seq_file);
		
	#gene_id	gene_name	transcript_id	status	biotype	no_codons ccds_id tsl aa_len
	#
	#fun_res
	#con_struct
	#vert_signal
	#dom_signal
	#tmh_signal
	#pep_mit_signal
	#mit_signal
	#u_evol
	#num_peptides
	#appris_score
	foreach my $input_line (@{$content})
	{
		$input_line=~s/\n*$//;
		my(@split_line)=split("\t", $input_line);
		next if(scalar(@split_line)<=7);

		my($gene_id)=$split_line[0];
		my($gene_name)=$split_line[1];
		my($transcript_id)=$split_line[2];
		my($translation_id)=$split_line[3];
		my($trans_translation)=$split_line[4];
		my($trans_biotype)=$split_line[5];
		my($no_codons)=$split_line[6];
		my($ccds_id)=$split_line[7];
		my($tsl)=$split_line[8];
		my($prot_len)=$split_line[9];

		my($firestar_annot)=$split_line[10];
		my($matador3d_annot)=$split_line[11];
		my($corsair_annot)=$split_line[12];
		my($spade_annot)=$split_line[13];
		my($thump_annot)=$split_line[14];
		my($crash_annot)=$split_line[15];
		my(@aux_crash_annot)=split(',', $crash_annot);
		my($crash_sp_annot)=$aux_crash_annot[0];
		my($crash_tp_annot)=$aux_crash_annot[1];
		my($inertia_annot)=$split_line[16];
		my($proteo_annot)=$split_line[17];
		my($appris_annot)=$split_line[18];
		my($appris_label)=$split_line[19];

		if(	defined $gene_id and defined $gene_name and defined $transcript_id and defined $trans_translation and 
			defined $trans_biotype and defined $no_codons and defined $ccds_id and defined $prot_len
		){
			#if ( exists $data_cont->{$gene_id} and exists $data_cont->{$gene_id}->{'external_id'} and defined $data_cont->{$gene_id}->{'external_id'} ) {
			#	$report->{$gene_id}->{'gene_name'}=$data_cont->{$gene_id}->{'external_id'};
			#}
			$report->{$gene_id}->{'gene_name'}=$gene_name;
			$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'biotype'}=$trans_biotype;
			unless (exists $report->{$gene_id}->{'num_trans'}) {
				$report->{$gene_id}->{'num_trans'}=0;
				$report->{$gene_id}->{'num_isof'}=0;
			}
			$report->{$gene_id}->{'num_trans'}++;
			if($no_codons ne '-') {
				$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'no_codons'}=$no_codons;
			}
			
			# add the list of proteins
			my ($seq) = undef;
			if ( exists $seq_report->{$gene_id} and exists $seq_report->{$gene_id}->{'transcripts'}->{$transcript_id} ) {
				$seq = $seq_report->{$gene_id}->{'transcripts'}->{$transcript_id};
				unless ( exists $report->{$gene_id}->{'isof_seqs'}->{$seq} ) {
					push(@{$report->{$gene_id}->{'isof_seqs'}->{$seq}}, $transcript_id); $report->{$gene_id}->{'num_isof'}++;
				}
			}

			# add the protein length
			$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'aa_length'}=$prot_len;

			# get ccds taking into account the pep sequence. We take into account the unique ccds,
			# then we add new ccds if the seq is unique.				
			if($ccds_id ne '-') {
				$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'ccds_id'}=$ccds_id;
				if ( defined $seq ) {
					if (exists $report->{$gene_id}->{'ccds_id'} and
						exists $report->{$gene_id}->{'ccds_id'}->{$ccds_id} and 
						defined $report->{$gene_id}->{'ccds_id'}->{$ccds_id}
					){
						my ($old_seq) = $report->{$gene_id}->{'ccds_id'}->{$ccds_id};
						if ($old_seq ne $seq) { # add new ccds if the seq is different
							my $key = $ccds_id.'-'.$transcript_id;
							$report->{$gene_id}->{'ccds_id'}->{$key} = $seq;
						}
					}
					else {
						$report->{$gene_id}->{'ccds_id'}->{$ccds_id} = $seq;
					}
				}
			}
			$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'tsl'}=$tsl;
			
			if(	defined $firestar_annot and defined $matador3d_annot and defined $corsair_annot and defined $spade_annot and 
				defined $thump_annot and defined $crash_sp_annot and defined $crash_sp_annot and
				defined $inertia_annot and 
				defined $proteo_annot and				 
				defined $appris_annot and defined $appris_label
			){						
				push(@{$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'annotations'}},
					$firestar_annot,
					$matador3d_annot,
					$corsair_annot,
					$spade_annot,
					$thump_annot,
					$crash_sp_annot,
					$crash_tp_annot,
					$inertia_annot,
					$proteo_annot,
					$appris_annot,
					$appris_label,
				);
			}			
		}		
	}
	return ($report, $seq_report);
	
} # end get_main_report
						
=head2 get_label_report

  Arg[1]      : String $label_file - File path of label annotations
  Arg[2]      : String $seq_file  - File path of translation sequences
  Arg[2]      : String $data_file - File path of data (Optional)
  Description : Retrieves report of label annotations.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub get_label_report($$;$)
{
	my ($label_file, $seq_file, $data_file) = @_;
	my ($report);
	
	# get content of appris annots
	my ($content) = getTotalStringFromFile($label_file);
	
	# get report from gencode
	my ($data_cont);
	if ( defined $data_file ) {
		$data_cont = APPRIS::Parser::_parse_indata($data_file);
	}
	
	# create sequence report
	my ($seq_report) = get_seq_report($seq_file);
			
	#gene_id	gene_name	transcript_id	status	biotype	no_codons ccds_id aa_len
	#
	# firestar_annot
	# matador3d_annot
	# corsair_annot
	# spade_annot
	# thump_annot
	# crash_sp_annot
	# crash_tp_annot
	# inertia_annot
	# proteo_annot
	# appris_annot
	# appris_relia
	foreach my $input_line (@{$content})
	{
		$input_line=~s/\n*$//;
		my(@split_line)=split("\t", $input_line);
		next if(scalar(@split_line)<=7);
		
		my($gene_id)=$split_line[0];
		my($gene_name)=$split_line[1];
		my($transcript_id)=$split_line[2];
		my($translation_id)=$split_line[3];
		my($trans_translation)=$split_line[4];
		my($trans_biotype)=$split_line[5];
		my($no_codons)=$split_line[6];
		my($ccds_id)=$split_line[7];
		my($tsl)=$split_line[8];
		my($prot_len)=$split_line[9];

		my($firestar_annot)=$split_line[10];
		my($matador3d_annot)=$split_line[11];
		my($corsair_annot)=$split_line[12];
		my($spade_annot)=$split_line[13];
		my($thump_annot)=$split_line[14];
		my($crash_annot)=$split_line[15];
		my(@aux_crash_annot)=split(',', $crash_annot);
		my($crash_sp_annot)=$aux_crash_annot[0];
		my($crash_tp_annot)=$aux_crash_annot[1];
		my($inertia_annot)=$split_line[16];		
		my($proteo_annot)=$split_line[17];
		my($appris_annot)=$split_line[18];
		my($appris_relia)=$split_line[19];

		if(	defined $gene_id and defined $gene_name and defined $transcript_id and defined $trans_translation and 
			defined $trans_biotype and defined $no_codons and defined $ccds_id and defined $prot_len
		){
			#if ( exists $data_cont->{$gene_id} and exists $data_cont->{$gene_id}->{'external_id'} and defined $data_cont->{$gene_id}->{'external_id'} ) {
			#	$report->{$gene_id}->{'gene_name'}=$data_cont->{$gene_id}->{'external_id'};
			#}
			$report->{$gene_id}->{'gene_name'}=$gene_name;
			unless (exists $report->{$gene_id}->{'num_trans'}) {
				$report->{$gene_id}->{'num_trans'}=0;
				$report->{$gene_id}->{'num_isof'}=0;
			}
			$report->{$gene_id}->{'num_trans'}++;
			if($no_codons ne '-') {
				$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'no_codons'}=$no_codons;
			}
			
			# add the list of proteins
			my ($seq) = undef;
			if ( exists $seq_report->{$gene_id} and exists $seq_report->{$gene_id}->{'transcripts'}->{$transcript_id} ) {
				$seq = $seq_report->{$gene_id}->{'transcripts'}->{$transcript_id};
				unless ( exists $report->{$gene_id}->{'isof_seqs'}->{$seq} ) {
					push(@{$report->{$gene_id}->{'isof_seqs'}->{$seq}}, $transcript_id); $report->{$gene_id}->{'num_isof'}++;
				}
			}
			
			# add the protein length
			$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'aa_length'}=$prot_len;			

			# get ccds taking into account the pep sequence. We take into account the unique ccds,
			# then we add new ccds if the seq is unique.				
			if($ccds_id ne '-') {
				$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'ccds_id'}=$ccds_id;
				if ( defined $seq ) {
					if (exists $report->{$gene_id}->{'ccds_id'} and
						exists $report->{$gene_id}->{'ccds_id'}->{$ccds_id} and 
						defined $report->{$gene_id}->{'ccds_id'}->{$ccds_id}
					){
						my ($old_seq) = $report->{$gene_id}->{'ccds_id'}->{$ccds_id};
						if ($old_seq ne $seq) { # add new ccds if the seq is different
							my $key = $ccds_id.'-'.$transcript_id;
							$report->{$gene_id}->{'ccds_id'}->{$key} = $seq;
						}
					}
					else {
						$report->{$gene_id}->{'ccds_id'}->{$ccds_id} = $seq;
					}
				}
			}			
			$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'tsl'}=$tsl;
			
			if(	defined $firestar_annot and defined $matador3d_annot and defined $corsair_annot and defined $spade_annot and 
				defined $thump_annot and defined $crash_sp_annot and defined $crash_sp_annot and
				defined $inertia_annot and defined $proteo_annot and
				defined $appris_annot and defined $appris_relia 
			){						
				push(@{$report->{$gene_id}->{'transcripts'}->{$transcript_id}->{'annotations'}},
					$firestar_annot,
					$matador3d_annot,
					$corsair_annot,
					$spade_annot,
					$thump_annot,
					$crash_sp_annot,
					$crash_tp_annot,
					$inertia_annot,
					$proteo_annot,
					$appris_annot,
					$appris_relia
				);
			}			
		}		
	}
	return ($report, $seq_report);
	
} # end get_label_report

sub get_prin_report($)
{
	my ($in_file) = @_;
	my ($report);
	
	# get content of appris annots
	my ($content) = getTotalStringFromFile($in_file);
	
	#gene_name	gene_id	transcript_id	ccds_id appris_label
	foreach my $input_line (@{$content})
	{
		$input_line=~s/\n*$//;
		my(@split_line)=split("\t", $input_line);
		next if(scalar(@split_line)<=0);

		my($gene_name)=$split_line[0];
		my($gene_id)=$split_line[1];
		my($transcript_id)=$split_line[2];
		my($ccds_id)=$split_line[3];
		my($appris_label)=$split_line[4];

		if(	defined $gene_name and defined $gene_id and defined $transcript_id and defined $ccds_id and defined $appris_label ){
			unless (exists $report->{$gene_id}->{'num_trans'}) {
				$report->{$gene_id}->{'num_trans'}=0;		
			}
			$report->{$gene_id}->{'num_trans'}++;
			$report->{$gene_id}->{'gene_name'} = $gene_name;
			$report->{$gene_id}->{'varsplic'}->{$transcript_id} = $appris_label;
			push( @{$report->{$gene_id}->{'appris'}->{$appris_label}}, $transcript_id);
			unless (exists $report->{$gene_id}->{'appris_num'}->{$appris_label}) {
				$report->{$gene_id}->{'appris_num'}->{$appris_label} = 1;
			}
			else {
				$report->{$gene_id}->{'appris_num'}->{$appris_label}++;
			}
		}
	}
	return $report;
	
} # end get_prin_report

=head2 src_appris_decision

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

sub src_appris_decision($$$$$;$)
{
	my ($name, $ids, $seq_file, $main_file, $label_file, $ext_data_ccds) = @_;
	
	# Declare appris cutoffs	
	my ($cfg) = new Config::IniFiles( -file =>  $CONFIG_FILE );
	$main::APPRIS_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'cutoff');
	$main::APPRIS_METHODS			= $cfg->val( 'APPRIS_VARS', 'methods');
	$main::FIRESTAR_MINRES			= $cfg->val( 'APPRIS_VARS', 'firestar_minres');
	$main::FIRESTAR_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'firestar_cutoff');
	$main::MATADOR3D_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'matador3d_cutoff');
	$main::SPADE_CUTOFF				= $cfg->val( 'APPRIS_VARS', 'spade_cutoff');
	$main::CORSAIR_AA_LEN_CUTOFF	= $cfg->val( 'APPRIS_VARS', 'corsair_aa_cutoff');
	$main::CORSAIR_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'corsair_cutoff');
	$main::THUMP_CUTOFF				= $cfg->val( 'APPRIS_VARS', 'thump_cutoff');
	$main::TRIFID_DB				= $cfg->val( 'TRIFID_VARS', 'db');
	
	# determine the methods involved in the final decision
	# By default, we use Matador3D2 but If we have genome coordinates, then we use Matador3D.
	# In this case, we don't have genome coordinates.
	my ($involved_methods) = $main::APPRIS_METHODS;
		
	# create temporal sequence file for the given ids
	my ($seq_content) = '';
	my ($seq_tmpfile) = File::Temp->new( UNLINK => 1, SUFFIX => '.appris.seq' );
	my ($seq_tmpfilename) = $seq_tmpfile->filename;
	my ($seq_found);
	my (%exist_ids) = map { $_ => 1 } @{$ids};
	my ($in) = Bio::SeqIO->new(
						-file => $seq_file,
						-format => 'Fasta'
	);
	while ( my $seq = $in->next_seq() ) {
		if ( $seq->id=~/([^|]*)/ ) {
			my ($transc_id, $transl_id, $gene_id, $gene_name, $ccds_id, $len, $s) = (undef,undef,undef,undef,undef,undef);						

			my (@ids) = split('\|', $seq->id);
			if ( scalar(@ids) > 4 ) {
				$transc_id = $ids[0];
				$transl_id = $ids[1];
				$gene_id = $ids[2];
				$gene_name = $ids[3];
				$ccds_id = $ids[4];
				$len = $ids[5];
				$s = $seq->seq;
				
				# if not exists CCDS, add from external
				if ( $ccds_id eq '-' ) {
					if ( defined $ext_data_ccds and exists $ext_data_ccds->{$s} and defined $ext_data_ccds->{$s} ) {
						$ccds_id = $ext_data_ccds->{$s};
					}
				}
				
				if ( exists $exist_ids{$gene_id} ) {
					$seq_content .= ">".$transc_id."|".$transl_id."|".$name."|".$gene_name."|".$ccds_id."|".$len."\n".$s."\n";
					$seq_found->{$gene_id} = 1;
				}				
			}
		}
	}
	if ( scalar(keys(%{$seq_found})) != scalar(keys(%exist_ids)) ) {
		die("NOT_FOUND\n".Dumper(%exist_ids)."\n");
	}
	if ( $seq_content ne '' ) {
		my ($p) = APPRIS::Utils::File::printStringIntoFile($seq_content, $seq_tmpfilename);		
	}
		
	# get sequence data
	my ($gene);
	if ( -e $seq_tmpfilename ) {
		my ($gencode_data) = APPRIS::Parser::parse_transl_data($seq_tmpfilename);
		if ( defined $gencode_data and UNIVERSAL::isa($gencode_data, 'ARRAY') and (scalar(@{$gencode_data}) > 0) ) {
				$gene = $gencode_data->[0];
		}		
	}
	
	# get object of reports from the given ids
	my ($ids_str) = join("\\|^", @{$ids});
	$ids_str = '^'.$ids_str;
	my ($main_result);
	eval {
		my ($cmd) = "grep \"$ids_str\" $main_file";
		#info("** script: $cmd\n");
		my (@cmd_out) = `$cmd`;
		$main_result = join('', @cmd_out);
	};
	my ($label_result);
	eval {
		my ($cmd) = "grep \"$ids_str\" $label_file";
		#info("** script: $cmd\n");
		my (@cmd_out) = `$cmd`;
		$label_result = join('', @cmd_out);
	};
	my ($reports) = parse_appris_methods($gene, undef, undef, undef, undef, undef, undef, undef, undef, undef, $main_result, $label_result);
		
	# get scores of methods for each transcript
	my ($scores,$s_scores) = appris::get_method_scores_from_appris_rst($gene, $reports, $involved_methods);
	
	# get list of metrics from list of involved methods
	my $involved_metrics = appris::get_method_metric_names($involved_methods);

	# get annots of methods for each transcript
	my ($annots) = appris::get_method_annots($gene, $s_scores, $involved_metrics);

	# get scores/annots of appris for each transcript
	my ($nscores) = appris::get_normalized_method_scores($gene, $annots, $involved_metrics, $scores, $s_scores);

	# get trifid report, if available
	my ($trifid_report);
	my ($gene_id) = $gene->stable_id;
	my ($cmd) = "tabix -h $main::TRIFID_DB $gene_id";
	my (@trifid_lines) = `$cmd`;
	if ( scalar(@trifid_lines) > 1 ) {
		$trifid_report = appris::parse_trifid_rst(join("\n", @trifid_lines));
	}

	# get annotations indexing each transcript
	appris::get_final_annotations($gene, $scores, $s_scores, $nscores, $annots,
	                              $involved_metrics, $trifid_report);
	
	# print outputs
	my ($score_content) = appris::get_score_output($gene, $scores, $annots);
	my ($nscore_content) = appris::get_nscore_output($gene, $nscores);
	my ($label_content) = appris::get_label_output($gene, $annots);

	return ($seq_content, $score_content, $nscore_content, $label_content);
}

1;