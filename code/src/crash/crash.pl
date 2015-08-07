#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile );

####################
# Global variables #
####################
use vars qw(
	$LOCAL_PWD
	$PROG_IN_SUFFIX
	$PROG_OUT_SUFFIX
	$WSPACE_BASE
	$WSPACE_CACHE
	$RUN_PROGRAM_1
	$RUN_PROGRAM_2
	$PROG_IN_SUFFIX
	$PROG1_OUT_SUFFIX
	$PROG2_OUT_SUFFIX
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL
);

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($input_file) = undef;
my ($output_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'conf=s'			=> \$config_file,
	'input=s'			=> \$input_file,
	'output=s'			=> \$output_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);
unless ( defined $config_file and defined $input_file and defined $output_file )
{
	print `perldoc $0`;
	exit 1;
}

# Get conf vars
my ($cfg) = new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$WSPACE_BASE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('CRASH_VARS', 'name').'/';
$WSPACE_CACHE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('CACHE_VARS', 'name').'/';
$RUN_PROGRAM_1		= $cfg->val( 'CRASH_VARS', 'program1');
$RUN_PROGRAM_2		= $cfg->val( 'CRASH_VARS', 'program2');
$PROG_IN_SUFFIX		= 'faa';
$PROG1_OUT_SUFFIX	= 'signalp';
$PROG2_OUT_SUFFIX	= 'targetp';
$OK_LABEL			= 'YES';
$UNKNOWN_LABEL		= 'UNKNOWN';
$NO_LABEL			= 'NO';

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);

#####################
# Method prototypes #
#####################
sub get_tmp_file($$);
sub run_sp($$$);
sub run_tp($$$);
sub get_sp_cutoffs($);
sub get_tp_cutoffs($);
sub get_sp_annotations($$);
sub get_tp_annotations($$$);
sub print_annotations($$$$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Get id of output file ----------
	$logger->info("##Get id ---------------\n");
	my ($id);
	my (@out_file_paths) = split('/',$output_file);
	if ( scalar(@out_file_paths) > 0 ) {
		$id = $out_file_paths[scalar(@out_file_paths)-1];
		$id =~ s/\.[^\$]*$//g;
	}
	
	# Get list transcripts ----------
	$logger->info("##Get list of transcripts ---------------\n");
	my ($trans_list, $tmp_in_file) = get_tmp_file($id, $input_file);
	$logger->debug("\n".Dumper($trans_list)."\n");
            
	# Run SignalP ----------
	$logger->info("##Run signalp ---------------\n");
	my ($sp_result) = run_sp($id, $tmp_in_file, $output_file);
	
	# Run TargetP ----------
	$logger->info("##Run targetp ---------------\n");
	my ($tp_result) = run_tp($id, $tmp_in_file, $output_file);
	
	# Parse the results ----------
	$logger->info("##Get sp cutoffs ---------------\n");
	my ($sp_cutoffs) = get_sp_cutoffs($sp_result);
	$logger->debug("\n".Dumper($sp_cutoffs)."\n");
	$logger->info("##Get tp cutoffs ---------------\n");
	my ($tp_cutoffs) = get_tp_cutoffs($tp_result);
	$logger->debug("\n".Dumper($tp_cutoffs)."\n");

	# Get annotations ----------
	$logger->info("##Get sp annots ---------------\n");
	my ($sp_annotations) = get_sp_annotations($trans_list, $sp_cutoffs);
	$logger->debug("\n".Dumper($sp_annotations)."\n");
	$logger->info("##Get tp annots ---------------\n");
	my ($tp_annotations) = get_tp_annotations($trans_list, $tp_cutoffs, $sp_cutoffs);
	$logger->debug("\n".Dumper($tp_annotations)."\n");
	
	# Print annotations ----------
	$logger->info("##Print annots ---------------\n");
	my ($crash_annotations) = print_annotations($trans_list, $sp_cutoffs, $sp_annotations, $tp_cutoffs, $tp_annotations);
	my ($print_out) = printStringIntoFile($crash_annotations, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}
	
	$logger->finish_log();
	
	exit 0;
}

# Get list of trans
sub get_tmp_file($$)
{
	my ($id,$input_file) = @_;
	my ($list);
	my ($tmp_in_file);
	my ($tmp_in_cont) = '';
    my $fasta_object = Bio::SeqIO->new(
                        -file => $input_file,
                        -format => 'Fasta'
    );
	while ( my $seq = $fasta_object->next_seq() )
	{
		if ( $seq->id=~/([^|]*)/ )
		{			
			my ($sequence_id) = $1;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
			my ($sequence) = $seq->seq;
			if ( length($sequence) > 2 ) {
				$tmp_in_cont .= ">$sequence_id\n$sequence\n";
			}
            push(@{$list}, $sequence_id);
		}
	}
	if ( $tmp_in_cont ne '' ) {
		$tmp_in_file = $WSPACE_CACHE.'/'.$id.'.'.$PROG_IN_SUFFIX;
		my ($print_out) = printStringIntoFile($tmp_in_cont, $tmp_in_file);
		unless( defined $print_out ) {
			$logger->error("Can not create output file: $!\n");
		}		
	}
		
	return ($list,$tmp_in_file);
}

# Run SignalP
sub run_sp($$$)
{
	my ($id, $input_file, $output_file) = @_;
	
	# Declare variables
	my ($sp_input_file) = $input_file;
	my ($sp_output_file) = $WSPACE_CACHE.'/'.$id.'.'.$PROG1_OUT_SUFFIX;
	my ($sp_err_file) = '/dev/null';
	
	# Execute program
	unless ( -e $sp_output_file and (-s $sp_output_file > 0 ) ) {
		eval {
			$logger->debug("$RUN_PROGRAM_1 $sp_input_file 1> $sp_output_file 2> $sp_err_file\n");		
			system("$RUN_PROGRAM_1 $sp_input_file 1> $sp_output_file 2> $sp_err_file") == 0 or $logger->error("ERROR: running signalp program");
		};
		$logger->error("ERROR: running signalp program") if($@);		
	}
	
	# Get file result
	local (*FILE);
	open (FILE, $sp_output_file) or die ("ERROR: openning signalp result");
	my (@sp_result) = <FILE>;
	close(FILE);
	
	return \@sp_result;
}

# Run TargetP
sub run_tp($$$)
{
	my ($id, $input_file, $output_file) = @_;

	# Declare variables	
	my ($tp_input_file) = $input_file;
	my ($tp_output_file) = $WSPACE_CACHE.'/'.$id.'.'.$PROG2_OUT_SUFFIX;	
	my ($tp_err_file) = '/dev/null';

	# Execute program
	unless ( -e $tp_output_file and (-s $tp_output_file > 0 ) ) {
		eval {
			$logger->debug("$RUN_PROGRAM_2 $tp_input_file 1> $tp_output_file 2> $tp_err_file\n");		
			system("$RUN_PROGRAM_2 $tp_input_file 1> $tp_output_file 2> $tp_err_file") == 0 or $logger->error("ERROR: running targetp program");
		};
		$logger->error("ERROR: running targetp program") if($@);
	}
		
	# Get file result
	local (*FILE);
	open (FILE, $tp_output_file) or die ("ERROR: openning targetp result");
	my (@tp_result) = <FILE>;
	close(FILE);
	
	return \@tp_result;
}

# Get the cutoffs from SignalP
sub get_sp_cutoffs($)
{
	my ($result) = @_;

	# Declare variables
	my ($cutoffs);
	my ($id);
	
	foreach my $line (@{$result}) {
		
		# Get the identifier of protein
		if ($line =~ /^>([^\s]*)\s*/) {
			$id = $1;
		}
		
		# SignalP-NN result: # Measure  Position  Value  Cutoff  signal peptide?
		if (defined $id and ($line =~ /^\s*mean S\s*([^-]+)\-([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\n]+)\n*/mg)) {
			my ($measure) = 'Smean';
			$cutoffs->{$id}->{$measure}->{'start'}=$1 if (defined $1);
			$cutoffs->{$id}->{$measure}->{'end'}=$2  if (defined $2);
			$cutoffs->{$id}->{$measure}->{'score'}=$3 if (defined $3);
			$cutoffs->{$id}->{$measure}->{'cutoff'}=$4 if (defined $4);
			$cutoffs->{$id}->{$measure}->{'peptide'}=$5 if (defined $5);				

			# Init the range of peptide signal
			if($1 > 0) {
				$cutoffs->{$id}->{'start'}=$1;							
			} else {
				$cutoffs->{$id}->{'start'}=1;
			}
			if($2 > 0) {
				$cutoffs->{$id}->{'end'}=$2;							
			} else {
				$cutoffs->{$id}->{'end'}=1;
			}	
		}
		if (defined $id and ($line =~ /^\s*D\s*([^-]+)\-([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\n]+)\n*/mg)) {
			my ($measure) = 'Dscore';
			$cutoffs->{$id}->{$measure}->{'start'}=$1 if (defined $1);
			$cutoffs->{$id}->{$measure}->{'end'}=$2 if (defined $2);
			$cutoffs->{$id}->{$measure}->{'score'}=$3 if (defined $3);
			$cutoffs->{$id}->{$measure}->{'cutoff'}=$4 if (defined $4);
			$cutoffs->{$id}->{$measure}->{'peptide'}=$5 if (defined $5);
			
			# Get the biggest range of peptide signal
			if (exists $cutoffs->{$id}->{'start'} and defined $cutoffs->{$id}->{'start'} and
				($cutoffs->{$id}->{'start'} > 0) and ($1 < $cutoffs->{$id}->{'start'})) {
					$cutoffs->{$id}->{'start'}=$1;
			}
			if (exists $cutoffs->{$id}->{'end'} and defined $cutoffs->{$id}->{'end'} and
				($cutoffs->{$id}->{'end'} > 0) and ($2 > $cutoffs->{$id}->{'end'})) {
					$cutoffs->{$id}->{'end'}=$2;
			}				
		}
		
		# SignalP-HMM result:
		if (defined $id and ($line =~ /^\s*Signal peptide probability:\s*([^\n]+)\n*/mg)) {
			my ($measure) = 'Sprob';
			$cutoffs->{$id}->{$measure}->{'score'}=$1 if(defined($1));
		}
		if(defined $id and ($line =~ /^\s*Max cleavage site probability:\s*([^\s]+)\s*between\s*pos\.\s*([^\s]+)\s*and\s*([^\n]+)\n*/mg)) {
			my ($measure) = 'Cmax';
			$cutoffs->{$id}->{$measure}->{'score'}=$1 if (defined $1);
			$cutoffs->{$id}->{$measure}->{'start'}=$2 if (defined $2);
			$cutoffs->{$id}->{$measure}->{'end'}=$3 if (defined $3);

			# Get the biggest range of peptide section
			if (exists $cutoffs->{$id}->{'start'} and defined $cutoffs->{$id}->{'start'} and
				($2 > 0) and ($2 < $cutoffs->{$id}->{'start'})) {
				$cutoffs->{$id}->{'start'}=$2;
			}
			if (exists $cutoffs->{$id}->{'end'} and defined $cutoffs->{$id}->{'end'} and
				($3 > 0) and ($3 > $cutoffs->{$id}->{'end'})) {
					$cutoffs->{$id}->{'end'}=$3;
			}
		}
	}
	
	return $cutoffs;
}

# Get the cutoffs from TargetP
sub get_tp_cutoffs($)
{
	my ($result) = @_;

	# Declare variables
	my ($cutoffs);
	my ($id);	
	my ($start);
	
	foreach my $line (@{$result}) {
		
		# Say when start and end the result
		if (defined $start and $line =~ /^-/) {
			$start = undef;
		} elsif ($line =~ /^-/) {
			$start = 1;
		} 
		#Name                  Len            mTP     SP  other  Loc  RC
		elsif (defined $start and ($line =~ /^([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*/mg)) {
			my ($id) = $1;
			$cutoffs->{$id}->{'localization'}->{'score'} = $6 if (defined $6);
			$cutoffs->{$id}->{'reliability'}->{'score'} = $7 if (defined $7);
		}
	}
	
	return $cutoffs;	
}

# Get the annotations of SignalP
sub get_sp_annotations($$)
{
	my ($list,$cutoffs) = @_;

	# Declare variables
	my ($annots);

	# Scan transcripts
	foreach my $id (@{$list}) {
		if ( exists $cutoffs->{$id} and defined $cutoffs->{$id} ) {
			my ($features) = $cutoffs->{$id};

			if(	exists $features->{'Smean'} and defined $features->{'Smean'} and
				exists $features->{'Dscore'} and defined $features->{'Dscore'} and
				exists $features->{'Sprob'} and defined $features->{'Sprob'} and
				exists $features->{'Cmax'} and defined $features->{'Cmax'}
			)
			{
				my ($score_num) = 0;
				my ($score) = 0;
				
				if (exists $features->{'Smean'}->{'score'} and defined $features->{'Smean'}->{'score'}) {
					if ($features->{'Smean'}->{'score'} > 0.85) {
						$score += 1;
					} elsif ($features->{'Smean'}->{'score'} < 0.68) {
						$score += -1;
					} else {
						$score += 0;
					}
					$score_num++;
				}
				
				if (exists $features->{'Dscore'}->{'score'} and defined $features->{'Dscore'}->{'score'}) {
					if ($features->{'Dscore'}->{'score'} > 0.75) {
						$score += 1;
					} elsif ($features->{'Dscore'}->{'score'} < 0.53) {
						$score += -1;
					} else {
						$score += 0;
					}
					$score_num++;
				}
				if (exists $features->{'Sprob'}->{'score'} and defined $features->{'Sprob'}->{'score'}) {
					if ($features->{'Sprob'}->{'score'} > 0.99) {
						$score += 1;
					} elsif ($features->{'Sprob'}->{'score'} < 0.8) {
						$score += -1;
					} else {
						$score += 0;
					}
					$score_num++;
				}
				if (exists $features->{'Cmax'}->{'score'} and defined $features->{'Cmax'}->{'score'}) {
					if ($features->{'Cmax'}->{'score'} > 0.85) {
						$score += 1;
					} elsif ($features->{'Cmax'}->{'score'} < 0.40) {
						$score += -1;
					} else {
						$score += 0;
					}
					$score_num++;
				}
				# Get the annotation from every score
				if ($score_num == 4) {
					$annots->{$id}->{'score'} = $score;
					my ($label);
					if ($score >= 2) {
						$label =$OK_LABEL;
					} elsif (($score == 0) or ($score == 1)) {
						$label = $UNKNOWN_LABEL;
					} elsif ($score <= -1) {
						$label = $NO_LABEL;
					}
					if (defined $label) {
						$annots->{$id}->{'peptide_signal'} = $label;
					}
				}
			}
		}
		else {
			$annots->{$id}->{'peptide_signal'} = $NO_LABEL;
		}
	}

	return $annots;	
}

# Get the annotations of TargetP
sub get_tp_annotations($$$)
{
	my ($list, $cutoffs, $sp_cuttoffs) = @_;

	# Declare variables
	my ($annots);

	# Scan transcripts
	foreach my $id (@{$list}) {
		if ( exists $cutoffs->{$id} and defined $cutoffs->{$id} ) {
			my ($features) = $cutoffs->{$id};
	
			if(	exists $features->{'localization'} and defined $features->{'localization'} and
				exists $features->{'reliability'} and defined $features->{'reliability'}
			)
			{
				my ($score_num) = 0;
				my ($score) = 0;
				if(exists $features->{'localization'}->{'score'} and defined $features->{'localization'}->{'score'}) {
					if ($features->{'localization'}->{'score'} eq 'M') {
						$score += 1;
					} else {
						$score += -4;
					}
					$score_num++;
				}
				if (exists $features->{'reliability'}->{'score'} and defined $features->{'reliability'}->{'score'}) {
					if ($features->{'reliability'}->{'score'} <= 2) {
						$score += 1;
					} elsif ($features->{'reliability'}->{'score'} >= 4) {
						$score += -1;
					} else {
						$score += 0;
					}
					$score_num++;
				}
				# Get signalp info
				if(exists $sp_cuttoffs->{$id}->{'Sprob'} and defined $sp_cuttoffs->{$id}->{'Sprob'}) {
					if ($sp_cuttoffs->{$id}->{'Sprob'}->{'score'} > 0.99) {
						$score += 1;
					} elsif ($sp_cuttoffs->{$id}->{'Sprob'}->{'score'} < 0.8) {
						$score += -1;
					} else {
						$score += 0;
					}
					$score_num++;
				}
				if (exists $sp_cuttoffs->{$id}->{'Cmax'} and defined $sp_cuttoffs->{$id}->{'Cmax'}) {
					if($sp_cuttoffs->{$id}->{'Cmax'}->{'score'} > 0.85) {
						$score += 1;
					} elsif ($sp_cuttoffs->{$id}->{'Cmax'}->{'score'} < 0.40) {
						$score += -1;
					} else {
						$score += 0;
					}
					$score_num++;
				}					
	
				if ($score_num == 4) {
					$annots->{$id}->{'score'} = $score;
					my ($label);
					if ($score >= 2) {
						$label = $OK_LABEL;
					} elsif (($score == 0) or ($score == 1)) {
						$label = $UNKNOWN_LABEL;
					} elsif ($score <=- 1) {
						$label = $NO_LABEL;
					}
					if (defined $label) {
						$annots->{$id}->{'mitochondrial_signal'} = $label;
					}
				}
			}
		}
		else {
			$annots->{$id}->{'mitochondrial_signal'} = $NO_LABEL;
		}
	}

	return $annots;	
}

# Print annotations of CRASH
sub print_annotations($$$$$)
{
	my ($list, $sp_cutoffs, $sp_annotations, $tp_cutoffs, $tp_annotations) = @_;
	
	# Declare variables
	my ($annots) = '';

	# Scan transcripts
	$annots .= "### crash: signalp 3.0 and targetp v1.1 prediction results ##################################";
	$annots .= "\n";
	$annots .= "----------------------------------------------------------------------\n";
	$annots .= 'id'."\t".
				'start'."\t".
				'end'."\t".
				's_mean'."\t".
				'd_score'."\t".
				'c_max'."\t".
				's_prob'."\t".
				'sp_score'."\t".
				'peptide_signal'."\t".
				'localization'."\t".
				'reliability'."\t".
				'tp_score'."\t".
				'mitochondrial_signal'."\n";
	foreach my $id (@{$list}) {
		if ( exists $sp_cutoffs->{$id} and defined $sp_cutoffs->{$id} ) {
			my ($sp_features) = $sp_cutoffs->{$id};
				
	#while ( my ($id, $sp_features) = each(%{$sp_cutoffs}) ) {
			$annots .= 	">".$id."\t".
							$sp_features->{'start'}."\t".
							$sp_features->{'end'}."\t".
							$sp_features->{'Smean'}->{'score'}."\t".
							$sp_features->{'Dscore'}->{'score'}."\t".
							$sp_features->{'Sprob'}->{'score'}."\t".
							$sp_features->{'Cmax'}->{'score'}."\t".
							$sp_annotations->{$id}->{'score'}."\t".
							$sp_annotations->{$id}->{'peptide_signal'}."\t".
							$tp_cutoffs->{$id}->{'localization'}->{'score'}."\t".
							$tp_cutoffs->{$id}->{'reliability'}->{'score'}."\t".
							$tp_annotations->{$id}->{'score'}."\t".
							$tp_annotations->{$id}->{'mitochondrial_signal'}."\n";
		}
		else {
			$annots .= 	">".$id."\t".
							'-'."\t".
							'-'."\t".
							'-'."\t".
							'-'."\t".
							'-'."\t".
							'-'."\t".
							'NULL'."\t".
							$sp_annotations->{$id}->{'peptide_signal'}."\t".
							'-'."\t".
							'-'."\t".
							'NULL'."\t".
							$tp_annotations->{$id}->{'mitochondrial_signal'}."\n";
		}
	}
	
		
	return $annots;	
}

main();

1;

__END__

=head1 NAME

crash

=head1 DESCRIPTION

Run web services of CRASH, which executes SignalP and TargetP programs.

* SignalP predicts the presence and location of signal peptide cleavage sites 
in amino acid sequences from different organisms: Gram-positive bacteria,
Gram-negative bacteria, and eukaryotes.

For more information:
  http://www.cbs.dtu.dk/services/SignalP/

* TargetP predicts the subcellular location of eukaryotic protein sequences.

For more information:
  http://www.cbs.dtu.dk/services/TargetP/

=head1 VERSION

0.1

=head1 PARAMETERS

--input    Amino acid sequences file as FASTA format

--output   Outut file where method reponse will be saving

=head2 Required arguments:

Specific parameters by service:

* SignalP:

    --format    Output format [full, summary, short] (full, by default)
  
    --truncate  Truncate each sequence to max. number of residues (70, by default)

    --type      Organism group: (euk, by default)
                     euk   => Eukaryutes
                     gram- => Gram-negative bacteria
                     gram+ => Gram-positive bacteria

    --method    Neural networks, Hidden Markov models, Both [nn, hmm, nn+hmm] (nn+hmm, by default)
                     
* TargetP: No specific parameters. Only default values:
    
      Non-plant Organism group
      No prediction scope
      No cutoffs

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
	
	
=head1 EXAMPLE

perl crash.pl
								
    --input=examples/OTTHUMG00000086781.faa
    
    --output=examples/OTTHUMG00000086781.out

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut