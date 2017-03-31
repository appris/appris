#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile prepare_workspace );

####################
# Global variables #
####################
use vars qw(
	$LOCAL_PWD
	$PROG_IN_SUFFIX
	$PROG_OUT_SUFFIX
	$WSPACE_TMP
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
my ($cfg)			= new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$WSPACE_TMP			= $ENV{APPRIS_TMP_DIR};
$WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
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
sub run_sp($$);
sub run_tp($$);
sub get_sp_cutoffs($);
sub get_tp_cutoffs($);
sub get_sp_annotations($);
sub get_tp_annotations($$);
sub print_annotations($$$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# For every sequence run the method ---------------
    my ($sp_cutoffs);
    my ($tp_cutoffs);
    my ($fasta_object) = Bio::SeqIO->new(
                        -file => $input_file,
                        -format => 'Fasta'
    );
	while ( my $seq = $fasta_object->next_seq() )
	{
		if ( $seq->id=~/^([^|]*)\|([^|]*)/ )
		{			
			my ($sequence_id) = $2;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
            my ($sequence) = $seq->seq;
            $logger->info("-- $sequence_id\n");
			
			# Create cache obj
			my ($cache) = APPRIS::Utils::CacheMD5->new(
				-dat => $sequence,
				-ws  => $WSPACE_CACHE			
			);		
			my ($seq_idx) = $cache->idx;
			my ($seq_sidx) = $cache->sidx;
			my ($seq_dir) = $cache->dir;
					
			# prepare cache dir
			my ($ws_cache) = $cache->c_dir();
			# prepare tmp dir
			my ($ws_tmp) = $WSPACE_TMP.'/'.$seq_idx;
			prepare_workspace($ws_tmp);
			
			# Cached fasta
			my ($seq_file) = $ws_cache.'/seq.faa';
			unless(-e $seq_file and (-s $seq_file > 0) ) {				
				my ($seq_cont) = ">Query\n$seq\n";			
				open (SEQ_FILE,">$seq_file");
				print SEQ_FILE $seq_cont;
				close (SEQ_FILE);
			}
						
			# Run SignalP ----------
			$logger->info("##Run signalp ---------------\n");
			my ($sp_result) = run_sp($seq_file, $ws_cache);
			
			# Run TargetP ----------
			$logger->info("##Run targetp ---------------\n");
			my ($tp_result) = run_tp($seq_file, $ws_cache);
			
			# Parse the results ----------
			$logger->info("##Get sp cutoffs ---------------\n");
			$sp_cutoffs->{$sequence_id} = get_sp_cutoffs($sp_result);
			$logger->debug("\n".Dumper($sp_cutoffs)."\n");
			$logger->info("##Get tp cutoffs ---------------\n");
			$tp_cutoffs->{$sequence_id} = get_tp_cutoffs($tp_result);
			$logger->debug("\n".Dumper($tp_cutoffs)."\n");		
		}
	}	

	# Get annotations ----------
	$logger->info("##Get sp annots ---------------\n");
	my ($sp_annotations) = get_sp_annotations($sp_cutoffs);
	$logger->debug("\n".Dumper($sp_annotations)."\n");
	$logger->info("##Get tp annots ---------------\n");
	my ($tp_annotations) = get_tp_annotations($tp_cutoffs, $sp_cutoffs);
	$logger->debug("\n".Dumper($tp_annotations)."\n");
	
	# Print annotations ----------
	$logger->info("##Print annots ---------------\n");
	my ($crash_annotations) = print_annotations($sp_cutoffs, $sp_annotations, $tp_cutoffs, $tp_annotations);
	my ($print_out) = printStringIntoFile($crash_annotations, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}	
	
	$logger->finish_log();
	
	exit 0;
}

# Run SignalP
sub run_sp($$)
{
	my ($input_file, $ws_cache) = @_;
	
	# Declare variables
	my ($output_file) = $ws_cache.'/seq.signalp';	
	my ($err_file) = '/dev/null';
	
	# Execute program
	unless ( -e $output_file and (-s $output_file > 0 ) ) {
		eval {
			my ($cmd) = "$RUN_PROGRAM_1 $input_file 1> $output_file 2> $err_file";
			$logger->debug("\n** script: $cmd\n");		
			system($cmd) == 0 or $logger->error("ERROR: running signalp program");
		};
		$logger->error("ERROR: running signalp program") if($@);		
	}
	
	# Get file result
	local (*FILE);
	open (FILE, $output_file) or die ("ERROR: openning signalp result");
	my (@sp_result) = <FILE>;
	close(FILE);
	
	return \@sp_result;
}

# Run TargetP
sub run_tp($$)
{
	my ($input_file, $ws_cache) = @_;

	# Declare variables	
	my ($output_file) = $ws_cache.'/seq.targetp';	
	my ($err_file) = '/dev/null';

	# Execute program
	unless ( -e $output_file and (-s $output_file > 0 ) ) {
		eval {
			my ($cmd) = "$RUN_PROGRAM_2 $input_file 1> $output_file 2> $err_file";
			$logger->debug("\n** script: $cmd\n");		
			system($cmd) == 0 or $logger->error("ERROR: running targetp program");
		};
		$logger->error("ERROR: running targetp program") if($@);
	}
		
	# Get file result
	local (*FILE);
	open (FILE, $output_file) or die ("ERROR: openning targetp result");
	my (@tp_result) = <FILE>;
	close(FILE);
	
	return \@tp_result;
}

# Get the cutoffs from SignalP
sub get_sp_cutoffs($)
{
	my ($result) = @_;
	my ($cutoffs);
	foreach my $line (@{$result})
	{		
		# SignalP-NN result: # Measure  Position  Value  Cutoff  signal peptide?
		if ( $line =~ /^\s*mean S\s*([^-]+)\-([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\n]+)\n*/mg ) {
			my ($measure) = 'Smean';
			$cutoffs->{$measure}->{'start'}=$1 if (defined $1);
			$cutoffs->{$measure}->{'end'}=$2  if (defined $2);
			$cutoffs->{$measure}->{'score'}=$3 if (defined $3);
			$cutoffs->{$measure}->{'cutoff'}=$4 if (defined $4);
			$cutoffs->{$measure}->{'peptide'}=$5 if (defined $5);				

			# Init the range of peptide signal
			if($1 > 0) {
				$cutoffs->{'start'}=$1;							
			} else {
				$cutoffs->{'start'}=1;
			}
			if($2 > 0) {
				$cutoffs->{'end'}=$2;							
			} else {
				$cutoffs->{'end'}=1;
			}	
		}
		if ( $line =~ /^\s*D\s*([^-]+)\-([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\n]+)\n*/mg ) {
			my ($measure) = 'Dscore';
			$cutoffs->{$measure}->{'start'}=$1 if (defined $1);
			$cutoffs->{$measure}->{'end'}=$2 if (defined $2);
			$cutoffs->{$measure}->{'score'}=$3 if (defined $3);
			$cutoffs->{$measure}->{'cutoff'}=$4 if (defined $4);
			$cutoffs->{$measure}->{'peptide'}=$5 if (defined $5);
			
			# Get the biggest range of peptide signal
			if (exists $cutoffs->{'start'} and defined $cutoffs->{'start'} and
				($cutoffs->{'start'} > 0) and ($1 < $cutoffs->{'start'})) {
					$cutoffs->{'start'}=$1;
			}
			if (exists $cutoffs->{'end'} and defined $cutoffs->{'end'} and
				($cutoffs->{'end'} > 0) and ($2 > $cutoffs->{'end'})) {
					$cutoffs->{'end'}=$2;
			}				
		}
		
		# SignalP-HMM result:
		if ( $line =~ /^\s*Signal peptide probability:\s*([^\n]+)\n*/mg ) {
			my ($measure) = 'Sprob';
			$cutoffs->{$measure}->{'score'}=$1 if(defined($1));
		}
		if( $line =~ /^\s*Max cleavage site probability:\s*([^\s]+)\s*between\s*pos\.\s*([^\s]+)\s*and\s*([^\n]+)\n*/mg ) {
			my ($measure) = 'Cmax';
			$cutoffs->{$measure}->{'score'}=$1 if (defined $1);
			$cutoffs->{$measure}->{'start'}=$2 if (defined $2);
			$cutoffs->{$measure}->{'end'}=$3 if (defined $3);

			# Get the biggest range of peptide section
			if (exists $cutoffs->{'start'} and defined $cutoffs->{'start'} and
				($2 > 0) and ($2 < $cutoffs->{'start'})) {
				$cutoffs->{'start'}=$2;
			}
			if (exists $cutoffs->{'end'} and defined $cutoffs->{'end'} and
				($3 > 0) and ($3 > $cutoffs->{'end'})) {
					$cutoffs->{'end'}=$3;
			}
		}
	}
	
	return $cutoffs;
}

# Get the cutoffs from TargetP
sub get_tp_cutoffs($)
{
	my ($result) = @_;
	my ($cutoffs);
	my ($start);
	foreach my $line (@{$result})
	{		
		# Say when start and end the result
		if (defined $start and $line =~ /^-/) {
			$start = undef;
		} elsif ($line =~ /^-/) {
			$start = 1;
		} 
		#Name                  Len            mTP     SP  other  Loc  RC
		elsif (defined $start and ($line =~ /^([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*/mg)) {
			$cutoffs->{'localization'}->{'score'} = $6 if (defined $6);
			$cutoffs->{'reliability'}->{'score'} = $7 if (defined $7);
		}
	}	
	return $cutoffs;	
}

# Get the annotations of SignalP
sub get_sp_annotations($)
{
	my ($cutoffs) = @_;

	# Declare variables
	my ($annots);

	# Scan transcripts
	while ( my ($id,$features) = each(%{$cutoffs}) ) {
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
	return $annots;	
}

# Get the annotations of TargetP
sub get_tp_annotations($$)
{
	my ($cutoffs, $sp_cuttoffs) = @_;

	# Declare variables
	my ($annots);

	# Scan transcripts
	while ( my ($id,$features) = each(%{$cutoffs}) ) {
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
	return $annots;	
}

# Print annotations of CRASH
sub print_annotations($$$$)
{
	my ($sp_cutoffs, $sp_annotations, $tp_cutoffs, $tp_annotations) = @_;
	my ($annots) = '';
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
	foreach my $id (sort {$a cmp $b} keys(%{$sp_cutoffs})) {
		if ( exists $sp_cutoffs->{$id} and defined $sp_cutoffs->{$id} ) {
			my ($sp_features) = $sp_cutoffs->{$id};
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