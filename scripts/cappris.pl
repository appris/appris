#!/usr/bin/perl -W

use strict;
use warnings;
use threads;

use Getopt::Long;
use FindBin;
use Config::IniFiles;
use Bio::SeqIO;
use List::Util qw(sum);
use Data::Dumper;

use APPRIS::Utils::File qw( updateStringIntoFile );
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################

# TEMPORAL SOLUTION!!!!
$ENV{APPRIS_METHODS}="firestar,matador3d,spade,corsair,thump,crash,appris";
# TEMPORAL SOLUTION!!!!

# Input parameters
my ($id) = undef;
my ($inpath) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'id=s'				=> \$id,
	'inpath=s'			=> \$inpath,
	'outfile=s'			=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
	defined  $id and
	defined  $inpath and
	defined  $outfile 	
){
	print `perldoc $0`;
	exit 1;
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();

#################
# Method bodies #
#################
sub _avg_seq_length($);
sub check_appris($$);

# Main subroutine
sub main()
{	
	# Create gene report
	$logger->info("-- check time results\n");
	my ($runtime) = check_appris($id, $inpath);
	unless ( defined $runtime ) {
		$logger->error("checking time results");
	}

	# Create gene report
	$logger->info("-- print time results\n");
	my ($a) = updateStringIntoFile($runtime, $outfile);
	unless ( defined $a ) {
		$logger->error("printing time results");
	}
	
	$logger->finish_log();
	
	exit 0;	
}

sub check_appris($$)
{
	my ($gene_id, $workspace) = @_;
	my ($runtime) = undef;
	my ($avg_length) = 0;
	my ($output) = '';
		
	# input files
	my ($in_files) = {
			'data'			=> $workspace.'/'.$gene_id.'.annot.gtf',
			'pdata'			=> $workspace.'/'.$gene_id.'.pannot.gtf',
			'transc'		=> $workspace.'/'.$gene_id.'.transc.fa',
			'transl'		=> $workspace.'/'.$gene_id.'.transl.fa',
			'cdsseq'		=> $workspace.'/'.$gene_id.'.cdsseq.fa',			
			'firestar'		=> $workspace.'/'.$gene_id.'.firestar',
			'matador3d'		=> $workspace.'/'.$gene_id.'.matador3d',
			'spade'			=> $workspace.'/'.$gene_id.'.spade',
			'corsair'		=> $workspace.'/'.$gene_id.'.corsair',
			'thump'			=> $workspace.'/'.$gene_id.'.thump',
			'crash'			=> $workspace.'/'.$gene_id.'.crash',
			'appris'		=> $workspace.'/'.$gene_id.'.appris',
			'appris-label'	=> $workspace.'/'.$gene_id.'.appris.label',
			'appris-score'	=> $workspace.'/'.$gene_id.'.appris.score',
			'logfile'		=> $workspace.'/'.$gene_id.'.log',
			'errfile'		=> $workspace.'/'.$gene_id.'.err',
	};
	
	# loop until find log file
	#$logger->info("-- loop until find log file");
	#while ( !(-e $in_files->{'logfile'}) ) {
	#	sleep(20);
	#	$logger->info(".");		
	#}
	#sleep(2);
	#$logger->info("\n");
	
	# check if pipeline has got inputs
	if ( -e $in_files->{'data'} and (-s $in_files->{'data'} > 0) and
		 -e $in_files->{'pdata'} and (-s $in_files->{'pdata'} > 0) and 
		 -e $in_files->{'transc'} and (-s $in_files->{'transc'} > 0) and
		 -e $in_files->{'transl'} and (-s $in_files->{'transl'} > 0) 
		 #-e $in_files->{'cdsseq'} and (-s $in_files->{'cdsseq'} > 0)
	) {
		# get the avarage length of transcripts
		my ($avg_length) = _avg_seq_length($in_files->{'transl'});
		
		# check if pipeline run correctly
		$logger->info("-- check whether appris pipeline finished correctly\n");
		eval {
			my ($cmd) = 'sed -e \'/./{H;$!d;}\' -e \'x;/All done for/!d\' '.$in_files->{'logfile'}.' ';
			$logger->info("** $cmd\n");
			my (@tail_logfile) = `$cmd`;		
			my ($logfile_cont) = join "",@tail_logfile;		
			my ($pattern2) = 'All done for ([^\.]*)\.pl\.';
			my ($pattern3) = '\d* warnings\.';
			my ($pattern4) = 'Runtime\: (\d*)h (\d*)min (\d*)sec \[([^\,]*)\, mem ([^\]]*)\]';
			while ( $logfile_cont =~ /$pattern2\n*$pattern3\s*$pattern4/mg ) {			
				# if a time already exists then save the bigger one
				if ( exists $runtime->{$1} ) {
					my ($passtime) = $runtime->{$1}->{'hour'}.$runtime->{$1}->{'min'}.$runtime->{$1}->{'sec'};
					my ($currenttime) = $2.$3.$4;
					if ( $currenttime - $passtime > 0 ) { # Warning: this is not a good way to compare date-times
						$runtime->{$1} = {
							'hour'	=> $2,
							'min'	=> $3,
							'sec'	=> $4,
							'date'	=> $5,
							'mem'	=> $6
						};					
					}
				}
				else {
					$runtime->{$1} = {
						'hour'	=> $2,
						'min'	=> $3,
						'sec'	=> $4,
						'date'	=> $5,
						'mem'	=> $6
					};				
				}
			}
		};
		throw("checking whether pipeline executed correctly: $!\n") if($@);
		
		# print the file-runtime for each method
		if ( defined $runtime ) {
			$output .= $gene_id."\t".$avg_length."\t";
			foreach my $method ( split(',',$ENV{APPRIS_METHODS}) ) {
				# result files
				if ( -e $in_files->{$method} and (-s $in_files->{$method} > 0) ) {
					$output .= "OK,";
				}
				else {
					$output .= "NO,";
				}
				
				# runtimes
				my ($h) = 0;
				my ($m) = 0;
				my ($s) = 0;
				if ( exists $runtime->{$method} and defined $runtime->{$method} ) {
					my ($rtime) = $runtime->{$method};
					$h = $rtime->{'hour'};
					$m = $rtime->{'min'};
					$s = $rtime->{'sec'};
				}
				$output .= "$h:$m:$s"."\t";
			}
			$output =~ s/\t$/\n/;		
		}
	}
	else {
		# print empty results for each method
		if ( defined $runtime ) {
			$output .= $gene_id."\t".$avg_length."\t";
			foreach my $method ( split(',',$ENV{APPRIS_METHODS}) ) {
				# result files
				$output .= "-,";
				
				# runtimes
				my ($h) = 0;
				my ($m) = 0;
				my ($s) = 0;
				$output .= "$h:$m:$s"."\t";
			}
			$output =~ s/\t$/\n/;		
		}
	}
	
	return $output;
} # end check_appris

sub _avg_seq_length($)
{
	my ($file) = @_;
	my ($data) = 0;
	if (-e $file and (-s $file > 0) ) {
		my (@lengths);
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() ) {
			if ( $seq->id=~/([^|]*)\|([^|]*)/ ) {
				push(@lengths, length($seq->seq));
			}
		}		
		if ( scalar(@lengths) > 0 ) {
			sub mean { return @_ ? sum(@_) / @_ : 0 }
			$data = int(mean(@lengths));
		}
	}
	return $data;
} # End _avg_seq_length

main();


1;

__END__

=head1 NAME

cappris

=head1 DESCRIPTION

script that checks the results of APPRIS 

=head1 SYNOPSIS

run_appris

=head2 Input arguments:

	--id= <Ensembl gene identifier>
	
	--inpath= <Acquire input files from PATH>
	
	--outfile= <Output file>
	
=head2 Optional arguments (log arguments):
	
		--loglevel=LEVEL <define log level (default: NONE)>	
	
		--logfile=FILE <Log to FILE (default: *STDOUT)>
		
		--logpath=PATH <Write logfile to PATH (default: .)>
		
		--logappend= <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl cappris.pl

	--id=ENSG00000093072.11
	
	--inpath=/home/jmrodriguez/projects/Encode/gencode15/annotations/chr22/ENSG00000093072.11
	
	--outfile=/home/jmrodriguez/projects/Encode/gencode15/annotations/chr22/appris.runtime.txt
	
	--loglevel=INFO
	
	--logappend
	
	--logpath=/home/jmrodriguez/projects/Encode/gencode15/logs/
			
	--logfile=cappris.log


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
