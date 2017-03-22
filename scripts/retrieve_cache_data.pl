#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::File qw(
	prepare_workspace
	printStringIntoFile
	getTotalStringFromFile
);

use lib "$FindBin::Bin/lib";
use common qw( get_main_report get_label_report );

###################
# Global variable #
###################

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($input_main_file) = undef;
my ($input_seq_file) = undef;
my ($methods) = undef;
my ($outdir) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'input-main=s' 			=> \$input_main_file,
	'input-seq=s'    		=> \$input_seq_file,
	'methods=s'				=> \$methods,
	'outdir=s'				=> \$outdir,
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);
unless ( defined $input_main_file and defined $input_seq_file and defined $methods and defined $outdir )
{
	print `perldoc $0`;
	exit 1;
}

# Optional arguments

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
sub get_cache_data($$$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Get the method report
	$logger->info("-- get method report -------\n");
	my (%met_report) = map { $_ => 1 } split(',', $methods);

	# Prepare workspaces
	$logger->info("-- prepare workspaces -------\n");
	$outdir = prepare_workspace($outdir);
	foreach my $met ( keys(%met_report) ) { prepare_workspace($outdir.$met) }

	# Get data from file
	$logger->info("-- get main data from files -------\n");
	my ($main_report, $seq_report) = common::get_main_report($input_main_file, $input_seq_file);	
	#$logger->debug("MAIN_REPORT:\n".Dumper($main_report)."\n");
	#$logger->debug("SEQ_REPORT:\n".Dumper($seq_report)."\n");
	
	# Get annots
	$logger->info("-- get cache data -------\n");
	get_cache_data($main_report, $seq_report, \%met_report);

	$logger->finish_log();
	
	exit 0;	
	
}
sub get_cache_data($$$)
{
	my ($main_report, $seq_report, $met_report) = @_;
	
	# get appris annot
	while (my ($gene_id, $g_report) = each(%{$main_report}) )
	{
		while (my ($transc_id, $t_report) = each(%{$g_report->{'transcripts'}}) ) {
			# get sequence
			my ($translation_seq) = '';
			if (exists $seq_report->{$gene_id}->{'transcripts'}->{$transc_id}) {
				$translation_seq = $seq_report->{$gene_id}->{'transcripts'}->{$transc_id};
			}
			
			# get cache id
			if ( $translation_seq ne '' ) {
				my ($cache) = APPRIS::Utils::CacheMD5->new(
					-dat => $translation_seq,
					-ws  => $ENV{APPRIS_PROGRAMS_CACHE_DIR}			
				);		
				my ($seq_dir) = $cache->dir;
				
				# copy cache files
				my ($method) = '';
				$method = 'firestar';
				if ( exists $met_report->{$method} ) {
					my ($outdir_met) = $outdir.'/'.$method;
					eval {
						my ($cmd) = "cp -r $seq_dir/seq.firestar $outdir_met/$transc_id.firestar && ".
									"cp -r $seq_dir/seq.psi      $outdir_met/$transc_id.psi && ".
									"cp -r $seq_dir/seq.hhr      $outdir_met/$transc_id.hhr ";
						$logger->debug("\n** script: $cmd\n");
						system ($cmd);
					};
					throw("copying cache files of $method") if($@);
				}
				$method = 'matador3d';
				if ( exists $met_report->{$method} ) {
					my ($outdir_met) = $outdir.'/'.$method;
					eval {
						my ($cmd) = "cp -r $seq_dir/seq.pdb       $outdir_met/$transc_id.pdb && ".
									"cp -r $seq_dir/seq.pdb70.hmm $outdir_met/$transc_id.pdb70.hmm && ".
									"cp -r $seq_dir/seq.matador3d $outdir_met/$transc_id.matador3d ";
						$logger->debug("\n** script: $cmd\n");
						system ($cmd);
					};
					throw("copying cache files of $method") if($@);
				}
				$method = 'spade';
				if ( exists $met_report->{$method} ) {
					my ($outdir_met) = $outdir.'/'.$method;
					eval {
						my ($cmd) = "cp -r $seq_dir/seq.pfam $outdir_met/$transc_id.pfam ";
						$logger->debug("\n** script: $cmd\n");
						system ($cmd);
					};
					throw("copying cache files of $method") if($@);
				}
				$method = 'corsair';
				if ( exists $met_report->{$method} ) {
					my ($outdir_met) = $outdir.'/'.$method;
					eval {
						my ($cmd) = "cp -r $seq_dir/seq.refseq $outdir_met/$transc_id.refseq ";
						$logger->debug("\n** script: $cmd\n");
						system ($cmd);
					};
					throw("copying cache files of $method") if($@);
				}
				$method = 'thump';
				if ( exists $met_report->{$method} ) {
					my ($outdir_met) = $outdir.'/'.$method;
					eval {
						my ($cmd) = "cp -r $seq_dir/seq.phobius    $outdir_met/$transc_id.phobius && ".
									"cp -r $seq_dir/seq.prodiv     $outdir_met/$transc_id.prodiv && ".
									"cp -r $seq_dir/seq.memsat     $outdir_met/$transc_id.memsat && ".
									"cp -r $seq_dir/seq.chk_swtr90 $outdir_met/$transc_id.chk_swtr90 && ".
									"cp -r $seq_dir/seq.swtr90     $outdir_met/$transc_id.swtr90 && ".
									"cp -r $seq_dir/seq.mtx        $outdir_met/$transc_id.mtx && ".
									"cp -r $seq_dir/seq.memsat_aln $outdir_met/$transc_id.memsat_aln && ".
									"cp -r $seq_dir/seq.kalign     $outdir_met/$transc_id.kalign ";
						$logger->debug("\n** script: $cmd\n");
						system ($cmd);
					};
					throw("copying cache files of $method") if($@);
				}
				$method = 'crash';
				if ( exists $met_report->{$method} ) {
					my ($outdir_met) = $outdir.'/'.$method;
					eval {
						my ($cmd) = "cp -r $seq_dir/seq.signalp    $outdir_met/$transc_id.signalp && ".
									"cp -r $seq_dir/seq.targetp    $outdir_met/$transc_id.targetp ";
						$logger->debug("\n** script: $cmd\n");
						system ($cmd);
					};
					throw("copying cache files of $method") if($@);
				}				
			}
						
		} # end transc-loop
	}
	return undef;
}

main();


1;

__END__

=head1 NAME

retrieve_cache_data

=head1 DESCRIPTION

Get the cache data from given method.

=head1 SYNOPSIS

retrieve_cache_data

=head2 Required arguments:

	--input-main <Result file of APPRIS's scores>

	--input-seq <Result file of APPRIS protein sequences>
	
	--methods= <List of APPRIS's methods ('firestar,matador3d,spade,corsair,thump,crash,appris')>	

	--outdir    <Output dir>
	
=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl retrieve_cache_data.pl
	
		--input-main=data/appris_data.appris.txt
		
		--input-seq=data/appris_data.transl.fa
		
		--methods='spade,corsair'

		--outdir=cache_data/
		
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
