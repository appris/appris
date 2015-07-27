#!/usr/bin/perl -W

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Registry;
use APPRIS::CDS;
use APPRIS::Parser qw( parse_gencode );
use APPRIS::Utils::File qw( printStringIntoFile );
use APPRIS::Utils::Logger;


###################
# Global variable #
###################
use vars qw(
	$CONFIG_INI_APPRIS_DIFF_FILE
	$CFG
	$OLD_LIST_GENES
);

$CONFIG_INI_APPRIS_DIFF_FILE	= $FindBin::Bin.'/conf/apprisdiff.ini';

# Input parameters
my ($old_data_file) = undef;
my ($new_data_file) = undef;
my ($new_score_file) = undef;
my ($apprisdiff_conf_file) = undef;
my ($output_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'old-data=s'			=> \$old_data_file,
	'new-data=s'			=> \$new_data_file,
	'new-score=s'			=> \$new_score_file,
	'apprisdiff-conf=s'		=> \$apprisdiff_conf_file,	
	'output=s'				=> \$output_file,
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);

# Required arguments
unless ( defined $old_data_file and defined $new_data_file and defined $new_score_file and defined $output_file )
{
	print `perldoc $0`;
	exit 1;
}

# Get old version of gencode to be compared
unless ( defined $apprisdiff_conf_file ) {
	$apprisdiff_conf_file = $CONFIG_INI_APPRIS_DIFF_FILE;
}
$CFG = new Config::IniFiles( -file =>  $apprisdiff_conf_file );

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();


#####################
# Method prototypes #
#####################
sub _changes_between_releases($$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	my ($gene_chages) = '';
	my ($new_output) = '';
	my ($chan_output) = '';
	
	# create gencode object (old)
	my ($old_genes);
	my ($ocmd) = "awk '{if (\$3 ==\"gene\") {print \$10\"\t\"\$18} }' $old_data_file | sed -e 's/[\"|\;]//g'";
	$logger->info("** script; $ocmd\n");
	my (@old_gene_list) = `$ocmd`;
	for my $line (@old_gene_list) {
		if ( $line =~ /([^\.]*)\.([^\t]*)\t([^\n]*)\n/ ) {
			$old_genes->{$1} = {
				'version'	=> $2,
				'name'		=> $3
			};			
		}
	}
$logger->debug("-- OLD_GENES\n".Dumper($old_genes)."\n");	

	# get the list of genes (new)
	my ($new_genes);
	my ($ncmd) = "awk '{if (\$3 ==\"gene\") {print \$10\"\t\"\$18} }' $new_data_file | sed -e 's/[\"|\;]//g'";
	$logger->info("** script: $ncmd\n");
	my (@new_gene_list) = `$ncmd`;
	for my $line (@new_gene_list) {
		if ( $line =~ /([^\.]*)\.([^\t]*)\t([^\n]*)\n/ ) {
			$new_genes->{$1} = {
				'version'	=> $2,
				'name'		=> $3
			};			
		}
	}
$logger->debug("-- NEW_GENES\n".Dumper($new_genes)."\n");
	
	
	# acquire differences from new gene list
	while ( my ($new_gene_id, $new_gene_rep) = each(%{$new_genes}) ) {
		my ($change) = 0;
		my ($new_gene_name) = $new_gene_rep->{'name'};
		if ( exists $old_genes->{$new_gene_id} ) {
			my ($old_gene_rep) = $old_genes->{$new_gene_id};
			if ( $new_gene_rep->{'version'} ne $old_gene_rep->{'version'} ) {
				$change = 2; # Changed gene based on ensembl version
			}
		}
		else {
			$change = 1; # New gene
		}
	$logger->debug("-- CHANGE: $new_gene_id > $change\n");
		if ( ($change == 1) or ($change == 2) ) {
			my ($max_fire, $max_mat, $max_cor, $max_spa, $max_thump) = ('-', '-', '-', '-', '-');
			my ($cmd) = "grep $new_gene_id $new_score_file";
			my (@appris_scores_list) = `$cmd`;
$logger->debug("-- APPRIS_SCORE: $new_gene_id\n".Dumper(@appris_scores_list)."\n");
			for my $line (@appris_scores_list) {
				my (@lines) = split ("\t", $line);
				my ($g_id) = $lines[0];
				my ($t_id) = $lines[1];
				my ($fire) = $lines[7];
				my ($mat) = $lines[8];
				my ($cor) = $lines[9];
				my ($spa) = $lines[10];
				my ($thump) = $lines[13];
				my ($c_sp) = $lines[14];
				my ($c_tp) = $lines[15];

				if ( $max_fire eq '-' ) { $max_fire = $fire; }
				else { $max_fire = $fire if ( $fire > $max_fire ); }				

				if ( $max_mat eq '-' ) { $max_mat = $mat; }
				else { $max_mat = $mat if ( $mat > $max_mat ); }				

				if ( $max_cor eq '-' ) { $max_cor = $cor; }
				else { $max_cor = $cor if ( $cor > $max_cor ); }				

				if ( $max_spa eq '-' ) { $max_spa = $spa; }
				else { $max_spa = $spa if ( $spa > $max_spa ); }

				if ( $max_thump eq '-' ) { $max_thump = $thump; }
				else { $max_thump = $thump if ( $thump > $max_thump ); }
			}
			if ( $change == 1 ) {
				my ($t) = $new_gene_id."\t".$new_gene_name."\t".$max_fire."\t".$max_mat."\t".$max_cor."\t".$max_spa."\t".$max_thump;
$logger->debug("-- OUTPUT: $t\n");
				$new_output .= $t."\n";				
			}
			elsif ( $change == 2 ) {
				my ($t) = $new_gene_id."\t".$new_gene_name."\t".$max_fire."\t".$max_mat."\t".$max_cor."\t".$max_spa."\t".$max_thump;
$logger->debug("-- OUTPUT_2: $t\n");
				$chan_output .= $t."\n";				
			}			
		}
		
	}
	
	# Print list of genes with their change code
	$logger->info("-- print list of genes\n");	
	my ($printing_file_log) = printStringIntoFile($new_output, $output_file.'.new.txt');
	$logger->error("printing output") unless ( defined $printing_file_log );		
	my ($printing_file_log2) = printStringIntoFile($chan_output, $output_file.'.changed.txt');
	$logger->error("printing output") unless ( defined $printing_file_log2 );
		
	$logger->finish_log();
	
	exit 0;	
}


main();


1;

__END__

=head1 NAME

diff_genes

=head1 DESCRIPTION

Get the list of genes which are different among versions of GENCODE.

=head2 Required arguments:

	--old-data=  <Gencode data file>
	
	--new-data=  <Gencode data file>
	
	--new-scores=  <Score file>
	
	--output= <File where genes with differences will save>
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl diff_genes.pl
		
	--old-data /home/jmrodriguez/projects/Encode/appris_rel12/features/gencode.v12.annotation.gtf
	
	--new-data /home/jmrodriguez/projects/Encode/gencode19/features/gencode.v19.annotation.gtf
	
	--new-score /home/jmrodriguez/projects/Encode/gencode19/data/gen19.v2.18Dec2013/appris_data.prin_scores.txt
	
	--output /home/jmrodriguez/projects/Encode/gencode19/data/gen19.v2.18Dec2013/diff_genes.gen12_vs_gen19	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
