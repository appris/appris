#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Config::IniFiles;
use FindBin;

use APPRIS::Registry;
use APPRIS::Exporter;
use APPRIS::Utils::File qw( updateStringIntoFile printStringIntoFile );
use APPRIS::Utils::Logger;

use Data::Dumper;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_APPRIS_DB_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_APPRIS_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/apprisdb.ini';

# Input parameters
my ($species) = undef;
my ($chr) = undef;
my ($methods) = undef;
my ($format) = undef;
my ($output_file) = undef;
my ($apprisdb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'species=s'			=> \$species,
	'chr=s'				=> \$chr,
	'methods=s'			=> \$methods,
	'format=s'			=> \$format,
	'output=s'			=> \$output_file,
	'apprisdb-conf=s'	=> \$apprisdb_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);
unless(defined $species and defined $format and defined $output_file)
{
	print `perldoc $0`;
	exit 1;
}

# Optional arguments
# get vars of appris db
unless ( defined $apprisdb_conf_file ) {
	$apprisdb_conf_file = $CONFIG_INI_APPRIS_DB_FILE;
}

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
sub get_data_by_chr($$);
sub get_data_by_method($$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Check inputs
	if (defined $format and ($format ne '')) {
		$format = lc($format);
		unless (($format ne '') and ( ($format eq 'gtf') or ($format eq 'gff3') or ($format eq 'bed') or ($format eq 'bed12') or ($format eq 'json') ) ) {
			print `perldoc $0`;
			exit 1;
		}
	}	
	
	# APPRIS Registry from given specie
	my ($cfg) = new Config::IniFiles( -file => $apprisdb_conf_file );
	my ($spe) = $species; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	my ($param) = {
			'-dbhost'       => $cfg->val('APPRIS_DATABASES', 'host'),
			'-dbuser'       => $cfg->val('APPRIS_DATABASES', 'user'),
			'-dbpass'       => $cfg->val('APPRIS_DATABASES', 'pass'),
			'-dbport'       => $cfg->val('APPRIS_DATABASES', 'port'),
			'-dbname'       => $cfg->val($specie_db, 'db'),
	};
	$logger->debug(Dumper($param)."\n");
	my ($registry) = APPRIS::Registry->new();
	$registry->load_registry_from_db(%{$param});
		
	# Get annotations from given chr
	$logger->info("-- get annotations from given chromosome: $chr\n");
	my ($output_report) = get_data_by_chr($registry, $chr);
	$logger->debug("REPORT\n".Dumper($output_report)."\n");
	
	# Print annotations
	if ( defined $output_report ) {
		while ( my ($met,$met_out) = each(%{$output_report}) ) {
			if ( $met_out ne '' ) {
				my ($met_file) = $output_file.'.'.$met.'.'.$format;
				my ($printing_file_log) = printStringIntoFile($met_out,$met_file);
				die ("ERROR: printing $chr:$met ") unless(defined $printing_file_log);						
			}				
		}
	}
	
	$logger->finish_log();
	
	exit 0;	
	
}
sub get_data_by_chr($$)
{
	my ($registry, $chr) = @_;
	my ($output_report);
	my ($chromosome) = $registry->fetch_by_region($chr, undef, undef, 'all');
	#$logger->debug(Dumper($chromosome)."\n");

	foreach my $method ( split(',',$methods) ) {
		$logger->debug("#method $method -------\n");
		$output_report->{$method} = get_data_by_method($chromosome,$method);
	}
	return $output_report;
}
sub get_data_by_method($$)
{
	my ($chromosome, $method) = @_;
	my ($output) = '';

	foreach my $gene (@{$chromosome}) {
		my ($gene_id) = $gene->stable_id;
		my ($exporter) = APPRIS::Exporter->new();
		if ($gene and ($format eq 'bed' or $format eq 'bed12' ) ){
			$output .= $exporter->get_bed_annotations($gene, $method, undef, $format);
		}
		elsif ($gene and $format eq 'gtf') {
			$output .= $exporter->get_gtf_annotations($gene, $method);
	    }
		elsif ($gene and $format eq 'gff3') {
			$output .= $exporter->get_gff3_annotations($gene, $method);
	    }
		elsif ($gene and $format eq 'json') {
			$output .= $exporter->get_json_annotations($gene, $method);
    	}
	}
	return $output;
}

main();


1;

__END__

=head1 NAME

retrieve_method_data

=head1 DESCRIPTION

Get the annotations of methods in several formats

=head1 SYNOPSIS

retrieve_method_data

=head2 Required arguments:

	--format <Output format: 'gtf', 'bed', or 'json'>
	
	--output <Output file that has method's annotations>
	
=head2 Optional arguments:

	--chr  <Genomic region>
	
	--methods= <List of APPRIS's methods ('firestar,matador3d,spade,corsair,thump,crash,appris'. Default: ALL)>

	--apprisdb-conf <Config file of APPRIS database>
				
=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl retrieve_method_data.pl

		--chr=chr21
		
		--methods=firestar,matador3d,spade,corsair,thump,crash,appris

		--format=gtf
		
		--output=../features/data/appris.results.rel7.v1.chr21
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
