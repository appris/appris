#!/usr/bin/perl -W

use strict;
use warnings;
use Getopt::Long;
use FindBin;

use APPRIS::Registry;
use APPRIS::Utils::File qw( printStringIntoFile getTotalStringFromFile );
use APPRIS::Utils::Logger;

use Data::Dumper;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_APPRIS_DB_FILE
	
	$CFG
	$SPECIE_LIST
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_APPRIS_DB_FILE	= $LOCAL_PWD.'/conf/apprisdb.ALL.ini';

# Input parameters
my ($methods) = undef;
my ($output_file) = undef;
my ($apprisdb_in_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'methods=s'			=> \$methods,
	'output=s'			=> \$output_file,
	'apprisdb-conf=s'	=> \$apprisdb_in_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);
unless (
	defined $methods and
	defined $output_file
) {
	print `perldoc $0`;
	exit 1;
}

# Optional arguments
# get vars of appris db
unless ( defined $apprisdb_in_conf_file ) {
	$apprisdb_in_conf_file = $CONFIG_INI_APPRIS_DB_FILE;
}

# get config vars
$CFG = new Config::IniFiles( -file => $apprisdb_in_conf_file );
$SPECIE_LIST = [ split( ',', $CFG->val('APPRIS_DATABASES', 'species') ) ];

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
sub get_princ_for_genename($$\$);
sub get_princ_for_genename_datafile($$\$);
sub get_annot_for_genename($$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# get principal isoform for each gene (sorting by gene name) 	
	$logger->info("-- get principal isoform for each gene (sorting by gene name) -------\n");
	my ($princ_gname);
	foreach my $specie ( @{$SPECIE_LIST} )
	{
		# get principal file for specie
		$logger->info("-- -- get data files for: $specie\n");
		my ($specie_files) = uc($specie.'_files');
		my ($princ_data_file) = $CFG->val($specie_files, 'appris_data_dir').'/'.$CFG->val($specie_files, 'appris_data_principal');

		# get the genes for every specie and save them by gene name
		$logger->info("-- -- get_princ_for_genename: $princ_data_file\n");
		get_princ_for_genename_datafile($princ_data_file, $specie, $princ_gname);
		$logger->debug(Dumper($princ_gname)."\n");
				
#		# APPRIS registry
#		$logger->info("-- -- get appris registry for: $specie\n");		
#		my ($specie_db) = uc($specie.'_db');
#		my ($registry) = APPRIS::Registry->new();		
#		$registry->load_registry_from_db(
#								-dbhost	=> $CFG->val($specie_db, 'host'),
#								-dbname	=> $CFG->val($specie_db, 'db'),
#								-dbuser	=> $CFG->val($specie_db, 'user'),
#								-dbpass	=> $CFG->val($specie_db, 'pass'),
#								-dbport	=> $CFG->val($specie_db, 'port'),
#		);
#		$logger->debug(Dumper($registry)."\n");
#		
#		# get the genes for every specie and save them by gene name
#		$logger->info("-- -- get_princ_for_genename\n");
#		get_princ_for_genename($registry, $specie, $princ_gname);
#		$logger->debug(Dumper($princ_gname)."\n");
	}

	# get appris methods for each principal (sorting by gene name)	
	$logger->info("-- get principal isoform for each gene (sorting by gene name) -------\n");
	my ($met_gname);
	while ( my ($gene_name, $princ_rep) = each(%{$princ_gname}) )
	{
		$logger->debug("-- -- all: ".scalar(keys(%{$princ_rep}))." >= ".(scalar(@{$SPECIE_LIST})-1)."\n");
		# get use only common gene names (at least in two species)
		if ( scalar(keys(%{$princ_rep})) >= (scalar(@{$SPECIE_LIST})-1)  )
		{

			$logger->debug("-- -- $gene_name\n");
			while ( my ($specie, $g_report) = each(%{$princ_rep}) ) {
				# APPRIS registry
				$logger->info("-- -- get appris registry for: $specie\n");
				my ($specie_db) = uc($specie.'_db');
				my ($registry) = APPRIS::Registry->new();		
				$registry->load_registry_from_db(
										-dbhost	=> $CFG->val($specie_db, 'host'),
										-dbname	=> $CFG->val($specie_db, 'db'),
										-dbuser	=> $CFG->val($specie_db, 'user'),
										-dbpass	=> $CFG->val($specie_db, 'pass'),
										-dbport	=> $CFG->val($specie_db, 'port'),
				);
				$logger->debug(Dumper($registry)."\n");				
				while ( my ($gene_id, $g_rep) = each(%{$g_report}) ) {
					my ($transc_id) = $g_rep->{'transcript_id'};
					my ($met_report) = get_annot_for_genename($registry, $transc_id);
					if ( defined $met_report ) {
						$met_report->{'transcript_id'} = $transc_id;
						$met_gname->{$gene_name}->{$specie}->{$gene_id} = $met_report;
					}
				}				
			}
			
		}
	}
	$logger->debug(Dumper($met_gname)."\n");	

	# print annots sorting by gene name
	my ($output_content) = '';
	while ( my ($gene_name, $met_rep) = each(%{$met_gname}) ) {
		my ($g_out_cont);
		foreach my $specie ( @{$SPECIE_LIST} ) {
			if ( exists $met_rep->{$specie} ) {
				my ($g_report) = $met_rep->{$specie};
				while ( my ($gene_id, $g_rep) = each(%{$g_report}) ) {
					my ($g_output) =	$gene_id."\t";
					foreach my $met ( split(',',$methods) ) {
						if ( defined $met ) {
							$g_output .=	$g_rep->{$met}."\t";						
						}
					}
					push(@{$g_out_cont},$g_output);
				}				
			}
			else {
				my ($g_output) = '-'."\t";
				foreach my $met ( split(',',$methods) ) {
					$g_output .= '-'."\t";
				}
				push(@{$g_out_cont},$g_output);				
			}
		}
		my ($out_content) = '';
		foreach my $out (@{$g_out_cont}) {
			$out_content .=	$out;
		}
		$out_content =~ s/\t$//;
		$output_content .=	$gene_name."\t".$out_content."\n";
	}	
	if ($output_content ne '') {
		my ($printing_file_log) = printStringIntoFile($output_content,$output_file);
		die ("ERROR: printing ") unless(defined $printing_file_log);		
	}

	$logger->finish_log();
	
	exit 0;	
	
}
sub get_princ_for_genename_datafile($$\$)
{
	my ($data_file, $specie, $report) = @_;
	my ($data_cont) = getTotalStringFromFile($data_file);
	
	foreach my $input_line (@{$data_cont}) {
		$input_line=~s/\n*$//;
		my(@split_line)=split("\t", $input_line);		
		next if(scalar(@split_line)<=0);
		
		my ($chr) = $split_line[0];
		my ($gene_id) = $split_line[1];
		my ($transc_id) = $split_line[2];
		my ($transc_name) = $split_line[3];
		my ($ccds_id) = $split_line[4];
		my ($appris_annot) = $split_line[5];
		
		my ($gene_name) = uc($transc_name);
		$gene_name =~ s/\-(\d*)$//;
		$logger->debug("-- $gene_id: $gene_name\n");
				
		my ($g_report) = {
					'chr'					=> $chr,
					'gene_id'				=> $gene_id,
					'transcript_id'			=> $transc_id,
					'transcript_name'		=> $transc_name,
					'ccds_id'				=> $ccds_id,
					'appris_annot'			=> $appris_annot,				
		};
		${$report}->{$gene_name}->{$specie}->{$gene_id} = $g_report;
	}
}

sub get_princ_for_genename($$\$)
{
	my ($registry, $specie, $report) = @_;
	my ($chromosome) = $registry->fetch_basic_by_region(undef);
	
	foreach my $gene (@{$chromosome}) {
		my ($gene_id) = $gene->stable_id;
		my ($gene_name) = $gene->external_name;
		my ($chr) = $gene->chromosome;
		$logger->debug("-- $gene_id: ");
		
		my ($g_report);
		foreach my $transcript (@{$gene->transcripts}) {		
			my ($transcript_id) = $transcript->stable_id;			
			my ($transcript_name) = '-';
			if ($transcript->external_name) {
				$transcript_name = $transcript->external_name;	
			}

			if ($transcript->translate) {
				my ($translation_length) = length($transcript->translate->sequence);
				my ($ccds_id) = '-';				
				if ( $transcript->xref_identify ) {
					foreach my $xref_identify (@{$transcript->xref_identify}) {								
						if ($xref_identify->dbname eq 'CCDS') {
							$ccds_id = $xref_identify->id;
						}
					}					
				}
				
				my ($analysis) = $registry->fetch_analysis_by_stable_id($transcript_id,'appris');				
				if ( $analysis and $analysis->appris and $analysis->appris->principal_isoform_signal ) {
					my ($appris_annot) = $analysis->appris->principal_isoform_signal;
					if ( $appris_annot eq 'YES' ) {
						$g_report->{'chr'}					= $chr;
						$g_report->{'gene_id'}				= $gene_id;
						$g_report->{'transcript_id'}		= $transcript_id;
						$g_report->{'transcript_name'}		= $transcript_name;
						$g_report->{'ccds_id'}				= $ccds_id;
						$g_report->{'translation_length'}	= $translation_length;
						$g_report->{'appris_annot'}			= 'PRINCIPAL';
					}
					elsif ( $appris_annot eq 'UNKNOWN' ) {
						unless ( defined $g_report ) {
							$g_report->{'chr'}					= $chr;
							$g_report->{'gene_id'}				= $gene_id;
							$g_report->{'transcript_id'}		= $transcript_id;
							$g_report->{'transcript_name'}		= $transcript_name;
							$g_report->{'ccds_id'}				= $ccds_id;
							$g_report->{'translation_length'}	= $translation_length;
							$g_report->{'appris_annot'}			= 'LONGEST';
						}
						else {
							if ( $translation_length > $g_report->{'translation_length'} ) { # maximum length
								$g_report->{'chr'}					= $chr;
								$g_report->{'gene_id'}				= $gene_id;
								$g_report->{'transcript_id'}		= $transcript_id;
								$g_report->{'transcript_name'}		= $transcript_name;
								$g_report->{'ccds_id'}				= $ccds_id;
								$g_report->{'translation_length'}	= $translation_length;
								$g_report->{'appris_annot'}			= 'LONGEST';
							}
						}
					}
				}
			}
		}
		if ( defined $gene_name and defined $g_report ) {
			${$report}->{$gene_name}->{$specie}->{$g_report->{'gene_id'}} = $g_report;
		}
		$logger->debug("\n");
	}
}

sub get_annot_for_genename($$)
{
	my ($registry, $transc_id) = @_;
	my ($met_report);
	my ($analysis) = $registry->fetch_analysis_by_stable_id($transc_id,'appris');

	if ( $analysis and $analysis->appris ) {
		foreach my $method ( split(',',$methods) ) {
			my ($annot);
			if ( $method eq 'firestar' and defined $analysis->appris->functional_residues_score ) {
				$annot = $analysis->appris->functional_residues_score;					
			}
			elsif ( $method eq 'matador3d' and defined $analysis->appris->homologous_structure_score ) {
				$annot = $analysis->appris->homologous_structure_score;					
			}
			elsif ( $method eq 'corsair' and defined $analysis->appris->vertebrate_conservation_score ) {
				$annot = $analysis->appris->vertebrate_conservation_score;					
			}
			elsif ( $method eq 'spade' and defined $analysis->appris->domain_score ) {
				$annot = $analysis->appris->domain_score;					
			}
			if ( defined $annot ) {
				$met_report->{$method} = $annot;				
			}
		}
	}
	
	return $met_report;
}

main();


1;

__END__

=head1 NAME

retrieve_homolog_data

=head1 DESCRIPTION

Get the list of genes with the same name between species.

=head1 SYNOPSIS

retrieve_homolog_data

=head2 Required arguments:

	--methods= <List of APPRIS's methods (default: ALL)>

	--output <Output file that has the main isoforms>	
	
=head2 Optional arguments:

	--apprisdb-conf <Config file of APPRIS database>
				
=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl retrieve_homolog_data.pl
	
		--apprisdb-conf=conf/apprisdb.ALL.ini
		
		--methods=firestar,matador3d		
		
		--output=stats/appris_data.same_genename.g15-v4_Mmus70-v4_Rnor70-v4.txt
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
