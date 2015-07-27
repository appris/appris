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
	$OLD_CFG
);

$CONFIG_INI_APPRIS_DIFF_FILE	= $FindBin::Bin.'/conf/apprisdiff.ini';

# Input parameters
my ($data_file) = undef;
my ($translations_file) = undef;
my ($apprisdiff_conf_file) = undef;
my ($output_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'data=s'				=> \$data_file,
	'translations=s'		=> \$translations_file,
	'apprisdiff-conf=s'		=> \$apprisdiff_conf_file,	
	'output=s'				=> \$output_file,
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);

# Required arguments
unless ( defined $data_file and defined $translations_file and defined $output_file )
{
	print `perldoc $0`;
	exit 1;
}

# Get old version of gencode to be compared
unless ( defined $apprisdiff_conf_file ) {
	$apprisdiff_conf_file = $CONFIG_INI_APPRIS_DIFF_FILE;
}
$OLD_CFG = new Config::IniFiles( -file =>  $apprisdiff_conf_file );

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
sub _copy_last_release($);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	my ($gene_chages) = '';
	my ($output) = '';
		
	# create gencode object
	$logger->info("-- create gencode object\n");
	my ($gencode_data) = parse_gencode($data_file, undef, $translations_file);
	unless ( defined $gencode_data ) {
		$logger->error("can not create gencode object\n");
	}
	
	# create old version of gencode object
	$logger->info("-- create old gencode object\n");
	my ($old_registry) = APPRIS::Registry->new(
										-dbhost  => $OLD_CFG->val( 'OLD_REGISTRY_DB', 'host'),
										-dbname  => $OLD_CFG->val( 'OLD_REGISTRY_DB', 'db' ),
										-dbuser  => $OLD_CFG->val( 'OLD_REGISTRY_DB', 'user'),
	);
	
	foreach my $gene (@{$gencode_data}) {
		my ($gene_id) = $gene->stable_id;		
		my ($old_gene_id) = $gene_id;
		$old_gene_id =~ s/\.\d*$//; # delete suffix because older versions don't use it
		my ($change,$translate,$trans)= (0,0,'');		
		$logger->debug("-- $gene_id:");
		my ($old_gene);
		eval {
			$old_gene = $old_registry->fetch_basic_by_stable_id('gene', $old_gene_id);				
		};
		$change = 1 if ($@); # New gene
		if ( $change == 0 ) {
			($change, $translate,$trans) = _changes_between_releases($gene, $old_gene);
		}
		$logger->info("$gene_id:$change:$translate:$trans\n");
		if ( $change != 0 ) {
			$output .= $gene_id."\t".$change."\n";			
		}
		
		# Copy the results of gene if there is any change
		if ( ($change == 0) and ($translate == 1) ) {
			# genes that have changed
		}
	}
	
	# Print list of genes with their change code
	$logger->info("-- print list of genes\n");	
	my ($printing_file_log) = printStringIntoFile($output, $output_file);
	$logger->error("printing output") unless ( defined $printing_file_log );		
		
	$logger->finish_log();
	
	exit 0;	
}

# Get the transcript that are different between releases
sub _changes_between_releases($$)
{
	my ($gene, $old_gene) = @_;
	my ($change) = 0;
	my ($translate) = 0;
	my ($trans) = '';
	my ($num_translations) = 0;
	my ($old_num_translations) = 0;
	
#$logger->debug("GENE:\n".Dumper($gene)."\n");
#$logger->debug("OLD_GENE:\n".Dumper($old_gene)."\n");
	
	foreach my $transcript (@{$gene->transcripts}) {			
		my ($transcript_id) = $transcript->stable_id;
		my ($old_transcript_id) = $transcript_id;
		$old_transcript_id =~ s/\.\d*$//; # delete suffix because older versions don't use it
		$logger->debug("$transcript_id: ");
		
		if ( $transcript->translate ) {
			$translate = 1; # Gene has got one translation at least. This case is not discriminate			
			my ($translation) = $transcript->translate;
									 
				if ( exists $old_gene->{'_index_transcripts'}->{$old_transcript_id} ) {
					my ($index) = $old_gene->{'_index_transcripts'}->{$old_transcript_id};
					my ($old_transcript) = $old_gene->transcripts->[$index];

					if ( $old_transcript->translate ) {
						my ($old_translation) = $old_transcript->translate;						

						if ( $translation->cds and $old_translation->cds and 
							(scalar(@{$translation->cds}) == scalar(@{$old_translation->cds}))
						) {
							my ($cds_list) = $translation->cds;
							my ($old_cds_list) = $old_translation->cds;
							
							for ( my $i = 0; $i < scalar(@{$cds_list}); $i++ ) {						
								if ( exists $cds_list->[$i] and $old_cds_list->[$i]) {
									unless ( ($cds_list->[$i]->start == $old_cds_list->[$i]->start) or
											 ($cds_list->[$i]->end == $old_cds_list->[$i]->end)
									) {
										$change = 2;  # Changes when cds coord of a transcript has changed
										$trans = $old_transcript_id;
										last;
									}
								}
								else {
									$change = 3; # Change when a transcript has new cds
									$trans = $old_transcript_id;
									last;
								}
							}
						}
						else {
							$change = 4; # Change when a transcript has new cds
							$trans = $old_transcript_id;
							last;
						}
				}			
				else {
					$change = 5; # Change when a transcript with no protein
					$trans = $old_transcript_id;
					last;
				}
			}			
			else {
				$change = 6; # Change when there is a new transcript
				$trans = $old_transcript_id;
				last;
			}
			$num_translations++;
		}
	}
	
	if ( ($change == 0) and ($translate == 1) ) {
		foreach my $old_transcript (@{$old_gene->transcripts}) {			
			$old_num_translations++ if ( $old_transcript->translate );
		}
		if ( $num_translations != $old_num_translations ) {
			$change = 7; # Change when new gene lost proteins
		}
	}

	$logger->debug("\n");
	return ($change,$translate,$trans);
}

main();


1;

__END__

=head1 NAME

get_diff_genes

=head1 DESCRIPTION

Get the list of genes which are different among versions of GENCODE.

=head2 Required arguments:

	--data=  <Gencode data file>
	
	--translations=  <Gencode translations file>
	
	--apprisdiff-conf <Config file that has the old databases info>
		
	--output= <File where genes with differences will save>
		$change = 1 # new gene
		$change = 2 # when cds coord of a transcript has changed
		$change = 3 # when a transcript has new cds
		$change = 4 # when a transcript has new cds
		$change = 5 # when a transcript with no protein
		$change = 6 # when there is a new transcript
		$change = 7 # when new gene lost proteins	

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl get_diff_genes.pl
		
	--output=../features/rel7_vs_rel11.genes.txt
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
