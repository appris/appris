#!/usr/bin/perl -w

use strict;
use warnings;
use threads;

use Getopt::Long qw(:config no_auto_abbrev);
use FindBin;
use Data::Dumper;

use lib "$FindBin::Bin/lib";
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_ENSEMBL_DB_FILE
	$CONFIG_INI_APPRIS_DB_FILE
	$LOGLEVEL
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_ENSEMBL_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/ensembldb.ini';
$LOGLEVEL					= 'INFO';

# Input parameters
my ($conf_species) = undef;
my ($methods) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'conf|c=s'			=> \$conf_species,	
	'methods|m=s'		=> \$methods,
	'outfile|o=s'		=> \$outfile,
	'loglevel|l=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $conf_species and defined $methods ) {
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


# Main subroutine
sub main()
{
	# get the list of input genes with translation
	my (@data_genes);
	eval {
		my ($cmd) = "appris_check_ls_report -c $conf_species -g";
		$logger->debug("\n** script: $cmd\n");
		@data_genes = `$cmd`; map { s/\s+$// } @data_genes;
	};
	$logger->error("checking data genes") if($@);
	
	# get the list of input transcripts with translation
	my (@data_transc);
	eval {		
		my ($cmd) = "appris_check_ls_report -c $conf_species -t";
		$logger->debug("\n** script: $cmd\n");
		@data_transc = `$cmd`; map { s/\s+$// } @data_transc;
	};
	$logger->error("checking data transcripts") if($@);	

	# get the list of annotation files
	my (@annot_list);
	my ($annot_appris);
	eval {
		my ($cmd) = "appris_check_ls_annots -c $conf_species ";
		$logger->debug("\n** script: $cmd\n");
		@annot_list = `$cmd`; map { s/\s+$// } @annot_list;
	};
	$logger->error("deleting annot files") if($@);
	my ($methods_patt) = $methods; ($methods_patt =~ s/\,/\|/g);
	foreach my $annot_file (@annot_list) {
		if ( $annot_file =~ /^(.*)\.($methods_patt)$/ ) {
			my ($gene_id) = $1;
			my ($method) = $2;
			$annot_appris->{$method}->{$gene_id} = 1;
		}
	}	
	
	# compare the data gen set with the list of genes with annotatios
	my ($output) = '';
	foreach my $method ( split(',', $methods) ) {
		my (@diffs) = grep(!defined($annot_appris->{$method}->{$_}), @data_genes);		
		map( $output .= "[TRACELOG]\t". uc($method) . "\t" . $_ . "\n" , @diffs);		
	}
	print STDOUT $output;

$logger->debug("DATA_GENES:\n".Dumper(@data_genes));
#$logger->debug("DATA_TRANS:\n".Dumper(@data_transc));
#$logger->debug("ANNOT_LIST:\n".Dumper(@annot_list));
#$logger->debug("ANNOT_APPRIS:\n".Dumper($annot_appris));

	$logger->finish_log();
	
	exit 0;	
}

main();


1;

__END__

=head1 NAME

check_files

=head1 DESCRIPTION

Compare number the APPRIS result files and GENCODE/Ensembl data geneset  

=head1 SYNOPSIS

check_files

=head2 Required arguments (data input):

  -s, --species {string} <Name of species>
	* Hsap  - Homo sapiens -\n
	* Mmus  - Mus musculus -\n
	* Rnor  - Rattus norvegicus -\n
	* Drer  - Danio rerio -\n
	* Sscr  - Sus scrofa -\n
	* Dmel  - Drosophila melanogaster -\n
	* Cele  - Caenorhabditis elegans -\n
	* Lpar  - Lynx pardinus -\n  
  
 Or
 
  -c, --conf {file} <Config file name>
  
  -m, --methods= <List of APPRIS's methods ('firestar,matador3d,spade,corsair,thump,crash,appris')>  

=head2 Optional arguments (log arguments):
	
  --loglevel=LEVEL <define log level (default: NONE)>
  
  --logfile=FILE <Log to FILE (default: *STDOUT)>
  
  --logpath=PATH <Write logfile to PATH (default: .)>
  
  --logappend= <Append to logfile (default: truncate)>

=head1 EXAMPLE

	check_files.pl -s Hsap

=head1 EXAMPLE2

	check_files.pl -c conf/script/apprisrc.XXX
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
