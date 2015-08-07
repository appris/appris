#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Utils::File qw( getTotalStringFromFile getStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$READTHROUGH_TRANSL_GENCODE
	$NUM_TRANSL_GENES_GENCODE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
my ($file_readthrough) = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/features/homo_sapiens/e76_g20/readthrough_transcript.txt');
while ( $file_readthrough =~  /([^\.]*)\.[0-9]*/g ) {
	my ($t_id) = $1;
	$t_id =~ s/\s*//g;
	$READTHROUGH_TRANSL_GENCODE->{$t_id} = 1;
	
}
my ($file_num_transl) = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/features/homo_sapiens/e79_g22/num_transl_per_gene.txt');
while ( $file_num_transl =~  /(\d+)\t*([^\n]*)/g ) {
	my ($g_id) = $2;
	my ($num) = $1;
	$g_id =~ s/\s*//g;
	$num =~ s/\s*//g;
	$NUM_TRANSL_GENES_GENCODE->{$g_id} = $num;
	
}

# Input parameters
my ($appris_infile) = undef;
my ($cano_infile) = undef;
my ($uni_infile) = undef;
my ($tsl_infile) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'in-appris=s'		=> \$appris_infile,
	'in-cano=s'			=> \$cano_infile,
	'in-uniprot=s'		=> \$uni_infile,
	'in-tsl=s'			=> \$tsl_infile,
	'outfile=s'			=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $appris_infile and defined $cano_infile and defined $uni_infile and defined $tsl_infile and defined $outfile )
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
$logger->init_log();


#####################
# Method prototypes #
#####################
sub get_appris_uniprot_tsl_report($$$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	$logger->debug("-- get appris/uniprot/tsl report -------\n");
	my ($report) = get_appris_uniprot_tsl_report($appris_infile, $cano_infile, $uni_infile, $tsl_infile);
	$logger->debug("REPORT:\n".Dumper($report)."\n");

	$logger->info("-- print appris/uniprot/tsl report -------\n");
	my ($output_content) = get_report_content($report);	
	if ($output_content ne '') {
		my ($printing_file_log) = printStringIntoFile($output_content, $outfile);
		$logger->error("printing") unless(defined $printing_file_log);		
	}

}

sub get_appris_uniprot_tsl_report($$$$) {
	my ($afile, $cfile, $ufile, $tfile) = @_;
	my ($report);
	my ($t_g_report);
	my ($c_report);
	
	my ($acontent) = getTotalStringFromFile($afile);	
	foreach my $line (@{$acontent}) {
		if ( defined $line and ($line ne '') ) {
			my (@cols) = split("\t", $line);
			my ($g_id) = $cols[1];
			my ($t_id) = $cols[2];
			my ($c_id) = $cols[3];
			my ($label) = $cols[4];
			$g_id =~ s/\s*//mg; $g_id =~ s/\.\d*$//;
			$t_id =~ s/\s*//mg; $t_id =~ s/\.\d*$//;
			$label =~ s/\s*//mg;
			$report->{$g_id}->{$t_id}->{'appris'} = $label;
			$report->{$g_id}->{$t_id}->{'ccds_id'} = $c_id;
			$t_g_report->{$t_id} = $g_id;
		}
	}
	my ($ccontent) = getTotalStringFromFile($cfile);
	foreach my $c_id (@{$ccontent}) {
		$c_id =~ s/\s*//mg;
		$c_report->{$c_id} = 1;
	}	
	my ($ucontent) = getTotalStringFromFile($ufile);
	foreach my $line (@{$ucontent}) {
		if ( defined $line and ($line ne '') ) {
			my (@cols) = split("\t", $line);
			my ($u_id) = $cols[0];
			my ($t_id) = $cols[1];
			$u_id =~ s/\s*//mg;
			$t_id =~ s/\s*//mg; $t_id =~ s/\.\d*$//;
			if ( exists $t_g_report->{$t_id} ) {
				my ($g_id) = $t_g_report->{$t_id};
				$report->{$g_id}->{$t_id}->{'uniprot'} = $u_id;
				if ( exists $c_report->{$u_id} ) {
					$report->{$g_id}->{$t_id}->{'canonical'} = 1;
				}				
			}
		}
	}
	my ($tcontent) = getTotalStringFromFile($tfile);	
	for (my $i = 1; $i <= scalar(@{$tcontent}); $i++) {
		my ($line) = $tcontent->[$i];
		if ( defined $line and ($line ne '') ) {
			my (@cols) = split("\t", $line);
			my ($g_nm) = $cols[0];
			my ($g_id) = $cols[1];
			my ($t_id) = $cols[2];
			my ($label) = $cols[3];
			$g_id =~ s/\s*//mg; $g_id =~ s/\.\d*$//;
			$t_id =~ s/\s*//mg; $t_id =~ s/\.\d*$//;
			if ( $label =~ s/([^\s]*)\s// ) {
				$report->{$g_id}->{$t_id}->{'tsl'} = $1;				
			}
		}
	}
	return $report;
}

sub get_report_content($) {
	my ($report) = @_;
	my ($output) = '';	
	foreach my $g_id ( sort( {$a cmp $b} keys(%{$report})) ) {
		# skip genes with unique transcript
		if ( exists $NUM_TRANSL_GENES_GENCODE->{$g_id} and ($NUM_TRANSL_GENES_GENCODE->{$g_id} == 1) ) {
			next;
		}		
		my ($g_rep) = $report->{$g_id};
		foreach my $t_id ( sort( {$a cmp $b} keys(%{$g_rep})) ) {			
			# skip read-through transcripts
			if ( exists $READTHROUGH_TRANSL_GENCODE->{$t_id} ) {
				next;
			}
			my ($t_rep) = $g_rep->{$t_id};
			if ( exists $t_rep->{'appris'} or exists $t_rep->{'uniprot'} ) {
				my ($a_label) = '-';
				my ($ccds_id) = '-';
				my ($u_label) = '-';
				my ($c_label) = '-';
				my ($t_label) = '-';		
				$a_label = $t_rep->{'appris'} if ( exists $t_rep->{'appris'} );
				$ccds_id = $t_rep->{'ccds_id'} if ( exists $t_rep->{'ccds_id'} );
				$u_label = $t_rep->{'uniprot'} if ( exists $t_rep->{'uniprot'} );
				$c_label = $t_rep->{'canonical'} if ( exists $t_rep->{'canonical'} );
				$t_label = $t_rep->{'tsl'} if ( exists $t_rep->{'tsl'} );
				$output .= $g_id."\t".$t_id."\t".$a_label."\t".$ccds_id."\t".$u_label."\t".$c_label."\t".$t_label."\n";
			}
		}
	}
	return $output;
}

main();


1;

__END__

=head1 NAME

cmpUniProtvsAPPRIS

=head1 DESCRIPTION

Compare UniProt anntos agains APPRIS 

=head1 SYNOPSIS

retrieve_exon_data

=head2 Required arguments:

	--in-appris= <APPRIS Principal Isoforms>

	--in-tsl=  <TSL data>
		
	--outfile= <Output file>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl cmpUniProtvsAPPRIS.pl
		
	--in-appris=/Users/jmrodriguez/projects/APPRIS/data/homo_sapiens/ens79.v8.2Apr2015/appris_data.principal.txt
	
	--in-cano=/Users/jmrodriguez/projects/APPRIS/workspaces/uniprot/uniprot_sprot_human.canonical.2015_06.dat
	
	--in-uniprot=/Users/jmrodriguez/projects/APPRIS/workspaces/uniprot/UniProt_Ensembl_idmapping.2015_06.dat
	
	--in-tsl=/Users/jmrodriguez/projects/APPRIS/workspaces/tsl/TSL.annots.e79_g21.txt
	
	--outfile=/Users/jmrodriguez/projects/APPRIS/workspaces/uniprot/cmpUniProtvsAPPRIS.ens79.txt
	
	--loglevel=debug --logfile=/Users/jmrodriguez/projects/APPRIS/workspaces/uniprot/cmpUniProtvsAPPRIS.log
		

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
