#!/usr/bin/perl -W

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
my ($tsl_infile) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'in-appris=s'		=> \$appris_infile,
	'in-tsl=s'			=> \$tsl_infile,
	'outfile=s'			=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $appris_infile and defined $tsl_infile and defined $outfile )
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
sub get_appris_tsl_report($$);
#sub get_appris_report($);
#sub get_tsl_report($);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	$logger->debug("-- get appris/tsl report -------\n");
	my ($report) = get_appris_tsl_report($appris_infile, $tsl_infile);
	$logger->debug("REPORT:\n".Dumper($report)."\n");

	$logger->info("-- print appris/tsl report -------\n");
	my ($output_content) = get_report_content($report);	
	if ($output_content ne '') {
		my ($printing_file_log) = printStringIntoFile($output_content, $outfile);
		$logger->error("printing") unless(defined $printing_file_log);		
	}

}

sub get_appris_tsl_report($$) {
	my ($afile,$tfile) = @_;
	my ($report);
	
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
	my ($num_total_principals) = {
		'PRINCIPAL:1' 	=> 0,
		'PRINCIPAL:2' 	=> 0,
		'PRINCIPAL:3' 	=> 0,
		'PRINCIPAL:4' 	=> 0,
		'PRINCIPAL:5' 	=> 0,
		'ALTERNATIVE:1' => 0,
		'ALTERNATIVE:2' => 0,
	};
	my ($num_total_tsl) = {			
			'tsl1'	=> 0,
			'tsl2' 	=> 0,
			'tsl3' 	=> 0,
			'tsl4' 	=> 0,
			'tsl5' 	=> 0,
			'tslNA'	=> 0,
	};	
	my ($num_principals);
	foreach my $p_id ( keys(%{$num_total_principals}) ) {
		$num_principals->{$p_id} = {			
			'tsl1'	=> 0,
			'tsl2' 	=> 0,
			'tsl3' 	=> 0,
			'tsl4' 	=> 0,
			'tsl5' 	=> 0,
			'tslNA'	=> 0,
		};
	}	
	
	foreach my $g_id ( sort( {$a cmp $b} keys(%{$report})) ) {
		# skip genes with unique transcript
		if ( exists $NUM_TRANSL_GENES_GENCODE->{$g_id} and ($NUM_TRANSL_GENES_GENCODE->{$g_id} == 1) ) {
			next;
		}
		my ($print_flag) = 0;
		my ($g_output) = '';
		my ($g_rep) = $report->{$g_id};
		my ($num_transc_tsl) = {			
				'tsl1'	=> 0,
				'tsl2' 	=> 0,
				'tsl3' 	=> 0,
				'tsl4' 	=> 0,
				'tsl5' 	=> 0,
				'tslNA'	=> 0,
		};		
		foreach my $t_id ( sort( {$a cmp $b} keys(%{$g_rep})) ) {			
			# skip read-through transcripts
			if ( exists $READTHROUGH_TRANSL_GENCODE->{$t_id} ) {
				next;
			}
			my ($t_rep) = $g_rep->{$t_id};
			my ($a_label) = '-';
			my ($t_label) = '-';		
			my ($ccds_id) = '-';
			if ( exists $t_rep->{'appris'} ) {
				my ($a) = $t_rep->{'appris'};
				$num_total_principals->{$a} += 1;
				if ( exists $t_rep->{'tsl'} ) {
					my ($t) = $t_rep->{'tsl'}; $t = ($t eq "" )? 'NA' : $t;
					#$num_total_tsl->{$t} += 1;
					$num_principals->{$a}->{$t} += 1;
				}				
			}
			if ( exists $t_rep->{'tsl'} ) {
				my ($t) = $t_rep->{'tsl'}; $t = ($t eq "" )? 'NA' : $t;
				$num_transc_tsl->{$t} += 1;
				$num_total_tsl->{$t} += 1;
			}			
			#if ( exists $t_rep->{'appris'} and exists $t_rep->{'tsl'} ) {
			#	if ( $t_rep->{'tsl'} eq "tsl1" ) {
			#		$t_label = $t_rep->{'tsl'};
			#		$print_flag = 1;
			#	}
			#	$a_label = $t_rep->{'appris'} if ( exists $t_rep->{'appris'} );
			#	$ccds_id = $t_rep->{'ccds_id'} if ( exists $t_rep->{'ccds_id'} );
			#	$g_output .= $g_id."\t".$t_id."\t".$a_label."\t".$ccds_id."\t".$t_label."\n";
			#}
			if ( exists $t_rep->{'tsl'} and ($t_rep->{'tsl'} eq "tsl1") ) {
				$t_label = $t_rep->{'tsl'};
				$a_label = $t_rep->{'appris'} if ( exists $t_rep->{'appris'} );
				$ccds_id = $t_rep->{'ccds_id'} if ( exists $t_rep->{'ccds_id'} );
				$print_flag = 1;
				$g_output .= $g_id."\t".$t_id."\t".$a_label."\t".$ccds_id."\t".$t_label."\n";
			}
		}
		# print TSL1 when is unique within a gene
		if ( ($print_flag == 1) and ($num_transc_tsl->{'tsl1'} >= 1) ) { $output .= $g_output }
	}
	
	# print SUMMARY
	my ($summary) = '';
	my ($s_head) .= "\t\t\t";
	my ($s_body) .= "\t\t\t";
	foreach my $t_id ( sort {$a cmp $b} keys(%{$num_total_tsl}) ) {
		$s_head .= $t_id."\t";
		$s_body .= $num_total_tsl->{$t_id}."\t";
	}
	$s_head =~ s/\t$//g; $s_body =~ s/\t$//g;
	$summary .= $s_head."\n".$s_body."\n\n";
	foreach my $p_id ( sort {$a cmp $b} keys(%{$num_principals}) ) {
		my ($p_report) = $num_principals->{$p_id};
		my ($total) = $num_total_principals->{$p_id};
		my ($s) = '';
		foreach my $t_id ( sort {$a cmp $b} keys(%{$p_report}) ) {
			my ($per) = sprintf( "%.2f", ($p_report->{$t_id}/$total)*100 );
			$s .= $per."\t";
		}
		$s =~ s/\t$/\n/g;
		$summary .= $p_id."\t".$total."\t".$s;
	}
	print STDOUT "\n\n$summary\n";



	return $output;
}

main();


1;

__END__

=head1 NAME

cmpTSLvsAPPRIS

=head1 DESCRIPTION

Compare TSL anntos agains APPRIS 

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

perl cmpTSLvsAPPRIS.pl
		
	--in-appris=/Users/jmrodriguez/projects/APPRIS/data/homo_sapiens/ens79.v8.2Apr2015/appris_data.principal.txt
	
	--in-tsl=/Users/jmrodriguez/projects/APPRIS/workspaces/tsl/TSL.annots.e79_g22.txt
	
	--outfile=/Users/jmrodriguez/projects/APPRIS/workspaces/tsl/cmpTSLvsAPPRIS.ens79.txt
	
	--loglevel=debug --logfile=/Users/jmrodriguez/projects/APPRIS/workspaces/tsl/cmpTSLvsAPPRIS.log
		

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
