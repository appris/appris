#!/usr/bin/perl -w
# _________________________________________________________________
# $Id$
# $Revision$
# Developed by:
#		Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es-
# _________________________________________________________________

use strict;
use warnings;
use FindBin;

use Digest::MD5;
use Data::Dumper;
use APPRIS::Utils::WSpace;
use APPRIS::Parser;

#
# This script is used when we have not time to execute APPRIS and we need the results.
# For example, when they give you a new assembly.
# It creates preliminar results from the last assembly.
#
# It works as...
# 1) When we have already the gene with at least one isoform with the same protein sequences, 
# the annotation of APPRIS is migrated.
# 
# 2) When we don't have APPRIS annotations for a gene => We get the longest sequences (PRINCIPAL:5)
#
# Note: It checks is there is a new isoform with the same protein sequences as older isoform.
#
# Run Exampl:
# perl createPreliminarAPPRIS.pl
#	~/projects/APPRIS/apprisws/data/danio_rerio/ens77.v7.9Feb2015/appris_data.principal.txt
#	~/projects/APPRIS/apprisws/features/danio_rerio/e77/Danio_rerio.Zv9.pep.all.fa
#	~/projects/APPRIS/apprisws/features/danio_rerio/e80/genes_translation.txt
#	~/projects/APPRIS/apprisws/features/danio_rerio/e80/Danio_rerio.GRCz10.pep.all.fa
#
#	1> /Users/jmrodriguez/projects/APPRIS/apprisws/data/danio_rerio/ens80.v8.16Apr2015/Danio_rerio.GRCz10.e80.appris_data.principal.txt
#	2> createPreliminarAPPRIS.log
#

my ($file1) = $ARGV[0];
my ($file2) = $ARGV[1];
my ($file3) = $ARGV[2];
my ($file4) = $ARGV[3];

local(*FILE);
open(FILE,$file1) or return undef;
my(@string1)=<FILE>;
close(FILE);

local(*FILE);
open(FILE,$file2) or return undef;
my(@string2)=<FILE>;
close(FILE);

local(*FILE);
open(FILE,$file3) or return undef;
my(@string3)=<FILE>;
close(FILE);

local(*FILE);
open(FILE,$file4) or return undef;
my(@string4)=<FILE>;
close(FILE);

# list of APPRIS genes for Zv9
my ($report);
foreach my $line (@string1) {
	if ( defined $line and ($line ne '') ) {
		my (@cols) = split("\t", $line);
		my ($g_id) = $cols[1];
		my ($t_id) = $cols[2];
		my ($label) = $cols[4];
		$g_id =~ s/\s*//mg; $g_id =~ s/\.\d*$//;
		$t_id =~ s/\s*//mg; $t_id =~ s/\.\d*$//;
		$label =~ s/\s*//mg;
		if ( $label =~ /PRINCIPAL/ ) {			
#			push(@{$report->{$g_id}->{'appris'}}, {
#				'label' 	=> $label,
#				't_id'		=> $t_id,
#			});
			$report->{$g_id}->{'appris'}->{$t_id} = {
				'label' 	=> $label,
			};
		}
	}
}
# Create Index using proteim sequences.
my ($sequences) = APPRIS::Parser::_parse_inseq_transl($file2);
while (my ($g_id,$g_rep) = each(%{$report}) ) {
	while (my ($t_id,$t_rep) = each(%{ $g_rep->{'appris'} }) ) {
		if ( exists $sequences->{$t_id} ) {		
			my ($seq) = $sequences->{$t_id};
			my ($len) = length($seq);
			# create id (md5)
			my ($idx);
			eval {
				my ($ctx) = Digest::MD5->new;
				$ctx->add($seq);
				$idx = $ctx->hexdigest;		
			};
			throw('Creating md5') if ($@);
			if ( defined $idx ) {
				$report->{$g_id}->{'appris'}->{$t_id}->{'idx'} = $idx;
				$report->{$g_id}->{'appris'}->{$t_id}->{'len'} = $len;
			}
		}		
	}
}
print STDERR "REPORT\n".Dumper($report)."\n";

# list of genes for Zv10
my ($report2);
foreach my $line (@string3) {
	if ( defined $line and ($line ne '') ) {
		my (@cols) = split(" ", $line);
		my ($g_id) = $cols[0];
		my ($t_id) = $cols[1];
		$g_id =~ s/\s*//mg; $g_id =~ s/\.\d*$//;
		$t_id =~ s/\s*//mg; $t_id =~ s/\.\d*$//;
		$report2->{$g_id}->{'appris'} = undef;
		$report2->{$g_id}->{'transc'}->{$t_id} = undef;		
	}
}
# Create Index using proteim sequences.
my ($sequences2) = APPRIS::Parser::_parse_inseq_transl($file4);
while (my ($g_id,$g_rep) = each(%{$report2}) ) {
	while (my ($t_id,$t_rep) = each(%{$g_rep->{'transc'}}) ) {
		if ( exists $sequences2->{$t_id} ) {		
			my ($seq) = $sequences2->{$t_id};
			my ($len) = length($seq);
			# create id (md5)
			my ($idx);
			eval {
				my ($ctx) = Digest::MD5->new;
				$ctx->add($seq);
				$idx = $ctx->hexdigest;		
			};
			throw('Creating md5') if ($@);
			if ( defined $idx ) {
				$report2->{$g_id}->{'transc'}->{$t_id} = {
					'idx'	 	=> $idx,
					'len'		=> $len,
				};
			}
		}		
	}
}
print STDERR "REPORT2\n".Dumper($report2)."\n";

# Check the same genes comparing Zv10 and Zv9
# If the translation exists, we use this isoform.
# Otherwise, we take the longest => PRINCIPAL:5
while (my ($g_id2,$g_rep2) = each(%{$report2}) ) {
	my ($a_tid_idx, $a_label);
	my ($t_rep2) = $g_rep2->{'transc'};
	my (@s_transc) = sort {$t_rep2->{$b}->{'len'} <=> $t_rep2->{$a}->{'len'}} keys(%{$t_rep2});
	my ($max_len) = $t_rep2->{ $s_transc[0] }->{'len'};
	
	# get annot from last results
	if ( exists $report->{$g_id2} and exists $report->{$g_id2}->{'appris'} and defined $report->{$g_id2}->{'appris'} ) {
		my (@a_tids) = keys( %{$report->{$g_id2}->{'appris'}} );
		my ($a_tid) = $a_tids[0];
		$a_tid_idx = $report->{$g_id2}->{'appris'}->{$a_tid}->{'idx'};
		$a_label = $report->{$g_id2}->{'appris'}->{$a_tid}->{'label'};
	}
	
	while (my ($t_id2,$rep2) = each(%{$t_rep2}) ) {
		my ($idx2) = $rep2->{'idx'};
		my ($len2) = $rep2->{'len'};
		if ( defined $a_tid_idx and defined $a_label and ($idx2 eq $a_tid_idx) ) {
			$report2->{$g_id2}->{'appris'}->{$t_id2} = {
					'label' 	=> $a_label,
					't_id'		=> $t_id2,
			};
		}
		elsif ( defined $max_len and ($len2 == $max_len) ) {
			$report2->{$g_id2}->{'appris'}->{$t_id2} = {
					'label' 	=> 'PRINCIPAL:5',
					't_id'		=> $t_id2,
			};			
		}	
	}
}
print STDERR "CONSENSUS\n".Dumper($report2)."\n";

# Print out
while (my ($g_id2,$g_rep2) = each(%{$report2}) ) {
	while (my ($t_id2,$t_rep2) = each(%{$g_rep2->{'appris'}}) ) {
		print STDOUT $g_id2."\t".$t_id2."\t".$t_rep2->{'label'}."\n";		
	}
}
