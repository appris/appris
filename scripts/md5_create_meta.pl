#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Digest::MD5;
use Bio::SeqIO;
use Data::Dumper;

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::File qw( updateStringIntoLockFile printStringIntoFile prepare_workspace );
use APPRIS::Utils::Exception qw(throw info warning deprecate);

my ($fasta_dir) = '/Users/jmrodriguez/projects/APPRIS/data/';
#my ($fasta_dir) = '/home/jmrodriguez/projects/APPRIS/data/';
my ($ALL_FILES) = [{
# ALREADY DONE:
	'transl'	=> $fasta_dir.'homo_sapiens/gen19.v17/appris_data.transl.fa',
	'species'	=> 'homo_sapiens',
	'ds'		=> 'gencode19'
},{
	'transl'	=> $fasta_dir.'homo_sapiens/ens84.v17/appris_data.transl.fa',
	'species'	=> 'homo_sapiens',
	'ds'		=> 'e84_g24'	
},{
	'transl'	=> $fasta_dir.'homo_sapiens/rs105.v17/appris_data.transl.fa',
	'species'	=> 'homo_sapiens',
	'ds'		=> 'rs105'	
},{
	'transl'	=> $fasta_dir.'homo_sapiens/rs107.v17/appris_data.transl.fa',
	'species'	=> 'homo_sapiens',
	'ds'		=> 'rs107'	
},{
	'transl'	=> $fasta_dir.'homo_sapiens/up201606.v18/appris_data.transl.fa',
	'species'	=> 'homo_sapiens',
	'ds'		=> 'up201606'
},{
	'transl'	=> $fasta_dir.'mus_musculus/ens84.v17/appris_data.transl.fa',
	'species'	=> 'mus_musculus',
	'ds'		=> 'e84_gM9'
},{
	'transl'	=> $fasta_dir.'danio_rerio/ens84.v17/appris_data.transl.fa',
	'species'	=> 'danio_rerio',
	'ds'		=> 'e84'
},{
	'transl'	=> $fasta_dir.'danio_rerio/ens77.v17/appris_data.transl.fa',
	'species'	=> 'danio_rerio',
	'ds'		=> 'e77'
},{
## NEW:
	'transl'	=> $fasta_dir.'rattus_norvegicus/ens84.v17/appris_data.transl.fa',
	'species'	=> 'rattus_norvegicus',
	'ds'		=> 'e84'
},{
	'transl'	=> $fasta_dir.'rattus_norvegicus/ens77.v17/appris_data.transl.fa',
	'species'	=> 'rattus_norvegicus',
	'ds'		=> 'e77'
},{
	'transl'	=> $fasta_dir.'sus_scrofa/ens84.v17/appris_data.transl.fa',
	'species'	=> 'sus_scrofa',
	'ds'		=> 'e84'
},{
	'transl'	=> $fasta_dir.'pan_troglodytes/ens84.v17/appris_data.transl.fa',
	'species'	=> 'pan_troglodytes',
	'ds'		=> 'e84'
},{
	'transl'	=> $fasta_dir.'caenorhabditis_elegans/ens84.v17/appris_data.transl.fa',
	'species'	=> 'caenorhabditis_elegans',
	'ds'		=> 'e84'
},{
	'transl'	=> $fasta_dir.'drosophila_melanogaster/ens84.v17/appris_data.transl.fa',
	'species'	=> 'drosophila_melanogaster',
	'ds'		=> 'e84'
}];

#$ENV{APPRIS_PROGRAMS_CACHE_DIR} = '/home/jmrodriguez/projects/APPRIS/appriscache/md5';
$ENV{APPRIS_PROGRAMS_CACHE_DIR} = '/Users/jmrodriguez/projects/APPRIS/appris/cache/md5';

sub main() {	
	my ($report);
	#my ($all_metadata) = '';
	foreach my $inrep (@{$ALL_FILES}) {
		my ($file) = $inrep->{'transl'};
		my ($species) = $inrep->{'species'};
		my ($ds) = $inrep->{'ds'};
		my ($ws_name) = "$species/$ds";		
		info("-- $ws_name\n");
		
		my ($in) = Bio::SeqIO->new(
					-file => $file,
					-format => 'Fasta'
		);
		while ( my $seqObj = $in->next_seq() ) {
			# create cache obj
			my ($seq_id) = $seqObj->id;
			my ($seq_s) = $seqObj->seq;		
			my ($cache) = APPRIS::Utils::CacheMD5->new(
				-dat => $seq_s,
				-ws  => $ENV{APPRIS_PROGRAMS_CACHE_DIR}			
			);		
			my ($seq_idx) = $cache->idx;
			my ($seq_sidx) = $cache->sidx;
			my ($seq_dir) = $cache->dir;
					
			# prepare cache dir
			#my ($idx_dir) = $cache->c_dir();
			
			# create metadata and seq
			my ($metadata) = $seq_sidx."\t".$ws_name."\t".$seq_id."\n";
			#$all_metadata .= $metadata;
			#my ($add_meta) = $cache->add_data($metadata,'meta');
print STDOUT "$metadata";
								
		}
	}
	
}

main();
