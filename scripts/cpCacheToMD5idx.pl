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

#my ($fasta_dir) = '/Users/jmrodriguez/projects/APPRIS/data/';
my ($fasta_dir) = '/home/jmrodriguez/projects/APPRIS/data/';
my ($ALL_FILES) = [{
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
	'ds'		=> 'e84'
}];

my ($OLDCACHEDIR) = '/home/jmrodriguez/projects/APPRIS/appriscache';

sub create_muldir($);

my ($report);
foreach my $inrep (@{$ALL_FILES}) {
	my ($file) = $inrep->{'transl'};
	my ($species) = $inrep->{'species'};
	my ($ds) = $inrep->{'ds'};
	
print STDERR "-- $species:$ds\n";

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
		my ($idx_dir) = $cache->c_dir();
		
		# create metadata and seq
		my ($metadata) = $ENV{APPRIS_WS_NAME}."\t".$seq_id."\n";
		my ($add_meta) = $cache->add_data($metadata,'meta');
		
		# extract seq identifiers
		my (@ids) = split('\|', $seq_id);
		my ($tranc_id) = $ids[0];
		my ($tranl_id) = $ids[1];
		my ($gene_id) = $ids[2];
		
		
		# copy old cache files, if thet do not exits
		my (@files) = (
			$idx_dir.'/seq.faa',
			$idx_dir.'/seq.psi',
			$idx_dir.'/seq.firestar',
			$idx_dir.'/seq.pfam',
			$idx_dir.'/seq.refseq',
			$idx_dir.'/seq.phobius',
			$idx_dir.'/seq.prodiv',
			$idx_dir.'/seq.memsat',
		);
		my ($shortage);
		foreach my $f ( @files ) {
			unless ( -e $f ) { $shortage = 1; }			
		}
		if ( defined $shortage ) {
			my ($old_cache_gdir) = $OLDCACHEDIR.'/'.$species.'/'.$ds.'/'.$gene_id;
			eval {
				my ($cmd) = "cd $old_cache_gdir/cache && tar -cf - $tranc_id.* | (cd $idx_dir; tar -xf -)";
				info($cmd);
				system($cmd);
			};
			throw('copying the global files of cache') if ($@);
			eval {
				my ($cmd) = "cd $old_cache_gdir/firestar && tar -cf - $tranc_id.out | (cd $idx_dir; tar -xf -)";
				info($cmd);
				system($cmd);
			};
			throw('copying cache file of firestar') if ($@);
		}
		
		# Change the names: transcript_id -> seq
		eval {
			my ($cmd) = "mv $idx_dir/$tranc_id.faa      $idx_dir/seq.faa & ".
						"mv $idx_dir/$tranc_id.psi      $idx_dir/seq.psi & ".
						"mv $idx_dir/$tranc_id.hhr      $idx_dir/seq.hhr & ".
						"mv $idx_dir/$tranc_id.out      $idx_dir/seq.firestar & ".
						"mv $idx_dir/$tranc_id.pdb      $idx_dir/seq.pdb & ".
						"mv $idx_dir/$tranc_id.pfam     $idx_dir/seq.pfam & ".
						"mv $idx_dir/$tranc_id.refseq   $idx_dir/seq.refseq & ".
						"mv $idx_dir/$tranc_id.phobius  $idx_dir/seq.phobius & ".
						"mv $idx_dir/$tranc_id.prodiv   $idx_dir/seq.prodiv & ".
						"mv $idx_dir/$tranc_id.memsat   $idx_dir/seq.memsat  ";
						
			info($cmd);
			system($cmd);
		};
		throw('changing the names') if ($@);						
		
		# Replace the transcript id to Query name
		eval {
			my ($cmd) = "cd $idx_dir && sed -i -- 's/$tranc_id/Query/g' seq.faa seq.psi seq.hhr seq.firestar seq.pdb seq.pfam seq.refseq seq.phobius seq.prodiv seq.memsat ";
			info($cmd);
			system($cmd);
		};
		throw('replacing the id') if ($@);
		
	}	
}
