#!/usr/bin/perl -w
# $Id: THUMPvX.pl (formerly known as TMHv13.pl)
# Developed by: Paolo Maietta -pmaietta@cnio.es-
# Created: 	20-May-2009
# Updated: 	20-Oct-2009
# reUpdated:	20-Aug-2010
# This version was developed with a parallel performing of Memsat and Sequence collecting
# _________________________________________________________________

use strict;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile getStringFromFile prepare_workspace );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$BIN_DIR
	$MEMSAT3_DIR
	$PHOBIUS_DIR
	$PRODIV_DIR
	
	$WSPACE_TMP
	$WSPACE_CACHE
	$PROG_DB
	
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL
	
	$LOGGER_CONF
);

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($input_file) = undef;
my ($output_file) = undef;
my ($appris) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'conf=s'			=> \$config_file,
	'input=s'			=> \$input_file,
	'output=s'			=> \$output_file,
	'appris'			=> \$appris,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $config_file and defined $input_file and defined $output_file )
{
    print `perldoc $0`;
    exit 1;
}

# Get conf vars
my ($cfg)			= new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$BIN_DIR			= $LOCAL_PWD.'/bin/';
$MEMSAT3_DIR		= $ENV{APPRIS_PROGRAMS_OPT_DIR}.'/memsat3/';
$PHOBIUS_DIR		= $ENV{APPRIS_PROGRAMS_OPT_DIR}.'/phobius/';
$PRODIV_DIR			= $ENV{APPRIS_PROGRAMS_OPT_DIR}.'/prodiv/';
$WSPACE_TMP			= $ENV{APPRIS_TMP_DIR};
$WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
$PROG_DB			= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('THUMP_VARS', 'db');

$OK_LABEL			= 'YES';
$UNKNOWN_LABEL		= 'UNKNOWN';
$NO_LABEL			= 'NO';

$LOGGER_CONF		= '';

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);
$LOGGER_CONF .= " --loglevel=$loglevel " if ( defined $loglevel );
$LOGGER_CONF .= " --logpath=$logpath " if ( defined $logpath );
$LOGGER_CONF .= " --logfile=$logfile " if ( defined $logfile );
$LOGGER_CONF .= " --logappend " if ( defined $logappend );
my ($logfilename) = $logger->logpath().'/'.$logger->logfile();

#####################
# Method prototypes #
#####################
sub run_phobius($$$);
sub run_memsat($$$);
sub run_kalign($$$);
sub convert_mod($$$);
sub run_prodiv($$$);
sub park_consensus($$$$$$);
sub park_consensus2($$$$$$$$$$$);
sub filter_length($$$);
sub _get_appris_annotations($);
sub _delete_tmp_dir($);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Run methods
	$logger->info("-- run thump pipeline for every variant\n");
	eval {
		my ($in) = Bio::SeqIO->new(
					-file => $input_file,
					-format => 'Fasta'
		);
		while ( my $seqObj = $in->next_seq() ) {			
			if ( $seqObj->id =~ /([^|]*)/ )
			{			
				my ($seq_id) = $1;
				if ( $seq_id =~ /^ENS/ ) { $seq_id =~ s/\.\d*$// }	
				my ($seq) = $seqObj->seq;
				$logger->info("-- $seq_id\n");				
				
				# Create cache obj
				my ($cache) = APPRIS::Utils::CacheMD5->new( -dat => $seq );
				my ($seq_idx) = $cache->idx;
				my ($seq_idx_s) = $cache->idx_s;
				
				my ($ws_tmp) = $WSPACE_TMP.'/'.$seq_idx;
				my ($ws_cache) = $WSPACE_CACHE.'/'.$seq_idx_s;
				prepare_workspace($ws_tmp);
				prepare_workspace($ws_cache);
				
				# Cached fasta
				my ($seq_file) = $ws_cache.'/seq.faa';
				unless(-e $seq_file and (-s $seq_file > 0) ) {				
					my ($seq_cont) = ">Query\n$seq\n";			
					open (SEQ_FILE,">$seq_file");
					print SEQ_FILE $seq_cont;
					close (SEQ_FILE);
				}
			
				# Run methods
				$logger->info("-- run phobius\n");
				my ($output_phobius) = run_phobius($seq_file, $ws_cache, $ws_tmp);
				
				$logger->info("-- run memsat (with psiblast)\n");
				my ($output_memsat, $output_blast, $output_align) = run_memsat($seq_file, $ws_cache, $ws_tmp);
				
				$logger->info("-- run kalign\n");
				my ($output_kalign) = run_kalign($output_align, $ws_cache, $ws_tmp);
	
				$logger->info("-- convert mod\n");
				my ($output_mod) = convert_mod($output_kalign, $ws_cache, $ws_tmp);
			
				$logger->info("-- run prodiv\n");
				my ($output_prodiv) = run_prodiv($output_mod, $ws_cache, $ws_tmp);
			}
		}
	};
	$logger->error("scaning input") if($@);
	
	# Run Consensus for each seq (if not cached)
	$logger->info("-- run consensus method for every variant\n");
	my ($num_helices_file, $consensus_file, $thump_file);
	my ($in) = Bio::SeqIO->new(
				-file => $input_file,
				-format => 'Fasta'
	);
	while ( my $seqObj = $in->next_seq() ) {
		if ( $seqObj->id =~ /([^|]*)/ )
		{			
			my ($seq_id) = $1;
			if ( $seq_id =~ /^ENS/ ) { $seq_id =~ s/\.\d*$// }
			my ($seq) = $seqObj->seq;
			$logger->info("-- $seq_id\n");			
			
			# Create cache obj
			my ($cache) = APPRIS::Utils::CacheMD5->new( -dat => $seq );
			my ($seq_idx) = $cache->idx;
			my ($seq_idx_s) = $cache->idx_s;
			
			my ($ws_tmp) = $WSPACE_TMP.'/'.$seq_idx;
			my ($ws_cache) = $WSPACE_CACHE.'/'.$seq_idx_s;
			prepare_workspace($ws_tmp);
			prepare_workspace($ws_cache);
			
			# Cached fasta
			my ($seq_file) = $ws_cache.'/seq.faa';
			unless(-e $seq_file and (-s $seq_file > 0) ) {				
				my ($seq_cont) = ">Query\n$seq\n";			
				open (SEQ_FILE,">$seq_file");
				print SEQ_FILE $seq_cont;
				close (SEQ_FILE);
			}
			
			# Output files
			my ($phobius_file) = $ws_cache."/seq.phobius";
			my ($memsat_file) = $ws_cache."/seq.memsat";
			my ($prodiv_file) = $ws_cache."/seq.prodiv";

			$logger->info("-- parse consensus\n");
			my ($park_phobius_file, $park_prodiv_file, $park_memsat_file, $phobius, $prodiv, $memsat, $phobius2, $prodiv2, $memsat2);		
			eval {
				($park_phobius_file, $park_prodiv_file, $park_memsat_file, $phobius, $prodiv, $memsat, $phobius2, $prodiv2, $memsat2) = park_consensus($seq_id, $seq_file, $phobius_file, $memsat_file, $prodiv_file, $ws_tmp);
				$logger->debug("PHOBIUS:\n".Dumper($phobius)."\n");
				$logger->debug("PRODIV:\n".Dumper($prodiv)."\n");
				$logger->debug("MEMSAT:\n".Dumper($memsat)."\n");
				$logger->debug("PHOBIUS2:\n".Dumper($phobius2)."\n");
				$logger->debug("PRODIV2:\n".Dumper($prodiv2)."\n");
				$logger->debug("MEMSAT2:\n".Dumper($memsat2)."\n");
			};
			$logger->error("parsing consensus") if($@);
			
			my ($aux_consensus_file);				
			if ( defined $park_phobius_file and 
				defined $park_prodiv_file and
				defined $park_memsat_file and
				defined $phobius and
				defined $prodiv and
				defined $memsat and
				defined $phobius2 and
				defined $prodiv2 and
				defined $memsat2 )
			{
				eval {
					$aux_consensus_file = park_consensus2($park_phobius_file, $park_prodiv_file, $park_memsat_file, $phobius, $prodiv, $memsat, $phobius2, $prodiv2, $memsat2, $ws_cache, $ws_tmp);
				};
				$logger->error("parsing consensus2") if($@);					
			}
			
			if ( defined $aux_consensus_file )
			{		
				eval {
					($num_helices_file, $consensus_file, $thump_file) = filter_length($aux_consensus_file, $ws_cache, $ws_tmp);			
				};
				$logger->error("running filter_length") if($@);
			}
		}
	}		
	
	# Print records by transcript ---------------
	my ($output_content) = getStringFromFile($thump_file);
	unless( defined $output_content ) {
		$logger->error("Can not obtain output file: $!\n");
	}
	
	# Get the annotations for the main isoform /* APPRIS */ ----------------
	if ( defined $appris )
	{		
		$output_content .= _get_appris_annotations($output_content);
	}
	
	# Print records by transcript ---------------
	my ($print_out) = printStringIntoFile($output_content, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}
	
	$logger->finish_log();
	
	exit 0;
}

sub run_phobius($$$)
{
	my ($input, $ws_cache, $ws_tmp) = @_;
	my ($output) = $ws_cache.'/seq.phobius';	
	unless ( -e $output and (-s $output > 0) )
	{
		eval {
			my ($cmd) = "perl $BIN_DIR/phobvX.pl --input=$input --output=$output $LOGGER_CONF";
			$logger->debug("\n** script: $cmd\n");
			my (@out) = `$cmd`;
		};
		$logger->error("run_phobius") if($@);
	}
	return $output;
}

sub run_memsat($$$)
{
	my ($input, $ws_cache, $ws_tmp) = @_;
	my ($out_memsat) = $ws_cache.'/seq.memsat';
	my ($out_chk) = $ws_cache."/seq.chk_swtr90";
	my ($out_blast) = $ws_cache."/seq.swtr90";
	my ($out_align) = $ws_cache."/seq.memsat_aln";	
	unless ( -e $out_memsat and (-s $out_memsat > 0) and 
			 -e $out_blast and (-s $out_blast > 0) and
			 -e $out_align and (-s $out_align > 0) )
	{
		eval {
			my ($cmd) = "perl $BIN_DIR/memsatvX.pl ".
								" --db=$PROG_DB ".
								" --name=Query ".
								" --input=$input ".
								" --out-memsat=$out_memsat ".
								" --out-chk=$out_chk ".
								" --out-blast=$out_blast ".
								" --out-align=$out_align ".
								" --tmp-dir=$ws_tmp ".
								" --cache-dir=$ws_cache ".
								" $LOGGER_CONF ";
			$logger->debug("\n** script: $cmd\n");
			my (@out) = `$cmd`;
		};
		$logger->error("run_memsat") if($@);
	}
	return ($out_memsat, $out_blast, $out_align);
}

sub run_kalign($$$)
{
	my ($input, $ws_cache, $ws_tmp) = @_;	
	my ($out_align) = $ws_cache.'/seq.kalign';
	unless ( -e $out_align and (-s $out_align > 0) )
	{
		eval {
			my ($cmd) = "kalign -i $input -o $out_align -f clu -c input -d wu ";
			$logger->debug("\n** script: $cmd\n");
			my (@out) = `$cmd`;
		};
		$logger->error("run_kalign") if($@);
	}
	return ($out_align);
}

sub convert_mod($$$)
{
	my ($input, $ws_cache, $ws_tmp) = @_;
	my ($out_mod) = $ws_cache.'/seq.mod';
	eval {
		my (%seen);
		my (@trash);
		my ($k) = 0;
		open (FILE, $input);
		open (MOD, ">$out_mod");
		my ($switch1) = "on";
		my ($switch2) = "off";
		while (<FILE>){
			if ($_=~/^\n/ && $switch1 eq "off") {
				$switch2="on";
				$k=0;
			}
			if ($_=~/^.+\s+[A-Z-]+\n/) {
				$switch1 = "off";
				@trash= split(/\s+/,$_);
				if ($switch2 eq "on"){
					$k++;
					$seen{$k}=$seen{$k}.$trash[1];
				}
				else {
					$k++;
					$seen{$k}=$trash[1];
				}
			}
		}
		@trash=keys(%seen);
		for($k=1;$k<scalar@trash;$k++){
			my (@parking) = split(//, $seen{$k});
			my ($sequence) = "<";
			foreach my $w (@parking){
				$sequence = $sequence.$w.";";
			}
			$sequence = $sequence.">\n";
			print MOD $sequence;
		}
		close (MOD);
	};
	$logger->error("run_kalign") if($@);
	return ($out_mod);
}

sub run_prodiv($$$)
{
	my ($input, $ws_cache, $ws_tmp) = @_;
	my ($out_align) = $ws_cache.'/seq.prodiv';
	my ($in_tmp_prodiv) = $ws_tmp.'/seq.prodiv.res';	
	unless ( -e $out_align and (-s $out_align > 0) )
	{		
		eval {
			#my ($cmd) = "perl $PRODIV_DIR/all_tmhmm_runner.pl $id $WS_MOD $WS_PREDICTION $WS_TMP &>> $logfilename";
			#my ($cmd) = "perl $PRODIV_DIR/all_tmhmm_runner.pl $id $WS_MOD $WS_PREDICTION $WS_TMP 1> /dev/null 2> /dev/null";
			my ($cmd) = "perl $PRODIV_DIR/all_tmhmm_runner.pl seq $ws_cache $ws_tmp 1> /dev/null";
			$logger->debug("\n** script: $cmd\n");
			my (@out) = `$cmd`;
		};
		$logger->error("run_prodiv") if($@);

		$logger->debug("PASE_PRODIV\n");	
		eval {

			my ($output_prodiv) = parse_prodiv($in_tmp_prodiv,$out_align);		
		};
		$logger->error("parse_prodiv") if($@);
	}	
	return ($out_align);
}

sub parse_prodiv($$)
{
	my ($input, $output) = @_;
	
	my ($fin) = '';
	my ($caracteristic) = '';
	
	my (@listado);
	if (-e $input and (-s $input > 0) ) {
		open (FILE, $input) or die "Can not open $input: $!\n";
		@listado = <FILE>;
		close(FILE);			
	}
	for ( my $i = 0; $i <= scalar(@listado); $i++ ) {
		if ( defined $listado[$i] && ($listado[$i] =~ /^(.+).mod/) ){
			open (PRE,">$output");
			print PRE "ID   $1\n";
		}
		if ( defined $listado[$i] && ($listado[$i] =~ /^Labeling:/) ){
			$i++;
			my ($cadena) = '';
			while ( $listado[$i] !~ /Posterior probabilities:/ ) {
				chomp ($listado[$i]);
				$cadena.=$listado[$i];
				$i++;
			}
			my (@secuencias) = split(//, $cadena);

			my ($contador) = 1;
			my ($control) = '';
			my ($ops) = 0;
			foreach my $k (@secuencias){
				if ( $k eq 'i' ) {
					if ( $control ne 'C' && $ops==1 ) {
						$fin = $contador-1;
						print PRE "$fin\t$caracteristic\n";
						$control='C';
						$caracteristic='CYTOPLASMIC.';
						print PRE "FT   DOMAIN\t$contador\t";
					}
					elsif ($control ne 'C' && $ops==0){
						print PRE "FT   DOMAIN\t$contador\t";
						$control='C';
						$caracteristic='CYTOPLASMIC.';
						$ops=1;
					}
				}
				elsif ($k eq 'M'){
					if ($control ne 'T' && $ops==1) {
						$fin=$contador-1;
						print PRE "$fin\t$caracteristic\n";
						$control='T';
						$caracteristic='';
						print PRE "FT   TRANSMEM\t$contador\t";
					}
					elsif ($control ne 'T' && $ops==0) {
						print PRE "FT   TRANSMEM\t$contador\t";
						$control='T';
						$caracteristic='';
						$ops=1;
					}
				}
				elsif ($k eq 'o') {
					if ($control ne 'NC' && $ops==1) {
						$fin=$contador-1;
						print PRE "$fin\t$caracteristic\n";
						$control='NC';
						$caracteristic='NON CYTOPLASMIC.';
						print PRE "FT   DOMAIN\t$contador\t";
					}
					elsif ($control ne 'NC' && $ops==0){
						print PRE "FT   DOMAIN\t$contador\t";
						$control='NC';
						$caracteristic='NON CYTOPLASMIC.';
						$ops=1;
					}
				}
				$contador++;
			}
			$fin = $contador-1;
			print PRE "$fin\t$caracteristic\n";
			print PRE "//";
			close(PRE);
		}
	}
}

sub park_consensus($$$$$$) 
{
	my ($id, $seq_file, $phobius_file, $memsat_file, $prodiv_file, $ws_tmp) = @_;
	
	my ($phobius);
	my ($prodiv);
	my ($memsat);
	my ($phobius2);
	my ($prodiv2);
	my ($memsat2);
	
	my (@fichero);
	my ($contador);
	my (@sequence);
	my ($w);
	my ($control);
	
	
	local(*SEQ_FILE);
	open (SEQ_FILE, $seq_file);
	my (@parking) = <SEQ_FILE>;
	close(SEQ_FILE);

	my ($park_phobius_file) = $ws_tmp."/parking.phobius";
	my ($park_prodiv_file) = $ws_tmp."/parking.prodiv";
	my ($park_memsat_file) = $ws_tmp."/parking.memsat";
	
	local(*PARKPH);
	open (PARKPH, ">>$park_phobius_file");
	local(*PARKPR);
	open (PARKPR, ">>$park_prodiv_file");
	local(*PARKME);
	open (PARKME, ">>$park_memsat_file");

	@fichero = undef;
	$contador = 0;
	@sequence = undef;
	$w = 0;
	$control = 'n';	
	if (-e $phobius_file and (-s $phobius_file > 0) ) {
		local(*PHOB_FILE);
		open (PHOB_FILE, $phobius_file);
		@fichero = <PHOB_FILE>;
		close(PHOB_FILE);
	}
	
	foreach my $z (@fichero) {
		if ( $z =~ /^FT\s+TRANSMEM\s+(\d+)\s+(\d+)/ ) {
			push(@{$phobius2},$1);
			push(@{$phobius2},$2);
			$contador++;
		}
	}		
	$phobius->{$id} = $contador;
	chomp($parking[1]);
	@sequence = split(//,$parking[1]);
	print PARKPH ">$id\t",scalar@sequence ," a.a.\n";
	for (my $i = 1; $i < (scalar@sequence+1); $i++ ) {
		if ( defined $phobius2->[$w] && $i == $phobius2->[$w] && $control eq 'n' ) {
			print PARKPH "X";
			$control='H';
			$w++;
		}
		elsif ( defined $phobius2->[$w] && ($i != $phobius2->[$w]) && ($control eq 'n') ) {
			print PARKPH "-";
		}
		elsif ( defined $phobius2->[$w] && ($i == $phobius2->[$w]) && ($control eq 'H') ){
			print PARKPH "X";
			$control='n';
			$w++;
		}
		elsif ( defined $phobius2->[$w] && ($i != $phobius2->[$w]) && ($control eq 'H') ) {
			print PARKPH "X";
		}
		else {
			print PARKPH "-";
		}
	}
	print PARKPH "\n";

	@fichero = undef;
	$contador = 0;
	@sequence = undef;
	$w = 0;
	$control = 'n';		
	if (-e $prodiv_file and (-s $prodiv_file > 0) ) {		
		local(*PROD_FILE);
		open (PROD_FILE, $prodiv_file);
		@fichero = <PROD_FILE>;
		close(PROD_FILE);
	}
			
	foreach my $z (@fichero) {
		if ( defined $z && ($z =~ /^FT\s+TRANSMEM\t+(\d+)\t+(\d+)/) ){
			push(@{$prodiv2},$1);
			push(@{$prodiv2},$2);
			$contador++;
		}
	}
	$prodiv->{$id}=$contador;
	chomp($parking[1]);
	@sequence = split(//,$parking[1]);
	print PARKPR ">$id\t",scalar@sequence ," a.a.\n";
	for (my $i = 1; $i <(scalar@sequence+1); $i++ ) {
		if ( defined $prodiv2->[$w] && ($i == $prodiv2->[$w]) && ($control eq 'n') ) {
			print PARKPR "X";
			$control='H';
			$w++;
		}
		elsif ( defined $prodiv2->[$w] && ($i != $prodiv2->[$w]) && ($control eq 'n') ) {
			print PARKPR "-";
		}
		elsif ( defined $prodiv2->[$w] && ($i == $prodiv2->[$w]) && ($control eq 'H') ) {
			print PARKPR "X";
			$control='n';
			$w++;
		}
		elsif ( defined $prodiv2->[$w] && ($i != $prodiv2->[$w]) && ($control eq 'H') ) {
			print PARKPR "X";
		}
		else {
			print PARKPR "-";
		}			
	}
	print PARKPR "\n";
	
	@fichero = undef;
	$contador = 0;
	@sequence = undef;
	$w = 0;
	$control = 'n';
	if (-e $memsat_file and (-s $memsat_file > 0) ) {
		local(*MEM_FILE);
		open (MEM_FILE, $memsat_file);
		@fichero = <MEM_FILE>;
		close(MEM_FILE);
	}

	for (my $z=0 ; $z < scalar(@fichero); $z++ ) {
		if ( $fichero[$z] =~ /^FINAL PREDICTION/ ) {
			$z++;
			while ( $z < scalar(@fichero) ) {
				if ( $fichero[$z] =~ /^\d+:\s+\(?([a-z]+)?\)?\s?(\d+)-(\d+)\t\((-?\d+)\.\d+\)/ ) {
					if ($4 > 1.9){
						push(@{$memsat2},$2);
						push(@{$memsat2},$3);
						$contador++;
					}
				}					
				$z++;
			}
		}
	}
	$memsat->{$id}=$contador;
	chomp($parking[1]);
	@sequence = split(//,$parking[1]);
	print PARKME ">$id\t",scalar@sequence ," a.a.\n";
	for (my $i = 1 ; $i <(scalar@sequence+1); $i++ ) {
		if ( defined $memsat2->[$w] && ($i == $memsat2->[$w]) && ($control eq 'n') ) {
			print PARKME "X";
			$control='H';
			$w++;
		}
		elsif ( defined $memsat2->[$w] && ($i != $memsat2->[$w]) && ($control eq 'n') ) {
			print PARKME "-";
		}
		elsif ( defined $memsat2->[$w] && ($i == $memsat2->[$w]) && ($control eq 'H') ) {
			print PARKME "X";
			$control='n';
			$w++;
		}
		elsif ( defined $memsat2->[$w] && ($i != $memsat2->[$w]) && ($control eq 'H') ) {
			print PARKME "X";
		}
		else {
			print PARKME "-";
		}			
	}
	print PARKME "\n";
	
	close (PARKPH);
	close (PARKPR);
	close (PARKME);
	
	return ($park_phobius_file, $park_prodiv_file, $park_memsat_file, $phobius, $prodiv, $memsat, $phobius2, $prodiv2, $memsat2);
}

sub park_consensus2($$$$$$$$$$$)
{
	my ($id, $park_phobius_file, $park_prodiv_file, $park_memsat_file, $phobius, $prodiv, $memsat, $phobius2, $prodiv2, $memsat2, $ws_tmp) = @_;
	
	my ($consensus_file) = $ws_tmp."/consensus.txt";
	my ($no_predict_file) = $ws_tmp."/no_predict.txt";
	
	local(*RESUL);
	open (RESUL, ">>$consensus_file");
	local(*RESUL2);
	open (RESUL2, ">>$no_predict_file");
	
	if ( exists $phobius->{$id} && exists $memsat->{$id} && exists $prodiv->{$id} ) {
		my (@phobius3) = ();
		my (@prodiv3) = ();
		my (@memsat3) = ();
		my ($longitud);		

		local(*PHOB_FILE);
		open (PHOB_FILE, $park_phobius_file);
		my (@phob_cont) = <PHOB_FILE>;
		close(PHOB_FILE);
		for ( my $z = 0; $z < scalar(@phob_cont); $z++ ) {
			if ( $phob_cont[$z] =~ />$id\t(\d+) a.a./ ) {
				$longitud=$1;
				chomp($phob_cont[$z+1]);
				@phobius3=split(//,$phob_cont[$z+1]);
			}
		}
		
		local(*PROD_FILE);
		open (PROD_FILE, $park_prodiv_file);
		my (@prod_cont) = <PROD_FILE>;
		close(PROD_FILE);
		for (my $z = 0; $z < scalar(@prod_cont); $z++ ) {
			if ( $prod_cont[$z] =~ />$id\t(\d+) a.a./ ) {
				chomp($prod_cont[$z+1]);
				@prodiv3=split(//,$prod_cont[$z+1]);
			}
		}
		
		local(*MEM_FILE);
		open (MEM_FILE, $park_memsat_file);
		my (@mem_cont) = <MEM_FILE>;
		close(MEM_FILE);
		for (my $z = 0; $z < scalar(@mem_cont); $z++ ) {
			if ( $mem_cont[$z] =~ />$id\t(\d+) a.a./ ) {
				chomp($mem_cont[$z+1]);
				@memsat3=split(//,$mem_cont[$z+1]);
			}
		}
			
		print RESUL ">$id\tlength $longitud a.a.\n";
		for (my $z = 0; $z < $longitud; $z++ ) {
			if ( defined $phobius3[$z] && defined $prodiv3[$z] && defined $memsat3[$z] ) {				
				if ( ($phobius3[$z] eq $prodiv3[$z]) && ($phobius3[$z] eq $memsat3[$z]) && ($phobius3[$z] eq 'X') ) {
					print RESUL "X";
				}
				else {
					print RESUL "-";
				}				
			}
			else {
				print RESUL "-";
			}
		}
		print RESUL "\n";
	}
	else {
		print RESUL2 "$id\tdoesn't have all the predictions.\n";
	}

	close (RESUL);
	close (RESUL2);
	
	return $consensus_file;
}

#filtrado por numero minimo de  9 a.a. en la helice
sub filter_length($$$)
{
	my ($tmp_consensus_file, $ws_cache, $ws_tmp) = @_;
	
	local(*RESUL);
	open (RESUL, $tmp_consensus_file);
	my (@content) = <RESUL>;
	close(RESUL);
	
	my ($longitud) = 0;
	my ($presencia);
	
	my ($num_helices_file) = $ws_tmp."/num_helices.txt";	
	my ($consensus_file) = $ws_tmp."/consensus.txt";
	my ($thump_file) = $ws_tmp."/thump.txt";
	
	open (RE,">$consensus_file");
	open (COORD,">$thump_file");
	
	for ( my $i = 0; $i < scalar(@content); $i++ ) {
		if ( $content[$i] =~ /^>.+/ ) {
			print RE $content[$i];
			print COORD $content[$i];
		}
		else {
			chomp($content[$i]);
			my (@parking) = ();
			my (@sequence) = ();
			my (@listado) = split(//,$content[$i]);
			for (my $w = 0; $w <= scalar(@listado); $w++ ) {
				if ( $w > 0 && defined $listado[$w] && defined $listado[$w-1] && ($listado[$w] eq '-') && ($listado[$w-1] eq '-') ) {
					push(@sequence,'-');
				}
				if ( $w == 0 && defined $listado[$w] && $listado[$w] eq '-' ) {
					push(@sequence,'-');
				}
				elsif ( defined $listado[$w] && defined $listado[$w-1] && ($listado[$w] eq '-') && ($listado[$w-1] eq 'X') ) {
					if ( exists $presencia->{scalar(@parking)} ) {
						$presencia->{scalar(@parking)}++;
					}
					else {
						$presencia->{scalar(@parking)}=1;
					}
					if ( scalar(@parking) > 9 ) {
						foreach my $z (@parking){
							push(@sequence,'X');
						}
					}
					else {
						foreach my $z (@parking){
							push(@sequence,'-');
						}
					}
					push(@sequence,'-');
					@parking = ();
				}  
				elsif ( defined $listado[$w] && ($listado[$w] eq 'X') ) {
					push(@parking,'X');
				}
			}
			my ($k) = 0;
			my ($num_helix) = 0;
			foreach my $z (@sequence) {
				$k++;
				print RE $z;
				if ( defined $sequence[$k-2] && ($z eq 'X') && ($sequence[$k-2] eq '-') ) {
					$num_helix++;
					$longitud=$k;
					print COORD "helix number $num_helix start: $k\t";
				}
				elsif ( defined $sequence[$k-2] && ($z eq '-') && ($sequence[$k-2] eq 'X') ) {
					if ( ($k - $longitud) > 14) {
						print COORD "end: ",$k-1,"\n";
					}
					else {
						print COORD "end: ",$k-1,"\tdamaged\n";
					}
				}						
			}
			print RE "\n";
		}
	}
	close(RE);
	close(COORD);
	
	local(*NUM);
	open (NUM,">$num_helices_file");	
	my (@key) = keys(%{$presencia});
	foreach my $z (@key){
		print NUM "$z\t",$presencia->{$z},"\n";
	}
	close(NUM);
	
	return ($num_helices_file, $consensus_file, $thump_file);
}

# This program localizes and opens the files THUMP.txt from all the analysis directories. 
# Then evaluates the output and tries to extrapolate information (when it exists) about the 
# principal isoform. There are a some different cases and the information about the damaged "helices" is considered too.
# The actual evaluationn schema is:
# if there's a unique protein that has the highest TMH, this will be tagged as "YES";
# if there are more than one protein that have the same highest TMH number for their parental gene, these are tagged as "UNKNOWN";
#	*** How we calculate the highest TMH number ?? Taking in account damaged helices too.
#	-> e.g. If a gene has 2 protein products, both with 3 helices, and 1 of the TMH of a transcript is tagged as "damaged",
#		both of them will be tagged as UNKNOWN, because if the damaged helix will be confirmed as good-helix, 
#		they could have the same max. TMH number
#	-> e.g. Imagine now that one gene has 2 protein products. The first has 5 helix, one tagged as damaged. The second one has
#		4 helices, none of them damaged. These proteins will be tagged as UNKNOWN, because if if the damaged helix will be not
#		confirmed as good-helix, they could have the same max. TMH number.
# if a protein couldn't have the highest TMH, even considering damaged helices too, this will be tagged as "NO".
#
# So, in general, if we have in a gene a protein product tagged as YES, all the others will be tagged as NO.
# We can have in the same gene more than one product tagged as UNKNOWN or NO.
# YES and UNKNOWN can't coexist in the same gene's products annotation.
sub _get_appris_annotations($)
{
	my ($aux_rst_cont) = @_;
	
	my ($cutoffs);
	my ($output_content) = '';
		$output_content .= "\n";
		$output_content .= "# ==================================== #\n";
		$output_content .= "# Prediction of trans-membrane helices #\n";
		$output_content .= "# ==================================== #\n";	
	
	my ($transcripts);
	my ($trans_num_helices);
	my ($trans_num_dam_helices);
	
	my (@trans_rsts) = split('>', $aux_rst_cont);
	foreach my $trans_rst (@trans_rsts) {
		if ( $trans_rst =~/^([^\t]*)\tlength \d+ a\.a\./ ) {
			my ($trans_id) = $1;
			my ($num_damaged_helices) = 0;
			my ($num_helices) = 0;
			while ( $trans_rst =~ /^(helix number.+[^\n]*)/mg ) {
				$num_helices++;
				$num_damaged_helices++ if ( $1 =~ /damaged/ );
			}
			$transcripts->{$trans_id} = [$num_helices, $num_damaged_helices];
			push(@{$trans_num_helices->{$num_helices}}, $trans_id);
			push(@{$trans_num_dam_helices->{$num_damaged_helices}}, $trans_id); 
		}		
	}
	
	if ( defined $trans_num_helices ) {
		my (@trans_num_helices_list) = sort { $a <= $b } keys(%{$trans_num_helices} );
		
		if ( scalar(@trans_num_helices_list) > 0 )
		{
			# We tag the transcript as UNKOWN whose num helices are biggest
			my (@trans_biggest_num_helices);
			my ($biggest_num_helices) = $trans_num_helices_list[0];
			my ($biggest_num_helices2) = $trans_num_helices_list[0];
			foreach my $trans_id (@{$trans_num_helices->{$biggest_num_helices}})
			{
				push(@trans_biggest_num_helices, $trans_id);
			}
			if ( scalar(@trans_biggest_num_helices) == 1 ) {
				my ($trans_id2) = $trans_biggest_num_helices[0];
				my ($b_n_helices) = $transcripts->{$trans_id2}[0];
				my ($b_n_d_helices) = $transcripts->{$trans_id2}[1];
				$biggest_num_helices2 = $biggest_num_helices - $b_n_d_helices;		
			}
			
			my ($unique) = 1;
			for (my $i = 1; $i < scalar(@trans_num_helices_list); $i++)
			{
				my ($current_score) = $trans_num_helices_list[$i];
				
				# If the biggest score is bigger than current => the transcripts are rejected
				if ( ($biggest_num_helices2 - $current_score) > 0 )
				{
					foreach my $trans_id (@{$trans_num_helices->{$current_score}})
					{
						$cutoffs->{$trans_id} = $NO_LABEL;
					}
				}
				else
				{
					$unique=0;
					foreach my $trans_id (@{$trans_num_helices->{$current_score}})
					{
						$cutoffs->{$trans_id} = $UNKNOWN_LABEL;
					}					
				}
			}
			# There is one transcript with the biggest score
			if (($unique == 1) and scalar(@trans_biggest_num_helices) == 1)
			{
				foreach my $trans_id (@trans_biggest_num_helices)
				{
					$cutoffs->{$trans_id} = $OK_LABEL;				
				}
			}
			else
			{
				foreach my $trans_id (@trans_biggest_num_helices)
				{
					$cutoffs->{$trans_id} = $UNKNOWN_LABEL;				
				}
			}
		}
		
	}
	
	# Get appris output
	if ( defined $cutoffs ) {
		while ( my ($trans_id,$annot) = each (%{$cutoffs}) ) {
			$output_content .= ">".$trans_id."\t".$annot."\n";
		}		
	}
	return $output_content;	
}
sub _delete_tmp_dir($)
{
	my ($dir) = @_;
	eval {
		my ($cmd) = "rm -rf $dir/*";
		$logger->debug("\n** script: $cmd\n");
		my (@out) = `$cmd`;
	};
	$logger->error("_delete_tmp_dir") if($@);
	return undef;
}

main();

__END__

=head1 NAME

thump

=head1 DESCRIPTION

This is a consensus program that combines the predictions of 3 algorithms:

	Memsat3 [ref], Phobius [ref]  and Prodiv [ref]
	
for predict alpha transmembrane helices in proteins using a strict filter schema, 
in order to minimize the False Positive Ratio.

=head1 ARGUMENTS

=head2 Required arguments:

	--conf <Config file>

    --input <Fasta sequence file>

    --output <Annotation output file>

=head2 Optional arguments:

	--appris <Flag that enables the output for APPRIS (default: NONE)>

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
    

=head1 EXAMPLE

perl thump.pl

	--conf=../conf/pipeline.ini
	
	--input=examples/ENSG00000154639.faa
	
	--output=examples/ENSG00000154639.output


=head1 AUTHOR

Created by

	Paolo Maietta -pmaietta@cnio.es-

Updated and mainteined by

	Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut