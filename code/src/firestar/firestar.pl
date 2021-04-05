#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin;
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
	$DEFAULT_CONFIG_FILE
	$DEFAULT_FIRESTAR_CONFIG_FILE
	$APPRIS_HOME
	$WSPACE_TMP
	$WSPACE_CACHE
	$CACHE_FLAG
	$NAME_DIR
	$PROG_EVALUE
	$PROG_CUTOFF
	$PROG_CSA
	$PROG_COG
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
my ($cfg) 				= new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD				= $FindBin::Bin;
$DEFAULT_FIRESTAR_CONFIG_FILE	= $ENV{APPRIS_CODE_CONF_DIR}.'/firestar.ini';
$APPRIS_HOME			= $ENV{APPRIS_HOME};
$WSPACE_TMP				= $ENV{APPRIS_TMP_DIR};
$WSPACE_CACHE			= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
$CACHE_FLAG				= $cfg->val('FIRESTAR_VARS', 'cache');
$NAME_DIR				= $cfg->val('FIRESTAR_VARS', 'name');
$PROG_EVALUE			= $cfg->val('FIRESTAR_VARS', 'evalue');
$PROG_CUTOFF			= $cfg->val('FIRESTAR_VARS', 'cutoff');
$PROG_CSA				= $cfg->val('FIRESTAR_VARS', 'csa');
$PROG_COG				= $cfg->val('FIRESTAR_VARS', 'cog');

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);

#####################
# Method prototypes #
#####################
sub exist_motif($$);

#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	# Declare vars
	my ($output_content) = '';
	my (%gene_vars);
	my (%motifs) = ();
	my (%res_count);
	my (%pre_count);	
	my (%motif_count);
	my (%motif_resnum);
	my (%motif_freq);
	my (%motif_score);
	my (%motif_rel);
	my (%motif_sort);
	my (@motif_sort);
	my (%motif_ids);
	my (%var_sumas);
	my (%order) = ();
	my (@winners);
	my (@loosers);
		
	# Declare and init the local variables
	$logger->info("-- declare and init the local variables\n");
		
	$output_content .= "/*\n";
	$output_content .= " * firestar\n";
	$output_content .= " *   prediction of functionally important residues using structural templates and alignment reliability.\n";
	$output_content .= " * Gonzalo Lopez; A. Valencia; M. Tress Nucleic Acids Research, doi:10.1093/nar/gkm297\n";
	$output_content .= " * Date: Jan 1, 2011\n";
	$output_content .= " */\n";

	
	# Get sequences
	$logger->info("-- get the sequences of variants\n");
    my $fasta_object = Bio::SeqIO->new(
                        -file => $input_file,
                        -format => 'Fasta'
    );
	while ( my $seq = $fasta_object->next_seq() )
	{
		if ( $seq->id=~/^([^|]*)\|([^|]*)/ )
		{			
			my ($sequence_id) = $2;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
            my ($sequence) = $seq->seq;
            $gene_vars{$sequence_id} = $sequence;
		}
	}	
	#$logger->debug("-- variant sequences:\n".Dumper(%gene_vars)."\n---------------\n");
	
	# Run firestar_pipeline for every sequence. ----------------
	# For every sequence save within data structure
	$logger->info("-- run firestar pipeline for every variant\n");	
	
	$output_content .= "\n";
	$output_content .= "# ============================================ #\n";
	$output_content .= "# Functional residues and structural templates #\n";
	$output_content .= "# ============================================ #\n";
	foreach my $varname (keys %gene_vars)
	{
		$logger->info("\t-- $varname ");
		my $seq = $gene_vars{$varname};
				
		# create cache obj
		my ($cache) = APPRIS::Utils::CacheMD5->new(
			-dat => $seq,
			-ws  => $WSPACE_CACHE			
		);		
		my ($seq_idx) = $cache->idx;
		my ($seq_sidx) = $cache->sidx;
		my ($seq_dir) = $cache->dir;
				
		# prepare cache dir
		my ($ws_cache) = $cache->c_dir();
		# prepare tmp dir
		my ($ws_tmp) = $WSPACE_TMP.'/'.$seq_idx;
		prepare_workspace($ws_tmp);		
		
		# Setup configure file of firestar
		my ($firestar_config_cont) = getStringFromFile($DEFAULT_FIRESTAR_CONFIG_FILE);
		my $subs_template = sub {
			my ($cont, $old, $new) = @_;	
			$cont =~ s/$old/$new/g;		
			return $cont;		
		};
-
		$firestar_config_cont = $subs_template->($firestar_config_cont, 'APPRIS__CACHE__WORKSPACE', $ws_cache);
		$firestar_config_cont = $subs_template->($firestar_config_cont, 'APPRIS__HOME', $APPRIS_HOME);
		my ($firestar_config_file) = $ws_tmp.'/firestar.ini';		
		my ($firePredText_log) = $ws_tmp.'/'.'firestar.log';

		my ($firestar_config) = new Config::IniFiles( -file => \$firestar_config_cont );
		my ($PROG_DB_UID) = $firestar_config->val('DATABASES', 'release');
		my ($firePredText_file_path) = $ws_cache.'/'."seq.firestar_$PROG_DB_UID";

		my ($print_log) = printStringIntoFile($firestar_config_cont, $firestar_config_file);
		$logger->error("-- printing firestar config annot") unless ( defined $print_log );

		# If output is not cached
		unless ( -e $firePredText_file_path and (-s $firePredText_file_path > 0) and ($CACHE_FLAG eq 'yes') ) {
			my ($cmd) = "perl $LOCAL_PWD/source/perl/firestar.pl -opt appris -q seq -e $PROG_EVALUE -cut $PROG_CUTOFF -csa $PROG_CSA -cog $PROG_COG -s $seq -o $firePredText_file_path -conf $firestar_config_file 2> $firePredText_log";
			$logger->debug("\n** script: $cmd\n");			
			my (@firePredText_out) = `$cmd`;
			$logger->error("Empty output of firestar: $varname") unless ( -e $firePredText_file_path && -s $firePredText_file_path > 0 );
		}
			
		my ($firePredText_cont) = getStringFromFile($firePredText_file_path);
		$firePredText_cont =~ s/^>>>Query/>>>$varname/mg; # changet the 'Query' identifier from firestar to the current identifier 
		my (@firePredText_cont2) = split(/\n/,$firePredText_cont);
		foreach my $line (@firePredText_cont2)
		{
			$output_content .= $line."\n";
			if ( $line =~ /^\d/ )
			{
				my @spl = split(/\t/,$line);
				my $resnum = $spl[0];
				my $key = $spl[1];
				my $motif = $spl[1];
				my $motif_ids = $spl[2];
				my $motif_list_all = $spl[2];
				# add motif reliability, square, and frequency
				my @motif_list_all2 = split(/\|/,$motif_list_all);
				foreach my $motif_list (@motif_list_all2) {
					$motif_list =~ s/\s*$//; $motif_list =~ s/^\s*//;
					if ( $motif_list =~ /^([^\[]*)\[([^\,]*)\,([^\,]*)\,([^\]]*)\]$/ ) {
						my $ligand = $1;
						my $freq_score = $2;
						my $square_score = $3;
						my $rel_score = $4;
						# add motif frequency
						if ( exists $motif_freq{$key}{'main'} ) {
							$motif_freq{$key}{'main'} = $freq_score if ( $motif_freq{$key}{'main'} < $freq_score );
						}
						else {
							$motif_freq{$key}{'main'} = $freq_score;
						}
						push(@{$motif_freq{$key}{'list'}}, $freq_score);
						# add motif score
						if ( exists $motif_score{$key}{'main'} ) {
							$motif_score{$key}{'main'} = $square_score if ( $motif_score{$key}{'main'} < $square_score );
						}
						else {
							$motif_score{$key}{'main'} = $square_score;
						}
						push(@{$motif_score{$key}{'list'}}, $square_score);
						# add motif reliability
						if ( exists $motif_rel{$key}{'main'} ) {
							$motif_rel{$key}{'main'} = $rel_score if ( $motif_rel{$key}{'main'} < $rel_score );
						}
						else {
							$motif_rel{$key}{'main'} = $rel_score;
						}
						push(@{$motif_rel{$key}{'list'}}, $rel_score);
					}
					
				}
				
				# add motif info
				push(@{$motifs{$key}}, $varname);
				unless ( exists $motif_sort{$key} ) 
				{
					$motif_sort{$motif} = "";
					push(@motif_sort, $key);
				}
				$motif_count{$key} = 0;
				$motif_ids{$varname}{$key} = $motif_ids;
				
				# Important: take into account that repeated motifs can exit.
				# 	For that reason, we save the list of residues for one motif
				push(@{$motif_resnum{$varname}{$key}}, $resnum);
			}
		}
	}

	$logger->debug("-- motifs\n".Dumper(%motifs)."\n---------------\n");
	$logger->debug("-- motif residue positions\n".Dumper(%motif_resnum)."\n---------------\n");
	$logger->debug("-- motif freq scores\n".Dumper(%motif_freq)."\n---------------\n");
	$logger->debug("-- motif square scores\n".Dumper(%motif_score)."\n---------------\n");
	$logger->debug("-- motif reliability scores\n".Dumper(%motif_rel)."\n---------------\n");
	$logger->debug("-- motif id list\n".Dumper(%motif_ids)."\n---------------\n");

	# Count the motifs ----------------
	foreach my $varname (keys %gene_vars)
	{	
		# BEGIN: DEPRECATED
		#$res_count{$varname}{'seis'} = 0;
		#$res_count{$varname}{'cinco'} = 0;
		#$res_count{$varname}{'cuatro'} = 0;
		#$res_count{$varname}{'tres'} = 0;
		# END: DEPRECATED
		$res_count{$varname}{'total'} = 0;

		if ( exists $motif_resnum{$varname} and defined $motif_resnum{$varname} )
		{
			foreach my $key (keys %{$motif_resnum{$varname}})
			{			
				if($key =~ /\w/)
				{
					# one motif can appear several times
					my ($count) = 1;
					if (scalar(@{$motif_resnum{$varname}{$key}}) > 1) {
						$count = scalar(@{$motif_resnum{$varname}{$key}});
					}					
					
					$motif_count{$key} += $count;
					$res_count{$varname}{'total'} += $count;
					
					# BEGIN: DEPRECATED
					#if ( $motif_score{$key}{'main'} == 6 )
					#{
					#	$res_count{$varname}{'seis'} += $count;
					#}
					#elsif ( $motif_score{$key}{'main'} == 5 )
					#{
					#	$res_count{$varname}{'cinco'} += $count;
					#}
					#elsif ( $motif_score{$key}{'main'} == 4 )
					#{
					#	$res_count{$varname}{'cuatro'} += $count;
					#}											
					#elsif ( $motif_score{$key}{'main'} == 3 )
					#{
					#	$res_count{$varname}{'tres'} += $count;
					#}
					# END: DEPRECATED
				}
			}			
			
		} 
	}
	$logger->debug("-- motif count\n".Dumper(%motif_count)."\n---------------\n");	
	$logger->debug("-- residues count\n".Dumper(%res_count)."\n---------------\n");

	
	# Check the repeated motifs ----------------
	$logger->debug("-- check the repeated motifs\n");
	my ($repeated_motifs);
	foreach my $varname (keys %gene_vars)
	{
		$logger->debug(">>$varname:\n");
		my $seq = $gene_vars{$varname};
		foreach my $key (keys %motifs)
		{
			$logger->debug("Motif: $key: ");
			# discard the motifs whose orignal length is less than 13			
			if ( ($key =~ /\w/) and (length $key == 13) ) 
			{
				# get the 10 residues beginning from left and rigth
				my ($motif) = $key;
				my ($motif_l) = $key;
				my ($motif_r) = $key;
				$motif_l =~ s/\w{3}$//g;
				$motif_r =~ s/^\w{3}//g;

				# scan for each
				my @new_motifs = ($motif, $motif_l, $motif_r);
				my ($index) = 0;
				my ($find_motif) = 0;
				while ( $index < scalar(@new_motifs) and ($find_motif == 0) )
				{
					my ($new_motif) = $new_motifs[$index]; 
					$logger->debug("\nNew_Motif:$index: $new_motif");

					# check if smaller motif match agains other motiff
					while ( $seq =~ m/($new_motif)/g )
					{
						$logger->debug("\nNew_Motif: $new_motif");

						# see if new motif does not exits within current var
						if ( !(exist_motif($motif_resnum{$varname}, $new_motif)) )
						{
							$logger->debug("\tnot_exist_new_motif:");

							if (exists $motifs{$key} and defined $motifs{$key} and (scalar(@{$motifs{$key}}) > 0) and 
								exists $motif_freq{$key} and defined $motif_freq{$key} and
								exists $motif_score{$key} and defined $motif_score{$key} and exists $motif_score{$key}{'main'}
							) {
								$logger->debug("ok1:");

								# find the correct position of center residue of current var
								my ($plus) = 7;
								if ($index == 1) { $plus = 7 } # left side
								elsif ($index == 2) { $plus = 4 } # right side
								my ($repeated_residue) = rindex($seq, $new_motif) + $plus;

								# add residue if not already exist
								unless (exists $repeated_motifs->{$varname}->{$repeated_residue})
								{
									$find_motif = 1;
									
									my ($repeated_motif) = $key;
									if ($index == 1) { $repeated_motif =~ s/\w{3}$/---/g; } # left side
									elsif ($index == 2) { $repeated_motif =~ s/^\w{3}/---/g; } # right side
									
									my ($repeated_varname) = $motifs{$key}->[0];
									my ($repeated_res_freq) = $motif_freq{$key}{'main'};
									my ($repeated_res_score) = $motif_score{$key}{'main'};
									my ($repeated_res_rel) = $motif_rel{$key}{'main'};
									my ($repeated_motif_ids) = $motif_ids{$repeated_varname}{$key};
									$repeated_motif_ids =~ s/\n*//g;
									$repeated_motifs->{$varname}->{$repeated_residue} = {
																					'motif'	=> $repeated_motif,
																					'freq'	=> $repeated_res_freq,
																					'score'	=> $repeated_res_score,
																					'rel'	=> $repeated_res_rel,
																					'ids'	=> $repeated_motif_ids,
									};
									# BEGIN: DEPRECATED
									#if ( $motif_score{$key}{'main'} == 6 )
									#{
									#	$res_count{$varname}{'seis'} += 1;
									#}
									#elsif ( $motif_score{$key}{'main'} == 5 )
									#{
									#	$res_count{$varname}{'cinco'} += 1;
									#}
									#elsif ( $motif_score{$key}{'main'} == 4 )
									#{
									#	$res_count{$varname}{'cuatro'} += 1;
									#}
									#elsif ( $motif_score{$key}{'main'} == 3 )
									#{
									#	$res_count{$varname}{'tres'} += 1;
									#}
									# END: DEPRECATED
									# count new motifs
									$motif_count{$key} += 1;
									$res_count{$varname}{'total'} += 1;
									
									$logger->debug("ok2\n");
									$logger->debug(Dumper($repeated_motifs->{$varname}->{$repeated_residue}));
								}
							}
						}
					}
					$index += 1;
					
					$logger->debug("\n");
				}
			}
		}
	}
	$logger->debug("-- motif count\n".Dumper(%motif_count)."\n---------------\n");
	$logger->debug("-- residues count 2\n".Dumper(%res_count)."\n---------------\n");
	$logger->debug("-- repeated motifs\n".Dumper($repeated_motifs)."\n---------------\n");
		

	
	
	# Select rejected and accepted sequences ----------------
	#if ( defined $appris )
	#{
		# get consensus residues
		$logger->info("-- print consensus residues\n");

		$output_content .= "\n";
		$output_content .= "# ============================== #\n";
		$output_content .= "# Prediction of consensus motifs #\n";
		$output_content .= "# ============================== #\n";
		foreach my $varname (keys %gene_vars)
		{
			$output_content .= "######\n";

			my ($num_residues) = '';
			my ($pos_list_residues) = '';

			if ( exists $repeated_motifs->{$varname} and defined $repeated_motifs->{$varname} )
			{
				my ($repeated_motifs_var) = $repeated_motifs->{$varname};
				my (@sort_residues) = sort {$a <=> $b} keys %{$repeated_motifs_var};

				foreach my $sort_residue (@sort_residues)
				{
					my ($repeated_res_score) = $repeated_motifs_var->{$sort_residue}->{'score'};
					my ($repeated_motifs) = $repeated_motifs_var->{$sort_residue}->{'motif'};
					my ($repeated_motif_ids) = $repeated_motifs_var->{$sort_residue}->{'ids'};
					$output_content .= "$sort_residue\t$repeated_motifs\t$repeated_motif_ids\n";
					$pos_list_residues .= "$sort_residue,";
				}

				$pos_list_residues =~ s/\,$//g;
				$num_residues = scalar(@sort_residues) if ( scalar(@sort_residues) > 0 );
			}

			$output_content .= "C>>\t$varname";
			$output_content .= "\t".$num_residues if ($num_residues ne '');
			$output_content .= "\t".$pos_list_residues if ($pos_list_residues ne '');
			$output_content .= "\n";
		}
	
		$logger->info("-- select which is the variant with more functional residues /* APPRIS */\n");	

		# get the score of consensus motifs		
		foreach my $varname (keys %res_count)
		{
			# BEGIN: DEPRECATED
			#my $sum = (6*$res_count{$varname}{'seis'} + 5*$res_count{$varname}{'cinco'} + 4*$res_count{$varname}{'cuatro'} + 0*$res_count{$varname}{'tres'});
			# END: DEPRECATED
			my $sum = $res_count{$varname}{'total'};
			$var_sumas{$varname} = $sum;
			if ( exists $order{$sum} )
			{	
				$order{$sum} = "$order{$sum} $varname";
			}
			else
			{
				$order{$sum} = $varname;
			}
		}
		#$logger->debug("-- motif scores\n".Dumper(%order)."\n---------------\n");
		
		# sort the order of motifs
		my @sort = sort {$b <=> $a} keys %order;
		my $win = shift @sort;
		push(@winners, $order{$win});
		
		# get the number of residues of winner sequence(s)
		my $num_residues_win = 0;
		foreach my $line (@winners) {
			foreach my $varname (split(/ /, $line)) {
				if ( $res_count{$varname}{'total'} > $num_residues_win )
					{ $num_residues_win = $res_count{$varname}{'total'}; }
			}
		}
		
#		# get the list of winners and loosers 
#		for ( my $i = 0; $i < $#sort+1; $i++ )
#		{
#			my $nextsum = $sort[$i];
#			
#			# We believe in firestar but when the number of functional residues are bigger than 2
#			if ( $num_residues_win <= $PROG_MIN_RESIDUES ) 
#			{		
#				my ($list_of_vars) = '';
#				my @vars = split(/ /, $order{$nextsum});
#				push(@winners, $order{$nextsum});				
#			}
#			else # we believe
#			{
#				if ( $win - $nextsum <= $PROG_DIFF_RESIDUES )
#				{
#					push(@winners, $order{$nextsum});
#				}
#				else
#				{
#					push(@loosers, $order{$nextsum});
#				}		
#			}
#		}
		# get the list of winners and loosers 
		for ( my $i = 0; $i < $#sort+1; $i++ )
		{
			my $nextsum = $sort[$i];
			push(@winners, $order{$nextsum});				
		}

		
		# print rejected and accepted sequences ----------------
		$output_content .= "\n";
		$output_content .= "# ========================================== #\n";
		$output_content .= "# Final annotations: no. functional residues #\n";
		$output_content .= "# ========================================== #\n";
		foreach my $line (@winners)
		{
			my @vars = split(/ /, $line);
			foreach my $varname (@vars)
			{
				#$output_content .= "ACCEPT: $varname\t$var_sumas{$varname}\t$res_count{$varname}{total}\n";
				$output_content .= "F>>\t$varname\t$var_sumas{$varname}\t$res_count{$varname}{total}\n";
			}
		}
		
		foreach my $line (@loosers)
		{
			my @vars = split(/ /, $line);
			foreach my $varname (@vars)
			{
				#$output_content .= "REJECT: $varname\t$var_sumas{$varname}\t$res_count{$varname}{total}\n";
				$output_content .= "F>>\t$varname\t$var_sumas{$varname}\t$res_count{$varname}{total}\n";
			}
		}	
	#}	

	# Print output ----------------	
	my ($print_out) = printStringIntoFile($output_content, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}
	
	$logger->finish_log();
		
	exit 0;	
}

sub exist_motif($$) 
{
	my ($var_motifs, $motif) = @_;
	
	my ($exits) = 0;
	foreach my $var_motif (keys(%{$var_motifs}))
	{
		if ( $var_motif =~ m/$motif/g )
		{
			$exits = $var_motif;
			return $exits;
		}
	}	
	return $exits;
}

main();

__END__

=head1 NAME

firestar

=head1 DESCRIPTION

Run the main script of firestar pipeline

=head1 SYNOPSIS

firestar

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

perl firestar.pl
	
	--conf=../conf/pipeline.ini
	
	--input=examples/ENSG00000198563.faa
	
	--output=examples/ENSG00000198563.output


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
