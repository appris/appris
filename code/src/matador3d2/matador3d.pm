#ÊBEGIN: CHANGE for APPRIS
##!/usr/bin/env perl
package matador3d;
#ÊEND: CHANGE for APPRIS

#ÊBEGIN: CHANGE for APPRIS
#use warnings;
#ÊEND: CHANGE for APPRIS
use strict;
use File::Basename;
use Getopt::Long;


my $input_fasta;
my $hmmscan_bin;
my $database;
my $tmp_dir = "/tmp";

#ÊBEGIN: CHANGE for APPRIS
### Test
##$input_fasta = "gencode.v23.pc_translations.fa";
#
#GetOptions ("i=s" 	=> \$input_fasta,
#			"t=s" 	=> \$tmp_dir,
#			"bin=s" => \$hmmscan_bin,
#			"db=s" 	=> \$database,
#			);
#
## Check mandatory parameters
#if(!$input_fasta || !-e $input_fasta){
#	print STDERR "ERROR: The input fasta file has not been specified or it does not exist\n";
#	usage();
#}
#if(!$hmmscan_bin || !-x $hmmscan_bin){
#	print STDERR "ERROR: The path to hmmscan binary has not been specified or cannot be executed\n";
#	usage();
#}
#if(!$database || !-e $database){
#	print STDERR "ERROR: The path to the HMM database has not been specified or it does not exist\n";
#	usage();
#}
#
## Load input file, $sequences is a hash, key1 is gene_id, key 2 is transcript_id, the value is the sequence of the correspondng transcript
#my $sequences = read_fasta_file($input_fasta);
#
## Run each hmmscan for each sequence and store results in temporary directory
#run_hmmscan($hmmscan_bin, $sequences, $tmp_dir, $database);
#
## Score isoforms
#my $scoring = run_scoring($sequences, $tmp_dir);
#
#
#exit;
#
#ÊEND: CHANGE for APPRIS


sub run_hmmscan{
	
	my $hmmscan_bin = shift;
	my $sequences 	= shift;
	my $tmp_dir 	= shift;
	my $database 	= shift;
	
	foreach my $gene_id (sort keys %{$sequences}){
		foreach my $transcript_id (sort keys %{$sequences->{$gene_id}}){
			
#			print "$transcript_id\n";
			
			#ÊBEGIN: CHANGE for APPRIS
			# TODO: Cambiar a directorio en appris
			# Comentario: Probablemente seria bastante mas rapido hacer una busqueda global con todas las secuencias y luego separar los ficheros
			# en lugar de hacer una busqueda por secuencia.
#			my $output_dir 				= "$tmp_dir/$transcript_id";
#			my $sequence_file 			= "$output_dir/seq.fasta";
#			my $output_ali_file 		= "$output_dir/matador3D.alis";
#			my $output_dom_file 		= "$output_dir/matador3D.domtblout";
#			my $output_dom_file_sorted 	= "$output_dir/matador3D.domtblout.sorted";
#
#			if(!-d "$output_dir"){
#				`mkdir -p $output_dir`;
#			}
			
			# get sequence
			my $sequence = $sequences->{$gene_id}{$transcript_id};
			
			# Create cache obj
			my ($cache) = APPRIS::Utils::CacheMD5->new(
				-dat => $sequence,
				-ws  => $tmp_dir			
			);		
			my ($seq_idx) = $cache->idx;
			my ($seq_sidx) = $cache->sidx;
			my ($seq_dir) = $cache->dir;
								
			# prepare cache dir
			my ($ws_cache) = $cache->c_dir();
			
			# var files
			my $sequence_file 			= "$ws_cache/seq.faa";
			#my $output_ali_file 		= "$ws_cache/seq.pdb70.hmm";
			my $output_ali_file 		= "$ws_cache/seq.".$database;
			my $output_dom_file 		= "$ws_cache/seq.matador3d.domtblout";
			my $output_dom_file_sorted 	= "$ws_cache/seq.matador3d";
			#ÊEND: CHANGE for APPRIS
			
			# Create input sequence file
			#ÊBEGIN: CHANGE for APPRIS
#			if(!-e "$output_dom_file_sorted"){
#				open(FH, ">", $sequence_file) || die "Error while opening file $sequence_file";
#				my $sequence = $sequences->{$gene_id}{$transcript_id};
#				print FH ">$transcript_id\n$sequence\n";
#				close FH;
			unless(-e $output_ali_file and (-s $output_ali_file > 0) and -e $output_dom_file_sorted and (-s $output_dom_file_sorted > 0) and ($main::CACHE_FLAG eq 'yes')) # Cached Blast
			{
				if(!-e "$sequence_file"){
					open(FH, ">", $sequence_file) || die "Error while opening file $sequence_file";
					my $sequence = $sequences->{$gene_id}{$transcript_id};
					print FH ">$transcript_id\n$sequence\n";
					close FH;
				}
				#ÊEND: CHANGE for APPRIS
				
				# Run hmmscan for this sequence against HMM database
				`$hmmscan_bin -o $output_ali_file --notextw --domtblout $output_dom_file $database $sequence_file`;
				`grep -v "^#" $output_dom_file | sort -grk14 > $output_dom_file_sorted`; # Remove superflous lines starting with # and sort by decreasing domain bitscore
				# Un checkeo de todo ha funcionado no estaria demas, habria que desactivar -o /dev/null
				
				# Delete temp files
				`rm $sequence_file $output_dom_file`;
			}
		}
	}
}


sub run_scoring{

	my $sequences 	= shift;
	my $tmp_dir 	= shift;
	
	#ÊBEGIN: CHANGE for APPRIS
	my $output = '';
	#ÊEND: CHANGE for APPRIS
	
	foreach my $gene_id (sort keys %$sequences){
	
		# Hashes where the data will be stored
		my %gene;
		my %domains;
		my %domains_sortable_list;
		
		# Load data from hmmscan output files
		load_data($gene_id, \%gene, \%domains, \%domains_sortable_list, $sequences, $tmp_dir);
		
		# Find the non-onverlapping combinations of domains with highest biscore
		find_highest_scoring_domains_combinations(\%gene, \%domains, \%domains_sortable_list);
		
		# Now, we will get the isoform with highest score. Then we will check if the templates of this isoforms 
		# can be transfer to other isoforms as a whole (all of them has to be transferable, otherwise it could be overlapping issues)
		# We tranfer them whenever it is possible. Therefore, all the sequence that have the same sequence in those pieces 
		# of sequence with PDB evidence will have the same score (as it should be).
		transfer_templates(\%gene, \%domains);
		
		# Report final scores for each isoform
		#ÊBEGIN: CHANGE for APPRIS
#		print_final_scoring($gene_id, \%gene, \%domains);

		$output .= print_final_scoring($gene_id, \%gene, \%domains);
		#ÊEND: CHANGE for APPRIS
	}
	
 	#ÊBEGIN: CHANGE for APPRIS
	return $output;
	#ÊEND: CHANGE for APPRIS
}


sub load_data{
	
	my $gene_id 				= shift;
	my $gene 					= shift;
	my $domains 				= shift;
	my $domains_sortable_list 	= shift;
	my $sequences 				= shift;
	my $tmp_dir 				= shift;
	
	foreach my $transcript_id (sort keys %{$sequences->{$gene_id}}){
		
		#ÊBEGIN: CHANGE for APPRIS
#		my $hmmscan_file = "$tmp_dir/$transcript_id/matador3D.domtblout.sorted";

		# get sequence
		my $seq = $sequences->{$gene_id}{$transcript_id};
		
		# Create cache obj
		my ($cache) = APPRIS::Utils::CacheMD5->new(
			-dat => $seq,
			-ws  => $tmp_dir			
		);		
		my ($seq_idx) = $cache->idx;
		my ($seq_dir) = $cache->dir;
						
		my $hmmscan_file = "$seq_dir/seq.matador3d";		
		#ÊEND: CHANGE for APPRIS
		
		if(! -e $hmmscan_file){
			print STDERR "Warning: No hmmscan file has been found for transcript $transcript_id\n";
			next;
		}
		
		my $sequence = $sequences->{$gene_id}{$transcript_id};
		
		if(!$sequence){
			print STDERR "Warning: No sequence found for transcript $transcript_id\n";
			next;
		}
		
		$gene->{$transcript_id}{seq} 			= $sequence;
		$gene->{$transcript_id}{cov} 			= [];
		$gene->{$transcript_id}{domains} 		= [];
		$gene->{$transcript_id}{global_score} 	= 0;
		
		open(FH, $hmmscan_file) || die "Error while opening file $hmmscan_file";
		while(my $line = <FH>){
			my @fields = split(/\s+/, $line);
			
			my $template 	= $fields[0];
			my $bitscore	= $fields[13];
			my $bias 		= $fields[14];
			my $start 		= $fields[17];
			my $end 		= $fields[18];
			my $num_domains	= $fields[10];
			
			if($bitscore < 1){
				next;
			}
			
			$domains->{$template}{$transcript_id}{$start}{$end}{bitscore} 		= $bitscore;
			$domains->{$template}{$transcript_id}{$start}{$end}{bias} 			= $bias;
			$domains->{$template}{$transcript_id}{$start}{$end}{num_domains} 	= $num_domains;
			
			$domains_sortable_list->{"$template|$transcript_id|$start|$end"} 	= $bitscore;
		}
		close FH;
	}
}



sub find_highest_scoring_domains_combinations{
	
	my $gene 					= shift;
	my $domains 				= shift;
	my $domains_sortable_list 	= shift;
	
	
	foreach my $domain (sort {$domains_sortable_list->{$b} <=> $domains_sortable_list->{$a}} keys %$domains_sortable_list){
		
		my $bitscore = $domains_sortable_list->{$domain};
		
		my @fields 		= split(/\|/, $domain);
		my $template 	= $fields[0];
		my $isoform 	= $fields[1];
		my $start 		= $fields[2];
		my $end 		= $fields[3];
		
		# If there fewer than ten position aligned, the domain is filtered
		if($end - $start + 1 < 10){
			next;
		}
		
		my $overlap = 0;
		# Check coverage on this isoform
		foreach my $range (@{$gene->{$isoform}{cov}}){
			my $cov_start 	= $range->[0]; 
			my $cov_end 	= $range->[1];
			
			if( !($cov_end < $start || $cov_start > $end) ){
				# There is overlap
				$overlap = 1;
				last;
			}		
		}
		
		if(!$overlap){
			$gene->{$isoform}{global_score} = $gene->{$isoform}{global_score} + $bitscore;
			push(@{$gene->{$isoform}{domains}}, $domain);
			push(@{$gene->{$isoform}{cov}}, [$start, $end]);
		}
		
		# Transfer to other isoforms with the same sequence if there no a other domain that covers this subsequence
		my $domain_seq = substr($gene->{$isoform}{seq}, $start - 1, $end - $start + 1);
		unless ( defined $domain_seq ) { next; }
		
		foreach my $other_isoform (keys %$gene){
			if($other_isoform eq $isoform){
				next;
			}
			
			my @sub_seqs = split(/$domain_seq/, $gene->{$other_isoform}{seq});
			
			my $num_subseqs = scalar @sub_seqs;
			

			
			# Possible enhancement: Take into account when the substring is found more than once
#			if($num_subseqs > 2){
#				print STDERR "More than two substring, gene $gene, isos $isoform $other_isoform, seq: $domain_seq\n";
#			}
			if($num_subseqs > 1){
				
				my $other_iso_start = length($sub_seqs[0]) + 1;
				my $other_iso_end 	= $other_iso_start + length($domain_seq) - 1;
				
				$overlap = 0;
				# Check coverage on the another isoform
				foreach my $range (@{$gene->{$other_isoform}{cov}}){
					my $cov_start 	= $range->[0]; 
					my $cov_end 	= $range->[1];
					
					if( !($cov_end < $other_iso_start || $cov_start > $other_iso_end) ){
						# There is overlap
						$overlap = 1;
						last;
					}		
				}
				
				if(!$overlap){
					my $other_domain = "$template|$other_isoform|$other_iso_start|$other_iso_end";
					$gene->{$other_isoform}{global_score} = $gene->{$other_isoform}{global_score} + $bitscore;
					$domains->{$template}{$other_isoform}{$other_iso_start}{$other_iso_end}{bitscore} 		= $bitscore;
					$domains->{$template}{$other_isoform}{$other_iso_start}{$other_iso_end}{bias} 			= $domains->{$template}{$isoform}{$start}{$end}{bias};
					$domains->{$template}{$other_isoform}{$other_iso_start}{$other_iso_end}{num_domains} 	= 0;
					push(@{$gene->{$other_isoform}{domains}}, $other_domain);
					push(@{$gene->{$other_isoform}{cov}}, [$other_iso_start, $other_iso_end]);				
				}
			}
		}
	}
}



sub transfer_templates{
	
	my $gene 	= shift;
	my $domains = shift;
	
	## Now, we will get the isoform with highest score. Then we will check if the templates of this isoforms can be transfer to other isoforms as a whole (all of them has to be transferable, otherwise it could be overlapping issues)
	# Tranfer them whenever it is possible. Therefore, all the sequence that have the same sequence in those pieces of sequence with PDB evidence will have the same score (as it shoudl be).
	
	# First, find the isoform (the first one if there are a draws):
	my $max_score 		= 0;
	my $max_isoform 	= undef;
	my $max_sequence 	= undef;
	foreach my $isoform (keys %$gene){
		my $global_score 	= $gene->{$isoform}{global_score};
		if($global_score > $max_score || !defined $max_isoform){
			$max_score 		= $global_score;
			$max_isoform 	= $isoform;
			$max_sequence 	= $gene->{$isoform}{seq};
		}
	}
	
	# We try to transfer the templates (domains) from the isoform with the highest score to other isoforms
	foreach my $isoform (keys %$gene){
		if($isoform eq $max_isoform){
			next;
		}
		
		my @cov;
		my @domains;
		my $transferables = 1; # Flag to check whether all templates are transferable or not
		foreach my $domain_to_transfer (@{$gene->{$max_isoform}{domains}}){
			# Find start and end of the domain
			my @fields 		= split(/\|/, $domain_to_transfer);
			my $template 	= $fields[0];
			my $start 		= $fields[2];
			my $end 		= $fields[3];
			
			# Find the sequence of the domain to transfer
			my $domain_seq 	= substr($max_sequence, $start - 1, $end - $start + 1);
			# Check if this sequence appear in the another isoform
			my @sub_seqs 	= split(/$domain_seq/, $gene->{$isoform}{seq});
	
			if(scalar @sub_seqs > 1){
				my $other_iso_start = length($sub_seqs[0]) + 1;
				my $other_iso_end 	= $other_iso_start + length($domain_seq) - 1;
				# This domain can be transfer, save info to make the transference
				push(@cov, [$other_iso_start, $other_iso_end]);
				push(@domains, "$template|$isoform|$other_iso_start|$other_iso_end");
				
				# Get the bias and bitscore of the max isoform
				my $bias 		= $domains->{$template}{$max_isoform}{$start}{$end}{bias};
				my $bitscore 	= $domains->{$template}{$max_isoform}{$start}{$end}{bitscore};
	
				# Save bias and bitscore in the another isoform			
				$domains->{$template}{$isoform}{$other_iso_start}{$other_iso_end}{bias} 	= $bias;
				$domains->{$template}{$isoform}{$other_iso_start}{$other_iso_end}{bitscore} = $bitscore;
			}
			else{
				# This domain can not be transfer. Stop this transference
				$transferables = 0;
				last;
			}
		}
		
		if($transferables){
			# Consolidate the transference
			$gene->{$isoform}{global_score} = $max_score;
			$gene->{$isoform}{domains} 		= \@domains;
			$gene->{$isoform}{cov} 			= \@cov;
		}
	}
}


sub print_final_scoring{
	
	my $gene_id = shift;
	my $gene 	= shift;
	my $domains = shift;
	
	#ÊBEGIN: CHANGE for APPRIS
	my $output = '';
	#ÊEND: CHANGE for APPRIS	
	
	foreach my $isoform (keys %$gene){
		
		my $global_score 	= $gene->{$isoform}{global_score};
		my $str_domains 	= join(';', @{$gene->{$isoform}{domains}});
		
	 	#ÊBEGIN: CHANGE for APPRIS
#		print "$gene_id\t$isoform\t$global_score";
		$output .= "$gene_id\t$isoform\t$global_score";
		#ÊEND: CHANGE for APPRIS
		
		
		foreach my $domain (@{$gene->{$isoform}{domains}}){
			
			my @fields 		= split(/\|/, $domain);
			my $template 	= $fields[0];
			my $isoform 	= $fields[1];
			my $start 		= $fields[2];
			my $end 		= $fields[3];
			
			my $bias 		= $domains->{$template}{$isoform}{$start}{$end}{bias};
			my $bitscore 	= $domains->{$template}{$isoform}{$start}{$end}{bitscore};
			
		 	#ÊBEGIN: CHANGE for APPRIS
#			print "\t$template;$bitscore;$bias;$start-$end;" . ($end - $start + 1);
			$output .= "\t$template;$bitscore;$bias;$start-$end;" . ($end - $start + 1);
			#ÊEND: CHANGE for APPRIS
		}
		
	 	#ÊBEGIN: CHANGE for APPRIS
#		print "\n";
		$output .= "\n";
		#ÊEND: CHANGE for APPRIS
	}
	
 	#ÊBEGIN: CHANGE for APPRIS
	return $output;
	#ÊEND: CHANGE for APPRIS	
}



# Read the input file. It returns a hash where the key is transcript_id, {gene} = gene_id, {sequence} = protein_sequence
# The expected format of the header is "> Protein_Id | Transcript_id | Gene_id | Gene_name | Seq_length". Whitepaces are ignored, only Transcript_id and Gene_id fields are used
sub read_fasta_file {
	
	my $fasta_file = shift;
		
	my %ret;
	my $new_header 	= "";
	my $old_header 	= "";
	my $sequence 	= "";
	open (IN, $fasta_file) or die "can't open $fasta_file: $!\n";
	while(my $line = <IN>){
		chomp $line;
		if($line =~ /^>(.+)/){
			my $new_header = $1;
			# Remove whitespaces
			if(!$old_header){
				$old_header = $new_header;
				next;
			}
			
			$old_header =~ s/\s+//g;
			
			my @fields 		= split(/\|/, $old_header);
			
			if(scalar @fields < 2){
				die "The header $line in file $fasta_file does not contain | characters. The expected format is \"> Protein_Id | Transcript_id | Gene_id | Gene_name | Seq_length\". Whitepaces are ignored, only Transcript_id and Gene_id fields are used\n";
			}
			
			#ÊBEGIN: CHANGE for APPRIS
			my $protein_id 		= $fields[0];
			my $transcript_id 	= $fields[1];
			my $gene_id 		= $fields[2];
			#ÊEND: CHANGE for APPRIS
			
			# Check transcript_id is not repeated, otherwise print warning and discard current sequence
			if($sequence){
				if(defined $ret{$gene_id}{$transcript_id}){
					print STDERR "Warning: Two sequences have the same transcript_id: $old_header. Only the first one will be processed\n";
					$sequence = "";
				}
				else{
					$ret{$gene_id}{$transcript_id} = $sequence;
				}
			}              
			
			$old_header = $new_header;
			$sequence 	= ""; # clear out old sequence
		}         
		else{    
			$line =~ s/\s+//g; # remove whitespace
			$sequence .= $line; # add sequence
		}         
	}    
	
	close IN;
	
	if($sequence){ # handle last sequence
		my @fields 			= split(/\|/, $old_header);
		#ÊBEGIN: CHANGE for APPRIS
		my $protein_id 		= $fields[0];
		my $transcript_id 	= $fields[1];
		my $gene_id 		= $fields[2];
		#ÊEND: CHANGE for APPRIS
		
		# Check transcript_id is not repeated, otherwise print warning and discard current sequence
		if(defined $ret{$gene_id}{$transcript_id}){
			print STDERR "Warning: Two sequences have the same transcript_id: $old_header. Only the first one will be processed\n";
			$sequence = "";
		}
		else{
			$ret{$gene_id}{$transcript_id} = $sequence;
		}
	}
	
	return \%ret;
}





sub usage {
	
	print STDERR "Usage: \nmatador3d.pl -i input_fasta_file -bin hmmscan_bin -db database [-t temp_dir]\n";
	print STDERR "Description:\nGiven a input fasta file, it estimates how well each sequence is covered by a 3D structure in the PDB. To that aim, each sequence is searched against a database of sequence profiles built from proteins with 3D structure.\n";
	print STDERR "Parameters:\n";
	print STDERR "-i -> Input fasta file. The format of the header is expected as \"> Protein_Id | Transcript_id | Gene_id | Gene_name | Seq_length\". Whitepaces are ignored, only Transcript_id and Gene_id fields are used within the program\n";
	print STDERR "-bin -> Path to hmmscan binary.\n";
	print STDERR "-db -> Path to HMM database.\n";
	print STDERR "-t -> Temporary directory. Temporary files will be generated and stored in this directory. Optional, default is /tmp.\n";
	
	exit;
}




#ÊBEGIN: CHANGE for APPRIS
#
# OUTPUT DESCRIPTION:
# Description of columns
#
# 1.  Gene identifier
# 2.  Transcript/Protein identifier
# 3.  Total bitscore from every HMMER alignment.
# 4.  Result from HMMER Alignment:
# 4.1 PDB identifier
# 4.2 Bitscore from alignment
# 4.3 Bias of aligment: bias of the distribution of amino acid compositions calculated by HMMER
# 4.4 Start-End of query position
# 4.5 Absolute value of aligmnet
#

1;

#ÊEND: CHANGE for APPRIS