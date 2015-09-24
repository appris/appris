=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::MotifTools

=head1 DESCRIPTION

This module collates a variety of useful methods for dealing with transcription
factor binding motifs (PFMs/PWMs).


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw(list your required methods here);

=cut

###############################################################################

package Bio::EnsEMBL::Funcgen::Sequencing::MotifTools;

use warnings;
use strict;

use File::Basename qw( basename fileparse );
use DBI            qw( :sql_types );
use base           qw( Exporter );
use vars           qw( @EXPORT );

use Bio::EnsEMBL::Utils::SqlHelper;
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd
                                               run_backtick_cmd
                                               check_file
                                               open_file
                                               which_path );
@EXPORT = qw(
  filter_pwm_mappings
  get_revcomp_file_path
  matrix_min_max
  parse_matrix_line
  read_matrix_file
  reverse_complement_matrix
  revcomp_matrix_file
  run_moods
  sprint_matrix
  write_matrix_file
);


#Make %swap  global ReadOnly
#so we don't have to instantiate it for every call
#ReadOnly is not a core module though
my %swap = (0 => 3, 1 => 2, 2 => 1, 3 => 0);
# assumes a jaspar matrix pfm file with >header and 4 rows A C G T

#Dispatch table
my %method_refs = ('parse_find_pssm_dna_to_bed' => \&parse_find_pssm_dna_to_bed);

#integrate matrix_min_max in here


=head2 read_matrix_file

  Arg [1]     : String - The input file path i.e. a Japar pmf file.
  Arg [2]     : Arrayref (optional) - A list of matrix IDs to read from file
  Description : Reads a matrix file and return each frequency or weight as 
                an element of an array. Can process multiple matrices within 
                a headered file, or single matrix from a header less file.
                In the case of a headerless file, the file name is used to 
                generate the hash key.
  Returntype  : Hashref - ID line/File prefix hash keys. Values are 3d array. 
                The first element being an arrayref of a 2d array representing
                the matrix. If the boolean argument is passed, the second 
                element will represent the reverse complement  
  Example     : my $matrix_info   = read_matrix_file('/path/to/MA0012.1.pfm', 1);
                my ($matrix_ref, $rev_comp_matrix_ref) = @{$matrix_info{'MA0012.1'};
  Exceptions  : Throws if any matrix does not have 4 rows or no matrices exist.
                Throws if headers are absent from a multi-matrix file.
  Status      : At risk

=cut

# Also parse out version?

sub read_matrix_file{
  my $file             = shift;
  my $matrix_ids_array = shift;
  my ($num_matrices, $matrix_ids);
  
  if (defined $matrix_ids_array){
    assert_ref($matrix_ids_array, 'ARRAY', '$matrix_ids_array');
    $num_matrices = scalar(@$matrix_ids_array);
    #convert to hash for speed
    $matrix_ids = {map {$_ => undef} @$matrix_ids_array};
  } 

  my $in_file = open_file($file, '<');
  my ($line, $header, $matrix_id, %matrix_cache, @matrix_tmp);

  while( ($line = $in_file->getline) && 
         defined $line){
    chomp $line;

    if($line =~ /^>/){ #e.g. >MA0012.1 br_Z3 

      if(@matrix_tmp && ! $header){ #Process and cache the matrix
        throw("Failed to find header for matrix before:\n\t$line");
      }
      elsif($header && (scalar(@matrix_tmp) != 4)){
        throw('Found incorrect number of matrix rows('.scalar(@matrix_tmp).
          ") for $header in file:\t$file");
      }
      elsif(@matrix_tmp){
        ($matrix_id = $header) =~ s/^>([^ ]+)\s*.*$/$1/o;
        #Strip off any other header info e.g. > and protein name

        if((! $matrix_ids) ||
           exists $matrix_ids->{$matrix_id}){
          #last before we populate $matrix_cache
          #as last one is done below
          if ($num_matrices){
            $num_matrices--;
            last if ! $num_matrices;
          }          

          #warn "setting $matrix_id matrix to ".\@matrix_tmp."\n\t".
          # join("\n\t", (map { join(' ', @$_)} @matrix_tmp))."\n";

          #Deref here as all records were getting the same 
          #matrix, i.e. the last one.
          $matrix_cache{$matrix_id} = {matrix => [@matrix_tmp],
                                       header => $header};        
        }  
      }
  
      undef @matrix_tmp;
      $header = $line;
    }
    else{
      #eval thsi and write $id, to avoid having to pass $id, as parse_matrix_line
      #is used in web code
      push @matrix_tmp, [parse_matrix_line($line)]; #, $matrix_id);
    }
  }

  close($in_file);

  if(scalar(@matrix_tmp) != 4){
    throw('Found incorrect number of rows('.scalar(@matrix_tmp).
      ") for $header in file:\t$file") ;
  }


  if((! $matrix_ids) ||
           exists $matrix_ids->{$matrix_id}){

    if(! defined $header){
      #%matrix_cache must be empty
      #Must be a single matrix in a headerless file
      #Get the matrix name from the file
      ($header = basename($file)) =~ s/\.p[fw]m$//io;
      $matrix_id = $header;
    }

    $matrix_cache{$matrix_id} = {matrix => \@matrix_tmp,
                                 header => $header};
  }

  return \%matrix_cache;
}

# Also used in BindingMatrix

sub parse_matrix_line{
  my $line = shift;
  #my $id   = shift;
  # remove leading whitespace or A|C|G|T [ ] wrapping
  (my $clean_line = $line) =~ s/^\s*[ACGTacgt]?\s*\[?\s*([\.0-9 ]+[\.0-9])/$1/o;
  $clean_line =~ s/\s*\]?\s*$//o;
  #(my $clean_line = $line) =~ s/^\s*[ACGTacgt]?\s*\[?\s*([0-9 ]+)\s*\]?\s*$/$1/o;

  if($clean_line =~ /[^\s0-9\.]/){
    throw("Found invalid characters in matrix line:\n$line");
  }

  return split(/\s+/, $clean_line);
}


sub reverse_complement_matrix{
  my $matrix = shift;
  assert_ref($matrix, 'ARRAY', '$matrix');
  #Skipping size check as probably already done in read_matrix_file;
  my @rc_matrix;

  foreach my $row(0..$#{$matrix}){
    $rc_matrix[$swap{$row}] = [reverse(@{$matrix->[$row]})]; 
  }

  return \@rc_matrix;
}


#This does not handle revcomp
#although we just pass that as the 'matrix'
#would likely have to alter header too.

sub write_matrix_file{
  my $matrix_hashes = shift;
  my $out_path      = shift;
  my $write_header  = shift;
  assert_ref($matrix_hashes, 'ARRAY', '$matrix_hashes');
  $write_header     = 1 if scalar(@$matrix_hashes) > 1;
  my $out_file      = open_file($out_path, '>');

  foreach my $mhash(@$matrix_hashes){
    print $out_file $mhash->{header}."\n" if $write_header;
    print $out_file sprint_matrix($mhash->{matrix});
  }

  close($out_file);
}


sub sprint_matrix{
  my $matrix = shift;
  assert_ref($matrix, 'ARRAY', 'Matrix array');
  my $mstring;

  for my $row(@$matrix){
    assert_ref($row, 'ARRAY', 'Matrix row array');
    $mstring .= join(" ", @{$row})."\n";
  }

  return $mstring;
}



sub get_revcomp_file_path{
  my $pfm_file_path = shift;
  my $out_dir       = shift;

  my ($pfm_file_name, $pfm_dir) = fileparse($pfm_file_path);  
  $out_dir ||= $pfm_dir;
  my $rc_file  = $out_dir.'/revcomp.'.$pfm_file_name;  
}

sub revcomp_matrix_file{
  my $pfm_file_path = shift;
  my $out_dir       = shift;
  my $write_header  = shift;

  my $matrix_info = read_matrix_file($pfm_file_path);
  $write_header   = 1 if scalar(keys %$matrix_info) > 1;
  my $rc_file     = get_revcomp_file_path($pfm_file_path, $out_dir);
  #add matrix_min_max functionality to read_matrix if still required.

  my $out_file = open_file($rc_file, '>');

  foreach my $mhash(values %{$matrix_info}){
    print $out_file $mhash->{header}."\n" if $write_header;

    for my $row(reverse_complement_matrix($mhash->{matrix})){
      print $out_file join(" ", @{$row})."\n";
    }
  }

  close($out_file);
  return $rc_file;
}


#To do
#1 Investigate usage of perl interface? This would remove intermediate output/parsing
# and deliver results direct to process enabling custom/standard output format
# Would also remove the need to write the rc matrix files, as we could do that here


# Currently doesn't allow over-ride of find_pssm_dna path, relying on $PATH

#Must be single chr fasta file! #As target fasta header not reported in output
#Validate we only have one header in the file?



sub run_moods{
  my ($p_thresh, $target, $queries, $out_prefix, $format, $mapper_params, $debug) =
   rearrange( [qw(p_thresh target queries out_prefix format mapper_params debug) ], @_);

  assert_ref($queries, 'ARRAY', 'Mandatory param -queries');    
  throw("Mandatory param -out_prefix not defined")                 if ! defined $out_prefix;
  #throw("Mandatory param -mapper not defined or does not exist")   if ! -f $target;
  throw("Mandatory param -queries must has at least 1 query file") if ! @$queries;
  map { throw("Query file is not defined or does not exist")       if ! -f $_ } @$queries;
  $mapper_params ||= [];
  assert_ref($mapper_params, 'ARRAY', 'Mapper params');
  $p_thresh = 0.001 if ! defined $p_thresh;

  my $mapper = 'find_pssm_dna';
  my $parse_method;

  if(defined $format){
    $parse_method = 'parse_find_pssm_dna_to_'.lc($format);

    if(! exists $method_refs{$parse_method}){
      throw("Unsupport output format:\t$format");
    }

    $parse_method = $method_refs{$parse_method};
  }


  #Validate we only have 1 header in fasta file.
  #As find_dna_pssm does not report fasta haeder in hit info
  if(run_backtick_cmd("grep '^>' $target | wc -l") > 1){
    throw('Found more than one header in target fasta file. Moods find_pssm_dna requires single sequence target fasta files.')
  }

  #Assume the file naming is correct
  #Could validate vs actual returned headers?
  (my $chr = basename $target) =~ s/\.fa(sta)*$//;


  #Or use validate_package_path
  #This is overkill, although maybe essentail if find_pssm_dna isn't in top level $PATH dir 
  #Will also ensure full path of mapper is used in cmdline
  #only if not already a path?
  #my $mapper_path = ($mapper =~ '/') ? $mapper : which_path($mapper);
  #if(! defined $mapper_path){
  #  throw("Failed to find path to mapper:\t$mapper\nMapper is incorrect or \$PATH needs configuring");
  #}

  #Enable running generic mapper with @mapper_params only?
  #No point as the utility here is to provide inline post-processing
  #of known mapper
  #Otherwise we could just run the mapper directly

  #Currently not using classes for PWMMappers, but could move towards same model as 
  #for Aligners

  #any point in calling run_find_dna_pssm
  #or can we just run directly from here
  #The reason is to handle the translation between params an cmdline options
  #Hence run_mapper is a little redundant
  #we should just have run_$mapper?
  #These would in effect replace the Mapper::run methods

  #But we probably want to maintain a similar interfact to the SeqlTools::run_aligner
  #method. That does some pre-validation, creates and runs the Aligner

  my $out_file = $out_prefix.'.'.$chr.'.out';


#iterate over queries here

  #Also need to generate jaspar $dbh if we are parsing the data and integrating the matric name.
  my $cmd = "$mapper -f $p_thresh ".join(' ', @$mapper_params)." $target ".
   join(' ', @$queries)." > ${out_file}.tmp";

  run_system_cmd($cmd);
  run_system_cmd("mv ${out_file}.tmp $out_file");

  if(defined $parse_method){
    my $reformated_out_file = $parse_method->($out_prefix.'.'.$chr, $chr);#, $jdb);
    run_system_cmd("rm -f $out_file");
    $out_file = $reformated_out_file;
  }


  #Validate we have some records?
  #STDERR will have some output, but find_pssm_dna
  #will not die if it fails to read a matrix or returns no mapping
  #throw is we have 0 output?

  return $out_file;
}

#This is jaspar specific at present due to access of Jaspar DB to get matrix name

sub parse_find_pssm_dna_to_bed{
  my $prefix = shift;
  my $chr    = shift;
  #my $jdb    = shift;

  if(! ($prefix && $chr)){
    throw('Must pass an output prefix and seq_region name arguments');
  }

  #assert_ref($jdb, 'Bio::EnsEMBL::DBSQL::DBConnection', '$jdb');

  my $in_fh  = open_file($prefix.'.out');
  my $out_fh = open_file("${prefix}.bed.tmp", '>');
  
  #Individual input pfm files do not have headers!!!!!
  #The results report the input pfm files name in the header
  #Also doesn't report target fasta header, so need to be run on 1 chr at a time!
  #such that find_dna_pssm will list hits based on the file name rather than the 
  #matrix header line?
  #This adds more to the argument for moving to the Moods.pm interface
  #Example find_pssm_dna output

  #>MA0376.1.pfm
  #Hits:2649 Length:20
  #3002127 5.18373
  #3003071 6.17587
  #3003490 4.34327
  #3004276 6.03918
  #3009721 4.39983
  #3012854 4.82056
  #3015039 7.59674
  #3015235 5.16521

  my ($strand, $id, $pfm_length, $tf_name, $line);
  #my $score_thresh; #my $cnt = 0;
  #my $sql = "SELECT NAME from MATRIX where BASE_ID=?";
  #my $sth = $jdb->prepare($sql);

  while(($line = $in_fh->getline) && 
        defined $line){
    chomp $line;

    if($line eq ''){ #often empty lines at EOF?
      next;
    }
    elsif($line =~ /^>/){ #Parse pfm hits header
      #$cnt ++;
      #Get the matrix ID
      ($id = basename($line)) =~ s/>//; #In case we have specified a path
      $id =~ s/\.pfm$//;  
      #$sth->execute($id) or die $sth->errstr;
      #($tf_name) = $sth->fetchrow_array;
      #This assumes that all rows with the same base id share the same TF name.
      #Although there appears to 1 entry where this is not true, probably a typo?
      #> select * from MATRIX m join MATRIX m1 using(BASE_ID) where m.name != m1.NAME;
      #+---------+-------+------------+---------+--------+-------+------------+---------+--------+
      #| BASE_ID | ID    | COLLECTION | VERSION | NAME   | ID    | COLLECTION | VERSION | NAME   |
      #+---------+-------+------------+---------+--------+-------+------------+---------+--------+
      #| MA0110  |  9339 | CORE       |       1 | ATHB-5 | 10733 | CORE       |       2 | ATHB5  |
      #| MA0110  | 10733 | CORE       |       2 | ATHB5  |  9339 | CORE       |       1 | ATHB-5 |
      #+---------+-------+------------+---------+--------+-------+------------+---------+--------+
      
      #Get the TF name
      #if(! defined $tf_name){ warn "Failed to get TF name from the DB:\t$id\n"; }
  
      #Set the strand
      $strand = '1';

      if($id =~ /^revcomp\.(.+$)/){
        $strand = '-1';
        $id = $1;
      }
      
      #$score_thresh = $max_scores_ref->{$id.'.pfm'} * $score_thresh_perc; # 70% of max
      #$score_thresh = 0; #???!!!
      #commentary( $max_scores_ref->{$id.'.pfm'} ." = max score for $id $tf_name\n") if $verbose;
      
      #Parse the pfm length
      $line = $in_fh->getline;
      ($pfm_length) = $line =~ /.*Length:([0-9]+)/;
    }
    else{ #Parse hits for this pfm
      my ($start, $score) = split("\t", $line);
      #if($score > $score_thresh){
      #why $start + 1?
      #why do we need the TF name here?
      #It's totally redundant, and we can grab it later anyway
      #This removes the need for the jdb
      #Is this just here to pad out to bed format?
      #Shouldn't this be bed rather than tab?
      #we're using bedtools in pwm_filter_mappings.pl
      #$id is also redundant

      #TODO Create cache here to buffer output

      #This is already half open coords
      #Maintain this as we are using bedtools for merge and overlap
      print $out_fh join("\t", 
        ($chr, $start, ($start + $pfm_length), $id, $score, $strand))."\n";

    }
  }

  $in_fh->close;;
  $out_fh->close;
  run_system_cmd("mv ${prefix}.bed.tmp ${prefix}.bed");

  #Could return $cnt here, but we don't do that for out file in run_moods
  return "${prefix}.bed";
}


#The reason we were using $perc_of_max here is because 
#the background rate is also a percentage?
#But this may not scale if we simply re-run this filtering on some 
#new mappings with higher scores. 

#todo this could be done quicker in awk?

sub filter_pwm_mappings{
  #my $perc_of_max  = shift;
  my $threshold    = shift;
  my $in_file      = shift;
  my $out_file     = shift;
  
  #my $max = run_backtick_cmd("cut -f 5 $in_file | sort -g -u | tail -1");
  #chop $max;
  #my $threshold = $max * $perc_of_max/100;
  #print "$mat $max $thresh\n" if $verbose;
  my $ifh = open_file($in_file);
  my $ofh = open_file($out_file, '>');
  
  while(my $line = <$ifh>){
    my (undef, undef, undef, undef, $score) = split("\t",$line);
    
    if($score >= $threshold){
      print $ofh $line;
    }
  }
      
  close($ifh);
  close($ofh);

  return;
}



1;