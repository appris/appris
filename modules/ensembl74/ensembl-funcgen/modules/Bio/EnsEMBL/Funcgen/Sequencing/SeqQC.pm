=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::Sequencing::SeqQC

=head1 DESCRIPTION

This module collates a variety of quality control utilities for sequence alignments
and peaks calling


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Utils::Sequencing qw(list your required methods here);



=cut

###############################################################################

package Bio::EnsEMBL::Funcgen::Sequencing::SeqQC;

use warnings;
use strict;

use File::Find; 
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw      );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd 
                                               open_file 
                                               which_path);
use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw( 
  run_QC
  run_spp_phantom
  );



#todo
#implement so infrastructure to for running all default QC for given stage
#will likely need a hash with some label keys:
# fastq, raw_bam, filtered_bam, peaks etc.
#Write wrapper methods for these if required, such that we don't need to handle option
#setting

#Return these methods as a list, such that they can either be run in series or submitted to the farm to run in parallel

#Input params should be standardised across methods with same qc type/label


sub run_QC{
    
  
  
}







=head2 run_spp_phantom

  Arg [1]     : String  - input.bam file path
  Arg [2]     : Hashref - Params: 
                  {-out_prefix => <String>,   #Output file prefix. Default = input.bam
                   -quality    => <Int>,      #Acceptable quality threshold [-2..2]
                   -out_dir    => <String>,   #Output directory(-odir), default is 
                                              directory of input bam file. 
                   -no_dups    => <Boolean>,  #For bam files where duplicates have been removed.
                                              Runs run_spp_nodups.R instead of run_spp.R.
                   -force      => <Boolean>}  #Force flag. Over-writes previous output
                                              by spcifying -rf
  Arg [3]     : Hashref - Raw SPP params. These can be any valid SPP command line
                options, with optional values. These will over-ride any of the above 
                params. The '=' does not need to be specified for key value pairs 
                e.g.
                  -odir=/your/out/put/directory
                  Should be specified as follows:
                  {'-odir' => '/your/out/put/directory'}             
  Description : A wrapper method to run the phantompeakqualtools of the SPP package.
                Defaults outputs will be:
                  $output_dir/$output_prefix.spp.txt
                  $output_dir/$output_prefix.spp.pdf
                
                Note: SPP phantom uses quite a lot of tmp space, roughly the same 
                amount of the input bam file. Hence it may be necessary to specify this 
                in any resource requirements. 
                  e.g. For LSF: bsub 
                  
                For more details about command line options and requirements see:
                  https://code.google.com/p/phantompeakqualtools/
  Returntype  : Int - SPP Quality tag [-2..2]
  Example     : my $qc_passed_results = run_spp_phantom($bam_file, {quality => 1}); 
  Exceptions  : Throws if input.bam not valid.
                Throws if acceptable quality threshold is defined and not met.
  Status      : At risk

=cut


#unix,bash,R-2.10 and above,awk,samtools,boost C++ libraries
#R packages: SPP, caTools, snow
#NOTE: The current package does not run on a MAC or WINDOWS.

#run post filtering with a threshold of 1 or higher
#if we get failures then we can always go back and look at the unfitlered bams if need be
#to look at pure library quality.


#8GB? Probably a lot less than that.

#Max 
#2MB for 459M

#We may want to define an RScript path to override the default
#This can be done outside of perl by appropriately order the $PATH
#But this may not always be possible

sub run_spp_phantom{
  my $in_bam     = shift;
  my $params     = shift || {};
  my $raw_params = shift || {}; #overrides $params
  my $result = 1;
  
  #VALIDATE PARAMS
  if(! defined $in_bam ){
    throw('Bam file argument is not defined');      
  }
  elsif(! -f $in_bam){
    throw("Bam file argument is not a file:\t".$in_bam);
  }
  
  assert_ref($params, 'HASH');
  assert_ref($raw_params, 'HASH');
  
  #PROCESS PARAMS
  #validate all param key start with -, 
  #otherwise rearrange will simply return them (hangover from CGI origins)
  #reported to core team
  
  if(my @non_hyphenated = grep {/^[^-]/} keys %$params){
    throw("Params must be prefixed with a hyphen:\t@non_hyphenated");  
  }
  
  my ($out_pfx, $req_qual, $out_dir, $no_dups, $force, $debug) = 
    rearrange([qw(out_prefix quality out_dir no_dups force debug )], %$params);
  
  #Generate defaults, spp would use these anyway, 
  #but we need them in the error output for clarity
  ($out_pfx = $in_bam) =~ s/.*\///o     if ! defined $out_pfx;
  ($out_dir = $in_bam) =~ s/$out_pfx//  if ! defined $out_dir;
  $out_dir ||= '.';
  
  #SET DEFAULT OPTIONS
  if(! keys(%$raw_params)){ #Set default save option
    $raw_params->{'-savp'} = "${out_dir}/${out_pfx}.spp.pdf"; 
  }
  
  if(! exists $raw_params->{'-odir'}){
    #This maybe an empty string if we have passed a bam file in the cwd
    $raw_params->{'-odir'} = $out_dir if $out_dir;   
  }

  if( ! exists $raw_params->{'-out'}){
    $raw_params->{'-out'} = "${out_dir}/${out_pfx}.spp.txt";      
  }
  
  $raw_params->{'-rf'} = undef if $force;
    
  #BUILD THE CMD 
  my $cmd = $no_dups ? 'run_spp_nodups.R' : 'run_spp.R';
  $cmd    = which_path($cmd);
  $cmd    = "Rscript $cmd -c=${in_bam}";
  
  foreach my $opt(keys(%$raw_params)){
    my $param = ' '.$opt;
    $param   .= '='.$raw_params->{$opt} if defined $raw_params->{$opt};
    $cmd     .= $param;
  }
  
  #It's a good idea to test the outputs will be writeable before wasting time
  
  #RUN THE CMD RETURN RESULTS
  warn "cmd is $cmd" if $debug;
  run_system_cmd($cmd);
  #could really do with passing error string back if this fails
  
  
  my $fh = open_file($raw_params->{'-out'});
  my $output = <$fh>;
  $fh->close;   
  chomp $output;
  my @output = split("\t", $output);
  #$output[0] is filename  
  my %results = (numReads         => $output[1],
                 estFragLen       => $output[2],
                 corr_estFragLen  => $output[3],
                 phantomPeak      => $output[4],
                 corr_phantomPeak => $output[5],
                 argmin_corr      => $output[6],
                 min_corr         => $output[7],
                 NSC              => $output[8],
                 RSC              => $output[9],
                 QualityTag       => $output[10]);
  
  $out_dir = $raw_params->{'-out'}.'/' || '';
  
  #todo change this so we don't throw, but we return a fail status
  
           
  if(defined $req_qual && 
     ($results{QualityTag} < $req_qual)){
    $result = 0;   
    warn("SPP:PhantomPeakQualTools failed on QualityTag:\tRequired $req_qual\tObserved ".
      $results{QualityTag}."\nInput file:\t$in_bam\nResults:\t${out_dir}${out_pfx}.spp.txt\n".
      "Plot:\t${out_dir}/${out_pfx}.spp.pdf\n");     
  }
  
  #return like as in the scalar context the last value of the array is returned
  #meaning we can do something like this in the caller
  #if(! run_spp_phantom()){ #do something }
  
  return \%results, $result;
}


#fastqc online docs seem focused on gui usage only, 
#cmdline options are not available on line only via man page?
#


#This can run fastqs (given a tracking DB) and/or alignment files
#get_alignment_files_by_formats_ResultSet is a Hive method
#can we move get_alignment_file_prefix_by_ResultSet outside of the hive and pass the relevant alignment dir?
#also need to test ALIGNED states

#where is out_dir going to be, it can't be in the fastq dir? 
#as this may be of unkown structure and inputs may be disparately located
#pre_db_QC dir in output?
#actually, we need non-DB specific dirs to capture all QC output and persistant files?
#this may clashe if we have different DBs processing the same peak set
#we just accomodate this clash in the alignments
#do we need to make alignment output DB specific too?

#we need align_dir to create bam files
#and tracking DB to get fastqs

sub run_fastqc_for_ResultSet{
  my $rset   = shift;
  my $params = shift || {};
  assert_ref($params, 'HASH');
  assert_ref($rset, 'Bio::EnsEMBL::Funcgen::ResultSet');
  my ($tracking_adaptor, $align_root, $debug) = 
    rearrange([qw(tracking_adaptor alignment_root_dir debug throw_on_fail)], %$params);
  my (%results, $result);  
  
  #fastq first
  if($tracking_adaptor){
    assert_ref($tracking_adaptor, 'Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor');
    $rset->is_stored($tracking_adaptor->db);
  
  
  }
  
  #Now bam (or bam_mapped? to be inline with fastqc?)
  if($align_root){
    
  }
  
  #return like this as in the scalar context the last value of the array is returned

   #meaning we can do something like this in the caller
   #if(! run_spp_phantom()){ #do something }
   
  
  
   return \%results, $result;
}
 


sub run_fastqc{
  my $in_files = shift;
  my $params   = shift || {};
  assert_ref($params, 'HASH');

  my ($out_dir,  $format, $debug) = 
    rearrange([qw(out_dir format debug throw_on_fail)], %$params);
  my (%results, $result);  
    

  if(! defined $in_files){
    throw('Input file(s) argument not defined, mus be a single file path or an array ref of file paths');  
  }
  elsif(ref($in_files)){
    assert_ref($in_files, 'ARRAY', 'fastqc input files');  
  }
  else{
    $in_files = [$in_files];  
  }
  


  #fastqc isn't normally installed into the top of the path 
  #a symlink is required, so let's do some whiching to see whether we have it in the $PATH
  #else call _which_path?
    
  #my $cmd = ;
  #warn "DEACTIVATED FASTQC FOR NOW:\nfastqc -f fastq -o ".$self->output_dir." @fastqs";
  #run_system_cmd('fastqc -o '.$self->output_dir." @fastqs"); 
  
  #run_system_cmd('fastqc -o '.$self->output_dir." @fastqs"); 
  #-f seems to be optional, and only acts to over-ride auto-detection

  
  #-t --threads    Specifies the number of files which can be processed
  #                simultaneously.  Each thread will be allocated 250MB of
  #                memory so you shouldn't run more threads than your
  #                available memory will cope with, and not more than
  #                6 threads on a 32 bit machine

  run_system_cmd('fastqc -extract -o '.$out_dir.' '.join(' ', @$in_files)); 
  
   
  #grep FAIL, WARN, PASS jobs


  #return like this as in the scalar context the last value of the array is returned
  #meaning we can do something like this in the caller
  #if(! run_spp_phantom()){ #do something }
  
  #change this to hashref keyed on test name?
  #or test status?
  
  return \%results, $result;
}



#This actually dumps data from the DB!
#We need to change this to be able to run off bed files
#as we need a peaks report for rep peaks, which are not loaded into the DB!

sub run_peaks_report{
  my $peak_file = shift;
    
  #  name=peaks_${DB_NAME}.report
  #report_file=${WORK_DIR}/${name}.pdf
  #echo ": Generating Report: $report_file"  
  #BackUpFile $report_file

  #cmd="perl $EFG_SRC/scripts/peaks_report.pl $DB_READ_SCRIPT_ARGS $DNADB_SCRIPT_ARGS -R -compare -feature_table annotated -no_outliers -name $name -outdir $WORK_DIR $sets"
  #echo $cmd

=pod
  
my ($species, $help, $R, $nodump, $compare, $regstats, $all_seq_regions, $no_outliers, $name, $outdir);
my ($feature_table, $host, $port, $user, $pass, $dbname, $dnadbhost, $dnadbport, $dnadbuser, $dnadbname, $dnadbpass);
my ($inset_main, $inset_compare);

#Default values
$inset_main=0.1;
$inset_compare=0.5;
$user = 'ensro';
$name = 'peaks_report_'.$$;#Add PID to avoid overwriting previous reports
$outdir = $ENV{'EFG_DATA'}."/output/".$ENV{DB_NAME};

#get command line options

my (@fset_names_tmp);

print "peaks_report.pl @ARGV\n";

GetOptions 
 (
      'species=s'          => \$species,
      'dnadb_host=s'       => \$dnadbhost,
      'dnadb_user=s'       => \$dnadbuser,
      'dnadb_port=i'       => \$dnadbport,
      'dnadb_pass=s'       => \$dnadbpass,
      'dnadb_name=s'       => \$dnadbname,
      'dbhost=s'           => \$host,
      'dbuser=s'           => \$user,
      'dbport=i'           => \$port,
      'dbpass=s'           => \$pass,
      'dbname=s'           => \$dbname,
      'outdir=s'           => \$outdir,
      "help|h"             => \$help,
      "R"                  => \$R,
      "nodump"             => \$nodump,
      "no_outliers"        => \$no_outliers,
      "compare"            => \$compare,
      "regstats"           => \$regstats,
      "all_seq_regions"    => \$all_seq_regions,
      "feature_sets=s{,}"  => \@fset_names_tmp,
            "feature_table=s",   => \$feature_table,
                  "inset_main=s",      => \$inset_main,
                  "inset_compare=s",   => \$inset_compare,
      "name=s"             => \$name,
       )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help);

#Reset to undef so we don't try with empty string
#$pass      ||= undef;
#$dnadbpass ||= undef;

# Sould be failing a little nicer now... 
if(!$feature_table) { print "Missing Type of Feature: annotated or regulatory (use -h for help)\n"; exit 0; }
if(!$host || !$port || !$user || !$dbname )  {  print "Missing connection parameters (use -h for help)\n"; exit 0; }
if(!$outdir )  {  print "\$EFG_DATA not defined and -outdir not specified\n"; exit 0; }

if(! $feature_tables{$feature_table}){
  die("You have specified an invalid -feature_table. Must be one of:\t".join("\t", (keys %feature_tables)));
}

if(! -d $outdir){
  die("Error: $outdir is not a valid output folder");
} else{
  print "Setting default output directory to:\t".$outdir;
  if(! -d $outdir){
    system("mkdir -p $outdir") == 0 or 
      die("Could not create output directory:\t".$outdir);
  }
}

#warn "dbpass is x${pass}x";
#warn "dnadbpass is x${dnadbpass}";

#Check database connections
my ($coredba, $efgdba);
if($dnadbname){
  
  my $coredba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $dnadbhost,
     -port => $dnadbport,
     -user => $dnadbuser,
     -dbname => $dnadbname,
     -species => $species,
     -group   => 'core',
     -pass    => $dnadbpass,
    );
}

 
$efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host    => $host,
   -port    => $port,
   -user    => $user,
   -dbname  => $dbname,
   -species => $species, #Not strictly necessary
   -dnadb   => $coredba, 
   -pass    => $pass,
  );


#Test connections
$efgdba->dbc->db_handle;
$efgdba->dnadb->dbc->db_handle;

my $fsa = $efgdba->get_FeatureSetAdaptor();

my @fset_names;

if(scalar(@fset_names_tmp)==0){
  map { push @fset_names, $_->name } @{$fsa->fetch_all_by_type($feature_table)};
}
else{#Validate fset names
  
  foreach my $fsname(@fset_names_tmp){
    my $fset = $fsa->fetch_by_name($fsname);
  
    if($fset){
      push @fset_names, $fset->name;
    } else {
      die("Could not fetch FeatureSet:\t$fsname\n");
    }
  }
}

#print the data of each set to individual files (maybe put all in one file??)
#give other saving options (e.g. clean-up, backup?)


my @sr_types = ('chromosome');
push @sr_types, 'non_chromosome' if $all_seq_regions;


my %sr_type_clauses = (
           'chromosome' => , "='chromosome'",
           'non_chromosome' => , "!='chromosome'",
          );

if(!$nodump){
  print "\n\n:: Dumping Datasets\n";
  
  #This was not accounting for nr sr_ids
  foreach my $sr_type(@sr_types){
    
    #Save to only one file... though it may be big...
    my $query = 'SELECT fs.name as 'name', s.name as 'region', (f.seq_region_end - f.seq_region_start) '.
      "as 'length' FROM ${feature_table}_feature f, (select distinct(seq_region_id), sr.name from seq_region".
      ' sr, coord_system cs where sr.coord_system_id=cs.coord_system_id and cs.name'.
      $sr_type_clauses{$sr_type}.' and cs.is_current is TRUE) s, feature_set fs '.
      'WHERE f.feature_set_id=fs.feature_set_id AND f.seq_region_id=s.seq_region_id AND fs.name IN ("'.
      join('","',@fset_names).'");';
    
    my $cmd = "mysql -e \"".$query."\" -quick -h$host -P$port -u$user ".
      ( $pass ? "-p$pass" : "")." $dbname >${outdir}/${name}.data.${sr_type}.txt";
    print $cmd."\n";
    system($cmd) == 0 || die('Failed to dump data');
  }
  
  if($regstats){
    my $query = "SELECT * from rf_stats";
    system("mysql -e \"".$query."\" -quick -h$host -P$port -u$user ".(($pass)? "-p$pass" : "")." $dbname > ${outdir}/regstats.txt");
  }
}

my $dir = getcwd;
#print R file with analysis
#Check which parameters one might want...
if (defined $R) {
  
  print "::Generating the R plots\n";
  open(FO,">${outdir}/${name}.R");
  print FO "require(graphics)\n";
  print FO "pdf(file=\"${outdir}/${name}.pdf\")\n";
  print FO "par(cex=0.7,cex.main=1)\n";
  
  if($regstats){ 
     print FO "require(gplots);\n";
     print FO "regstats <- read.table(\"${outdir}/regstats.txt\",row.names=1,header=TRUE);\n"; 
     print FO "textplot(regstats,halign=\"left\")\n";
  }
  
  foreach my $sr_type(@sr_types){
    
    #Load the data
    print FO "data_${sr_type} <- read.table(\"${outdir}/${name}.data.${sr_type}.txt\",header=TRUE,sep=\"\\t\")\n";
    
    #Give a little space for the legend... outside the graph... (test how much space... and the size of text in lengend)
    print FO "par(xpd=T, mar=par()\$mar+c(0,0,0,5))\n";
      
    #Print individual graphs
    print FO "for (subset in split(data_${sr_type},data_${sr_type}\$name)){\n";
    print FO "    subdata <- lapply(split(subset, subset\$region),function(x) x\$length)\n";
    print FO "    barplot(sapply(subdata, function(x) length(x)), main=subset\$name[1], xlab=\"region\",ylab=\"Number of Peaks\", col=rainbow(length(subdata)))\n";
    print FO "    legend('topright', inset=c(-".$inset_main.",0), legend=levels(subset\$region),fill=rainbow(length(subdata)), cex=0.8)\n";
    print FO "    boxplot(subdata,main=subset\$name[1],xlab='region',ylab='Peak length', col=rainbow(length(subdata))";
    if($no_outliers){ print FO ",outline=FALSE"; }
    print FO ")\n"; 
    print FO "    legend('topright',inset=c(-".$inset_main.",0),legend=levels(subset\$region), fill=rainbow(length(subdata)), cex=0.8)\n";    
    print FO "}\n";

   }

  if($compare){ 
    
    #For brevity only compare true chromosomes
    print FO "par(xpd=T, mar=par()\$mar+c(0,0,0,10))\n";

    #Global overview comparison
    print FO "data_region <-lapply(split(data_chromosome, data_chromosome\$name), function(x) x\$length)\n";
    print FO "barplot(sapply(data_region, function(x) length(x)),main='Number of Peaks per Set', xlab='Set',ylab='Number of Peaks', col=rainbow(length(data_region)), xaxt='n')\n";
    print FO "legend('topright', inset=c(-".$inset_compare.",0),legend=levels(data_chromosome\$name),fill=rainbow(length(data_region)), cex=0.8)\n";  
    #print FO "axis(1, labels=FALSE, at=1:length(data_region), tick=TRUE)\n";
    print FO "boxplot(data_region,main='Peak Length per Dataset',xlab=\"Set\",ylab=\"Peaks Length\", col=rainbow(length(data_region)), xaxt='n'";
    if($no_outliers){ print FO ",outline=FALSE"; }
    print FO ")\n";
    print FO "legend('topright',inset=c(-".$inset_compare.",0),legend=levels(data_chromosome\$name),fill=rainbow(length(data_region)), cex=0.8)\n";    
    print FO "axis(1, labels=FALSE, at=1:length(data_region), tick=TRUE)\n";

    #Print Comparative graphs by Region
    print FO "for (subset in split(data_chromosome,data_chromosome\$region)){\n";
    print FO "    subdata <- lapply(split(subset, subset\$name),function(x) x\$length)\n";
    print FO "    barplot(unlist(lapply(subdata, function(x) length(x))), main=subset\$region[1], xlab=\"region\",ylab=\"Number of Peaks\", col=rainbow(length(subdata)), xaxt='n')\n";
    print FO "    legend('topright',inset=c(-".$inset_compare.",0),legend=levels(subset\$name),fill=rainbow(length(subdata)), cex=0.8)\n";
    #print FO "    axis(1, labels=FALSE, at=1:length(subdata), tick=TRUE)\n";

    print FO "    boxplot(subdata,main=subset\$region[1],xlab='Set',ylab='Peak length', col=rainbow(length(subdata)), xaxt='n'";
    if($no_outliers){ print FO ",outline=FALSE"; }

    print FO ")\n"; 
    print FO "    legend('topright',inset=c(-".$inset_compare.",0),legend=levels(subset\$name),fill=rainbow(length(subdata)), cex=0.8)\n";    

    print FO "    axis(1, labels=FALSE, at=1:length(subdata), tick=TRUE)\n";
    print FO "}\n";
    
  }
    
   
  close FO;
  chdir($outdir);
  #This submits to the yesterday queue by default
  system "R CMD BATCH --slave ${name}.R";
}

chdir($dir);
  
  
=cut
  
   
}
 

1;