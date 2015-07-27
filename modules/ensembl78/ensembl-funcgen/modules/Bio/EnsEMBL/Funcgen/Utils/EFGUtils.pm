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

Bio::EnsEMBL::Funcgen::Utils::EFGUtils

=head1 DESCRIPTION

This module collates a variety of miscellaneous methods.


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(list your required methods here);



=cut

###############################################################################

package Bio::EnsEMBL::Funcgen::Utils::EFGUtils;

use warnings;
use strict;

use File::Find; #qw( find )
use Digest::MD5;
use Bio::EnsEMBL::Utils::Exception qw( throw      );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Data::Dumper                   qw( Dumper );
use Scalar::Util                   qw( blessed    );
use File::Path                     qw( make_path  );
use File::Basename                 qw( dirname fileparse );
use File::Spec;
use Time::Local;
use FileHandle;
use Carp qw( confess );

use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  add_DB_url_to_meta
  add_external_db
  assert_ref_do
  assert_refs
  assert_refs_do
  backup_file
  check_file
  create_Storable_clone
  convert_strand
  dump_data
  file_suffix_parse
  generate_checksum
  generate_slices_from_names
  get_alignment_file_prefix_by_ResultSet
  get_set_prefix_from_Set
  get_study_name_from_Set
  get_current_regulatory_input_names
  get_date
  get_file_format
  get_month_number
  gunzip_file
  is_bed
  is_gzipped
  is_sam
  mean
  median
  open_file
  parse_DB_url
  path_to_namespace
  run_backtick_cmd
  run_system_cmd
  scalars_to_objects
  set_attributes_by_my_scalar_names
  species_chr_num
  species_name
  split_CamelCase
  strip_param_args
  strip_param_flags
  url_from_DB_params
  validate_checksum
  validate_package_path
  validate_path
  validate_sam_header
  which_path
  write_checksum
 );




#Split out methods into FileUtils.pm?


my %strand_syns = ('+'  =>  1,
                   '-'  => -1,
                   '.'  =>  0,
                   '1'  => 1,
                   '0'  => 0,
                   '-1' => '-1');



#Handy function to add an external_db entry when needed
sub add_external_db{
  my ($efg_db, $db_name, $db_release, $db_display_name) = @_;
  my $sql = "select external_db_id from external_db where db_name='$db_name' and db_release='$db_release'";
  my ($db_id) =  $efg_db->dbc->db_handle->selectrow_array($sql);
  
  if($db_id){
    warn "External DB $db_name $db_release already exists in db with db_id $db_id\n";
  } else {
    #TODO check if it there was a failure
    $efg_db->dbc->do('insert into external_db (db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) '.
      "values('$db_name', '$db_release', 'KNOWNXREF', '1', '5', '$db_display_name', 'MISC')");
  }
  
  return;
}


=head2 assert_ref_do

  Arg [1]     : The reference to check
  Arg [2]     : The method we want to call
  Arg [3]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given object/ref is
                implements the given method. Will throw exceptions.

  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_ref_can($gene, 'dbID');
  Exceptions  : If the reference is not defined, if the object does not
                implement the given method and if no method was given to check
  Status      : Stable

=cut

#Ideally want to make assert_ref support can and do
#then we can make all these wrapper to that method passing a $do_can arg

sub assert_ref_do {
  my ($ref, $expected, $method, $attribute_name) = @_;
  #Can't do this as we expect a method return value
  #return 1 unless $ASSERTIONS;

  $attribute_name ||= '-Unknown-';
  throw('No method given') if ! defined $method;
  throw('No expected type given') if ! defined $expected;
  throw "The given reference $attribute_name is not defined" unless defined $ref;
  throw "The given reference $attribute_name is not blessed" unless blessed($ref);
  throw("${attribute_name}'s type '${ref}' is not an ISA of '${expected}'") if ! $ref->isa($expected);

  if(! $ref->can($method)) {
    my $str_ref = ref($ref);
    throw sprintf(q{The given blessed reference '%s' for attribute '%s' does not implement the method '%s'}, $str_ref, $attribute_name, $method);
  }


  return $ref->$method;
}


#Return type isn't actually boolean as it is never 0

sub assert_refs {
  my ($refs, $expected, $attribute_name) = @_;

  my $array_attr_name = (defined $attribute_name) ?
    $attribute_name.' array' : undef;

  assert_ref($refs, 'ARRAY', $array_attr_name);
  map { assert_ref($_, $expected, $attribute_name) } @$refs;

  return 1;
}


=head2 assert_refs_do

  Arg [1]     : An Arrayref of references to check
  Arg [2]     : The type we expect
  Arg [3]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given objects/refs are
                what you expect. This behaves in an identical manner as
                C<assert_ref> does works on an array ref of references

                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_refs([[],[],[]], 'ARRAY');
  Exceptions  : Throws is references argument is not an Arrayref, also
                if the expected type was not set and if the given reference
                was not assignable to the expected value.
  Status      : Stable

=cut


sub assert_refs_do {
  my ($refs, $expected, $method, $attribute_name) = @_;

  my $array_attr_name = (defined $attribute_name) ?
    $attribute_name.' array' : undef;

  assert_ref($refs, 'ARRAY', $array_attr_name);
  my @values;
  map { push @values, assert_ref_do($_, $expected, $method, $attribute_name) } @$refs;

  return \@values;
}


#can't do check_ref_do as this would make the return type ambiguous


sub backup_file{
  my $file_path = shift;

  throw("Must define a file path to backup") if(! $file_path);

  if (-f $file_path) {
    #$self->log("Backing up:\t$file_path");
    my $date = run_backtick_cmd('date \'+%T\'');
    system ("mv ${file_path} ${file_path}.${date}") == 0 || return 0;
  }

  return 1;
}


=head2 create_Storable_clone

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable to clone
  Arg [2]    : Hash containing list of object parameters linked to
               the clone and to be reset
  Example    :
  Description: Blesses an object and replaces all linked objects with the ones
               passed with the $params_hash
  Returntype : cloned object
  Exceptions : Throws if object is not a Storable
  Caller     : general
  Status     : At risk - not tested

=cut

#Currently not exposing the no db reset flag from reset_relational_attributes
#as we would never want to do this here

sub create_Storable_clone {
  my ($obj, $params_hash) = @_;

  assert_ref($obj, 'Bio::EnsEMBL::Funcgen::Storable', 'object');

  if(! $obj->can('reset_relational_attributes')){
     throw('You must pass a Bio::EnsEMBL::Funcgen::Storable which can '.
      'call the reset_relational_attributes method');
  }

  my $clone = bless({%{$obj}}, ref($obj));
  $clone->reset_relational_attributes($params_hash);
  return $clone;
}

sub dump_data {
  my $data   = shift;
  my $indent = shift;
  my $terse  = shift;
  
  if((defined $indent) and $indent !~ /[0-3]/o){
    throw("Indent must be 0,1,2,3 not $indent");
  }
  
  $indent = 2 if ! defined $indent;

  my $dumper = Data::Dumper->new([$data]);
  $dumper->Indent($indent);
  $dumper->Terse($terse) if $terse;
  
  return $dumper->Dump;
}


#Do this as File::Basename doesn't quite get there.

sub file_suffix_parse{
  my $filepath = shift;
  return (fileparse($filepath) =~ /(.*)\.([^\.]+)$/);
}


#Generates slices from names or optionally alll default top level nonref
#slice ref args can be array ref (inc empty) or undef

sub generate_slices_from_names{
  my $slice_adaptor = shift;
  my $slice_names   = shift;
  my $skip_slices   = shift;
  my $level         = shift;
  my $non_ref       = shift;
  my $inc_dups      = shift;
  my $assembly      = shift;
  my (@slices, $slice, $sr_name, $have_slice_names, $have_skip_slices);

  #Test if $assembly is old?

  #Validate array ref skip/slice name args

  if(defined $slice_names){
    if(ref($slice_names) ne 'ARRAY'){
      throw('Slice names argument must be an ARRAYREF');
    }

    $have_slice_names = 1 if @$slice_names;
  }

  if(defined $skip_slices){
    if(ref($skip_slices) ne 'ARRAY'){
      throw('Skip slices argument must be an ARRAYREF');
    }

    $have_skip_slices = 1 if @$skip_slices;
  }

  #Generate slices

  if($have_slice_names){

    foreach my $name(@$slice_names){
      $slice = $slice_adaptor->fetch_by_region(undef, $name, undef, undef, undef, $assembly);
      #Why is this failing for hap regions?

      if(! $slice){

        if(! (eval { $slice = $slice_adaptor->fetch_by_name($name); 1 } &&
              defined $slice) ){
          #Need to eval this as it will break with incorrect formating
          throw("Could not fetch slice by region or name:\t".$name);
        }
      }

      $sr_name = $slice->seq_region_name;

      next if(grep { /^${sr_name}$/ } @$skip_slices);
      push @slices, $slice;
    }
  }
  elsif($level){

    #Can't guarantee what level will have an assembly version
    #Can only throw if we have set toplevel and assembly
    #as this will alway default to the current assembly
    #if($assembly &&
    #   ($level eq 'toplevel') ){
    #  throw('You cannot specify an assembly version with the toplevel coordinate system');
    #}
    #Let core API handle this?


    my @tmp_slices = @{$slice_adaptor->fetch_all($level, $assembly, $non_ref, $inc_dups)};

    if($have_skip_slices){

      foreach my $slice(@tmp_slices){
        $sr_name = $slice->seq_region_name;
        push @slices, $slice if ! grep { /^${sr_name}$/ } @$skip_slices;
      }
    }
    else{
      @slices = @tmp_slices;
    }
  }
  else{
    throw('You must either pass an arrayref of slice names or specify the toplevel flag');
  }


  if(! @slices){
    throw("You have specified slice_names and skip_slices paramters which have generated no slices.\nslice_names:\t".join(' ',@$slice_names)."\nskip_slices:\t".join(' ', @$skip_slices));
  }

  return \@slices;
}


# Tracking DB methods
# Move to DBAdaptor? Can we add this as a separate package in the same module?

sub get_current_regulatory_input_names{
  my ($tdb, $efg_db, $focus) = @_;

  warn "Move get_current_regulatory_input_names to the TrackingAdaptor";

  #Validate is production?
  my $sql;



  if($focus){
    $focus = 'Focus';
    $sql   = 'SELECT efgdb_set_name from dataset where is_focus=true and is_current=true and species="'.$efg_db->species.'"';
  }
  else{
    $focus = 'Non-focus';
    #0 rather than false so we don't get NULLs
    $sql = 'SELECT efgdb_set_name from dataset where is_focus=0 and is_current=true and species="'.$efg_db->species.'"';
  }



  #Currently efgdb_set_name can either be data_set or feature_set name!
  #Need to standardise this

  my @prd_names = @{$tdb->db_handle->selectcol_arrayref($sql)};
  my @names;
  my @failed_sets;

  foreach my $prd_name(@prd_names){

    $sql = "SELECT name from feature_set where name like '${prd_name}%'";
    my @tmp_names =  @{$efg_db->dbc->db_handle->selectcol_arrayref($sql)};

    #This is causing problems with multiple feature sets with differing analyses

    #Do this via InputSets(using query extension?) instead of using like?

    #This is very hacky right now to get it to work
    #Need to standardise and review tracking db data.

    if(scalar(@tmp_names) > 1){

      $sql = "SELECT name from feature_set where name ='${prd_name}_ccat_histone'";
      @tmp_names =  @{$efg_db->dbc->db_handle->selectcol_arrayref($sql)};

      if(scalar(@tmp_names) == 1){
        push @names, $tmp_names[0];
      }else{
        push @failed_sets, $prd_name;
      }
    }
    elsif(scalar(@tmp_names) == 0){
      push @failed_sets, $prd_name;
    }
    else{
      push @names, $tmp_names[0];
    }
  }

  if(@failed_sets){
    throw("Failed to find unique $focus FeatureSets for production dataset names:\n\t".
          join("\n\t", @failed_sets)."\n");
  }

  return @names;
}


sub get_date{
	my ($format, $file) = @_;

	my ($time, $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst);


	throw("File does not exist or is not a regular file:\t$file") if $file && ! -f $file;


	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = (defined $file) ?
	  localtime((stat($file))[9]) : localtime();

	#print "	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst)\n";

	if((! defined $format && ! defined $file) || $format eq "date"){
		$time = ($year+1900)."-".$mday."-".($mon+1);
	}
	elsif($format eq "time"){#not working!
		$time = "${hour}:${min}:${sec}";
	}
	elsif($format eq "timedate"){#
	  $time = localtime();
	}
	else{#add mysql formats here, datetime etc...
		croak("get_date does not handle format:\t$format");
	}

	return $time;
}


sub write_checksum{
  my $file_path = shift;
  my $params    = shift || {};
  assert_ref($params, 'HASH');
  
  my $signature_file = (exists $params->{signature_file}) ? $params->{signature_file} : undef;
  my $digest_method  = (exists $params->{digest_method})  ? $params->{digest_method}  : undef;
  my $debug          = (exists $params->{debug})          ? $params->{debug}          : 0;
  $digest_method   ||= 'hexdigest';
  my $md5_sig        = generate_checksum($file_path, $digest_method);
  my $file_name      = fileparse($file_path);

  if(! defined $signature_file){
    $signature_file = $file_path.'.CHECKSUM';
    warn "Defining default signature file as $signature_file\n";
  }
  
  warn "Writing checksum:\t$md5_sig\t$file_name\nTo signature file:\t$signature_file\n" if $debug;
  my $checksum_row = $md5_sig."\t".$file_name."\t".$digest_method;

  #Update or create entry in signature_file
  my $sigfile_row;
  
  if(-f $signature_file){
    warn "Signature file is $signature_file\n" if $debug;
    
    $sigfile_row = run_backtick_cmd("grep '$file_name' $signature_file", 1);
    #no exit flag as grep will return 1 if no results are returned
    
    if($?){
      warn("$sigfile_row\nFailed to grep old checksum from existing signature file:\t$signature_file");
      $sigfile_row = '';
    }
  }

  if($sigfile_row){ #Update entry
    my $sed_string = "sed -r 's/^.*[[:space:]]$file_name";

    #handle absent digest method here
    if($sigfile_row =~ /${digest_method}$/){
      $sed_string .= "[[:space:]]$digest_method";
    }

    run_system_cmd("${sed_string}\$/$checksum_row/' $signature_file > ${signature_file}.tmp; ".
                   "mv ${signature_file}.tmp $signature_file");
  }
  else{#create/append to file
    run_system_cmd("echo '$checksum_row' >> $signature_file");
  }

  return $md5_sig; #or return sig file?
}


sub generate_checksum{
  my $file          = shift;
  my $digest_method = shift;
  #if($file =~ /\.gz$/){
    #warn("It is unsafe to generate checksums for a compressed file:\n\t$file");  
  #}
  
  $digest_method ||= 'hexdigest';

  my $ctx = Digest::MD5->new;

  if(! $ctx->can($digest_method)){
    throw("Digest::MD5 cannot call method $digest_method, ".
     'please choose a valid digest method or omit for default hexdigest method');
  }

  #Don't use bareword (FILE) for descriptor as is stored in symbol table for this package
  #meaning potential interference if FILE is used elsewhere. (PBP 203)
  open(my $CHK_FILE, '<', $file) or throw("Cannot open file for md5 digest:\t$file\n$!");
  binmode $CHK_FILE;
  $ctx->addfile($CHK_FILE);#eval this?
  my $md5_sig = $ctx->$digest_method;
  close($CHK_FILE);

  return $md5_sig;
}

#assume the format of the file is:
#checksum_sig\tfilename\tdigestmethod

sub validate_checksum{
  my $file_path = shift;
  my $params    = shift || {};
  assert_ref($params, 'HASH');
  my $signature_file = (exists $params->{signature_file})    ? $params->{signature_file}    : undef;
  my $digest_method  = (exists $params->{digest_method})     ? $params->{digest_method}     : undef;
  my $md5_sig        = (exists $params->{checksum})          ? $params->{checksum}          : undef; 
  my $md5_optional   = (exists $params->{checksum_optional}) ? $params->{checksum_optional} : undef;
  my $debug          = (exists $params->{debug})             ? $params->{debug}             : 0;
  
  if((defined $signature_file) && (defined $md5_sig)){
    throw("Params 'signature_file' and 'checksum' are mutually exclusive, please restrict to one");
  }

  my $file_name = fileparse($file_path);

  if(! defined $md5_sig){
    
    if(! defined $signature_file){ 
      $signature_file = $file_path.'.CHECKSUM';
    }
   
    if((! -f $signature_file) && 
       ($signature_file =~ /\.bam.CHECKSUM/o)){
      warn "!!!!! REMOVE THIS !!!!\nTEMPORARILY CREATING missing MD5 file:\t$signature_file";
      write_checksum($file_path);         
    }
      
    if(-f $signature_file){
      my $qtd_file_name = quotemeta($file_name);
      my $cmd = "grep -E '[[:space:]]+$qtd_file_name\[[:space:]]*.*\$' $signature_file";
      #warn "Checksum file grep:\t$cmd" if $debug;
      my $checksum_row = run_backtick_cmd($cmd);
      #warn "Checksum file row:\t$checksum_row" if $debug;
      my ($sig_file_name, $sig_digest_method);
      ($md5_sig, $sig_file_name, $sig_digest_method) = split(/\s+/, $checksum_row);
  
      if((! defined $sig_file_name) ||
         ($sig_file_name ne $file_name)){
        throw("Failed to find $file_name entry in checksum signature file:".
          "\n\t$signature_file\nUsing:\t$cmd");
      }
    
      #default to digest method in file
      $digest_method     ||= $sig_digest_method;
      #Also need to account for absent $sig_digest
      $sig_digest_method ||= $digest_method;
  
    
      if(defined $sig_digest_method){
        if($digest_method ne $sig_digest_method){
        throw("Specified digest method($digest_method) does not match method found in ".
          "checksum signature file($sig_digest_method):\n\t$file_path\n\t$signature_file");
        }
      }
      else{
        warn "Could not find digest method in checksum signature file, assuming $digest_method";
      }
    }
  }

  if(defined $md5_sig){
    my $new_md5_sig = generate_checksum($file_path, $digest_method);
    
    if($md5_sig ne $new_md5_sig){
      #This could be due to mismatched digest methods?
      my $msg = 'MD5 checksums do not match:'.
        "\n\tExpected\t$md5_sig\n\tFound\t\t$new_md5_sig".
        "\n\tFile:\t\t$file_path";
      $msg .= "\n\tSignature file:\t$signature_file" if $signature_file;
      throw($msg);
    }
  }
  elsif(! $md5_optional){
    throw("Failed to find checksum file:\t$signature_file\n".
          'Please specify one as an argument, create a default file or pass the checksum_optional param');
  }

  return $md5_sig;
}

#todo purge_checksum from file, as we are deleting the file, or replacing it with a gz file


#could try and get the file suffix first to be a bit more directed here
#especially if the format list grows

sub get_file_format{
  my $file = shift;

  my $format = &is_bed($file);

  if(! $format){
    $format =  &is_sam($file);

    #Add more tests here
  }


  return $format;
}


sub get_month_number{
  my($mon) = @_;
  my %month_nos =(
          "jan", "01",
          "feb", "02",
          "mar", "03",
          "apr", "04",
          "may", "05",
          "jun", "06",
          "jul", "07",
          "aug", "08",
          "sep", "09",
          "oct", "10",
          "nov", "11",
          "dec", "12",
         );
  return (exists $month_nos{lc($mon)}) ? $month_nos{lc($mon)} : undef;
}


sub gunzip_file {
  my $filepath = shift;
  my $was_gzipped = 0;
 

  if( is_gzipped($filepath) ){
    run_system_cmd("gunzip $filepath");
    $filepath =~ s/\.(?:t){0,1}gz$//;
    $was_gzipped = 1;
  }

  return($filepath, $was_gzipped);
}


#TODO Change this to use named pipes for the zcat so pipe exit can be caught properly
#do this in open file, by passing an array of open modes

sub is_bed {
  my $file = shift;
  my ($BED_FILE, @line);

  if(&is_gzipped($file, 1)){
    open($BED_FILE, "zcat $file 2>&1 |") or throw("Can't open file via zcat:\t$file");
  }
  else{
    open($BED_FILE, '<', $file) or throw("Can't open file:\t$file");
  }

  while (<$BED_FILE>) {
    chomp;
    @line = split("\t", $_);
    last;
  }
  
  close($BED_FILE);

  if (scalar @line < 6) {
    #warn "Infile '$file' does not have 6 or more columns. We expect bed format:\t".
    #  "CHROM START END NAME SCORE STRAND.\n";
    return 0;
    #} elsif ($line[0] !~ m/^((chr)?[MTXYNT_\d]+)$/) {
    #    warn ("1st column must contain name of seq_region (e.g. chr1 or 1) in '$file'");
    #    return 0;
    #Commented this out for now due to HSCHR_RANDOM seqs
    #How does the webcode handle this?
  }
  elsif ($line[1] !~ m/^\d+$/ && $line[2] =~ m/^\d+$/) {
    warn "2nd and 3rd column must contain start and end respectively in '$file'\n";
    return 0;
  }
  elsif ($line[5] !~ m/^[+-\.]$/) {
    warn "6th column must define strand (either +, - or .) in '$file'\n";
    return 0;
  }

  return 'bed';
}


#Have individual format methods, as this is likely to hit performance
#would be nive to have a no thro mode, which would return the original value

sub convert_strand {
  my $strand = shift;

  if(! defined $strand){
   throw('Strand argument is not defined');
  }

  if(! exists $strand_syns{$strand}){
    throw('Strand argument is not a supported strand synonym, should be one of: '.
      join(' ', keys %strand_syns));
  }

  return $strand_syns{$strand};
}


#This can return true for some file formats which use gzip intenrally
#but are not valid gzip files e.g. bam
#Hence we now test for known gzip suffixes
#Probably need an override flag here
#changed return type to file suffix, for file name subbing

sub is_gzipped {
  my $file               = shift; 
  my $fail_if_compressed = shift;
  my $gzip               = 0;

  throw ("File does not exist:\t$file") if ! -e $file;

  open(my $GZ_FILE, "file -L $file |")
    or throw("Can't execute command 'file' on '$file'");
  my $file_info = <$GZ_FILE>;
  close($GZ_FILE);


  if($file =~ /\.(?:t){0,1}gz$/){

    if($file_info =~ /compressed data/){

      if($file_info =~ /gzip/){
        $gzip = 1;
      }
      else{
        throw("File is compressed but not with gzip, please unzip or gzip:\t$file_info");
      }
    }
  }

  return $gzip;
}



#Change these to also return the gz status

sub is_sam{
  my $file = shift;

  #warn "Only checking file suffix for is_sam";
  #Could check for header here altho this is not mandatory!
  #Can we use web format guessing code?

  my $gz = (&is_gzipped($file, 1)) ? '.gz' : '';

  return ($file =~ /.sam${gz}/) ? 'sam' : 0;
}

#need is bam here too!


sub mean{
  my $scores = shift;
  my $total = 0;
  map { $total+= $_ } @$scores;
  return $total / scalar(@$scores);
}


#todo, add no sort flag, and sort by default?
#will need to deref to avoidaffecting the the array in the caller

sub median{
  my $scores = shift;
  return if ! @$scores;

  my ($median);
  my $count = scalar(@$scores);
  my $index = $count-1;
  #need to deal with undefs
  #deal with one score fastest
  return  $scores->[0] if ($count == 1);

  #taken from Statistics::Descriptive
  #remember we're dealing with size starting with 1 but indices starting at 0

  if ($count % 2) { #odd number of scores
    $median = $scores->[($index+1)/2];
  }
  else { #even, get mean of flanks
    $median = ($scores->[($index)/2] + $scores->[($index/2)+1] ) / 2;
  }

  return $median;
}


=head2 path_to_namespace

  Arg [1]    : String - Module path
  Example    : my $path = '/my/module.pm';
               eval { require $path };

               if(! $@){
                 my $namespace = path_to_namespace($path);
                 my $obj       = $namespace->new();
               }

  Description: Simply converts a path to a namespace. Useful when you have used
               required a module dynamically rather than 'using' at compile time.
  Returntype : String
  Exceptions : Throws if arguments not defined
  Caller     : General
  Status     : At risk

=cut

sub path_to_namespace {
  my $path = shift;

  if(! defined $path){
    throw('No path argument defined');
  }

  $path =~ s/[\/]+/::/g;
  $path =~ s/\.pmc*$//;
  return $path;
}


#todo test file exists and is readable/writable first to avoid zombie?
#todo Use IPC::Open open2/3 for piping? Take list of operators (and list of files e.g. in/out)

sub open_file{
  my $file             = shift;
  my $operator         = shift || '<';
  my $file_permissions = shift;

  if(! defined $file){
    throw('Must provide a file argument');
  }

  #Get dir here and create if not exists
  my $dir              = dirname($file);
  my $mkpath_opts      = {verbose => 1};
  $mkpath_opts->{mode} = $file_permissions if defined $file_permissions;

  if((! -d $dir) && 
     ($operator eq '>')){

    if(! eval { make_path($dir, $mkpath_opts); 1 }){
      throw("Failed to make_path:\t$dir\n$@");
    }
  }

  if ($operator !~ /%/) {
    $operator = "$operator $file";
  } else { #We have some piping to do
    $operator = sprintf($operator, $file);
  }

  my $fh = FileHandle->new($operator);

  #This does not catch errors when piping then redirecting
  #as the progam will have opened sucessfully
  #i.e. gzip -c | to file we don't have permission to create!
  #This simply carries on with the reset of the script
  #won't be caught until we close the file descriptor
  #Which is a little late, and we don't catch error codes on close yet

  #it maybe better to use PerlIO libs here e.g. PerlIO::gzip
  confess("Failed to open $operator") if ! defined $fh;

  #Have to chmod here as umask will over-ride permissions passed to FileHandle
  if (defined $file_permissions) {

    #Catch non-numeric here as chmod still returns true
    if ($file_permissions =~ /[^0-9]/) {
      confess("Failed to change $file permissions using:\t$file_permissions");
    }

    #chmod requires a literal octal number e.g. 0775 not '0775'
    #should catch numbers as strings here, but perl makes this very hard to test
    #Can't even system this as if we build the cmd line with an octal it will be converted to a decimal
    #There is still no way of testing for non-octal number or string
    #eval/sprintf will also not fail if there are non-octal digits i.e. 1999

    #eval will treat octal number and string as true octal number
    #else will pass non-octal string/number which we can't catch
    chmod(eval{$file_permissions}, $file);
  }

  return $fh;
}



sub parse_DB_url {
  my $url = shift;

  my ($dbtype, $user, $pass, $host, $dbname, $port);

  if($url =~ /(.*):\/\/(.*):(.*)@(.*):(.*)\/(.*)/){
    ($dbtype, $user, $pass, $host, $port, $dbname) = ($1, $2, $3, $4, $5, $6);

    if(! ($user && $host && $dbname)){
      die("The DB url must contain at least user, host and dbname elements i.e.\n\t".
        'DBTYPE://USER:[PASS]@DBHOST:[PORT]/DBNAME');
    }
  }
  else{
    throw("Unrecognized DB url format:\t$url\nShould be:\n\t".
      'DBTYPE://USER:[PASS]@DBHOST:[PORT]/DBNAME');
  }

  return [$user, $pass, $host, $dbname, $port, $dbtype];
}



################################################################################

=head2 run_system_cmd

 Description : Method to control the execution of the standard system() command

 ReturnType  : none

 Example     : run_system_cmd($my_cmd);

 Exceptions  : Throws exception if system command returns none zero and no exit boolean not set

=cut

################################################################################

#Allow no_exit as some programs(tab2mage) give successful non-zero exit codes!
#Would be nice to catch warn output also, but system doesn not handle this well and would
#to redirect STDERR to a file to read back in. 


#These are currently not capturing STDERR and so this is not reported in the
#throw. It is necessary to inspect the STDERR file to see these, which is a slight
#problem with the hive which will only list the API throw string, not the actual
#external error

#Need to use IPC::Open3
#http://learn.perl.org/faq/perlfaq8.html#How-can-I-capture-STDERR-from-an-external-command-

#Should reset $? here too? In case we want to handle that in the caller
#Although this is returned after bit shifting

sub run_system_cmd{
  my ($command, $no_exit) = @_;
  #$exit_status is same as $? i.e. $CHILD_ERROR
  #my $exit_status = system($command);#don't nest below, as this may cause problems with $!
  system($command);
  return _handle_exit_status($command, $?, $!, $no_exit);
}


################################################################################

=head2 run_backtick_cmd

 Description : Method to control the execution of the standard `backtick` command
               and return the results

 ReturnType  : Array or Scalar as required

 Example     : $Helper->debug(2,"dir=$dir file=$file");

 Exceptions  : throws exception if system command returns none zero

=cut

################################################################################


#todo Use IPC::Open3 instead of backtick here
#Although we are handling exit codes, this will allow us
#to capture STDOUT and STDERR sepatately
#http://search.cpan.org/dist/Perl-Critic/lib/Perl/Critic/Policy/InputOutput/ProhibitBacktickOperators.pm

#remove $@ from here, this is only for perl evals?
#remove $! from here as this is only set by system or library calls?

sub run_backtick_cmd{
  my $command = shift;
  my $no_exit = shift;
  my (@results, $result);
  
  if(wantarray){
    #perlcritic(3) suggests this 3 statement map should be subbed out
    #@results = map {my $line = $_; chomp $line; $line} `$command`;
    #This is substantially faster
    @results = `$command`;
  }
  else{
    $result  = `$command`;  
  }
  
  my $exit_status = $?; #$CHILD_ERROR
  #my $errno       = $!; 
  _handle_exit_status($command, $exit_status, $!, $no_exit); 
  $? = $exit_status if $no_exit;
  #Non-local 'Magic' $? assignment is a perlcritic severity 4 warning. 
  #But we actually want this 'global' behaviour 
   
  if(wantarray){
    chomp(@results); 
    return @results ; 
  }
  else{
    chomp($result);
    return $result;
  }
}


#This is currently not capturing error output very well.
#If we are running a script, the error message from the script 
#(and whatever binary it is running) is not captured
#Meaning Running scripts from script can mean errors go missing.
#use IPC::Run instead of system? But this is not a core module!
#IPC::Cmd is available in 5.12 but not 5.10 Grr!

sub _handle_exit_status{
  my ($cmd, $exit_status, $errno, $no_exit) = @_;
  my ($exit_code, $err_string);
  
  if ($cmd =~ /\|/){
    warn "Failed piped commands may not be caught:\n$cmd\n";    
  }

  $exit_code = $exit_status >> 8; #get the true exit code
  
  if ($exit_status == -1) {
    $err_string = "Failed to execute:\t$cmd\nError:\t$errno\n";
  }
  elsif ($exit_status & 127) {
    $err_string = sprintf("Child process died with signal %d, %s coredump\n$cmd\nError:\t$errno\n",
                    ($exit_status & 127),
                    ($exit_status & 128) ? 'with' : 'without');
  }
  elsif($exit_status != 0) {
    #$cmd may contain sprintf symbols here we need to escape
    $err_string = sprintf("Child process exited with value %d:", $exit_code)."\t$cmd\nError:\t$errno\n";
  }
  
  if(defined $err_string){
    if(! $no_exit){
      throw($err_string);
    }   
    else{
      warn $err_string;  
    }
  }
  
  return $exit_code;
}


=head2 scalars_to_objects

  Arg [1]    : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Arg [2]    : String - class name of object to retrieve
  Arg [3]    : String - method name to use
  Arg [4]    : Arrayref - Scalar arguments to use iteratively with the fetch method.
               A single scalar can also be passed in a scalar context.
  Example    : my @cell_types = @{scalars_to_object($db, 'CellType',
                                                    'fetch_by_name',
                                                    [ qw ( GM06990 HUVEC H1ESC ) ])};
  Description: Convenience method to fetch objects iteratively based on an
               array of scalar arguments. Useful for processing command line arguments.
  Returntype : Arrayref of Objects
  Exceptions : Throws if arguments not specified or valid
               Throws if cannot get relevant adaptor based on the class name
               Throws if Object Adaptor cannot call fetch method
               Throws if Object not retrieved
  Caller     : General
  Status     : At risk

=cut

#We could be clever here and support incorporation of this sub
#into a class by assuming that $db might also be $self if called as a method.
# i.e. could test if can('db') and use that?

sub scalars_to_objects{
  my ($db, $class_name, $fetch_method, $scalars) = @_;
  assert_ref($db, 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor', 'db');

  #Be kind and handle single scalars in scalar context
  if(defined $scalars && 
     (! ref($scalars)) ){
    $scalars = [$scalars];     
  }

  if(! ((defined $scalars && (ref($scalars) eq 'ARRAY')) &&
         defined $class_name &&
         defined $fetch_method) ){
    throw('You need to specify a class name and a fetch method and an Arrayref of scalar method arguments');
  }

  my @objs;
  my $adaptor_method = 'get_'.$class_name.'Adaptor';

  my $adaptor = $db->$adaptor_method;

  if(! ((defined $adaptor) &&
        $adaptor->can($fetch_method))){
    throw("Could not $adaptor_method or ${class_name}Adaptor cannot call method $fetch_method");
  }

  foreach my $str(@$scalars){

    my $obj = $adaptor->$fetch_method($str);

    if(! defined $obj){
      throw("Could not fetch object using ${class_name}Adaptor->${fetch_method}('$str')");
    }

    push @objs, $obj;
  }

  return \@objs;
}


sub species_chr_num{
    my ($species, $val) = @_;

    ($species = lc($species)) =~ s/ /_/;

    my %species_chrs = (
                        homo_sapiens => {(
                                          'x' => 23,
                                          'y' => 24,
                                          'mt' => 25,
                                         )},

                        mus_musculus => {(
                                          'x'  => 20,
                                          'y'  => 21,
                                          'mt' => 22,
                                           )},

                        rattus_norvegicus =>  {(
                                                'x'  => 21,
                                                'y'  => 22,
                                                'mt' => 23,
                                               )},
                       );

    die("species not defined in chromosome hash") if(! exists $species_chrs{$species});

    return (exists $species_chrs{$species}{lc($val)}) ? $species_chrs{$species}{lc($val)} : $val;
}



#migrate this data to defs file!!??
#must contain all E! species and any other species which are used in local DB extractions
#NEED TO ADD FLY!!

sub species_name{
  my($species) = @_;
  my %species_names = (
		       "HOMO_SAPIENS", "human",
		       "MUS_MUSCULUS", "mouse",
		       "RATTUS_NORVEGICUS", "rat",
		       "CANIS_FAMILIARIS", "dog",
		       "PAN_TROGOLODYTES", "chimp",
		       "GALLUS_GALLUS", "chicken",
		       "SACCHAROMYCES_CEREVISIAE", "yeast",
		       "HUMAN",  "HOMO_SAPIENS",
		       "MOUSE", "MUS_MUSCULUS",
		       "RAT","RATTUS_NORVEGICUS",
		       "DOG", "CANIS_FAMILIARIS",
		       "CHIMP", "PAN_TROGOLODYTES",
		       "CHICKEN", "GALLUS_GALLUS",
		       "YEAST", "SACCHAROMYCES_CEREVISIAE",
		      );

  return $species_names{uc($species)};
}

#returns array

sub split_CamelCase{
  my $string = shift || die('Must provide a CamelCase string to split');
  #I'd like to thank google and salva  
  #This will simply drop any non-letter characters
  return $string =~ /[A-Z](?:[A-Z]+|[a-z]*)(?=$|[A-Z])/g;
}


#These subs are useful for implementing
#a farm mode in a run script, where a script can
#submit itself to the farm as slice based jobs

#strip cmd line params and associated arguments from a list
#should not be used to remove flag options i.e. no following args
#as this may cause removal of any following @ARGV;
#Can this be used on flattened args hash?

sub strip_param_args{
  my ($args, @strip_params) = @_;

  my $param_name;
  my $seen_opt = 0;

  foreach my $i(0..$#{$args}){

	if($args->[$i] =~ /^[-]+/){
	  $seen_opt = 0;#Reset seen opt if we seen a new one

	  ($param_name = $args->[$i]) =~ s/^[-]+//;

	  if( grep { /^${param_name}$/ } @strip_params){
		  $seen_opt = 1;
	  }
	}

	#$args->[$i] = '' if $args->[$i] =~ /^[-]+farm/;#Only remove current flag
	#$seen_opt = 1 if $args->[$i] =~ /^[-]+skip_slices/;
	#$seen_opt = 1 if $args->[$i] =~ /^[-]+slice/;#Don't have full param name incase we have just specified -slice

	$args->[$i] = '' if $seen_opt;#Remove option and args following option
  }

  return $args;
}


sub strip_param_flags{
  my ($args, @strip_params) = @_;
  my @args;

  foreach my $flag(@strip_params){
	  push @args, grep { !/[-]+${flag}$/ } @$args;
  }

  return \@args;
}


sub add_DB_url_to_meta {
  my $url_type = shift;
  my $url      = shift;
  my $db       = shift;
  
  if(! defined $url_type){
    throw('Must provide URL_TYPE, URL and DBAdaptor arguments');
  }
  parse_DB_url($url);
  assert_ref($db, 'Bio::EnsEMBL::DBSQL::DBAdaptor', 'db');#restrict this to Funcgen?

  my $mc         = $db->get_MetaContainer;
  my $meta_key = $url_type.'_url';
  my $meta_value = $mc->single_value_by_key($meta_key);

  if( ! defined $meta_value){ #Store the new meta entry
    #Add key via API to store with appropriate species_id

    if(! eval { $mc->store_key_value($meta_key, $url); 1}){
      throw("Failed to store hive url meta entry:\t$meta_key\t$url\n$@");
    }
  }
  elsif($meta_value ne $url){
    throw("Could not set $meta_key meta entry:\t".$url."\nAs ".$db->dbc->dbname.
      " is currently locked to a different hive DB:\t$meta_value\n".
      'Please use that hive DB or drop that pipeline, or remove the meta entry');
  }
  #else meta values match and all is good.

  return;
}



sub url_from_DB_params {
  my $db_params = shift;

  my $host   = (exists $db_params->{'-host'})   ? $db_params->{'-host'} :
    throw('Missing mandatory -host paramter in DB parameters hash');
  my $user   = (exists $db_params->{'-user'})   ? $db_params->{'-user'} :
    throw('Missing mandatory -user paramter in DB parameters hash');
  my $dbname = (exists $db_params->{'-dbname'})   ? $db_params->{'-dbname'} :
    throw('Missing mandatory -dbname paramter in DB parameters hash');
  my $driver = (exists $db_params->{'-driver'}) ? $db_params->{'-driver'} : 'mysql';
  my $port   = (exists $db_params->{'-port'}) ? $db_params->{'-port'} : '';
  my $pass   = (exists $db_params->{'-pass'}) ? $db_params->{'-pass'} : '';

  return $driver.'://'.$user.':'.$pass.'@'.$host.':'.$port.'/'.$dbname;
}

#Very simple method to check for a file and it's compressed variants
#This may cause problems if the $file_path is already .gz

#todo tidy up suffix handling, will it ever be anything other than .gz?


# To run on checksum validation simply define the checksum param
#to the expected checksum or undef if we want to read from a sig file

sub check_file{
  my $file_path = shift;
  my $suffix    = shift;
  my $params    = shift;
  my $found_path;

  if(-f $file_path){
    $found_path = $file_path;
  }
  elsif(defined $suffix){
    
    if($file_path =~ /\.${suffix}/o){ #check without first
      #Just incase it was already in the $file_path
      $file_path =~ s/\.${suffix}$//;
      
      if(-f $file_path){
        $found_path = $file_path;  
      }
    }
    elsif(-f $file_path.'.'.$suffix){
      $found_path =  $file_path.'.'.$suffix
    }
  }

  if($found_path){

    if($params->{gunzip} && 
       ($found_path =~ /\.gz$/o)){          
      $found_path = gunzip_file($found_path);  
    }

    if(defined $params){
      assert_ref($params, 'HASH');
      
      if(exists $params->{checksum}){
        validate_checksum($found_path, $params);  
      }
    }
  }

  return $found_path;
}


sub validate_path{
  my ($path, $create, $is_dir, $alt_root) = @_;

  #how will alt_root and create interact
  #should they be mutually exclusive
  #or should we try both, but only create on the first root?
  #alt root will only work if path is a ref!
  #Is this the right way to do this,
  #or should we handle root subbing in hive where we have access to both roots?

  #We should always test the dir first?
  #i.e. we don't want to mix files from two different dirs?

  #How will thi swork when we have set teh alignment dir
  #but we can't find something in there
  #we would also need to set alt alignment dir?
  #Whenever we set a dir, we shoudl also test and set an alt dir
  #if the first dir is absent, but we find the alt dir, then we should
  #set the alt as the main dir, then undef the alt dir

  #when ever we validate a path, we will then have access to an alt dir for
  #which ever root we are using, or not if it doesn't exist

  #Still a little tricky about using files from two different dirs
  #we definately don't want to do this!
  #hard to manage this
  #maybe we just want to do this for fastqs and other 'reference' files etc?
  #all other intermediate piepline outptu/input should use the same root!


  if(! defined $path){
    throw('A path argument must be defined');
  }

  #$label = ($label) ? ' '.$label : '';
  my $path_ref = ref($path);

  ### Build path
  if($path_ref ne ''){

    if($path_ref ne 'ARRAY'){
      throw("Path argument is not valid:\t$path\n".
        'It must be a String or an Arrayref of string to pass to File::Spec');
    }

    #This works for files too instead of catfile
    #Will not prepend / if absent
    $path = File::Spec->catdir(@$path);
  }
  else{
    #Tidy up paths from config
    #This won't tidy up paths already stored in meta from init_pipeline.pl
    $path = File::Spec->canonpath($path);
  }

  ### Validate/create file/dir
  if(-e $path){

    if($is_dir){
      if( ! -d $path){
        throw("Found path but is not a directory:\n\t$path");
      }
      elsif($create &&
          ! -w $path){ #This chsould be writeable if we want to create it
        throw("Found path but is not writeable:\n\t$path");
      }
    }
    elsif(! -f $path){
      throw("Found path but is not a file:\n\t$path");
    }
    #could handle file create here?
  }
  elsif($is_dir && $create){
      #print STDOUT "Creating${label}:\n\t$path\n";

    if(! eval { make_path($path); 1 }){
      throw("Unable to create:\n\t$path\n$@");
    }
  }
  else{
    throw("Could not find:\n\t$path");
  }

  return $path;
}



# Move some of these to SetUtils?



#Is this really needed, could this not just be replaced with
#my @controls = grep { $_->is_control == 1 } @{$set->get_support};
#and test in caller, which is probably better?

#This also assumes InputSubset input if not ResultSet?!
#This arbitrarily returns the first, no all


#If this is just used for getting the control feature type, then we can probably to that all at once?
#and make this private?

#control experiment is now being handles by ResultSet itself
#should probably move feature_type method to ResultSet too

sub _get_a_control_InputSubset{
  my $set  = shift; 
  my @is_sets;
  
  if( $set->isa('Bio::EnsEMBL::Funcgen::ResultSet') ){
    #Enforce input_subset table_name here
    
    @is_sets = @{$set->get_support};
  
    if(! @is_sets){
      throw("Failed to identify control InputSubset support for ResultSet:\t".$set->name);  
    }
  }
  else{
    @is_sets = ($set);
  }
  
  my @ctrls = grep { $_->is_control == 1 } @is_sets;

  if(! @ctrls){
    throw('Could not identify a control InputSubset from '.ref($set).":\t".$set->name);  
  }

  return $ctrls[0];
}
  
#This may cause uncaught exception if there are no controls and 
#the caller does not test the return  
#Assumes is ResultSet with table_name input_subset
sub _get_control_InputSubsets{
  my $set = shift;
  return [ grep { $_->is_control == 1 } @{$set->get_support} ];  
}


#Should this conditionally append the analysis?
#This is only useful if your want to omit the _TRN 
#suffix.


sub get_set_prefix_from_Set{
  my $set        = shift;
  my $control    = shift; 
  my $study_name = get_study_name_from_Set($set, $control);
 
  if(! (ref($set) && 
        ($set->isa('Bio::EnsEMBL::Funcgen::InputSubset') ||
         $set->isa('Bio::EnsEMBL::Funcgen::ResultSet') ))){
    throw('Must pass a valid InputSet or ResultSet');         
  }
  
  my $ftype;
  
  if($control){
    $ftype = _get_a_control_InputSubset($set)->feature_type->name;
  }
  else{
     $ftype = $set->feature_type->name;
  }
 
  return $set->cell_type->name.'_'.$ftype.'_'.$study_name; 
}

#This currently only works for Experiments
#which have controls mixed in
    
#The pipeline does not currently expect 
#stand alone control experiments
#and will try and process these like a mixed signal/control set.
#This should be fine until after the alignment
#at which point we want to stop their processing
    
sub get_study_name_from_Set {
  my $set     = shift;
  my $control = shift;
  
  if(! (ref($set) && 
        ($set->isa('Bio::EnsEMBL::Funcgen::InputSubset') ||
         $set->isa('Bio::EnsEMBL::Funcgen::ResultSet') ))){
    throw('Must pass a valid InputSet or ResultSet');         
  }
  
  my ($exp_name, $ftype);
  my $ctype = $set->cell_type->name;
  
  if($control){
    $control = _get_a_control_InputSubset($set);
    #This is based on the assumption that all non-ctrl subsets are associated with a unique Experiment.
    my $exp   = $control->experiment;
    $exp_name = $exp->name;
 
    foreach my $isset(@{$exp->get_InputSubsets}){
      
      if(! $isset->is_control){
        $ftype = $isset->feature_type->name;  
        last;
      }  
    }
    
    if(! $ftype){   #We have a pure control experiment i.e. no signal InputSubsets
     $ftype = $exp->feature_type->name;
    }      
  }
  else{  
    $exp_name = $set->experiment->name ||
      throw('Cannot find unique experiment name for '.ref($set).":\t".$set->name);
    $ftype = $set->feature_type->name;
  }
  
  #\Q..\E escape mata characters in string variables 
  (my $study_name = $exp_name) =~ s/\Q${ctype}_${ftype}\E_(.*)/$1/i;
  
  if($study_name eq $exp_name){
    throw("Failed to create study name for Experiment $exp_name with cell type $ctype and feature type $ftype");  
  }
  
  return $study_name;
}

#This will also take a namespace

sub validate_package_path{
  my $pkg_path = shift ||  throw('Must defined a package path to validate');  
 
  #Quote so eval treats $aln_pkg as BAREWORD and converts :: to /
  if(! eval "{require $pkg_path; 1}" ){ 
    throw( "Failed to require:\t$pkg_path\n$@" );
  }
  
  #This might not always be correct if the file contains >1 package
  return path_to_namespace($pkg_path);
}

#This is to get around the problem that some scripts are not available in the top level of the $PATH dirs
#and creating a toplevel link is not appropriate as we need there full path.

#Shouldn't call this on a full path as will fail and take a long time about it
#Also,  we already know the location!

sub which_path{
  my $filename  = shift;
  my @env_paths = split(/:/, $ENV{PATH});
  my $path;
   
  find(sub{ if ($_ eq $filename){$path=$File::Find::name; return; }}, @env_paths);
  
  if(! defined $path){
    throw("Unable to find file in \$PATH:\t$filename\n\$PATH contains:\t".join("\t", @env_paths));
  }
  
  return $path;
}


1;

