#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor
#

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -user       => 'user',
   -pass       => 'pass',
   -dbname     => 'dbname',
   -host       => 'host',
   -species    => 'species',
   -DNADB_USER => 'dnadb_user',
   -DNADB_HOST => 'dnadb_host',
   -DNADB_NAME => 'dnadb_name',
   -DNADB_PORT => dnadb_port
);



my $rs_aadaptor   = $efgdba->get_ResultSetAdaptor;
my $dmf_adaptor   = $efgdba->get_DNAMethylationFeatureAdaptor;
my $slice_adaptor = $efgdba->dnadb->get_SliceAdaptor;

my ($rset) = @{ $rs_adaptor->fetch_all_by_name('ES_5mC_Stadler2011_PMID22170606') };
my $slice  = $slice_adaptor->fetch_by_region( 'chromosome', 1, 3010493, 3011550 );

my $dmf_ref = $dmf_adaptor->fetch_all_by_Slice_ResultSet($slice, $rset);

#print result_set cell type, feature_type and analysis here
#or leave to display label?

foreach my $dmf ( @{$dmf_ref} ) {
    print "Display label:\t"      . $dmf->display_label."\n";
    print "Location:\t"           . $dmf->feature_Slice->name."\n";
    print "Methylated reads:\t"   . $dmf->methylated_reads."\n";
    print "Total reads:\t"        . $dmf->total_reads."\n";
    print "Percent Methylation:\t". $dmf->percent_methylation."\n";
    print "Context:\t"            . $dmf->context."\n";
}


=head1 DESCRIPTION

The DNAMethylationFeatureAdaptor is a file based adaptor and inherits from a format
specific wrapper class, which provides generic methods to handle format specific parsers.
It provides an interface to retrieve DNAMethylationFeature objects.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DNAMethylationFeature
Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor
Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw );
#use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor; 
#Only need BaseFeatureAdaptor for non-DB based methods, e.g. projection, new_fast etc.
use Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor;

use base qw( Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor ); #Should also inherit from Funcgen BigBedAdaptor

#This 'mix-in' inheritance should be well defined, similar to a Moose 'has' (role) style relation ship.
#It should provide a well defined interface, which does not corrupt the interface of the inheriting class
#Any method naming clashes as these will never be executed in the mix-in, for example, should not contain 
#a new method.  But this behaviour could be employed to redefine a mix-in method in the inheritor

# To Do
#
# 1 Add file test and eval to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor::bigbed_open
#
# 2 Migrate fetch_features (with coderef support) back to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
#
# 3 Create funcgen BigBedAdaptor/BigBedCache|Handler which caches multiple Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptors
#
# 4 Add support for collection files above a certain slice length. Does this mean the web should just use out API from the start
#   Or do we need to make the method which decides what file to use available?
#   This depends on what we use for the methylation density features i.e. bigBed or .cols?
#
# Caveats
#
# This implementation adds extra dependancies on:
#    Bio::EnsEMBL::ExternalData
#     which in turn depends on
#    Bio::EnsEMBL::Web !!!
#    Bio::DB::BigFile


=head2 fetch_all_by_Slice_ResultSet

  Arg[1]     : Bio::EnsEMBL::Slice
  Arg[2]     : Bio::EnsEMBL::Funcgen::ResultSet
  Arg[3]     : HASHREF (optional) - valid contraint params e.g.
                 {
                  min_read_depth  => 5,
                  context         => 'CG',
                  min_methylation => 25,     # Float percentage
                  max_methylation => '75.5', # Float percentage
                 }
  Example    : my @dna_meth_feats = @{$dna_mf_adaptor->fetch_all_by_Slice_ResultSet($slice, $result_set)};
  Description: Fetches DNAMethylationFeatures for a given Slice and ResultSet
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DNAMethylationFeature objects
  Exceptions : Throws if args not valid and if ResultSet feature class is not dna_methylation
  Caller     : General
  Status     : Stable

=cut

#All calls in this method should be generic wrt file format

sub fetch_all_by_Slice_ResultSet{
  my ($self, $slice, $rset, $params) = @_;

  # VALIDATE ARGS
  if(! (defined $slice && $slice->isa('Bio::EnsEMBL::Slice') )){
    throw("Slice argument is not a valid Bio::EnsEMBL::Slice:\t$slice");
  }
  
  if(! (defined $rset && $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet') )){
    #Doesn't need to be stored
    throw("ResultSet argument is not a valid Bio::EnsEMBL::Funcgen::ResultSet:\t$rset");
  }

  if($rset->feature_class ne 'dna_methylation'){
    throw("ResultSet ".$rset->name." feature class is not dna_methylation:\t".$rset->feature_class);
  }

  #Here we need to decide which file/window_size to use based on the slice length
  #Need to lift code from ResultFeatureAdaptor::set_collection_config_by_Slice_ResultSets{
  #and available in generic CollectionAdaptor

  my $constructor_args = {
                          new_args => {
                                       slice => $slice,
                                       set   => $rset
                                      }
                         };
  

  #Validate constraints

  if(defined $params){

    if( ref($params) ne 'HASH'){
      throw("'params' argument must be a valid HASHREF of contraint key value pairs");
    }

    my $validate_method;

    foreach my $constraint(keys %{$params}){
      #No this all needs to be moved to the fetch method
      $validate_method  = '_validate_'.$constraint;

      if($self->can($validate_method)){
        $self->$validate_method($params->{$constraint})
      }
      else{
        #Could have these defined in a hash to print helpfully here
        throw("$constraint is not a valid DNAMethylationFeatureAdaptor constraint");
      }
    }
    
    #define contraints as separate key to avoid ne 'new_args' when constraining on each feature
    $constructor_args->{constraints} = $params;
  }
 

  my $constructor_wrapper = (defined $params) ? 
    '_create_DNAMethylationFeature_from_file_data_constraints' :
      '_create_DNAMethylationFeature_from_file_data';
  
  #calling fetch_features on self but should move back to
  #Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
  
  return $self->fetch_features
    (
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
     $rset->dbfile_data_dir,
     $self->can($constructor_wrapper),#return a code ref
     $constructor_args,
     $self
    );
}




#All the following methods should ultimately be moved elsewhere

#Move the following to funcgen BigBedAdaptor? (BigBedCache?)
#We shouldn't really have the same name, as it is not doing the same thing
#It's more like a cache
#What should the namespace of this be in the funcgen API?
#Bio::EnsEMBL::Funcgen::ExternalData::BigFile::BigBedAdaptor
#Bio::EnsEMBL::Funcgen::DBFile::BigFile::BigBedAdaptor
#which caches and delegates to separate Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptors



#Document bed format here

#Currently uses a 6 column format (conforming to bed standard) defined as follows:
#seq_id  start(0 based)  end    context/total_reads  score(1-1000)  strand(-/+/.)
#
#e.g.
#1       567890          567891 CG/120               950            +



my %strands = (
               '.' => 0,
               '+' => 1,
               '-' => -1,
              );

=head2 get_DNAMethylationFeature_params_from_file_data

  Arg[1]     : ARRAYREF - Split but unprocessed bed row (6 columns)
  Arg[2]     : Bio::EnsEMBL::Slice
  Example    : my ($chr, $start, $end, $strand, $methylated_reads,
                   $total_reads, $percent_methylation, $context) = 
               @{$dna_mf_adaptor->get_DNAMethylationFeature_params_from_file_data(\@bed_row)};
  Description: Processes a line of bed data to return valid parameters required
               for DNAMethylationFeature contructor.
  Returntype : ARRAYREF
  Exceptions : Throws if args not valid
  Caller     : General
  Status     : At risk

=cut

sub get_DNAMethylationFeature_params_from_file_data {
  my ($self, $bed_data, $slice) = @_;

  if(! (defined $bed_data && ref($bed_data) ne 'ARRAYREF') ){
    throw("bed_data argument in not an ARRAYREF:\t$bed_data");
  }
  
  #Most likely already validated in fetch method
  #Remove here for speed?
  if(! (defined $slice && $slice->isa('Bio::EnsEMBL::Slice') )){
    throw("slice argument is not a valid Bio::EnsEMBL::Slice:\t$slice");
  }
  
  my $slice_start = $slice->start;
  my $slice_end   = $slice->end;
  my $total_reads;
  my ($chr, $start, $end, $context, $percent_methylation, $strand) = @$bed_data;
  #$chr is actually never used in the Ensembl code as we already have the query slice

  $start       = $start - $slice_start;       #half open
  $end         = $end   -  $slice_start + 1 ;
  
  ( $context, $total_reads ) = split /\//, $context;
  $percent_methylation /= 10;  
  my $methylated_reads  = sprintf "%.0f", ( ( $percent_methylation / 100 ) * $total_reads );
  $strand = $strands{$strand} if defined $strand; #undef strand is not valid bed!

  return [ $chr, $start, $end, $strand, $methylated_reads,
           $total_reads, $percent_methylation, $context ];
}


=head2 get_file_adaptor

  Arg[1]     : String - Path to bigBed file
  Example    : my $bb_file_adaptor = $bb_cache_adaptor->get_file_adaptor($bb_path)
  Description: Generates and caches a BigBedAdaptor given a file path
  Returntype : Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
  Exceptions : Warns if cannot define BigBedAdaptor
  Caller     : Self
  Status     : At risk

=cut

#Make private?

# This will ultimately move to funcgen BigBedAdaptor
# This should use Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor


#This one should probably be in a base class

sub get_file_adaptor{
  my ($self, $path) = @_;

  #This doesn't need to be generically named as it will always be called from
  #funcgen BigBedAdaptor, not a feature adaptor


  #This is similar to some to core FileAdaptor::get_filehandle code

  #my $bb_path = $rset->dbfile_data_dir;
  #my $rset_id = $rset->dbID;


  if(! exists $self->{result_set_file_adaptors}{$path}){

    #Removed a lot of this as it should be in external BigBedAdaptor::bigbed_open
    #or whatever new IO core module is created

    #my $bb_path = $rset->dbfile_data_dir;

    #unless ( ( -e $bb_path ) && 
    #         ( -r $bb_path ) 
    #       ) {
    #  throw("The ResultSet dbfile ' $bb_path ' is not a valid path or the file is not readable");
    #}

    #eval to catch croak from BigFile
    #eval { $self->{result_set_file_adaptors}{$rset_id} 
    #         = Bio::DB::BigFile->bigBedFileOpen($bb_path) };

    #if($@){
    #  throw("Failed to open bigBedFile:\t$bb_path\n$@");
    #}

    $self->{result_set_file_adaptors}{$path} = 
      Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor->new($path);

  }
  
  if(! defined $self->{result_set_file_adaptors}{$path}){
    warn "Encountered error when opening bigBedFile\t".$path;
  }
  
  return $self->{result_set_file_adaptors}{$path};
}

=head2 _create_DNAMethylationFeature_from_file_data

  Arg[1]     : ARRAYREF - Split but unprocessed bed row (6 columns)
  Arg[2]     : HASHREF  - Containing additional required 'new_args' parameters for 
               DNAMethylationFeature::new_fast constructor:
                 {
                  new_args => {
                               slice => $slice,
                               set   => $result_set,
                              }
                 }
  Example    : my $constructor_wrapper_ref = $self->can('_create_DNAMethylationFeature_from_file_data');
               #Then passed as arg to BigBedHandler/Adaptor::fetch_features and called as follows
               my $feature = $constructor_wrapper_ref->
                               (
                                $constructor_wrapper_obj, #explicitly pass containing object when using coderef
                                $bed_data_ref,
                                $constructor_wrapper_args
                               );
  Description: Internal method which is passed as a code reference to the 
               fetch_features method of the specific file format adaptor
  Returntype : Bio::EnsEMBL::Funcgen::DNAMethylationFeature
  Exceptions : None
  Caller     : fetch_all_by_Slice_ResultSet and BigBedHandler/Adaptor::fetch_features
  Status     : At risk

=cut


sub _create_DNAMethylationFeature_from_file_data{
  my ($self, $bed_data, $args) = @_;

  #Could validate $args is HAHREF here
  #But omit for speed as is private method

  $bed_data = $self->get_DNAMethylationFeature_params_from_file_data
    (
     $bed_data,
     $args->{new_args}{slice}
    );

  return Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new_fast
    ({
      start               => $bed_data->[1],
      end                 => $bed_data->[2],
      strand              => $bed_data->[3],
      methylated_reads    => $bed_data->[4],
      percent_methylation => $bed_data->[5],
      total_reads         => $bed_data->[6],
      context             => $bed_data->[7],
      adaptor             => $self,
      %{$args->{new_args}}                     # args are slice and result_set
     });
}

=head2 _create_DNAMethylationFeature_from_file_data

  Arg[1]     : ARRAYREF - Split but unprocessed bed row (6 columns)
  Arg[2]     : HASHREF  - Containing additional required 'new_args' parameters for 
               DNAMethylationFeature::new_fast constructor and contraint definitions :
                 {
                  new_args    => {
                                  slice => $slice,
                                  set   => $result_set,
                                 }
                  constraints => {
                                  min_read_depth      => 10,
                                  context             => 'CG', #
                                  percent_methylation => 80,   #Minimum
                                 }
                 }
  Example    : my $constructor_wrapper_ref = $self->can('_create_DNAMethylationFeature_from_file_data');
               #Then passed as arg to BigBedHandler/Adaptor::fetch_features and called as follows
               my $feature = $constructor_wrapper_ref->
                               (
                                $constructor_wrapper_obj, #explicitly pass containing object when using coderef
                                $bed_data_ref,
                                $constructor_wrapper_args
                               );
  Description: Internal method which is passed as a code reference to the 
               fetch_features method of the specific file format adaptor
  Returntype : Bio::EnsEMBL::Funcgen::DNAMethylationFeature
  Exceptions : None
  Caller     : fetch_all_by_Slice_ResultSet and BigBedHandler/Adaptor::fetch_features
  Status     : At risk

=cut

sub _create_DNAMethylationFeature_from_file_data_constraints{
  my ($self, $bed_data, $args) = @_;
  
  #loop through the contraints
  #calling methods based on the contraint key, passing the constraint value.
  #This will mean we won't have to do this check here for every single feature if we don't want to filter
  #just once in the caller!

  #$bed_data is not yet processed, so get_DNAMethylationFeature_params_from_file_data
  #and pass attrs directly here

  $bed_data = $self->get_DNAMethylationFeature_params_from_file_data
    (
     $bed_data,
     $args->{new_args}{slice}
    );
  

  my $constrain_method;
  
  foreach my $constraint(keys %{$args->{constraints}} ){
    $constrain_method = '_constrain_'.$constraint;
       
    if( $self->$constrain_method($args->{constraints}{$constraint}, $bed_data) ){
      #warn "constraining $constraint";
      return;
    }
  }
  
  return Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new_fast
    ({
      start               => $bed_data->[1],
      end                 => $bed_data->[2],
      strand              => $bed_data->[3],
      methylated_reads    => $bed_data->[4],
      percent_methylation => $bed_data->[5],
      total_reads         => $bed_data->[6],
      context             => $bed_data->[7],
      adaptor             => $self,
      %{$args->{new_args}}                     # args are slice and result_set
     });
}

#Have to be class methods rather than subs
#when using strict refs. subs would probably be a little faster

### Private _contraint and _validate methods
### _validate methods simply validate specified contraint param
### _constraint methods return true if contraint is to be applied

# Some of these are generic and would sit better in the BigBedAdaptor itself
# or integrated into the fetch_features query directly
# This may break the generic implementation?

sub _constrain_min_read_depth{
  my ($self, $min_read_depth, $bed_data) = @_;
  return ($bed_data->[6] < $min_read_depth) ? 1 : 0;
}

#Validate methods could return defaults if param not defined?

sub _validate_min_read_depth{
  my ($self, $min_read_depth) = @_;

  if(! defined $min_read_depth ||
     $min_read_depth !~ /^\d+$/){
    throw($min_read_depth.' is not a valid integer to contrain using min_read_depth')
  }
  return;
}


my %valid_contexts = (
                      CG  => undef,
                      CGG => undef,
                      CAG => undef,
                      CTG => undef,
                     );

##ambiguity codes
#CNG => {N => [A, C, G]},
#CHG => {H => [A, C, T, G]},
#CHH => {H => [A, C, T, G]},
#CHN? - would need separate %ambiguity_codes
#CH
#CN
#C[^GATC]
#Non-cytosine methylation? i.e. EG support for plants!
#Will the broad range of possible contexts
#make this almost redundant?


sub _constrain_context{
  my ($self, $context, $bed_data) = @_;
  return ($bed_data->[7] ne $context) ? 1 : 0;
}

sub _validate_context{
  my ($self, $context) = @_;

  if(! defined $context ||
     ! exists $valid_contexts{$context} ){
    throw($context." is not a valid context to contrain, please use one of the following:\n\t"
          .join("\t", keys %valid_contexts));
  }

  return;
}

#probably want min and max operations here

sub _constrain_min_methylation{
  my ($self, $perc, $bed_data) = @_;
  return ($bed_data->[5] < $perc) ? 1 : 0;
}

sub _constrain_max_methylation{
  my ($self, $perc, $bed_data) = @_;
  return ($bed_data->[5] > $perc) ? 1 : 0;
}


#Could warn if erroneously set to 0 or 100 for max and min respectively

sub _validate_min_methylation{
  my ($self, $perc_float) = @_;

  if( (! defined $perc_float)              ||
      ($perc_float !~ /^[0-9]*\.?[0-9]+$/) ||
      ($perc_float > 100) ){
    throw($perc_float.' is not a valid positive percentage to contrain using min_methylation');
  }
  return;
}

sub _validate_max_methylation{
  my ($self, $perc_float) = @_;

  if( (! defined $perc_float)              ||
      ($perc_float !~ /^[0-9]*\.?[0-9]+$/) ||
      ($perc_float > 100) ){
    throw($perc_float.' is not a valid positive percentage to contrain using max_methylation');
  }
  return;
}

#other contraints
#strand?


=head2 fetch_features

  Arg[1]     : String   - sequence id or name
  Arg[2]     : Int      - Genomic start (1 based)
  Arg[3]     : Int      - Genomic end
  Arg[4]     : String   - Path to bigBed file
  Arg[5]     : CODEREF  - Reference to Feature constructor wrapper method/sub
  Arg[6]     : HASHREF|ARRAYREF|SCALAR (optional) - constructor wrapper args
  Arg[7]     : Object (optional)  - Object containing constructor

  Example    : my @feats = @{$bb_adaptor->fetch_features
                             (
                              $chr_name,
                              $start,
                              $end,
                              $bigBed_path,
                              $contructor_wrapper_coderef,
                              $constructor_wrapper_args_ref,
                              $constructor_wrapper_obj,
                             )};
  Description: 
  Returntype : ARRAYREF likely, but depends on return type of constructor wrapper reference
  Exceptions : Warns if cannot open bigBeg file
               Throws is constructor wrapper argument is not a CODEREF
  Caller     : General
  Status     : At risk

=cut


#This should be completely agnostic to calling adaptor
#And should really be migrated back to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
#as it is generic due to the passing of a coderef

#Will still need a wrapper method in the Funcgen BigBed adaptor
#which should probably be renamed and standardised across all
#file adaptor wrapper

sub fetch_features {
  my ($self, $chr_id, $start, $end, $bb_path, $constructor_ref, $constructor_args, $constructor_obj) = @_;
  #Have to pass CODEREF rather than method name string
  #to have standard valid way of calling both subs and class methods
    
  my $bb_adaptor = $self->get_file_adaptor($bb_path);
  my $bb         = $bb_adaptor->bigbed_open;

  warn "Failed to open BigBed file" . $bb_path unless $bb;
  return [] unless $bb;
  
  #  Maybe need to add 'chr' 
  my $seq_id = $bb_adaptor->munge_chr_id($chr_id);
  
  return [] if ! defined $seq_id;

  if(ref($constructor_ref) ne 'CODE'){
    throw($constructor_ref.' is not a valid CODEREF required from the constructor method');
  }

  my @constructor_args = ($constructor_obj, undef, $constructor_args);
  my $bed_data_index   = 1;
  
  if(! defined $constructor_obj){
    #validate ref and is Object?
    $bed_data_index = 0;
    @constructor_args = (undef, $constructor_args);
  }


  # Remember this method takes half-open coords (subtract 1 from start)
  my $list_head = $bb->bigBedIntervalQuery($seq_id, $start-1, $end);
  my @features;
 

  for ( my $i=$list_head->head; $i; $i=$i->next ) {

    #my @bedline = ($chr_id,$i->start,$i->end,split(/\t/,$i->rest));
    #my $bed = EnsEMBL::Web::Text::Feature::BED->new(\@bedline);
    #$bed->coords([$chr_id,$i->start,$i->end]);
    ### Set score to undef if missing to distinguish it from a genuine present but zero score
    #$bed->score(undef) if @bedline < 5;
    
    #Passing $chr_id is redundant as we already know this from the query slice

    $constructor_args[$bed_data_index] = [ $chr_id, $i->start, $i->end, split(/\t/,$i->rest) ];
    my $feature = $constructor_ref->(@constructor_args);

    #do we need to handle undef score here?

    #This defined check is to allow the constructor to return undef
    #if there is any sort of filtering going on
    push @features, $feature if defined $feature;
  }
  
  return \@features;
}

#Add _true_tables required by BaseAdaptor
#This can go once the file and DB BaseAdaptors are split

sub _true_tables {
  return ([ ]);
}


1;
