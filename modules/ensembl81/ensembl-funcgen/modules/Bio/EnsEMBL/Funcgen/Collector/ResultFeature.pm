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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::Funcgen::Collector::ResultFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( open_file );
use Bio::EnsEMBL::Funcgen::Collection::ResultFeature;
use Bio::EnsEMBL::Funcgen::ProbeFeature; #Only used for _pre_storing slice/seq_region details
#use POSIX;#ceil

use base qw( Bio::EnsEMBL::Utils::Collector Bio::EnsEMBL::DBFile::CollectionAdaptor );


### Global config variables
# See associated Collector methods for more info
$Bio::EnsEMBL::Utils::Collector::bin_model  = 'SIMPLE';


#WARNING Need to change this package var as it is persisting and growing between instances
#within the same process! i.e. workers in the hive!
$Bio::EnsEMBL::Utils::Collector::window_sizes = [30, 65, 130, 260, 450, 648, 950, 1296];


# Tuning 900(Most used display size) - 172(non-drawable area) = 772 pixels/bins
# Optimal window sizes over zoom levels
# 1kb 5kb 10kb 50kb 100kb 200kb 500kb 1mb
# 2   7   13   65   130   260   648   1296

#Theoretical max for default 16MB max_allowed_packet_size is 30(~23kb) - sufficient min for read coverage.
#0 will be used for array intensies until 65 is more appropriate(see adaptor)
#450 and 950 added to handle middle ground, but are these really needed?
#New windows will use more memory due to smaller window sizes
#But we have dropped a window size from the mid-upper range
#Largest size is essential so we aren't retrieving/drawing more features than pixels
#Smaller sizes are desirable as more resolution required when zoomed in.

#count - v small(4) integers(little endian)
#RPKM is now a f(loat).  Which is native, not option to force endianess with packed floats
#So may cause problems if using a host machine with a different endian order
$Bio::EnsEMBL::Utils::Collector::bin_method   = 'RPKM'; #only used by collector
$Bio::EnsEMBL::Utils::Collector::packed_size   = 4;     #per score
$Bio::EnsEMBL::Utils::Collector::pack_template = 'f';   #per score



### Mandatory methods required by the base Collector
# You must have a methods like this to perform the store

sub store_window_bins_by_Slice_ResultSet {
  my ($self, $slice, $rset, %config) = @_;

  $self->source_set_type('result');#required by get_Feature_by_Slice   
  $self->set_collection_defs_by_ResultSet($rset);

  

  $self->set_config(%config);


  #Need to test for existing collection

  $self->store_window_bins_by_Slice($slice);
  
  return;
}



sub store_window_bins_by_Slice_Parser{
  my ($self, $slice, $imp, %config) = @_;

  
  #warn " store_window_bins_by_Slice_Parser ".$slice->name;
  #warn "wsizes before set ".join(', ', @{$self->window_sizes});
  #warn "wsizes before set ".join(', ', @{$Bio::EnsEMBL::Utils::Collector::window_sizes});
  
  
  #We need to test for new_assm if we have a reads result set and fail
  #currently don't have access to new assm for validation?
  #Can project to new_assm as this would require two passes storing initially on the 0 window level
  #which we don't want to do
  #throw if new_assm as we need to remap before running this
  #This is done in Collector so why are we parsing it here?
  #skip_zero window will always be set, so this would fail!!!

  #test parse config here, will this strip it out of the hash before
  #we pass it to the super method?
  #Also need to test skip_zero_window or window_sizes dependant on result set type(sequencing, array)
  #This needs to be defined when craeting the ResultSet in the Importer
  #NEED TO CHANGE ALL ResultSet generation to add type!!!!!
  
  my ($skip_zero_window, $force) = rearrange( [ 'SKIP_ZERO_WINDOW',  'FORCE'], %config );



  #There is no way of setting IMPORTED status for slice based jobs
  #We need the accumulator to set the IMPORTED/RESULT_FEATURE_SET status
  #Would need to turn parse_and_import into RunnableDB?
  #To run these as slice jobs we would need to code to configure some slice based input ids


  $self->source_set_type('input');#required by get_Feature_by_Slice
  $self->set_collection_defs_by_ResultSet($imp->result_set);

  $self->parser($imp);
 
  #For safety, set skip_zero window if we are using SEQUENCING data
    
  if(! $skip_zero_window){
	  #Assume we only have one set here(enforced in define_and_validate_sets)
	  my ($iset) = @{$imp->result_set->get_support};

    if(! $iset->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip')){
	    $config{'-skip_zero_window'} = 1;
	  }
  }
 
  warn "wsizes ".join(', ', @{$self->window_sizes});

  $self->set_config(%config, (-method_config => 
                              {
                               #RPKM method config
                               -dnadb          => $imp->db->dnadb,
                               -total_features => $imp->total_features,
                               -gender         => $imp->result_set->cell_type->gender,
                               -window_sizes   => $self->window_sizes,
                               #pass here as Collector is currently 'readding' 0 wsize
                              }
                             ));

 
  ### Check for existing dbfile dir/files
  my $rset_dir        = $self->result_set->dbfile_data_dir;
  my $dbfile_data_dir = $imp->get_dir('output');

  #This is causing problems with mismatch between
  #output dir and staging dir generated from meta dbfile_data_dir
  #Need to use force!

  if((! $force ) &&
	 $rset_dir   &&
	 ($rset_dir ne $dbfile_data_dir)){
	throw("ResultSet dbfile_data_dir($rset_dir) and -output_dir($dbfile_data_dir) do not match. Please rectify or specify -force to update");
	
  }
  
  #Update/set ResultSet dbfile_data_dir
  if( (! $rset_dir) ||
	  ( $rset_dir                       &&
		($rset_dir ne $dbfile_data_dir) &&
		$force ) ){
	$self->result_set->dbfile_data_dir($dbfile_data_dir);
	$self->result_set->adaptor->store_dbfile_data_dir($self->result_set);
  }

  if((! defined $dbfile_data_dir) ||
	 (! -d $dbfile_data_dir)){
	throw('Your -dbfile_data_dir is either not set or not a valid directory');
  }


  #Check and open each window_size file
  for my $wsize(@{$self->window_sizes}){

	if($wsize == 0){
	  #Would not normally happen
	  warn ("Need to fix 0bp window_size re-adding in Collector?");
	  next;
	}
	else{
	  my $col_file_name = $self->result_set->get_dbfile_path_by_window_size($wsize, $slice);

	  if(-e $col_file_name && $force){
		unlink($col_file_name) or warn("Failed to remove exisiting col file:\t$col_file_name\n$!");
		#no throw as file will be over-written anyway
	   #This was removing open file an existing fh when there are duplicate wsizes!
		#Do we need to some sort of cache check first?
	  }

      warn "Getting filehandle for:\t".$col_file_name;

	  #Generate and cache filhandle here to avoid failing after data has been generated
	  my $fh = $self->get_filehandle($col_file_name, {-file_operator => '>'});
	  
	  #if(! defined $fh){
	  #throw("Could not get_filehandle for >${col_file_name}");
	  #}
	  
	  #Set AUTOFLUSH to enable validate_file_length in store_collection
	  $fh->autoflush;
	}
  }

  $self->store_window_bins_by_Slice($slice);

  #Index creation and merging done in post-processing script
  return;
}

# add rollback method here

=head2 get_Features_by_Slice

  Args[0]    : Bio::EnsEMBL::Slice
  Example    : $collector->get_Feature_by_Slice($slice);
  Description: Wrapper method to fetch input features for building the Collections
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Collection::ResultFeatures
  Exceptions : None
  Caller     : Bio::EnsEMBL::Utils::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

#Can remove this method if we pass a method ref
#Default for normal FeatureAdaptors would fetch_all_by_Slice
#However cannot do this as we don't get self passed via a code ref
#we could however manually pass this when we call the method?
#self will either be a collector or a parser
#Just maintain this wrapper for now for simplicity

#CR Put stub in Colelctor.pm

sub get_Features_by_Slice{
  my ($self, $slice) = @_;

  #Can we pass this as method ref to Collector, this would prvent the need to write
  #this wrapper for standard fetch_all_by_Slice based access
  my $features;
  my $source_set_type = $self->source_set_type;

  if ($source_set_type eq 'result'){
	#Add more args here? status?
	#Add default window size 0, if ResultSet is alread a RESULT_FEATURE_SET
	#This may occur if you are projecting features from an old assembly to 
	#0 window size before generating the other windows in a second pass
	$features = $self->result_set->get_ResultFeatures_by_Slice($slice, undef, undef, undef, 0);
  }
  elsif($source_set_type eq 'input'){
    $features = $self->parser->parse_Features_by_Slice($slice);

	#This method assumes a sorted file handle
	#Can we set markers for disk seeking on a sorted handle?
	#Either we make it handle an slice passed
	#Or we assume the next query slice will be after the last
	#Just restrict to one slice at a time for now
	
	#Also needs to hold cache of long features
	#Kind of reinventing the BaseFeatureAdaptor wheel here?


  }
  else{
	#We should have already validated this by now
	throw('get_Features_by_Slice only support input and result set_types');
  }

  return $features;

}



=head2 write_collection

  Args[0]    : int : window_size
  Args[1]    : Bio::EnsEMBL::Slice
  Args[2]    : optional int : Feature end
  Args[3]    : optional int : Feature strand
  Args[4]    : optional ARRAYREF of scores
  Example    : $self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  Description: Writes Collections by caching scores, defining collection end and strand and 
               storing when end of collection seen or max_packed_size is exceeded (currently does 
               not support fragmented collections). Stores last Collection if only window_size
               and Slice are passed.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Utils::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

sub write_collection{
  my ($self, $wsize, $slice, $slice_end, $strand, $scores) = @_;

  if(defined $slice_end){
	#collection_start is defined in Collector for first bin in collection
	$self->collection_end($wsize, $slice_end);
	$self->collection_strand($wsize, $strand);
  }
  

  if(defined $scores){
	my $new_cps = $self->current_packed_size($wsize) + 
	                                         (scalar(@$scores) * $self->packed_size);	
	#This loop is very similar to the _obj_from_sth loop
	#We need to set the slice and modify the start end appropriately
	
	
	#This is no longer relevant as we are storing in files, not in a DB table
	#if($new_cps >= $self->max_data_type_size){
	#  warn("Have found $wsize collection larger that max_data_type_size(16MB). Need to implement cross collection querying");
	#  #Via a union of two substr queries
	#  #This is no set to 64MB, but this does not directly translate to the maximum size allowed in a single insert
	#}
	  


	if($wsize == 0){
	  #We need to 0 wsize collections to be projected as single features
	  #Other wise the assembly projection will not work
	  #Then use the 0 wsize ResultFeatures to generate

	  #if($new_cps >= $self->max_data_type_size){
	  #	warn('Have found collection larger that max_data_type_size(16MB). Need to implement cross collection querying');
	  #	#Via a union of two substr queries
	  #  }
	  

	  #Need to cache current score for 0 window
	  $self->score_cache($wsize, $scores) if $wsize == 0;
	  $self->store_collection($wsize, $slice, $self->collection_start($wsize), $slice_end, $strand);

	}
	
	$self->score_cache($wsize, $scores) if $wsize != 0; 
  }
  else{#No score info, so we store the remaining last collections for each wsize
	$self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  }

  return;

}


#To get this to print to file, we need to set up a slice/window specific file handle
#Doing this on a slice basis would also reduce the initial seek time
#at the expense of the number of open file handles (similar to partitions)
#Hence need to cat together slices for same window in defined order in a post processing step


=head2 get_score_by_Feature


=cut 

#Not required for count/density based collections

sub get_score_by_Feature{
  my ($self, $feature) = @_;

  #For speed, assume we have a valid Bio::EnsEMBL::Funcgen::Collection::ResultFeature
  #This will always be a 0 wsize collection, hence we always want the first score

  return $feature->scores->[0];

}



=head2 reinitialise_input 


=cut 

#optional method - resets file handle for iterative slice parsing for different chunk length sets

sub reinitialise_input{
  my $self = shift;

  #No need to do this for DB queries
  return if $self->source_set_type eq 'result';

  #This will make all the parser counts go screwy
  #as we are potentially counting the same features again


  #Move this to parser method
  #change to seek fh, 0, 0
  $self->parser->file_handle->close;
  $self->parser->file_handle( open_file($self->parser->input_file, $self->parser->input_file_operator) );

  #Reset caches 
  $self->parser->{'last_slice'} = undef;
  $self->parser->{'overhang_features'} = [];
  $self->parser->last_line('');
  return;
}


##############################################################
### Following methods are not mandatory for a Collector
### but support above methods for this ResultFeature Collector
##############################################################




=head2 store_collection

  Args[0]    : int : window_size
  Args[1]    : Bio::EnsEMBL::Slice
  Args[2]    : int : Feature start
  Args[3]    : int : Feature end
  Args[4]    : int : Feature strand
  Example    : $self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  Description: Collection storage method. Resets seq_regios_start/end and scores appropriately
               if collection exceeds storage slice. Resets score cache and next collection_start.
  Returntype : None
  Exceptions : None
  Caller     : write_collection and Bio::EnsEMBL::Utils::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

sub store_collection{
  my ($self, $wsize, $full_slice, $slice_start, $slice_end, $strand) =  @_;
	  
  #warn "Storing collection $wsize, $full_slice, $slice_start, $slice_end, $strand";

  my $sr_start = $slice_start;
  my $sr_end   = $slice_end;
  $strand = 0 if ! defined $strand;
  #Overwriting strand value with 0 
  #As we have collected features from both strands?
  $strand = 0 if $wsize != 0;


  #This happens if the last collection was already written in the loop
  #Handle this here so we don't have to set/test collection_start for every record
  return if($sr_start == ($sr_end + 1));


  #Set store slice to start at 1 (PARs, test slices)
  #Can we remove this store_slice now?
  #How are we going to handle storing/fetching on PARs?
  my $store_slice = $full_slice;

  
  if($store_slice->start != 1){
	#Alter sr_start/end
	$sr_start  = $slice_start + $full_slice->start - 1;
	$sr_end = $slice_start + $full_slice->start - 1;
	$store_slice = $store_slice->adaptor->fetch_by_region(undef, $store_slice->seq_region_name);
  }

  
  ### Splice scores if collection overhang store slice
  #We reassign this below to avoid any weird ref updating behaviour
  my $scores_ref = $self->score_cache($wsize);

  if($sr_end > $full_slice->end){
	#Can happen for 0 wsize if we are running with small test slice
	#But we never want to splice scores for 0 wsize

	if($wsize == 0){  #Throw away any which are not at least partially on this slice
	  return if $sr_start > $full_slice->end;
	}
	else{

	  #Trim the scores and sr_end to the nearest bin
	  my $tmp_end = int($full_slice->end/$wsize) * $wsize;
	  $tmp_end += $wsize if $tmp_end < $full_slice->end;
	  my $overhang_length = ($sr_end - $tmp_end)/$wsize;
	  #This should always be an int as $slice_end is always a valid bin end
	  $sr_end = $tmp_end;
	  
	  #Now remove the surplus scores from the cache
	  splice(@$scores_ref, ($#{$scores_ref} - $overhang_length), $overhang_length);
	}

  }
 
  print 'Storing '.scalar(@$scores_ref)." bins(window_size=$wsize) for:\t".$store_slice->name."\n";


  if($wsize == 0){

	$self->store([Bio::EnsEMBL::Funcgen::Collection::ResultFeature->new_fast
				  ({
					start  => $sr_start,
					end    => $sr_end,
					strand => $strand,
					scores => $scores_ref,
					#probe  => undef,
					result_set_id => $self->result_set->dbID,
					window_size   => $wsize,
					slice         => $store_slice,#Full slice
				   }
				  )], $self->result_set, $self->new_assembly);
  }
  else{  #Write to col file!
   	#This needs moving to ResultFeatureAdaptor::write_to_file?

	#store slice should always be the full length slice here
	#and sr_start should always ==1
	#and sr_end should always >= slice->end/length

	my $col_file_name = $self->result_set->get_dbfile_path_by_window_size($wsize, $store_slice);
	my $fh = $self->get_filehandle($col_file_name);

	
	#, {-file_operator => '>'});#not strictly needed here as it should be open for writing?
	#Is this true? We are getting some file not exist errors on validate_file_length
	#this should fail before then anyway
	#get_filehandle just warn if it fails to acquire $fh
	

	# STORE SEQ_REGION
	if(! $self->get_seq_region_id_by_Slice($store_slice)){
	  $self->_pre_store(Bio::EnsEMBL::Funcgen::ProbeFeature->new
						(
						 -slice => $store_slice,
						 -start => 1,
						 -end   => 1,
						 -strand => 0,
						)
					   );	
	}
	
	# PACK AND PRINT
	my $pack_template = '('.$self->pack_template.')'.scalar(@{$scores_ref});
	#warn "Pack template is $pack_template";
	
	print "Writing to file:\t".$col_file_name."\n";
	print $fh pack($pack_template, @{$scores_ref});
	
	# VALIDATE
   	my $total_packed_size = $self->get_total_packed_size($sr_start, $sr_end, $wsize);
	#This will only work due to autoflush set when opening in store_window_bins_by_Slice_Parser{
	#This will never be true is the collection end does not match the slice end?
	#This effectively validates the pack template
	$self->validate_file_length($col_file_name, $total_packed_size, 1);#1 is binmode;
  }
  
  #Reassign/clean score_cache reference to avoid any reference updating problems
 
  $self->{'score_cache'}{$wsize} = [];

  #Set collection start, this is reset to the actual feature start for wsize == 0
  $self->collection_start($wsize, ($slice_end+1));
 
  return;
}




sub get_total_packed_size{
  my ($self, $sr_start, $sr_end, $wsize) = @_;

  #Can afford to validate here as this is not for fetch purposes
  my $start_bin = ($sr_start + $wsize - 1) / $wsize;
  my $end_bin   = $sr_end / $wsize;

  if(($start_bin != int($start_bin)) ||
	 ($end_bin != int($end_bin))){
	throw("The seq_region_start/end($sr_start/$sr_end) are not valid bin bounds for window_size $wsize");
  }

  my $bin_count = ($sr_end - $sr_start +1) / $wsize;
  #Could use POSIX::ceil here
  my $tmp_bc      = int($bin_count);
  $bin_count    = $tmp_bc + 1  if $tmp_bc != $bin_count;
  
  return ($bin_count * $self->packed_size);
}




=head2 set_collection_defs_by_ResultSet

  Args[0]    : Bio::EnsEMBL::Funcgen::ResultSet
  Example    : $self->set_collection_defs_by_ResultSet($rset);
  Description: Sets the packed_size and pack_template and bin_method dependant on the 
               ResultSet type(reads/array)
  Returntype : None
  Exceptions : throws if supporting InputSet is not of type result (i.e. short reads import)
               throws if ResultSet is not and input_set or experimental_chip based ResultSet (i.e. channel etc)
  Caller     : store_window_bins_by_Slice_ResultSet
               Bio::EnsEMBL::Funcgen:DBSQL::ResultFeatureAdaptor::fetch_all_by_Slice_ResultSet
  Status     : At Risk

=cut

#This needs merging with set_collection_config_by_ResultSets

sub set_collection_defs_by_ResultSet{
  my ($self, $rset) = @_;
  $self->result_set($rset);
  warn "set_collection_defs_by_ResultSet ".join(', ', @{$Bio::EnsEMBL::Utils::Collector::window_sizes});

  if(defined $rset){
    #This is most likely already done in the caller
    #$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
	
    if($rset->table_name eq 'experimental_chip'){ #Array Intensities i.e. single float
  	  #perl only offers native endian order for floats (http://www.perlmonks.org/?node_id=629530)
  	  #$Bio::EnsEMBL::Utils::Collector::packed_size       = 4;#per score
  	  #$Bio::EnsEMBL::Utils::Collector::pack_template     = 'f';#per score
  	  #Use defaults for these now
  
      $Bio::EnsEMBL::Utils::Collector::bin_method        = 'max_magnitude';#only used by collector
      $Bio::EnsEMBL::Utils::Collector::window_sizes->[0] = 0;#Can have natural resolution for low density array data
    }
    elsif($rset->table_name eq 'input_subset'){
  	  #Need to reset this as we may be doing serial queries.
      $Bio::EnsEMBL::Utils::Collector::window_sizes->[0] = 30;
  
  	  #Currently only expecting int from InputSet
  	  #my @isets = @{$rset->get_support};
  	  #my @tmp_isets = grep(!/result/, (map $_->feature_class, @isets));
  	  
  	  #if(@tmp_isets){
  		#throw("Bio::EnsEMBL::Funcgen::Collector::ResultFeature only supports result type InputSets, not @tmp_isets types");
  	  #}
  
  	  #We still have no way of encoding pack_type for result_feature InputSets
  	}
    else{
      throw('Bio::EnsEMBL::Funcgen::Collector:ResultFeature does not support ResultSets of type '.
        $rset->table_name);
    }
  }

  #Do we need to validate the smallest non-0 window size
  #against the max pack size?
  #This should be done in the Collector
  warn "set_collection_defs_by_ResultSet ".join(', ', @{$Bio::EnsEMBL::Utils::Collector::window_sizes});
  return;
}



=head2 result_set

  Args[0]    : optional Bio::EnsEMBL::Funcgen::ResultSet
  Example    : $self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  Description: result_set attribute method required for fetch wrapper method
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : throws if arg is not valid
  Caller     : get_Features_by_Slice
  Status     : At Risk

=cut

sub result_set{
  my ($self, $rset) = @_;

 
  #Can't use is_stored_and_valid here

  if($rset && ! (ref($rset) && $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet'))){
	throw('You must pass a valid Bio::EnsEMBL::Funcgen::ResultSet');
  }
  elsif($rset){

	#Now done in the caller
	#if($rset->has_status('RESULT_FEATURE_SET')){
	#throw('ResultSet('.$rset->name.') already has precomputed ResultFeatures stored, please rollback ResultFeature first');
	  #Retrieving ResultFeatures here would end up retrieving them from the result_feature table
	  #which is where we want to store them
	  #They are only retrived from the probe/result/probe_feature table if they do not have this status.
	  #REMEMBER TO SET THIS IN THE CALLING SCRIPT!?
	#}
  
	$self->{'result_set'} = $rset;
  }
  
  return $self->{'result_set'};
}


=head2 parser

  Args[0]    : optional Bio::EnsEMBL::Funcgen::Parsers::InputSet
  Example    : $self->parser($parser);
  Description: Getter/Setter for parser attribute if this ResultFeature Collector
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::InputSet
  Exceptions : throws if arg is not valid
  Caller     : general
  Status     : At Risk

=cut

sub parser{
  my ($self, $parser) = @_;
 
  #Can't use is_stored_and_valid here

  if($parser && ! (ref($parser) && $parser->isa('Bio::EnsEMBL::Funcgen::Parsers::InputSet'))){
	throw('You must pass a valid Bio::EnsEMBL::Funcgen::Parsers::InputSet');
  }
  elsif($parser){
  
	$self->{'parser'} = $parser;
  }
  
  return $self->{'parser'};
}




sub source_set_type{
  my ($self, $type) = @_;

  $self->{source_set_type} = $type if $type;

  return $self->{source_set_type};
}


sub merge_and_index_slice_collections {
  my ($self, $rset, $slices, $data_dir, $force_overwrite) = @_;  
  
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
  

  if(! -d $data_dir){
    die("data_dir argument is not valid:'t$data_dir\n".
      "This should contains the slice/seq_region .col files");
  }
  

  
  if(! (defined $slices &&
        ref($slices)    &&
        (ref($slices) eq 'ARRAY') &&
        (scalar(@$slices) > 0)  ) ){
    throw("Slices argument must be a valid Arrayref of Bio::EnsEMBL::Slice objects:\t$slices"); 
  }
  
  
  
  $self->set_collection_defs_by_ResultSet($rset);
  my $rset_name = $rset->name;
  my $packed_size = $self->packed_size;
  my @wsizes      = @{$self->window_sizes};
  my $first_wsize = 1;
  my (%index, $fh, @rm_files);
  
  
  foreach my $wsize(@wsizes){
  
    #Check for target file first
    my $merged_fname = join('.', ('result_features',  $rset_name, $wsize, 'col'));
    my $merged_file =  "${data_dir}/${merged_fname}";
  
    if(-f $merged_file){
  	  
  	  if(! $force_overwrite){
    	throw("Failed to over-write existing output file:\t$merged_file\n$!");
      }else{
    	  #print STDOUT "Over-writing existing output file:\t$merged_file\n";
  	  }
    }
  
    my $total_length = 0;
    my $off_set      = 0;
    my $last_length  = 0;
    my $seen_y       = 0;
    my (@efg_sr_ids, @file_list);
     
    for my $slice(@$slices){

      my $sr_name    = $slice->seq_region_name;
      my $slice_file = $data_dir.'/'.join('.', ('result_features',  $rset_name, $wsize, $sr_name, 'col'));
  
      if(! -f $slice_file){

        if(($sr_name =~ /^[0-9]+$/) ||
           ($sr_name eq 'X') ||
           ($sr_name eq 'Y')) {
          throw("Found absent col file for mandatory seq_region $sr_name:\t$merged_file");
        }
        
        next;
	  }
	   
	     
      my $efg_sr_id = $self->get_seq_region_id_by_Slice($slice);
      push @efg_sr_ids, $efg_sr_id;
      push @file_list, $slice_file; 
  
      #Set new byte off set for this slice
      $off_set           = $off_set + ($last_length);
      $index{$efg_sr_id} = $off_set;
  
      #Calculate last length
      $last_length     = $slice->length / $wsize;#rather than end, just in case we don't start at 1
      my $tmp_last_length = int($last_length);
      $last_length     = $tmp_last_length + 1 if ($last_length != $tmp_last_length);
      $last_length    *= $packed_size; 
      $total_length    = $off_set + $last_length;

      #Now validate file length
      #This should have already been done in the Collector
      $self->validate_file_length($slice_file, $last_length, 1);#1 is binmode;
    }
 
    #Could probably do this via the dump by ordering by seq_region_name/id
    #But mysql sort is not reliably the same as unix/perl sort
    #So let's do it here to be sure.
    
    print "Total $wsize window_size data length:\t\t$total_length\n";
  
    #Encode release version in the file name for visibility.
    #How big does the index need to be? I (32 bit unsigned) endian order depends on arch
    #v for size of index
    #v for index key i.e. sr_id
    #V for offset (4bytes?, max unsigned is 4,294,967,295)
    #Max offset is??
    #=>Nchrs * (max chr length) / (min window size) * packed_size
    #=> ~23(ignore non-chr for now) * 249250621 / 30 * 4 (also ignore index length)
    #=>764368568 < 4294967295
    #              404840556
    #Plenty of room for ~ 6* more data!
    
    #Encode index as key(sr_id - v 2 bytes) value(offset - V 4 bytes) pairs
  
    
    
    my $num_keys      = scalar(keys %index);
    my $index_size    = scalar(keys %index)*(6);
    #Full index is index_size + key-values pairs
    #Could add file length here to enable validation
    #i.e. seek to full length and make sure we don't get and EOF until we seek one byte more
  
    #perl -V:{short,int,long{,long}}size
    #shortsize='2';
    #intsize='4';
    #longsize='8';
    #longlongsize='8';
    #Actual size can change due to compiler
    #	   use Config;
    #       print $Config{shortsize},    "\n";
    #       print $Config{intsize},      "\n";
    #       print $Config{longsize},     "\n";
    #       print $Config{longlongsize}, "\n";
  
    
    #Does this also mean that changing data between perl 32bit -> 64bit will screw the reads?   
    
    my $total_index_size = 2 + ($num_keys * (2 + 4));
    
    #Adjust the offsets to account for the total index size
    foreach my $key(keys %index){
    	#warn "Adjusting $key index value from $index{$key}\t";
    	$index{$key} += $total_index_size;
    	#warn "to ".$index{$key};
    }
  
  
    my $pack_template = 'v(vV)'.$num_keys;
    #warn "pack template is $pack_template";
    my @index = %index;
    my $packed_index  = pack($pack_template, ($index_size, @index));
    
      
    #Let's just validate the index here
    #my @unpacked_index = unpack($pack_template, $packed_index);
    #my $unpacked_index_size = shift(@unpacked_index);
    #warn "Index size($index_size) unpacked is $unpacked_index_size";
    #my %unpacked_index = @unpacked_index;
    #warn "Keys ".$num_keys.' unpacked is '.scalar(keys(%unpacked_index));
    #for my $key(keys(%unpacked_index)){
    #  warn "key $key value ".$index{$key}.' unpacked is '.$unpacked_index{$key};
    #}
  
  
    #Dumps the index to file and add it to the start of the cat list
    my $idx_name   = join('.', ('result_features',  $rset->name, $wsize, 'idx')); 
    my $index_file = "${data_dir}/${idx_name}";
    print "Writing index:\t\t\t$idx_name\n";
    open($fh, '>', $index_file) or die("Cannot open $index_file");
    binmode($fh);
    print $fh $packed_index;
    close($fh);
  
  
    #Validate index file
    $fh           = $self->get_filehandle($index_file, {-binmode => 1});
    #Do this directly here as we have little need for grabbing the whole index normally
    my $index_ref = $self->{file_cache}{$index_file}{off_sets};
  
    #### Now validate written index
  
    #Test num keys is same
    if($num_keys != scalar(keys %{$index_ref})){
      throw("Original index has $num_keys whilst read and unpacked index has ".scalar(keys %{$index_ref}));
    }
    
    for my $key(keys(%{$index_ref})){ #Test values match
  	
      if($index{$key} != $index_ref->{$key}){
        throw("Original key $key value ".$index{$key}.
         ' does not match read and unpacked value '.$index_ref->{$key});
      }
    }
    
  
    #Generate merged file
    my $cat = "cat $index_file @file_list > $merged_file";
    print "Generating indexed file:\t$merged_file\n";
    system($cat) == 0 or die ("Cannot cat files:\n$cat");
  
    #Validate length
    $total_length += $total_index_size;
    $self->validate_file_length($merged_file, $total_length, 1);#1 is binmode;
  
    push @rm_files, ($index_file, @file_list);
  }
  
  		
  #Remove index and slice files
  #print STDOUTT "Removing index and sr_name col files\n";
  unlink(@rm_files) || throw("Failed to remove tmp files:\n\t".join("\n\t", @rm_files)."\n$!");
    
  return 1;
}

1;
