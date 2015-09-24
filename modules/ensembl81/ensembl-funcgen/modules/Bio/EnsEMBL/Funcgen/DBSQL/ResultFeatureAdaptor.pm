#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ResultFeatureAdaptor
#
#

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

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::ResultFeatureAdaptor - adaptor for fetching and storing ResultFeature objects

=head1 SYNOPSIS

my $rfeature_adaptor = $db->get_ResultFeatureAdaptor();

my @result_features = @{$rfeature_adaptor->fetch_all_by_ResultSet_Slice($rset, $slice)};


=head1 DESCRIPTION

The ResultFeatureAdaptor is a database adaptor for storing and retrieving
ResultFeature objects.

This will automatically query the web optimised result_feature
table if a data is present, else it will query the underlying raw data tables.


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


#How are we going to track the association between raw experimental_chip and collection result set?  
#ExperimentalChip/Channel ResultSets will probably go away, so ignore this association problem for now.
#Association between FeatureSet and Input/ResultSet to be handled in import pipeline

#ResultFeature ResultSets now have a RESULT_FEATURE_SET status entry.
#Add some details about file based collections


package Bio::EnsEMBL::Funcgen::DBSQL::ResultFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception  qw( throw deprecate );
use Bio::EnsEMBL::Utils::Scalar     qw( assert_ref );
use Bio::EnsEMBL::Funcgen::ResultFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor; #DBI sql_types import

use base qw( Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor
			       Bio::EnsEMBL::Funcgen::Collector::ResultFeature  );

#Private vars to used to maintain simple implementation of Collector
#Should be set in each method to enable trimmingof  the start and end bins. 
#Cannot depend on $dest_slice_start/end in _objs_from_sth
#As _collection_start/end are adjusted to the nearest bin
my ($_scores_field, $_collection_start, $_collection_end);
#my ($_window_size);#Need to be a method as it is used by the BaseFeatureAdaptor. or our?


=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  my $self = shift;

  return ([ 'result_feature', 'rf' ]);
}


=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  my $self = shift;

  return  ('rf.seq_region_start', 'rf.seq_region_end', 'rf.seq_region_strand', 
		   "$_scores_field", 'rf.result_set_id');
}




=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk - Moving to DBFile implementation

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;
 
  if(! $dest_slice){
	throw('ResultFeatureAdaptor always requires a dest_slice argument');
	#Is this correct?
	#ExperimentalChip based normalisation?
	#Currently does not use this method
	#Never have non-Slice based fetchs, so will always have dest_slice and seq_region info.
  }
  
  my (%rfeats, $start, $end, $strand, $scores, $rset_id);
  my $window_size = 0; #Always natural resolution from table

  #Could dynamically define simple obj hash dependant on whether feature is stranded and new_fast?
  #We never call _obj_from_sth for extended queries
  #This is only for result_feature table queries i.e. standard/new queries
  $sth->bind_columns(\$start, \$end, \$strand, \$scores, \$rset_id);


  #Test slice is loaded in eFG?
  my $slice_adaptor = $self->db->get_SliceAdaptor;
  
  if(! $slice_adaptor->get_seq_region_id($dest_slice)){
	warn "Cannot get eFG slice for:".$dest_slice->name.
	  "\nThe region you are using is not present in the current dna DB";
	return;
  }

  if($mapper) {
	throw('Cannot dynamically assembly map Collections yet');
	#See GeneAdaptor for mapping code if required
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;

  $dest_slice_start  = $dest_slice->start();
  $dest_slice_end    = $dest_slice->end();
  $dest_slice_strand = $dest_slice->strand();
  $dest_slice_length = $dest_slice->length();
  $dest_slice_sr_name = $dest_slice->seq_region_name();
  $dest_slice_sr_id = $dest_slice->get_seq_region_id();
    
  
  my (@scores, $slice, $start_pad, $end_pad);
  #Set up the %rfeats arrays here to prevent having to test in loop
  #This will speed up 0 wsize, but most likely slow others?
  
 FEATURE: while ( $sth->fetch() ) {

	if($window_size == 0){
	  warn "0bp window size array based result_features are no longer supported";
	  #Remove completely or re-instate?
	}
	else{

	  #Cannot have a collection which does not start at 1
	  #As we cannot compute the bin start/ends correctly?
	  #Actually we can do so long as they have been stored correctly
	  #i.e. start and end are valid bin bounds(extending past the end of the seq_region if needed)
	  #Let's keep it simple for now
	  throw("Collections with a window size > 0 must start at 1, not ($start)") if $start !=1;
	  #Can remove this if we test start and end are valid bin bounds for the given wsize

	  #Account for oversized slices
	  #This is if the slice seq_region_start/end are outside of the range of the record
	  #As collections should really represent a complete seq_region
	  #This should only happen if a slice is defined outside the the bounds of a seq_region
	  #i.e. seq_region_start < collection_start or seq_region_end > slice length
	  #if a test slice has been stored which does not represent the complete seq_region
	  	  #We don't need to pad at all, just adjust the $_collection_start/ends!!!
	  #Don't need to account for slice start < 1
	  #Start and end should always be valid bin bounds
	  #These could be removed if we force only full length seq_region collections
	  $_collection_start = $start if($_collection_start < $start);
	  $_collection_end   = $end   if($_collection_end   > $end);
	  
	  #warn "col start now $_collection_start";
	  #warn "col end   now $_collection_end";
	  
	  # If the dest_slice starts at 1 and is foward strand, nothing needs doing
	  # else convert coords
	  # These need to use $_collection_start/end rather than dest_slice_start/end

	  if($dest_slice_start != 1 || $dest_slice_strand != 1) {
		if($dest_slice_strand == 1) {
		  $start = $_collection_start - $dest_slice_start + 1;
		  $end   = $_collection_end   - $dest_slice_start + 1;
		} else {
		  my $tmp_seq_region_start = $_collection_start;
		  $start = $dest_slice_end - $_collection_end + 1;
		  $end   = $dest_slice_end - $tmp_seq_region_start + 1;
		  $strand *= -1;
		}	  
	  }
	  #What about 0 strand slices?
	  
	  #throw away features off the end of the requested slice or on different seq_region
	  if($end < 1 || $start > $dest_slice_length){# ||
		#( $dest_slice_sr_id ne $seq_region_id )) {
		#This would only happen if assembly mapper had placed it on a different seq_region
		#Dynamically mapped features are not guaranteed to come back in correct order?
		
		next FEATURE;
	  }
	
	  @scores = unpack('('.$self->pack_template.')'.(($_collection_end - $_collection_start + 1)/$window_size ), $scores);
	}


	push @{$rfeats{$rset_id}}, Bio::EnsEMBL::Funcgen::Collection::ResultFeature->new_fast({
																			  start  => $start,
																			  end    => $end, 
																			  strand =>$strand, 
																			  scores => [@scores], 
																			  #undef, 
																			  #undef, 
																			  window_size => $window_size,
																			  slice       => $dest_slice,
																			 });

  }
  

  
  return \%rfeats;
}
  

=head2 store

  Args[0]    : List of Bio::EnsEMBL::Funcgen::ResultFeature objects
  Args[1]    : Bio::EnsEMBL::Funcgen::ResultSet
  Args[2]    : Optional - Assembly to project to e.g. GRCh37
  Example    : $rfa->store(@rfeats);
  Description: Stores ResultFeature objects in the result_feature table.
  Returntype : None
  Exceptions : Throws if a List of ResultFeature objects is not provided or if
               any of the attributes are not set or valid.
  Caller     : General
  Status     : At Risk - Moving to DBFile implementation

=cut

sub store{
  my ($self, $rfeats, $rset, $new_assembly) = @_;

  #We can't project collections to a new assembly after they have been generated 
  #as this will mess up the standardised bin bounds.
  #Just project the 0 window size and then rebuild other window_sizes form there

  $self->set_collection_defs_by_ResultSet($rset);
  throw("Must provide a list of ResultFeature objects") if(scalar(@$rfeats == 0));
 
  #These are in the order of the ResultFeature attr array(excluding probe_id, which is the result/probe_feature query only attr))
  my $sth = $self->prepare('INSERT INTO result_feature (result_set_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, scores) VALUES (?, ?, ?, ?, ?, ?)');  
  my $db = $self->db();
  my ($pack_template, $packed_string);


 
  #my @max_allowed_packet = $self->dbc->db_handle->selectrow_array('show variables like "max_allowed_packet"'); 
  #warn "@max_allowed_packet";

 FEATURE: foreach my $rfeat (@$rfeats) {
    
    if( ! (ref($rfeat) && $rfeat->isa('Bio::EnsEMBL::Funcgen::Collection::ResultFeature'))) {
      throw('Must be a Bio::EnsEMBL::Funcgen::Collection::ResultFeature object to store');
    }
    
	if($rfeat->window_size == 0){
	  throw('Non 0bp window_size ResultFeatures cannot be stored in the result_feature table, write a col file instead');
	}


	#This is the only validation! So all the validation must be done in the caller as we are simply dealing with ints?
	#Remove result_feature_set from result_set and set as status?
	
	my $seq_region_id;
	($rfeat, $seq_region_id) = $self->_pre_store($rfeat, $new_assembly);

	next if ! $rfeat;#No projection to new assembly
	#Is there a way of logging which ones don't make it?

	#This captures non full length collections at end of seq_region
	$pack_template = '('.$self->pack_template.')'.scalar(@{$rfeat->scores});


	#Check that we have non-0 values in compressed collections
	if($rfeat->window_size != 0){
	
	  if(! grep { /[^0]/ } @{$rfeat->scores} ){
		warn('Collection contains no non-0 scores. Skipping store for '.
			 $rfeat->slice->name.' '.$rfeat->window_size." window_size\n");
		next;
	  }
	}

	
	$packed_string = pack($pack_template, @{$rfeat->scores});	
	
	$sth->bind_param(1, $rfeat->result_set_id, SQL_INTEGER);
	$sth->bind_param(2, $seq_region_id,        SQL_INTEGER);
    $sth->bind_param(3, $rfeat->start,         SQL_INTEGER);
    $sth->bind_param(4, $rfeat->end,           SQL_INTEGER);
	$sth->bind_param(5, $rfeat->strand,        SQL_INTEGER);
	$sth->bind_param(7, $packed_string,        SQL_BLOB);
	$sth->execute();
  }

  return $rfeats;
}



=head2 _list_dbIDs

  Description: Re-implementation of parent method, as we now store in flat files
               and do not have internal DB IDs ResultFeature 
  Returntype : None
  Exceptions : Warns
  Caller     : General
  Status     : At risk

=cut

sub _list_dbIDs {
  warn('_list_dbIDs is not appropriate for the ResultFeatureAdaptor as it is likely all data is now stored in flat files');
  return $_[0]->_list_dbIDs;
}


=head2 _window_size

  Args       : None
  Example    : my $wsize = $self->_window_size
  Description: Gets the window_size of the current ResultFeature query.
               This needs to be a method rather than just a private variable
               as it is used by the BaseFeatureAdaptor.
  Returntype : int
  Exceptions : None
  Caller     : Bio::EnsEMBL::BaseFeatureAdaptor
  Status     : At risk - ??? Is this required anymore

=cut


sub _window_size{
  my $self = shift;

  return $self->{'window_size'};
}




=head2 set_collection_config_by_Slice_ResultSets

  Args[0]    : Bio::EnsEMBL::Slice
  Args[1]    : ARRAYREF of Bio::EnsEMBL::Funcgen::ResultSet object
  Args[2]    : int - Maximum number of bins required i.e. number of pixels in drawable region
  Example    : $self->set_collection_defs_by_ResultSet([$rset]);
  Description: Similar to set_collection_defs_by_ResultSet, but used
               to set a config hash used for multi-ResultSet fetches.
  Returntype : None
  Exceptions : throws is args are not valid
               throws if supporting InputSet is not of type result (i.e. short reads import)
               throws if supporting InputSet format is not SEQUENCING (i.e. short reads import)
               throws if ResultSet is not and input_subset or experimental_chip based ResultSet (i.e. channel etc)
  Caller     : ResultFeatureAdaptor::fetch_all_by_Slice_ResultSets
  Status     : At Risk

=cut

sub set_collection_config_by_Slice_ResultSets{
  my ($self, $slice, $rsets, $max_bins, $window_size) = @_;
  assert_ref($slice, 'Bio::EnsEMBL::Slice');
  assert_ref($rsets, 'ARRAY', 'ResultSets');
  
  if(scalar(@$rsets) == 0){
	  throw('You must pass an ARRAYREF of Bio::EnsEMBL::ResultSet objects');
  }
  
  my ($wsize, $window_element, @rsets, %wsize_config);
  my ($collection_start, $collection_end);

  if($window_size && $max_bins){
	  warn "Over-riding max_bins with specific window_size, omit window_size to calculate window_size using max_bins";
  }

  my ($is_rf_set);
  my $rf_source = 'file';

  foreach my $rset(@{$rsets}){
  	$is_rf_set = $rset->has_status('RESULT_FEATURE_SET') || 0;
  	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
    my $table_name = $rset->table_name;
   
  	if(! ($table_name eq 'input_subset' || $is_rf_set)){
  	  warn("Skipping non-ResultFeature $table_name ResultSet:\t".$rset->name);
  	  next;
  	}
  
  	#Eventually will want to completely remove DB sets?
  	#Depends on how we handle 0bp window sets, BigWig?
  	#Leave 0bp sets in DB for v62 which means we need an extra status
  	#Assume all input_subset rsets are file based
  
  	#NOTE:
  	#Don't need to check for packed_size and packed_template differences as they are now the
  	#same for all collections i.e. float. These would need to be separate queries otherwise.
  	#This method makes assumptions about the window_sizes array structure
  	#If this is to change more the the 0 - 30 bp change below then the config hash generation
  	#needs to be reviewed
  	
   	if($table_name eq 'experimental_chip'){ #Array Intensities i.e. single float
  	  $Bio::EnsEMBL::Utils::Collector::window_sizes->[0] = 0;#Can have natural resolution for low density array data
  	}
  	elsif($table_name eq 'input_subset'){  
  	  $Bio::EnsEMBL::Utils::Collector::window_sizes->[0] = 30;
  	}
  	else{
      throw('Bio::EnsEMBL::Funcgen::Collector:ResultFeature does not support ResultSets of type'.$rset->table_name);
  	}

  	### SET ResultSet CONFIG BASED ON OPTIMAL WINDOW SIZE
  	# This sets all based on $rf_source=file
  	# For wsize==0 (experimental_chip), the rf_source is fixed after this loop
  	# This is done to prevent cyclical dependancy and improve cache checking speed
   
  	if( (defined $window_element ) &&
  		exists $wsize_config{$rf_source}{$Bio::EnsEMBL::Utils::Collector::window_sizes->[$window_element]} ){
  	  $wsize_config{$rf_source}{$Bio::EnsEMBL::Utils::Collector::window_sizes->[$window_element]}{'result_sets'}{$rset->dbID} = $rset;
  	}
  	else{ #We have not seen this wsize before
  
  
  	  #This is currently entirely based on the position of the first wsize.
  	  #We can't strictly rely that the same window_element will be optimal for each collection
  	  #However this will work if they are size ordered and only the first element changes e.g. 0 || 30
  	  #Will need to do for each set if window_sizes change
  	  
  	  if(! defined $window_element){
  			
    		my @sizes = @{$self->window_sizes};
    		#we need to remove wsize 0 if ResultSet was generated from high density seq reads
    		#0 should always be first
    
    		shift @sizes if ($sizes[0] == 0 && ($rset->table_name eq 'input_subset'));
    		$max_bins ||= 700;#This is default size of display?
            #CR Change to true drawable pixel width
    		
    		#The speed of this track is directly proportional
    		#to the display size, unlike other tracks!
    		#e.g
    		#let's say we have 300000bp
    		#700  pixels will use 450 wsize > Faster but lower resolution
    		#2000 pixels will use 150 wsize > Slower but higher resolution
    		
    		if(defined $window_size){
    		
    		  if(! grep { /^${window_size}$/ } @sizes){
      			warn "The ResultFeature window_size specifed($window_size) is not valid, the next largest will be chosen from:\t".join(', ', @sizes);
    		  }
    		  else{
      			$wsize = $window_size;
    		  }
    		}
    		else{#! defined $window_size
    		
    		  #Work out window size here based on Slice length
    		  #Select 0 wsize if slice is small enough
    		  #As loop will never pick 0
    		  #probably half 150 max length for current wsize
    		  #Will also be proportional to display size
    		  #This depends on size ordered window sizes arrays
    		  
    		  $window_size = ($slice->length)/$max_bins;
    		
    		  if($Bio::EnsEMBL::Utils::Collector::window_sizes->[0] == 0){
      			my $zero_wsize_limit = ($max_bins * $sizes[1])/2;
      			
      			if($slice->length <= $zero_wsize_limit){
      			  $wsize = 0;
      			  $window_element = 0;
      			}
    		  }
    		}
      
      
    		#Let's try and avoid this loop if we have already grep'd or set to 0
    		#In the browser this is only ever likely to speed up the 0 window
    		
    		if (! defined  $wsize) {
    		  #default is maximum
    		  $wsize = $sizes[-1];  #Last element
    		  
    		  #Try and find the next biggest window
    		  #As we don't want more bins than there are pixels
    
    		  for (my $i = 0; $i <= $#sizes; $i++) {
      			#We have problems here if we want to define just one window size
      			#In the store methods, this resets the wsizes so we can only pick from those
      			#specified, hence we cannot force the use of 0
      			#@sizes needs to always be the full range of valid windows sizes
      			#Need to always add 0 and skip_zero window if 0 not defined in window_sizes?
      		  
      			if ($window_size <= $sizes[$i]) {
      			  $window_element = $i;
      			  $wsize = $sizes[$i];
      			  last;    
      			}
    		  }
    		}
      }
      else{   #ASSUME the same window_element has the optimal window_size
    		$wsize = $Bio::EnsEMBL::Utils::Collector::window_sizes->[$window_element];
      }
  	  
  	 
  	  #Set BLOB access & collection_start/end config
  
  	  if ( $wsize == 0) {
    		$wsize_config{$rf_source}{$wsize}{scores_field} = 'rf.scores';
    		#No need to set start end config here as these are normal features
  	  } else {
    		#We want a substring of a whole seq_region collection
    	  
    		#Correct to the nearest bin bounds	  
    		#int rounds towards 0, not always down!
    		#down if +ve or up if -ve
    		#This causes problems with setting start as we round up to zero
    		
    
    		#Sub this lot as we will use it in the Collector for building indexes?
    
    
    		my $start_bin      = $slice->start/$wsize;
    		$collection_start = int($start_bin);
    		
    		#Add 1 here so we avoid multiply by 0 below
    		if ($collection_start < $start_bin) {
    		  $collection_start +=1;	#Add 1 to the bin due to int rounding down
    		}
    		
    		$collection_start = ($collection_start * $wsize) - $wsize + 1 ; #seq_region
    		$collection_end   = int($slice->end/$wsize) * $wsize; #This will be <= $slice->end
    		
    		#Add another window if the end doesn't meet the end of the slice
    		if (($collection_end > 0) &&
    			($collection_end < $slice->end)) {
    		  $collection_end += $wsize;
    		}
  	
    		#Now correct for packed size
    		#Substring on a blob returns bytes not 2byte ascii chars!
    		#start at the first char of the first bin
    		my $sub_start = (((($collection_start - 1)/$wsize) * $self->packed_size) + 1); #add first char
    		#Default to 1 as mysql substring starts < 1 do funny things
    		$sub_start = 1 if $sub_start < 1;
    		my $sub_end   = (($collection_end/$wsize) * ($self->packed_size));
    		
    		#Set local start end config for collections
    		$wsize_config{$rf_source}{$wsize}{collection_start} = $collection_start;
    		$wsize_config{$rf_source}{$wsize}{collection_end}   = $collection_end;
    
    
    		if($rf_source eq 'file'){  #file BLOB access config
    		  $wsize_config{$rf_source}{$wsize}{'byte_offset'} = $sub_start -1;   #offset is always start - 1
    		  $wsize_config{$rf_source}{$wsize}{'byte_length'} = ($sub_end - $sub_start + 1);
    		  #warn "byte_offset = $sub_start";
    		  #warn "byte_length = ($sub_end - $sub_start + 1) => ". ($sub_end - $sub_start + 1);
    		}
    		else{ #DB BLOB access config
    		  #Finally set scores column for fetch
    		   $wsize_config{$rf_source}{$wsize}{scores_field} = "substring(rf.scores, $sub_start, ".
    			 ($sub_end - $sub_start + 1).')';
    		}
  	  }
  
  	  #Set the result_set and scores field config
  	  #Would also need to set pack template here if this 
  	  #were to change between collections
  	  $wsize_config{$rf_source}{$wsize}{result_sets}{$rset->dbID} = $rset;
  	}
  	
  }

  #Fix the 0 wsize source as will be set to 'file' but should be 'db'
  
  if(exists $wsize_config{'file'}{0}){
  	$wsize_config{'db'}{0} = $wsize_config{'file'}{0};
  	delete $wsize_config{'file'}{0};
  
  	if(keys(%{$wsize_config{'file'}}) == 0){
  	  delete $wsize_config{'file'};
  	}
  }

  $self->{'_collection_config'} = \%wsize_config;
  return $self->{'_collection_config'};
}





=head2 fetch_all_by_Slice_ResultSets

  Arg[1]     : Bio::EnsEMBL::Slice - Slice to retrieve results from
  Arg[2]     : ARRAYREF of Bio::EnsEMBL::Funcgen::ResultSets - ResultSet to retrieve results from
  Arg[3]     : OPTIONAL int    - max bins, maximum number of scores required
  Arg[4]     : OPTIONAL int    - window size, size of bin/window size in base pirs
  Arg[5]     : OPTIONAL string - sql contraint for use only with DB collections
  Arg[6]     : OPTIONAL hasref - config hash
  Example    : my %rfeatures = %{$rsa->fetch_all_by_Slice_ResultSets($slice, [@rsets])};
  Description: Gets a list of lightweight ResultFeature collection(s) for the ResultSets and Slice passed.
  Returntype : HASHREF of ResultSet dbID keys with a LISTREF of ResultFeature collection values
  Exceptions : Warns and skips ResultSets which are not RESULT_FEATURE_SETS.
  Caller     : general
  Status     : At risk

=cut

sub fetch_all_by_Slice_ResultSets{
  my ($self, $slice, $rsets, $max_bins, $window_size, $orig_constraint) = @_;
 
  $orig_constraint .= (defined $orig_constraint) ? ' AND ' : '';
  #currently does not accomodate raw experimental_chip ResultSets!
  my $conf = $self->set_collection_config_by_Slice_ResultSets($slice, $rsets, $max_bins, $window_size);

  #Loop through each wsize build constraint set private vars and query
  my (%rset_rfs, $constraint);

  #Can remove this DB conf when we remove support for 0 wsize experimental_chip
  #result_features

  my $rf_conf = $conf->{db};
  
  foreach my $wsize(keys(%{$rf_conf})){
    $self->{'window_size'} = $wsize;
    $_scores_field          = $rf_conf->{$wsize}->{scores_field};
    $constraint = 'rf.result_set_id IN ('.join(', ', keys(%{$rf_conf->{$wsize}->{result_sets}})).')'.
      " AND rf.window_size=$wsize";

    if ($wsize != 0){
      $_collection_start = $rf_conf->{$wsize}->{collection_start};
      $_collection_end   = $rf_conf->{$wsize}->{collection_end};
    }
	
    my ($rset_results) = @{$self->fetch_all_by_Slice_constraint($slice, $orig_constraint.$constraint)};
    	  	  
    #This maybe undef if the slice is not present
    $rset_results ||= {};
    
    #Will this work for wsize split queries?
    #We are not setting values for empty queries
    #Which is causing errors on deref
    	
    %rset_rfs = (%rset_rfs, %{$rset_results});
    
    #Account for DB rsets which return no features
    foreach my $rset(values %{$rf_conf->{$wsize}{result_sets}}){
    
      if(! exists $rset_rfs{$rset->dbID}){
        $rset_rfs{$rset->dbID} = [];
      }
    }
  }


  #Two blocks to prevent extra if-conditional for each rset.
  #Can easily remove one block if/when we deprecate DB collections
  #Should always be the same wsize if we are using a file i.e. no 0bp
  #Unless we have differing window_size ranges
  my $rf;
  $rf_conf = $conf->{file};

  foreach my $wsize(keys(%{$rf_conf})){

    foreach my $file_rset(values %{$rf_conf->{$wsize}{result_sets}}){
      $rf = $self->_fetch_from_file_by_Slice_ResultSet($slice, $file_rset, $wsize, $rf_conf);
      $rset_rfs{$file_rset->dbID} = defined($rf) ? [$rf] : [];
    }
  }

  return \%rset_rfs;
}



=head2 _fetch_from_file_by_Slice_ResultSet

  Arg[1]     : Bio::EnsEMBL::Slice - Slice to retrieve results from
  Arg[2]     : Bio::EnsEMBL::Funcgen::ResultSets - ResultSet to retrieve results from
  Arg[3]     : int     - valid window size defined by set_collection_config_by_ResultSets
  Arg[4]     : HASHREF - Window size config for file based ResultSets
  Example    : my $rfeat = $self->_fetch_from_file_by_Slice_ResultSet($slice, $rset, $wsize, $conf);
  Description: Generates ResultFeature collection reading packed scores from a 'col' file.
  Returntype : Bio::EnsEMBL::Funcgen::Collection::ResultFeature
  Exceptions : 
  Caller     : fetch_all_by_Slice_ResultSets
  Status     : At risk

=cut


#Set/store filepath in ResultSet to avoid having to regenerate?
#TODO change the start/end to local coords wrt query slice. How does the webcode handle this at present?

sub _fetch_from_file_by_Slice_ResultSet{
  my ($self, $slice, $rset, $window_size, $conf) = @_;
  #private as window_size needs to ba valid i.e. generated by set_collection_config_by_ResultSets
  #and no class tests
  
  #Cache this as  ResultSet::get_dbfile_path_prefix?
  #Is there any point if the rsets aren't cached?
  #Data is cached anyway, so no redundant calls
  #ResultSets are always regenerated
  #Key would either have to be query or dbID
  #Former would be hard?(constraint key, is this already done for features?)
  #Later would be object cache, so we would still do sql query but skip the object generation 
  #given the dbID is enough to pull a valid object from the cache

  #How can we cache result/feature sets?
  #Would this really speed anything up?
  
  #if(! exists $conf->{$window_size}){
  #    throw("Specified window_size($window_size) is not present in config.\n".
  #		  "Have you called this private method directly?\n".
  #		  "Try using the fetch_by_Slice_ResultSets warpapper method\n".
  #		  "Or set the RESULT_FEATURE_FILE_SET status and window_size correctly.");
  #  }

  my $rf;
  my $efg_sr_id = $self->get_seq_region_id_by_Slice($slice);

  if($efg_sr_id){
    my $packed_scores =  $self->read_collection_blob
	    ($rset->get_dbfile_path_by_window_size($window_size),
	     #Would be in analysis object for unique analysis tracks/data
	     $efg_sr_id,
	     $conf->{$window_size}{'byte_offset'},
	     $conf->{$window_size}{'byte_length'});

	my ($start, $end, @scores);
	

	if(defined $packed_scores){
	  ($start, $end) = ($conf->{$window_size}{collection_start}, 
						$conf->{$window_size}{collection_end});
	  
	  #Need to capture unpack failure here and undef the fh?
	  #i.e. pack/unpack repeat count overflow 
	  
	  @scores = unpack('('.$self->pack_template.')'.(($end - $start + 1)/$window_size),
					   $packed_scores);
	  
	  #could validate scores size here
	  #will only ever be <= than expected value
	  #as unpack will discard excess
	  #better validate length of $packed_scores
	  
	  $rf =  Bio::EnsEMBL::Funcgen::Collection::ResultFeature->new_fast
		({
		  start  => $start,
		  end    => $end,
		  strand => 0,           #These are strandless features
		  scores => [@scores],
		  window_size => $window_size,
		 });
	}
  }
	
  return $rf;
}



=head2 fetch_all_by_Slice_ResultSet

  Arg[1]     : Bio::EnsEMBL::Slice - Slice to retrieve results from
  Arg[2]     : Bio::EnsEMBL::Funcgen::ResultSet - ResultSet to retrieve results from
  Arg[3]     : OPTIONAL int    - max bins, maximum number of scores required
  Arg[4]     : OPTIONAL int    - window size, size of bin/window size in base pirs
  Arg[5]     : OPTIONAL string - sql contraint for use only with DB collections
  Arg[6]     : OPTIONAL hasref - config hash
  Example    : my %rfeatures = %{$rsa->fetch_all_by_Slice_ResultSets($slice, [@rsets])};
  Description: Gets a list of lightweight Collection of ResultFeatures  for the ResultSet and Slice passed.
               NOTE: ExperimentalChip/Channel based ResultFeatures was removed in version 63.
  Returntype : Bio::EnsEMBL::Funcgen::Collection::ResultFeature
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut



#To do 
#remove Bio::EnsEMBL::Funcgen::ResultFeature in favour of Collection::ResultFeature?

sub fetch_all_by_Slice_ResultSet{
  my ($self, $slice, $rset, $max_bins, $window_size, $orig_constraint) = @_;
  #Do this first to avoid double validation of rset
  my $rf_hashref = $self->fetch_all_by_Slice_ResultSets($slice, [$rset], $max_bins, $window_size, $orig_constraint);

  return $rf_hashref->{$rset->dbID};
}



sub fetch_Iterator_by_Slice_ResultSet{
  my ($self, $slice, $rset, $max_bins, $window_size, $constraint, $chunk_size) = @_;

  return $self->fetch_collection_Iterator_by_Slice_method
	($self->can('fetch_all_by_Slice_ResultSet'),
	 [$slice, $rset, $max_bins, $window_size, $constraint],
	 0,#Slice idx
	 $chunk_size #Iterator chunk length
	);

}

=head2 fetch_collection_Iterator_by_Slice_method

  Arg [1]    : CODE ref of Slice fetch method
  Arg [2]    : ARRAY ref of parameters for Slice fetch method
  Arg [3]    : Optional int: Slice index in parameters array
  Arg [4]    : Optional int: Slice chunk size. Default=500000
  Example    : my $slice_iter = $feature_adaptor->fetch_Iterator_by_Slice_method
                               	      ($feature_adaptor->can('fetch_all_by_Slice_Arrays'),
	                                   \@fetch_method_params,
	                                   0,#Slice idx
	                                   #500 #chunk length
	                                  );

               while(my $feature = $slice_iter->next && defined $feature){
                 #Do something here
               }

  Description: Creates an Iterator which chunks the query Slice to facilitate
               large Slice queries which would have previously run out of memory
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : Throws if mandatory params not valid
  Caller     : general
  Status     : at risk - move to core BaseFeatureAdaptor

=cut


#Essentially the only difference here we have one feature with an array of 'scores'

sub fetch_collection_Iterator_by_Slice_method{
  my ($self, $slice_method_ref, $params_ref, $slice_idx, $chunk_size) = @_;

  if(! ( defined $slice_method_ref &&
		 ref($slice_method_ref) eq 'CODE')
	){
	throw('Must pass a valid Slice fetch method CODE ref');
  }

  if (! ($params_ref && 
		 ref($params_ref) eq 'ARRAY')) {
	#Don't need to check size here so long as we have valid Slice
	throw('You must pass a method params ARRAYREF');
  }
  
  $slice_idx    = 0 if(! defined $slice_idx);
  my $slice     = $params_ref->[$slice_idx];
  $chunk_size ||= 1000000;
		
  my $collection;
  my $finished     = 0;
  my $start        = 1;	#local coord for sub slice
  my $end          = $slice->length;
  my $overlap      = 0;
  
  my $coderef = 
	sub {
	  my $collection;

	  if(! $finished) {
		
		my $new_end = ($start + $chunk_size - 1);
		
		if ($new_end >= $end) {
		  # this is our last chunk
		  $new_end = $end;
		  $finished = 1;  
		}
		
		#Chunk by sub slicing
		my $sub_slice             = $slice->sub_Slice($start, $new_end);
		$params_ref->[$slice_idx] = $sub_slice;
		($collection) = @{ $slice_method_ref->($self, @$params_ref)};
		
		
		if(! $collection){
		  $finished = 1;
		}
		else{
		  
		  #Trim score and start if overlapping found
		  #Will only ever be on overlapping 'score'
		  if($overlap){
			shift @{$collection->scores};
			$collection->{start} = $collection->start +  $collection->window_size;
			$overlap = 0;
		  }
		  
		  
		  if ( $collection->end > $sub_slice->end ){
			$overlap = 1;
		  }
		  
		  $start = $new_end + 1;
		}
	  }
	  
	  return $collection;
	};

  return Bio::EnsEMBL::Utils::Iterator->new($coderef);
}



# Over-ride/deprecate generic methods
# which do not work with ResultFeature Collections

sub fetch_all{
  deprecate('The fetch_all method has been disabled as it is not appropriate for the ResultFeatureAdaptor');
  return;
}

sub fetch_by_dbID{
  warn 'The fetch_by_dbID method has been disabled as it is not appropriate for the ResultFeatureAdaptor';
  #Could use it for 0 wsize DB based data, but not useful.
  return;
}

sub fetch_all_by_dbID_list {
  warn 'The fetch_all_by_dbID_list method has been disabled as it is not appropriate for the ResultFeatureAdaptor';
  #Could use it for 0 wsize DB based data, but not useful.
  return;
}

sub fetch_all_by_logic_name {
  warn 'The fetch_all_by_logic_name method has been disabled as it is not appropriate for the ResultFeatureAdaptor';
  #Could use it for 0 wsize DB based data, but not useful.
  return;
}

sub _list_seq_region_ids{
 warn 'The _list_seq_region_ids method has been disabled as it is not appropriate for the ResultFeatureAdaptor';
 #Could use it for 0 wsize DB based data, but not useful.
 return
}

#Over-ride fetch_all_by_display_label? Or move this to the individual FeatureAdaptors?
#Same with fetch_all_by_stable_Storable_FeatureSEts and wrappers (fetch_all_by_external_name)?
#Basically this is not a DBAdaptor anymore so should inherit from somewhere else. 
#Need to separate the common utility methods and have co-inheritance e.g.
#DBFile::Adaptor    Utils::FeatureAdaptor
#DBAdaptor::Adaptor Utils::FeatureAdaptor



1;


