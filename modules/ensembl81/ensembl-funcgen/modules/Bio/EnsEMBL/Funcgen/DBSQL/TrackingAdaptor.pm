#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::TrackingAdaptor
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::TrackingAdaptor - Handles interactions with the funcgen tracking DB

=head1 SYNOPSIS

my $tdb = Bio::EnsEMBL::DBSQL::Funcgen::DBSQL::TrackingAdaptor->new
  (
   -user     => $tdb_user,
   -pass     => $tdb_pass,
   -dbname   => $tdb_name,
   -host     => $tdb_host,
   -species  => $species,
  );


=head1 DESCRIPTION

Handles interaction with the funcgen tracking DB is used to store meta data and tracking information
about data sets. This provides a staging area where analysis and QC can be performed prior to
migration to the DB destined for a particular release. This adaptor methods to store and access
these information given the appropriately related objects e.g.
    InputSet, InputSubset, ResultSet, FeatureSet & DataSet

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;

use strict;
use warnings;
use feature qw(say);

use DateTime;
#use POSIX                                  qw( strftime );
use Bio::EnsEMBL::Utils::Exception         qw( throw warning );
use Bio::EnsEMBL::Utils::Scalar            qw( check_ref assert_ref );
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( get_month_number split_CamelCase);

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );
#Currently don't use and Funcgen BaseAdaptor methods in here
#and this would require _true_tables to be defined

#values are mandatory booleans
my %tracking_columns =
 (input_subset => {availability_date => 1,
                   download_url      => 1,
                   download_date     => 0,
                   local_url         => 0,
                   md5sum            => 0,
                   notes             => 0},
  result_set   => {idr_max_peaks        => 0,
                   idr_peak_analysis_id => 0},
  experiment   => {experiment_id => 0,
                   notes             => 0},
 );


=head2 new

  Arg [-DNADB_HOST]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_USER]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_PASS]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_PORT]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_NAME]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_ASSEMBLY] : String - Overrides defaults (ensembldb or registry)

  Arg [...]         : Other args are passed to superclass Bio::EnsEMBL::DBSQL::DBAdaptor
  Example1          : $db = new Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
                              (
						       -user   => 'readonly', #No password
						       -dbname => 'pog',
						       -host   => 'caldy',
						    );

                      #If dnadb is not defined in registry, this will automatically
                      #set it from the default dnadb_host (e.g. ensembldb)

  Exmaple2          : $db = new Bio::EnsEMBL::DBSQL::Funcgen::DBAdaptor
                              (
						       -user           => 'write',
                               -pass           => 'password',
						       -dbname         => 'pog',
						       -host           => 'caldy',
                               -dnadb_assmebly => '36',
						      );
                       #This will specifically look for a dnadb with assembly version 36
                       #on the default dnadb_host

   Exmaple2          : $db = new Bio::EnsEMBL::DBSQL::Funcgen::DBAdaptor
                              (
						       -user           => 'write',
                               -pass           => 'password',
						       -dbname         => 'pog',
						       -host           => 'caldy',
                               #The following will over-ride the default dnadb setting
                               -dnadb_name     => 'my_homo_sapiens_core_67_37',
                               -dnadb_host     => 'my_host',
                               #Can add more dnadb params if required
                               #-dnadb_host
                               #-dnadb_user
                               #-dnadb_port
						      );

  Description: Constructor for DBAdaptor. Will automatically set the dnadb based on dnadb params.
               This makes some assumptions about how the DBs name are defined i.e. the last part must
               conform to the _RELEASE_ASSEMBLY format e.g. 67_37 in homo_sapiens_funcgen_67_37
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions : Throws if conflicting dnadb params found
  Caller     : general
  Status     : Stable

=cut


# We don't want this to be available via $efgdb->get_TrackingAdaptor
# So have to handle setting db
# alternatively just pass DBAdaptor here?

sub new {
  my ($class, @args) = @_;
  return $class->SUPER::new(Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(@args));
}

sub repository_path{
  my ($self, $rpath) = shift;

  #Use MetaContainer instead
  my $sql;

  if(! defined $self->{repository_path}){
    $sql = "SELECT meta_value from meta where meta_key='repository';";
    ($self->{repository_path}) = $self->db->dbc->db_handle->selectrow_array($sql);
  }

  if(defined $rpath){

	if(defined $self->{repository_path} &&
	   ($self->{repository_path} ne $rpath ) ){
	  throw("Supplied repository path does not match meta_value: $rpath ne ".$self->{repository_path});
	}
	elsif(! defined $self->{repository_path}){
	  $sql = "INSERT into meta(meta_key, meta_value) values('repository_path', '${rpath}'";
	  $self->db->dbc->do($sql);
	}
  }

  return $self->{repository_path};
}


#Can take InputSet or InputSubset
#Change this to take sets rather than one?

#todo update this wrt local_url
#what are force_download, data, and skip_beds for
#surely that should be handled in a wrapper script?

#download_input_set_data will need revising


#Remove InputSet support from this to simplify it
#why are we still getting two records for a single InputSubset?
#we are still having to contend with merged input_subsets here
#we really need to fix this?

sub _columns{
  my $self       = shift;
  #unlikely table_name will be 0
  my $table_name = shift || throw('Must defined a table_name argument');

  if(! defined $table_name){
    throw('Must defined a table_name argument');
  }

  return keys %{$tracking_columns{$table_name}};
}

sub _mandatory_columns{
  my $self       = shift;
  #unlikely table_name will be 0
  my $table_name = shift || throw('Must defined a table_name argument');

  if(! defined $table_name){
    throw('Must defined a table_name argument');
  }

  #force list context so we don't just get the number of matches
  return (grep { $tracking_columns{$_} } keys %tracking_columns);
}

sub _is_mandatory_column{
  my $self       = shift;
  #unlikely table_name or col will be 0
  my $table_name = shift || throw('Must defined a table_name argument');
  my $col        = shift || throw('Must defined a columns argument');

  if(! exists $tracking_columns{$table_name}{$col}){
    throw($col." is not a valid column, must be one of:\t".join(' ', $self->_columns($table_name)));
  }

  return $tracking_columns{$table_name}{$col};
}


sub _is_column{
  my $self       = shift;
  #unlikely table_name or col will be 0
  my $table_name = shift || throw('Must defined a table_name argument');
  my $col        = shift || throw('Must defined a columns argument');

  return exists $tracking_columns{$table_name}{$col};
}


#optionally add in class specific subs like is_embargoed and set_dowload_status?
#As these methods are not simple getter/setters it maybe cleaner not to inject them
#But this leaves some methods available via the object, and some by the adaptor only :/
#This is a case to have separate class specific _inject methods, which call the generic
#_inject methods(for the getter/setters), then inject the additionally inject the
#class specific ones
#One downside to this is error handling, where line numbers will reference anonymous subroutines
#and not the true line number in this module. This is likely fine for the getter/setters, but might
#be problemtic with more complicated methods.

sub _inject_tracking_methods{
  my $self       = shift;
  my $storable   = shift || throw('Must pass a Storable argument');
  my $table_name = $self->get_valid_stored_table_name($storable);

  #We don't need to test for existing methods in the namespace here(CvGV_name_or_bust)
  #as we know exactly what methods are already present
  my @cols = $self->_columns($table_name);

  if(! $storable->can($cols[0])){
    my $ref = ref($storable);
    no strict 'refs';

    for my $col(@cols){
      *{$ref."::${col}"} =
        sub {
          my $self = shift;
          my $val  = undef;

          #to prevent auto vivification
          #which may prevent fetch
          if(exists $self->{tracking_info}){
            $val = $self->{tracking_info}{$col};
          }
          else{
            $self->throw("Cannot get $col attribute as tracking_info has not be set yet.".
            ' You need to call $tdb->fetch_tracking_info($storable) first.');
          }

          return $val;
        };
    }

    use strict;
    #avoids being able to work with symbolic reference
    #i.e. Your::Package::Name->$param_method
    #rather than $obj->$param_method
  }

  return;
}


#Todo
#Handle clashes between existing tracking_info attribute and fetched data

sub fetch_tracking_info{
  my $self       = shift;
  my $storable   = shift || throw('Must pass a Storable argument');
  my $table_name = $self->get_valid_stored_table_name($storable);
  my $db         = $self->db;

  $self->_inject_tracking_methods($storable);
  my $hashref;

  if(exists $storable->{tracking_info}){
    throw("$table_name Storable(".$storable->dbID.
      ") already has tracking_info attribute. Please store/update before fetching\n");

      #This may over-write existing data
  }else{

    my $sql = 'SELECT '.join(', ', ($self->_columns($table_name))).
              " FROM ${table_name}_tracking WHERE ${table_name}_id =".$storable->dbID;
    my $sth = $self->prepare($sql);
    $sth->execute;
    $hashref = $sth->fetchrow_hashref; #Will this be undef if there is no data?
    $storable->{tracking_info} = $hashref;
    $sth->finish;#otherwise we get disconnect warnings

    #my %column;
    #pseudo array/hash?
    #$sth->bind_columns( \( @column{ @{$sth->{NAME_lc} } } ));
    #while( $sth->fetch ){
    #  my $record = {%column}; #deref properly as %columns will be updated
    #  push @tracking_info, {%$record}; #Keep the dbID here
    #  #otherwise there will be now way to identify what the record refers to
    #  my $dbID = delete($record->{input_subset_id}); #Don't need this
    #  push @{$subset_cache{$dbID}->{tracking_info}}, $record;
    #}
  }

  return $hashref;
}


#TO DO
#This should also set the local_url?
#The is a question of validation of the tracking info before storing
#in this case, if we have set the download date, then we should also set the local_url
#is this a case to have a validate_input_subset_tracking_info method?
#This is currently of very limited utility, and can probably be handled in the caller (for now)
#Remove this in favour of store_tracking_info($storable, {allow_update =>1}) ?
#Although the is extra logic here: to_null and potential local_url check?

#This could be a wrapper to store_tracking_info. But we would need to get the stored data first
#and selectively overwrite just the date (and local_url) and set purge, to ensure that the NULL date
#is updated

sub set_download_status_by_input_subset_id{
  my ($self, $iss_id, $to_null, $lurl) = @_;



  my $date = ($to_null) ? 'NULL' : 'NOW()';
  my $sql = "UPDATE input_subset_tracking SET download_date=${date}";
  $sql .= ", local_url='${lurl}'" if $lurl;
  $sql .= " WHERE input_subset_id=${iss_id};";
  $self->db->dbc->do($sql);

  return;
}



#Careful this is american format for some reason!!!
#Can we inject a version this for input_subsets?

sub is_InputSubset_embargoed {
  my ($self, $isset, $yyyy_mm_dd) = @_;

  my $now_date = DateTime->today();
  my $rel_date;

  if(! defined $yyyy_mm_dd){
    $rel_date = $now_date;
  }
  else{

    my ($y, $m, $d) = split(/-/, $yyyy_mm_dd);

    if(! (defined $m && defined $d)){ #we have just passed a month number or abbrv
      $m = $y;
      undef $y;

      if($m !~ /^[0-9]{1,2}$/){
        my $month_number = get_month_number($m);

        if(! $month_number){
          throw("Unable to identify month number from:\t$m\n".
            'Please specify a two digit number (e.g. 01) rf a three letter code (e.g. Jan)');
        }

        $m = $month_number;
      }


      #we need to figure out whether the month is this year or next year
      $y = $now_date->year;

      if($m < $now_date->month){
        $y++;
      }

      $rel_date = DateTime->new(month => $m, year  => $y);
    }
  }

  #InputSubsets can have more than one tracking record at the moment
  #as they are currently merged at the InputSubset level rather than the InputSet level
  #What will happen here is we have other unstored tracking info
  #It will still barf.
  #we need to fix fetch_tracking_info to handle this, similar to how store does
  #fetching without method ssets means we should preferentially take the new data
  #what about NULLs? Should warn too
  #If the methods are available, it means we have already fetched
  #only if they are not available and tracking_info hash is present is it an issue
  #This should suffice for now
  $self->fetch_tracking_info($isset) if ! $isset->can('availability_date');

  #This was causing issues as DateTime cannot take the default null date of 0000-00-00
  my $avail_date = $isset->availability_date;

  if($avail_date eq '0000-00-00'){
    warn "Found $avail_date availability date for InputSubset:\t".$isset->name.
      "\nDefaulting to 0001-01-01\n";
    $avail_date = '0001-01-01';
  }

  my ($year, $month, $day) = split(/-/, $avail_date);
  my $isset_date           = DateTime->new(day   => $day,
                                           month => $month,
                                           year  => $year  );
  return ($isset_date > $rel_date) ? 1 : 0; #Nice operator overloading!
}



#Subroutine "store_tracking_info" with high complexity score (28) Consider refactoring.  (Severity: 3)

sub store_tracking_info{
  my $self         = shift;
  my $storable     = shift;
  my $params       = shift || {};
  my $table_name   = $self->get_valid_stored_table_name($storable);
  my $allow_update = (exists $params->{allow_update}) ? $params->{allow_update} : 0;
  my $purge        = (exists $params->{purge})        ? $params->{purge}        : 0;
  my $info         = (exists $params->{info})         ? $params->{info}         : undef;


  if(! defined $info){
    if(! exists $storable->{tracking_info}){
      throw('Found no tracking info to store. Please set info using tracking methods '.
      'if they are available, else pass the info hash parameters, but not both');
    }

    $info = $storable->{tracking_info};
  }
  elsif(exists $storable->{tracking_info}){
    #Does this risk unecessary failure if the tracking info
    #has has been autovivified for some reason with null values?
    #This should never happen.
    throw('An info hash parameter has been passed for a '.$table_name.
      " Storable which already has a tracking_info hash set\n".
      'Please use the existing tracking methods if available else pass the info hash parameter, but not both');
  }

  assert_ref($info, 'HASH', 'Tracking info HASH');

  #Grab the existing data to ensure we return/set the complete set of tracking
  #info after an update
  my $stored_info = $self->fetch_tracking_info($storable);

  #Compare hashes to avoid unecessary INSERT which may
  #fail if we have not specified allow_update
  my $info_is_identical = 0;

  if(keys %$stored_info){
     $info_is_identical = 1;
    foreach my $key(keys %$stored_info){ #This will have all the column keys
      if( ( (defined $stored_info->{$key}) &&
             ((! exists $info->{$key}) || ($stored_info->{$key} ne $info->{$key})) ) ||
          ( (! defined $stored_info->{$key}) &&
             (exists $info->{$key}) && (defined $info->{$key})) ){
        #This last test does not handle empty strings
        $info_is_identical = 0;;
      }
    }
  }
  #Test for unexpect info items
  my @valid_cols = $self->_columns($table_name);

  foreach my $col(keys %$info){

    if(! $self->_is_column($table_name, $col) ){
      #Could also grep this from @valid_cols?
      throw("Found unexpected parameter in tracking info hash:\t".$col.
        "Must be one of:\t@valid_cols");
    }
  }
  return if $info_is_identical ;

  #Test for mandatory info items, and build cols/values
  my @cols       = ($table_name.'_id');
  my @values     = ($storable->dbID);

  foreach my $col(@valid_cols){

    #Need to handle update here. There is no way of knowing whether the record has
    #already been stored, so we just skip this test if $allow_update is set
    #and let the eval handle the failure

    if($self->_is_mandatory_column($table_name, $col) &&
       ((! exists $info->{$col}) || (! defined $info->{$col}))){
      throw("Mandatory tracking column must be defined:\t$col\n$table_name:\t".$storable->name);
    }

    #This only updates those attributes which are set unless purge is set.
    #In which case undef attribs will be autovivified as undef(or NULL/default in the DB)

    #Hmm, this will not update current update NULLs!
    #e.g. when we want to reset a download_date to NULL
    #This makes the line between overwrite and update a bit murky
    #If we allow NULLs to be specified then we risk overwriting the data
    #especially if the tracking hash has been autovivified
    #This is only done after a fetch, so the data always matches the DB
    #The only risk is that someone may use the _columns method
    #to generate the info hash and leave some keys undef, which maybe defined in the DB
    #Even grabbing the store_info in set_download_status_by_input_subset_id
    #wont protect again this? What if we fetch the stored_info first
    #create the full valid hash, then specify purge in set_download_status_by_input_subset_id

    #Never update NULLs here unless $purge is set!

    if(defined $info->{$col} || $purge){
      push @cols,   $col;
      push @values, $info->{$col};
    }
  }

  my $sql = "INSERT IGNORE into ${table_name}_tracking(".join(', ', @cols).
    ') values("'.join('", "', @values).'")';

  if($allow_update){
    $sql .= 'ON DUPLICATE KEY UPDATE '.join(', ', (map {"$_=values($_)"} @cols));
  }

  #Use SQLHelper::execute_update here?
  #Although working with DBConnection here provides error handling
  my $row_cnt;
  my $db = $self->db;
  if(! eval { $row_cnt = $db->dbc->do($sql); 1 }){
    my $update_warn = '';

    if($allow_update){
      $update_warn = 'This maybe because allow_update is specified for an absent row '.
        "and the mandatory columns have not been specified:\n\t".
        join(' ', $self->_mandatory_columns)."\n";
    }

    throw("Failed to insert using sql:\t$sql\n${update_warn}$@");
  }

  #This shoudl always be greater than 1
  #$row_cnt = 0 if $row_cnt eq '0E0';

  #Overwrite set info KV pairs in stored_info
  #so we're sure have the complete set of tracking data
  #if allow_update was used with incomplete data.

  foreach my $key(keys %$info){
    $stored_info->{$key} = $info->{$key};
  }

  $self->_inject_tracking_methods($storable);
  $storable->{tracking_info} = $stored_info;

  return $storable;
}


sub get_valid_stored_table_name{
  my $self     = shift;
  my $storable = shift || throw('Must pass a Storable object argument');

  (my $class = ref($storable)) =~ s/.*:://;
  my $table_name = get_table_name_from_class($class);

  if(! exists $tracking_columns{$table_name}){
    throw("$class table_name ($table_name) does not correspond to a valid tracking class:\n\t".
      join("\n\t", keys %tracking_columns));
  }

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::'.$class, $storable);

  return $table_name;
}


#Move this to DBAdaptor?
#Not BaseAdaptor as this is not a BaseAdaptor (is multi table adaptor)
#but probably shouldn't be DBAdaptor either?

sub get_table_name_from_class{
  my $class = shift || throw('Must provide a class argument');
  return lc(join('_', split_CamelCase($class)));
}



1;
