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

Bio::EnsEMBL::Funcgen::Utils::Helper
 
=head1 SYNOPSIS


 e.g. 


 my $object = Bio::EnsEMBL::Object->new
 (
     logging     => 1,
     log_file    => "/tmp/Misc.log",
     debug_level => 2,
     debug_file  => "/tmp/Misc.dbg",
 );

 $object->log("This is a log message.");
 $object->debug(1,"This is a debug message.");
 $object->system("rmdir /tmp/test");


 ----------------------------------------------------------------------------


=head1 OPTIONS

=over 8


=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> performs several debugging and logging functions, aswell as providing several inheritable EFGUtils methods.

=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Utils::Helper;

use strict;
use warnings;
use Carp;    #? Can't use unless we can get it to redirect
use File::Basename;
use Data::Dumper                           qw( Dumper );
use Bio::EnsEMBL::Utils::Exception         qw( throw stack_trace );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref check_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( get_date assert_refs );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::ResultSet;

use base qw( Bio::Root::Root );

#todo remove or use Bio::Root::Root properly
#we get throw from Exception and we have re-implemented debug
#and verbose (in subclasses)


#default dump level
#set this to debug level + 1?
$Data::Dumper::Maxdepth = 2;

my @rollback_tables = qw( data_set   feature_set
                          result_set input_set
                          experiment array
                          array_chip experimental_chip );

#List of valid rollback levels, which are set in _velidate_Set_config
#based on set type string
#To be used in conjunction with -full_delete
#values here used in > ops below, to implement
#hierarchical rollback



my %rollback_modes = 
  (
   feature_set => 1,
   #deletes all associated SetFeatures (etc.) and revokes_states on FeatureSet
   #full delete also removes feature/data/supporting_set records
 
   #data_set => 2, #This is now handled by feature_set mode


   result_set => 3,
   #revokes_states and remove dbfile_registry records
   #full_delete also removes result_set/_input records
        
   #add final levels
   input_set  => 4,
   #This is moot as there isn't really a rollback mode for InputSets
   #apart from revoking the states which can be done directly
   #should have delete_InputSet instead?
   
   #input_subsets may be done independantly? as they may belong to >1 input_set/experiment?
   
   experiment => 5,
);

#Some local filevars to avoid assigning to package typeglobs
my ( $DBGFILE, $LOGFILE );

################################################################################

=head2 new

 Description : Constructor method to create a new object with passed or
               default attributes.

 Arg  [1]    : hash containing optional attributes :-
                 log_file    - name of log file (default = undef -> STDOUT)
                 debug_level - level of detail of debug message [1-3] (default = 0 = off)
                 debug_file  - name of debug file (default = undef -> STDERR)

 ReturnType  : Helper

 Example     : my $Helper = new Bio::EnsEMBL::Helper(
                                                      debug_level => 3,
                                                      debug_file  => "/tmp/efg.debug",
                                                      log_file    => "/tmp/efg.log",
                                                     );

 Exceptions  : throws exception if failed to open debug file
             : throws exception if failed to open log   file

=cut

################################################################################

#To do , change to rearrange

sub new {
  my ( $caller, %args ) = @_;

  my ( $self, %attrdata, $argname );
  my $class = ref($caller) || $caller;

  #Create object from parent class
  $self = $class->SUPER::new(%args);

  #we need to mirror ensembl behaviour here
  #use rearrange and set default afterwards if not defined

  # objects private data and default values
  #Not all of these need to be in main

  %attrdata = 
    (
    _tee         => $main::_tee,
    _debug_level => $main::_debug_level,
    _debug_file  => $main::_debug_file,
    _log_file    => $main::_log_file,   #default should be set in caller
    _no_log      => $main::_no_log,    #suppresses log file generation if log file not defined
    _default_log_dir => $main::_default_log_dir, 
    );

  # set each class attribute using passed value or default value
  foreach my $attrname ( keys %attrdata ) {
    ( $argname = $attrname ) =~ s/^_//;    # remove leading underscore
    $self->{$attrname} = ( exists $args{$argname} ) ? $args{$argname} :
      $attrdata{$attrname};
  }

  $self->{'_tee'} = 1 if $self->{'_no_log'};

  #should we undef log_file here too?
  #This currently only turns off default logging

  $self->{_default_log_dir} ||= $ENV{'HOME'} . '/logs';
  $self->{'_report'} = [];

  # DEBUG OUTPUT & STDERR

  #should default to lowest or highest debug level here!

  if ( defined $self->{_debug_level} && $self->{_debug_level} ) {
    $main::_debug_level = $self->{_debug_level};

    if ( defined $self->{_debug_file} ) {
      $main::_debug_file = $self->{_debug_file};

      open( $DBGFILE, '>>', $self->{_debug_file} ) or
        throw("Failed to open debug file : $!");

#open (DBGFILE, "<STDERR | tee -a ".$self->{_debug_file});#Mirrors STDERR to debug file
    }
    else {
      open( $DBGFILE, '>&STDERR' );
    }

    select $DBGFILE;
    $| = 1;    # make debug file unbuffered

    $self->debug( 1,
                  "Debugging started " .
                    localtime() .
                    " on $0 at level " . $self->{_debug_level} . " ..."
    );
  }

  my $log_file = $self->{_log_file};

  # LOG OUTPUT
  if ( defined $self->{_log_file} ) {

    #This causes print on unopened file as we try and log in the DESTROY
    throw(
'You have specified mutually exclusive parameters log_file and no_log'
    ) if ( $self->{'_no_log'} );
    $main::_log_file = $self->{_log_file};

    #we need to implment tee here
    if ( $self->{'_tee'} ) {
      open( $LOGFILE, ' | tee -a ' . $log_file );
    }
    else {
      open( $LOGFILE, '>>', $log_file ) or
        throw("Failed to open log file : $log_file\nError: $!");
    }
  }
  else {
#Change this to get the name of the control script and append with PID.out
#This is to ensure that we always capture output
#We need to also log params
#We will have to call this from the child class.

    #Only do this if we don't have supress default logs set
    #To avoid loads of loags during testing
    if ( !$self->{'_no_log'} ) {

      my @stack     = stack_trace();
      my $top_level = $stack[$#stack];
      my ( undef, $file ) = @{$top_level};
      $file =~ s/.*\///;

      $self->run_system_cmd( 'mkdir ' . $self->{_default_log_dir} )
        if ( !-e $self->{_default_log_dir} );
      $self->{'_log_file'} =
        $self->{_default_log_dir} . '/' . $file . '.' . $$ . '.log';
      warn "No log file defined, defaulting to:\t" .
        $self->{'_log_file'} . "\n";

      #we should still tee here
      if ( $self->{'_tee'} ) {
        open( $LOGFILE, '| tee -a ' . $self->{'_log_file'} );
      }
      else {
        open( $LOGFILE, '>', $self->{'_log_file'} ) or
          throw( 'Failed to open log file : ' .
                 $self->{'_log_file'} . "\nError: $!" );
      }

    } ## end if ( !$self->{'_no_log'...})
    else {
      #Have to include STD filehandles in operator
      open( $LOGFILE, '>&STDOUT' );
    }
  } ## end else [ if ( defined $self->{_log_file...})]

  select $LOGFILE;
  $| = 1;    # make log file unbuffered
  $self->log( "\n\nLogging started at " . localtime() . "..." );

  # RESET STDOUT TO DEFAULT
  select STDOUT;
  $| = 1;

  $self->debug( 2, "Helper class instance created." );

  return $self;
} ## end sub new

################################################################################

=head2 DESTROY

 Description : Called by gargbage collection to enable tidy up before object deleted

 ReturnType  : none

 Example     : none - should not be called directly

 Exceptions  : none

=cut

################################################################################

sub DESTROY {
  my ($self) = @_;

  #This prevents having to explicitly call report
  $self->report;

  if ( $self->{_log_file} ) {
    $self->log( "Logging complete " . localtime() . "." );
    $self->log( 'Virtual Memory ' . `ps -p $$ -o vsz |tail -1` );
    $self->log( 'Resident Memory ' . `ps -p $$ -o rss |tail -1` );

#       close LOGFILE;  # if inherited object then cannot close filehandle !!!
  }

  if ( $self->{_debug_level} ) {
    $self->debug( 1, "Debugging complete " . localtime() . "." );

#       close DBGFILE;  # if inherited object then cannot close filehandle !!!
  }

  if ( defined $self->{'_timer'} ) {
    $self->{'_timer'}->report();
  }

  $self->debug( 2, "Bio::EnsEMBL::Helper class instance destroyed." );

  return;
} ## end sub DESTROY

##Need generic method in here to get stack and line info
###Use Root.pm stack methods!
# and replace this with caller line method for logging

sub _get_stack {
  my ($self) = shift;

#need to resolve this method with that in debug, pass log or debug arg for different format

  my @prog =
    ( caller(2) ) ? caller(2) :
    ( caller(1) ) ?
    caller(1) :
    ( undef, "undef", 0 );

  return "[" .
    localtime() . " - " . basename( $prog[1] ) . ":$prog[2]]";
}

################################################################################

=head2 log

 Arg[0]      : string  - log message.
 Arg[1]      : boolean - memory usage, appends current process memory stats
 Description : Method to write messages to a previously set up log file.
 Return type : none
 Example     : $root->log("Processing file $filename ...", 1);
 Exceptions  : none

=cut

################################################################################

sub log {
  my ( $self, $message, $mem, $date, $no_return ) = @_;

  if ($mem) {
    $message .= " :: " . `ps -p $$ -o vsz |tail -1`;
    chomp $message;
    $message .= " KB";
  }

  if ($date) {
    my $time = localtime();
    chomp($time);
    $message .= ' - ' . localtime();
  }

  $message .= "\n" if ! $no_return;

  print $LOGFILE "::\t$message";

  if( ! $self->{'_no_log'}){
    $self->debug( 1, $message );
  }
  
} ## end sub log

################################################################################

=head2 report

 Arg[0]      : optional string  - log message.
 Arg[1]      : optional boolean - memory usage, appends current process memory stats
 Description : Wrapper method for log, which also stores message for summary reporting
 Return type : none
 Example     : $root->report("WARNING: You have not done this or that and want it reported at the end of a script");
 Exceptions  : none

=cut

################################################################################

sub report {
  my ( $self, $message, $mem ) = @_;

  if ( defined $message ) {

    $self->log( $message, $mem );

    push @{ $self->{'_report'} }, $message;
  }
  elsif ( scalar( @{ $self->{'_report'} } ) ) {
    my $report = "\n::\tSUMMARY REPORT\t::\n".
      join( "\n", @{ $self->{'_report'} } ) . "\n";
    print $LOGFILE $report;
    
    if($self->{report_fail}){
      die($report); 
    }

    $self->{'_report'} = [];
  }

  return;
}

sub report_fail{
  my ( $self, $msg, @args ) = @_;
  $self->{report_fail} = 1;
  return $self->report("FAIL:\t".$msg, @args); 
}


#add report_warning here
#should print to STDERR or $DBGFILE?


################################################################################

=head2 log_header

 Arg[0]      : string  - log message.
 Arg[1]      : boolean - memory usage, appends current process memory stats
 Description : Wrapper method to format a log as a header line
 Return type : none
 Example     : $root->log("Processing file $filename ...", 1);
 Exceptions  : none

=cut

################################################################################

sub log_header {
  my ( $self, $message, $mem, $date ) = @_;

  print $LOGFILE "\n\n";
  $self->log( "::\t$message\t::\t::", $mem, $date );
  print $LOGFILE "\n";
}

################################################################################

=head2 log_error

 Description : Wrapper method to log, warns to STDERR first
 Return type : none
 Example     : $root->log_error("could not find", ...);
 Exceptions  : none

=cut

################################################################################

sub log_error {
  my ( $self, @args ) = @_;

  warn $args[0] . "\n";

  return $self->log(@args);
}

################################################################################

=head2 debug

 Description : Method to write debug info to a previously set up debug file.
               Over-rides Root.pm on/off style debugging

 Args        : int: debug level and string: log message.

 ReturnType  : none

 Example     : $root->debug(2,"dir=$dir file=$file");

 Exceptions  : none

=cut

###############################################################################

sub debug {
  my ( $self, $level, $message, $data, $depth ) = @_;

  # if debug on at the requested level then output the passed message
  if ( defined $self->{_debug_level} &&
       $level <= $self->{_debug_level} ) {
  
    

    if($data){
      $message.= "\n".$self->dump($data, $depth); 
    }

    my ( @call, $cnt, $prog_name, $prog_line, $call_name, $call_line, $debug_caller );
    #$prog_name = $call_name = $debug_caller = "undef";
    #$prog_line = $call_line = $cnt = 0;
  
    ######Replace this with Carp method?
    
    while ( @call = caller( $cnt++ ) ) {
    #  warn "call($cnt) is @call";
      if ( $cnt ==  1 ) { ##This is the caller of debug
        #$debug_caller = $call[2];  
        $call_line = $call[2];  
      }

      if ( $cnt == 2 ) { #This is the caller of the method which debugs
        #$call_name = basename( $call[1] );
        #$call_line = $call[2];
        $call_name = $call[3];
      }

      #$prog_name = basename( $call[1] );
      #$prog_line = $call[2];
    }

    #This still attempts to print if file not opened
    print $DBGFILE "DEBUG :: $message\t: [ $call_name:$call_line ]\n";
      #"[$$ - $prog_name:$prog_line  $call_name:$call_line]\n";# $debug_caller]\n";

    #carp("carping $message");
  }
} ## end sub debug

################################################################################

#rename this so it doesn't clash with core dump (which causes a 134 SIGABRT)

sub dump {
  my ($self, $data, $depth) = @_;
 
  my $tmp_depth = $Data::Dumper::Maxdepth;
     
  if($depth &&
     ($depth != $Data::Dumper::Maxdepth) ){ 
      $Data::Dumper::Maxdepth = $depth;
  }
  
  my $dump = Dumper($data)."\n";    
  $Data::Dumper::Maxdepth = $tmp_depth;
  
  return $dump;
}



################################################################################

=head2 run_system_cmd

 Description : Method to control the execution of the standard system() command

 ReturnType  : none

 Example     : $Helper->debug(2,"dir=$dir file=$file");

 Exceptions  : throws exception if system command returns none zero

=cut

################################################################################

#Move most of this to EFGUtils.pm
#Maintain wrapper here with throws, only warn in EFGUtils

sub run_system_cmd {
  my ( $self, $command, $no_exit ) = @_;

  my $redirect = '';

  $self->debug( 3, "system($command)" );

  # decide where the command line output should be redirected

  #This should account for redirects
  #This just sends everything to 1 no?

  if ( defined $self->{_debug_level} && $self->{_debug_level} >= 3 ) {

    if ( defined $self->{_debug_file} ) {
      $redirect = " >>" . $self->{_debug_file} . " 2>&1";
    }
    else {
      $redirect = "";
    }
  }
  else {
    #$redirect = " > /dev/null 2>&1";
  }

  # execute the passed system command
  my $status    = system("$command $redirect");
  my $exit_code = $status >> 8;

  if ( $status == -1 ) {
    warn "Failed to execute: $!\n";
  }
  elsif ( $status & 127 ) {
    warn sprintf( "Child died with signal %d, %s coredump\nError:\t$!",
                  ( $status & 127 ),
                  ( $status & 128 ) ? 'with' : 'without' );
  }
  elsif ( $status != 0 ) {
    warn
      sprintf( "Child exited with value %d\nError:\t$!\n", $exit_code )
      ;    #get the true exit code
  }

  if ( $exit_code != 0 ) {

    if ( !$no_exit ) {
      throw(
        "System command failed:\t$command\nExit code:\t$exit_code\n$!");
    }
    else {
      warn("System command returned non-zero exit code:\t$command\n".
           "Exit code:\t$exit_code\n$!");
    }
  }

#reverse boolean logic for perl...can't do this anymore due to tab2mage successful non-zero exit codes :/

  return $exit_code;
} ## end sub run_system_cmd

#add sys_get method ehre to handle system calls which retrieve data?
#i.e.backtick commands `find . -name *fasta`
#or use want or flag with above method?
#should open pipe instead to capture error?

sub get_data {
  my ( $self, $data_type, $data_name ) = @_;

#This method is just to provide standard checking for specific get_data/config methods

  if ( defined $data_name ) {
    throw(
      "Defs data name $data_name for type '$data_type' does not exist\n"
    ) if ( !exists $self->{"${data_type}"}{$data_name} );
  }
  else {
    throw("Defs data type $data_type does not exist\n")
      if ( !exists $self->{"${data_type}"} );
  }

  return
    ( defined $data_name ) ? $self->{"${data_type}"}{$data_name} :
    $self->{"${data_type}"};
}

#sub Timer{
#	my ($self) = shift;

#	$self->{'_timer'} = new Devel::Timer()  if(! defined $self->{'_timer'});

#	return $self->{'_timer'};

#}

sub set_header_hash {
  my ( $self, $header_ref, $fields ) = @_;

  my %hpos;

  for my $x ( 0 .. $#{$header_ref} ) {
    $hpos{ $header_ref->[$x] } = $x;
  }

  if ($fields) {

    foreach my $field (@$fields) {

      if ( !exists $hpos{$field} ) {
        throw("Header does not contain mandatory field:\t${field}");
      }
    }
  }

  return \%hpos;
}


#This should move to Utils
#as it is a simple string manipulation

sub get_schema_and_build {
  my ( $self, $dbname ) = @_;
  my @dbname = split /_/, $dbname;
  return [ $dbname[ ( $#dbname - 1 ) ], $dbname[ ($#dbname) ] ];
}

=head2 get_regbuild_set_states

  Arg [1]    : Bio::EnsEMBL::DBAdaptor
  Example    : my ($dset_states, $rset_states, $fset_states) = $helper->get_regbuild_set_states($db);
  Description: Returns Array refs of appropriate states for sets use din the regulatory build
  Returntype : Array
  Exceptions : Warns if cannot find chromosome CoordSystem
  Caller     : HealthChecker & regulatory build code
  Status     : At risk

=cut

sub get_regbuild_set_states {
  my ( $self, $db ) = @_;

  my $cs_a = $db->get_CoordSystemAdaptor;

  #These states need to be mirrored in RegulatorySets.java

  my $chrom_cs = $cs_a->fetch_by_name('chromosome');
  my ( @dset_states, @rset_states, @fset_states );

  if ( !defined $chrom_cs ) {

    #This species most likely does not have a regbuild
    #really just need to get the 'highest' level here
    warn "Could not find Chromosome CoordSystem. " . $db->dbc->dbname .
      ". most likely does not contain a RegulatoryBuild";
  }
  else {
    my $imp_cs_status =
      'IMPORTED_' . $cs_a->fetch_by_name('chromosome')->version;

    #What about non-chromosome assemblies?
    #top level will not return version...why not?
    @dset_states = ('DISPLAYABLE');
    @rset_states = ( @dset_states, 'DAS_DISPLAYABLE', $imp_cs_status );
    @fset_states = ( @rset_states, 'MART_DISPLAYABLE' );
  }

  return ( \@dset_states, \@rset_states, \@fset_states );
} ## end sub get_regbuild_set_states

=head2 define_ResultSet

  Arg [1]    : Hash - set constructor parameters:
                   
  Example    : my $rset = $self->define_ResultSet(%params);
  Description: 
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : 
  Caller     : Importers and Parsers
  Status     : At risk

=cut

#There is some overlap between this and the new migration code, can it be reused?
#The migration code deals with two different DBs
#our comparisons are between an unstored object and one fetched from the DB.


#do we really need result_set_mode here?
#isn't this just recovery?
#no, this can be none also
#both of these are valid, but we need to manage them better
#as recovery is not being passed on to set_mode
#which in turn is use in the rollback methods

#do we even need result_set == none?
#This support creating a DataSet without a ResultSet
#so 'none' here handles validating it doesn't exists
#or deletion if the appropriate rollback level is set?
#This can only happen if we rollback to result_set
#hence, we would necessarily delete the feature_set and data_set
#which might not be what you want to do?
#This behaviour should probably be handled individually in a rollback script?
  
sub define_ResultSet {
  my $self = shift; #Need to modify the stack here before passing to _validate_Set_config
  my ($rollback_level, $ctype, $ftype, $db, 
      $slices, $recover, $full_delete) =
    $self->_validate_Set_config(@_);#this could also return the following if these are generic

  # feature_class can be infered from the feature_type from the inp sets
  my ($name, $anal, $rset_mode, $sets) = 
    rearrange( ['RESULT_SET_NAME','RESULT_SET_ANALYSIS', 
                'RESULT_SET_MODE', 'SUPPORTING_SETS'], @_ );

  throw("Mandatory parameter note defined:\t-RESULT_SET_NAME") if ! defined $name;
  assert_ref($sets, 'ARRAY'); 
  throw('Must pass at least one supporting set') if scalar(@$sets) == 0;
  $self->debug(1, "Defining ResultSet $name");
  
  foreach my $set(@$sets){
    assert_ref($set, 'Bio::EnsEMBL::Funcgen::InputSubset', 'InputSubset support');
        
    #This isn't fully robust as we are not testing if they are defined
    
    if($ctype->name ne $set->cell_type->name){
      throw('Found mismatch between InputSubset '.$set->name.
        " CellType and CellType specified:\n\t".$set->cell_type->name.
        "\t".$ctype->name);
    }
    
    if( (! $set->is_control) &&
        ($ftype->name ne $set->feature_type->name) ){
      throw('Found mismatch between InputSubset '.$set->name.
        " FeatureType and FeatureType specified:\n\t".$set->feature_type->name.
        "\t".$ftype->name);
    }             
  }

  
   
  #slight hack until we sort out 5mC class and db_file_registry.format
  #also overlaps with FeatureSet -feature_class
  #Move this to ResultSet new?
  my $fclass = ($sets->[0]->feature_type->name eq '5mC') ? 
    'dna_methylation' : 'result';
  
 
  
  #This is now mandatory to avoid erroneous naming when omitted
  #my $name ||= $sets->[0]->experiment->name.' '.$anal->logic_name;
  
  #This will catch mandatory params
  my $rset = Bio::EnsEMBL::Funcgen::ResultSet->new
    (-name          => $name,
     -feature_type  => $ftype,
     -cell_type     => $ctype,
     -support       => $sets,
     -analysis      => $anal,
     -feature_class => $fclass   );
 
  my $rset_adaptor   = $db->get_ResultSetAdaptor; 
  my ($stored_rset)  = @{$rset_adaptor->fetch_all_by_name($name)};
      
  return $self->_validate_rollback_Set($stored_rset, $rset, 'result_set', $rollback_level, 
                                       $rset_adaptor, $slices, $recover, $full_delete, 
                                       $rset_mode);
} ## end sub define_ResultSet





sub define_FeatureSet {
  my $self = shift; #Need to modify the stack here before passing to _validate_Set_config
  my ($rollback_level, $ctype, $ftype, $db, 
      $slices, $recover, $full_delete) =
    $self->_validate_Set_config(@_);
 
  my ($name, $desc, $fclass, $dlabel, $anal) = 
    rearrange(['NAME','DESCRIPTION', 'FEATURE_CLASS', 
               'DISPLAY_LABEL', 'FEATURE_SET_ANALYSIS'], @_ );
                
  #This will catch mandatory params
  my $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
    (
     -name          => $name,
     -feature_type  => $ftype,
     -cell_type     => $ctype,
     -analysis      => $anal,
     -feature_class => $fclass,
     -description   => $desc,
     -display_label => $dlabel,
     #-input_set     => $ssets->[0],
    );
 
  my $fset_adaptor = $db->get_FeatureSetAdaptor; 
  my $stored_fset  = $fset_adaptor->fetch_by_name($name);
  
  return $self->_validate_rollback_Set($stored_fset, $fset, 'feature_set', $rollback_level, 
                                       $fset_adaptor, $slices, $recover, $full_delete);
} ## end sub define_FeatureSet



#set modes
#recover = allows rollback of a ResultSet even if it has a dependant DataSet
#none    = Means we don't want the ResultSet to exist

sub _validate_rollback_Set {
  my ($self, $stored_set, $new_set, $set_type, 
      $rollback_level, $adaptor, $slices, $recover, $full_delete, $set_mode) = @_; 
  #set_type is not validated
  #assumes other vars are pre-validated by _validate_Set_config

  #this will fail if we don't have a stored set already
  #$self->debug(1, "Helper::_validate_rollback_set:".
  #                "\n\tStored set:\t".$stored_set.
  #                "\n\tNew set:\t".$new_set.
  #                "\n\tSet type:\t".$set_type.
  #                "\n\tRollback level:\t".$rollback_level.
  #                "\n\tAdaptor:\t".$adaptor.
  #                "\n\tSlices:\t".$slices.
  #                "\n\tSet mode:\t".$set_mode.
  #                "\n\tIMPORTED:\t".$stored_set->has_status('IMPORTED'));


  (my $rollback_method = ucfirst($set_type)) =~ s/_([a-z])/uc($1)/eg;
  $rollback_method =~ s/_//g;
  $rollback_method = 'rollback_'.$rollback_method;


  #move this to define_ResultSet

  if($set_mode){
    #This rollback shortcut is a bit result set specific, 
    #but could be generic hence why it is here
    
    if($set_type ne 'result_set'){
      throw('set_mode is currently only valid for result_sets'); 
    }
    elsif($set_mode eq 'none'){
      
      if($stored_set){
    
        if($rollback_level < $rollback_modes{result_set}){
          throw('Cannot import with -result_set_mode none with a pre-existing result set.'.
                ' Please specify -rollback result_set');
        }
        else{ #rollback directly here as we aren't concerned with diffs
          $self->$rollback_method($stored_set, undef, undef, 'full'); #full delete
        }
      }
   
      return;  #Bail out early as we don't want a ResultSet!
    }
    elsif($set_mode ne 'recover'){
     throw("Invalid -result_set_mode $set_mode. ".
           'Please omit or specify -result_set_mode none|recover'); 
    }
  }

  if(! defined $stored_set){
    ($stored_set) = @{ $adaptor->store($new_set) };
  }
  else{    
    #Do compare so we catch differences even if we are not rolling back
    #otherwise we may have specified a new set, but this will be silently replaced
    #with the stored set.
    #This fails if there are diffs present and rollback level is insufficient
    my $diffs = $self->_compare_set_for_rollback($new_set, $stored_set, $set_type, 
                                                 $rollback_level, $slices);
    
    #FINALLY DO ROLLBACK
    #Independant of whether there are diffs as there maybe some fault with the 
    #data which has not been caught by the diffs e.g. truncated input
         
    #we have some discrepancy between recovery and set_mode eq 'recover'  
    #Need to change this so we just use recovery instead of set_mode
    #set_mode is then only 'none', so can become a no_result_set flag instead
     
    #is it safe to use recovery for set_mode here?
    #this would only be unsafe if we ever want to treat sets differently wrt recovery
    #that will never happen  
     
    $recover ||=0; #to avoid undef in debug  
    $self->debug(1, "rollback: if ((\$rollback_level($rollback_level) >= \$rollback_modes{$set_type}(".$rollback_modes{$set_type}.")) ||\n".
      "(! \$stored_set->has_status('IMPORTED')(".$stored_set->has_status('IMPORTED').")) ||\n".
      " \$recover($recover))\n". 
      "   &&\n".
      "   ( ! (( \$stored_set->has_status('IMPORTED')(".$stored_set->has_status('IMPORTED').")      &&\n".
      "     (\$rollback_level($rollback_level ) < \$rollback_modes{$set_type}(".$rollback_modes{$set_type}.")) ))");   

    if( (($rollback_level >= $rollback_modes{$set_type}) ||
         (! $stored_set->has_status('IMPORTED')        ) ||
         $recover                                       ) 
         &&
         ( ! (( $stored_set->has_status('IMPORTED')      &&
              ($rollback_level < $rollback_modes{$set_type})) ))){ 
      
      #unset full delete if we do not have the sufficient rollback level
      #i.e. we are rolling back not because of the rollback level, but other criterion
      #full_delete is always specific to rollback level
      $full_delete  = ($rollback_level < $rollback_modes{$set_type}) ? 0 : $full_delete;
  
      if(%$diffs && ! $full_delete){
        throw('Cannot rollback '.ref($new_set).' '.$new_set->name.
          " without specifying full_delete as it has diffs (stored vs new):\n".$self->dump($diffs));  
        #Some times it maybe valid to have diffs
        #but only when rolling back, either:
        #1 With full delete
        #or
        #2 We want to replace some linked data e.g. result_set_inputs
        #The types of diffs allowed are highly dependant on the rollback mode (full/replace)
        #and the set type. So leave this to the individual rollback methods
        #This will be tricky due to traversing dependant sets and their states
      }
  
  
      #We just want to unset it if rollback level is insufficient
      #what about clash between fulldelete feature_set 
      #and recover result_set
      #Do we need recovery to have a level?
      #Just let this fail for now and sort out with rollback script
      
      
      
      #This is not right?
      #do we really want to full delete?
      #or should recover over-ride this?
      #do we need a -full_rollback? param, otherwise
      #-full_rollback and -recovery should be mutually exclusive?
      #Or can we full rollback feature_set, the recover ResultSet?
      
      #This will currently rollback an Unimported ResultSet
      #even if we only specified feature_set rollback level
      
      #Currently not necessary to have an IMPORTED ResultSet before we load a FeatureSet
      
      #So we want to:
      #1 rollback un-IMPORTED sets (with dependants if recovery)
      #2 rollback un/IMPORTED sets if rollback level sufficient (full if -full_delete)
      # slices?
  
  
      #Can't leave full delete to script if we have result_set_mode == none.
      #Do we even want to consider full delete here?
      #let's leave this to a script?
    
      #recover simply allows the above if there are dependants sets
      #although not for full delete
      
    
      
      #These are Importer methods
      #Where is the split between Importer and Helper methods here
      #rollback methods should definitely go in the Helper
      #as we may want to use those outside of the Importer
      #should define_Set methods go in the Importer?
      #probably
      
      #what else uses define_and_validate_sets?
      #but isn't an importer?
      #regbuild!!!!????
      
      
      
      #How do we allow slice rollback here?
      #This needs to be done for IMPORTED sets
      #Should not slice rollback on un IMPORTED set? NO this will screw import slice IMPORT
      #There is a grey area here, where one might perform a slice rollback
      #on an partially IMPORTED set. This may then get erroneously marked as IMPORTED?
      #Can we use batch job as a way to get around this?
      #The assumption being that batch_jobs will have already been rolled back prior to submission?
      #or can we specify a no_slice_rollback flag?
      
     
      #Handle full delete and slices clash as rollback_Data/ResultSet don't take slices arg
      #-recovery has to be set for a slice run
     
      
      #This depends on the rollback level!
      #We pass full delete or slices dependant on the rollback level!
      #Catch this here on in 
      
      #This is currently rollingback everything if recover is set
      
      #This cannot happen if the set is IMPORTED! and the rollback level is insufficient
      
      # DO THE ROLLBACK
      $self->$rollback_method($stored_set, $slices, $recover, $full_delete);
          
      # REDEFINE STORED SET AFTER FULL DELETE
      if($full_delete){
        ($stored_set) = @{ $adaptor->store($new_set) };
      }      
    }
    #else sets are identical and we have no rollback specified
  }
 
  return $stored_set;
}


#Currently sets FeatureSet and ResultSet f/ctypes to same value  

sub _validate_Set_config {
  my $self = shift;
  my ($db, $rollback, $slices, $ftype, $ctype, $fclass, $recover, $full_delete) = 
   rearrange( ['DBADAPTOR', 'ROLLBACK', 'SLICES', 'FEATURE_TYPE',  
               'CELL_TYPE', 'FEATURE_CLASS', 'RECOVER', 'FULL_DELETE'], @_ );
 
  my $rollback_level = 0;

  if ($rollback) {

    if ( ! exists $rollback_modes{$rollback} ) {
      throw("$rollback is not a valid rollback mode, please specify one of the following:\n\t"
          . join( ', ', keys %rollback_modes ) . "\n" );
    }
    
    $rollback_level = $rollback_modes{$rollback};
  }

  if ( $slices && ( ref($slices) ne 'ARRAY' ) ) {
    throw('-slices param must be an ARRAYREF of Bio::EnsEMBL::Slice objects');
    #Rest of slice validation done in other methods
  }

  #Check mandatory params
  if ( ! (ref($db) && $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))) {
    throw('Must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }

  #returns rearranged vars for convenience 
  return ( $rollback_level, $ctype, $ftype, $db, $slices, $recover, $full_delete);  
}


#This is a controlling method, which calls define_Result/FeatureSet
#can we extend this for use in the registration pipeline
#by adding define_InputSet & define_InputSubset etc.

sub define_DataSet {
  my $self = shift; 
  
  my ($fclass, $name, $ssets, $db, $rset_mode) = 
    rearrange(['FEATURE_CLASS' ,'NAME', 'SUPPORTING_SETS',
               'DBADAPTOR', 'RESULT_SET_MODE'], @_ );
  #No need to _validate_Set_config here as it will be done in the below methods
 
  
  #Always create the FeatureSet first
  #This will rollback_DataSet if required  
  my $fset = $self->define_FeatureSet(@_);
  
  #Always create the ResultSet before the DataSet
  #so we can pass it as support
  my $rset;
  
  if($fclass eq 'annotated'){
    #what about result_set_mode = none here?
    #This is handled via define_ResultSet
    #which will also handle a rollback/full delete
    #this is a bit odd that define_ResultSet will delete 
    #alternate is to fetch the stored set and call _validate_rollback_Set
    #directly from here
    
    if(defined $ssets){
      
      if(! ((scalar(@$ssets) == 1) &&
            check_ref($ssets->[0], 'Bio::EnsEMBL::Funcgen::ResultSet'))){ #We likely have InputSubsets
        $ssets = [$self->define_ResultSet(@_)];
      }#else we have a single ResultSet
    }
    else{
      throw("Mandatory parameter not defined:\t-SUPPORTING_SETS"); 
    }
  }   
    
 
  my $dset_adaptor = $db->get_DataSetAdaptor;

  #Generate new DataSet first to validate the parameters
  my $new_dset = Bio::EnsEMBL::Funcgen::DataSet->new
    (-name            => $name,
     -feature_set     => $fset,
     -supporting_sets => $ssets);
  my $stored_dset = $dset_adaptor->fetch_by_name($name);
 
  #Could we actually cascade a compare_to from DataSet to feature_set and result_set?
  #This would be much harder to manage rollbacks
  #We are effectively doing this above, without complicating the the DataSet compare_to method
  #If we did it all in one go, then we would have to define eveything here
  #meaning the define_ResultSet and define_FeatureSet would beome redundant?
  #No, we could remove the _validate_rollback_Set from the return value
  #and let that be called in a controlling method define_validate_rollback_DataSet ? 
  #it would also be very hard to cascade rollback of a DataSet withing _validate_rollback_Set?
  #Hence we would have to have a control method, which is essentially this method!
  #change the name of this method to define_validate_rollback_DataSet???
  #should we change this other define method, such that they don't call _validate_rollback_Set?
  #Would we have a use for the very simple define Set methods?

  if(defined $stored_dset){    
    $self->_compare_set_for_rollback($new_dset, $stored_dset, 'data_set');
    #No need to pass slices here as we can't rollback a data_slice based on a slice
  }
  else{
    ($stored_dset) = @{$dset_adaptor->store($new_dset)};
  }
  
 return $stored_dset;
} ## end sub define_DataSet


=head2 rollback_FeatureSet

  Arg [0]    : Bio::EnsEMBL::Funcgen::FeatureSet
               for another DataSet.
  Arg [1]    : Arrayref of Bio::EnsEMBL::Slice objects to rollback (Optional)
  Arg [2]    : Boolean flag to perform full delete of feature_set and data_set records. 
  Example    : $self->rollback_FeatureSet($fset);
  Description: Deletes all features for the passed FeatureSet, along with associated
               status and xref records. Checks whether FeatureSet is a supporting set 
               in any other DataSet.
  Returntype : None
  Exceptions : Throws if:
                Any deletes fails
                Dependant DataSets are found
                FeatureSet 'adaptor' method is not defined
                Defined Slices are not valid
                More than 1 sub-Slice is defined (i.e. not full length)
  Caller     : Importers and Parsers
  Status     : At risk

=cut

#Not safe to implement recover here, as we never want to rollback
#and supporting FeatureSet without re-running the depedant DataSet
#allow force tho, so we can re-call a set without having to rollback 
#the depedant dataset (i.e. regbuild) before hand?
#Can we we add a special status to the DataSet? SUPPORT_CHANGED?
#There should never be any SUPPORT_CHANGED states present before release?

#or should this be the recovery mode analagous to the ResultSet recovery mode?
#force is for full delete, should be force_delete?

#todo change recover to also take the 'known' data set
#and allow implicity recovery, or full delete.


sub rollback_FeatureSet {
  my ( $self, $fset,  $slices, $recovery, $full_delete, $force_delete ) = @_;
  my ( $sql, $slice_name );
  my $slice_join = '';
   
  #print "Full: $full_delete\n";
  if ( ! ($fset &&
         (ref($fset) eq 'Bio::EnsEMBL::Funcgen::FeatureSet') &&
          defined $fset->adaptor ) ){
    throw('Must provide a valid stored Bio::EnsEMBL::Funcgen::ResultSet');
  }
    
  my $db    = $fset->adaptor->db; #Assumes access to DBAdaptor
  $db->is_stored_and_valid( 'Bio::EnsEMBL::Funcgen::FeatureSet', $fset );
  my $table = $fset->feature_class . '_feature';
  
  $self->log_header('Rolling back '.$fset->feature_class.
                    " FeatureSet:\t" . $fset->name );

  #Check whether this is a supporting set for another data_set
  my @dsets =  @{ $db->get_DataSetAdaptor->fetch_all_by_supporting_set($fset) };

  if (@dsets) {
    my $txt = $fset->name." is a supporting set of the following DataSets:\t" .
              join( ', ', ( map { $_->name } @dsets ) );
    throw( $txt ."\nPlease resolve by deleting dependant Feature/DataSets" ); 
    #This will currently not allow rollback if we have already associated it 
    #with a regulatory build
  }
  
  
  #Validate all slices before we commence any rollback
  if ($slices) {

    if ($full_delete) {
      throw("Cannot do a \'full delete\' on a sub set of \'slices\' for FeatureSet:\t".
        $fset->name."\nPlease omit one");  
    }

    assert_refs($slices, 'Bio::EnsEMBL::Slice');
   
    $self->log( "Restricting to slices:\n\t\t" .
                join( "\n\t\t", ( map { $_->name } @$slices ) ) );
  }
  else{ #Set undef slice for no slice definition
    $slices ||= [undef];  
  }

  
  for my $slice(@$slices){
    my $slice_join = '';
    
    if(defined $slice){  
      my $efg_sr_id =
      $fset->get_FeatureAdaptor->get_seq_region_id_by_Slice($slice);

      if ( ! $efg_sr_id ) {
        $self->log("Slice is not present in eFG DB:\t" . $slice->name );
      }
      else {   
        #Test is not subslice
        my $full_slice = $slice->adaptor->fetch_by_region(undef,
                                                          $slice->seq_region_name);

        if ( ( $slice->start != 1 ) ||
             ( $full_slice->end != $slice->end ) ){
          $slice_join = " and f.seq_region_id=$efg_sr_id and f.seq_region_start<=" .
                        $slice->end.' and f.seq_region_end>=' . $slice->start;
        }
        else{
          $slice_join = ' and f.seq_region_id = '.$efg_sr_id; 
        }
      }
    }
           
    # Now do the rollback 
    $fset->adaptor->revoke_states($fset);
        
    if ( $fset->feature_class eq 'regulatory' ) {   #Rollback reg attributes
      $sql = "DELETE ra from regulatory_attribute ra, $table f where ".
        "f.${table}_id=ra.${table}_id and f.feature_set_id=".$fset->dbID.$slice_join;
          
      $db->rollback_table( $sql, 'regulatory_attribute');
    }
    elsif($fset->feature_class eq 'annotated'){   #Handle amfs
      $sql = 'DELETE amf from annotated_feature af, associated_motif_feature amf where '.
        'af.feature_set_id='.$fset->dbID.
        ' AND af.annotated_feature_id = amf.annotated_feature_id'.$slice_join;
      $db->rollback_table( $sql, 'associated_motif_feature' );
    } ## end if ( $fset->feature_class...)
      
    #Remove object_xref records (but not xref which may be used by soemthing else)
    $sql = "DELETE ox from object_xref ox, $table f where ox.ensembl_object_type='"
      .ucfirst( $fset->feature_class )."Feature' and ox.ensembl_id=f.${table}_id and ".
      "f.feature_set_id=". $fset->dbID . $slice_join;
    $db->rollback_table( $sql, 'object_xref', 'object_xref_id' );
      
    #Remove associated_feature_type records
    $sql = "DELETE aft from associated_feature_type aft, $table f where ".
      "f.feature_set_id=".$fset->dbID." and f.${table}_id=aft.table_id and ".
      "aft.table_name='".$fset->feature_class . "_feature'" . $slice_join;
    $db->rollback_table( $sql, 'associated_feature_type' );
      
    #Remove features
    $sql = "DELETE f from $table f where f.feature_set_id=" .
      $fset->dbID . $slice_join;
      
    # warn $sql;  
    $db->rollback_table( $sql, $table, "${table}_id" );
  }
  
  if ($full_delete) {  #Also delete feature/data/supporting_set records
    $self->log( "Deleting Feature/DataSet:\t" . $fset->name );
  
    #Delete regbuild strings first
    if ( $fset->feature_class eq 'regulatory' ) {  
      $sql = "DELETE from regbuild_string where name like 'regbuild." . $fset->cell_type->name.".%'";
      $db->rollback_table( $sql, 'regbuild_string', 'regbuild_string_id' );
      $self->log( "Deleted regbuild_string entries for:\t" . $fset->name );
    }  
  
    $sql = "DELETE from feature_set where feature_set_id=" . $fset->dbID;
    $db->rollback_table( $sql, 'feature_set', 'feature_set_id' );
    $self->log( "Deleted feature_set entry for:\t" . $fset->name );

    $sql = 'DELETE ss, ds from data_set ds, supporting_set ss where '.
      'ds.feature_set_id='.$fset->dbID.' AND ds.data_set_id=ss.data_set_id';
    $db->rollback_table( $sql, 'data_set', 'data_set_id' );
    $self->log("Deleted associated data/supporting_set entries for:\t" . $fset->name );
  }  
  
  return;
} ## end sub rollback_FeatureSet


=head2 rollback_ResultSet

  Arg[1]     : Bio::EnsEMBL::Funcgen::ResultSet
  Arg[2]     : Boolean - optional flag to roll back array results
  Example    : $self->rollback_ResultSet($rset);
  Description: Deletes all status. chip_channel and result_set entries for this ResultSet.
               Will also rollback_results sets if rollback_results specified.  This will also
               update or delete associated ResultSets where appropriate.
  Returntype : Arrayref containing the ResultSet and associated DataSet which have not been rolled back
  Exceptions : Throws if ResultSet not valid
               Throws is result_rollback flag specified but associated product FeatureSet found.
  Caller     : General
  Status     : At risk

=cut

#Need to change slice to slices ref here
#Need to add full rollback, which will specify to remove all sets
#as well as results and
#These params need clarifying as their nature changes between input_set and array rsets
#Don't we always want to rollback_results?
#force should only really be used to rollback InputSet ResultFeature sets
#i.e. Read collections which are not used as direct input for the linked product FeatureSet
#This should fail with array data associated with a product feature set


#todo update callers to add full_delete, recovery and force flags as appropriate
#force should really only be used via the script i.e. we don't really want to use force
#in the automated pipeline, as this will probably lead to errors


#todo add slices support here

#It would be better to take the known dataset here rather than 'recovery'

#There is a case here, for just being able to redefine the result_set_input entries
#This can be identified from _compare_set_for_rollback
#Can recover handle this, 
#maybe full_delete and force_delete should just be full_delete 1|2
#and recover/update should be recover 1|2
#force delete is never passed from _validate_rollback_Set at present

sub rollback_ResultSet {
  my ( $self, $rset, $slices, $recover, $full_delete, $force_delete ) = @_;

  $self->debug(1, "rollback_ResultSet:\t ResultSet = ".$rset->name.
    "\trecover = $recover\tfull_delete = $full_delete\tforce_delete = $force_delete");

  if($slices){
    warn "Slice rollback not currently supported, performing full rollback for ResultSet:\t".
      $rset->name; 
  }


  if($full_delete){
    
    if($recover){
      throw('The \'full delete\' and \'recover\' argument are mutually exclusive'.
        ', please omit one or both');   
    }
    
    if ($slices) {
      throw("Cannot do a \'full delete\' on a sub set of \'slices\' for ResultSet:\t".
        $rset->name."\nPlease omit one");  
    }          
  }

  if ( ! ($rset &&
         (ref($rset) eq 'Bio::EnsEMBL::Funcgen::ResultSet') &&
          defined $rset->adaptor ) ){
    throw('Must provide a valid stored Bio::EnsEMBL::Funcgen::ResultSet');
  }
  
  
  my $db = $rset->adaptor->db; #Assumes db is accessible  
  $db->is_stored_and_valid( 'Bio::EnsEMBL::Funcgen::ResultSet', $rset );

  #Just omit experimental_chip/channel/result_feature support for now as we probably 
  #won't ever use it again
  if($rset->table_name ne 'input_subset'){
    throw('rollback_ResultSet not longer support non-InputSet rollbacks');
    #This would need to check co-dependant ResultSets, see old version in cvs
  } 
  
  ### Check if this ResultSet is part of a DataSet with a product feature set
  my @dsets = @{ $db->get_DataSetAdaptor->fetch_all_by_supporting_set($rset) };

  if(@dsets) {
     
    #This will cause failures for slice jobs unless we specify recover
    #won't this also cause slice jobs to rollback in parallel and potentially
    #delete data from a parallel job
    
    #We need a flag to signify that the rollback has alreayd been done
    #how was this done previously?
    #probably as the recovery flag simply used the existing set without rollback
    #was this the distinction between clobber and recover? 
    #Or do the slice modes simply rollback slice only data!
    #There is no way to know whether an import will be done as one or in parallel
    #so we have to rollback in both instances!
    #This is fine
     
    if( ! ($recover || $full_delete)){
      
      if(scalar (@dsets) > 1){
        throw('ResultSet '.$rset->name." has associated DataSets, please specify".
          " -recover or use rollback script with -force_delete:\n\t".join(',', (map $_->name, @dsets)));
      }
      else{
        throw('ResultSet '.$rset->name." has an associated DataSet, please specify".
          " 'recover', 'full_delete' or use rollback script:\n\t".
          join(',', (map $_->name, @dsets)));
      }
    }
    elsif($recover){
      
      if(scalar (@dsets) > 1){
         warn("Recovering/redefining ResultSet which has more than one dependant DataSet:\t".
          $rset->name."\n:DataSets:\n\t".join("\n\t", (map { $_->name } @dsets))."\n");
        #Only warn here else we can never import more than 1 DataSet for a ResultSet
        
        #Can we we add a special status to the DataSet? SUPPORT_CHANGED?
        #There should never be any SUPPORT_CHANGED states present before release?
          
      }
      #else assume we know about this DataSet?     
    }
    elsif(! $force_delete){ # && full delete
      #full_delete may invalidate states for dependant Data/FeatureSets
      #this may cause problems in a pipeline, unless the dependant Data/FeatureSets
      #are handled first.
      #Cannot traverse set relationships as these maybe very complex
      #This is the same for result_set_input redefinition
      #This would require rolling back in reverse order to which we are currently
      #creating/comparing the setsq
      
      #Force delete is completely unsafe, as this will not only result in 
      #an absent supporting result_set, which will break the DataSet
      #but the rollback, may mean that this result_set_id get's re-used
      #by a completely different ResultSet!
    
      throw("Cannot 'full delete' ResultSet which has dependant DataSet(s):\t".
          $rset->name."\n:DataSets:\n\t".join("\n\t", (map { $_->name } @dsets))."\n".
          "Please use rollback script with -force_delete instead");
    }
  }


  #This currently also revokes pre-imported states
  #i.e. ALIGNED/ING
  $db->get_ResultSetAdaptor->revoke_imported_states_by_Set($rset);  
  
  
  if($rset->table_name ne 'input_subset'){
    throw('rollback_ResultSet does not support non-input_subset/dbfile rollbacks'); 
  }
  
  my $sql;
  
  if($rset->dbfile_data_dir){ #delete the dbfile_registry entry 
    $sql = 'DELETE from dbfile_registry where table_name="result_set" and table_id='.$rset->dbID;
     $self->debug(2, "rollback_ResultSet - deleting dbfile_registry enrty:\n\t$sql");
    $db->rollback_table($sql, 'dbfile_registry');
  }

  #We may want to just change the result_set_inputs, without changing the result_set record?
  #is this sane? as it may involve leaving an orphaned result_set

  if($full_delete){ #delete the result_set and result_set_input entries
  
    if(@dsets){
      $sql = 'DELETE ss from result_set rs, supporting_set ss WHERE '.
        'rs.result_set_id=ss.supporting_set_id AND ss.type="result"';
      $self->debug(2, "rollback_ResultSet - deleting supporting_set entries:\n\t$sql");
      $db->rollback_table($sql, 'supporting_set' ); 
    } 
  
    $db->get_ResultSetAdaptor->revoke_states($rset);
  
    $sql = 'DELETE rs, rsi from result_set rs, result_set_input rsi WHERE '.
      'rsi.result_set_id=rs.result_set_id AND rs.result_set_id='.$rset->dbID;
      
    $self->debug(2, "rollback_ResultSet - deleting ResultSet:\n\t$sql" );
    $db->rollback_table($sql, 'result_set', 'result_set_id');   
  }
  
 $self->debug(1, "rollback_Result - finished rolling back:\t" . $rset->name ); 
  
 return;
} ## end sub rollback_ResultSet


#todo add delete_DataSet (with flag to full_delete supporting_sets?
#todo add delete_InputSet


=head2 rollback_ArrayChips

  Arg[1]     : ARRAYREF: Bio::EnsEMBL::Funcgen::ArrayChip objects
  Example    : $self->rollback_ArrayChips([$achip1, $achip2]);
  Description: Deletes all Probes, ProbeSets, ProbeFeatures and 
               states associated with this ArrayChip
  Returntype : None
  Exceptions : Throws if ArrayChip not valid and stored
               Throws if ArrayChips are not of same class
  Caller     : General
  Status     : At risk

=cut

#This should be tied to a CS id!!!
#And analysis dependant?
#We may not want to delete alignment by different analyses?
#In practise the slice methods ignore analysis_id for this table
#So we currently never use this!
#So IMPORTED status should be tied to CS id and Analysis id?

sub rollback_ArrayChips {
  my ( $self, $acs, $mode, $force, $keep_xrefs, $no_clean_up,
       $force_clean_up ) = @_;

#no_clean_up and force_clean_up allow analyze/optimize to be skipped until the last rollback
#We could get around this by specifying all ArrayChips for all formats at the same time?
#Need to implement in RollbackArrays

  $mode ||= 'probe';

  if ( $mode &&
       ( $mode ne 'probe' &&
         $mode ne 'probe_feature'        &&
         $mode ne 'ProbeAlign'           &&
         $mode ne 'ProbeTranscriptAlign' &&
         $mode ne 'probe2transcript' ) )
  {
    throw(
"You have passed an invalid mode argument($mode), you must omit or specify either 'probe2transcript', 'probe', 'ProbeAlign, 'ProbeTranscriptAlign' or 'probe_feature' for all of the Align output"
    );
  }

  if ( $force && ( $force ne 'force' ) ) {
    throw(
"You have not specified a valid force argument($force), you must specify 'force' or omit"
    );
  }

  if ( $keep_xrefs && ( $keep_xrefs ne 'keep_xrefs' ) ) {
    throw(
"You have not specified a valid keep_xrefs argument($keep_xrefs), you must specify 'keep_xrefs' or omit"
    );
  }

  if ($keep_xrefs) {

    if ( $mode eq 'probe' || $mode eq 'probe2transcript' ) {
      throw(
"You cannot specify 'keep_xrefs' with mode $mode, you can only rollback features e.g. probe_feature, ProbeAlign or ProbeTranscriptAlign"
      );
    }

    if ($force) {
      throw(
"You cannot 'force' delete the probe2transcript xrefs and 'keep_xrefs' at the same time. Please specify just one."
      );
    }
  }

  my ( $adaptor, $db, %classes );

  foreach my $ac (@$acs) {
    $adaptor ||=
      $ac->adaptor || throw('ArrayChip must have an adaptor');
    $db ||= $adaptor->db;
    $db->is_stored_and_valid( 'Bio::EnsEMBL::Funcgen::ArrayChip', $ac );

    if ( !$ac->get_Array->class ) {
      throw(
'The ArrayChip you are trying to rollback does not have a class attribute'
      );
    }

    $classes{ $ac->get_Array->class } = undef;

#if($class && ($class ne $ac->get_Array->class)){
#  throw('You can only rollback_ArrayChips for ArrayChips with the same class');
#}
  }

#This is always the case as we register the association before we set the Import status
#Hence the 2nd stage of the import fails as we have an associated ExperimentalChip
#We need to make sure the ExperimentalChip and Channel have not been imported!!!
  warn
"NOTE: rollback_ArrayChips. Need to implement ExperimentlChip check, is the problem that ExperimentalChips are registered before ArrayChips imported?";

#Check for dependent ExperimentalChips
#if(my @echips = @{$db->get_ExperimentalChipAdaptor->fetch_all_by_ArrayChip($ac)}){
#	my %exps;
#	my $txt = "Experiment\t\t\t\tExperimentalChip Unique IDs\n";

  #	foreach my $ec(@echips){
  #	  $exps{$ec->get_Experiment->name} ||= '';

  #	  $exps{$ec->get_Experiment->name} .= "\t".$ec->unique_id;
  #	}

  #	map {$txt.= "\t".$_.":".$exps{$_}."\n"} keys %exps;

  #	throw("Cannot rollback ArrayChip:\t".$ac->name.
  #		  "\nFound Dependent Experimental Data:\n".$txt);
  #  }

  my $ac_names = join( ', ', ( map { $_->name } @$acs ) );
  my $ac_ids   = join( ', ', ( map { $_->dbID } @$acs ) );

  $self->log("Rolling back ArrayChips $mode entries:\t$ac_names");
  my ( $row_cnt, $probe_join, $sql );

#$ac->adaptor->revoke_states($ac);#This need to be more specific to the type of rollback
  my $species = $db->species;

  if ( !$species ) {
    throw(
'Cannot rollback probe2transcript level xrefs without specifying a species for the DBAdaptor'
    );
  }

  #Will from registry? this return Homo sapiens?
  #Or homo_sapiens
  ( $species = lc($species) ) =~ s/ /_/;

  my $transc_edb_name = "${species}_core_Transcript";
  my $genome_edb_name = "${species}_core_Genome";

#Maybe we want to rollback ProbeAlign and ProbeTranscriptAlign output separately so we
#can re-run just one part of the alignment step.

#We want this Probe(Transcript)Align rollback available in the environment
#So we can do it natively and before we get to the RunnableDB stage,
#where we would be trying multiple rollbacks in parallel
#Wrapper script?
#Or do we keep it simple here and maintain probe_feature wide rollback
#And just the ProbeAlign/ProbeTranscriptAlign roll back in the environment?

  #We can restrict the probe deletes using the ac_id
  #We should test for other ac_ids using the same probe_id
  #Then fail unless we have specified force delete

  #These should be deleted for all other modes but only if force is set?
  #This may delete xrefs for other ArrayChips

#The issues is if we need to specify force for one delete but don't want to delete something else?
#force should only be used to delete upto and including the mode specified
#no mode equates to probe mode
#if no force then we fail if previous levels/modes have xrefs etc...

#Let's grab the edb ids first and use them directly, this will avoid table locks on edb
#and should also speed query up?

  if ( $mode eq 'probe2transcript' || $force ) {

    #Delete ProbeFeature UnmappedObjects
    $self->log("Deleting probe2transcript ProbeFeature UnmappedObjects");
    $sql = "
      DELETE 
        uo 
      FROM 
        analysis a, 
        unmapped_object uo, 
        probe p, 
        probe_feature pf, 
        external_db e 
      WHERE 
        a.logic_name           = 'probe2transcript'   AND 
        a.analysis_id          =  uo.analysis_id      AND
        p.probe_id             =  pf.probe_id         AND 
        pf.probe_feature_id    =  uo.ensembl_id       AND 
        uo.ensembl_object_type = 'ProbeFeature'       AND 
        uo.external_db_id      =  e.external_db_id    AND
        e.db_name              = '${transc_edb_name}' AND 
        p.array_chip_id       IN ($ac_ids)
        ";
    $db->rollback_table( $sql, 'unmapped_object', 'unmapped_object_id', $no_clean_up );

    #Delete ProbeFeature Xrefs/DBEntries
    $self->log("Deleting probe2transcript ProbeFeature Xrefs");
    $sql = "
      DELETE 
        ox 
      FROM 
        xref          x, 
        object_xref   ox, 
        probe         p, 
        probe_feature pf, 
        external_db   e 
      WHERE 
        x.external_db_id        = e.external_db_id        AND 
        e.db_name               = '${transc_edb_name}'    AND 
        x.xref_id               = ox.xref_id              AND 
        ox.ensembl_object_type  ='ProbeFeature'           AND 
        ox.ensembl_id           = pf.probe_feature_id     AND 
        pf.probe_id             = p.probe_id              AND 
        ox.linkage_annotation  != 'ProbeTranscriptAlign'  AND 
        p.array_chip_id        IN($ac_ids)
      ";
    $db->rollback_table( $sql, 'object_xref', 'object_xref_id', $no_clean_up );

    #Probe/Set specific entries
    for my $xref_object ( 'Probe', 'ProbeSet' ) {
      $probe_join = ( $xref_object eq 'ProbeSet' ) ? 'p.probe_set_id' :
        'p.probe_id';

      #Delete Probe/Set UnmappedObjects

      $self->log("Deleting probe2transcript $xref_object UnmappedObjects");

      $sql = "DELETE uo FROM analysis a, unmapped_object uo, probe p, external_db e ".
              "WHERE a.logic_name='probe2transcript' AND a.analysis_id=uo.analysis_id AND ".
              "uo.ensembl_object_type='${xref_object}' AND $probe_join=uo.ensembl_id AND ".
              "uo.external_db_id=e.external_db_id AND e.db_name='${transc_edb_name}' ".
              "AND p.array_chip_id IN($ac_ids)";

      #.' and edb.db_release="'.$schema_build.'"';
      $db->rollback_table( $sql, 'unmapped_object', 'unmapped_object_id', $no_clean_up );

      #Delete Probe/Set Xrefs/DBEntries
      $sql = "
        DELETE 
          ox 
        FROM 
          xref        x, 
          object_xref ox, 
          external_db e, 
          probe       p 
        WHERE 
          x.xref_id               = ox.xref_id            AND 
          e.external_db_id        = x.external_db_id      AND
          e.db_name               = '${transc_edb_name}'  AND 
          ox.ensembl_object_type  = '${xref_object}'      AND 
          ox.ensembl_id           = ${probe_join}         AND 
          p.array_chip_id        IN ($ac_ids)
          ";
      $self->log("Deleting probe2transcript $xref_object xref records");
      $db->rollback_table( $sql, 'object_xref', 'object_xref_id', $no_clean_up );
    }
  } ## end if ( $mode eq 'probe2transcript'...)
  elsif ( !$keep_xrefs )
  {    #Need to check for existing xrefs if not force
        #we don't know whether this is on probe or probeset level
     #This is a little hacky as there's not way we can guarantee this xref will be from probe2transcript
     #until we get the analysis_id moved from identity_xref to xref
     #We are also using the Probe/Set Xrefs as a proxy for all other Xrefs and UnmappedObjects
     #Do we need to set a status here? Would have problem rolling back the states of associated ArrayChips

    for my $xref_object ( 'Probe', 'ProbeSet' ) {

      $probe_join = ( $xref_object eq 'ProbeSet' ) ? 'p.probe_set_id' :
        'p.probe_id';

      $row_cnt = $db->dbc->db_handle->selectrow_array("
        SELECT 
          COUNT(*) 
        FROM 
          xref        x, 
          object_xref ox, 
          external_db e, 
          probe       p 
        WHERE 
          x.xref_id               = ox.xref_id            AND 
          e.external_db_id        = x.external_db_id      AND 
          e.db_name               = '${transc_edb_name}'  AND 
          ox.ensembl_object_type  = '${xref_object}'      AND 
          ox.ensembl_id           = ${probe_join}         AND 
          p.array_chip_id        IN ($ac_ids)
        ");

      if ($row_cnt) {
        throw("
          Cannot rollback ArrayChips($ac_names), found $row_cnt $xref_object Xrefs. Pass 'force' argument or 'probe2transcript' mode to delete"
        );
      }
      else {
        #$self->log("Found $row_cnt $xref_object Xrefs");
      }
    }
  } ## end elsif ( !$keep_xrefs )  [ if ( $mode eq 'probe2transcript'...)]

  #ProbeFeatures inc ProbeTranscriptAlign xrefs

  if ( $mode ne 'probe2transcript' ) {

    if ( ( $mode eq 'probe' && $force ) ||
         $mode eq 'probe_feature' ||
         $mode eq 'ProbeAlign'    ||
         $mode eq 'ProbeTranscriptAlign' )
    {

      #Should really revoke some state here but we only have IMPORTED

      #ProbeTranscriptAlign Xref/DBEntries

#my (@anal_ids) = @{$db->get_AnalysisAdaptor->generic_fetch("a.module='ProbeAlign'")};
#Grrrr! AnalysisAdaptor is not a standard BaseAdaptor implementation
#my @anal_ids = @{$db->dbc->db_handle->selectall_arrayref('select analysis_id from analysis where module like "%ProbeAlign"')};
#@anal_ids = map {$_= "@$_"} @anal_ids;

      if ( $mode ne 'ProbeAlign' ) {
        my $lnames = join( ', ',
                           ( map { "'${_}_ProbeTranscriptAlign'" }
                               keys(%classes) ) );

        $sql = "
        DELETE 
          ox 
        FROM 
          object_xref   ox, 
          xref          x,
          probe         p, 
          probe_feature pf, 
          external_db   e 
        WHERE 
          ox.ensembl_object_type  = 'ProbeFeature'         AND 
          ox.linkage_annotation   = 'ProbeTranscriptAlign' AND 
          ox.xref_id              = x.xref_id              AND 
          e.external_db_id        = x.external_db_id       AND 
          e.db_name               = '${transc_edb_name}'   AND 
          ox.ensembl_id           = pf.probe_feature_id    AND 
          pf.probe_id             = p.probe_id             AND 
          p.array_chip_id        IN($ac_ids)
          ";

        $self->log(
            "Deleting ProbeFeature Xref/DBEntry records for:\t$lnames");
        $db->rollback_table( $sql, 'object_xref', 'object_xref_id', $no_clean_up );

#Can't include uo.type='ProbeTranscriptAlign' in these deletes yet as uo.type is enum'd to xref or probe2transcript
#will have to join to analysis and do a like "%ProbeTranscriptAlign" on the the logic name?
#or/and ur.summary_description='Promiscuous probe'?

        $sql = "
        DELETE
          uo 
        FROM 
          unmapped_object uo, 
          probe p, 
          external_db e, 
          analysis a 
        WHERE 
          uo.ensembl_object_type  = 'Probe'               AND 
          uo.analysis_id          = a.analysis_id         AND
          a.logic_name           IN (${lnames})           AND 
          e.external_db_id        = uo.external_db_id     AND
          e.db_name               = '${transc_edb_name}'  AND 
          uo.ensembl_id           = p.probe_id            AND 
          p.array_chip_id        IN ($ac_ids)
                ";

        $self->log("Deleting UnmappedObjects for:\t${lnames}");
        $db->rollback_table( $sql, 'unmapped_object', 'unmapped_object_id', $no_clean_up );

        #Now the actual ProbeFeatures
        $sql = "
        DELETE
          pf 
        FROM 
          probe_feature pf, 
          probe         p, 
          analysis      a 
        WHERE 
          a.logic_name    IN (${lnames})    AND 
          a.analysis_id    = pf.analysis_id AND 
          pf.probe_id      = p.probe_id     AND 
          p.array_chip_id IN ($ac_ids)
          ";
        $self->log("Deleting ProbeFeatures for:\t${lnames}");
        $db->rollback_table( $sql, 'probe_feature', 'probe_feature_id', $no_clean_up );
      } ## end if ( $mode ne 'ProbeAlign')

      if ( $mode ne 'ProbeTranscriptAlign' ) {

        my $lnames = join( ', ',
                           ( map { "'${_}_ProbeAlign'" } keys(%classes)
                           ) );

        $sql = "
        DELETE 
          uo 
        FROM 
          unmapped_object uo, 
          probe           p, 
          external_db     e, 
          analysis        a 
        WHERE 
          uo.ensembl_object_type = 'Probe'              AND 
          uo.analysis_id         = a.analysis_id        AND
          a.logic_name          IN (${lnames})          AND 
          e.external_db_id       = uo.external_db_id    AND
          e.db_name              = '${genome_edb_name}' AND 
          uo.ensembl_id          = p.probe_id           AND 
          p.array_chip_id       IN ($ac_ids)
          ";
        $self->log("Deleting UnmappedObjects for:\t${lnames}");
        $db->rollback_table( $sql, 'unmapped_object', 'unmapped_object_id', $no_clean_up );

        $sql = "
        DELETE 
          pf 
        FROM 
          probe_feature pf, 
          probe p, 
          analysis a 
        WHERE
          a.logic_name    IN(${lnames})       AND 
          a.analysis_id    = pf.analysis_id   AND 
          pf.probe_id      = p.probe_id       AND 
          p.array_chip_id IN ($ac_ids)
          ";
          $self->log("Deleting ProbeFeatures for:\t${lnames}");
        $db->rollback_table( $sql, 'probe_feature', 'probe_feature_id', $no_clean_up );
      }
    } ## end if ( ( $mode eq 'probe'...))
    else {
#Need to count to see if we can carry on with a unforced probe rollback?
#Do we need this level of control here
#Can't we assume that if you want probe you also want probe_feature?
#Leave for safety, at least until we get the dependant ExperimetnalChip test sorted
#What about if we only want to delete one array from an associated set?
#This would delete all the features from the rest?

      $sql =
"select count(*) from object_xref ox, xref x, probe p, external_db e WHERE ox.ensembl_object_type='ProbeFeature' AND ox.linkage_annotation='ProbeTranscriptAlign' AND ox.xref_id=x.xref_id AND e.external_db_id=x.external_db_id and e.db_name='${transc_edb_name}' AND ox.ensembl_id=p.probe_id AND p.array_chip_id IN($ac_ids)";
      $row_cnt = $db->dbc->db_handle->selectrow_array($sql);

      if ($row_cnt) {
        throw(
"Cannot rollback ArrayChips($ac_names), found $row_cnt ProbeFeatures. Pass 'force' argument or 'probe_feature' mode to delete"
        );
      }
      else {
        $self->log("Found $row_cnt ProbeFeatures");
      }
    }

    if ( $mode eq 'probe' ) {

#Don't need to rollback on a CS as we have no dependant EChips?
#Is this true?  Should we enforce a 3rd CoordSystem argument, 'all' string we delete all?

      foreach my $ac (@$acs) {
        $ac->adaptor->revoke_states($ac)
          ;    #Do we need to change this to revoke specific states?
         #Current states are only IMPORTED, so not just yet, but we could change this for safety?
      }

      #ProbeSets
      $sql =
"DELETE ps from probe p, probe_set ps where p.array_chip_id IN($ac_ids) and p.probe_set_id=ps.probe_set_id";
      $db->rollback_table( $sql, 'probe_set', 'probe_set_id', $no_clean_up );

      #Probes
      $sql = "DELETE from probe where array_chip_id IN($ac_ids)";
      $db->rollback_table( $sql, 'probe', 'probe_id', $no_clean_up );
    }
  } ## end if ( $mode ne 'probe2transcript')

  $self->log("Finished $mode roll back for ArrayChip:\t$ac_names");
  return;
} ## end sub rollback_ArrayChips

#This will just fail silently if the reset value
#Is less than the true autoinc value
#i.e. if there are parallel inserts going on
#So we can never assume that the $new_auto_inc will be used

#Todo Move these to DBAdaptor?
#although this will remove logging capability?


=head2 get_core_display_name_by_stable_id

  Args [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Args [2]   : stable ID from core DB.
  Args [3]   : stable feature type e.g. gene, transcript, translation
  Example    : $self->validate_and_store_feature_types;
  Description: Builds a cache of stable ID to display names.
  Returntype : string - display name
  Exceptions : Throws is type is not valid.
  Caller     : General
  Status     : At risk

=cut

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> display_name
# Return hashref keyed on {$type}{$stable_id}
#Need to update cache if we're doing more than one 'type' at a time
# as it will never get loaded for the new type!

sub get_core_display_name_by_stable_id {
  my ( $self, $cdb, $stable_id, $type ) = @_;

  $type = lc($type);

  if ( $type !~ /(gene|transcript|translation)/ ) {
    throw(
      "Cannot get display_name for stable_id $stable_id with type $type"
    );
  }

  if ( !exists $self->{'display_name_cache'}->{$stable_id} ) {
    ( $self->{'display_name_cache'}->{$stable_id} ) =
      $cdb->dbc->db_handle->selectrow_array(
"SELECT x.display_label FROM $type t, xref x where t.display_xref_id=x.xref_id and t.stable_id='${stable_id}'"
      );
  }

  return $self->{'display_name_cache'}->{$stable_id};
}

=head2 get_core_stable_id_by_display_name

  Args [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Args [2]   : display name (e.g. from core DB or GNC name)
  Example    : 
  Description: Builds a cache of stable ID to display names.
  Returntype : string - gene stable ID
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> display_name
# Return hashref keyed on {$type}{$stable_id}
#Need to update cache if we're doing more than one 'type' at a time
# as it will never get loaded for the new type!

sub get_core_stable_id_by_display_name {
  my ( $self, $cdb, $display_name ) = @_;

#if($type !~ /(gene|transcript|translation)/){
#	throw("Cannot get display_name for stable_id $stable_id with type $type");
#  }

  if ( !exists $self->{'stable_id_cache'}->{$display_name} ) {
    ( $self->{'stable_id_cache'}->{$display_name} ) =
      $cdb->dbc->db_handle->selectrow_array(
"SELECT g.stable_id FROM gene g, xref x where g.display_xref_id=x.xref_id and and x.display_label='${display_name}'"
      );
  }

  return $self->{'stable_id_cache'}->{$display_name};
}

#This could be a simple sub, if we passed $rollback_modes
#add flag for states comparison?

sub _compare_set_for_rollback {
  my ($self, $new_set, $stored_set, $set_type, $rollback_level, $slices) = @_;
  my $diffs = $stored_set->compare_to($new_set, undef, undef, undef, 1);
  #by default this tests nested objects via is_stored and dbID match
  #1 is skip states flag for data set

  #delete get_states_diffs as these are never set and will likely
  #always be different and we handle IMPORTED below
  delete $diffs->{'get_all_states'} if exists $diffs->{'get_all_states'};  
          
  if(%$diffs){
    $self->debug(1, 'Found diffs for stored vs new '.ref($new_set).' '.
      $new_set->name, $diffs);
    
    if($set_type eq 'data_set'){
      #data_set does not have a rollback_mode level as this is handled by feature_set
      throw('Found differences in specified and stored DataSet '.
            "support or product FeatureSet, please rectify manually\n".
            Dumper($diffs));      
    }
    elsif($rollback_level < $rollback_modes{$set_type}){
      throw("Found $set_type(".$new_set->name.') mismatch, please rectify manually or specify '.
            "-rollback $set_type\n".Dumper($diffs));
    }
    elsif($slices){
      #Should never have diffs and slices set, this indicates we are 
      #redefining the inputs but only rerunning a subset of the data 
      #There should be no diffs in parallel mode(single slice), as we should
      #have resolved this in the previous setup/submit analysis! 
      throw("It is unsafe to rollback $set_type with a sub set of slices ".
            'please do a full rollback i.e. omit -slices');
    }
  }
    
  return $diffs;
}



### DEPRECATED METHODS ###

1;
