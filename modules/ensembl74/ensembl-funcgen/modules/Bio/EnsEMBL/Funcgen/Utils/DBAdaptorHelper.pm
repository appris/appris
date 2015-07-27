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

Bio::EnsEMBL::Funcgen::Utils::DBAHelper

=head1 DESCRIPTION

This module provides a set of useful methods for creating Ensembl DBAdaptors.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Utils::DBAHelpers qw(list your required methods here);

=cut


# These methods were removed from EFGUtils do to a cyclical dependancy
# caused by the DBAdaptor using EFGUtils. Hence this is not to be used by
# the DBAdaptor itself, just as a way to generate them.

###############################################################################

package Bio::EnsEMBL::Funcgen::Utils::DBAdaptorHelper;

use strict;
use warnings;
use Getopt::Long                   qw( GetOptionsFromArray );
use Bio::EnsEMBL::Utils::Exception qw( throw               );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref          );
use Data::Dumper                   qw( Dumper              );
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; #This was causing cyclical requirements as it requires EFGUtils
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base qw( Exporter );
use vars   qw( @EXPORT_OK );

@EXPORT_OK = qw(               
                create_DBAdaptor_from_params
                create_DBAdaptor_from_options
                create_DBConnection_from_options
                create_Funcgen_DBAdaptor_from_options
                create_Funcgen_DBAdaptor_from_params
                get_MYSQL_args_from_DBAdaptor_params
                get_MYSQL_args_from_options
                process_DB_options
                get_DB_options_config
               );


#Keys here map to distinct DBAdaptors
my %db_type_options = ('funcgen' => ['user=s', 'pass=s', 'port=i', 'host=s', 'dbname=s'],   
                       'core'    => ['dnadb_user=s', 'dnadb_pass=s', 'dnadb_port=i', 'dnadb_host=s', 'dnadb_name=s'],
                       'hive'    => ['pdb_user=s', 'pdb_pass=s', 'pdb_port=i', 'pdb_host=s', 'pdb_name=s']);
#Should probably add in full support for other DBAdaptor params e.g. 
#driver, group, no_cache, 
#multi_species, species_id ??


#Add in some aliases
$db_type_options{dna}  = $db_type_options{core};
#These follwing two don't map to distinct DBAdaptors
#Hence will need manually will need manually re-assigning
#before use with create_DBAdaptors_from_params
$db_type_options{db}   = $db_type_options{funcgen};
$db_type_options{pipeline} = $db_type_options{hive}; 



#Split these out into:
#get_DB_options_config
#process_DB_options (to params or script/mysql options)
#get_MYSQL_args_from_options (this is now done by process_DB_options
#create_DBadaptors_from_params (already exists)

#Then have wrapper methods which
#create_Funcgen_DBAdaptor_from_options

=head2 get_DB_options_config

  Arg [1]    : Arrayref - DB types to handle e.g. funcgen|db, core|dna or 
               hive|pipeline
  Example    : my $db_opts = get_DB_options_config(['funcgen', 'core']);
               GetOptions(%$db_opts, ...);
               
               #OR
               
               my @db_opts_config = keys %{get_DB_options_config(['funcgen', 'core'])};
               my $db_opts = {...}; #references to vars/defaults/sub handlers in here
               GetOptions($db_opts, @db_opts_config);
               
               #Then get and adaptor with the options
               my $efg_db = create_Funcgen_DBAdaptor_from_options($db_opts);
  DESCRIPTION: Returns an options hash for use with GetOptions and other 'options'
               methods in this module
  Returntype : Hashref
  Exceptions : Throws if arguments are not valid.
               Throws if DB type is not valid.
  Caller     : General, scripts
  Status     : At risk

=cut 

#custom dbtypes cannot be used with process_DB_options (yet)
#shall we have an allow custom flag?


sub get_DB_options_config{
  my ($db_types, $allow_custom) = @_;  
  
  if(! defined $db_types){
    $db_types = ['funcgen', 'pipeline', 'core'];
  }
  elsif( (ref($db_types) ne 'ARRAY') ||
          (scalar(@$db_types) < 1) ){
    throw("DB types arg must be an arrayref of valid cmdline DB types:\t". 
      join("\t", keys(%db_type_options)) );            
  }
     
  my %db_opts = ();

  foreach my $db_type(@$db_types){
    
    if(! exists $db_type_options{$db_type}){
      
      if(! $allow_custom){
        throw("$db_type is not valid. Valid DB types are:\t". 
          join("\t", keys(%db_type_options)) ); 
      }
   
      map {my $param; $db_opts{$_} = \$param } 
        ( map {$db_type.$_ } ('_user=s', '_pass=s', '_port=i', '_host=s', '_name=s') );  
    }
    else{
      map {my $param; $db_opts{$_} = \$param } @{$db_type_options{$db_type}};  
    } 
  }
  
  
  
  
  return \%db_opts;
}


=head2 process_DB_options

  Arg [1]    : Hashref  - Options hash after GetOptions processing
  Arg [2]    : Arrayref - DB types to handle e.g. funcgen|db, core|dna or 
               hive|pipeline
  Arg [3]    : String - Validation mode
                 pass     - Adds 'pass' to the list of parameters to check
                 optional - Does not validate any 
               Omitting will validate user, host and dbname are defined.  
  Arg [4]    : String - option type:
                 dnadb  - Retains original options names and does not check mandatory 
                          params i.e. for use with Funcgen DBAdaptor.
                 script - Retains original option names and returns a String
                 mysql  - Converts to MySQL cmdline options and returns a String
  Example    : my $db_params = process_DB_iptions($db_opts, ['funcgen', 'core']);
  Description: Processes the hash of options returned by GetOptions, either into DBAdaptor
               parameters or maintains the original option names for use with scripts.
  Returntype : Hashref - DB type keys. Values can be:
                  DBAdaptor constructor parameter hashref
                  Script arguments string 
                  Mysql arguments string
  Exceptions : Throws if arguments are not valid.
               Throws if DB type is not valid.
               Throws if optional boolean is not specified and mandatory params
               for each given DB type are not found.
  Caller     : General, scripts.
  Status     : At risk

=cut

#Now supports both way of getting options hash
#change optional to support write(pass defined) and full optional param sets?
#validation level - optional|pass required (default no pass required)

sub process_DB_options {
  my ($db_opts, $db_types, $vlevel, $option_type) = @_;  
  
  if(! defined $db_types){
    $db_types = ['funcgen', 'pipeline', 'core'];
  }
  elsif( (ref($db_types) ne 'ARRAY') ||
          (scalar(@$db_types) < 1) ){
    throw("DB types arg must be an arrayref of valid cmdline DB types:\t". 
      join("\t", keys(%db_type_options)) );            
  }
     
  assert_ref($db_opts, 'HASH', 'db_opts');
                 
  if(defined $vlevel && 
     (($vlevel ne 'optional') && ($vlevel ne 'pass')) ){
    throw('Validation level argument is not valid, please specify \'optional\' or \'pass\'');     
  } 
  
  my %valid_option_types = (dnadb  => 1,
                            script => 1,
                            mysql  => {'host'   => '-h',
                                       'user'   => '-u',
                                       'pass'   => '-p',
                                       'port'   => '-P',
                                       'dbname' => ' '});
 
  if(defined $option_type && 
     (! exists $valid_option_types{$option_type}) ){
     throw("$option_type is not valid. Valid option types:\t".
      join("\t", keys %valid_option_types) ); 
  }
  
  my (%db_params, $param, $check_param);
  
  foreach my $db_type(@$db_types){
    
    if(exists $db_params{$db_type}){
      next;
      #We have a duplicate db_type and don't want to over-write
      #processed options
    }
    
    #Caller will expect a hashref even if empty 
    $db_params{$db_type} = {};
    
    
    if(! exists $db_type_options{$db_type}){
      throw("$db_type is not valid. Valid DB types are:\t". 
        join("\t", keys(%db_type_options)) ); 
    }
    
    #cache sub'd check params in case we are mainting options
    my %check_params = ('user'   => 0,
                        'dbname' => 0,
                        'host'   => 0);   
       
    if (defined $vlevel && ($vlevel eq 'pass')){                    
      $check_params{pass} = 0;   
    }
      
    my $dbtype_prefix = '';  
          
    foreach my $opt(@{$db_type_options{$db_type}}){
      my $deref = 1;
      my $value;
      ($param = '-'.$opt) =~ s/\=.*//o;  
      #This conditional sub/deref is to handle the two 
      #different types of opts hash we can take
      #$opt entry should already 'exist' from the first 
      #call of process_DB_options
      if(exists $db_opts->{$opt}){
        
        if(defined $db_opts->{$opt}){
          $value = ${$db_opts->{$opt}};      
        } 
      }
      elsif(exists $db_opts->{$param}){
        $value = $db_opts->{$param};
      }
       
                          
      if(defined $value){
        $param =~ /-(.*_)?(.*)/o;
        $check_param   = $2;
        $dbtype_prefix = $1 || ''; 
        
        #($check_param = $param) =~ s/-(.*_)*//o;           #Strip off -(pdb_|dnadb_) prefixes
        $check_param = 'dbname' if $check_param eq 'name'; #handle dnadb_name oddity  
        
        if (exists $check_params{$check_param}){
          $check_params{$check_param} = 1;  
        }
        
        #Do some more substitution or not
        if(! defined $option_type){
          $param =~ s/.*_/-/o;                     #Strip off -pdb_ or -dnadb_ prefixes       
          $param = '-dbname' if $param eq '-name'; #handle dbname oddity
        }
        elsif($option_type eq 'mysql'){
          #We always need to use the DBADaptor param here
          #rather any maintained script/dnadb option
          $param = $valid_option_types{$option_type}->{$check_param};
        }
    
        $db_params{$db_type}{$param} = $value;
      }
    }
   
    # CHECK PARAMS       
    if((! defined $vlevel) ||
       ($vlevel eq 'pass') ){
        
      if(grep(/0/, values %check_params) ){
        #This lists -dnadb_name option as -dnadb_dbname
        throw("Failed to validate $db_type DB options. Could not find mandatory (".
          join(' ', (map {"-${dbtype_prefix}".$_ } keys %check_params)).
          ") options:\n".Dumper($db_opts));        
      }
    }
    
    # BUILD STRING for cmdline
    if(defined $option_type &&
      ($option_type ne 'dnadb') ){
      $db_params{$db_type} = join(' ', (map {$_.' '.$db_params{$db_type}->{$_}} 
                                        keys %{$db_params{$db_type}}) );
    }
    
  }
  
 
  return \%db_params;
}  


=head2 create_Funcgen_DBAdaptor_from_options

  Arg [1]    : Hashref - DB options generated by passing the get_DB_options_config hash 
               to GetOptions
  Arg [2]    : String - Validation mode, 'pass' or 'optional'.
               Omitting will validate user, host and dbname are defined.
  Arg [3]    : Boolean - DNA DB validation, conditional on the Funcgen DB.          
  Example    : my $db = create_Funcgen_DBAdaptor_from_params($db_opts);
  DESCRIPTION: Creates a Funcgen DBAdaptor given the options hash returned from GetOptions.
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor 
  Exceptions : None
  Caller     : General, scripts
  Status     : At risk

=cut

#This currently doesn't support pass validation for dnadb

sub create_Funcgen_DBAdaptor_from_options {
  my ($db_opts, $vlevel, $validate_dnadb) = @_;
 
  my $dnadb_vlevel   = 'optional'; 
  my $funcgen_params = ${&process_DB_options($db_opts, ['funcgen'], $vlevel)}{funcgen};

  if(%{$funcgen_params} && $validate_dnadb){
    $dnadb_vlevel = undef;
  }

  my $core_params    = ${&process_DB_options($db_opts, ['core'], $dnadb_vlevel, 'dnadb')}{core};
  
  return create_DBAdaptor_from_params({%{$funcgen_params},
                                       %{$core_params}},
                                      'funcgen');  
}

=head2 create_DBAdaptor_from_options

  Arg [1]    : Hashref - DB options generated by passing the get_DB_options_config hash 
               to GetOptions
  Example    : my $db = create_Funcgen_DBAdaptor_from_params($db_opts);
  DESCRIPTION: Creates a Funcgen DBAdaptor given the options hash returned from GetOptions.
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor 
  Exceptions : None
  Caller     : General, scripts
  Status     : At risk

=cut


sub create_DBAdaptor_from_options {
  my ($db_opts, $db_type, $vlevel) = @_;
  $db_type ||= 'core'; #Enable creation of core adaptor for any db type options
  
  return create_DBAdaptor_from_params(${&process_DB_options($db_opts, [$db_type], $vlevel)}{$db_type},
                                      'core');  
}




=head2 create_Funcgen_DBAdaptor_from_params

  Arg [1]    : Hashref - DB type keys, DB contructor param hash values i.e. what 
               is returned from process_DB_options.
  Arg [2]    : Boolean - Flag to validate a -pass parameter has been set in the DB parameters   
  Example    : my $db = create_Funcgen_DBAdaptor_from_params($db_params);
  DESCRIPTION: Wrapper method to create a Funcgen DBAdaptor from a DB type => parameters hash.
               This will try and do the write thing if the core DB parameters are set i.e. 
               create a core DB directly or pass the dnadb parameters onto the Funcgen DBAdaptor.
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions : Throws if arguments not valid 
               Throws if Hash argument does not contain a funcgen key.
               Throws if 'core' Hash exists but does not contain mandatory DB parameters or
               and dnadb paramters
  Caller     : General
  Status     : At risk

=cut

sub create_Funcgen_DBAdaptor_from_params {
  my ($db_params, $pass_required) = @_;
  my ($dnadb, %dnadb);
  
  if(! exists $db_params->{funcgen}){
    throw('DB params argument does not contain a \'funcgen\' key');
  }
  
  if($pass_required && 
     (! exists $db_params->{funcgen}{-pass})){
    throw('The -pass boolean has been set but no -pass paramter was'.
      ' found in the DB parameters');     
  }
  
  if(exists $db_params->{core}){  
    #Check whether we have complete DBADaptor param first   
    eval {$dnadb = create_DBAdaptor_from_params($db_params->{core}, 'core'); };
    my $error = $@;
        
    #These checks are not exhaustive, but should catch most behaviour
    
    if($@){
      if( (exists $db_params->{core}->{'-user'}) ||
          (exists $db_params->{core}->{'-host'}) ||
          (exists $db_params->{core}->{'-dbname'}) ){
        throw("Some mandatory parameters were omited when creating the core DBAdaptor.\n$@");            
      }
            
      if(! grep(/-dnadb/, keys %{$db_params->{core}}) ){
        throw('DB parameter does not contain either valid core DBAdaptor params'.
          " or funcgen DBAdaptor dnadb params.\n".Dumper($db_params->{core}) );         
      }
    }
  }
   
  if(defined $dnadb){
    $dnadb{-dnadb} = $dnadb;  
  }
  elsif(exists $db_params->{core}){
    %dnadb = %{$db_params->{core}};
  }     
    
  return create_DBAdaptor_from_params({%{$db_params->{funcgen}}, %dnadb}, 'funcgen');
}


=head2 create_DBAdaptor_from_params

  Arg [1]    : Hashref - DB contructor params
  Arg [2]    : String  - DB type e.g. funcgen, core or hive
  Example    : my $db = create_DBAdaptor_from_params($funcgen_params, 'funcgen');
  DESCRIPTION: Creates and tests the relevant DBAdaptors from the passed parameters.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor (Or Funcgen or Hive DBAdaptor)
  Exceptions : Throws if arguments not valid 
               Throws if DB type is not supported
  Caller     : General
  Status     : At risk

=cut
    
sub create_DBAdaptor_from_params {
  my ($db_params, $db_type) = @_;
   
  assert_ref($db_params, 'HASH', 'db_params');
  my %dba_modules = (
                     funcgen => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
                     core    => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
                     hive    => undef, #'Bio::EnsEMBL::Hive::DBSQL::DBAdaptor',
                     #This adds a unecessary requirement for ensembl-hive package
                     );
  
             
  if( scalar(keys %$db_params) < 1 ){
    throw('You must pass a hashref DB parameters.');
  }  
  
  if(! (defined $db_type && 
        exists $dba_modules{$db_type}) ){
    throw("DB type argument is invalid, valid DB types are:\t".join("\t", (keys %dba_modules)) );        
  }
  
  if(! defined $dba_modules{$db_type}){
    my $dba_class = 'Bio::EnsEMBL::'.ucfirst($db_type).'::DBSQL::DBAdaptor';
    eval "require $dba_class";
    
    if($@){
      throw($@."\nFailed to require DBAdaptor class:\t$dba_class");
    }    
    
    $dba_modules{$db_type} = $dba_class;
  }
  
  my $dba = $dba_modules{$db_type}->new(%$db_params);   
  
  #Test connections  
  $dba->dbc->db_handle;
         
  if($db_type eq 'funcgen'){
    $dba->dnadb->dbc->db_handle; 
  }
  
  return $dba; 
}    
 




sub create_DBConnection_from_options {
  my ($db_opts, $db_type, $vlevel) = @_;
  return create_DBConnection_from_params(${process_DB_options($db_opts, [$db_type], $vlevel)}{$db_type});
}


sub create_DBConnection_from_params {
  my $db_params    = shift;
  my $extra_params = shift;
  $extra_params  ||= {};
  assert_ref($extra_params, 'HASH', 'extra_params');
  assert_ref($db_params, 'HASH', 'db_params');
          
  if( scalar(keys %$db_params) < 1 ){
    throw('You must pass a hashref DB parameters.');
  }  

  my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$db_params, %$extra_params);   
  
  #Test connection  
  $dbc->db_handle;
         
  return $dbc; 
}   


=head2 get_MYSQL_args_from_options

  Arg [1]    : Hashref - DB options generated by passing the get_DB_options_config hash 
               to GetOptions
  Arg [2]    : String - DB type e.g. funcgen, core  
  Example    : my $mysql_args = get_MYSQL_args_from_options($db_opts, 'funcgen');
  DESCRIPTION: Wrapper method to process_DB_options.
  Returntype : String
  Exceptions : None
  Caller     : General, scripts
  Status     : At risk

=cut 

sub get_MYSQL_args_from_options {
  my ($db_opts, $db_type, $optional) = @_;

  return ${&process_DB_options($db_opts, [$db_type], $optional, 'mysql')}{$db_type};
}


=head2 get_MYSQL_args_from_DBAdaptor_params

  Arg [1]    : Hashref - DBAdaptor constructor params
  Example    : my $mysql_args = get_MYSQL_args_from_options(\%funcgen_db_params);
  DESCRIPTION: Builds MYSQL cmdline argument string based on DBAdaptor parameters
  Returntype : String
  Exceptions : None
  Caller     : General, scripts
  Status     : At risk

=cut 


sub get_MYSQL_args_from_DBAdaptor_params {
  my $db_params = $_[0];    
  assert_ref($db_params, 'HASH', 'db_params');
  
  my @db_mysql_params = (['-host',  '-h', 1],
                         ['-user',  '-u', 1],
                         ['-pass',  '-p'],
                         ['-port',  '-P'],
                         ['-dbname', ' ', 1]);
  
  my $mysql_params= '';
  
  foreach my $param_info(@db_mysql_params){
    my ($db_param, $mysql_param, $mandatory) = @$param_info;
    
    my $param = '';
    
    if(! exists ${$db_params}{$db_param}) {
      if($mandatory) {
        throw($db_param.' is a mandatory DB param');  
      }
    }
    else{
      $param .= ' '.$mysql_param.$db_params->{$db_param};  
    }
    
    $mysql_params .= $param;
  }  
  
  return $mysql_params;
}
   


1;
