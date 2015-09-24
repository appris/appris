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

Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor -
MySQL Database queries to load and store external object references.

=head1 SYNOPSIS

$db_entry_adaptor = $db_adaptor->get_DBEntryAdaptor();
$db_entry = $db_entry_adaptor->fetch_by_dbID($id);

my $gene = $db_adaptor->get_GeneAdaptor->fetch_by_stable_id('ENSG00000101367');
@db_entries = @{$db_entry_adaptor->fetch_all_by_Gene($gene)};
@gene_ids = $db_entry_adaptor->list_gene_ids_by_extids('BAB15482');

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use base qw( Bio::EnsEMBL::DBSQL::DBEntryAdaptor Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType $feature_type 
               (The feature type to retrieve DBEntries for)
  Arg [2]    : optional external database name
  Arg [3]    : optional external_db type 
  Example    : @db_entries = @{$db_entry_adaptor->fetch_by_FeatureType($feature_type)};
  Description: This returns a list of DBEntries associated with this feature type.
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or GoXref if there is linkage data.
  Exceptions : throws if feature type object not passed
  Caller     : Bio::EnsEMBL::Funcgen::FeatureType
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType {
  my ( $self, $feature_type, $ex_db_reg, $exdb_type ) = @_;
  assert_ref($feature_type, 'Bio::EnsEMBL::Funcgen::FeatureType');
  return $self->_fetch_by_object_type($feature_type->dbID, 'FeatureType', $ex_db_reg, $exdb_type);
}

# TODO - Add more fetch_all_by methods for other classes(e.g. Probe & ProbeSet etc.),
# or make the above generic as this is very rarely used. This could utilise a class
# validation hash which could be generated dynamically in new via the available adaptors
# and a little bit of logic. This would avoid missing support for new classes.

# NOTE: There is a case to remove all these _ids_by_extids wrappers
# There origins are in the core API where is is necessary to have wrappers
# to call related methods i.e. for Gene methods the API also calls Transcript 
# and Translation methods.
# There are no obvious relationships like this in Funcgen, with the potential
# exception of Probe & ProbeSet.
# The only usage of these wrappers was in BaseAdaptor::fetch_all_by_external_name
# which has now been optimised to use the _type_by_external_ids method directly


=head2 list_feature_type_ids_by_extids

  Arg [1]    : ARRAYREF of external name strings
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_feature_type_ids_by_extids(['ENST00012398371', ...]);
  Description: Gets a list of feature type IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_feature_type_ids_by_extids {
  my ( $self, $external_names, $external_db_name ) = @_;
  return $self->_type_by_external_ids( $external_names, 'FeatureType', 
									                     undef, $external_db_name );
}


=head2 list_regulatory_feature_ids_by_extids

  Arg [1]    : ARRAYREF of external name strings
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_regulatory_feature_ids_by_extids(['ENSG00283757289', ...]);
  Description: Gets a list of regulatory_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_regulatory_feature_ids_by_extids {
  my ( $self, $external_names, $external_db_name ) = @_;
  return $self->_type_by_external_ids( $external_names, 'RegulatoryFeature', 
									                     undef, $external_db_name );
}


=head2 list_external_feature_ids_by_extids

  Arg [1]    : ARRAYREF of external name strings
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_external_feature_ids_by_extids('ENSG00085672387',...]);
  Description: Gets a list of external_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_external_feature_ids_by_extids {
  my ( $self, $external_names, $external_db_name ) = @_;
  return $self->_type_by_external_ids( $external_names, 'ExternalFeature', 
                                       undef, $external_db_name );
}


=head2 list_annotated_feature_ids_by_extids

  Arg [1]    : ARRAYREF of external name strings
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_annotated_feature_ids_by_extids('ENSG00023847582', ...]);
  Description: Gets a list of annotated_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut


sub list_annotated_feature_ids_by_extids{
  my ( $self, $external_names, $external_db_name ) = @_;
  return $self->_type_by_external_ids( $external_names, 'AnnotatedFeature', 
                                  undef, $external_db_name );
}


=head2 list_probe_feature_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_annotated_feature_ids_by_extid('ENST000000000001');
  Description: Gets a list of annotated_feature IDs by external display ID
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_probe_feature_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;
  return $self->_type_by_external_id( $external_name, 'ProbeFeature', 
                                      undef, $external_db_name );
}


=head2 list_probe_feature_ids_by_extids

  Arg [1]    : ARRAYREF of external name strings
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_annotated_feature_ids_by_extids(['ENST000000000001', ...]);
  Description: Gets a list of annotated_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_probe_feature_ids_by_extids {
  my ( $self, $external_names, $external_db_name ) = @_;
  return $self->_type_by_external_ids( $external_names, 'ProbeFeature', 
                                      undef, $external_db_name );
}


=head2 list_probe_ids_by_extids

  Arg [1]    : ARRAYREF of external name strings
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_probe_id_by_extids(['ENST000000000001'], ...);
  Description: Gets a list of probe IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_probe_ids_by_extids {
  my ( $self, $external_names, $external_db_name ) = @_;
  return $self->_type_by_external_ids( $external_names, 'Probe', 
                                       undef, $external_db_name );
}


=head2 list_probeset_ids_by_extids

  Arg [1]    : ARRAYREF of external name strings
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_probeset_ids_by_extids(['ENST000000000001'], ...);
  Description: Gets a list of probeset IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_probeset_ids_by_extids {
  my ( $self, $external_names, $external_db_name ) = @_;
  return $self->_type_by_external_ids( $external_names, 'ProbeSet', 
                                       undef, $external_db_name );
}


=head2 _type_by_external_id

  Arg [1]    : string $name - dbprimary_acc
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  Arg [4]    : (optional) string $external_db_name
  	       other object type to be returned
  Example    : $self->_type_by_external_id($name, 'regulatory_feature');
  Description: Gets
  Returntype : list of dbIDs (regulatory_feature, external_feature )
  Exceptions : none
  Caller     : list_regulatory/external_feature_ids_by_extid
  Status     : Stable

=cut

#NOTE: Can't deprecate/remove this as it is part of the core API interface

sub _type_by_external_id {
  my ( $self, $name, $ensType, $extraType, $external_db_name ) = @_;
  return $self->_type_by_external_ids([$name], $ensType, $extraType, $external_db_name);
} ## end sub _type_by_external_id


=head2 _type_by_external_ids

  Arg [1]    : ARRAYREF of external names strings
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  Arg [4]    : (optional) string $external_db_name
  	            other object type to be returned
  Example    : $self->_type_by_external_ids([$name, $name2] 'regulatory_feature');
  Description: Gets
  Returntype : list of dbIDs (regulatory_feature, external_feature )
  Exceptions : none
  Caller     : list_regulatory/external_feature_ids_by_extid
  Status     : Stable

=cut

sub _type_by_external_ids {
  my ( $self, $names, $ensType, $extraType, $external_db_name ) = @_;

  my $from_sql  = '';
  my $where_sql = '';
  my $ID_sql    = "oxr.ensembl_id";

  if ( defined $extraType ) {
    #This was DBLinks query? We could do this for ProbeSet->Probe->ProbeFeature?
    #See core method for missing code
    throw('Extra types not accomodated in Funcgen xref API');
  }

  if ( defined($external_db_name) ) {
    # Involve the 'external_db' table to limit the hits to a particular
    # external database.
    $from_sql .= 'external_db xdb, ';
    $where_sql .= 'xdb.db_name LIKE '.$self->dbc->db_handle->quote( $external_db_name . '%' )
     .' AND xdb.external_db_id = x.external_db_id AND';
   }


   my $in_clause = '('.join(', ', (map $self->db->dbc->db_handle->quote($_), @$names)).')';
  #For use with use selectcol_arrayref 

  #my $in_clause = '('.join(', ', (('?') x scalar(@$names))).')';
  #For use with fetchall_hashref

  my @queries = 
   ("SELECT $ID_sql FROM $from_sql xref x, object_xref oxr ".
    "WHERE $where_sql x.dbprimary_acc IN $in_clause AND ".
    "x.xref_id = oxr.xref_id AND oxr.ensembl_object_type= '${ensType}'",
    "SELECT $ID_sql FROM $from_sql xref x, object_xref oxr ".
    "WHERE $where_sql x.display_label IN $in_clause AND ".
    "x.xref_id = oxr.xref_id AND oxr.ensembl_object_type= '${ensType}'");

  if ( defined $external_db_name ) {
    # If we are given the name of an external database, we need to join
    # between the 'xref' and the 'object_xref' tables on 'xref_id'.

    push @queries, "SELECT $ID_sql FROM $from_sql xref x, object_xref oxr, external_synonym syn ".
     "WHERE $where_sql syn.synonym IN $in_clause AND x.xref_id = oxr.xref_id AND ".
     "oxr.ensembl_object_type= '${ensType}' AND syn.xref_id = oxr.xref_id";
  } else {
    # If we weren't given an external database name, we can get away
    # with less joins here.
    push @queries, "SELECT $ID_sql FROM $from_sql object_xref oxr, external_synonym syn ".
     "WHERE $where_sql syn.synonym IN $in_clause AND ".
     "oxr.ensembl_object_type= '${ensType}' AND syn.xref_id = oxr.xref_id";
   }

  # Increase speed of query by splitting the OR in query into three
  # separate queries.  This is because the 'or' statments renders the
  # index useless because MySQL can't use any fields in it.
  # Changed this to a UNION and grab the col arrayref directly
  #UNION will return NR row but ensembl_ids may have >1 ox row
  #hence selectcol won't  necessarily return NR ID list, 
  #as caller will use fetch_all_by_dbID_list which will make it NR
  #change to selectall_hashref if this changes

  return @{$self->db->dbc->db_handle->selectcol_arrayref(join(' UNION ', @queries))};
} ## end sub _type_by_external_ids


=head2 _type_by_external_db_id

  Arg [1]    : string $type - external_db type
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  	       other object type to be returned
  Example    : $self->_type_by_external_id('1030', 'Translation');
  Description: Gets
  Returntype : list of dbIDs (gene_id, transcript_id, etc.)
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  	       geneids_by_extids
  Status     : Stable

=cut

#Use split_CamelCase here and a validation hash to remove this elsif trail
#This could be generated in new from the available adaptor list and a little logic
#is there any way we can dynamically build this validation list?
#and could be used to validate on store too

#Don't actaully use this, but we need to redefine it to over-ride the core specific method

sub _type_by_external_db_id{
  my ($self, $external_db_id, $ensType, $extraType) = @_;
  my $from_sql  = '';
  my $where_sql = '';
  my $ID_sql    = "oxr.ensembl_id";

  if (defined $extraType) { #See core method for missing code
    throw('Extra types not accomodated in eFG xref schema');
  }

  if(lc($ensType) eq 'regulatoryfeature'){
    $from_sql  = 'regulatory_feature rf, ';
    $where_sql = qq( rf.regulatory_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'externalfeature'){
    $from_sql  = 'external_feature ef, ';
    $where_sql = qq( ef.external_feature_id = oxr.ensembl_id AND );
  } 
  elsif(lc($ensType) eq 'annotatedfeature'){
    $from_sql  = 'annotated_feature af, ';
    $where_sql = qq( af.annotated_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'featuretype'){
    $from_sql  = 'featuretype ft, ';
    $where_sql = qq( ft.feature_type_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probefeature'){
    $from_sql  = 'probe_feature pf, ';
    $where_sql = qq( pf.probe_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probe'){
    $from_sql  = 'probe p, ';
    $where_sql = qq( p.probe_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probeset'){
    $from_sql  = 'probe_set ps, ';
    $where_sql = qq( pf.probe_set_id = oxr.ensembl_id AND );
  }
  else{
    throw("ensembl_object_type $ensType is not accommodated");
  }

  my $query = "SELECT $ID_sql FROM $from_sql xref x, object_xref oxr ".
   "WHERE $where_sql x.external_db_id = ? AND x.xref_id = oxr.xref_id ".
   "AND oxr.ensembl_object_type= ?";

  # Increase speed of query by splitting the OR in query into three separate 
  # queries. This is because the 'or' statments render the index useless 
  # because MySQL can't use any fields in the index.
  my (%hash, @result);
  my $sth = $self->prepare( $query );
  $sth->bind_param(1, "$external_db_id", SQL_VARCHAR);
  $sth->bind_param(2, $ensType, SQL_VARCHAR);
  $sth->execute();

  while( my $r = $sth->fetchrow_array() ) {
    if( !exists $hash{$r} ) {
      $hash{$r} = 1;
      push( @result, $r );
    }
  }

  return @result;
}


# NOTE: list_regulatory_feature_ids_by_external_db_id is never 
# used and never likely to be but there as an example of the
# _type_by_external_db_id implementation above

sub list_regulatory_feature_ids_by_external_db_id{
  my ($self,$external_db_id) = @_;
 
  my %T = map { ($_, 1) }
   $self->_type_by_external_db_id( $external_db_id, 'RegulatoryFeature' );
  
  return keys %T;
}


#Placeholders to catch error from inherited methods
#These now work in reverse as the Gene/Transcript/Translation 
#is the xref not the ensembl_object as with the core code

sub fetch_all_by_Gene {  throw('Not implemented for Funcgen API, maybe you want the core DBEntryAdaptor?'); }
sub fetch_all_by_Transcript { throw('Not implemented for Funcgen API, maybe you want the core DBEntryAdaptor?'); }
sub fetch_all_by_Translation {  throw('Not implemented for Funcgen API, maybe you want the core DBEntryAdaptor?'); }
  

### DEPRECATED ###

#These wrappers are not really need at all and have been replaced by the _ids_by_extids
#wrappers
#They were inherited from core where is makes more sense as the
#Gene and Transcript methods also assitionally call the Transcript and 
#Transclation method accordingly.
#Here they simply harcode the type once, instead of multiple times
#and also allow callers to validate the type by the presence of the
#method. 
#Down sides of removing them is potential for typos in caller
#This impact is minimal in the context of the test suite, as caller are
#generally internal methods

#Likewise the list_TABLE_NAME_ids_by_extids could also be deprecated
#unless there is a use case similar to the core for a Gene method calling
#related Transcript and Translation methods 
#i.e. Probe calling ProbeSet methods or vice-versa

sub list_regulatory_feature_ids_by_extid{ #DEPRECATED in v79
  my ( $self, $external_name, $external_db_name ) = @_;
  deprecate('Please use _type_by_external_ids directly for internal methods');
  return $self->_type_by_external_id( $external_name, 'RegulatoryFeature', 
                    undef, $external_db_name );
}

sub list_probeset_ids_by_extid{ #DEPRECATED in v79
  my ( $self, $external_name, $external_db_name ) = @_;
 deprecate('Please use _type_by_external_ids directly for internal methods');
  return
    $self->_type_by_external_id( $external_name, 'ProbeSet', undef,
                                 $external_db_name );
}

sub list_feature_type_ids_by_extid{ #DEPRECATED in v79
  my ( $self, $external_name, $external_db_name ) = @_;
 deprecate('Please use _type_by_external_ids directly for internal methods');
  return $self->_type_by_external_id( $external_name, 'FeatureType', 
                    undef, $external_db_name );
}

sub list_external_feature_ids_by_extid{ #DEPRECATED in v79
  my ( $self, $external_name, $external_db_name ) = @_;
 deprecate('Please use _type_by_external_ids directly for internal methods');
  return
    $self->_type_by_external_id( $external_name, 'ExternalFeature', undef,
                                 $external_db_name );
}

sub list_annotated_feature_ids_by_extid { #DEPRECATED in v79
  my ( $self, $external_name, $external_db_name ) = @_;
 deprecate('Please use _type_by_external_ids directly for internal methods');
  return
    $self->_type_by_external_id( $external_name, 'AnnotatedFeature', undef,
                                 $external_db_name );
}

sub list_probe_ids_by_extid { #DEPRECATED in v79
  my ( $self, $external_name, $external_db_name ) = @_;
 deprecate('Please use _type_by_external_ids directly for internal methods');
  return
    $self->_type_by_external_id( $external_name, 'Probe', undef,
                                 $external_db_name );
}





1;

