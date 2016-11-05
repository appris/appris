=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Registry

=head1 SYNOPSIS

  $registry = APPRIS::Registry->new(
    -dbhost  => 'localhost',
    -dbname  => 'homo_sapiens_encode_3c',
    -dbuser  => 'jmrodriguez'
    );

  $gene = $registry->fetch_by_stable_id($stable_id);

  @genes = @{ $registry->fetch_by_region('X', 1, 10000) };

=head1 DESCRIPTION

All Adaptors are stored/registered using this module.
This module should then be used to get the adaptors needed.

The registry can be loaded from a configuration file or from database info.

=head1 METHODS

=cut

package APPRIS::Registry;

use strict;
use warnings;
use Data::Dumper;
use Config::IniFiles;

use APPRIS::Gene;
use APPRIS::Transcript;
use APPRIS::Translation;
use APPRIS::Exon;
use APPRIS::CDS;
use APPRIS::Codon;
use APPRIS::XrefEntry;
use APPRIS::Analysis;
use APPRIS::Analysis::Region;
use APPRIS::Analysis::Firestar;
use APPRIS::Analysis::Matador3D;
use APPRIS::Analysis::SPADE;
use APPRIS::Analysis::INERTIA;
use APPRIS::Analysis::CRASH;
use APPRIS::Analysis::THUMP;
use APPRIS::Analysis::CORSAIR;
use APPRIS::Analysis::PROTEO;
use APPRIS::Analysis::APPRIS;
use APPRIS::DBSQL::DBAdaptor;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

my ($API_VERSION) = 'rel15';

{
    # Encapsulated class data
    #___________________________________________________________
    my %_attr_data = # DEFAULT
		(
			dbhost			=>  undef,
			dbport			=>  '3306',
			dbname			=>  undef,
			dbuser			=>  undef,
			dbpass			=>  '',
			dbadaptor		=>  undef,
		);
    #_____________________________________________________________

	# Classwide default value for a specified object attribute
	sub _default_for {
		my ($self, $attr) = @_;
		$_attr_data{$attr};
	}

	# List of names of all specified object attributes
	sub _standard_keys {
		keys %_attr_data;
	}
	
	sub dbhost {
		my ($self, $arg) = @_;
		$self->{dbhost} = $arg if defined $arg;
		return $self->{dbhost};
	}

	sub dbport {
		my ($self, $arg) = @_;
		$self->{dbport} = $arg if defined $arg;
		return $self->{dbport};
	}

	sub dbname {
		my ($self, $arg) = @_;
		$self->{dbname} = $arg if defined $arg;
		return $self->{dbname};
	}

	sub dbuser {
		my ($self, $arg) = @_;
		$self->{dbuser} = $arg if defined $arg;
		return $self->{dbuser};
	}

	sub dbpass {
		my ($self, $arg) = @_;
		$self->{dbpass} = $arg if defined $arg;
		return $self->{dbpass};
	}
	
	sub dbadaptor {
		my ($self, $arg) = @_;
		$self->{dbadaptor} = $arg if defined $arg;
		return $self->{dbadaptor};
	}
}

=head2 new

  Example :

    APPRIS::Registry->new()

  Description: Will load the correct versions of the appris
               databases for the software release it can find on a
               database instance into the registry.

  Exceptions : None.
  Status     : Stable

=cut

sub new {
	my ($caller, %args) = @_;
	
	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;

	foreach my $attrname ($self->_standard_keys) {
		my ($attr) = "-".$attrname;
		if (exists $args{$attr} && defined $args{$attr}) {
			$self->{$attrname} = $args{$attr};
		} else {
			$self->{$attrname} = $self->_default_for($attrname);
		}
	}

	return $self;
}

=head2 software_version
  
  get the software version.
  
  Args       : none
  ReturnType : int
  Status     : At Risk
  
=cut
  
sub software_version{
	my ($self) = @_;
	return $API_VERSION;
}

=head2 load_registry_from_db

  Arg [dbhost] : string
                The domain name of the database host to connect to.

  Arg [dbport] : string (Optional)
                The port to use when connecting to the database.

  Arg [dbname] : string
                The name of database to connect to.

  Arg [dbuser] : string
                The name of the database user to connect with.

  Arg [dbpass] : string
                The password to be used to connect to the database.

  Example :

    $registry->load_registry_from_db(
       -dbhost  => 'localhost',
       -dbname  => 'homo_sapiens_encode_3c',
       -dbuser  => 'jmrodriguez'
       -dbpass  => ''
    );

  Description: Will load the correct versions of the appris
               databases for the software release it can find on a
               database instance into the registry.

  Exceptions : None.
  Status     : Stable

=cut

sub load_registry_from_db {
	my ($self) = shift;
	
	my (
		$dbhost,	$dbport,
		$dbname,	$dbuser,
		$dbpass
    )
    = rearrange( [
		'dbhost',	'dbport',
		'dbname',	'dbuser',
		'dbpass'
	],
	@_
	);
	return undef
		unless (
			defined $dbhost and defined $dbname and
			defined $dbuser	and defined $dbpass 
	);
	
 	$self->dbhost($dbhost);
 	$self->dbname($dbname);
 	$self->dbuser($dbuser);
 	$self->dbpass($dbpass);

	return $self;
}

=head2 load_registry_from_ini

  Arg [file] : string
                Config ini file to connect to database.

  Example :

    $registry->load_registry_from_ini(
       -file  => 'config.ini'
    );

  Description: Will load the correct versions of the appris
               databases for the software release it can find on a
               database instance into the registry.
               The ini file has to be like:
                  [APPRIS_DB]
                    host=localhost
                    db=
                    user=
                    pass=
                    port=

  Exceptions : None.
  Status     : Stable

=cut

sub load_registry_from_ini {
	my ($self) = shift;

	my (
		$file
    )
    = rearrange( [
		'file'
	],
	@_
	);
	return undef
		unless (
			defined $file 
	);

	my ($cfg) = new Config::IniFiles(
		-file =>  $file
	); 
	my ($dbhost) = $cfg->val( 'APPRIS_DB', 'host');
	my ($dbname) = $cfg->val( 'APPRIS_DB', 'db' );
	my ($dbport) = $cfg->val( 'APPRIS_DB', 'port');
	my ($dbuser) = $cfg->val( 'APPRIS_DB', 'user');
	my ($dbpass) = $cfg->val( 'APPRIS_DB', 'pass');
	return undef
		unless (
			defined $dbhost and defined $dbname and
			defined $dbuser	and defined $dbpass 
	);
	
 	$self->dbhost($dbhost);
 	$self->dbname($dbname);
 	$self->dbuser($dbuser);
 	$self->dbpass($dbpass);

	return $self;
}

=head2 fetch_basic_by_region

  Arg [1]    : string $chr
               The name of the sequence region that the slice will be
               created on.
  Arg [2]    : int $start (optional)
               The start of the slice on the sequence region
  Arg [3]    : int $end (optional)
               The end of the slice on the sequence region
  Example    : $gene = $registry->fetch_basic_by_region('M', 1, 1000);
  Description: Retrieves a listref of gene object 
               from the database via its stable id.
               If the regions is not found undef is returned instead.
  Returntype : Listref of APPRIS::Gene or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_basic_by_region {
	my ($self, $chr, $start, $end) = @_;

	my ($entity) = $self->fetch_entity_by_region($chr, $start, $end);
	return $entity;
}

=head2 fetch_by_region

  Arg [1]    : string $chr
               The name of the sequence region that the slice will be
               created on.
  Arg [2]    : int $start (optional)
               The start of the slice on the sequence region
  Arg [3]    : int $end (optional)
               The end of the slice on the sequence region
  Example    : $gene = $registry->fetch_by_region('M', 1, 1000);
  Description: Retrieves a listref of gene object 
               from the database via its stable id.
               If the regions is not found undef is returned instead.
               Note: Entity object contains the results of analysis methods.
  Returntype : Listref of APPRIS::Gene or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_region {
	my ($self, $chr, $start, $end, $sources) = @_;
	
	# Flag of 'get_analysis' is active. By default 'all'
	unless ( defined $sources ) { $sources = 'all' }

	my ($entity) = $self->fetch_entity_by_region($chr, $start, $end, $sources);
	return $entity;
}

=head2 fetch_entity_by_region

  Arg [1]    : string $chr
               The name of the sequence region that the slice will be
               created on.
  Arg [2]    : int $start (optional)
               The start of the slice on the sequence region
  Arg [3]    : int $end (optional)
               The end of the slice on the sequence region
  Arg [4]    : Int $analysis (optional)
               Flag that controls if the analysis are reported
  Example    : $gene = $registry->fetch_entity_by_region('M', 1, 1000);
  Description: Retrieves a listref of gene object 
               from the database via its stable id.
               If the regions is not found undef is returned instead.
               Note: Entity object contains the results of analysis methods.
  Returntype : Listref of APPRIS::Gene or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_entity_by_region {
	my ($self, $chr, $start, $end, $sources) = @_;

	my ($entity_list);
	
	# Connection to DBAdaptor
	unless ( defined $self->dbadaptor ) {
		my ($dbadaptor) = APPRIS::DBSQL::DBAdaptor->new
		(
			dbhost => $self->dbhost,
			dbport => $self->dbport,
			dbname => $self->dbname,
			dbuser => $self->dbuser,
			dbpass => $self->dbpass,
		);
		$self->dbadaptor($dbadaptor);	
	}

	# Get gene list
	my ($gene_list);
	eval {
		my ($datasource) = $self->dbadaptor->query_datasource(name => 'Gene_Id')->[0];
		my (@args);
		push(@args, {
			datasource_id=> $datasource->{'datasource_id'}
		});
		if (defined $chr) {
			push(@args, {
				chromosome	=> $chr
			});			
		}
		if (defined $start and defined $end) {
			push(@args, [
				{
					start	=> $start,
					end		=> $end,
				},
				{
					start	=> $start,
					end		=> $end,
				},			
			]);			
		}
		$gene_list = $self->dbadaptor->query_entity_coordinate5(@args);
		warning('Argument must be a correct region') unless (defined $gene_list and scalar(@{$gene_list}) > 0);
	};
	warning('Argument must be a correct region') if ($@);
	
	# Do commit of everything
	$self->dbadaptor->commit();

	# Get every gene
	foreach my $gene (@{$gene_list})
	{
		my ($region);
		if (defined $start and defined $end) {
			$region = {
				start	=> $start,
				end		=> $end,
			};
		}
		my ($entity) = $self->fetch_by_gene_stable_id($gene->{'identifier'},$sources,$region);
		push(@{$entity_list},$entity) if (defined $entity);
	}
	
	return $entity_list;
}

=head2 fetch_basic_by_position

  Arg [1]    : string $chr
               The name of the sequence region that the slice will be
               created on.
  Arg [2]    : int $start
               The start of the slice on the sequence region
  Arg [3]    : int $end (optional)
               The end of the slice on the sequence region
               If is not defined, the position starts and ends in the 
               given 'start' position
  Arg [4]    : string $type (optional)
               Type of entity: gene or transcript. By default are both
  Example    : $entity = $registry->fetch_basic_by_position('M', 1000);
  Description: Retrieves a listref of entity object 
               from the database via its stable id.
               If the regions is not found undef is returned instead.
  Returntype : Listref of APPRIS::Gene and APPRIS::Transcript or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_basic_by_position {
	my ($self, $chr, $start, $end, $type) = @_;

	my ($entity) = $self->fetch_entity_by_position($chr, $start, $end);
	return $entity;
}

=head2 fetch_by_position

  Arg [1]    : string $chr
               The name of the sequence region that the slice will be
               created on.
  Arg [2]    : int $start
               The start of the slice on the sequence region
  Arg [3]    : int $end (optional)
               The end of the slice on the sequence region
               If is not defined, the position starts and ends in the 
               given 'start' position
  Arg [4]    : string $type (optional)
               Type of entity: gene or transcript. By default are both               
  Example    : $gene = $registry->fetch_by_position('M', 1000);
  Description: Retrieves a listref of entity objects 
               from the database via its stable id.
               If the regions is not found undef is returned instead.
               Note: Entity object contains the results of analysis methods.
  Returntype : Listref of APPRIS::Gene and APPRIS::Transcript or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_position {
	my ($self, $chr, $start, $end, $sources) = @_;

	# Flag of 'get_analysis' is active. By default 'all'
	unless ( defined $sources ) { $sources = 'all' }
	
	my ($entity) = $self->fetch_entity_by_position($chr, $start, $end, $sources);
	return $entity;
}

=head2 fetch_entity_by_position

  Arg [1]    : string $chr
               The name of the sequence region that the slice will be
               created on.
  Arg [2]    : int $start
               The start of the slice on the sequence region
  Arg [3]    : int $end (optional)
               The end of the slice on the sequence region
               If is not defined, the position starts and ends in the 
               given 'start' position
  Arg [4]    : Int $analysis (optional)
               Flag that controls if the analysis are reported               
  Arg [5]    : string $type (optional)
               Type of entity: gene or transcript. By default are both               
  Example    : $gene = $registry->fetch_entity_by_position('M', 1000);
  Description: Retrieves a listref of entity objects 
               from the database via its stable id.
               If the regions is not found undef is returned instead.
               Note: Entity object contains the results of analysis methods.
  Returntype : Listref of APPRIS::Gene and APPRIS::Transcript or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_entity_by_position {
	my ($self, $chr, $start, $end, $sources, $type) = @_;

	my ($entity_list);
	
	# Connection to DBAdaptor
	unless ( defined $self->dbadaptor ) {
		my ($dbadaptor) = APPRIS::DBSQL::DBAdaptor->new
		(
			dbhost => $self->dbhost,
			dbport => $self->dbport,
			dbname => $self->dbname,
			dbuser => $self->dbuser,
			dbpass => $self->dbpass,
		);
		$self->dbadaptor($dbadaptor);	
	}

	# Get list of several entities
	my ($aux_entity_list);
	eval {
		my (%region) = (
				chromosome		=> $chr
		);
		$region{'start'} = $start if (defined $start);
		if (defined $end) {
			$region{'end'} = $end;
		} else {
			$region{'end'} = $start;
		}
		
		$aux_entity_list = $self->dbadaptor->query_entity_coordinate4(%region);
		return undef unless (defined $aux_entity_list and scalar(@{$aux_entity_list}) > 0);
	};
	throw('Argument must be a correct position') if ($@);

	# Get gene and trascript list
	my ($gene_datasource_id);
	my ($transcript_datasource_id);
	eval {
		my ($datasource) = $self->dbadaptor->query_datasource(name => 'Gene_Id')->[0];
		$gene_datasource_id = $datasource->{'datasource_id'} if (defined $datasource);
		my ($datasource2) = $self->dbadaptor->query_datasource(name => 'Transcript_Id')->[0];
		$transcript_datasource_id = $datasource2->{'datasource_id'} if (defined $datasource2);
	};
	throw('Wrong datasource') if ($@);

	# Do commit of everything
	$self->dbadaptor->commit();
	
	foreach my $aux_entity (@{$aux_entity_list})
	{
		# If the type of entity is not defined, we return the list of gene and transcript
		if (defined $type) {
			if ( ($type eq 'gene') and ($aux_entity->{'datasource_id'} eq $gene_datasource_id) ) {
				my ($entity) = $self->fetch_by_gene_stable_id($aux_entity->{'identifier'},$sources);
				push(@{$entity_list},$entity) if (defined $entity);
			}
			elsif ( ($type eq 'transcript') and ($aux_entity->{'datasource_id'} eq $transcript_datasource_id) ) {
				my ($entity) = $self->fetch_by_transc_stable_id($aux_entity->{'identifier'},$sources);
				push(@{$entity_list},$entity) if (defined $entity);			
			}			
		}
		else {
			if ($aux_entity->{'datasource_id'} eq $gene_datasource_id) {
				my ($entity) = $self->fetch_by_gene_stable_id($aux_entity->{'identifier'},$sources);
				push(@{$entity_list},$entity) if (defined $entity);
			}
			elsif ($aux_entity->{'datasource_id'} eq $transcript_datasource_id) {
				my ($entity) = $self->fetch_by_transc_stable_id($aux_entity->{'identifier'},$sources);
				push(@{$entity_list},$entity) if (defined $entity);			
			}			
		}
	}	

	return $entity_list;
}

=head2 fetch_basic_by_xref_entry

  Arg [1]    : String $xref_entry 
               An external identifier of the entity to be obtained
  Example    : $gene = $registry->fetch_basic_by_xref_entry('TRMT2A');
  Description: Retrieves a entity object (gene or transcript) 
               from the database via external reference entry.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Gene or APPRIS::Transcript or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_basic_by_xref_entry {
	my ($self, $xref_entry) = @_;

	my ($entity) = $self->fetch_entity_by_xref_entry($xref_entry);
	return $entity;
}

=head2 fetch_by_xref_entry

  Arg [1]    : String $xref_entry 
               An external identifier of the entity to be obtained
  Example    : $gene = $registry->fetch_by_xref_entry('TRMT2A');
  Description: Retrieves a entity object (gene or transcript) 
               from the database via external reference entry.
               If the gene or transcript is not found
               undef is returned instead.
               Note: Entity object contains the results of analysis methods.
  Returntype : APPRIS::Gene or APPRIS::Transcript or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_xref_entry {
	my ($self, $xref_entry, $sources) = @_;

	# Flag of 'get_analysis' is active. By default 'all'
	unless ( defined $sources ) { $sources = 'all' }
	
	my ($entity) = $self->fetch_entity_by_xref_entry($xref_entry, $sources);
	return $entity;
}

=head2 fetch_basic_by_stable_id

  Arg [1]    : name of the type of entity in the registry.
  Arg [2]    : String $id 
               The stable ID of the entity to retrieve
  Example    : $gene = $registry->fetch_basic_by_stable_id('gene','ENSG00000148944');
  Description: Retrieves a entity object (gene or transcript) 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Gene or APPRIS::Transcript or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_basic_by_stable_id {
	my ($self, $type, $stable_id) = @_;

	my ($entity);

	if ($type eq 'gene') {
		$entity = $self->fetch_by_gene_stable_id($stable_id);
	} elsif ($type eq 'transcript') {
		$entity = $self->fetch_by_transc_stable_id($stable_id);
	}
	return $entity;
}

=head2 fetch_by_stable_id

  Arg [1]    : name of the type of entity in the registry.
  Arg [2]    : String $id 
               The stable ID of the entity to retrieve
  Example    : $gene = $registry->fetch_by_stable_id('gene','ENSG00000148944');
  Description: Retrieves a entity object (gene or transcript) 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
               Note: Entity object contains the results of analysis methods.
  Returntype : APPRIS::Gene or APPRIS::Transcript or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
	my ($self, $type, $stable_id, $sources) = @_;
	my ($entity);
	
	# Flag of 'get_analysis' is active. By default 'all'
	unless ( defined $sources ) { $sources = 'all' }
			
	if ($type eq 'gene') {
		$entity = $self->fetch_by_gene_stable_id($stable_id, $sources); 
	} elsif ($type eq 'transcript') {
		$entity = $self->fetch_by_transc_stable_id($stable_id, $sources);
	}
	return $entity;
}

=head2 fetch_entity_by_xref_entry

  Arg [1]    : String $xref_entry 
               An external identifier of the entity to be obtained
  Arg [2]    : Int $analysis (optional)
               Flag that controls if the analysis are reported
  Example    : $gene = $registry->fetch_entity_by_xref_entry('TRMT2A');
  Description: Retrieves a gene object 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Gene or APPRIS::Transcript or undef
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_entity_by_xref_entry {
	my ($self, $xref_entry, $sources) = @_;
	
	my ($entity_list);
	
	# Connection to DBAdaptor
	unless ( defined $self->dbadaptor ) {
		my ($dbadaptor) = APPRIS::DBSQL::DBAdaptor->new
		(
			dbhost => $self->dbhost,
			dbport => $self->dbport,
			dbname => $self->dbname,
			dbuser => $self->dbuser,
			dbpass => $self->dbpass,
		);
		$self->dbadaptor($dbadaptor);	
	}

	# Get list of datasource
	my ($datasource_list);
	my ($datasource_list2);
	eval {
		$datasource_list = $self->dbadaptor->query_datasource();
		foreach my $ds (@{$datasource_list})
		{
			$datasource_list2->{$ds->{'datasource_id'}} = $ds;
		}
	};
	
	# Get entity report
	# - the query is not correct.
	# - the query has to be alphanumeric values.
	# - the query is too short (more than three characters).
	# - the query matched with the combination of theses values: ENS, ENST, ENSG, ENSP, CCDS, OTTHUM, OTTHUMG, OTTHUMT.
	# For this case, you have to give the exact identifier.	
	my ($aux_query_entity_list);
	my ($query_entity_list);
	# Extract "Gene ID" entity
	eval {
		$query_entity_list = $self->dbadaptor->query_entity(identifier => $xref_entry);
	};
	warning('Argument must be a correct entry') if ($@);
	# Extract "Gene Name" entity
	unless ( defined $query_entity_list and (scalar(@{$query_entity_list}) > 0) ) {
		eval {
			$aux_query_entity_list = $self->dbadaptor->query_xref_identify(identifier => $xref_entry);		
			foreach my $xref_report (@{$aux_query_entity_list})
			{
				my ($xref_identifier_list) = $self->dbadaptor->query_entity(entity_id => $xref_report->{'entity_id'});
				foreach my $xref_identifier (@{$xref_identifier_list})
				{
					push(@{$query_entity_list},{
										'entity_id'		=> $xref_identifier->{'entity_id'},
										'identifier'	=> $xref_identifier->{'identifier'},
										'datasource_id'	=> $xref_identifier->{'datasource_id'},
										'biotype'		=> $xref_identifier->{'biotype'},
										'status'		=> $xref_identifier->{'status'},
										'source'		=> $xref_identifier->{'source'},
										'level'			=> $xref_identifier->{'level'},
										'version'		=> $xref_identifier->{'version'},
										'tsl'			=> $xref_identifier->{'tsl'},
										'tag'			=> $xref_identifier->{'tag'}
					});
				}		
			}
		};		
		warning('Argument must be a correct entry') if ($@);
	}	
	unless ( defined $query_entity_list and (scalar(@{$query_entity_list}) > 0) ) {
		warning('Argument must be a correct entry');
	}
	
	# For each mapped entity it prints the list of objets
	foreach my $query_entity (@{$query_entity_list})
	{
		my ($entity_report);
		my ($entity_id) = $query_entity->{'entity_id'};
		my ($datasource_id) = $query_entity->{'datasource_id'};
		my ($namespace) = $datasource_list2->{$datasource_id}->{'name'};
		
		$entity_report->{'identifier'} = $query_entity->{'identifier'};		
		$entity_report->{'namespace'} = $namespace;
		$entity_report->{'biotype'} = $query_entity->{'biotype'};
		$entity_report->{'status'} = $query_entity->{'status'};
		$entity_report->{'source'} = $query_entity->{'source'};
		$entity_report->{'level'} = $query_entity->{'level'};
		$entity_report->{'version'} = $query_entity->{'version'};
		$entity_report->{'tsl'} = $query_entity->{'tsl'} if ( defined $query_entity->{'tsl'} );
		$entity_report->{'tag'} = $query_entity->{'tag'} if ( defined $query_entity->{'tag'} );

		# Get coordinate of entity
		my ($coordinate);
		eval {
			my ($coordinate_list) = $self->dbadaptor->query_coordinate(entity_id => $entity_id);
			if (defined $coordinate_list and defined $coordinate_list->[0]) {
				$coordinate = $coordinate_list->[0];
			}
			else {
				#throw('Wrong coordinate');
				warning('No coordinate');
			}
		};
		throw('Wrong coordinate') if ($@);
					
		# Get external name, xref identify
		my ($external_name);
		my ($xref_identifies);
		my ($trans_id_list);
		foreach my $datasource (@{$datasource_list})
		{
			eval {
				my ($datasource_name) = $datasource->{'name'};
				my ($datasource_id) = $datasource->{'datasource_id'};
				my ($datasource_desc) = $datasource->{'description'};
				my ($xref_id_list) = $self->dbadaptor->query_xref_identify
				(
					entity_id => $entity_id,
					datasource_id => $datasource_id
				);
				my ($xref_id);
				if (defined $xref_id_list and scalar(@{$xref_id_list}) > 0) {
					$xref_id = $xref_id_list->[0]->{'identifier'};
				}
				# External name			
				if ($datasource_name eq 'External_Id' and defined $xref_id) {
					$external_name = $xref_id;
				}
				# Xref identifiers
				elsif ((
						$datasource_name eq 'CCDS' or
						$datasource_name eq 'UniProtKB_SwissProt' or 
						$datasource_name eq 'Gene_Id' or
						$datasource_name eq 'Protein_Id'
						)
						and defined $xref_id) {
					push(@{$xref_identifies},
						APPRIS::XrefEntry->new
						(
							-id				=> $xref_id,
							-dbname			=> $datasource_name,
							-description	=> $datasource_desc
						)
					);
				}
				# Xref identifiers
				elsif (($datasource_name eq 'Transcript_Id') and defined $xref_id_list and (scalar(@{$xref_id_list}) > 0)) {
					foreach my $xref (@{$xref_id_list}) {
						push(@{$trans_id_list}, $xref->{'identifier'});
					}
				}							
			};
			throw('Wrong xref identify') if ($@);	
		}
		
		# Get transcript list
		my ($transcripts);
		my ($index_transcripts);
		my ($index) = 0;
		foreach my $transcript_id (@{$trans_id_list}) {
			my ($transcript) = $self->fetch_by_transc_stable_id($transcript_id, $sources);
			if (defined $transcript_id) {
				push(@{$transcripts}, $transcript);
				$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
			}
			else {
				throw('Wrong transcript');	
			}
		}
		
		# Create object
		my ($entity);
		if ($namespace eq 'Gene_Id')
		{
			$entity = APPRIS::Gene->new();
			$entity->transcripts($transcripts, $index_transcripts);
		}
		elsif ($namespace eq 'Transcript_Id')
		{
			$entity = APPRIS::Transcript->new();			
		}

		$entity->chromosome($coordinate->{'chromosome'});
		$entity->start($coordinate->{'start'});
		$entity->end($coordinate->{'end'});
		$entity->strand($coordinate->{'strand'});
		$entity->stable_id($entity_report->{'identifier'});
		$entity->biotype($entity_report->{'biotype'});
		$entity->status($entity_report->{'status'});
		$entity->source($entity_report->{'source'});
		$entity->level($entity_report->{'level'});
		$entity->version($entity_report->{'version'});
		$entity->tsl($entity_report->{'tsl'}) if (defined $entity_report->{'tsl'});
		$entity->tag($entity_report->{'tag'}) if (defined $entity_report->{'tag'});
		$entity->external_name($external_name) if (defined $external_name);
		$entity->xref_identify($xref_identifies) if (defined $xref_identifies);
		
		push(@{$entity_list}, $entity) if (defined $entity);
	}	

	# Do commit of everything
	$self->dbadaptor->commit();

	return $entity_list;
}

=head2 fetch_by_gene_stable_id

  Arg [1]    : String $id 
               The stable ID of the gene to retrieve
  Arg [2]    : Int $analysis (optional)
               Flag that controls if the analysis are reported
  Example    : $gene = $registry->fetch_by_gene_stable_id('ENSG00000148944');
  Description: Retrieves a gene object 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Gene or undef
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_gene_stable_id {
	my ($self, $stable_id, $sources, $region) = @_;
	
	# Connection to DBAdaptor
	unless ( defined $self->dbadaptor ) {
		my ($dbadaptor) = APPRIS::DBSQL::DBAdaptor->new
		(
			dbhost => $self->dbhost,
			dbport => $self->dbport,
			dbname => $self->dbname,
			dbuser => $self->dbuser,
			dbpass => $self->dbpass,
		);
		$self->dbadaptor($dbadaptor);	
	}

	# Get entity report
	my ($entity);
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $stable_id);
		if (defined $entity_list and scalar(@{$entity_list}) > 0) {
			$entity = $entity_list->[0];
		}
		else {
			warning('Argument must be a correct stable id');
		}
	};
	return undef if ( $@ or !(defined $entity) );

	# Get coordinate of entity
	my ($coordinate);
	eval {
		my ($coordinate_list) = $self->dbadaptor->query_coordinate(entity_id => $entity->{'entity_id'});
		if (defined $coordinate_list and defined $coordinate_list->[0]) {
			$coordinate = $coordinate_list->[0];
		}
		else {
			throw('Wrong coordinate');
		}
	};
	#throw('Wrong coordinate') if ($@);
	warning('Wrong coordinate') if ($@);

	# Get the list of datasources for this gene
	my ($datasource_list) = $self->dbadaptor->query_datasource();
	
	# Get external name, xref identifiers, and transcript id list
	my ($external_name);
	my ($xref_identifies);
	my ($trans_id_list);	
	foreach my $datasource (@{$datasource_list}) {
		eval {
			my ($datasource_name) = $datasource->{'name'};
			my ($datasource_id) = $datasource->{'datasource_id'};
			my ($datasource_desc) = $datasource->{'description'};
			
			my ($xref_id_list) = $self->dbadaptor->query_xref_identify
			(
				entity_id => $entity->{'entity_id'},
				datasource_id => $datasource_id
			);
			
			# External name
			if (($datasource_name eq 'External_Id') and defined $xref_id_list and (scalar(@{$xref_id_list}) > 0)) {
				my ($xref_id) = $xref_id_list->[0]->{'identifier'};
				$external_name = $xref_id if (defined $xref_id);
			}
			# Xref identifiers
			elsif (($datasource_name eq 'UniProtKB_SwissProt') and defined $xref_id_list and (scalar(@{$xref_id_list}) > 0)) {
				my ($xref_id) = $xref_id_list->[0]->{'identifier'};
				if (defined $xref_id) {
					push(@{$xref_identifies},
						APPRIS::XrefEntry->new
						(
							-id				=> $xref_id,
							-dbname			=> $datasource_name,
							-description	=> $datasource_desc
						)
					);					
				}
			}
			# Xref identifiers
			elsif (($datasource_name eq 'Transcript_Id') and defined $xref_id_list and (scalar(@{$xref_id_list}) > 0)) {
				foreach my $xref (@{$xref_id_list}) {
					push(@{$trans_id_list}, $xref->{'identifier'});
				}
			}			
		};
		throw('Wrong xref identify') if ($@);	
	}

	# Get transcript list
	my ($transcripts);
	my ($index_transcripts);
	my ($index) = 0;
	foreach my $transcript_id (@{$trans_id_list}) {
		my ($transcript) = $self->fetch_by_transc_stable_id($transcript_id, $sources, $region);
		if (defined $transcript_id) {
			push(@{$transcripts}, $transcript);
			$index_transcripts->{$transcript_id} = $index; $index++; # Index the list of transcripts
		}
		else {
			throw('Wrong transcript');	
		}
	}
	
	# Do commit of everything
	$self->dbadaptor->commit();

	# Create object
	my ($gene) = APPRIS::Gene->new
	(
		-chr		=> $coordinate->{'chromosome'},
		-start		=> $coordinate->{'start'},
		-end		=> $coordinate->{'end'},
		-strand		=> $coordinate->{'strand'},
		-stable_id	=> $entity->{'identifier'},
		-biotype	=> $entity->{'biotype'},
		-status		=> $entity->{'status'},
		-source		=> $entity->{'source'},
		-level		=> $entity->{'level'},
		-version	=> $entity->{'version'}
	);
	$gene->external_name($external_name) if (defined $external_name);
	$gene->xref_identify($xref_identifies) if (defined $xref_identifies);
	$gene->transcripts($transcripts, $index_transcripts);
	
	return $gene;
}

=head2 fetch_by_transc_stable_id

  Arg [1]    : String $id 
               The stable ID of the transcript to retrieve
  Arg [2]    : Int $analysis (optional)
               Flag that controls if the analysis are reported
  Example    : $transcript = $registry->fetch_trans_by_stable_id('ENST00000148944');
  Description: Retrieves a transcript object 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Transcript or undef
  Exceptions : if we cant get the transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_transc_stable_id {
	my ($self, $stable_id, $sources, $region) = @_;
	
	# Connection to DBAdaptor
	unless ( defined $self->dbadaptor ) {
		my ($dbadaptor) = APPRIS::DBSQL::DBAdaptor->new
		(
			dbhost => $self->dbhost,
			dbport => $self->dbport,
			dbname => $self->dbname,
			dbuser => $self->dbuser,
			dbpass => $self->dbpass,
		);
		$self->dbadaptor($dbadaptor);	
	}
	
	# Get gene and trascript list
	my ($gene_datasource_id);
	my ($transcript_datasource_id);
	eval {
		my ($datasource) = $self->dbadaptor->query_datasource(name => 'Gene_Id')->[0];
		$gene_datasource_id = $datasource->{'datasource_id'} if (defined $datasource);
		my ($datasource2) = $self->dbadaptor->query_datasource(name => 'Transcript_Id')->[0];
		$transcript_datasource_id = $datasource2->{'datasource_id'} if (defined $datasource2);
	};
	throw('Wrong datasource') if ($@);	

	# Get entity report
	my ($entity);
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $stable_id, datasource_id => $transcript_datasource_id);
		if (defined $entity_list and scalar(@{$entity_list}) > 0) {
			$entity = $entity_list->[0];
		}
		else {
			warning('Argument must be a correct stable id');
		}
	};
	return undef if ( $@ or !(defined $entity) );
	
	# Get coordinate of entity
	my ($coordinate);
	eval {
		my ($coordinate_list) = $self->dbadaptor->query_coordinate(entity_id => $entity->{'entity_id'});
		if (defined $coordinate_list and defined $coordinate_list->[0]) {
			$coordinate = $coordinate_list->[0];
		}
		else {
			#throw('Wrong coordinate');
			warning('No coordinate');
		}
	};
	throw('Wrong coordinate') if ($@);

	# Get sequence
	my($sequence);
	eval {
		my ($type_list) = $self->dbadaptor->query_type(name => 'transcript');
		my ($transcript_type_id) = $type_list->[0]->{'type_id'};			
		
		my($sequence_list)=$self->dbadaptor->query_sequence
		(
			entity_id => $entity->{'entity_id'},
			type_id => $transcript_type_id
		);
		if(defined $sequence_list and scalar(@{$sequence_list}) > 0) {
			$sequence=$sequence_list->[0]->{'sequence'};		
		}
	};
	throw('Wrong nucleotide sequence') if ($@);	
	
	# Get the list of datasource for this transcript
	my ($datasource_list) = $self->dbadaptor->query_datasource();

	# Get external name, xref identify
	my ($external_name);
	my ($xref_identifies);
	foreach my $datasource (@{$datasource_list}) {
		eval {
			my ($datasource_name) = $datasource->{'name'};
			my ($datasource_id) = $datasource->{'datasource_id'};
			my ($datasource_desc) = $datasource->{'description'};
			my ($xref_id_list) = $self->dbadaptor->query_xref_identify
			(
				entity_id => $entity->{'entity_id'},
				datasource_id => $datasource_id
			);
			my ($xref_id);
			if (defined $xref_id_list and scalar(@{$xref_id_list}) > 0) {
				$xref_id = $xref_id_list->[0]->{'identifier'};
			}

			# External name			
			if ($datasource_name eq 'External_Id' and defined $xref_id) {
				$external_name = $xref_id;
			}
			# Xref identifiers
			elsif ((
					$datasource_name eq 'CCDS' or
					$datasource_name eq 'UniProtKB_SwissProt' or 
					$datasource_name eq 'Gene_Id' or
					$datasource_name eq 'Protein_Id'
					)
					and defined $xref_id) {
				push(@{$xref_identifies},
					APPRIS::XrefEntry->new
					(
						-id				=> $xref_id,
						-dbname			=> $datasource_name,
						-description	=> $datasource_desc
					)
				);
			}
		};
		throw('Wrong xref identify') if ($@);	
	}

	# Get translation
	my ($translate) = $self->fetch_by_transl_stable_id($stable_id);

	# Get exons
	my ($exons);
	eval {
		my ($exon_list) = $self->dbadaptor->query_exon(entity_id => $entity->{'entity_id'}); 
		foreach my $exon (@{$exon_list}){
			my($exon_report);
			$exon_report->{'exon_id'}=$exon->{'exon_id'} if(defined $exon->{'exon_id'});
			$exon_report->{'identifier'}=$exon->{'identifier'} if(defined $exon->{'identifier'});						
			$exon_report->{'strand'}=$exon->{'strand'} if(defined $exon->{'strand'});
			$exon_report->{'start'}=$exon->{'start'} if(defined $exon->{'start'});
			$exon_report->{'end'}=$exon->{'end'} if(defined $exon->{'end'});						
			
			push(@{$exons},
				APPRIS::Exon->new
				(
					-start		=> $exon->{'start'},
					-end		=> $exon->{'end'},
					-strand		=> $exon->{'strand'},
					-stable_id	=> $exon->{'identifier'}
				)
			);
		}
	};
	throw('Wrong exons') if ($@);	

	# Get Analysis
	my ($analysis);
	if (defined $sources) {
		$analysis = $self->fetch_analysis_by_stable_id($stable_id, $sources, $region);
	}
	
	# Do commit of everything
	$self->dbadaptor->commit();

	# Create object
	my ($transcript) = APPRIS::Transcript->new
	(
		-chr		=> $coordinate->{'chromosome'},
		-start		=> $coordinate->{'start'},
		-end		=> $coordinate->{'end'},
		-strand		=> $coordinate->{'strand'},
		-stable_id	=> $entity->{'identifier'},
		-biotype	=> $entity->{'biotype'},
		-status		=> $entity->{'status'},
		-source		=> $entity->{'source'},
		-level		=> $entity->{'level'},
		-version	=> $entity->{'version'},
		-tsl		=> $entity->{'tsl'},
		-tag		=> $entity->{'tag'},
	);
	$transcript->sequence($sequence) if (defined $sequence);
	$transcript->external_name($external_name) if (defined $external_name);
	$transcript->xref_identify($xref_identifies) if (defined $xref_identifies);
	$transcript->translate($translate) if (defined $translate);
	$transcript->exons($exons) if (defined $exons);
	$transcript->analysis($analysis) if (defined $analysis);

	return $transcript;	
}

=head2 fetch_by_transl_stable_id

  Arg [1]    : String $id 
               The stable ID of the transcript to retrieve
  Example    : $translate = $registry->fetch_by_transl_stable_id('ENST000001');
  Description: Retrieves a transcript object 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Translation or undef
  Exceptions : if we cant get the transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_transl_stable_id {
	my ($self, $stable_id) = @_;
	
	# Connection to DBAdaptor
	unless ( defined $self->dbadaptor ) {
		my ($dbadaptor) = APPRIS::DBSQL::DBAdaptor->new
		(
			dbhost => $self->dbhost,
			dbport => $self->dbport,
			dbname => $self->dbname,
			dbuser => $self->dbuser,
			dbpass => $self->dbpass,
		);
		$self->dbadaptor($dbadaptor);	
	}
	
	# Get gene and trascript list
	my ($gene_datasource_id);
	my ($transcript_datasource_id);
	eval {
		my ($datasource) = $self->dbadaptor->query_datasource(name => 'Gene_Id')->[0];
		$gene_datasource_id = $datasource->{'datasource_id'} if (defined $datasource);
		my ($datasource2) = $self->dbadaptor->query_datasource(name => 'Transcript_Id')->[0];
		$transcript_datasource_id = $datasource2->{'datasource_id'} if (defined $datasource2);
	};
	throw('Wrong datasource') if ($@);	

	# Get entity report
	my ($entity);
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $stable_id, datasource_id => $transcript_datasource_id);
		if (defined $entity_list and scalar(@{$entity_list}) > 0) {
			$entity = $entity_list->[0];
		}
		else {
			throw('Argument must be a correct stable id');
		}
	};
	throw('Argument must be a correct stable id') if ($@);
	
	# Get sequence
	my ($sequence);
	eval {
		my ($type_list) = $self->dbadaptor->query_type(name => 'peptide');
		my ($transcript_type_id) = $type_list->[0]->{'type_id'};			
		
		my($sequence_list)=$self->dbadaptor->query_sequence
		(
			entity_id => $entity->{'entity_id'},
			type_id => $transcript_type_id
		);
		if(defined $sequence_list and scalar(@{$sequence_list}) > 0) {
			$sequence=$sequence_list->[0]->{'sequence'};		
		}
	};
	throw('Wrong aminoacid sequence') if ($@);

	# Get translate id
	my ($translate_id);
	eval {
		my ($datasource) = $self->dbadaptor->query_datasource(name => 'Protein_Id')->[0];
		if (defined $datasource) {
			my ($datasource_name) = $datasource->{'name'};
			my ($datasource_id) = $datasource->{'datasource_id'};
			my ($datasource_desc) = $datasource->{'description'};
			my ($xref_id_list) = $self->dbadaptor->query_xref_identify
			(
				entity_id => $entity->{'entity_id'},
				datasource_id => $datasource_id
			);
			if (defined $xref_id_list and scalar(@{$xref_id_list}) > 0) {
				$translate_id = $xref_id_list->[0]->{'identifier'};
			}	
		}
	};
	throw('Wrong translate identify') if ($@);	

	# Get CDS
	my ($cds);
	eval {
		my ($cds_list) = $self->dbadaptor->query_cds(entity_id => $entity->{'entity_id'});
		foreach my $cds_rep (@{$cds_list}){
			push(@{$cds},
				APPRIS::CDS->new
				(
					-start		=> $cds_rep->{'start'},
					-end		=> $cds_rep->{'end'},
					-strand		=> $cds_rep->{'strand'},
					-phase		=> $cds_rep->{'phase'}
				)
			);
		}
	};
	throw('Wrong cds') if ($@);	
	
	# Get Codon
	my ($codons);
	eval {
		my ($codon_list) = $self->dbadaptor->query_codon(entity_id => $entity->{'entity_id'});  		
		foreach my $codon_rep (@{$codon_list}){
			push(@{$codons},
				APPRIS::Codon->new
				(
					-start		=> $codon_rep->{'start'},
					-end		=> $codon_rep->{'end'},
					-strand		=> $codon_rep->{'strand'},
					-phase		=> $codon_rep->{'phase'},
					-type		=> $codon_rep->{'type'}
				)
			);
		}
	};
	throw('Wrong codon') if ($@);	
	
	# Do commit of everything
	$self->dbadaptor->commit();

	my ($translate);
	if (defined $sequence) {
		$translate = APPRIS::Translation->new();
		$translate->stable_id($translate_id) if (defined $translate_id);	
		$translate->sequence($sequence) if (defined $sequence);	
		$translate->cds($cds) if (defined $cds);
		$translate->codons($codons) if (defined $codons);
	}
		
	return $translate;
}

=head2 fetch_analysis_by_xref_entry

  Arg [1]    : String $xref_entry 
               An external identifier of the entity to be obtained
  Arg [2]    : String $type (optional)
               The type of analysis of the transcript to retrieve
  Example    : $analysis = $registry->fetch_analysis_by_xref_entry('CRISP3-001','firestar');
  Description: Retrieves an analysis object 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Analysis or undef
  Exceptions : if we cant get the transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_analysis_by_xref_entry {
	my ($self, $xref_entry, $sources) = @_;
	my ($entity);
	my ($xref_list) = $self->fetch_basic_by_xref_entry($xref_entry);
	if (defined $xref_list and scalar(@{$xref_list}) > 0) {
		my ($xref) = $xref_list->[0];
		if ( $xref->stable_id ) {
			$entity = $self->fetch_analysis_by_stable_id($xref->stable_id, $sources);
		}
	}
	
	return $entity;	
}

=head2 fetch_analysis_by_stable_id

  Arg [1]    : String $id 
               The stable ID of the transcript to retrieve
  Arg [2]    : String $sources (optional)
               The type of analysis of the transcript to retrieve
  Example    : $analysis = $registry->fetch_analysis_by_stable_id('ENST000001','firestar,corsair');
  Description: Retrieves an analysis object 
               from the database via its stable id.
               If the gene or transcript is not found
               undef is returned instead.
  Returntype : APPRIS::Analysis or undef
  Exceptions : if we cant get the transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_analysis_by_stable_id {
	my ($self, $stable_id, $sources, $region) = @_;
	my ($analysis);
	my ($firestar);
	my ($matador3d);
	my ($spade);
	my ($inertia);
	my ($crash);
	my ($thump);
	my ($corsair);
	my ($proteo);
	my ($appris);
		
	# Connection to DBAdaptor
	unless ( defined $self->dbadaptor ) {
		my ($dbadaptor) = APPRIS::DBSQL::DBAdaptor->new
		(
			dbhost => $self->dbhost,
			dbport => $self->dbport,
			dbname => $self->dbname,
			dbuser => $self->dbuser,
			dbpass => $self->dbpass,
		);
		$self->dbadaptor($dbadaptor);	
	}
	
	# Get gene and trascript list
	my ($gene_datasource_id);
	my ($transcript_datasource_id);
	eval {
		my ($datasource) = $self->dbadaptor->query_datasource(name => 'Gene_Id')->[0];
		$gene_datasource_id = $datasource->{'datasource_id'} if (defined $datasource);
		my ($datasource2) = $self->dbadaptor->query_datasource(name => 'Transcript_Id')->[0];
		$transcript_datasource_id = $datasource2->{'datasource_id'} if (defined $datasource2);
	};
	throw('Wrong datasource') if ($@);	

	# Get entity report
	my ($entity);
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $stable_id, datasource_id => $transcript_datasource_id);
		if (defined $entity_list and scalar(@{$entity_list}) > 0) {
			$entity = $entity_list->[0];
		}
		else {
			throw('Argument must be a correct stable id');
		}
	};
	throw('Argument must be a correct stable id') if ($@);

	foreach my $source ( split(',', $sources) ) {
		
		# Get FIRESTAR analysis -----------------
		if (defined $source and ($source eq 'firestar' or $source eq 'all')) {
			my ($method);
			my ($regions);		
			eval {	
				my ($list) = $self->dbadaptor->query_firestar(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No firestar analysis') if ($@);
			eval {
				my ($list);
				if (defined $region) {
					$list = $self->dbadaptor->query_firestar_residues(entity_id => $entity->{'entity_id'}, region => {start => $region->{'start'}, end => $region->{'end'}})
				}
				else {
					$list = $self->dbadaptor->query_firestar_residues(entity_id => $entity->{'entity_id'})
				}
				foreach my $residue (@{$list}) {
					push(@{$regions},
						APPRIS::Analysis::FirestarRegion->new
						(
							-start		=> $residue->{'trans_start'},
							-end		=> $residue->{'trans_end'},
							-strand		=> $residue->{'trans_strand'},
							-residue	=> $residue->{'peptide_position'},
							-domain		=> $residue->{'domain'},
							-ligands	=> $residue->{'ligands'},
						)				
					);
				}
			};
			throw('No firestar analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$firestar = APPRIS::Analysis::Firestar->new
					(
						-result					=> $method->{'result'},
						-functional_residue		=> $method->{'functional_residue'},
						-num_residues			=> $method->{'num_residues'}
					);
					$firestar->residues($regions) if (defined $regions);				
				};
				throw('No firestar object') if ($@);
			}	
		}
		
		# Get MATADOR3D analysis -----------------
		if (defined $source and ($source eq 'matador3d' or $source eq 'all')) {
			my ($method);
			my ($regions);
			my ($num_alignments) = 0; # TODO: add this field into database
			eval {	
				my ($list) = $self->dbadaptor->query_matador3d(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No matador analysis') if ($@);
			eval {	
				my ($list);
				if (defined $region) {
					$list = $self->dbadaptor->query_matador3d_alignments(entity_id => $entity->{'entity_id'}, region => {start => $region->{'start'}, end => $region->{'end'}})
				}
				else {
					$list = $self->dbadaptor->query_matador3d_alignments(entity_id => $entity->{'entity_id'})
				}				
				foreach my $residue (@{$list}) {
					$num_alignments++ if ($residue->{'alignment_score'} ne '0');
					push(@{$regions},
						APPRIS::Analysis::Matador3DRegion->new
						(
							-start	=> $residue->{'trans_start'},
							-end	=> $residue->{'trans_end'},
							-strand	=> $residue->{'trans_strand'},
							-cds_id	=> $residue->{'cds_id'},
							-pstart	=> $residue->{'start'},
							-pend	=> $residue->{'end'},
							-type	=> $residue->{'type'},
							-alignment_start	=> $residue->{'alignment_start'},
							-alignment_end	=> $residue->{'alignment_end'},
							-score	=> $residue->{'alignment_score'},
							-pdb_id	=> $residue->{'pdb_id'},
							-identity	=> $residue->{'identity'},
							-external_id	=> $residue->{'external_id'}
						)				
					);
				}
			};
			throw('No matador analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$matador3d = APPRIS::Analysis::Matador3D->new
					(
						-result						=> $method->{'result'},
						-conservation_structure		=> $method->{'conservation_structure'},
						-score						=> $method->{'score'},
						-num_alignments				=> $num_alignments
					);
					$matador3d->alignments($regions) if (defined $regions);				
				};
				throw('No matador object') if ($@);
			}
		}
		
		# Get SPADE analysis -----------------
		if (defined $source and ($source eq 'spade' or $source eq 'all')) {
			my ($method);
			my ($regions);		
			eval {	
				my ($list) = $self->dbadaptor->query_spade(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No spade analysis') if ($@);
			eval {
				my ($list);
				if (defined $region) {
					$list = $self->dbadaptor->query_spade_alignments(entity_id => $entity->{'entity_id'}, region => {start => $region->{'start'}, end => $region->{'end'}})
				}
				else {
					$list = $self->dbadaptor->query_spade_alignments(entity_id => $entity->{'entity_id'})
				}				
				foreach my $residue (@{$list}) {
					push(@{$regions},
						APPRIS::Analysis::SPADERegion->new
						(
							-start								=> $residue->{'trans_start'},
							-end								=> $residue->{'trans_end'},
							-strand								=> $residue->{'trans_strand'},						
							-score								=> $residue->{'score'},						
							-type_domain						=> $residue->{'type_domain'},						
							-alignment_start					=> $residue->{'alignment_start'},					
							-alignment_end						=> $residue->{'alignment_end'},
							-envelope_start						=> $residue->{'envelope_start'},					
							-envelope_end						=> $residue->{'envelope_end'},
							-hmm_start							=> $residue->{'hmm_start'},					
							-hmm_end							=> $residue->{'hmm_end'},
							-hmm_length							=> $residue->{'hmm_length'},					
							-hmm_acc							=> $residue->{'hmm_acc'},
							-hmm_name							=> $residue->{'hmm_name'},					
							-hmm_type							=> $residue->{'hmm_type'},
							-bit_score							=> $residue->{'bit_score'},
							-evalue								=> $residue->{'evalue'},
							-significance						=> $residue->{'significance'},
							-clan								=> $residue->{'clan'},
							-predicted_active_site_residues		=> $residue->{'predicted_active_site_residues'},
							-external_id						=> $residue->{'external_id'}						
						)				
					);
				}
			};
			throw('No spade analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$spade = APPRIS::Analysis::SPADE->new
					(
						-result							=> $method->{'result'},
						-domain_signal					=> $method->{'domain_signal'},
						-num_domains					=> $method->{'num_domains'},
						-num_possibly_damaged_domains	=> $method->{'num_possibly_damaged_domains'},
						-num_damaged_domains			=> $method->{'num_damaged_domains'},
						-num_wrong_domains				=> $method->{'num_wrong_domains'}
					);
					$spade->regions($regions) if (defined $regions);				
				};
				throw('No spade object') if ($@);
			}
		}
		
		# Get INERTIA analysis -----------------
		if (defined $source and ($source eq 'inertia' or $source eq 'all')) {
	
			my ($method);
			my ($method2);
			my ($method3);
			my ($method4);
			my ($regions);
			my ($regions2);		
			my ($regions3);
			my ($regions4);
			eval {	
				my ($list) = $self->dbadaptor->query_inertia(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];				
				}
			};
			throw('No inertia analysis') if ($@);
			eval {	
				my ($list2) = $self->dbadaptor->query_inertia_residues(entity_id => $entity->{'entity_id'});
				foreach my $residue (@{$list2}) {
					push(@{$regions},
						APPRIS::Analysis::INERTIARegion->new
						(
							-start				=> $residue->{'trans_start'},
							-end				=> $residue->{'trans_end'},
							-strand				=> $residue->{'trans_strand'},
							-unusual_evolution	=> $residue->{'inertia_residue_unusual_evolution'}
						)
					);
				}
			};
			throw('No inertia analysis') if ($@);	
			
			eval { # Get mafft alignment 			
				my ($type_list) = $self->dbadaptor->query_slr_type(name => 'filter');
				my ($type_id) = $type_list->[0]->{'slr_type_id'};
	
				my ($list) = $self->dbadaptor->query_omega(
									entity_id => $entity->{'entity_id'},
									slr_type_id => $type_id
				);
				if (defined $list and scalar(@{$list}) > 0) {
					$method2 = $list->[0];				
				}
				my ($list2) = $self->dbadaptor->query_omega_residues(
									entity_id => $entity->{'entity_id'},
									slr_type_id => $type_id
				);
				foreach my $residue (@{$list2}) {
					push(@{$regions2},
						APPRIS::Analysis::OmegaRegion->new
						(
							-start				=> $residue->{'trans_start'},
							-end				=> $residue->{'trans_end'},
							-omega_mean			=> $residue->{'omega_mean'},
							-st_deviation		=> $residue->{'st_deviation'},
							-p_value			=> $residue->{'p_value'},
							-difference_value	=> $residue->{'difference_value'},
							-unusual_evolution	=> $residue->{'unusual_evolution'}
						)
					);
				}
			};
			throw('No inertia-omega-slr-mafft analysis') if ($@);				
	
			eval { # Get prank alignment		
				my ($type_list) = $self->dbadaptor->query_slr_type(name => 'prank');
				my ($type_id) = $type_list->[0]->{'slr_type_id'};
	
				my ($list) = $self->dbadaptor->query_omega(
									entity_id	=> $entity->{'entity_id'},
									slr_type_id	=> $type_id
				);
				if (defined $list and scalar(@{$list}) > 0) {
					$method3 = $list->[0];				
				}
				my ($list2) = $self->dbadaptor->query_omega_residues(
									entity_id	=> $entity->{'entity_id'},
									slr_type_id	=> $type_id
				);
				foreach my $residue (@{$list2}) {
					push(@{$regions3},
						APPRIS::Analysis::OmegaRegion->new
						(
							-start				=> $residue->{'trans_start'},
							-end				=> $residue->{'trans_end'},
							-omega_mean			=> $residue->{'omega_mean'},
							-st_deviation		=> $residue->{'st_deviation'},
							-p_value			=> $residue->{'p_value'},
							-difference_value	=> $residue->{'difference_value'},
							-unusual_evolution	=> $residue->{'unusual_evolution'}
						)
					);
				}
			};
			throw('No inertia-omega-slr-prank analysis') if ($@);	
	
			eval { # Get kalign alignment		
				my ($type_list) = $self->dbadaptor->query_slr_type(name => 'kalign');
				my ($type_id) = $type_list->[0]->{'slr_type_id'};
	
				my ($list) = $self->dbadaptor->query_omega(
									entity_id	=> $entity->{'entity_id'},
									slr_type_id	=> $type_id
				);
				if (defined $list and scalar(@{$list}) > 0) {
					$method4 = $list->[0];				
				}
				my ($list2) = $self->dbadaptor->query_omega_residues(
									entity_id	=> $entity->{'entity_id'},
									slr_type_id	=> $type_id
				);
				foreach my $residue (@{$list2}) {
					push(@{$regions4},
						APPRIS::Analysis::OmegaRegion->new
						(
							-start				=> $residue->{'trans_start'},
							-end				=> $residue->{'trans_end'},
							-omega_mean			=> $residue->{'omega_mean'},
							-st_deviation		=> $residue->{'st_deviation'},
							-p_value			=> $residue->{'p_value'},
							-difference_value	=> $residue->{'difference_value'},
							-unusual_evolution	=> $residue->{'unusual_evolution'}
						)
					);
				}
			};
			throw('No inertia-omega-slr-kalign analysis') if ($@);	
	
			if (defined $method) { # Create object
				eval {
					$inertia = APPRIS::Analysis::INERTIA->new
					(
						-result					=> $method->{'result'},
						-unusual_evolution		=> $method->{'unusual_evolution'}
					);
					$inertia->regions($regions) if (defined $regions);				
				};
				throw('No inertia object') if ($@);
				
				if (defined $method2) { # Create object
					eval {		
						my ($mafft) = APPRIS::Analysis::Omega->new
						(
							-average			=> $method2->{'omega_average'},
							-st_desviation		=> $method2->{'omega_st_desviation'},
							-result				=> $method2->{'result'},
							-unusual_evolution	=> $method2->{'unusual_evolution'}
						);
						$mafft->regions($regions2) if (defined $regions2);
						$inertia->mafft_alignment($mafft);
					};
					throw('No omega-maf object') if ($@);
				}
		
				if (defined $method3) { # Create object
					eval {		
						my ($prank) = APPRIS::Analysis::Omega->new
						(
							-average			=> $method3->{'omega_average'},
							-st_desviation		=> $method3->{'omega_st_desviation'},
							-result				=> $method3->{'result'},
							-unusual_evolution	=> $method3->{'unusual_evolution'}
						);
						$prank->regions($regions3) if (defined $regions3);
						$inertia->prank_alignment($prank);
					};
					throw('No omega-prank object') if ($@);
				}
		
				if (defined $method4) { # Create object
					eval {				
						my ($kalign) = APPRIS::Analysis::Omega->new
						(
							-average			=> $method4->{'omega_average'},
							-st_desviation		=> $method4->{'omega_st_desviation'},
							-result				=> $method4->{'result'},
							-unusual_evolution	=> $method4->{'unusual_evolution'}
						);
						$kalign->regions($regions4) if (defined $regions4);
						$inertia->kalign_alignment($kalign);
					};
					throw('No inertia object') if ($@);
				}
			}
	
		}
	
		# Get CRASH analysis -----------------
		if (defined $source and ($source eq 'crash' or $source eq 'all')) {
			my ($method);
			my ($regions);		
			eval {	
				my ($list) = $self->dbadaptor->query_crash(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No crash analysis') if ($@);
			eval {
				my ($list);
				if (defined $region) {
					$list = $self->dbadaptor->query_crash_residues(entity_id => $entity->{'entity_id'}, region => {start => $region->{'start'}, end => $region->{'end'}})
				}
				else {
					$list = $self->dbadaptor->query_crash_residues(entity_id => $entity->{'entity_id'})
				}			
				foreach my $residue (@{$list}) {
					push(@{$regions},
						APPRIS::Analysis::CRASHRegion->new
						(
							-start						=> $residue->{'trans_start'},
							-end						=> $residue->{'trans_end'},
							-strand						=> $residue->{'trans_strand'},						
							-pstart						=> $residue->{'start'},
							-pend						=> $residue->{'end'},
							-s_mean						=> $residue->{'s_mean'},						
							-s_prob						=> $residue->{'s_prob'},						
							-d_score					=> $residue->{'d_score'},					
							-c_max						=> $residue->{'c_max'},
							-reliability				=> $residue->{'reliability'},					
							-localization				=> $residue->{'localization'},
						)				
					);
				}
			};
			throw('No crash analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$crash = APPRIS::Analysis::CRASH->new
					(					
						-result						=> $method->{'result'},
						-sp_score					=> $method->{'sp_score'},
						-tp_score					=> $method->{'tp_score'},
						-peptide_signal				=> $method->{'peptide_signal'},
						-mitochondrial_signal		=> $method->{'mitochondrial_signal'}
					);
					$crash->regions($regions) if (defined $regions);				
				};
				throw('No crash object') if ($@);
			}
		}
	
		# Get THUMP analysis -----------------
		if (defined $source and ($source eq 'thump' or $source eq 'all')) {
			my ($method);
			my ($regions);		
			eval {	
				my ($list) = $self->dbadaptor->query_thump(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No thump analysis') if ($@);
			eval {
				my ($list);
				if (defined $region) {
					$list = $self->dbadaptor->query_thump_helixes(entity_id => $entity->{'entity_id'}, region => {start => $region->{'start'}, end => $region->{'end'}})
				}
				else {
					$list = $self->dbadaptor->query_thump_helixes(entity_id => $entity->{'entity_id'})
				}
				foreach my $residue (@{$list}) {
					push(@{$regions},
						APPRIS::Analysis::THUMPRegion->new
						(
							-start			=> $residue->{'trans_start'},
							-end			=> $residue->{'trans_end'},
							-strand			=> $residue->{'trans_strand'},						
							-pstart			=> $residue->{'start'},						
							-pend			=> $residue->{'end'},						
							-damaged		=> $residue->{'damaged'}					
						)				
					);
				}
			};
			throw('No thump analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$thump = APPRIS::Analysis::THUMP->new
					(
						-result						=> $method->{'result'},
						-transmembrane_signal		=> $method->{'transmembrane_signal'},
						-num_tmh					=> $method->{'num_tmh'},
						-num_damaged_tmh			=> $method->{'num_damaged_tmh'}
					);
					$thump->regions($regions) if (defined $regions);				
				};
				throw('No thump object') if ($@);
			}
		}
		
		# Get CORSAIR analysis -----------------
		if (defined $source and ($source eq 'corsair' or $source eq 'all')) {
			my ($method);
			my ($regions);		
			eval {	
				my ($list) = $self->dbadaptor->query_corsair(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No corsair analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$corsair = APPRIS::Analysis::CORSAIR->new
					(
						-vertebrate_signal		=> $method->{'vertebrate_signal'},
						-score					=> $method->{'score'},
						-result					=> $method->{'result'}
					);
				};
				throw('No corsair object') if ($@);
			}
		}	

		# Get PROTEO analysis -----------------
		if (defined $source and ($source eq 'proteo' or $source eq 'all')) {
			my ($method);
			my ($regions);		
			eval {	
				my ($list) = $self->dbadaptor->query_proteo(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No proteo analysis') if ($@);
			eval {	
				my ($list);
				if (defined $region) {
					$list = $self->dbadaptor->query_proteo_peptides(entity_id => $entity->{'entity_id'}, region => {start => $region->{'start'}, end => $region->{'end'}})
				}
				else {
					$list = $self->dbadaptor->query_proteo_peptides(entity_id => $entity->{'entity_id'})
				}				
				foreach my $residue (@{$list}) {
					push(@{$regions},
						APPRIS::Analysis::PROTEORegion->new
						(
							-start				=> $residue->{'trans_start'},
							-end				=> $residue->{'trans_end'},
							-strand				=> $residue->{'trans_strand'},
							-pstart				=> $residue->{'start'},						
							-pend				=> $residue->{'end'},							
							-peptide_id			=> $residue->{'peptide_id'},
							-sequence			=> $residue->{'sequence'},
							-num_experiments	=> $residue->{'num_experiments'},
							-experiments		=> $residue->{'experiments'},
						)				
					);
				}
			};
			throw('No proteo analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$proteo = APPRIS::Analysis::PROTEO->new
					(
						-result					=> $method->{'result'},
						-peptide_evidence		=> $method->{'peptide_evidence'},
						-num_peptides			=> $method->{'num_peptides'}
					);
					$proteo->peptides($regions) if (defined $regions);				
				};
				throw('No proteo object') if ($@);
			}	
		}
			
		# Get APPRIS analysis -----------------
		if (defined $source and ($source eq 'appris' or $source eq 'all')) {
			my ($method);
			my ($regions);		
			eval {	
				my ($list) = $self->dbadaptor->query_appris(entity_id => $entity->{'entity_id'});
				if (defined $list and scalar(@{$list}) > 0) {
					$method = $list->[0];
				}				
			};
			throw('No appris analysis') if ($@);
			if (defined $method) { # Create object
				eval {
					$appris = APPRIS::Analysis::APPRIS->new
					(
						-functional_residues_score			=> $method->{'functional_residues_score'},
						-homologous_structure_score			=> $method->{'homologous_structure_score'},
						-vertebrate_conservation_score		=> $method->{'vertebrate_conservation_score'},
						-domain_score						=> $method->{'domain_score'},
						-transmembrane_helices_score		=> $method->{'transmembrane_helices_score'},
						-peptide_score						=> $method->{'peptide_score'},
						-mitochondrial_score				=> $method->{'mitochondrial_score'},
						-unusual_evolution_score			=> $method->{'unusual_evolution_score'},
						-peptide_evidence_score				=> $method->{'peptide_evidence_score'},
						-principal_isoform_score			=> $method->{'principal_isoform_score'},
						-functional_residues_signal			=> $method->{'functional_residues_signal'},
						-homologous_structure_signal		=> $method->{'homologous_structure_signal'},
						-vertebrate_conservation_signal		=> $method->{'vertebrate_conservation_signal'},
						-domain_signal						=> $method->{'domain_signal'},
						-transmembrane_helices_signal		=> $method->{'transmembrane_helices_signal'},
						-peptide_signal						=> $method->{'peptide_signal'},
						-mitochondrial_signal				=> $method->{'mitochondrial_signal'},
						-unusual_evolution_signal			=> $method->{'unusual_evolution_signal'},
						-peptide_evidence_signal			=> $method->{'peptide_evidence_signal'},
						-principal_isoform_signal			=> $method->{'principal_isoform_signal'},
						-reliability						=> $method->{'reliability'},
						-result								=> $method->{'result'}						
					);
				};
				throw('No appris object') if ($@);
			}
		}	
		
		# Do commit of everything
		$self->dbadaptor->commit();
				
		# Create object
		if (defined $source) {
			$analysis = APPRIS::Analysis->new();
			if (defined $firestar) {
				$analysis->firestar($firestar);
				$analysis->number($analysis->number+1);
			}
			if (defined $matador3d) {
				$analysis->matador3d($matador3d);
				$analysis->number($analysis->number+1);
			}
			if (defined $spade) {
				$analysis->spade($spade);
				$analysis->number($analysis->number+1);
			}
			if (defined $inertia) {
				$analysis->inertia($inertia);
				$analysis->number($analysis->number+1);
			}
			if (defined $crash) {
				$analysis->crash($crash);
				$analysis->number($analysis->number+1);
			}		
			if (defined $thump) {
				$analysis->thump($thump);
				$analysis->number($analysis->number+1);
			}
			if (defined $corsair) {
				$analysis->corsair($corsair);
				$analysis->number($analysis->number+1);
			}
			if (defined $proteo) {
				$analysis->proteo($proteo);
				$analysis->number($analysis->number+1);
			}
			if (defined $appris) {
				$analysis->appris($appris);
				$analysis->number($analysis->number+1);
			}
		}
	}
	
	return $analysis;
}

sub DESTROY {
	my ($self) = @_;
	
	$self->dbadaptor->disconnect() if ( defined $self->dbadaptor );
	foreach my $attrname ($self->_standard_keys) {
		$self->{$attrname} = $self->_default_for($attrname);
	}
}

1;
