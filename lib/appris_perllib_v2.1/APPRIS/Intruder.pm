=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Intruder

=head1 SYNOPSIS

  $intruder = APPRIS::Intruder->new(
    -dbhost  => 'localhost',
    -dbname  => 'homo_sapiens_encode_3c',
    -dbuser  => 'jmrodriguez'
    );

  $intruder->feed_by_analysis_from_stable_id($gene);

=head1 DESCRIPTION

All Adaptors are stored/registered using this module.
This module should then be used to get the adaptors needed.

The intruder can be loaded from a configuration file or from database info.


=head1 METHODS

=cut

package APPRIS::Intruder;

use strict;
use warnings;
use Data::Dumper;
use Config::IniFiles;

use APPRIS::Gene;
use APPRIS::Transcript;
use APPRIS::Translation;
use APPRIS::Exon;
use APPRIS::CDS;
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
use APPRIS::Utils::Exception qw(throw warning deprecate info);
use APPRIS::Utils::Constant qw(
        $API_VERSION
        $VERSION
        $DATE
);

my ($API_VERSION) = $APPRIS::Utils::Constant::API_VERSION;
my ($VERSION) = $APPRIS::Utils::Constant::VERSION;
my ($DATE) = $APPRIS::Utils::Constant::DATE;

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

    APPRIS::Intruder->new()

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
  
sub software_version {
	my ($self) = @_;
	return $API_VERSION;
}

=head2 date
  
  get the date of exported data.
  
  Args       : none
  ReturnType : string
  Status     : At Risk
  
=cut
  
sub date {
	my ($self) = @_;
	return $DATE;
}

=head2 version
  
  get the version of exported data.
  
  Args       : none
  ReturnType : string
  Status     : At Risk
  
=cut
  
sub version {
	my ($self) = @_;
	return $VERSION;
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

    $intruder->load_registry_from_db(
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

    $intruder->load_registry_from_ini(
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

=head2 _fetch_datasources

  Example    : $intruder->_fetch_datasources();
  Description: Retrieves a listref of all datasources from database.
  Returntype : Listref of hash
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub _fetch_datasources {
	my ($self) = @_;
	my ($ids);
	
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
		
	eval {
		my ($list) = $self->dbadaptor->query_datasource();
		foreach my $data (@{$list}) {
			my ($name) = $data->{'name'};
			my ($id) = $data->{'datasource_id'};
			$ids->{$name} = $id;
		}
	};
	throw('Getting datasource identifiers') if ($@);
	
	# Do commit of everything
	$self->dbadaptor->commit();

	return $ids;
}

=head2 _fetch_type_seq

  Example    : $intruder->_fetch_type_seq();
  Description: Retrieves a listref of all type of sequences.
  Returntype : Listref of hash
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub _fetch_type_seq {
	my ($self) = @_;
	my ($ids);
	
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
		
	eval {
		my ($list) = $self->dbadaptor->query_type();
		foreach my $data (@{$list}) {
			my ($name) = $data->{'name'};
			my ($id) = $data->{'type_id'};
			$ids->{$name} = $id;
		}
	};
	throw('Getting sequence types') if ($@);
	
	# Do commit of everything
	$self->dbadaptor->commit();

	return $ids;
}

=head2 feed_by_gencode

  Arg [1]    : Listref of APPRIS::Gene or undef 
               The listred of gene objects
  Example    : $intruder->feed_by_gencode($genes);
  Description: Add a listref of gene objects into the database.
  Returntype : None
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub feed_by_gencode {
	my ($self, $genes) = @_;

	# Get the Identifiers for entity datasources
	my ($datasources) = $self->_fetch_datasources();

	# Get the types of sequences
	my ($types) = $self->_fetch_type_seq();

	# Scan genes
	if ( $genes and (ref($genes) eq 'ARRAY') ) {
		foreach my $gene (@{$genes}) {
			my ($gene_entity_id) = $self->feed_gene_by_gencode($gene, $datasources, $types);
		}
	}
	else {
		if ( $genes->isa("APPRIS::Gene") ) {
			my ($gene) = $genes;
			my ($gene_entity_id) = $self->feed_gene_by_gencode($gene, $datasources, $types);
		}
	}
	
	return undef;
}

=head2 feed_gene_by_gencode

  Arg [1]    : APPRIS::Gene $entity
  Arg [2]    : Listref of datasources (optional)
  Arg [3]    : Listref of sequence types (optional)
  Example    : $intruder->feed_gene_by_gencode($entity);
  Description: Add a listref of gene objects into the database.
  Returntype : Int - Entity id
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub feed_gene_by_gencode {
	my ($self, $entity, $datasources, $types) = @_;
	my ($entity_id);
	my ($stable_id) = $entity->stable_id;

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
	
	# Get datasource ids
	unless (defined $datasources) {
		$datasources = $self->_fetch_datasources();
	}

	# Get the types of sequences
	unless (defined $types) {
		$types = $self->_fetch_type_seq();
	}
	
	# Check if exits => ERROR
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $stable_id);
		throw('Gene stable id already exists') if(defined $entity_list and scalar(@{$entity_list})>0);
	};
	throw('Gene stable id already exists') if ($@);

	# Insert Gene Entity info
	eval {
		my (%parameters) = (
								datasource_id	=> $datasources->{'Gene_Id'},
								identifier		=> $stable_id,
								source			=> $entity->source,
								biotype			=> $entity->biotype,
								status			=> $entity->status,
								level			=> $entity->level,
								version			=> $entity->version
		);
		$entity_id = $self->dbadaptor->insert_entity(%parameters);
		throw('Inserting gene entity') unless ( defined $entity_id );			
	};
	throw('Inserting gene entity') if ($@);

	# Insert Gene Coordinate
	if ( $entity->chromosome and $entity->start and $entity->end and $entity->strand ) {
		eval {
			my (%parameters) = (
									entity_id		=> $entity_id,
									chromosome		=> $entity->chromosome,
									start			=> $entity->start,
									end				=> $entity->end,
									strand			=> $entity->strand
			);
			$self->dbadaptor->insert_coordinate(%parameters);
		};
		throw('Inserting gene coordinates') if ($@);		
	}
		
	# Insert Xrefs
	eval {
		if ( $entity->xref_identify ) {
			foreach my $xref (@{$entity->xref_identify}) {
				my (%parameters) = (
								datasource_id	=> $datasources->{$xref->dbname},
								entity_id		=> $entity_id,
								identifier		=> $xref->id
				);
				$self->dbadaptor->insert_xref_identify(%parameters);					
			}
		}
	};
	throw('Inserting gene xrefs') if ($@);
		
	# Insert Transc entity
	if ( $entity->transcripts ) {
		foreach my $transcript (@{$entity->transcripts})
		{
			# insert xref transc id
			eval {			
				my (%parameters) = (
								datasource_id	=> $datasources->{'Transcript_Id'},
								entity_id		=> $entity_id,
								identifier		=> $transcript->stable_id
				);
				$self->dbadaptor->insert_xref_identify(%parameters);
			};
			throw('Inserting gene xrefs (transc)') if ($@);
			
			my ($transc_entity_id) = $self->feed_transc_by_gencode($transcript, $datasources, $types);
		}			
	}
	
	# Do commit of everything
	$self->dbadaptor->commit();

	return $entity_id;
}

=head2 feed_transc_by_gencode

  Arg [1]    : APPRIS::Transcript $entity
  Arg [2]    : Listref of datasources (optional)
  Arg [3]    : Listref of sequence types (optional)
  Example    : $intruder->feed_transc_by_gencode($entity);
  Description: Add a listref of gene objects into the database.
  Returntype : Int - Entity id
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub feed_transc_by_gencode {
	my ($self, $entity, $datasources, $types) = @_;
	my ($entity_id);
	my ($stable_id) = $entity->stable_id;

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
	
	# Get datasource ids
	unless (defined $datasources) {
		$datasources = $self->_fetch_datasources();
	}
	
	# Get the types of sequences
	unless (defined $types) {
		$types = $self->_fetch_type_seq();
	}
	
	# Check if exits => ERROR
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $stable_id, datasource_id	=> $datasources->{'Transcript_Id'});
		throw('Transcript stable id already exists') if(defined $entity_list and scalar(@{$entity_list})>0);
	};
	throw('Transcript stable id already exists') if ($@);

	# Insert Entity info
	eval {
		my (%parameters) = (
								datasource_id	=> $datasources->{'Transcript_Id'},
								identifier		=> $stable_id,
								source			=> $entity->source,
								biotype			=> $entity->biotype,
								status			=> $entity->status,
								level			=> $entity->level,
								version			=> $entity->version,
								tsl				=> $entity->tsl,
								tag				=> $entity->tag
		);
		$entity_id = $self->dbadaptor->insert_entity(%parameters);
		throw('Inserting transcript entity') unless ( defined $entity_id );			
	};
	throw('Inserting transcript entity') if ($@);
	
	# Insert Coordinate
	if ( $entity->chromosome and $entity->start and $entity->end and $entity->strand ) {
		eval {
			my (%parameters) = (
									entity_id		=> $entity_id,
									chromosome		=> $entity->chromosome,
									start			=> $entity->start,
									end				=> $entity->end,
									strand			=> $entity->strand
			);
			$self->dbadaptor->insert_coordinate(%parameters);
		};
		throw('Inserting transcript coordinates') if ($@);
	}
		
	# Insert Xrefs
	eval {
		if ( $entity->xref_identify ) {
			foreach my $xref (@{$entity->xref_identify}) {
				my (%parameters) = (
								datasource_id	=> $datasources->{$xref->dbname},
								entity_id		=> $entity_id,
								identifier		=> $xref->id
				);
				$self->dbadaptor->insert_xref_identify(%parameters);					
			}
		}
	};
	throw('Inserting transcript xrefs') if ($@);
	
	# Insert Exons
	eval {
		if ( $entity->exons ) {
			for (my $i = 0; $i < scalar(@{$entity->exons}); $i++) {
				my ($exon) = $entity->exons->[$i];
				my (%parameters) = (
								entity_id		=> $entity_id,
								exon_id			=> $i+1,
								identifier		=> $exon->stable_id,
								start			=> $exon->start,
								end				=> $exon->end,
								strand			=> $exon->strand
				);
				$self->dbadaptor->insert_exon(%parameters);
			}			
		}
	};
	throw('Inserting transcript exons') if ($@);

	# Insert Sequence
	eval {
		if ( $entity->sequence ) {
			my (%parameters) = (
								entity_id		=> $entity_id,
								type_id			=> $types->{'transcript'},
								length			=> length($entity->sequence),
								sequence		=> $entity->sequence
			);
			$self->dbadaptor->insert_sequence(%parameters);
		}
	};
	throw('Inserting transcript sequence') if ($@);
	
	# Insert Translation
	if ( $entity->translate ) {
		my ($transl_entity_id) = $self->_feed_transl_by_gencode($entity_id, $entity->translate, $datasources, $types);
	}
	
	# Do commit of everything
	$self->dbadaptor->commit();

	return $entity_id;
}

=head2 _feed_transl_by_gencode

  Arg [1]    : DBAdaptor
  Arg [2]    : Int - Transcript entity identifier $entity_id
  Arg [3]    : APPRIS::Translation $entity
  Arg [4]    : Listref of datasources
  Arg [5]    : Listref of sequence types
  Example    : $intruder->_feed_transl_by_gencode($entity);
  Description: Add a listref of gene objects into the database.
  Returntype : Int - Entity id
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub _feed_transl_by_gencode {
	my ($self, $entity_id, $entity, $datasources, $types) = @_;
	my ($stable_id) = $entity->stable_id;

	# Insert Entity id
	eval {
		if ( $stable_id =~ /^ENSP/ ) {
			my (%parameters) = (
								datasource_id	=> $datasources->{'Protein_Id'},
								entity_id		=> $entity_id,
								identifier		=> $stable_id
				);
			$self->dbadaptor->insert_xref_identify(%parameters);			
		}
	};
	throw('Inserting translation id') if ($@);

	# Insert CDS
	eval {
		if ( $entity->cds ) {
			for (my $i = 0; $i < scalar(@{$entity->cds}); $i++) {
				my ($cds) = $entity->cds->[$i];
				my (%parameters) = (
								entity_id		=> $entity_id,
								cds_id			=> $i+1,
								start			=> $cds->start,
								end				=> $cds->end,
								strand			=> $cds->strand,
								phase			=> $cds->phase
				);
				$self->dbadaptor->insert_cds(%parameters);
			}			
		}
	};
	throw('Inserting translation cds') if ($@);
	
	# Insert Condons
	eval {
		if ( $entity->codons ) {
			foreach my $codon (@{$entity->codons}) {
				my (%parameters) = (
								entity_id		=> $entity_id,
								type			=> $codon->type,
								start			=> $codon->start,
								end				=> $codon->end,
								strand			=> $codon->strand,
								phase			=> $codon->phase
				);
				$self->dbadaptor->insert_codon(%parameters);
			}			
		}
	};
	throw('Inserting translation codons') if ($@);
	
	# Insert Sequence
	eval {
		if ( $entity->sequence ) {
			my (%parameters) = (
								entity_id		=> $entity_id,
								type_id			=> $types->{'peptide'},
								length			=> length($entity->sequence),
								sequence		=> $entity->sequence
			);
			$self->dbadaptor->insert_sequence(%parameters);
		}
	};
	throw('Inserting translation sequence') if ($@);
	
		
	return $entity_id;
}

=head2 feed_by_analysis

  Arg [1]    : APPRIS::Gene $gene 
               The analysis of the result to insert
  Arg [2]    : String $type (optional)
               The type of analysis of the transcript to insert
  Example    : $intruder->feed_by_analysis($gene,'inertia');
  Description: Add an analysis object into the database.
  Returntype : None
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub feed_by_analysis {
	my ($self, $gene, $type) = @_;
	
	$self->feed_gene_by_analysis($gene, $type);
	
	return undef;
}

=head2 feed_gene_by_analysis

  Arg [1]    : APPRIS::Gene $entity 
               The analysis of the result to insert
  Arg [2]    : String $type (optional)
               The type of analysis of the transcript to insert
  Example    : $intruder->feed_by_analysis($gene,'inertia');
  Description: Add an analysis object into the database.
  Returntype : None
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub feed_gene_by_analysis {
	my ($self, $entity, $type) = @_;
	
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
	
	# Get datasource ids
	my ($datasources) = $self->_fetch_datasources();
	
	# Get Gene Entity info
	my ($entity_id) = $entity->stable_id;		
	my($internal_entity_id);
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $entity_id, datasource_id	=> $datasources->{'Gene_Id'});		
		unless (defined $entity_list and (scalar(@{$entity_list}) > 0) ) {
			throw('Argument must be a correct stable id');
		}			
		$internal_entity_id = $entity_list->[0]->{'entity_id'};			
	};
	throw('Argument must be a correct stable id') if ($@);
	throw('No a correct stable id') unless (defined $internal_entity_id);

	# If not defined 'type' of analysis, we retrieve all the results
	$type = 'all' unless (defined $type);

	if ($entity->analysis) {
		my ($analysis) = $entity->analysis;
			
		# Insert FIRESTAR analysis -----------------
		if (defined $type and ($type eq 'firestar' or $type eq 'all') and $analysis->firestar) {
			my ($method) = $analysis->firestar;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;					
			}
			else {
				throw('No firestar analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_firestar(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No firestar analysis') if ($@);
		}
		
		# Insert MATADOR3D analysis -----------------
		if (defined $type and ($type eq 'matador3d' or $type eq 'all') and $analysis->matador3d) {
			my ($method) = $analysis->matador3d;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;					
			}
			else {
				throw('No matador3d analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_matador3d(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No matador3d analysis') if ($@);
		}

		# Insert SPADE analysis -----------------
		if (defined $type and ($type eq 'spade' or $type eq 'all') and $analysis->spade) {
			my ($method) = $analysis->spade;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;					
			}
			else {
				throw('No spade analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_spade(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No spade analysis') if ($@);
		}
		
		# Insert INERTIA analysis -----------------
		if (defined $type and ($type eq 'inertia' or $type eq 'all') and $analysis->inertia) {
			my ($method) = $analysis->inertia;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;					
			}
			else {
				throw('No inertia analysis');			
			}
			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_inertia(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No inertia analysis') if ($@);

			# Insert omega
			my (@ALIGMENT_TYPES);
			eval {	
				my ($type_list) = $self->dbadaptor->query_slr_type();
				for (my $i = 0; $i < scalar(@{$type_list}); $i++) {
					push(@ALIGMENT_TYPES, $type_list->[$i]->{'name'});						
				}
			};
			throw('No alignment types') if ($@);
			
			foreach my $type (@ALIGMENT_TYPES) {				
				# get the type of alignment
				my ($alignment);
				if ( ($type eq 'filter') and $method->mafft_alignment) {
					$alignment = $method->mafft_alignment;
				}
				elsif ( ($type eq 'prank') and $method->prank_alignment) {
					$alignment = $method->prank_alignment;
				}
				elsif ( ($type eq 'kalign') and $method->kalign_alignment) {
					$alignment = $method->kalign_alignment;
				}
				else {
					next; # we don't have aligment
				}

				# insert main
				my ($type_id);
				eval {	
					my ($type_list) = $self->dbadaptor->query_slr_type(name => $type);
					$type_id = $type_list->[0]->{'slr_type_id'};
				};
				throw("No inertia-omega $type analysis") if ($@);
						
				if ( defined $alignment->result ) {
					eval {
						my ($global_id2) = $self->dbadaptor->insert_omega(
										inertia_id			=> $global_id,
										slr_type_id			=> $type_id,
										result				=> $alignment->result,
						);
					};
					throw("No inertia-omega $type analysis") if ($@);
				}
			}
			
		}
				
		# Insert CRASH analysis -----------------
		if (defined $type and ($type eq 'crash' or $type eq 'all') and $analysis->crash) {
			my ($method) = $analysis->crash;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;
			}
			else {
				throw('No crash analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_crash(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No crash analysis') if ($@);
		}

		# Insert THUMP analysis -----------------
		if (defined $type and ($type eq 'thump' or $type eq 'all') and $analysis->thump) {
			my ($method) = $analysis->thump;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;
			}
			else {
				throw('No thump analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_thump(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No thump analysis') if ($@);
		}
				
		# Insert CORSAIR analysis -----------------
		if (defined $type and ($type eq 'corsair' or $type eq 'all') and $analysis->corsair) {
			my ($method) = $analysis->corsair;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;
			}
			else {
				throw('No corsair analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_corsair(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No corsair analysis') if ($@);
		}

		# Insert PROTEO analysis -----------------
		if (defined $type and ($type eq 'proteo' or $type eq 'all') and $analysis->proteo) {
			my ($method) = $analysis->proteo;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;					
			}
			else {
				throw('No proteo analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_proteo(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No proteo analysis') if ($@);
		}
		
		# Insert APPRIS analysis -----------------
		if (defined $type and ($type eq 'appris' or $type eq 'all') and $analysis->appris) {
			my ($method) = $analysis->appris;
				
			# Insert method annotation
			my ($result);
			if ($method->result) {
				$result = $method->result;
			}
			else {
				throw('No appris analysis');			
			}

			my($global_id);
			eval {
				$global_id = $self->dbadaptor->insert_appris(
								entity_id		=> $internal_entity_id,
								result			=> $result
				);
			};
			throw('No appris analysis') if ($@);
		}
				
	}
	
	# Insert Method analysis by transcript
	if ($entity->transcripts) {
		foreach my $transc (@{$entity->transcripts}) {
			$self->feed_transc_by_analysis($transc, $type);
		}	
	}
	
	# Do commit of everything
	$self->dbadaptor->commit();

	return undef;	
}
		
=head2 feed_transc_by_analysis

  Arg [1]    : APPRIS::Transcript $entity
               The analysis of the result to insert
  Arg [2]    : String $type (optional)
               The type of analysis of the transcript to insert
  Example    : $intruder->feed_transc_by_analysis($entity,'inertia');
  Description: Add an analysis object into the database.
  Returntype : None
  Exceptions : if we find some problems
  Caller     : general
  Status     : Stable

=cut

sub feed_transc_by_analysis {
	my ($self, $entity, $type) = @_;
	
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
	
	# Get datasource ids
	my ($datasources) = $self->_fetch_datasources();	
		
	# Get Transcript Entity info
	my ($entity_id) = $entity->stable_id;
	my($internal_entity_id);
	eval {
		my ($entity_list) = $self->dbadaptor->query_entity(identifier => $entity_id, datasource_id	=> $datasources->{'Transcript_Id'});		
		unless (defined $entity_list and (scalar(@{$entity_list}) > 0) ) {
			throw('Argument must be a correct stable id');
		}			
		$internal_entity_id = $entity_list->[0]->{'entity_id'};			
	};
	throw('Argument must be a correct stable id') if ($@);
	throw('No a correct stable id') unless (defined $internal_entity_id);
	
	# If not defined 'type' of analysis, we retrieve all the results
	$type = 'all' unless (defined $type);

	if ($entity->analysis) {
		my ($analysis) = $entity->analysis;
		
		# Insert FIRESTAR analysis -----------------
		if (defined $type and ($type eq 'firestar' or $type eq 'all') and $analysis->firestar) {
			my ($method) = $analysis->firestar;

			# insert annotation
			my($global_id);
			if ( defined $method->result ) {
				eval {
					$global_id = $self->dbadaptor->insert_firestar(
									entity_id				=> $internal_entity_id,
									num_residues			=> $method->num_residues,
									result					=> $method->result
					);
				};
				throw('No firestar analysis') if ($@);
			}
			else {
				throw('No firestar analysis');			
			}

			# insert regions
			if ( defined $method->residues ) {
				foreach my $region (@{$method->residues}) {
					if ( defined $region->residue ) {
						eval {
							my (%parameters) = (
											firestar_id			=> $global_id,
											peptide_position	=> $region->residue
							);
							$parameters{trans_start} = $region->start if ( defined $region->start );
							$parameters{trans_end} = $region->end if ( defined $region->end );
							$parameters{trans_strand} = $region->strand if ( defined $region->strand );
							
							$parameters{domain} = $region->domain if ( $region->domain );
							$parameters{ligands} = $region->ligands if ( $region->ligands );
							my ($method_residues_id) = $self->dbadaptor->insert_firestar_residues(%parameters);
						};
						throw('No firestar residue analysis') if ($@);					
					}
					else {
						throw('No firestar residue analysis');						
					}
				}
			}
		}
			
		# Insert MATADOR3D analysis -----------------
		if (defined $type and ($type eq 'matador3d' or $type eq 'all') and $analysis->matador3d) {
			my ($method) = $analysis->matador3d;

			# insert annotation
			my($global_id);
			if ( defined $method->score and defined $method->result ) {
				eval {
					$global_id = $self->dbadaptor->insert_matador3d(
									entity_id				=> $internal_entity_id,
									score					=> $method->score,
									result					=> $method->result
					);
				};
				throw('No matador3d analysis') if ($@);
			}
			else {
				throw('No matador3d analysis');			
			}

			# insert regions
			if ( defined $method->alignments ) {
				foreach my $region (@{$method->alignments}) {
					if ( defined $region->cds_id and defined $region->pstart and defined $region->pend and defined $region->score ) {
						eval {
							my (%parameters) = (
											matador3d_id		=> $global_id,
											cds_id				=> $region->cds_id,
											start				=> $region->pstart,
											end					=> $region->pend,
											score				=> $region->score
							);
							$parameters{trans_start} = $region->start if ( defined $region->start );
							$parameters{trans_end} = $region->end if ( defined $region->end );
							$parameters{trans_strand} = $region->strand if ( defined $region->strand );							
							
							$parameters{type} = $region->type if ( $region->type );
							$parameters{alignment_start} = $region->alignment_start if ( $region->alignment_start );
							$parameters{alignment_end} = $region->alignment_end if ( $region->alignment_end );
							$parameters{pdb_id} = $region->pdb_id if ( $region->pdb_id );
							$parameters{identity} = $region->identity if ( $region->identity );
							$parameters{external_id} = $region->external_id if ( $region->external_id );
							my ($method_residues_id) = $self->dbadaptor->insert_matador3d_alignments(%parameters);
						};
						throw('No matador3d residue analysis') if ($@);					
					}
					else {
						throw('No matador3d residue analysis');						
					}
				}
			}
		}
				
		# Insert SPADE analysis -----------------
		if (defined $type and ($type eq 'spade' or $type eq 'all') and $analysis->spade) {
			my ($method) = $analysis->spade;

			# insert annotation
			my($global_id);
			if ( defined $method->result and 
				 defined $method->num_domains and defined $method->num_possibly_damaged_domains and
				 defined $method->num_damaged_domains and defined $method->num_wrong_domains ) {
				eval {
					$global_id = $self->dbadaptor->insert_spade(
									entity_id						=> $internal_entity_id,
									result							=> $method->result,
									num_domains						=> $method->num_domains,
									num_possibly_damaged_domains	=> $method->num_possibly_damaged_domains,
									num_damaged_domains				=> $method->num_damaged_domains,
									num_wrong_domains				=> $method->num_wrong_domains,
					);
				};
				throw('No spade analysis') if ($@);
			}
			else {
				throw('No spade analysis');			
			}

			# insert regions
			if ( defined $method->regions ) {
				foreach my $region (@{$method->regions}) {

					if ( defined $region->type_domain and 
						 defined $region->alignment_start and defined $region->alignment_end and
						 defined $region->envelope_start and defined $region->envelope_end and
						 defined $region->hmm_start and defined $region->hmm_end and defined $region->hmm_length and
						 defined $region->hmm_acc and defined $region->hmm_name and defined $region->hmm_type and
						 defined $region->bit_score and defined $region->evalue ) {						 	
						eval {
							my (%parameters) = (
											spade_id			=> $global_id,
											type_domain			=> $region->type_domain,
											alignment_start		=> $region->alignment_start,
											alignment_end		=> $region->alignment_end,
											envelope_start		=> $region->envelope_start,
											envelope_end		=> $region->envelope_end,
											hmm_start			=> $region->hmm_start,
											hmm_end				=> $region->hmm_end,
											hmm_length			=> $region->hmm_length,
											hmm_acc				=> $region->hmm_acc,
											hmm_name			=> $region->hmm_name,
											hmm_type			=> $region->hmm_type,
											bit_score			=> $region->bit_score,
											evalue				=> $region->evalue
							);
							$parameters{trans_start} = $region->start if ( defined $region->start );
							$parameters{trans_end} = $region->end if ( defined $region->end );
							$parameters{trans_strand} = $region->strand if ( defined $region->strand );
							
							$parameters{significance} = $region->significance if ( $region->significance );
							$parameters{clan} = $region->clan if ( $region->clan );
							$parameters{predicted_active_site_residues} = $region->predicted_active_site_residues if ( $region->predicted_active_site_residues );
							$parameters{external_id} = $region->external_id if ( $region->external_id );
							$parameters{discarded} = $region->discarded if ( defined $region->discarded );

							my ($method_residues_id) = $self->dbadaptor->insert_spade_alignments(%parameters);
						};
						throw('No spade residue analysis') if ($@);					
					}
					else {
						throw('No spade residue analysis');						
					}
				}
			}
		}

		# Insert INERTIA analysis -----------------
		if (defined $type and ($type eq 'inertia' or $type eq 'all') and $analysis->inertia) {
			my ($method) = $analysis->inertia;
			
			# insert annotation
			my($global_id);
			if ( defined  $method->unusual_evolution ) {
				eval {
					$global_id = $self->dbadaptor->insert_inertia(
									entity_id			=> $internal_entity_id,
									unusual_evolution	=> $method->unusual_evolution
					);
				};
				throw('No inertia analysis') if ($@);
			}
			else {
				throw('No inertia analysis');			
			}
			
			# insert region
			if ( defined $method->regions ) {
				for (my $icds = 0; $icds < scalar(@{$method->regions}); $icds++) {
					my ($region) = $method->regions->[$icds];
					if ( defined $region->start and defined $region->end and defined $region->strand and defined $region->unusual_evolution ) {
						eval {
							my ($method_residues_id) = $self->dbadaptor->insert_inertia_residues(
																	inertia_id			=> $global_id,
																	inertia_exon_id		=> $icds+1,
																	trans_start			=> $region->start,
																	trans_end			=> $region->end,
																	trans_strand		=> $region->strand,
																	unusual_evolution	=> $region->unusual_evolution
							);
						};
						throw('No inertia residue analysis') if ($@);					
					}
					else {
						throw('No inertia residue analysis');						
					}
				}
			}
	
			# Insert omega
			my (@ALIGMENT_TYPES);
			eval {	
				my ($type_list) = $self->dbadaptor->query_slr_type();
				for (my $i = 0; $i < scalar(@{$type_list}); $i++) {
					push(@ALIGMENT_TYPES, $type_list->[$i]->{'name'});						
				}
			};
			throw('No alignment types') if ($@);				
					
			foreach my $type (@ALIGMENT_TYPES) {
				
				# get the type of alignment
				my ($alignment);
				if ( ($type eq 'filter') and $method->mafft_alignment) {
					$alignment = $method->mafft_alignment;
				}
				elsif ( ($type eq 'prank') and $method->prank_alignment) {
					$alignment = $method->prank_alignment;
				}
				elsif ( ($type eq 'kalign') and $method->kalign_alignment) {
					$alignment = $method->kalign_alignment;
				}
				else {
					next; # we don't have aligment
				}

				# insert main
				my ($type_id2);
				eval {	
					my ($type_list) = $self->dbadaptor->query_slr_type(name => $type);
					$type_id2 = $type_list->[0]->{'slr_type_id'};
				};
				throw("No inertia-omega $type analysis") if ($@);
						
				if ( defined $alignment->result and defined $alignment->unusual_evolution ) {
					my($global_id2);
					eval {
						$global_id2 = $self->dbadaptor->insert_omega(
										inertia_id			=> $global_id,
										slr_type_id			=> $type_id2,
										omega_average		=> $alignment->average,
										omega_st_desviation	=> $alignment->st_desviation,
										result				=> $alignment->result,
										unusual_evolution	=> $alignment->unusual_evolution
						);
					};
					throw("No inertia-omega $type analysis") if ($@);
							
					# insert residues
					if ( defined $alignment->regions ) {
						for (my $icds = 0; $icds < scalar(@{$alignment->regions}); $icds++) {
							my ($region2) = $alignment->regions->[$icds];
							if ( defined $region2->start and defined $region2->end and defined $region2->unusual_evolution ) {
								eval {
									my ($method_residues_id) = $self->dbadaptor->insert_omega_residues(
																			omega_id			=> $global_id2,
																			omega_exon_id		=> $icds+1,
																			trans_start			=> $region2->start,
																			trans_end			=> $region2->end,
																			omega_mean			=> $region2->omega_mean,
																			st_deviation		=> $region2->st_deviation,
																			p_value				=> $region2->p_value,
																			difference_value	=> $region2->difference_value,
																			unusual_evolution	=> $region2->unusual_evolution
																			
									);
								};
								throw("No inertia-omega $type residues") if ($@);
							}
							else {
								throw("No inertia-omega $type residues");
								
							}
						}
					}			
				}
				else {
					throw("No inertia-omega $type analysis");
				}
			}
		}
		
		# Insert CRASH analysis -----------------
		if (defined $type and ($type eq 'crash' or $type eq 'all') and $analysis->crash) {
			my ($method) = $analysis->crash;

			# insert annotation
			my($global_id);
			if ( defined $method->result and 
				 defined $method->sp_score and defined $method->tp_score and
				 defined $method->peptide_signal and defined $method->mitochondrial_signal ) {
				eval {
					$global_id = $self->dbadaptor->insert_crash(
									entity_id					=> $internal_entity_id,
									result						=> $method->result,
									sp_score					=> $method->sp_score,
									tp_score					=> $method->tp_score,
									peptide_signal				=> $method->peptide_signal,
									mitochondrial_signal		=> $method->mitochondrial_signal,
					);
				};
				throw('No crash analysis') if ($@);
			}
			else {
				throw('No crash analysis');			
			}

			# insert regions
			if ( defined $method->regions ) {
				foreach my $region (@{$method->regions}) {

					if ( defined $region->s_mean and defined $region->s_prob and 
						 defined $region->d_score and defined $region->c_max and
						 defined $region->reliability and defined $region->localization and
						 defined $region->pstart and defined $region->pend ) {
						eval {
							my (%parameters) = (
											crash_id			=> $global_id,
											s_mean				=> $region->s_mean,
											s_prob				=> $region->s_prob,
											d_score				=> $region->d_score,
											c_max				=> $region->c_max,
											reliability			=> $region->reliability,
											localization		=> $region->localization,
											start				=> $region->pstart,
											end					=> $region->pend
							);
							$parameters{trans_start} = $region->start if ( defined $region->start );
							$parameters{trans_end} = $region->end if ( defined $region->end );
							$parameters{trans_strand} = $region->strand if ( defined $region->strand );
							
							my ($method_residues_id) = $self->dbadaptor->insert_crash_residues(%parameters);
						};
						throw('No crash residue analysis') if ($@);					
					}
					else {
						throw('No crash residue analysis');						
					}
				}
			}
		}
		
		# Insert THUMP analysis -----------------
		if (defined $type and ($type eq 'thump' or $type eq 'all') and $analysis->thump) {
			my ($method) = $analysis->thump;

			# insert annotation
			my($global_id);
			if ( defined $method->transmembrane_signal and defined $method->result and 
				 defined $method->num_tmh and defined $method->num_damaged_tmh ) {
				eval {
					$global_id = $self->dbadaptor->insert_thump(
									entity_id						=> $internal_entity_id,
									transmembrane_signal			=> $method->transmembrane_signal,
									result							=> $method->result,
									num_tmh							=> $method->num_tmh,
									num_damaged_tmh					=> $method->num_damaged_tmh
					);
				};
				throw('No thump analysis') if ($@);
			}
			else {
				throw('No thump analysis');			
			}

			# insert regions
			if ( defined $method->regions ) {
				foreach my $region (@{$method->regions}) {
					if ( defined $region->pstart and defined $region->pend ) {
						eval {
							my (%parameters) = (
											thump_id			=> $global_id,
											start				=> $region->pstart,
											end					=> $region->pend
							);
							$parameters{trans_start} = $region->start if ( defined $region->start );
							$parameters{trans_end} = $region->end if ( defined $region->end );
							$parameters{trans_strand} = $region->strand if ( defined $region->strand );							
							
							$parameters{damaged} = $region->damaged if ( $region->damaged );

							my ($method_residues_id) = $self->dbadaptor->insert_thump_helixes(%parameters);
						};
						throw('No thump residue analysis') if ($@);					
					}
					else {
						throw('No thump residue analysis');						
					}
				}
			}
		}
		
		# Insert CORSAIR analysis -----------------
		if (defined $type and ($type eq 'corsair' or $type eq 'all') and $analysis->corsair) {
			my ($method) = $analysis->corsair;

			# insert annotation
			my($global_id);
			if ( defined $method->vertebrate_signal and defined $method->result and defined $method->score ) {
				eval {
					$global_id = $self->dbadaptor->insert_corsair(
									entity_id						=> $internal_entity_id,
									vertebrate_signal				=> $method->vertebrate_signal,
									result							=> $method->result,
									score							=> $method->score,
					);
				};
				throw('No corsair analysis') if ($@);
			}
			else {
				throw('No corsair analysis');			
			}
			
			# insert regions
			if ( defined $method->alignments ) {
				foreach my $region (@{$method->alignments}) {
					if ( defined $region->cds_id and defined $region->pstart and defined $region->pend and defined $region->score ) {
						eval {
							my (%parameters) = (
											corsair_id			=> $global_id,
											cds_id				=> $region->cds_id,
											start				=> $region->pstart,
											end					=> $region->pend,
											score				=> $region->score
							);
							$parameters{trans_start} = $region->start if ( defined $region->start );
							$parameters{trans_end} = $region->end if ( defined $region->end );
							$parameters{trans_strand} = $region->strand if ( defined $region->strand );
							
							$parameters{type} = $region->type if ( $region->type );
							$parameters{maxscore} = $region->maxscore if ( $region->maxscore );
							$parameters{sp_report} = $region->sp_report if ( $region->sp_report );
							my ($method_residues_id) = $self->dbadaptor->insert_corsair_alignments(%parameters);
						};
						throw('No corsair residue analysis') if ($@);					
					}
					else {
						throw('No corsair residue analysis');						
					}
				}
			}
		}
		
		# Insert PROTEO analysis -----------------
		if (defined $type and ($type eq 'proteo' or $type eq 'all') and $analysis->proteo) {
			my ($method) = $analysis->proteo;

			# insert annotation
			my($global_id);
			if ( defined $method->peptide_evidence ) {
				eval {
					$global_id = $self->dbadaptor->insert_proteo(
									entity_id				=> $internal_entity_id,
									peptide_evidence		=> $method->peptide_evidence,
									num_peptides			=> $method->num_peptides,
									result					=> $method->result
					);
				};
				throw('No proteo analysis') if ($@);
			}
			else {
				throw('No proteo analysis');			
			}
			
			# insert regions
			if ( defined $method->peptides ) {
				foreach my $region (@{$method->peptides}) {
					if ( defined $region->sequence and defined $region->num_experiments and defined $region->pstart and defined $region->pend ) {
						eval {
							my (%parameters) = (
											proteo_id			=> $global_id,
											#peptide_id			=> $region->peptide_id,
											sequence			=> $region->sequence,
											num_experiments		=> $region->num_experiments,
											experiments			=> $region->experiments,
											start				=> $region->pstart,
											end					=> $region->pend
							);
							$parameters{trans_start} = $region->start if ( defined $region->start );
							$parameters{trans_end} = $region->end if ( defined $region->end );
							$parameters{trans_strand} = $region->strand if ( defined $region->strand );
							
							my ($method_peptides_id) = $self->dbadaptor->insert_proteo_peptides(%parameters);
						};
						throw('No proteo residue analysis') if ($@);					
					}
					else {
						throw('No proteo residue analysis');						
					}
				}
			}
		}
		
		# Insert APPRIS analysis -----------------
		if (defined $type and ($type eq 'appris' or $type eq 'all') and $analysis->appris) {
			my ($method) = $analysis->appris;

			# insert annotation
			my($global_id);
			if ( defined $method->principal_isoform_signal and defined $method->result ) {				
				eval {
					$global_id = $self->dbadaptor->insert_appris(
									entity_id							=> $internal_entity_id,									
									functional_residues_score			=> $method->functional_residues_score,
									homologous_structure_score			=> $method->homologous_structure_score,
									vertebrate_conservation_score		=> $method->vertebrate_conservation_score,
									domain_score						=> $method->domain_score,
									transmembrane_helices_score			=> $method->transmembrane_helices_score,
									peptide_score						=> $method->peptide_score,
									mitochondrial_score					=> $method->mitochondrial_score,
									unusual_evolution_score				=> $method->unusual_evolution_score,
									peptide_evidence_score				=> $method->peptide_evidence_score,
									principal_isoform_score				=> $method->principal_isoform_score,									
									functional_residues_signal			=> $method->functional_residues_signal,
									homologous_structure_signal			=> $method->homologous_structure_signal,						
									vertebrate_conservation_signal		=> $method->vertebrate_conservation_signal,
									domain_signal						=> $method->domain_signal,
									transmembrane_helices_signal		=> $method->transmembrane_helices_signal,
									peptide_signal						=> $method->peptide_signal,
									mitochondrial_signal				=> $method->mitochondrial_signal,
									unusual_evolution_signal			=> $method->unusual_evolution_signal,
									peptide_evidence_signal				=> $method->peptide_evidence_signal,
									principal_isoform_signal			=> $method->principal_isoform_signal,
									reliability							=> $method->reliability,
									result								=> $method->result,								
					);
				};
				throw('No appris analysis') if ($@);
			}
			else {
				throw('No appris analysis');			
			}
		}

		
		
	}
	else {
		info('No analysis: '.$entity_id);			
	}
	
	# Do commit of everything
	$self->dbadaptor->commit();

	return undef;	
}

sub DESTROY {
	my ($self) = @_;
	
	$self->dbadaptor->disconnect();
	foreach my $attrname ($self->_standard_keys) {
		$self->{$attrname} = $self->_default_for($attrname);
	}
}


1;
