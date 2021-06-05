=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

DBRetriever

=head1 DESCRIPTION

Package the retrieves the results of methods

=head1 SYNOPSIS

	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = # DEFAULT
		(
			conf			=>  undef,
			specie			=>  undef,
		);

	#_____________________________________________________________
    
  my $features = DBRetriever->new(
					-specie => 'mouse',
				    -conf   => 'conf/config.ini',
  );
  
  
=head1 METHODS

=cut

package DBRetriever;

use strict;
use warnings;
use FindBin;
use JSON;
use Config::IniFiles;
use File::Temp;
use Data::Dumper;

# these modules are loaded in load_registry function
use APPRIS::Registry;
use APPRIS::Exporter;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(info throw warning deprecate);
use APPRIS::Utils::File qw(getStringFromFile);

use APPRIS::Analysis;
use APPRIS::Parser qw(parse_trifid);
use APPRIS::Utils::TRIFID qw(get_trifid_pred_file_path);
#use WSRetriever;

{
	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = # DEFAULT
		(
			conf			=>  undef,
			species			=>  undef,
			dataset			=>  '',
			assembly		=>  undef,
			source			=>  undef,

			type			=>  undef,
			input			=>  undef,
			trifid_base_dir	=>  undef,
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
	
	# Substitute values within template
	sub dbconf {
		my ($self, $arg) = @_;
		$self->{'dbconf'} = $arg if defined $arg;
		return $self->{'dbconf'};
	}
	
	sub srvconf {
		my ($self, $arg) = @_;
		$self->{'srvconf'} = $arg if defined $arg;
		return $self->{'srvconf'};
	}

	sub species {
		my ($self, $arg) = @_;
		$self->{'species'} = $arg if defined $arg;
		return $self->{'species'};
	}
	
	sub assembly {
		my ($self, $arg) = @_;
		$self->{'assembly'} = $arg if defined $arg;
		return $self->{'assembly'};
	}
	
	sub source {
		my ($self, $arg) = @_;
		$self->{'source'} = $arg if defined $arg;
		return $self->{'source'};
	}

	sub dataset {
		my ($self, $arg) = @_;
		$self->{'dataset'} = $arg if defined $arg;
		return $self->{'dataset'};
	}

	sub type {
		my ($self, $arg) = @_;
		$self->{'type'} = $arg if defined $arg;
		return $self->{'type'};
	}
	
	sub input {
		my ($self, $arg) = @_;
		$self->{'input'} = $arg if defined $arg;
		return $self->{'input'};
	}

	sub trifid_base_dir {
		my ($self, $arg) = @_;
		$self->{'trifid_base_dir'} = $arg if defined $arg;
		return $self->{'trifid_base_dir'};
	}
}

=head2 new

  Arg [-dbconf]:
        string - config file of Database host
  Arg [-srvconf]:
        string - config file of Server
  Arg [-species]:
        string - specie name (human, mouse, rat)
  Arg [-assembly]:
        string - Assembly version
  Arg [-source]:
        string - source name (ensembl,refseq,uniprot)
  Arg [-dataset]:
        string - Dataset version
  Arg [-type]:
        string - type of input (id, name, position)
  Arg [-input]:
        string - input identifier/name/genomic region
  Example    : $features = DBRetriever->new(...);
  Description: Creates a new retrieve object
  Returntype : DBRetriever
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
	my ($caller) = shift;
	
	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;

	foreach my $attrname ($self->_standard_keys) {
		$self->{$attrname} = $self->_default_for($attrname);
	}

	my (
		$dbconf,
		$srvconf,
		$species,
		$assembly,
		$source,
		$dataset,
		$type,
		$input,
		$trifid_base_dir
	)
	= rearrange( [
		'dbconf',
		'srvconf',
		'species',
		'assembly',
		'source',
		'dataset',
		'type',
		'input',
		'trifid_base_dir'
	],
	@_
	);

	# require paramater
	if ( defined $dbconf and -e $dbconf ) { $self->dbconf($dbconf); }
	else { return undef; }
	if ( defined $srvconf and -e $srvconf ) { $self->srvconf($srvconf); }
	else { return undef; }

	# optional parameter
	if ( defined $species ) { $self->species($species); }
	if ( defined $assembly ) { $self->assembly($assembly); }
	if ( defined $source ) { $self->source($source); }
	if ( defined $dataset ) { $self->dataset($dataset); }
	if ( defined $type ) { $self->type($type); }
	if ( defined $input ) { $self->input($input); }
	if ( defined $trifid_base_dir ) { $self->trifid_base_dir($trifid_base_dir); }
	
	return $self;
}

=head2 load_registry

  Arg [1]    : String $type
               type of datase (current, archives)
  Example    : use DBRetriever qw(load_registry);
               my $rst = load_registry();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub load_registry {
	my ($self) = shift;
	my ($species)  = $self->species  || undef;
	my ($assembly) = $self->assembly || undef;
	my ($source)   = $self->source   || undef;
	my ($dataset)  = $self->dataset  || undef;
	my ($registries);
	
	# config ini
	my ($db_cfg) = new Config::IniFiles( -file => $self->dbconf );

	# Get Server config file
	my ($server_json) = JSON->new();
	my ($server_cfg) = $server_json->decode( getStringFromFile($self->srvconf) );
		
	# For each species. By default, all species
	my (@all_species);
	if ( defined $species ) { push(@all_species, $species) }
	else { @all_species = keys(%{$server_cfg->{'species'}}) }
	foreach my $species_id ( @all_species ) {
		my ($cfg_species) = $server_cfg->{'species'}->{$species_id};
		
		# For each assembly. By default, all assemblies
		foreach my $cfg_assembly (@{$cfg_species->{'assemblies'}}) {
			my ($found_as) = 0;
			my ($cfg_datasets) = [  # by default, scan all queryable datasets
				grep { ! exists($_->{'queryable'}) || $_->{'queryable'} }
				@{$cfg_assembly->{'datasets'}}
			];
			if ( defined $assembly ) {				
				if ( lc($assembly) eq lc($cfg_assembly->{'id'}) ) { $found_as = 1 }
				# NOTE!!! scan only one (the first, the newest queryable dataset)
				if ( !defined $dataset and !(defined $source) ) { $cfg_datasets = [$cfg_datasets->[0]]; }
			} else {
				#if ( lc($cfg_species->{'official'}) eq lc($cfg_assembly->{'id'}) ) { $found_as = 1 }
				$found_as = 1;
			}
			if ( $found_as == 1 ) {
				foreach my $cfg_dataset (@{$cfg_datasets}) {
					my ($found_sc,$found_ds) = (0,0);
					if ( defined $source ) {
						if ( lc($source) eq lc($cfg_dataset->{'source'}->{'name'}) ) { $found_sc = 1 }
					} else { $found_sc = 1 }

					if ( defined $dataset ) {
						if ( lc($cfg_dataset->{'id'}) =~ /^$dataset/ ) { $found_ds = 1 }				
					} else { $found_ds = 1 }

					if ( $found_sc == 1 and $found_ds == 1 and exists $cfg_dataset->{'database'} and exists $cfg_dataset->{'database'}->{'name'} ) {
						my ($db_host, $db_user, $db_pass, $db_port) = (undef, undef, undef, undef);
						my ($db_name) = $cfg_dataset->{'database'}->{'name'}.'_'.$cfg_dataset->{'id'};
						my ($trifid_release);
						if ( exists($cfg_dataset->{'trifid'}) && exists($cfg_dataset->{'trifid'}{'release'}) ) {
							$trifid_release = $cfg_dataset->{'trifid'}{'release'};
						}
						my ($registry) = APPRIS::Registry->new();
						if ( exists $cfg_dataset->{'type'} and lc($cfg_dataset->{'type'}) eq 'archive' ) {
							$db_host  = $db_cfg->val('APPRIS_ARCHIVES', 'host');
							$db_user  = $db_cfg->val('APPRIS_ARCHIVES', 'user');
							$db_pass  = $db_cfg->val('APPRIS_ARCHIVES', 'pass');
							$db_port  = $db_cfg->val('APPRIS_ARCHIVES', 'port');
						}
						else {
							$db_host  = $db_cfg->val('APPRIS_DATABASES', 'host');
							$db_user  = $db_cfg->val('APPRIS_DATABASES', 'user');
							$db_pass  = $db_cfg->val('APPRIS_DATABASES', 'pass');
							$db_port  = $db_cfg->val('APPRIS_DATABASES', 'port');							
						}
						$registry->load_registry_from_db(
												-dbhost	=> $db_host,
												-dbuser	=> $db_user,
												-dbpass	=> $db_pass,
												-dbport	=> $db_port,
												-dbname	=> $db_name
						);
						push(@{$registries}, {
								'species'	=> $species_id,
								'assembly'	=> $cfg_assembly->{'id'},
								'source'	=> $cfg_dataset->{'source'},
								'dataset'	=> $cfg_dataset->{'id'},
								'trifid_release' => $trifid_release,
								'type'		=> lc($cfg_dataset->{'type'}),
								'registry'	=> $registry								
						});
					}
				}
						
			}				
		}
	}
		
	return $registries;
	
} # end load_registry

=head2 get_features

  Arg [1]    : Hash $report
  Arg [2]    : String $method
               method list  
  Example    : use WSRetriever qw(get_features);
               my $rst = get_features($reports);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_features
{
	my ($self) = shift;
	my ($methods) = shift if (@_);
	my ($t_inputs) = shift if (@_);
	my ($type) = $self->type;
	my ($inputs) = $self->input;
	my ($species)  = $self->species;
	my ($assembly) = $self->assembly || undef;
	my ($source)   = $self->source   || undef;
	my ($dataset)  = $self->dataset  || undef;	
	my ($features);
	
	# We the list of methods will be able to contain the type of track: Eg. spade:domain|damaged_domain,firestar
	my ($methods_str) = '';
	foreach my $method (split(',', $methods)) {
		if ( $method =~ /^([^\-|\$]*)/ ) { $methods_str .= $1 . ',' }
	}
	$methods_str =~ s/\,$//;

	# Get the list of APPRIS::Registry
	my ($registries) = $self->load_registry();

	# For each registry, extract the input query
	REGISTRY: foreach my $registry (@{$registries}) {
		if ( defined $type ) {
			my (@registry_features);
			if ( ($type eq 'id') and $inputs ) {
				my ($query_type) = 'gene';
				if ( defined $t_inputs ) { $query_type = 'transcript'; $inputs = $t_inputs; }
				foreach my $input (split(/[;|]/, $inputs)) {
					my ($feat) = $self->get_feat_by_stable_id($registry, $query_type, $input, $methods_str);
					if ( defined $feat ) {
						push(@registry_features, $feat);
					}
				}
			}
			elsif ( ($type eq 'name') and $inputs ) {
				foreach my $input (split(/[;|]/, $inputs)) {
					my ($feat) = $self->get_feat_by_xref_entry($registry, $input, $methods_str);
					if ( defined $feat and scalar(@{$feat}) > 0 ) {
						foreach my $f (@{$feat}) { push(@registry_features, $f); }
					}
				}
			}
			elsif ( ($type eq 'position') and $inputs ) {
				foreach my $input (split(/[;|]/, $inputs)) {
					my ($feat) = $self->get_feat_by_region($registry, $input, $methods_str);
					if ( defined $feat and scalar(@{$feat}) > 0 ) {
						foreach my $f (@{$feat}) { push(@registry_features, $f); }
					}
				}
			}
			else {
				next REGISTRY;
			}

			# At the last possible moment, fetch TRIFID analyses for each gene, if available.
			if ( $methods =~ /trifid/ ) {
				foreach my $feature (@registry_features) {
					my ($trifid_release) = $registry->{'trifid_release'};
					if ( defined($trifid_release) ) {
						my ($pred_file_path) = get_trifid_pred_file_path($self->trifid_base_dir,
						                                                 $registry->{'species'},
						                                                 $trifid_release);
						my ($trifid_report) = _fetch_trifid_report($feature, $pred_file_path);
						if ( defined($trifid_report) ) {
							_add_trifid_analyses($feature, $trifid_report);
						}
					}
				}
			}

			push(@{$features}, @registry_features)
		}
	}

	return $features;
	
} # end get_features

=head2 get_seek_features

  Arg [1]    : Hash $report
  Arg [2]    : String $method
               method list  
  Example    : use WSRetriever qw(get_features);
               my $rst = get_features($reports);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_seek_features
{
	my ($self) = shift;
	my ($methods) = shift;
	my ($inputs) = $self->input;
	my ($species)  = $self->species  || undef;
	my ($assembly) = $self->assembly || undef;
	my ($source)   = $self->source   || undef;
	my ($dataset)  = $self->dataset  || undef;	
	my ($features);
	
	# Get the list of APPRIS::Registry
	my ($registries) = $self->load_registry();
	
	# For each registry, extract the input query
	foreach my $registry (@{$registries}) {
		foreach my $input (split(',', $inputs)) {
			my ($feat) = $self->get_feat_by_xref_entry($registry, $input, $methods);
			foreach my $f (@{$feat}) {
				push(@{$features}, {
					'species'	=> $registry->{'species'},
					'assembly'	=> $registry->{'assembly'},
					'source'	=> $registry->{'source'},
					'dataset'	=> $registry->{'dataset'},
					'type'		=> $registry->{'type'},
					'entity'	=> $f
				});
			}
		}
	}
				
	return $features;
	
} # end get_seek_features

sub get_feat_by_stable_id {
	my ($self) = shift;
	my ($aregistry) = shift;
	my ($type) = shift;
	my ($input) = shift;
	my ($methods) = shift;
	my ($registry) = $aregistry->{'registry'};
	my ($feat);
		
#	if (
#	(
#		lc($id) =~ /^ens(\w\w\w)?g([\d+|\.]+)$/ or
#		lc($id) =~ /^ens(\w\w\w)?gr([\d+|\.]+)$/ or
#		lc($id) =~ /^ensgr([\d+|\.]+)$/ or
#		lc($id) =~ /^fbgn([\d+|\.]+)$/ or
#		lc($id) =~ /^wbgene([\d+|\.]+)$/
#	)
#	and
#		($aregistry->{'source'} eq 'ensembl')
#	) {
#		$type = 'gene';
#	}
#	elsif (	(lc($id) =~ /^(\d+)$/) and ($aregistry->{'source'} eq 'refseq') ) {
#		$type = 'gene';
#	}
#	elsif (	(lc($id) !~ /-/) and ($aregistry->{'source'} eq 'uniprot') ) {
#		$type = 'gene';
#	}
#	elsif (
#	(
#		lc($id) =~ /^ens(\w\w\w)?t([\d+|\.]+)$/ or
#		lc($id) =~ /^ens(\w\w\w)?tr([\d+|\.]+)$/ or
#		lc($id) =~ /^enstr([\d+|\.]+)$/ or
#		lc($id) =~ /^fbtr([\d+|\.]+)$/
#	)
#	and
#		($aregistry->{'source'} eq 'ensembl')
#	) {
#		$type = 'transcript';
#	}
#	elsif (  (lc($id) =~ /^xm\_([\d+|\.]+)$/ or lc($id) =~ /^nm\_([\d+|\.]+)$/) and ($aregistry->{'source'} eq 'refseq') ) {
#		$type = 'transcript';
#	}
#	else {
#		$type = 'transcript';
#	}

	if ( defined $input and defined $type ) {
		$feat = $registry->fetch_by_stable_id($type, $input, $methods);
	}
	
	return $feat;	
}

sub get_feat_by_xref_entry {
	my ($self) = shift;
	my ($aregistry) = shift;
	my ($name) = shift;
	my ($methods) = shift;
	my ($registry) = $aregistry->{'registry'};
	
	my ($feat) = $registry->fetch_by_xref_entry($name, $methods);
	
	return $feat;	
}

sub get_feat_by_region {
	my ($self) = shift;
	my ($aregistry) = shift;
	my ($input) = shift;
	my ($methods) = shift;
	my ($chr,$start,$end) = (undef,undef,undef);
	my ($registry) = $aregistry->{'registry'};
	my ($feat);
	
	if ( $input =~ /^([^\:]*)\:(\d*)-(\d*)$/ ) {
		my ($chr) = $1;
		my ($start) = $2;
		my ($end) = $3;
		$chr = lc($chr); $chr =~ s/chr//;
		
		$feat = $registry->fetch_by_region($chr, $start, $end, $methods);
	}		
	elsif ( $input =~ /^chr([^\$]*)$/ ) {
		my ($chr) = $1;
		$chr = lc($chr); $chr =~ s/chr//;
		$feat = $registry->fetch_by_region($chr, $start, $end, $methods);
	}		
	
	return $feat;	
}

##
## EXPORTER
##

=head2 export_features

  Arg [1]    : String $type
               type of input: 'id' or 'name'
  Arg [2]    : String $input
               gene or transcript identifier (or name)
  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(export_features);
               my $rst = export_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub export_features
{
	my ($self) = shift;
	my ($format) = shift;
	my ($methods) = shift;
	my ($ids) = shift if (@_);
	my ($res) = shift if (@_);
	my ($result) = '';
	my ($features) = $self->get_features($methods, $ids); # if methods is not defined => retrieve all features

	if ( defined $features ) {
		my ($exporter) = APPRIS::Exporter->new();
		if ($format eq 'tsv') {
			$result = $exporter->get_tsv_annotations($features, $methods, $res);
		}
		elsif ( ($format eq 'bed') or ($format eq 'bed12') ) {
			$result = $exporter->get_bed_annotations($features, $methods, undef, $format, $self->source);
	    }
		elsif ($format eq 'json') {
			$result = $exporter->get_json_annotations($features, $methods);
	    }
		elsif ($format eq 'gtf') {
			$result = $exporter->get_gtf_annotations($features, $methods);
	    }
		elsif ($format eq 'gff3') {
			$result = $exporter->get_gff3_annotations($features, $methods);
	    }
		elsif ($format eq 'raw') {
			$result = $exporter->get_raw_annotations($features, $methods);
	    }
	}
	
	return $result;
	
} # end export_features


##
## SEQUENCER
##

=head2 export_seq_features

  Arg [1]    : String $type
               type of input: 'id' or 'name'
  Arg [2]    : String $input
               gene or transcript identifier (or name)
  Arg [3]    : String $operation
               method name
  Example    : use DBRetriever qw(export_seq_features);
               my $rst = export_seq_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub export_seq_features
{
	my ($self) = shift;
	my ($operation) = shift;
	my ($format) = shift;
	my ($type) = shift if (@_);
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);
	my ($res) = shift if (@_);
	my ($result) = '';	
	my ($features) = $self->get_features($methods, $ids); # if methods is not defined => retrieve all features
	
	if ( defined $features ) {
		my ($exporter) = APPRIS::Exporter->new();
		if ( $operation eq 'sequences' ) {
			$result = $exporter->get_seq_annotations($features, $type, $format);
		}
		elsif ( $operation eq 'align' ) { # "A minimum of 2 sequences is required\n"
			my ($seq_result) = $exporter->get_seq_annotations($features, $type, 'fasta');
			my ($num_seq) = $seq_result =~ tr/\>//;
			if ( $num_seq >= 2 ) {
				require WSRetriever;
				my ($wsretriever) = new WSRetriever();
				$result = $wsretriever->get_aln_annotations($seq_result, $format, $ids);
			}
			else {
				$result = $exporter->get_seq_annotations($features, $type, 'json');
			}			
		}
		elsif ($operation eq 'residues') {
			$result = $exporter->get_res_annotations($features, $methods, $res);
	    }
		elsif ($operation eq 'cds') {
			$result = $exporter->get_cds_annotations($features, $methods, $res);
	    }	    
		elsif ($operation eq 'genome') {
			$result = $self->get_gen_annotations($self->source, $methods, $ids);
	    }
	}
	
	return $result;
	
} # end export_seq_features

=head2 export_seq_aln_features

  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(export_seq_features);
               my $rst = export_seq_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

#sub export_seq_aln_features
#{
#	my ($self) = shift;
#	my ($operation) = shift;
#	my ($format) = shift;
#	my ($methods) = shift if (@_);
#	my ($ids) = shift if (@_);	
#	my ($result) = '';	
#	
#	# get protein fasta
#	my ($seq_result) = $self->export_seq_features('sequences', 'fasta', 'aa', $methods, $ids);
#	my ($num_seq) = $seq_result =~ tr/\>//;
#		
#	# get annotated residues
#	my ($res_result) = $self->export_seq_features('residues', 'fasta', 'aa', $methods, $ids);
#	
#	# retrieve align from REST service of Muscle(EBI)
#	my ($wsretriever) = new WSRetriever();
#	if ( $operation eq 'align' ) {
#		if ( $num_seq <= 1 ) { # "A minimum of 2 sequences is required\n"
#			$operation = 'sequences'; # sequence operator
#		}
#		else { # overwrite seq report by align report
#			if ( defined $seq_result ) {
#				$seq_result = $wsretriever->get_tool_result($seq_result, 'muscle', 'clwstrict');		
#			}			
#		}	
#	}
#	
#	# get sequence/align annotations
#	if ( defined $res_result and defined $seq_result ) {
#		$result = $wsretriever->get_seq_aln_features($operation, $seq_result, $res_result, $ids, $format);
#	}
#				
#	return $result;
#	
#} # end export_seq_aln_features

=head2 get_gen_annotations

  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(export_seq_features);
               my $rst = export_seq_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_gen_annotations
{
	my ($self) = shift;
	my ($source) = shift if (@_);
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);
	my ($result) = '';
	
	require WSRetriever;
	my ($wsretriever) = new WSRetriever();
	my ($query_id) = 'exporter/' . $self->type . '/' . $self->species . '/' . $self->input;
	$result = $wsretriever->get_gen_features($query_id, $self->species, $self->assembly, $source, $self->dataset, $methods, $ids);
	unless ( defined $result ) {
		$result = "Job query needs genome information\n";
	}
			
	return $result;
	
} # end get_gen_annotations

##
## SEEKER
##

=head2 seek_features

  Arg [1]    : String $type
               type of input: 'id' or 'name'
  Arg [2]    : String $input
               gene or transcript identifier (or name)
  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(seek_features);
               my $rst = seek_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub seek_features
{
	my ($self) = shift;
	my ($format) = shift;
	my ($result) = '';	
	my ($features) = $self->get_seek_features('none');
	
	if ( defined $features ) {
		my ($wsretriever) = new WSRetriever();		
		my ($report) = $wsretriever->create_seeker_report($self->input, $features);
		
		if ( defined $report ) {
			if ($format eq 'json') {
				require JSON;
				my ($json) = new JSON;
				$result = $json->encode($report);
		    }
			elsif ($format eq 'xml') {
				require XML::LibXML;
				my ($xml_doc) = $wsretriever->get_xml_report($report);
				if ( defined $xml_doc ) {
					$result = $xml_doc->toString(1);
				}
		    }
		}		
	}
	
	return $result;
	
} # end seek_features








#sub seek_features
#{
#	my ($self) = shift;
#	my ($format) = shift;
#	my ($result) = '';	
#	my ($features) = $self->get_seek_features('none');
#	
#	if ( defined $features ) {
#		
#		my ($report) = $self->create_seeker_report($features);
#		if ( defined $report ) {
#			if ($format eq 'json') {
#				require JSON;
#				my ($json) = new JSON;
#				$result = $json->encode($report);
#		    }
#			elsif ($format eq 'xml') {
#				require XML::LibXML;
#				my ($xml_doc) = $self->get_xml_report($report);
#				if ( defined $xml_doc ) {
#					$result = $xml_doc->toString(1);
#				}
#		    }
#		}		
#	}
#	
#	return $result;
#	
#} # end seek_features

#sub create_seeker_report($)
#{
#	my ($self) = shift;
#	my ($features) = shift;
#	my ($report);
#	$report->{'query'} = $self->input;
#	
#	while ( my ($specie, $q_report) = each(%{$features}) )
#	{
#		foreach my $entity (@{$q_report})
#		{
#			if ( $entity and ref($entity) ) {
#				my ($match) = {
#					'specie'	=> $specie,
#					'id'		=> $entity->stable_id,
#					'version'	=> $entity->version,
#					'label'		=> $entity->stable_id,
#					'chr'		=> $entity->chromosome,
#					'start'		=> $entity->start,
#					'end'		=> $entity->end,
#					'biotype'	=> $entity->biotype,
#					'status'	=> $entity->status,
#				};
#				if ( $entity->isa("APPRIS::Gene") ) {
#					$match->{'namespace'} = 'Gene_Id';	
#				}
#				elsif ( $entity->isa("APPRIS::Transcript") ) {
#					$match->{'namespace'} = 'Transcript_Id';
#				}
#				if ( $entity->external_name ) {
#					my ($dblink) = {
#						'id'		=> $entity->external_name,
#						'namespace'	=> 'External_Id'
#					};
#					push(@{$match->{'dblink'}}, $dblink);
#				}
#				if ( $entity->xref_identify ) {
#					foreach my $xref_identify (@{$entity->xref_identify}) {
#						my ($dblink) = {
#							'id'		=> $xref_identify->id,
#							'namespace'	=> $xref_identify->dbname
#						};
#						push(@{$match->{'dblink'}}, $dblink);
#					}			
#				}
#				if ( ($entity->isa("APPRIS::Gene")) and $entity->transcripts ) {
#					foreach my $transcript (@{$entity->transcripts}) {
#						my ($dblink) = {
#							'id'		=> $transcript->stable_id,
#							'namespace'	=> 'Transcript_Id'
#						};
#						push(@{$match->{'dblink'}}, $dblink);
#					}
#				}
#				push(@{$report->{'match'}}, $match);
#			}
#			else {
#				return $report; # Empty return
#			}
#		}
#	}
#	
#	return $report;
#	
#} # end create_seeker_report
#
#sub get_xml_report($)
#{
#	my ($self) = shift;
#	my ($report) = shift;	
#	my ($xml_doc) = XML::LibXML::Document->new('1.0','UTF-8');
#	
#	my ($e_query) = $xml_doc->createElement('query');
#	$e_query->setNamespace("https://appris-tools.org", "appris", 0);
#	$xml_doc->setDocumentElement($e_query);
#	my ($a_query) = $xml_doc->createAttribute('query', $report->{'query'});
#	$e_query->setAttributeNode($a_query);
#
#	foreach my $match (@{$report->{'match'}})
#	{
#		my ($e_match) = $xml_doc->createElement('match');
#		my ($a_match) = $xml_doc->createAttribute('label', $match->{'label'});
#		my ($a_match2) = $xml_doc->createAttribute('namespace', $match->{'namespace'});
#		my ($a_match3) = $xml_doc->createAttribute('chr', $match->{'chr'});
#		my ($a_match4) = $xml_doc->createAttribute('start', $match->{'start'});
#		my ($a_match5) = $xml_doc->createAttribute('end', $match->{'end'});
#		my ($a_match6) = $xml_doc->createAttribute('specie', $match->{'specie'});
#		$e_match->setAttributeNode($a_match);
#		$e_match->setAttributeNode($a_match2);
#		$e_match->setAttributeNode($a_match3);
#		$e_match->setAttributeNode($a_match4);
#		$e_match->setAttributeNode($a_match5);
#		$e_match->setAttributeNode($a_match6);
#
#		my ($e_class) = $xml_doc->createElement('biotype');
#		my ($e_class_txt) = $xml_doc->createCDATASection($match->{'biotype'});
#		$e_class->appendChild($e_class_txt);
#		$e_match->appendChild($e_class);
#
#		my ($e_status) = $xml_doc->createElement('status');
#		my ($e_status_txt) = $xml_doc->createCDATASection($match->{'status'});
#		$e_status->appendChild($e_status_txt);		
#		$e_match->appendChild($e_status);
#		
#		foreach my $dblink (@{$match->{'dblink'}}) {
#			my ($e_dblink) = $xml_doc->createElement('dblink');
#			my ($a_dblink) = $xml_doc->createAttribute('namespace', $dblink->{'namespace'});
#			my ($a_dblink2) = $xml_doc->createAttribute('id', $dblink->{'id'});
#			$e_dblink->setAttributeNode($a_dblink);
#			$e_dblink->setAttributeNode($a_dblink2);
#			$e_match->appendChild($e_dblink);
#		}
#		$e_query->appendChild($e_match);
#	}
#
#	return $xml_doc;
#}










#
#=head2 get_tool_result
#
#  Arg [1]    : String $type
#               type of input: 'id' or 'name'
#  Arg [2]    : String $input
#               gene or transcript identifier (or name)
#  Arg [3]    : String $method
#               method name
#  Example    : use DBRetriever qw(get_method_result);
#               my $rst = get_residues_result($type, $input, $method);
#  Description: Throws an exception which if not caught by an eval will
#               provide a stack trace to STDERR and die.  If the verbosity level
#               is lower than the level of the throw, then no error message is
#               displayed but the program will still die (unless the exception
#               is caught).
#  Returntype : none
#  Exceptions : thrown every time
#  Caller     : generally on error
#
#=cut
#sub get_tool_result {
#	my ($self) = shift;
#	my ($sequences, $type, $format) = @_;
#	my ($result) = undef;
#
#	my ($seq_tmpfile) = File::Temp->new( UNLINK => 0, SUFFIX => '.web_tool.appris' );
#	print $seq_tmpfile $sequences;
#
#	my ($cmd) = "perl $FindBin::Bin/";
#	if ( $type eq 'blastp' ) {
#		$cmd .= "ncbi_rest/web_blast.pl $type nr $seq_tmpfile";		
#	}
#	elsif ( $type eq 'muscle' ) {
#		$cmd .= "ebi_rest/muscle_lwp.pl --async --email appris\@cnio.es --format $format $seq_tmpfile ";		
#	}
#	else {
#		$result = undef;
#	}	
#	eval {
#		my (@cmd_stdout_list) = `$cmd`;
#		my ($cmd_stdout) = $cmd_stdout_list[0];
#		if ( defined $cmd_stdout and $cmd_stdout =~ /NCBI:([^\n]*)/ ) { # ncbi blast
#			my ($job_id) = $1;
#			$result = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=$job_id";
#		}
#		elsif ( defined $cmd_stdout and $cmd_stdout =~ /($type\-[^\n]*)/ ) { # ebi tools
#			my ($job_id) = $1;
#			$result = "https://www.ebi.ac.uk/Tools/services/web/toolresult.ebi?jobId=$job_id"
#		}
#	};
#	$seq_tmpfile->unlink_on_destroy(1);
#	
#	return $result;	
#}

sub _add_trifid_analyses($$) {
	my ($gene, $trifid_report) = @_;

	foreach my $transcript (@{$gene->transcripts})
	{
		my ($transc_id) = $transcript->stable_id;
		next unless( exists($trifid_report->{'_index_transcripts'}{$transc_id}) );
		my ($transc_report) = $trifid_report->transcript($transc_id);
		next unless( defined($transc_report->analysis) );

		my ($trifid_analysis) = $transc_report->analysis->trifid;
		$transcript->analysis(APPRIS::Analysis->new) if ( ! defined($transcript->analysis) );
		$transcript->analysis->trifid($trifid_analysis);
		$transcript->analysis->number($transcript->analysis->number+1);
	}
}

sub _fetch_trifid_report($$) {
	my ($gene, $pred_file_path) = @_;
	my ($report);

	my ($gene_id) = $gene->stable_id;
	my ($cmd) = "tabix -Dh $pred_file_path $gene_id";
	my (@lines) = `$cmd`;
	if ( scalar(@lines) > 1 ) {  # first line is header
		my $result = join('', @lines);
		$report = parse_trifid($gene, $result);
	}

	return $report;
}

1;
