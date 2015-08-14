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
use Config::IniFiles;
use File::Temp;
use Data::Dumper;

# these modules are loaded in load_registry function
#use APPRIS::Registry;
#use APPRIS::Exporter;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(info throw warning deprecate);

use WSRetriever;

{
	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = # DEFAULT
		(
			conf			=>  undef,
			specie			=>  undef,
			ens				=>  '', # emtpy for the last version of ensembl
			assembly		=>  undef,
			type			=>  undef,
			input			=>  undef,
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
	sub conf {
		my ($self, $arg) = @_;
		$self->{'conf'} = $arg if defined $arg;
		return $self->{'conf'};
	}
	
	sub specie {
		my ($self, $arg) = @_;
		$self->{'specie'} = $arg if defined $arg;
		return $self->{'specie'};
	}
	
	sub ens {
		my ($self, $arg) = @_;
		$self->{'ens'} = $arg if defined $arg;
		return $self->{'ens'};
	}

	sub assembly {
		my ($self, $arg) = @_;
		$self->{'assembly'} = $arg if defined $arg;
		return $self->{'assembly'};
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
}

=head2 new

  Arg [-conf]:
        string - config ini
  Arg [-specie]:
        string - specie name (human, mouse, rat)
  Arg [-ens]:
        integer - Ensembl version
  Arg [-assembly]:
        string - Assembly version
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

	my ( $conf, $specie, $ens, $assembly, $type, $input ) = rearrange( [ 'conf', 'specie', 'ens', 'assembly', 'type', 'input' ], @_ );

	# require paramater
	if ( defined $conf and -e $conf ) { $self->conf($conf); }
	else { return undef; }
	if ( defined $specie ) { $self->specie($specie); }
	#else { return undef; }

	# optional parameter
	if ( defined $ens ) { $self->ens($ens); }
	if ( defined $assembly ) { $self->assembly($assembly); }
	if ( defined $type ) { $self->type($type); }
	if ( defined $input ) { $self->input($input); }
	
	return $self;
}

=head2 is_registry

  Arg [1]    : String $specie
               specie name (humna, mouse)
  Arg [2]    : Integer $ens
               Ensembl version
  Example    : use DBRetriever qw(is_registry);
               my $rst = is_registry($specie);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut
sub is_registry {
	my ($self) = shift;
	my ($specie) = shift;
	my ($ens) = shift;
	my ($species);
	my ($is_registry);
	
	# config ini
	my ($cfg) = new Config::IniFiles( -file => $self->{'conf'} );
	
	# unless defined, then check all species
	if ( defined $specie ) {
		$species = [ $specie ];
	} else {
		$species = [ split( ',', $cfg->val('APPRIS_DATABASES', 'species') ) ];
	}
	
	foreach my $specie_ens (@{$species}) {
		if ( defined $ens and ($ens ne '') ) {
			$specie_ens .= '_ens'.$ens; 
		}
		else {
			# assign default version for each assembly/database. We get the first one
			my ($versions) = [ split( ',', $cfg->val($specie, 'versions') ) ];
			$ens = $versions->[0];
			$self->ens($ens);
			$specie_ens .= '_ens'.$ens; 		
		}
		
		my ($specie_db) = uc($specie_ens.'_db');
		my ($specie_dbname) = $cfg->val($specie_db, 'db');		
		if ( defined $specie_dbname ) {
			$is_registry = $specie_dbname;
		}		
	}
		
	return $is_registry;
}

=head2 load_registry

  Arg [1]    : String $specie
               specie name (humna, mouse)
  Arg [2]    : Integer $ens
               Ensembl version
  Example    : use DBRetriever qw(load_registry);
               my $rst = load_registry($specie);
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
	my ($specie) = shift;
	my ($ensembl) = shift;
	my ($assembly) = shift if (@_);
	my ($registry);
	
	# config ini
	my ($cfg) = new Config::IniFiles( -file => $self->{'conf'} );
	
	# if defined 'assembly' we overwrite the official version of ensembl
	my ($ens) = $ensembl;
	if ( defined $assembly ) {
		foreach my $version ( split( ',', $cfg->val($specie, 'versions') ) ) {
			my ($species_ens) = $specie.'_ens'.$version;
			my ($species_db) = uc($species_ens.'_db');
			my ($species_assembly) = $cfg->val($species_db, 'assembly');				
			if ( !defined $assembly or ( defined $assembly and (lc($species_assembly) eq $assembly) ) ) {
				$ens = $version;
			}
		}		
	}
		
	# get database name for available species
	my ($specie_ens) = $specie;
	if ( defined $ens and ($ens ne '') ) {
		$specie_ens .= '_ens'.$ens; 
	}
	else {
		# assign default version for each assembly/database. We get the first one
		my ($versions) = [ split( ',', $cfg->val($specie, 'versions') ) ];
		$ens = $versions->[0];
		$self->ens($ens);
		$specie_ens .= '_ens'.$ens; 		
	}
	
	# load database and specific libs
	my ($specie_db) = uc($specie_ens.'_db');
	my ($specie_dbname) = $cfg->val($specie_db, 'db');
	if ( defined $specie_dbname ) {
		
		# load perllib dynamically!!!!!
		my ($appris_perllib) = $cfg->val($specie_db, 'perllib');
		eval "use lib qw($FindBin::Bin/../lib/$appris_perllib);
		require APPRIS::Registry;
		require APPRIS::Exporter;";
		if ($@) { die "Error loading APPRIS modules:\n$@\n"; }
				
		# load registry
		$registry = APPRIS::Registry->new();		
		$registry->load_registry_from_db(
								-dbhost	=> $cfg->val('APPRIS_DATABASES', 'host'),
								-dbuser	=> $cfg->val('APPRIS_DATABASES', 'user'),
								-dbpass	=> $cfg->val('APPRIS_DATABASES', 'pass'),
								-dbport	=> $cfg->val('APPRIS_DATABASES', 'port'),
								-dbname	=> $cfg->val($specie_db, 'db'),
		);		
	}
	
	return $registry;
}	

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
	my ($ids) = shift if (@_);
	my ($type) = $self->type;
	my ($inputs) = $self->input;
	my ($specie) = $self->specie;
	my ($ensembl) = $self->ens if ( $self->ens );
	my ($assembly) = lc($self->assembly) if ( $self->assembly );	
	my ($features);
	
	# We the list of methods will be able to contain the type of track: Eg. spade:domain|damaged_domain,firestar
	my ($methods_str) = '';
	foreach my $method (split(',', $methods)) {
		if ( $method =~ /^([^\:|\$]*)/ ) { $methods_str .= $1 . ',' }
	}
	$methods_str =~ s/\,$//;
	 
	my ($registry) = $self->load_registry($specie,$ensembl,$assembly);
	return undef unless (defined $registry);
		
	if ( defined $type ) {
		if ( ($type eq 'id') and $inputs ) {
			if ( defined $ids ) { $inputs = $ids }
			foreach my $input (split(',', $inputs)) {
				my ($feat) = $self->get_feat_by_stable_id($registry, $input, $methods_str);
				if ( defined $feat ) {
					push(@{$features}, $feat) ;
				}
			}
		}
		elsif ( ($type eq 'name') and $inputs ) {
			if ( defined $ids ) { $inputs = $ids }
			foreach my $input (split(',', $inputs)) {
				my ($feat) = $self->get_feat_by_xref_entry($registry, $input, $methods_str);
				if ( defined $feat and scalar(@{$feat}) > 0 ) {
					foreach my $f (@{$feat}) { push(@{$features}, $f); }
				}
			}
		}
		elsif ( ($type eq 'position') and $inputs ) {
			if ( defined $ids ) { $inputs = $ids }
			foreach my $input (split(',', $inputs)) {
				my ($feat) = $self->get_feat_by_region($registry, $input, $methods_str);
				if ( defined $feat and scalar(@{$feat}) > 0 ) {
					foreach my $f (@{$feat}) { push(@{$features}, $f); }
				}
			}
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
	my ($ensembl) = $self->ens if ( $self->ens );
	my ($assembly) = lc($self->assembly) if ( $self->assembly );
	my ($features);
	
	# check all species using default ensembl version
	my ($cfg) = new Config::IniFiles( -file => $self->{'conf'} );
	my ($species) = [ split( ',', $cfg->val('APPRIS_DATABASES', 'species') ) ];	
	foreach my $spe (@{$species}) {
		my ($ens_dbs) = $cfg->val($spe, 'dbs');
		if ( defined $ensembl ) {
			$ens_dbs = $ensembl;
		}
		else {
			$ens_dbs = $cfg->val($spe, 'dbs');
		}
		foreach my $ens ( split(',',$ens_dbs) ) {
			my ($registry) = $self->load_registry($spe, $ens, $assembly);
			if ( defined $registry ) {
				my ($species_ens) = $spe.'_ens'.$ens;
				my ($species_db) = uc($species_ens.'_db');
				my ($species_assembly) = $cfg->val($species_db, 'assembly');				
				#if ( !defined $assembly or ( defined $assembly and (lc($species_assembly) =~ /^$assembly\|/ or lc($species_assembly) =~ /\|$assembly$/) ) ) {
					if ( !defined $assembly or ( defined $assembly and (lc($species_assembly) eq $assembly) ) ) {				
					foreach my $input (split(',', $inputs)) {
						my ($feat) = $self->get_feat_by_xref_entry($registry, $input, $methods);
						#foreach my $f (@{$feat}) { push(@{$features->{$spe}->{$species_assembly}}, $f); }
						foreach my $f (@{$feat}) {
							push(@{$features}, {
								'species'	=> $spe,
								'assembly'	=> $species_assembly,
								'dataset'	=> $ens,
								'entity'	=> $f
							});
						}
					}
				}			
			}			
		}
	}
		
	return $features;
	
} # end get_seek_features


sub get_feat_by_stable_id {
	my ($self) = shift;
	my ($registry) = shift;
	my ($id) = shift;
	my ($methods) = shift;
	my ($feat);

	if ( lc($id) =~ /^ens(\w\w\w)?g(\d){11,13}/ or
		 lc($id) =~ /^ens(\w\w\w)?gr(\d){10,13}/ or
		 lc($id) =~ /^ensgr(\d){10,13}/ or
		 lc($id) =~ /^fbgn(\d){7}/ or
		 lc($id) =~ /^wbgene(\d){8}/		 
	) {
		$feat = $registry->fetch_by_stable_id('gene', $id, $methods);
	}
	elsif (  lc($id) =~ /^ens(\w\w\w)?t(\d){11,13}/ or
			 lc($id) =~ /^ens(\w\w\w)?tr(\d){10,13}/ or
			 lc($id) =~ /^enstr(\d){10,13}/ or
			 lc($id) =~ /^fbtr(\d){7}/
	) {
		$feat = $registry->fetch_by_stable_id('transcript', $id, $methods);
	}
	
	return $feat;	
}

sub get_feat_by_xref_entry {
	my ($self) = shift;
	my ($registry) = shift;
	my ($name) = shift;
	my ($methods) = shift;
	
	my ($feat) = $registry->fetch_by_xref_entry($name, $methods);
	
	return $feat;	
}

sub get_feat_by_region {
	my ($self) = shift;
	my ($registry) = shift;
	my ($input) = shift;
	my ($methods) = shift;
	my ($chr,$start,$end) = (undef,undef,undef);	
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
	my ($typebed) = shift if (@_);
	my ($result) = '';
	my ($features) = $self->get_features($methods, $ids); # if methods is not defined => retrieve all features
	
	if ( defined $features ) {
		my ($exporter) = APPRIS::Exporter->new();
		if ($format eq 'tsv') {
			$result = $exporter->get_tsv_annotations($features, $methods, $res);
		}
		elsif ($format eq 'bed') {
			$result = $exporter->get_bed_annotations($features, $methods, undef, $typebed);
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
			$result = $self->get_gen_annotations($methods, $ids);
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
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);
	my ($result) = '';
	
	# assign default ensembl version from assembly. Otherwise, the first value.
	unless ( $self->ens ) {
		my ($cfg) = new Config::IniFiles( -file => $self->conf );
		my ($assembly) = $self->assembly;
		my ($specie) = $self->specie;
		my ($versions) = [ split( ',', $cfg->val($specie, 'versions') ) ];
		my ($ens) = $versions->[0];		
		if ( defined $assembly ) {
			foreach my $version ( @{$versions} ) {
				my ($species_ens) = $specie.'_ens'.$version;
				my ($species_db) = uc($species_ens.'_db');
				my ($species_assembly) = $cfg->val($species_db, 'assembly');				
				if ( !defined $assembly or ( defined $assembly and (lc($species_assembly) eq $assembly) ) ) {
					$ens = $version;
				}
			}		
		}
		$self->ens($ens);
	}
	
	# retrieve images from UCSC
	my ($wsretriever) = new WSRetriever();
	if ( $self->specie and $self->ens ) {
		my ($query_id) = 'exporter/' . $self->type . '/' . $self->specie . '/' . $self->input;
		$result = $wsretriever->get_gen_features($query_id, $self->specie, $self->ens, $methods, $ids);
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
#					$match->{'namespace'} = 'Ensembl_Gene_Id';	
#				}
#				elsif ( $entity->isa("APPRIS::Transcript") ) {
#					$match->{'namespace'} = 'Ensembl_Transcript_Id';
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
#							'namespace'	=> 'Ensembl_Transcript_Id'
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
#	$e_query->setNamespace("http://appris.bioinfo.cnio.es", "appris", 0);
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
#			$result = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=$job_id";
#		}
#		elsif ( defined $cmd_stdout and $cmd_stdout =~ /($type\-[^\n]*)/ ) { # ebi tools
#			my ($job_id) = $1;
#			$result = "http://www.ebi.ac.uk/Tools/services/web/toolresult.ebi?jobId=$job_id"
#		}
#	};
#	$seq_tmpfile->unlink_on_destroy(1);
#	
#	return $result;	
#}

1;
