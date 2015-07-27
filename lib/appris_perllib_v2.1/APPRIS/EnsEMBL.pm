=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::EnsEMBL

=head1 SYNOPSIS

  $registry = APPRIS::EnsEMBL->new(
    -dbhost  => 'localhost',
    -dbname  => 'homo_sapiens_encode_3c',
    -dbuser  => 'jmrodriguez'
    );

  $gene = $registry->fetch_by_stable_id($stable_id);

  @genes = @{ $registry->fetch_by_chr_start_end('X', 1, 10000) };

=head1 DESCRIPTION

All Adaptors are stored/registered using this module.
This module should then be used to get the adaptors needed.

The registry can be loaded from a configuration file or from database info.

=head1 METHODS

=cut

package APPRIS::EnsEMBL;

use strict;
use warnings;
use Data::Dumper;
use URI::Escape;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;

use APPRIS::Export::GTF qw( print_annotations );
use APPRIS::Utils::Exception qw( throw warning deprecate );

# Constants
my ($GEN_SOURCE) = 'ENSEMBL';
my ($PROT_SEQNAME) = 'SEQ';
my ($PROT_SOURCE) = 'Ensembl';
my ($PROT_TYPE) = 'Protein';

{
    # Encapsulated class data
    #___________________________________________________________
    my %_attr_data = # DEFAULT
		(
			db_version	=> undef,		
			host		=> 'ensembldb.cnio.es',	# 'ensembldb.ensembl.org'
			user		=> 'ensembl',			# 'anonymous'
			pass		=> '',			# 'anonymous'
			verbose		=> 0,
			species		=> 'Homo sapiens',	
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
}

=head2 new

  Arg [-host]  : 
       String - the name of the host where the Ensembl database lives
  Arg [-user]  : 
       String - the user name used to access the database
  Arg [-species]    : 
       String - name of the species to add the adaptor to in the registry
       
  Example    : $ensembl = APPRIS::EnsEMBL->new(...);
  Description: creates a new ensembl registry object
  Returntype : APPRIS::EnsEMBL
  Exceptions : thrown if the given MySQL database cannot be connected to
               or there is any error whilst querying the database.
  Status     : stable

=cut

sub new {
	my ($caller, %args) = @_;
	
	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;
	
	# get input or default values
	foreach my $attrname ($self->_standard_keys) {
		my ($attr) = "-".$attrname;
		if (exists $args{$attr} && defined $args{$attr}) {
			$self->{$attrname} = $args{$attr};
		} else {
			$self->{$attrname} = $self->_default_for($attrname);
		}
	}


	# Ensembl connection details
	my ($registry) = 'Bio::EnsEMBL::Registry';
	eval {
		$registry->load_registry_from_db(
				-db_version	=> $self->{'db_version'},
				-host		=> $self->{'host'},
				-user		=> $self->{'user'},
				-pass		=> $self->{'pass'},
				-verbose	=> $self->{'verbose'},
				-species	=> $self->{'species'},				
		);
	};
	throw('Ensembl connection') if $@;
	$self->_registry($registry);
	
	# create gene adaptor from specie
	$self->_get_gene_adaptor($self->{'species'});
	

	return $self;
}

###################
# PRIVATE methods #
###################

=head2 _registry

  Example    : $ent = $ens->_registry("21");
  Description: setter for Ensembl registry
  Returntype : registry adaptor
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub _registry {
	my ($self) = shift;
	$self->{'_registry'} = shift if(@_);
	return $self->{'_registry'};
}

=head2 _get_gene_adaptor

  Example    : $ent = $ens->_get_gene_adaptor();
  Description: setter object adaptor (gene)
  Returntype : gene adaptor or undef
  Exceptions : DBAdaptor instances which can be found in the registry
  Caller     : general
  Status     : stable

=cut

sub _get_gene_adaptor {
	my ($self) = shift;
	my ($specie) = shift;
	my ($adaptor) = $self->{'_registry'}->get_adaptor($specie,'core','gene');
	if ( defined $adaptor ) {
		$self->{'_gene_adaptor'} = $adaptor;
	}
	else { throw('Undefined adaptor'); }
	return $self->{'_gene_adaptor'};
}

=head2 _get_ccdsid_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : $ccds_id = _get_ccdsid_by_transc_adaptor($entity);
  Description: get the CCDS identifier from transcript object
  Returntype : String or undef
  Exceptions : none

=cut

sub _get_ccdsid_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($ccds_id);
	if ( $entity and (ref($entity) ne 'ARRAY') ) {			
	   	if ( $entity->isa("Bio::EnsEMBL::Transcript") ) {
			my (@dblinks) = @{ $entity->get_all_DBLinks };
			foreach my $entry (@dblinks) {
				if ( $entry->dbname eq 'CCDS' ) {
					$ccds_id = $entry->primary_id().'.'.$entry->version();
					last;
				}
			}
	   	}
	}
	return $ccds_id;
}

=head2 fetch_by_stable_id

  Arg [1]    : String - the stable identifier of gene
  Example    : $gene = $ens->fetch_by_stable_id(id);
  Description: getter Ensembl gene objet
  Returntype : Bio::EnsEMBL::Gene objects or undef
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub fetch_by_stable_id {
	my ($self) = shift;
	my ($id) = shift;	
	my ($gene) = $self->{'_gene_adaptor'}->fetch_by_stable_id($id);				
	if ( defined $gene ) {
		$self->{'gene'} = $gene;
	}
	else { throw('Undefined gene'); }
	return $self->{'gene'};
}

=head2 fetch_all_by_external_name

  Arg [1]    : String - the external name of gene
  Example    : $gene = $ens->fetch_all_by_external_name(id);
  Description: getter Ensembl gene objet
  Returntype : Bio::EnsEMBL::Gene objects or undef
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub fetch_all_by_external_name {
	my ($self) = shift;
	my ($id) = shift;
	my (@genes) = @{ $self->{'_gene_adaptor'}->fetch_all_by_external_name($id) };
	if ( scalar(@genes) > 0 ) {
		my ($gene) = $genes[0];
		$self->{'gene'} = $gene;
	}
	else { throw('Undefined gene'); }
	return $self->{'gene'};
}

=head2 get_all_transcripts

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $transcripts = $ens->get_all_transcripts(id);
  Description: getter list of Ensembl transcript objet
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : No transcripts
  Caller     : general
  Status     : stable

=cut

sub get_all_transcripts {
	my ($self) = shift;
	my ($gene) = shift;

	my (@transcripts) = @{ $gene->get_all_Transcripts };
	if ( scalar(@transcripts) > 0 ) {
		$self->{'transcripts'} = \@transcripts;
	}
	else { throw('No transcripts'); }
	return $self->{'transcripts'};
}

=head2 _get_start_codon_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : $content = $ens->_get_start_codon_by_transc_adaptor($entity);
  Description: get start codon coordinates, if applied
  Returntype : Coord object or empty
  Exceptions : None
  Caller     : general
  Status     : stable

=cut

sub _get_start_codon_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($data);
	
	if ( defined $entity and defined $entity->stable_id() ) {
		
		# check if cds start/stop has founded
		my ($hasstart, $hasend) = _check_start_and_stop($entity);
		
		# get transcript annotation remarks to check NOT FOUND tags		
		my ($remarks) = _transcript_remarks($entity);

		# get codon feature		
		if ( $entity->translation && $hasstart ) {
			if ( ! ($remarks =~ /(cds_start_NF)/) ) {
				my (@codons) =  _make_start_codon_features($entity, $entity->stable_id());
				if ( scalar(@codons) > 0 ) {					
					my ($codon) = $codons[0];
					my ($sliceoffset) = $entity->slice->start-1;
					my ($start) = ($codon->start+$sliceoffset);
					my ($end) =  ($codon->end+$sliceoffset);
					my ($strand) = '.';
					if ($entity->strand() == 1) { $strand = '+' }
					elsif ($entity->strand() == -1) { $strand = '-' }
					my ($phase) = $codon->phase;
					if ( $phase >= 3 ) { $phase = 0 }
					$data->{'start'}	= $start;
					$data->{'end'}		= $end;
					$data->{'strand'}	= $strand;
					$data->{'phase'}	= $phase;		
				}
			}		
		}
	}
	return $data;
}

=head2 _get_stop_codon_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : $content = $ens->_get_stop_codon_by_transc_adaptor($entity);
  Description: get stop codon coordinates, if applied
  Returntype : Coord object or empty
  Exceptions : None
  Caller     : general
  Status     : stable

=cut

sub _get_stop_codon_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($data);
	
	if ( defined $entity and defined $entity->stable_id() ) {
		
		# check if cds start/stop has founded
		my ($hasstart, $hasend) = _check_start_and_stop($entity);
		
		# get transcript annotation remarks to check NOT FOUND tags		
		my ($remarks) = _transcript_remarks($entity);

		# get codon feature
		if ( $entity->translation && $hasend ) {
			if ( ! ($remarks =~ /(cds_end_NF)/) ) {
				my (@codons) =  _make_stop_codon_features($entity, $entity->stable_id());
				if ( scalar(@codons) > 0 ) {					
					my ($codon) = $codons[0];
					my ($sliceoffset) = $entity->slice->start-1;
					my ($start) = ($codon->start+$sliceoffset);
					my ($end) =  ($codon->end+$sliceoffset);
					my ($strand) = '.';
					if ($entity->strand() == 1) { $strand = '+' }
					elsif ($entity->strand() == -1) { $strand = '-' }
					my ($phase) = $codon->phase;
					if ( $phase >= 3 ) { $phase = 0 }
					$data->{'start'}	= $start;
					$data->{'end'}		= $end;
					$data->{'strand'}	= $strand;
					$data->{'phase'}	= $phase;		
				}
			}		
		}
	}	
	return $data;
}

=head2 get_gtf_by_adaptor

  Arg [1]    : Bio::EnsEMBL::Gene or Bio::EnsEMBL::Transcript
  Example    : $content = $ens->get_gtf_by_adaptor($entity);
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_by_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($output) = '';
	if ( $entity and (ref($entity) ne 'ARRAY') ) {			
	   	if ( $entity->isa("Bio::EnsEMBL::Gene") ) {
   			$output = $self->get_gtf_by_gene_adaptor($entity);
	   	}
	   	elsif ( $entity->isa("Bio::EnsEMBL::Transcript") ) {
			$output = $self->get_gtf_by_transc_adaptor($entity);
	   	}
	}
	return $output;
}

=head2 get_gtf_by_gene_adaptor

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $content = $ens->get_gtf_by_gene_adaptor($entity);
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_by_gene_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		my ($seqname) = 'chr'.$entity->slice->seq_region_name();
		my ($source) = $GEN_SOURCE;
		my ($type) = 'gene';
		my ($start) = $entity->start();
		my ($end) = $entity->end();
		my ($strand) = '.';
		if ($entity->strand() == 1) { $strand = '+' }
		elsif ($entity->strand() == -1) { $strand = '-' }
		my ($score) = '.';
		my ($phase) = '.';
		my ($data) = {
			'seqname'	=> $seqname,
			'source'	=> $source,
			'type'		=> $type,
			'start'		=> $start,
			'end'		=> $end,
			'score'		=> $score,
			'strand'	=> $strand,
			'phase'		=> $phase			
		};
		my ($optional);
		my ($gene_id) =  $entity->stable_id().'.'. $entity->version();
		$optional->{'gene_id'}			= $gene_id;
		$optional->{'transcript_id'}	= $gene_id;
		$optional->{'gene_name'}		= $entity->external_name();
		$optional->{'gene_type'}		= $entity->biotype();
		$optional->{'gene_status'}		= $entity->status();
		
		$output = APPRIS::Export::GTF::print_annotations($data, $optional);
	}

	return $output;
}

=head2 get_gtf_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_gtf_by_transc_adaptor($entity,'ENG..');
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		my ($seqname) = 'chr'.$entity->slice->seq_region_name();
		my ($source) = $GEN_SOURCE;
		my ($type) = 'transcript';
		my ($start) = $entity->start();
		my ($end) = $entity->end();
		my ($strand) = '.';
		if ($entity->strand() == 1) { $strand = '+' }
		elsif ($entity->strand() == -1) { $strand = '-' }
		my ($score) = '.';
		my ($phase) = '.';
		my ($data) = {
			'seqname'	=> $seqname,
			'source'	=> $source,
			'type'		=> $type,
			'start'		=> $start,
			'end'		=> $end,
			'score'		=> $score,
			'strand'	=> $strand,
			'phase'		=> $phase			
		};
		my ($optional);
		unless ( defined $gene_id ) {
			$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
		}
		my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
		$optional->{'gene_id'}				= $gene_id;
		$optional->{'transcript_id'}		= $transcript_id;
		$optional->{'transcript_name'}		= $entity->external_name();
		$optional->{'transcript_type'}		= $entity->biotype();
		$optional->{'transcript_status'}	= $entity->status();
		my ($ccds_id) = $self->_get_ccdsid_by_transc_adaptor($entity);
		if ( defined $ccds_id ) {
			$optional->{'ccdsid'}	= $ccds_id;	
		}
		if ( defined $entity->translation() and defined $entity->translation()->stable_id() ) {
			$optional->{'peptide_id'} = $entity->translation()->stable_id().'.'.$entity->translation()->version();
		}
		
		$output = APPRIS::Export::GTF::print_annotations($data, $optional);
	}
	
	return $output;
}

=head2 get_gtf_start_codon_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_gtf_start_codon_by_transc_adaptor($entity,'ENG..');
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_start_codon_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		my ($codon) = $self->_get_start_codon_by_transc_adaptor($entity);		
		if ( defined $codon ) {
			my ($seqname) = 'chr'.$entity->slice->seq_region_name();
			my ($source) = $GEN_SOURCE;
			my ($type) = 'start_codon';
			my ($start) = $codon->{'start'};
			my ($end) =  $codon->{'end'};
			my ($strand) =  $codon->{'strand'};
			my ($score) = '.';
			my ($phase) = $codon->{'phase'};
			my ($data) = {
				'seqname'	=> $seqname,
				'source'	=> $source,
				'type'		=> $type,
				'start'		=> $start,
				'end'		=> $end,
				'score'		=> $score,
				'strand'	=> $strand,
				'phase'		=> $phase			
			};
			my ($optional);
			unless ( defined $gene_id ) {
				$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
			}
			my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
			$optional->{'gene_id'}				= $gene_id;
			$optional->{'transcript_id'}		= $transcript_id;
			$output .= APPRIS::Export::GTF::print_annotations($data, $optional);
		}
	}
	
	return $output;
}

=head2 get_gtf_stop_codon_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_gtf_stop_codon_by_transc_adaptor($entity,'ENG..');
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_stop_codon_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		my ($codon) = $self->_get_stop_codon_by_transc_adaptor($entity);		
		if ( defined $codon ) {
			my ($seqname) = 'chr'.$entity->slice->seq_region_name();
			my ($source) = $GEN_SOURCE;
			my ($type) = 'stop_codon';
			my ($start) = $codon->{'start'};
			my ($end) =  $codon->{'end'};
			my ($strand) =  $codon->{'strand'};
			my ($score) = '.';
			my ($phase) = $codon->{'phase'};
			my ($data) = {
				'seqname'	=> $seqname,
				'source'	=> $source,
				'type'		=> $type,
				'start'		=> $start,
				'end'		=> $end,
				'score'		=> $score,
				'strand'	=> $strand,
				'phase'		=> $phase			
			};
			my ($optional);
			unless ( defined $gene_id ) {
				$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
			}
			my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
			$optional->{'gene_id'}				= $gene_id;
			$optional->{'transcript_id'}		= $transcript_id;
			$output .= APPRIS::Export::GTF::print_annotations($data, $optional);
		}
	}
	
	return $output;
}

=head2 get_gtf_exon_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_gtf_exon_by_transc_adaptor($entity,'ENG..');
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_exon_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		
		my (@transc_exons) = @{ $entity->get_all_Exons() };
		for (my $i = 0; $i < scalar(@transc_exons); $i++ ) {
			my ($transc_exon) = $transc_exons[$i];				
			my ($seqname) = 'chr'.$entity->slice->seq_region_name();
			my ($source) = $GEN_SOURCE;
			my ($type) = 'exon';
			my ($start) = $transc_exon->start();
			my ($end) = $transc_exon->end();
			my ($strand) = '.';
			if ($transc_exon->strand() == 1) { $strand = '+' }
			elsif ($transc_exon->strand() == -1) { $strand = '-' }
			my ($score) = '.';
			my ($phase) = '.';
			
			my ($data) = {
				'seqname'	=> $seqname,
				'source'	=> $source,
				'type'		=> $type,
				'start'		=> $start,
				'end'		=> $end,
				'score'		=> $score,
				'strand'	=> $strand,
				'phase'		=> $phase			
			};
			my ($optional);
			unless ( defined $gene_id ) {
				$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
			}
			my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
			my ($exon_id) =  $transc_exon->stable_id().'.'. $transc_exon->version();
			$optional->{'gene_id'}				= $gene_id;
			$optional->{'transcript_id'}		= $transcript_id;
			$optional->{'exon_id'}				= $exon_id;
						
			$output .= APPRIS::Export::GTF::print_annotations($data, $optional);			
		}
	}
	
	return $output;
}

=head2 get_gtf_cds_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_gtf_cds_by_transc_adaptor($entity,'ENG..');
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_cds_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		
		if ( $entity->translation ) {
			my (@transl_exons) = @{ $entity->get_all_translateable_Exons() };
			for (my $i = 0; $i < scalar(@transl_exons); $i++ ) {
				my ($transl_exon) = $transl_exons[$i];				
				my ($seqname) = 'chr'.$entity->slice->seq_region_name();
				my ($source) = $GEN_SOURCE;
				my ($type) = 'CDS';
				my ($start) = $transl_exon->start();
				my ($end) = $transl_exon->end();
				my ($strand) = '.';
				if ($transl_exon->strand() == 1) { $strand = '+' }
				elsif ($transl_exon->strand() == -1) { $strand = '-' }
				my ($score) = '.';
				my ($phase) = 3 - $transl_exon->phase();
				if ( $phase >= 3 ) { $phase = 0 }
				
				# check if stop codon exists
				if ( $i == scalar(@transl_exons)-1 ) {
					my ($codon) = $self->_get_stop_codon_by_transc_adaptor($entity);
					if ( defined $codon ) {
						if ( $strand eq '+' ) {
							$end = $codon->{'start'} -1;
						}
						if ( $strand eq '-' ) {
							$start = $codon->{'end'} +1;
						}
					}
				}
				my ($data) = {
					'seqname'	=> $seqname,
					'source'	=> $source,
					'type'		=> $type,
					'start'		=> $start,
					'end'		=> $end,
					'score'		=> $score,
					'strand'	=> $strand,
					'phase'		=> $phase			
				};
				my ($optional);
				unless ( defined $gene_id ) {
					$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
				}
				my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
				$optional->{'gene_id'}				= $gene_id;
				$optional->{'transcript_id'}		= $transcript_id;
							
				$output .= APPRIS::Export::GTF::print_annotations($data, $optional);			
			}
		}
	}
	
	return $output;
}

=head2 get_gtf_prot_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_gtf_cds_by_transc_adaptor($entity,'ENG..');
  Description: getter text annotation as GTF format
  Returntype : String or empty
  Exceptions : Undefined gene
  Caller     : general
  Status     : stable

=cut

sub get_gtf_prot_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		
		if ( $entity->translation ) {
			my ($trmapper) = Bio::EnsEMBL::TranscriptMapper->new($entity);
			my (@transl_exons) = @{ $entity->get_all_translateable_Exons() };
			for (my $i = 0; $i < scalar(@transl_exons); $i++ ) {
				my ($transl_exon) = $transl_exons[$i];
				my ($exon_id) = $transl_exon->stable_id();
				my ($seqname) = $PROT_SEQNAME;
				my ($source) = $PROT_SOURCE;
				my ($type) = $PROT_TYPE;
				
				# get peptide coords
				my ($pep_coords) = $trmapper->genomic2pep($transl_exon->start(), $transl_exon->end(), $transl_exon->strand());
				my ($start) = $pep_coords->start();
				my ($end) = $pep_coords->end();
				my ($strand) = '.';
				my ($score) = '.';
				my ($phase) = '.';
				
				# get phases of peptide
				my ($start_phase) = 3 - $transl_exon->phase();
				if ( $start_phase >= 3 ) { $start_phase = 0 }				
				my ($end_phase) = 3 - $transl_exon->end_phase();
				if ( $end_phase >= 3 ) { $end_phase = 0 }
								
				my ($data) = {
					'seqname'	=> $seqname,
					'source'	=> $source,
					'type'		=> $type,
					'start'		=> $start,
					'end'		=> $end,
					'score'		=> $score,
					'strand'	=> $strand,
					'phase'		=> $phase			
				};
				my ($optional);
				unless ( defined $gene_id ) {
					$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
				}
				my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
				$optional->{'gene_id'}				= $gene_id;
				$optional->{'transcript_id'}		= $transcript_id;
				$optional->{'exon_id'}				= $exon_id;
				$optional->{'start_cds'}			= $transl_exon->start();
				$optional->{'end_cds'}				= $transl_exon->end();
				$optional->{'start_phase'}			= $start_phase;
				$optional->{'end_phase'}			= $end_phase;
							
				$output .= APPRIS::Export::GTF::print_annotations($data, $optional);			
			}
		}
	}
	
	return $output;
}

=head2 get_transc_seq_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_transc_seq_by_transc_adaptor($entity,'ENG..');
  Description: getter transcript sequence as FASTA format
  Returntype : String or empty
  Exceptions : None
  Caller     : general
  Status     : stable

=cut

sub get_transc_seq_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		if ( defined $entity->seq() ) {
			my ($seqio) = $entity->seq();
			unless ( defined $gene_id ) {
				$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
			}
			my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
			my ($trans_name) = $entity->external_name();
			$output .= ">".$transcript_id."|".$gene_id."|".$trans_name."|".$seqio->length."\n";
			$output .= $seqio->seq."\n";			
		}
	}
	
	return $output;
}

=head2 get_transl_seq_by_transc_adaptor

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : String - the stable identifier of gene - Optional -
  Example    : $content = $ens->get_transl_seq_by_transc_adaptor($entity,'ENG..');
  Description: getter translate sequence as FASTA format
  Returntype : String or empty
  Exceptions : None
  Caller     : general
  Status     : stable

=cut

sub get_transl_seq_by_transc_adaptor {
	my ($self) = shift;
	my ($entity) = shift;
	my ($gene_id) = shift;
	my ($output) = '';
	
	if ( defined $entity and defined $entity->stable_id() ) {
		if ( $entity->translation() ) {
			my ($translation) = $entity->translation();
			if ( defined $translation->seq() ) {
				my ($seq) = $translation->seq();
				my ($len) = length($seq);
				unless ( defined $gene_id ) {
					$gene_id = $self->{'gene'}->stable_id().'.'.$self->{'gene'}->version();
				}
				my ($transcript_id) =  $entity->stable_id().'.'. $entity->version();
				my ($trans_name) = $entity->external_name();
				$output .= ">".$transcript_id."|".$gene_id."|".$trans_name."|".$len."\n";
				$output .= $seq."\n";			
			}
			
		}
	}
	
	return $output;
}

sub _make_start_codon_features {
	my ($trans, $id) = @_;

	return (()) if (!$trans->translation);

	my (@translateable) = @{$trans->get_all_translateable_Exons};
	my (@pepgencoords) = $trans->pep2genomic(1,1);

	if ( scalar(@pepgencoords) > 2 ) {
		#throw("pep start does not map cleanly: $id");
		warning("pep start does not map cleanly: $id");
		return (());
	}
	elsif ( scalar(@pepgencoords) == 2 ) {
		#throw("WOW got a 2 feature start codon for $id strand " . $translateable[0]->strand);
		warning("WOW got a 2 feature start codon for $id strand " . $translateable[0]->strand);
		return (());
 	}

	unless ( $pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
		#throw("pep start maps to gap: ".$pepgencoords[0]->start."-".$pepgencoords[0]->end." - ".$id." ".$trans->biotype."\n");
		warning("pep start maps to gap: ".$pepgencoords[0]->start."-".$pepgencoords[0]->end." - ".$id." ".$trans->biotype."\n");
		return (());
	}
	unless ( $pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
		#throw("pep start (end of) maps to gap: ".$pepgencoords[$#pepgencoords]->start."-".$pepgencoords[$#pepgencoords]->end." -".$id." ".$trans->biotype."\n");
		warning("pep start (end of) maps to gap: ".$pepgencoords[$#pepgencoords]->start."-".$pepgencoords[$#pepgencoords]->end." -".$id." ".$trans->biotype."\n");
		return (());
	}

	@translateable = @{$trans->get_all_translateable_Exons};
	my (@startc_feat);
	my ($phase) = 0;
	foreach my $pepgencoord ( @pepgencoords ) {
		push @startc_feat, new Bio::EnsEMBL::SeqFeature(
					                             -seqname => $id,
					                             -source_tag => 'starttrans',
					                             -primary_tag => 'similarity',
					                             -start => $pepgencoord->start,
					                             -end   => $pepgencoord->end,
					                             -phase => $phase,
					                             -strand => $translateable[0]->strand);
		$phase = 3 - ($pepgencoord->end - $pepgencoord->start + 1);
	}
	if ( $translateable[0]->strand == 1 ) {
		@startc_feat = sort { $a->start <=> $b->start } @startc_feat;
	}
	else {
		@startc_feat = sort { $b->start <=> $a->start } @startc_feat;
	}
	return @startc_feat;
}

sub _make_stop_codon_features {
	my ($trans,$id) = @_;
	
	return (()) if (!$trans->translation);
	
	my (@translateable) = @{$trans->get_all_translateable_Exons};
	my ($cdna_endpos) = $trans->cdna_coding_end;
	my (@pepgencoords) = $trans->cdna2genomic($cdna_endpos-2,$cdna_endpos);

	if ( scalar(@pepgencoords) > 2 ) {
		#throw("pep end does not map cleanly: $id");
		warning("pep end does not map cleanly: $id");
		return (());
	}
	elsif ( scalar(@pepgencoords) == 2 ) {
		#throw("WOW got a 2 feature stop codon for $id strand " . $translateable[0]->strand);
		warning("WOW got a 2 feature stop codon for $id strand " . $translateable[0]->strand);
		return (());
	}
	
	unless ( $pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
		#throw("pep end maps to gap\n");
		warning("pep end maps to gap\n");
		return (());
	}
	unless ( $pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
		#throw("pep end (end of) maps to gap\n");
		warning("pep end (end of) maps to gap\n");
		return (());
	}

	my (@stopc_feat);
	my ($phase) = 0;
	foreach my $pepgencoord ( @pepgencoords ) {
		push @stopc_feat, new Bio::EnsEMBL::SeqFeature(
					                             -seqname => $id,
					                             -source_tag => 'endtrans',
					                             -primary_tag => 'similarity',
					                             -start => $pepgencoord->start,
					                             -end   => $pepgencoord->end,
					                             -phase => $phase,
					                             -strand => $translateable[0]->strand);
		$phase = 3 - ($pepgencoord->end-$pepgencoord->start+1);
	}
	
	if ( $translateable[0]->strand == 1 ) {
		@stopc_feat = sort {$a->start <=> $b->start } @stopc_feat;
	}
	else {
		@stopc_feat = sort {$b->start <=> $a->start } @stopc_feat;
	}
	return @stopc_feat;
}

sub _h_do_escape {
	my $data = shift;
	
	my $equal_esc     = uri_escape("=");
	my $comma_esc     = uri_escape(",");
	my $semicolon_esc = uri_escape(";");
	
	$data =~ s/,/$comma_esc/g;
	$data =~ s/;/$semicolon_esc/g;
	$data =~ s/=/$equal_esc/g;
	
	return $data;
}

sub _transcript_remarks {
	my ($transcript) = @_;

	#my (@notes) = ();
	my ($remarks) = '';

	foreach my $att ( @{$transcript->get_all_Attributes()} ) {
		my ($val)  = _h_do_escape($att->value);
		my ($code) = $att->code;
		if ( ($code =~ /cds_end_NF|cds_start_NF|mRNA_end_NF|mRNA_start_NF/) and ($val != 0) ) {
			#push(@notes, " tag \"$code\";");
			$remarks .= "$code|";
		}
		elsif ( ($code eq 'remark') and ($val =~ /non.ATG/) ) {
			#push(@notes, " tag \"non_ATG_start\";");
			$remarks .= "non_ATG_start|";
		}
		elsif ( ($code eq "remark") and ($val =~/^(alternative 3' UTR|alternative 5' UTR|readthrough|NMD exception|[Nn]o.+org.+sup.+ted|not best-in-genome evidence|non-submitted evidence|upstream ATG|downstream ATG|upstream uORF|overlapping uORF|NAGNAG splice site|non canonical conserved|non canonical genome sequence error|non canonical other|non canonical polymorphism|non canonical U12|non canonical TEC)$/) ){
			$val =~ s/ /_/g;
			$val =~ s/'//g;
			$val =~ s/-/_/g;
			$val =~ s/^[Nn]o.+org.+sup.+ted$/not_organism_supported/;
			#push(@notes, " tag \"$val\";");
			$remarks .= "$val|";
		}
	}
	#return (\@notes);
	return $remarks;
}

sub _check_start_and_stop {
	my ($trans) = @_;
	
	return (0,0) if (!defined($trans->translation));
	
	my ($coding_start) = $trans->cdna_coding_start;
	my ($coding_end)   = $trans->cdna_coding_end;
	my ($cdna_seq)     = uc($trans->spliced_seq);
	my ($startseq)     = substr($cdna_seq,$coding_start-1,3);
	my ($endseq)       = substr($cdna_seq,$coding_end-3,3);
	my ($has_start) = 1;
	my ($has_end) = 1;
	$has_start = 0  if ($startseq ne "ATG");
	$has_end = 0 if (($endseq ne "TAG") && ($endseq ne "TGA") && ($endseq ne "TAA"));
	
	return ($has_start, $has_end);
}

sub DESTROY {}

1;