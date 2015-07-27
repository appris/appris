=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Utils::EnsEMBL - Utility functions for protein handling

=head1 SYNOPSIS

  use APPRIS::Utils::EnsEMBL
    qw(
       get_id
     );

  or to get all methods just

  use APPRIS::Utils::EnsEMBL;

  eval { get_id(cds_list) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

The functions exported by this package provide a set of useful methods 
to handle files.

=head1 METHODS

=cut

package APPRIS::Utils::EnsEMBL;

use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	get_xref_identifiers
	get_ccds_id
	get_exons
);


###################
# Global variable #
###################
use APPRIS::Utils::Constant qw(
	$ENSEMBL_VERSION
	$ENSEMBL_HOST
	$ENSEMBL_USER
	$ENSEMBL_VERBOSE
	$ENSEMBL_SPECIES
);

#my ($ENSEMBL_VERSION) = $APPRIS::Utils::Constant::ENSEMBL_VERSION;
my ($ENSEMBL_VERSION) = $Bio::EnsEMBL::ApiVersion::software_version;
my ($ENSEMBL_HOST)    = $APPRIS::Utils::Constant::ENSEMBL_HOST;
my ($ENSEMBL_USER)    = $APPRIS::Utils::Constant::ENSEMBL_USER;
my ($ENSEMBL_VERBOSE) = $APPRIS::Utils::Constant::ENSEMBL_VERBOSE;
my ($ENSEMBL_SPECIES) = $APPRIS::Utils::Constant::ENSEMBL_SPECIES;
			

#####################
# Method prototypes #
#####################
sub get_xref_identifiers($$);
sub get_ccds_id($$);
sub get_exons($$);

=head2 get_xref_identifiers

  Arg [1]    : String Stable identifier (Havana, Ensembl) $id
  Example    : use APPRIS::Utils::EnsEMBL qw(get_id);
               get_xref_identifiers($id);
  Description: Get xref identifiers from id.
  Returntype : Hash of identifiers or undef
  Exceptions : none

=cut

sub get_xref_identifiers($$)
{
	my ($registry, $id) = @_;
	my ($report);
		
	# If Ensembl database is not defined
	unless (defined $registry ) {
		$registry = 'Bio::EnsEMBL::Registry';
		eval {
			$registry->load_registry_from_db(
					-db_version	=> $ENSEMBL_VERSION,
					-host		=> $ENSEMBL_HOST,
					-user		=> $ENSEMBL_USER,
					-species	=> $ENSEMBL_SPECIES,
					-verbose	=> $ENSEMBL_VERBOSE,
			);
		};
		return undef if $@;		
	}
		
	# Get Ensembl adaptors
	my ($gene, $transcript);	
	eval {
		if ( $id =~ /^OTTHUMG/ )
		{
			my ($gene_adaptor) = $registry->get_adaptor('human','vega','gene');
			my ($gene_aux) = @{$gene_adaptor->fetch_all_by_external_name($id)}; # ListRef
			$gene = $gene_aux; # We get the first one
		}
		elsif ( $id =~ /^ENSG/ )
		{
			my ($gene_adaptor) = $registry->get_adaptor('human','core','gene');
			$gene = $gene_adaptor->fetch_by_stable_id($id);
		}
		elsif ( $id =~ /^OTTHUMT/ )
		{
			my ($transcript_adaptor) = $registry->get_adaptor('human','vega','transcript');
			my ($transcript_aux) = @{$transcript_adaptor->fetch_all_by_external_name($id)}; # ListRef
			$transcript = $transcript_aux; # We get the first one
		}
		elsif ( $id=~/^ENST/ )
		{
			my ($transcript_adaptor) = $registry->get_adaptor('human','core','transcript');
			$transcript = $transcript_adaptor->fetch_by_stable_id($id);
		}
	};
	return undef if $@;
		
	# When is a gene object
 	if ( defined $gene and defined $gene->stable_id() )
 	{
 		if ( defined $gene->description() ) {	 			
 			my ($description) = $gene->description();
 			$report->{'description'} = $description;
 			if ( $description =~ /\[Source:UniProtKB\/Swiss-Prot;Acc:([^\]]*)\]/ ) {
 				my ($uniprot_id) = $1;
 				$report->{'uniprot_id'} = $uniprot_id if(defined $uniprot_id);	
 			}
 		}
		my(@db_entries) = @{$gene->get_all_DBEntries};
		foreach my $dbe (@db_entries)
		{
			if(defined $dbe->dbname() and defined $dbe->display_id())
			{
				if(($dbe->dbname() eq 'Ens_Hs_gene') and defined $dbe->display_id())
				{
					$report->{'ensembl_id'} = $dbe->display_id();
				}
				elsif(($dbe->dbname() eq 'HGNC') and defined $dbe->display_id())
				{
					$report->{'hgnc'} = $dbe->display_id();
				}
				elsif(($dbe->dbname() eq 'Uniprot/SWISSPROT') and defined $dbe->display_id())
				{
					$report->{'uniprot_id'} = $dbe->display_id();
				}
			}
		}
 	}
 	# When is a transcript object	
	elsif ( defined $transcript and defined $transcript->stable_id() )
	{
 		# Get peptide id
		if(defined $transcript->translation() and defined $transcript->translation()->stable_id())
		{
			my ($translation) = $transcript->translation();
			$report->{'peptide_id'} = $translation->stable_id();
		}	
	}
		
	return $report;
}

=head2 get_ccds_id

  Arg [1]    : String Transcript identifier $id
  Example    : use APPRIS::Utils::EnsEMBL qw(get_ccds_id);
               get_ccds_id($id);
  Description: Get the CCDS identifier from stable transcript id.
  Returntype : String or undef
  Exceptions : none

=cut

sub get_ccds_id($$)
{
	my ($registry, $id) = @_;
	my ($ccdsIDs);

	# If Ensembl database is not defined
	unless (defined $registry ) {
		$registry = 'Bio::EnsEMBL::Registry';
		eval {
			$registry->load_registry_from_db(
					-db_version	=> $ENSEMBL_VERSION,
					-host		=> $ENSEMBL_HOST,
					-user		=> $ENSEMBL_USER,
					-species	=> $ENSEMBL_SPECIES,
					-verbose	=> $ENSEMBL_VERBOSE,
			);
		};
		return undef if $@;		
	}
		
	eval {	
		my ($transcript_adaptor) = $registry->get_adaptor('human','core','transcript');
		my ($transcript) = $transcript_adaptor->fetch_by_stable_id($id);
		if($transcript) {
			my(@dblinks) = @{$transcript->get_all_DBLinks};
			
			foreach my $entry (@dblinks) {
				if ($entry->dbname eq 'CCDS') {
					push(@{$ccdsIDs}, $entry->primary_id());
				}
			}		
		}
	};
	return undef if $@;
		
	return $ccdsIDs;
}

=head2 get_exons

  Arg [1]    : String Transcript identifier $id
  Example    : use APPRIS::Utils::EnsEMBL qw(get_exons);
               get_exons($id);
  Description: Get the sorted list of exons from stable transcript id.
  Returntype : String or undef
  Exceptions : none

=cut

sub get_exons($$)
{
	my ($registry, $id) = @_;
	my ($exons);
	
	# If Ensembl database is not defined
	unless (defined $registry ) {
		$registry = 'Bio::EnsEMBL::Registry';
		eval {
			$registry->load_registry_from_db(
					-db_version	=> $ENSEMBL_VERSION,
					-host		=> $ENSEMBL_HOST,
					-user		=> $ENSEMBL_USER,
					-species	=> $ENSEMBL_SPECIES,
					-verbose	=> $ENSEMBL_VERBOSE,
			);
		};
		return undef if $@;		
	}
	
	eval {
		$id =~ s/\.\d*$//;		
		my ($transcript_adaptor) = $registry->get_adaptor('human','core','transcript');		
		my ($transcript) = $transcript_adaptor->fetch_by_stable_id($id);
		if($transcript) {
			my(@e_exons) = @{$transcript->get_all_Exons};
			foreach my $exon (@e_exons) {
				my ($o_exon) = {
					'id'		=> $exon->stable_id(),
					'start'		=> $exon->start(),
					'end'		=> $exon->end(),
					'strand'	=> $exon->strand(),
				};
				push(@{$exons}, $o_exon);
			}		
		}
	};
	return undef if $@;
	
	return $exons;
}

1;