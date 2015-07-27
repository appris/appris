=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Utils::Constant - Utility constants

=head1 SYNOPSIS

  use APPRIS::Utils::Constant
    qw(
       get_protein_cds_sequence
     );

  or to get all methods just

  use APPRIS::Utils::Constant;

  eval { get_protein_cds_sequence(cds_list) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

The functions exported by this package provide a set of useful methods 
to handle files.

=head1 METHODS

=cut

package APPRIS::Utils::Constant;

use strict;
use warnings;

use vars qw(
	$API_VERSION
	$VERSION
	$DATE
	$PI_METHOD
	$UNKNOWN_LABEL
	$OK_LABEL
	$NO_LABEL
	$PI_LABEL	
	$APPRIS_PRINC_LABEL
	$APPRIS_CANDI_LABEL
	$APPRIS_CANDI_LONG_LABEL
	$APPRIS_ALTER_LABEL
	$SIGNALP_METHOD
	$SP_TYPE
	$SP_FORMAT
	$SP_METHOD
	$SP_TRUNCATE
	$TARGETP_METHOD
	$FIRESTAR_METHOD
	$FIRESTAR_ACCEPT_LABEL
	$FIRESTAR_REJECT_LABEL
	$CEXONIC_METHOD
	$MATADOR3D_METHOD
	$THUMP_METHOD
	$SPADE_METHOD
	$CORSAIR_METHOD	
	$OMEGA_METHOD
	$OMEGA_THRESHOLD
	$OMEGA_D_VALUE_THRESHOLD
	$OMEGA_P_VALUE_THRESHOLD
	$INERTIA_METHOD
	$MOBY_NAMESPACE
	$MOBY_CENTRAL_URL
	$MOBY_CENTRAL_URI
	$HAVANA_SOURCE
	$ENSEMBL_SOURCE
	$METHOD_DESC
	$METHOD_SCORE_DESC
	$METHOD_LABEL_DESC
	$ENSEMBL_VERSION
	$ENSEMBL_HOST
	$ENSEMBL_USER
	$ENSEMBL_VERBOSE
	$ENSEMBL_SPECIES
);

# Version variables
$API_VERSION = 'rel15';
$VERSION = 'rel15_v1';
$DATE = '2-May-2013';

# PI annotations
$PI_METHOD = 'PI';
$UNKNOWN_LABEL = 'UNKNOWN';
$OK_LABEL = 'YES';
$NO_LABEL = 'NO';
$PI_LABEL = 'principal_isoform';
$APPRIS_PRINC_LABEL = "appris_principal";
$APPRIS_CANDI_LABEL = "appris_candidate";
$APPRIS_CANDI_LONG_LABEL = "appris_candidate_longest";
$APPRIS_ALTER_LABEL = "appris_alternative";


# Constant for 'SignalP' method
$SIGNALP_METHOD = 'SignalP';
$SP_TYPE = "euk"; 
$SP_FORMAT = "full";
$SP_METHOD = "nn+hmm";
$SP_TRUNCATE = 70;

# Constant for 'TargetP' method
$TARGETP_METHOD = 'TargetP';

# Constant for 'Firestar' method
$FIRESTAR_METHOD = 'Firestar';
$FIRESTAR_ACCEPT_LABEL = 'ACCEPT';
$FIRESTAR_REJECT_LABEL = 'REJECT';

# Constant for 'CExonic' method
$CEXONIC_METHOD = 'CExonic';

# Constant for 'Matador3D' method
$MATADOR3D_METHOD = 'Matador3D';

# Constant for 'THUMP' method
$THUMP_METHOD = 'THUMP';

# Constant for 'SPADE' method
$SPADE_METHOD = 'SPADE';

# Constant for 'CORSAIR' method
$CORSAIR_METHOD = 'CORSAIR';

# Constant for 'Omega' method
$OMEGA_METHOD = 'Omega';
$OMEGA_THRESHOLD = 0.25;
$OMEGA_D_VALUE_THRESHOLD = 0.35;
$OMEGA_P_VALUE_THRESHOLD = 0.025;

# Constant for 'Inertia' method
$INERTIA_METHOD = 'Inertia';

# Source
$HAVANA_SOURCE = 'HAVANA';
$ENSEMBL_SOURCE = 'ENSEMBL';

# Method description
$METHOD_DESC = {
	'appris'	=> 'principal_isoform',
	'firestar'	=> 'functional_residue',
	'matador3d'	=> 'homologous_structure',
	'corsair'	=> 'vertebrate_conservation',
	'spade'		=> 'functional_domain',
	'thump'		=> 'transmembrane_signal',
	'crash'		=> 'peptide_mitochondrial_signal',
	'crash_sp'	=> 'signal_peptide',
	'crash_tp'	=> 'mitochondrial_signal',
	'inertia'	=> 'neutral_evolution',
	'proteo'	=> 'proteomic_evidence',
};
$METHOD_SCORE_DESC = {
	'appris'	=> "Principal Isoform",
	'firestar'	=> "Num. Functional Residues",
	'matador3d'	=> "Tertiary Structure Score",
	'corsair'	=> "Conservation score",
	'spade'		=> "Whole Domains",
	'thump'		=> "Num. Transmembrane Helices",
	'crash'		=> "Peptide / Mitochondrial Signal",
	'crash_sp'	=> "Peptide Signal",
	'crash_tp'	=> "Mitochondrial Signal",
	'inertia'	=> "Unusual Evolution exons",
	'proteo'	=> "Num. Mapping Peptides",
};
$METHOD_LABEL_DESC = {
	'appris'	=> "APPRIS principal isoform",
	'firestar'	=> "Known functional residues",
	'matador3d'	=> "Regions with known 3D structure",
	'corsair'	=> "Isoform with most cross-species evidence",
	'spade'		=> ["Whole Pfam functional domains","Damaged Pfam functional domains", "Whole (and Damaged) Pfam functional domains"],
	'thump'		=> ["Transmembrane helices","Damaged transmembrane helices"],
	'crash'		=> ["Signal peptide sequences","Mitochondrial signal sequences"],
	'inertia'	=> ["Neutrally evolving exons","Unusually evolving exons"],
	'proteo'	=> "CNIO Proteomic Evidence",
};

# Constants (default values) for Ensembl API
$ENSEMBL_VERSION = '70';
$ENSEMBL_HOST    = 'ensembldb.cnio.es';
$ENSEMBL_USER    = 'ensembl';
$ENSEMBL_VERBOSE = 0;
$ENSEMBL_SPECIES = 'Homo sapiens';

1;