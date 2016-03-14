--
-- Table structure for table datasource
--
CREATE TABLE datasource (
  datasource_id INT(11) unsigned NOT NULL, 
  name VARCHAR(50) NOT NULL,
  description TEXT DEFAULT NULL,
  url TEXT DEFAULT NULL,
  PRIMARY KEY  (datasource_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table entity
--
CREATE TABLE entity (
  entity_id INT(11) unsigned NOT NULL auto_increment,
  datasource_id INT(11) unsigned NOT NULL,
  identifier VARCHAR(50) DEFAULT NULL,
  source ENUM('GENCODE','HAVANA','ENSEMBL') DEFAULT NULL,
  biotype VARCHAR(50) DEFAULT NULL,
  status VARCHAR(50) DEFAULT NULL,
  level INT(1) DEFAULT NULL,
  version INT(5) DEFAULT NULL,
  CONSTRAINT fk_entity_datasource FOREIGN KEY (datasource_id) REFERENCES datasource (datasource_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (entity_id),
  UNIQUE KEY unique_key_entity_entity_id (entity_id),
  KEY key_entity_identifier (identifier),
  INDEX index_xref_entity_id_identifier (entity_id,identifier)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table xref_identify
--
CREATE TABLE xref_identify (
  entity_id INT(11) unsigned NOT NULL,
  datasource_id INT(11) unsigned NOT NULL,
  identifier VARCHAR(50) NOT NULL,
  CONSTRAINT fk_identify_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT fk_identify_datasource FOREIGN KEY (datasource_id) REFERENCES datasource (datasource_id) ON DELETE CASCADE ON UPDATE CASCADE,
  KEY key_identify_entity_id (entity_id),
  KEY key_identify_datasource_id (datasource_id),
  KEY key_identify_identifier (identifier),
  INDEX index_xref_identify_entity_id_datasource_id (entity_id,datasource_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table coordinate
--
CREATE TABLE coordinate (
  entity_id INT(11) unsigned NOT NULL,
  chromosome VARCHAR(25) NOT NULL,
  start INT(20) unsigned NOT NULL,
  end INT(20) unsigned NOT NULL,
  strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value)
  CONSTRAINT fk_coordinate_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  KEY unique_key_coordinate_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table type
--
CREATE TABLE type (
  type_id INT(11) unsigned NOT NULL,
  name VARCHAR(50) NOT NULL,
  description TEXT DEFAULT NULL,
  PRIMARY KEY  (type_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table sequence
--
CREATE TABLE sequence (
  entity_id INT(11) unsigned NOT NULL,
  type_id INT(11) unsigned NOT NULL,
  length INT(10) unsigned DEFAULT NULL,
  sequence TEXT NOT NULL,
  CONSTRAINT fk_sequence_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT fk_sequence_type FOREIGN KEY (type_id) REFERENCES type (type_id) ON DELETE CASCADE ON UPDATE CASCADE,
  KEY key_sequence_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table alignment
--
CREATE TABLE alignment (
  entity_id INT(11) unsigned NOT NULL,
  type_id INT(11) unsigned NOT NULL,
  num_species INT(10) unsigned DEFAULT NULL,
  length INT(10) unsigned DEFAULT NULL,
  alignment TEXT NOT NULL,
  CONSTRAINT fk_alignment_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT fk_alignment_type FOREIGN KEY (type_id) REFERENCES type (type_id) ON DELETE CASCADE ON UPDATE CASCADE,
  KEY key_alignment_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table exon
--
CREATE TABLE exon (
  entity_id INT(11) unsigned NOT NULL,
  exon_id INT(11) unsigned NOT NULL,
  identifier VARCHAR(50) DEFAULT NULL,
  start INT(10) unsigned NOT NULL,
  end INT(10) unsigned NOT NULL,
  strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value) 
  CONSTRAINT fk_exon_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  KEY key_exon_id (exon_id),
  KEY key_exon_identifier (identifier),  
  KEY key_exon_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table cds
--
CREATE TABLE cds (
  entity_id INT(11) unsigned NOT NULL,
  cds_id INT(11) unsigned NOT NULL,
  start INT(10) unsigned NOT NULL,
  end INT(10) unsigned NOT NULL,
  strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value)
  phase ENUM('.','0','1','2') NOT NULL,  -- default value '.' (the first value)
  CONSTRAINT fk_cds_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  KEY key_cds_id (cds_id),
  KEY key_cds_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table codon
--
CREATE TABLE codon (
  entity_id INT(11) unsigned NOT NULL,
  type ENUM('start','stop') NOT NULL,
  start INT(10) unsigned NOT NULL,
  end INT(10) unsigned NOT NULL,
  strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value)
  phase ENUM('.','0','1','2') NOT NULL,  -- default value '.' (the first value)
  CONSTRAINT fk_codon_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  KEY key_codon_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table crash
--
CREATE TABLE crash (
  crash_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  result LONGTEXT DEFAULT NULL,
  sp_score INT(3) DEFAULT NULL,
  tp_score INT(3) DEFAULT NULL,
  peptide_signal ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  mitochondrial_signal ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  CONSTRAINT fk_crash_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (crash_id),
  UNIQUE KEY unique_key_crash_crash_id (crash_id),
  UNIQUE KEY key_crash_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table crash_residues
-- Save aminoacid residue relative to transcript coordinates
--
CREATE TABLE crash_residues (
  crash_residues_id INT(11) unsigned NOT NULL auto_increment,
  crash_id INT(11) unsigned NOT NULL,
  s_mean FLOAT(3,3) NOT NULL,
  s_prob FLOAT(3,3) NOT NULL,
  d_score FLOAT(3,3) NOT NULL,
  c_max FLOAT(3,3) NOT NULL,
  reliability INT(1) NOT NULL,
  localization CHAR(1) NOT NULL,
  start INT(11) unsigned NOT NULL,
  end INT(11) unsigned NOT NULL,  
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL,
  CONSTRAINT fk_crash_residues_crash FOREIGN KEY (crash_id) REFERENCES crash (crash_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (crash_residues_id),
  UNIQUE KEY unique_key_crash_residues_crash_residues_id (crash_residues_id),
  INDEX index_crash_residues_crash_residues_id_crash_id (crash_residues_id,crash_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table firestar
--
CREATE TABLE firestar (
  firestar_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  num_residues INT(5) DEFAULT NULL,
  result LONGTEXT DEFAULT NULL,
  functional_residue ENUM('ACCEPT','REJECT') DEFAULT NULL,
  CONSTRAINT fk_firestar_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (firestar_id),
  UNIQUE KEY unique_key_firestar_firestar_id (firestar_id),
  UNIQUE KEY key_firestar_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table firestar_residues
-- Save aminoacid residue relative to transcript coordinates
--
CREATE TABLE firestar_residues (
  firestar_residues_id INT(11) unsigned NOT NULL auto_increment,
  firestar_id INT(11) unsigned NOT NULL,
  peptide_position INT(11) unsigned NOT NULL,
  domain VARCHAR(13) DEFAULT NULL,
  ligands TEXT DEFAULT NULL,
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL,  
  CONSTRAINT fk_firestar_residues_firestar FOREIGN KEY (firestar_id) REFERENCES firestar (firestar_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (firestar_residues_id),
  UNIQUE KEY unique_key_firestar_residues_firestar_residues_id (firestar_residues_id),
  INDEX index_firestar_residues_firestar_residues_id_firestar_id (firestar_residues_id,firestar_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table matador3d
--
CREATE TABLE matador3d (
  matador3d_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  result LONGTEXT DEFAULT NULL,
  score FLOAT DEFAULT NULL,
  conservation_structure ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  CONSTRAINT fk_matador3d_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (matador3d_id),
  UNIQUE KEY unique_key_matador3d_matador3d_id (matador3d_id),
  UNIQUE KEY key_matador3d_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table matador3d_alignments
--
CREATE TABLE matador3d_alignments (
  matador3d_alignments_id INT(11) unsigned NOT NULL auto_increment,
  matador3d_id INT(11) unsigned NOT NULL,
  alignment_start INT(20) unsigned DEFAULT NULL,
  alignment_end INT(20) unsigned DEFAULT NULL,  
  pdb_id TEXT DEFAULT NULL,
  identity TEXT DEFAULT NULL,
  external_id TEXT DEFAULT NULL,
  cds_id INT(2) unsigned NOT NULL,  
  start INT(11) unsigned NOT NULL,
  end INT(11) unsigned NOT NULL,
  score FLOAT NOT NULL,
  type TEXT DEFAULT NULL,  
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value)
  CONSTRAINT fk_matador3d_alignments_matador3d FOREIGN KEY (matador3d_id) REFERENCES matador3d (matador3d_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (matador3d_alignments_id),
  UNIQUE KEY unique_key_matador3d_alignments_matador3d_alignments_id (matador3d_alignments_id),
  INDEX index_matador3d_alignments_matador3d_alignments_id_matador3d_id (matador3d_alignments_id,matador3d_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table thump
--
CREATE TABLE thump (
  thump_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  result LONGTEXT DEFAULT NULL,
  num_tmh INT(5) DEFAULT NULL,
  num_damaged_tmh INT(5) DEFAULT NULL,
  transmembrane_signal ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  CONSTRAINT fk_thump_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (thump_id),
  UNIQUE KEY unique_key_thump_thump_id (thump_id),
  UNIQUE KEY key_thump_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table thump_helixes
-- Save helix coordinates
--
CREATE TABLE thump_helixes (
  thump_helixes_id INT(11) unsigned NOT NULL auto_increment,
  thump_id INT(11) unsigned NOT NULL,
  start INT(20) unsigned NOT NULL,
  end INT(20) unsigned NOT NULL,
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value)  
  damaged TINYINT(1) DEFAULT 0,
  CONSTRAINT fk_thump_helixes_thump FOREIGN KEY (thump_id) REFERENCES thump (thump_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (thump_helixes_id),
  UNIQUE KEY unique_key_thump_helixes_thump_helixes_id (thump_helixes_id),
  INDEX index_thump_helixes_thump_helixes_id_thump_id (thump_helixes_id,thump_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table spade 
--
CREATE TABLE spade (
  spade_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  result LONGTEXT DEFAULT NULL,
  num_domains INT(5) DEFAULT NULL,
  num_possibly_damaged_domains INT(5) DEFAULT NULL,
  num_damaged_domains INT(5) DEFAULT NULL,
  num_wrong_domains INT(5) DEFAULT NULL,
  bitscore FLOAT NOT NULL,
  domain_signal ENUM('NO','YES','UNKNOWN') DEFAULT NULL,  
  CONSTRAINT fk_spade_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (spade_id),
  UNIQUE KEY unique_key_spade_spade_id (spade_id),
  UNIQUE KEY key_spade_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table spade_alignments
-- Save alignments, envelopes, and hmm's for each transcript
--
CREATE TABLE spade_alignments (
  spade_alignments_id INT(11) unsigned NOT NULL auto_increment,
  spade_id INT(11) unsigned NOT NULL,
  alignment_start INT(20) unsigned NOT NULL,
  alignment_end INT(20) unsigned NOT NULL,
  envelope_start INT(20) unsigned NOT NULL,
  envelope_end INT(20) unsigned NOT NULL,
  hmm_start INT(20) unsigned NOT NULL,
  hmm_end INT(20) unsigned NOT NULL,
  hmm_length INT(10) unsigned NOT NULL,
  hmm_acc VARCHAR(50) NOT NULL,
  hmm_name VARCHAR(50) NOT NULL,
  hmm_type VARCHAR(50) NOT NULL,
  bit_score FLOAT NOT NULL,
  evalue DOUBLE NOT NULL,
  significance VARCHAR(50) DEFAULT NULL,
  clan VARCHAR(50) DEFAULT NULL,
  predicted_active_site_residues TEXT DEFAULT NULL,
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value)  
  type_domain ENUM('domain','domain_possibly_damaged','domain_damaged','domain_wrong') DEFAULT NULL,  
  external_id VARCHAR(50) DEFAULT NULL,
  discarded INT(1) DEFAULT NULL,
  CONSTRAINT fk_spade_spade_alignments FOREIGN KEY (spade_id) REFERENCES spade (spade_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (spade_alignments_id),
  UNIQUE KEY unique_key_spade_alignments_spade_alignments_id (spade_alignments_id),
  INDEX index_spade_alignments_spade_alignments_id_spade_id (spade_alignments_id,spade_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table corsair
--
CREATE TABLE corsair (
  corsair_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  result LONGTEXT DEFAULT NULL,  
  score  FLOAT DEFAULT NULL,  
  vertebrate_signal ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  CONSTRAINT fk_corsair_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (corsair_id),
  UNIQUE KEY unique_key_corsair_corsair_id (corsair_id),
  KEY key_corsair_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table corsair_alignments
--
CREATE TABLE corsair_alignments (
  corsair_alignments_id INT(11) unsigned NOT NULL auto_increment,
  corsair_id INT(11) unsigned NOT NULL,
  cds_id INT(2) unsigned NOT NULL,  
  start INT(11) unsigned NOT NULL,
  end INT(11) unsigned NOT NULL,
  score FLOAT NOT NULL,
  maxscore FLOAT DEFAULT NULL,
  sp_report LONGTEXT DEFAULT NULL,
  type TEXT DEFAULT NULL,  
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL, -- default value '.' (the first value)
  CONSTRAINT fk_corsair_alignments_corsair FOREIGN KEY (corsair_id) REFERENCES corsair (corsair_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (corsair_alignments_id),
  UNIQUE KEY unique_key_corsair_alignments_corsair_alignments_id (corsair_alignments_id),
  INDEX index_corsair_alignments_corsair_alignments_id_corsair_id (corsair_alignments_id,corsair_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table inertia
--
CREATE TABLE inertia (
  inertia_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  result LONGTEXT DEFAULT NULL,  
  unusual_evolution ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  CONSTRAINT fk_inertia_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (inertia_id),
  UNIQUE KEY unique_key_inertia_inertia_id (inertia_id),
  KEY key_inertia_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table inertia_residues
--
CREATE TABLE inertia_residues (
  inertia_residues_id INT(11) unsigned NOT NULL auto_increment,
  inertia_id INT(11) unsigned NOT NULL,
  inertia_exon_id INT(11) unsigned NOT NULL,
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL,  
  unusual_evolution ENUM('NO','YES','UNKNOWN') DEFAULT NULL,  
  CONSTRAINT fk_inertia_residues_inertia FOREIGN KEY (inertia_id) REFERENCES inertia (inertia_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (inertia_residues_id),
  UNIQUE KEY unique_key_inertia_residues_inertia_residues_id (inertia_residues_id),
  KEY key_inertia_residues_inertia_id (inertia_id),
  KEY key_inertia_residues_exon_id (inertia_exon_id),
  INDEX index_inertia_residues_id_inertia_id_inertia_exon_id (inertia_residues_id,inertia_id,inertia_exon_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table slr_type
--
CREATE TABLE slr_type (
  slr_type_id INT(11) unsigned NOT NULL,
  name VARCHAR(50) NOT NULL,
  description TEXT DEFAULT NULL,
  PRIMARY KEY  (slr_type_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table omega
--
CREATE TABLE omega (
  omega_id INT(11) unsigned NOT NULL auto_increment,
  inertia_id INT(11) unsigned NOT NULL,
  slr_type_id INT(11) unsigned NOT NULL,
  omega_average FLOAT DEFAULT NULL,
  omega_st_desviation FLOAT DEFAULT NULL, 
  result LONGTEXT DEFAULT NULL,
  unusual_evolution ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  CONSTRAINT fk_omega_inertia FOREIGN KEY (inertia_id) REFERENCES inertia (inertia_id) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT fk_omega_slr_type FOREIGN KEY (slr_type_id) REFERENCES slr_type (slr_type_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (omega_id),
  UNIQUE KEY unique_key_omega_omega_id (omega_id),
  KEY key_omega_inertia_id (inertia_id),
  KEY key_omega_slr_type_id (slr_type_id),
  INDEX index_omega_id_inertia_id_slr_type_id (omega_id,inertia_id,slr_type_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table omega_residues
--
CREATE TABLE omega_residues (
  omega_residues_id INT(11) unsigned NOT NULL auto_increment,
  omega_id INT(11) unsigned NOT NULL,
  omega_exon_id INT(11) unsigned NOT NULL,
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  omega_mean DOUBLE NOT NULL,
  st_deviation DOUBLE NOT NULL,
  p_value DOUBLE NOT NULL,
  difference_value DOUBLE NOT NULL,
  unusual_evolution ENUM('NO','YES','UNKNOWN') DEFAULT NULL,  
  CONSTRAINT fk_omega_residues_omega FOREIGN KEY (omega_id) REFERENCES omega (omega_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (omega_residues_id),
  UNIQUE KEY unique_key_omega_residues_omega_residues_id (omega_residues_id),
  KEY key_omega_residues_omega_id (omega_id),
  KEY key_omega_residues_exon_id (omega_exon_id),
  INDEX index_omega_residues_id_omega_id_omega_exon_id (omega_residues_id,omega_id,omega_exon_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table proteo
--
CREATE TABLE proteo (
  proteo_id INT(11) unsigned NOT NULL auto_increment,
  entity_id INT(11) unsigned NOT NULL,
  num_peptides INT(5) DEFAULT NULL,
  result LONGTEXT DEFAULT NULL,
  peptide_evidence ENUM('NO','YES','UNKNOWN') DEFAULT NULL,
  CONSTRAINT fk_proteo_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (proteo_id),
  UNIQUE KEY unique_key_proteo_proteo_id (proteo_id),
  UNIQUE KEY key_proteo_entity_id (entity_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table proteo_peptides
-- Save peptide evidences relative to transcript coordinates
--
CREATE TABLE proteo_peptides (
  proteo_peptides_id INT(11) unsigned NOT NULL auto_increment,
  proteo_id INT(11) unsigned NOT NULL,
  peptide_id VARCHAR(50) DEFAULT NULL,
  sequence TEXT NOT NULL,
  num_experiments INT(5) unsigned NOT NULL,
  experiments VARCHAR(13) DEFAULT NULL,
  start INT(20) unsigned NOT NULL,
  end INT(20) unsigned NOT NULL,  
  trans_start INT(20) unsigned NOT NULL,
  trans_end INT(20) unsigned NOT NULL,
  trans_strand ENUM('.','+','-') NOT NULL,  
  CONSTRAINT fk_proteo_peptides_proteo FOREIGN KEY (proteo_id) REFERENCES proteo (proteo_id) ON DELETE CASCADE ON UPDATE CASCADE,
  PRIMARY KEY (proteo_peptides_id),
  UNIQUE KEY unique_key_proteo_peptides_proteo_peptides_id (proteo_peptides_id),
  INDEX index_proteo_peptides_proteo_peptides_id_proteo_id (proteo_peptides_id,proteo_id)
) ENGINE=InnoDB CHARSET=utf8;

--
-- Table structure for table appris
--
CREATE TABLE appris (
  entity_id INT(11) unsigned NOT NULL,  
  functional_residues_score VARCHAR(10) DEFAULT NULL, -- firestar
  homologous_structure_score VARCHAR(10) DEFAULT NULL, -- Matador3D
  vertebrate_conservation_score VARCHAR(10) DEFAULT NULL, -- CORSAIR
  domain_score VARCHAR(10) DEFAULT NULL, -- SPADE
  transmembrane_helices_score VARCHAR(10) DEFAULT NULL, -- THUMP
  peptide_score VARCHAR(10) DEFAULT NULL, -- CRASH - SignalP
  mitochondrial_score VARCHAR(10) DEFAULT NULL, -- CRASH - TargetP  
  unusual_evolution_score VARCHAR(10) DEFAULT NULL, -- INERTIA
  peptide_evidence_score VARCHAR(10) DEFAULT NULL, -- CRASH - PROTEO
  principal_isoform_score VARCHAR(10) DEFAULT NULL, -- APPRIS  
  functional_residues_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- firestar
  homologous_structure_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- Matador3D
  vertebrate_conservation_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- CORSAIR
  domain_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- SPADE
  transmembrane_helices_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- THUMP
  peptide_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- CRASH - SignalP
  mitochondrial_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- CRASH - TargetP
  unusual_evolution_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- INERTIA
  peptide_evidence_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- PROTEO
  principal_isoform_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-', -- APPRIS
  reliability VARCHAR(15) DEFAULT NULL,
  result LONGTEXT DEFAULT NULL,
  CONSTRAINT fk_appris_entity FOREIGN KEY (entity_id) REFERENCES entity (entity_id) ON DELETE CASCADE ON UPDATE CASCADE,
  UNIQUE KEY key_appris_entity_id (entity_id),
  INDEX index_appris_pi (principal_isoform_signal)
) ENGINE=InnoDB CHARSET=utf8;

-- COMMENT: ADD PROTEO COLUMNS 
-- ALTER TABLE appris ADD COLUMN peptide_evidence_signal ENUM('-','NO','YES','UNKNOWN') DEFAULT '-';
-- ALTER TABLE appris ADD COLUMN peptide_evidence_rscore CHAR(1) DEFAULT NULL;
-- ALTER TABLE appris ADD COLUMN peptide_evidence_score VARCHAR(10) DEFAULT NULL;



--
-- Insert default value into datasource table
--
INSERT INTO datasource SET 
  datasource_id='1',
  name='Havana_Gene_Id',
  description='Havana Gene Identifier',
  url='http://www.sanger.ac.uk';

INSERT INTO datasource SET
  datasource_id='2',
  name='Havana_Transcript_Id',
  description='Havana Transcript Identifier',
  url='http://www.sanger.ac.uk';

INSERT INTO datasource SET
  datasource_id='3',
  name='Havana_Peptide_Id',
  description='Havana Peptide Identifier',
  url='http://www.sanger.ac.uk';

INSERT INTO datasource SET
  datasource_id='4',
  name='External_Id',
  description='External Identifier',
  url='http://www.ensembl.org/index.html';

INSERT INTO datasource SET
  datasource_id='5',
  name='Ensembl_Gene_Id',
  description='Ensembl Gene Identifier',
  url='http://www.ensembl.org/index.html';

INSERT INTO datasource SET
  datasource_id='6',
  name='Ensembl_Transcript_Id',
  description='Ensembl Transcript Identifier',
  url='http://www.ensembl.org/index.html';

INSERT INTO datasource SET
  datasource_id='7',
  name='Ensembl_Peptide_Id',
  description='Ensembl Peptide Identifier',
  url='http://www.ensembl.org/index.html';

INSERT INTO datasource SET
  datasource_id='8',
  name='UniProtKB_SwissProt',
  description='UniProtKB/Swiss-Prot Accession',
  url='http://www.uniprot.org';

INSERT INTO datasource SET
  datasource_id='9',
  name='CCDS',
  description='Consensus CDS',
  url='http://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi';

INSERT INTO datasource SET
  datasource_id='10',
  name='TSL',
  description='Transcript Support Level',
  url='http://www.ensembl.org/Help/Glossary?id=492';

--
-- Insert default value into type table
--
INSERT INTO type SET 
  type_id='1',
  name='transcript',
  description='Nucleotide sequence corresponds to transcript';

INSERT INTO type SET 
  type_id='2',
  name='peptide',
  description='Aminoacid sequence corresponds to peptide';

INSERT INTO type SET  
  type_id='3',
  name='cds',
  description='Alignment of CDS (coming from Maf)';

INSERT INTO type SET  
  type_id='4',
  name='filtered_cds',
  description='Alignment of CDS that is filtered from the size of sequence (coming from Maf)';


--
-- Insert default value into slr_type table
--
INSERT INTO slr_type SET 
  slr_type_id='1',
  name='filter',
  description='Filtered Alignment coming directly from Maf';

INSERT INTO slr_type SET 
  slr_type_id='2',
  name='prank',
  description='Filtered Alignment coming from Prank software';

INSERT INTO slr_type SET 
  slr_type_id='3',
  name='kalign',
  description='Filtered Alignment coming from Kalign software';
