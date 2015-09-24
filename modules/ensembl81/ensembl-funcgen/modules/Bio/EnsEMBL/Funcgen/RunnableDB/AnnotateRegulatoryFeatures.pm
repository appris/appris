=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::AnnotateRegulatoryFeatures

=head1 DESCRIPTION

'AnnotateRegulatoryFeatures'
 implements all Damian's regulatory annotation scripts in one single process

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::AnnotateRegulatoryFeatures;

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Annotation');

use constant  NO_ROWS => '0E0';

sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  
  $self->SUPER::fetch_input();

  my $cell_type = $self->param('cell_type') || throw "No cell type defined";
  $self->_cell_type($cell_type);

  return 1;
}

sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift @_;

  #Create a database for this cell_type... maybe add this step to setup?
  my %wparams = %{$self->_workdb_params};
  my $dbparams = "-h".$wparams{'-host'}." -u".$wparams{'-user'}." -P".$wparams{'-port'}." -p".$wparams{'-pass'}; 
  #my $dbname = "annotation_".$self->_species."_".$self->_release."_".$self->_cell_type;
  my $ct = $self->_cell_type;
  $ct =~ s/\-/_/g;
  my $dbname = "annotation_".$self->_efgdba->dbc->dbname."_".$ct;

  system("mysql $dbparams -e 'drop database if exists ".$dbname."; create database ".$dbname."';") && throw "Could not create work database $dbname";
  $wparams{'-dbname'} = $dbname;

  eval{
       $self->_workdb(Bio::EnsEMBL::DBSQL::DBConnection->new( %wparams ));
  };
  if($@) { throw "Error creating the Work DB Connection: $@";  }    
  if(!$self->_workdb()){ throw "Could not create Work DB server connection"; }

  my $workdb = $self->_workdb;

  my $dump_dir = $self->_output_dir."/tmp_".$self->_cell_type."_".int(rand(100));

  open(REPORT,">".$self->_output_dir."/".$self->_cell_type.".report");

  _copy_data($self->_dnadba->dbc,$self->_efgdba->dbc,$workdb,$dump_dir);

  #Perhaps make a feature factory?... Each feature knows which data it needs...
  # ORDER MATTERS!!??
  _create_transcript_features($workdb);
  _create_repeat_features($workdb);
  _create_gene_features($workdb);
  _create_PolIII_features($workdb);
  _create_exon_features($workdb);
  _create_intergenic_features($workdb);
  _cage_ditag_transcript_tss($workdb);

  #Not in use at the moment...
  #&karyotype_features($dbh,$dbu);
  #&cpg_features($dbh,$dbu);

  #Process 
  _process_regulatory_features($self->_dnadba->dbc,$self->_efgdba->dbc,$workdb,$self->_cell_type);

  _generate_random_frags($workdb);

  #read from file or parameters...
  my @gen_feats;
  push @gen_feats,  "protein_coding_exon1_plus_enhancer";
  push @gen_feats,  "protein_coding_transcript_downstream_2500";
  push @gen_feats,  "protein_coding_single_exon_gene_plus_enhancer";
  push @gen_feats,  "protein_coding_intron1";
  push @gen_feats,  "protein_coding_gene_body";
  push @gen_feats,  "intergenic_2500";
  push @gen_feats,  "RNA_gene_single_exon_gene_plus_enhancer";
  # This one needs cage_ditag_transcript_tss
  push @gen_feats,  "tss_centred_5000";
  push @gen_feats,  "PolIII_transcribed_gene_plus_enhancer";

 _calculate_associations($workdb, $self->_output_dir."/".$self->_cell_type, \@gen_feats, $self->_species);

  #&clean_temp();

  close(REPORT);

  return 1;
}


sub write_output {  
  my $self = shift @_;

  return 1;

}

# Dumps data from core and other relevant databases and imports it into a working db
sub _copy_data {
  my ($dnadb, $efgdb, $workdb, $dump_dir) = (shift,shift,shift,shift);

  # Copy data from core
  system("mkdir -p ".$dump_dir) && throw "Error creating $dump_dir folder";

  my $dump_templ = "mysqldump --opt --skip-lock-tables -h ".$dnadb->host.
    " -u ".$dnadb->username.
      " -P ".$dnadb->port.
	" ".$dnadb->dbname.
	  ' %s '.
	    " > ${dump_dir}/".'%s'.".dump";
  
  my $load_templ = "mysql -h ".$workdb->host.
    " -u ".$workdb->username.
      " -P ".$workdb->port.
	" -p".$workdb->password.
	  " ".$workdb->dbname.
	    " < ${dump_dir}/".'%s'.".dump";
  
  my @core_tables = (
		     'transcript',
		     'transcript_stable_id',
		     'gene',
		     'gene_attrib',
		     'attrib_type',
		     'exon',
		     'exon_transcript',
		     'coord_system',
		     'meta_coord',
		     'seq_region',
		     'seq_region_attrib',
		     'assembly',
		     'seq_region_attrib',
		     'simple_feature',
		     'analysis',
		     'karyotype',
		     'xref',
		     'object_xref',
		     'external_db',
		     'translation',
		     'analysis_description'
		    );

  push @core_tables, ('repeat_feature','repeat_consensus');

  foreach my $table (@core_tables){
    
    #unless( $core_dbu->table_exists($table)){die ("ERROR: The core database does not contain table $table")}
    
    my $command = sprintf($dump_templ,$table,$table);
    system($command) && throw "Error running $command";
    $command = sprintf($load_templ,$table);
    system($command) && throw "Error running $command";
  }

  # copy data from funcgen
  $dump_templ = "mysqldump --opt --skip-lock-tables -h ".$efgdb->host.
    " -u ".$efgdb->username.
      " -P ".$efgdb->port;
  if($efgdb->password){
    $dump_templ .= " -p".$efgdb->password;
  }  
  $dump_templ .= " ".$efgdb->dbname.
    ' %s '.
      " > ${dump_dir}/".'%s'.".dump";
  
  my @funcgen_tables =( 'regulatory_attribute',
		     'feature_type',
		     'regulatory_feature',
		     'feature_set',
		     'data_set',
		     'supporting_set',
		     'external_feature',
		     'meta',
		     'annotated_feature'
		   );
  

  foreach my $table (@funcgen_tables){
    
    #Do some checking??
    #unless( $source_dbu->table_exists($table)){die ("ERROR: The funcgen database does not contain table $table")}
    
    my $command = sprintf($dump_templ,$table,$table);
    system($command) && throw "Error running $command";
    $command = sprintf($load_templ,$table);
    system($command) && throw "Error running $command";
  }
  
  system("rm -r -f ".$dump_dir) && throw "Error deleting $dump_dir folder";

}

sub _create_PolIII_features{
    my($dbc) = @_;

    my @sql;
    push @sql, "insert into PolIII_transcribed_gene select * from all_repeats where feature_type = 'PolIII_transcribed_repeat'";
    push @sql, "insert into PolIII_transcribed_gene select * from all_repeats where feature_type = 'tRNA_repeat'";

    push @sql, "drop table if exists PolIII_transcribed_gene_plus_enhancer";
    push @sql, "create table PolIII_transcribed_gene_plus_enhancer select  seq_region_id, seq_region_name,feature_id, 'PolIII_transcribed_gene_plus_enhancer' as feature_type,if(feature_strand = 1,feature_start-2500,feature_start) as feature_start, if(feature_strand = 1,feature_end,feature_end +2500) as feature_end, feature_strand from PolIII_transcribed_gene";

    push @sql, _col_types_and_indices('PolIII_transcribed_gene_plus_enhancer');

    _run_sql($dbc,@sql) or throw "cannot run sql";
}


sub _create_gene_features{
    my ($dbc) = @_;

    my @sql;
    push @sql,"drop table if exists protein_coding_gene";
    push @sql,"create table protein_coding_gene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,'protein_coding_gene' as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('protein_coding')";
    push @sql,"alter table protein_coding_gene add index(seq_region_name)";
    push @sql,"alter table protein_coding_gene add index(seq_region_id)";


    push @sql,"drop table if exists pseudogene";
    push @sql,"create table pseudogene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('pseudogene','repeat','retrotransposed')";
    push @sql,"alter table pseudogene add index(seq_region_name)";
    push @sql,"alter table pseudogene add index(seq_region_id)";


    push @sql,"drop table if exists RNA_gene";
    push @sql,"create table RNA_gene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype not like '%pseudogene%'";
    push @sql,"alter table RNA_gene add index(seq_region_name)";
    push @sql,"alter table RNA_gene add index(seq_region_id)";


    #This section is possibly not necessary...
    push @sql,"drop table if exists processed_transcript"; # used for QC
    push @sql,"create table processed_transcript select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('processed_transcript')";
    push @sql,"alter table processed_transcript add index(seq_region_name)";
    push @sql,"alter table processed_transcript add index(seq_region_id)";


    push @sql,"drop table if exists RNA_pseudogene";
    push @sql,"create table RNA_pseudogene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype like '%pseudogene%'";
    push @sql,"alter table RNA_pseudogene add index(seq_region_name)";
    push @sql,"alter table RNA_pseudogene add index(seq_region_id)";


    push @sql,"drop table if exists lincRNA_gene";
    push @sql,"create table lincRNA_gene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype = 'lincRNA' ";
    push @sql,"alter table lincRNA_gene add index(seq_region_name)";
    push @sql,"alter table lincRNA_gene add index(seq_region_id)";

    _run_sql($dbc,@sql) or throw "cannot run sql";

    @sql=();
    push @sql,"drop table if exists PolIII_transcribed_gene";
    push @sql,"create table PolIII_transcribed_gene select * from lincRNA_gene where 1=0"; # make empty table

    foreach my $desc ('Nuclear RNase P','RNase MRP','Vault','Y RNA','5S ribosomal RNA','7SK RNA',){
	push @sql,"insert into PolIII_transcribed_gene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id, '$desc' as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.description like '$desc".'%'."'";
    }
    push @sql,"alter table PolIII_transcribed_gene add index(seq_region_name)";
    push @sql,"alter table PolIII_transcribed_gene add index(seq_region_id)";

    _run_sql($dbc,@sql) or throw "cannot run sql";

}

# creates tables storing transcript related features
sub _create_transcript_features{
    my ($dbc) = @_;

    my @sql;
    push @sql,"drop table if exists protein_coding_transcript";
    push @sql,"create table protein_coding_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,'protein_coding_transcript' as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('protein_coding')";
    push @sql,"alter table protein_coding_transcript add index(seq_region_name)";
    push @sql,"alter table protein_coding_transcript add index(seq_region_id)";

    push @sql,"drop table if exists protein_coding_transcript_downstream_2500";
    push @sql,"create table protein_coding_transcript_downstream_2500 select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,'protein_coding_transcript_downstream_2500' as feature_type,if(seq_region_strand = 1, t.seq_region_end+1,t.seq_region_start -2500) as feature_start,if(seq_region_strand = 1,t.seq_region_end+2500,t.seq_region_start - 1) as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('protein_coding')";
    push @sql,"alter table protein_coding_transcript_downstream_2500 add index(seq_region_name)";
    push @sql,"alter table protein_coding_transcript_downstream_2500 add index(seq_region_id)";

    push @sql,"drop table if exists pseudogene_transcript";
    push @sql,"create table pseudogene_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,concat(biotype,'_transcript') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('pseudogene','repeat','retrotransposed')";
    push @sql,"alter table pseudogene_transcript add index(seq_region_name)";
    push @sql,"alter table pseudogene_transcript add index(seq_region_id)";

    push @sql,"drop table if exists RNA_transcript";
    push @sql,"create table RNA_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,concat(biotype,'_transcript') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype not like '%pseudogene%'";
    push @sql,"alter table RNA_transcript add index(seq_region_name)";
    push @sql,"alter table RNA_transcript add index(seq_region_id)";

    push @sql,"drop table if exists RNA_pseudogene_transcript";
    push @sql,"create table RNA_pseudogene_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,concat(biotype,'_transcript') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype like '%pseudogene%'";
    push @sql,"alter table RNA_pseudogene_transcript add index(seq_region_name)";
    push @sql,"alter table RNA_pseudogene_transcript add index(seq_region_id)";

    _run_sql($dbc,@sql) or throw "Couldn't create transcript features";

}

#for repeats on top level sequence
sub _create_repeat_features{
    my($dbc)=@_;

    my @sql;
    my $temp1 = "tmp_".int(rand(1000));
    # we lose a few mappings because some mappings are on parts of the contigs
    # which are not included in the assembly
    push @sql,"drop table if exists $temp1";
    
        my $q= "create table $temp1
        select f.*,
               f.seq_region_start as chr_start,
               f.seq_region_end as chr_end,
               f.seq_region_id as chromosome_id,
               asm_sr.name as chr_name 
from repeat_feature f,
     analysis an,
     seq_region asm_sr
where f.analysis_id = an.analysis_id
and an.logic_name = 'RepeatMask'
and f.seq_region_id = asm_sr.seq_region_id ";

	$q =~ tr/\n/ /;
	push @sql, $q;

	push @sql,"alter table $temp1 add index(repeat_consensus_id)";
        # now add the repeat type 
        push @sql,"drop table if exists all_repeats";
        push @sql,"create table all_repeats select t1.chromosome_id as seq_region_id,t1.chr_name as seq_region_name,repeat_feature_id as feature_id,concat(rc.repeat_class,'_repeat') as feature_type,chr_start as feature_start,chr_end as feature_end , 0 as feature_strand from $temp1 t1, repeat_consensus rc where rc.repeat_consensus_id = t1.repeat_consensus_id ";
        # we want the U6 snRNA, 5S rRNA, 7SLRNA etc separately for pol3 transcribed regions
        push @sql,"insert into all_repeats select t1.chromosome_id as seq_region_id,t1.chr_name as seq_region_name,repeat_feature_id as feature_id,'PolIII_transcribed_repeat' as feature_type,chr_start as feature_start,chr_end as feature_end , 0 as feature_strand from $temp1 t1, repeat_consensus rc where rc.repeat_consensus_id = t1.repeat_consensus_id and rc.repeat_name in ('U6','5S','7SLRNA','7SK') ";

	push @sql, "alter table all_repeats add index(feature_type)";

        push @sql,"drop table $temp1";

    
    _run_sql($dbc,@sql) or throw "cannot run sql";

    my $aref = ['Satellite_repeat','Satellite/centr_repeat','LTR/ERV1_repeat','LTR/MaLR_repeat','SINE/Alu_repeat' ];
    @sql = ();
    foreach my $type (@$aref){
        @sql = ();
	my $table = $type;
	$table =~ tr/\//_/;
	$table =~ tr/-/_/;
	$table =~ tr/?/_/;
	push @sql,"drop table if exists $table";
	push @sql,"create table $table select * from all_repeats where feature_type = '$type'";
	push @sql,"update $table set feature_type = '$table'";
        push @sql, _col_types_and_indices("$table");
        _run_sql($dbc,@sql) or throw "cannot run sql for $type";
	#print "$table\t$table\t$table\t0\t0\t0\n";
    }

}

sub _create_exon_features{
    my($dbc)=@_;
    my @sql;

    my $temp1 = "tmp_".int(rand(1000))."_1";
    # denormalise exon, transcript and exon_transcript
    push @sql, "drop table if exists $temp1";
    push @sql, "create table $temp1 select e.*,et.rank,et.transcript_id,t.biotype from exon e, exon_transcript et, transcript t where e.exon_id = et.exon_id and et.transcript_id = t.transcript_id order by et.transcript_id,e.seq_region_start";
    push @sql, "alter table $temp1 add index(transcript_id)";

    # get rank of last exon for each transcript
    my $temp2 = "tmp_".int(rand(1000))."_2";
    push @sql, "drop table if exists $temp2";
    push @sql, "create table $temp2 select transcript_id,max(rank) as max_rank from $temp1 group by transcript_id";
    push @sql, "alter table $temp2 add index(transcript_id)";
    my $temp3 = "tmp_".int(rand(1000))."_3";
    push @sql, "drop table if exists $temp3";
    push @sql, "create table $temp3 select t1.*,t2.max_rank from $temp1 t1, $temp2 t2 where t1.transcript_id = t2.transcript_id";
    push @sql, "drop table if exists $temp2";
    push @sql, "drop table if exists $temp1";
    push @sql, "alter table $temp3 add index(transcript_id)";
    push @sql, "alter table $temp3 add index(max_rank)";
    push @sql, "alter table $temp3 add index(rank)";
    push @sql, "alter table $temp3 add index(seq_region_id)";
    

    # get first intron for each transcript
    my $temp4 =  "tmp_".int(rand(1000))."_4";
    push @sql, "create table $temp4 select t3.seq_region_id,sr.name as seq_region_name,t3.transcript_id as feature_id, concat(t3.biotype,'_intron1') as feature_type,if(t3.seq_region_strand = 1,t3.seq_region_end+1,b.seq_region_end+1) as feature_start, if(t3.seq_region_strand = 1,b.seq_region_start-1,t3.seq_region_start -1) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, $temp3 b, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.transcript_id = b.transcript_id and t3.rank = 1 and b.rank = 2"; 
    push @sql,"drop table if exists intron1";
    push @sql,"alter table $temp4 rename as intron1";
    push @sql,_col_types_and_indices("intron1");
    push @sql,_split_by_biotype("intron1");
   _run_sql($dbc,@sql);

    # promoter defined as 500 bp upstream
    push @sql, "drop table if exists exon1_plus_promoter";
    push @sql, "create table exon1_plus_promoter select t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, concat(t3.biotype,'_exon1_plus_promoter') as feature_type,if(seq_region_strand = 1,t3.seq_region_start-500,seq_region_start) as feature_start, if(seq_region_strand = 1,seq_region_end,seq_region_end +500) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.rank = 1";
    push @sql,_col_types_and_indices("exon1_plus_promoter");
    push @sql,_split_by_biotype("exon1_plus_promoter");

    
    # enhancer defined as 2500 bp upstream
    push @sql, "drop table if exists exon1_plus_enhancer";
    push @sql, "create table exon1_plus_enhancer select t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, concat(t3.biotype,'_exon1_plus_enhancer') as feature_type,if(seq_region_strand = 1,t3.seq_region_start-2500,seq_region_start) as feature_start, if(seq_region_strand = 1,seq_region_end,seq_region_end +2500) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.rank = 1";
    push @sql,_col_types_and_indices("exon1_plus_enhancer");
    push @sql,_split_by_biotype("exon1_plus_enhancer");

    push @sql, "drop table if exists single_exon_gene";
    push @sql, "create table single_exon_gene select t3.seq_region_id,sr.name as seq_region_name,exon_id as feature_id, concat(t3.biotype, '_single_exon_gene') as feature_type,t3.seq_region_start as feature_start, t3.seq_region_end as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.max_rank = 1";
    push @sql,_col_types_and_indices("single_exon_gene");    
    push @sql,_split_by_biotype("single_exon_gene");

    push @sql, "drop table if exists single_exon_gene_plus_enhancer";
    push @sql, "create table single_exon_gene_plus_enhancer select  t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, concat(t3.biotype,'_single_exon_gene_plus_enhancer') as feature_type,if(seq_region_strand = 1,t3.seq_region_start-2500,seq_region_start) as feature_start, if(seq_region_strand = 1,seq_region_end,seq_region_end +2500) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.max_rank = 1";
    push @sql,_col_types_and_indices("single_exon_gene_plus_enhancer");    
    push @sql,_split_by_biotype("single_exon_gene_plus_enhancer");

    # gene body is all exons and introns except the first of each
    push @sql, "drop table if exists gene_body";
    push @sql, "create table gene_body select t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,  concat(t3.biotype,'_gene_body') as feature_type,min(t3.seq_region_start) as feature_start, max(t3.seq_region_end) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.max_rank != 1 and t3.rank != 1 group by transcript_id";
    push @sql,_col_types_and_indices("gene_body");  
    push @sql,_split_by_biotype("gene_body");


    push @sql, "drop table if exists exon_plus_flanks_500";
    push @sql, "create table exon_plus_flanks_500 select distinct t3.seq_region_id,sr.name as seq_region_name,exon_id as feature_id, concat(t3.biotype,'_exon_plus_flanks_500') as feature_type,t3.seq_region_start-500 as feature_start, seq_region_end + 500 as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id ";
    push @sql,_col_types_and_indices("exon_plus_flanks_500");    
    push @sql,_split_by_biotype("exon_plus_flanks_500");


    push @sql, "drop table if exists $temp3";

    _run_sql($dbc,@sql) or throw "cannot run sql";
    
}

sub _split_by_biotype{
    my($feat_name) = @_;

    my @sql;

    my $prots = "'protein_coding_$feat_name', 'IG_V_gene_$feat_name', 'IG_C_gene_$feat_name', 'IG_J_gene_$feat_name', 'IG_D_gene_$feat_name'";
    push @sql,"drop table if exists protein_coding_$feat_name";
    push @sql,"create table protein_coding_$feat_name select * from $feat_name where feature_type in ($prots)";
    push @sql,"update protein_coding_$feat_name set feature_type = 'protein_coding_$feat_name'"; 
    push @sql,_col_types_and_indices("protein_coding_$feat_name");

    push @sql,"drop table if exists RNA_gene_$feat_name";
    push @sql,"create table RNA_gene_$feat_name select * from $feat_name where feature_type like '%RNA_".$feat_name."'";
    push @sql,_col_types_and_indices("RNA_gene_$feat_name");

    if($feat_name eq 'exon1_plus_enhancer'){
	push @sql,"drop table if exists snRNA_gene_$feat_name";
	push @sql,"create table snRNA_gene_$feat_name select * from $feat_name where feature_type like '%snRNA_".$feat_name."'";
	push @sql,"update snRNA_gene_$feat_name set feature_type = 'snRNA_gene_$feat_name'";
	push @sql,_col_types_and_indices("snRNA_gene_$feat_name");

	push @sql,"drop table if exists snoRNA_gene_$feat_name";
	push @sql,"create table snoRNA_gene_$feat_name select * from $feat_name where feature_type like '%snoRNA_".$feat_name."'";
	push @sql,"update snoRNA_gene_$feat_name set feature_type = 'snoRNA_gene_$feat_name'";
	push @sql,_col_types_and_indices("snoRNA_gene_$feat_name");

	push @sql,"drop table if exists miRNA_gene_$feat_name";
	push @sql,"create table miRNA_gene_$feat_name select * from $feat_name where feature_type like '%miRNA_".$feat_name."'";
	push @sql,"update miRNA_gene_$feat_name set feature_type = 'miRNA_gene_$feat_name'";
	push @sql,_col_types_and_indices("miRNA_gene_$feat_name");

	push @sql,"drop table if exists miscRNA_gene_$feat_name";
	push @sql,"create table miscRNA_gene_$feat_name select * from $feat_name where feature_type like '%miscRNA_".$feat_name."'";
	push @sql,"update miscRNA_gene_$feat_name set feature_type = 'miscRNA_gene_$feat_name'";
	push @sql,_col_types_and_indices("miscRNA_gene_$feat_name");

	push @sql,"drop table if exists rRNA_gene_$feat_name";
	push @sql,"create table rRNA_gene_$feat_name select * from $feat_name where feature_type like '%rRNA_".$feat_name."'";
	push @sql,"update rRNA_gene_$feat_name set feature_type = 'rRNA_gene_$feat_name'";
	push @sql,_col_types_and_indices("rRNA_gene_$feat_name");
    }

    my $pseuds = "'pseudogene_$feat_name','repeat_$feat_name','retrotransposed_$feat_name'";
    push @sql,"drop table if exists pseudogene_$feat_name";
    push @sql,"create table pseudogene_$feat_name select * from $feat_name where feature_type in ($pseuds)";
    push @sql,_col_types_and_indices("pseudogene_$feat_name");

    return @sql;
}


sub _cage_ditag_transcript_tss{
    my($dbc)=@_;

    my $dbh = $dbc->db_handle;

    my @sql;
    my $temp1 = "tmp_".int(rand(1000))."_1";
    push @sql, "drop table if exists $temp1";
    push @sql, "create table $temp1 select f.seq_region_id,sr.name as seq_region_name,if(f.seq_region_strand = 1,f.seq_region_start,f.seq_region_end) as feature_start,f.seq_region_strand as feature_strand from transcript f , seq_region sr where sr.seq_region_id = f.seq_region_id and f.biotype not like '%pseudogene%' and f.biotype not in('repeat','retrotransposed')";

    if($dbh->selectrow_array("select count(*) from information_schema.tables where table_schema ='".$dbc->dbname."' and table_name='ditag_feature'") > 0){
    #if($dbu->table_exists("ditag_feature")){
        push @sql, "insert into $temp1  select f.seq_region_id,sr.name as seq_region_name,if(f.seq_region_strand = 1,f.seq_region_start,f.seq_region_end) as feature_start,f.seq_region_strand from ditag_feature f , seq_region sr where sr.seq_region_id = f.seq_region_id";
    }

    my $temp2 = "tmp_".int(rand(1000))."_2";
    push @sql, "drop table if exists $temp2";
    push @sql, "create table $temp2 select distinct * from $temp1";

    # we may want to do clustering here to reduce the total number of features

    push @sql, "drop table if exists tss_upstream_500";
    push @sql, "create table tss_upstream_500 select seq_region_id,seq_region_name,if(feature_strand = -1,feature_start,feature_start-500) as feature_start,if(feature_strand = -1,feature_start+500,feature_start) as feature_end,feature_strand,'tss_upstream_500' as feature_type from $temp2";
push @sql,"alter table tss_upstream_500 add column feature_id int(10) not null auto_increment primary key";
    push @sql,_col_types_and_indices("tss_upstream_500");


    push @sql, "drop table if exists tss_downstream_500";
    push @sql, "create table tss_downstream_500 select seq_region_id,seq_region_name,if(feature_strand = -1,feature_start-500,feature_start) as feature_start,if(feature_strand = -1,feature_start,feature_start+500) as feature_end,feature_strand,'tss_downstream_500' as feature_type from $temp2";
push @sql,"alter table tss_downstream_500 add column feature_id int(10) not null auto_increment primary key";
    push @sql,_col_types_and_indices("tss_downstream_500");


    push @sql, "drop table if exists tss_centred_500";
    push @sql, "create table tss_centred_500 select seq_region_id,seq_region_name,feature_start-250 as feature_start,feature_start+250 as feature_end,feature_strand,'tss_centred_500' as feature_type from $temp2";
push @sql,"alter table tss_centred_500 add column feature_id int(10) not null auto_increment primary key";
    push @sql,_col_types_and_indices("tss_centred_500");



    push @sql, "drop table if exists tss_centred_5000";
    push @sql, "create table tss_centred_5000 select seq_region_id,seq_region_name,feature_start-2500 as feature_start,feature_start+2500 as feature_end,feature_strand,'tss_centred_5000' as feature_type from $temp2";
push @sql,"alter table tss_centred_5000 add column feature_id int(10) not null auto_increment primary key";
    push @sql,_col_types_and_indices("tss_centred_5000");

    # 
    _run_sql($dbh,@sql) or die;

    @sql = ();
    push @sql,"drop table if exists $temp1";
    push @sql,"drop table if exists $temp2";
    _run_sql($dbh,@sql) or die;

}

sub _create_intergenic_features{
    my ($dbc) = @_;

    my $dbh = $dbc->db_handle;

    # we process one seq_region ata a time
    # if we want to include more biotypes of gene then we should select from 
    # the gene table
    my $q = "select distinct seq_region_id from protein_coding_gene";
    my $regions_aref = $dbh->selectcol_arrayref($q);
    unless(@$regions_aref > 0){ throw "no data returned by :\n$q\n" }

    my @sql;
    my $temp1 = "tmp_".int(rand(1000));
    push @sql, "drop table if exists $temp1";
    # start and end are signed so we can use -ve numbers as flags
    push @sql, "create table $temp1 (seq_region_id int(10) unsigned,feature_start int(10), feature_end int(10))";
    _run_sql($dbh,@sql) or throw "error running sql";

    # first we get the regions which lie around/between the genes
    my @intergenics;
    foreach my $id (@$regions_aref){

	@intergenics = ();

	my $region_length = $dbh->selectrow_array("select length from seq_region where seq_region_id = $id");

	$q = "select feature_start,feature_end from protein_coding_gene where seq_region_id = $id order by feature_start";
	my $genes_aaref = $dbh->selectall_arrayref($q);
        unless(defined $genes_aaref){
            throw "query failed on:\n$q\n".$dbh->errstr;
	}

	my $last_end;
        if($genes_aaref->[0]->[0] == 1){
	    $last_end = $genes_aaref->[0]->[1];
	}else{
	    my @intergenic;
            $intergenic[0] = -1; # use -ve to flag chromosome end
            $intergenic[1] = $genes_aaref->[0]->[0] -1;
	    push @intergenics,\@intergenic;
            $last_end = $genes_aaref->[0]->[1];
	}

        
        for(my $i=1;$i< @$genes_aaref;$i++){

	    #print "gene $i ".join("\t",@{$genes_aaref->[$i]})."\n";
	    my @intergenic;
	    if($genes_aaref->[$i]->[0] > $last_end+1){
                $intergenic[0] = $last_end+1;
                $intergenic[1] = $genes_aaref->[$i]->[0] -1;
 	        push @intergenics,\@intergenic;
	    }

            if($genes_aaref->[$i]->[1] > $last_end){
                $last_end = $genes_aaref->[$i]->[1];
	    }

	}

        if($last_end < $region_length){
	    my @intergenic;
            $intergenic[0] = $last_end+1;
            $intergenic[1] = - $region_length; # use -ve to flag chromosome end
 	    push @intergenics,\@intergenic;
	}

	@sql = ();
        foreach my $aref (@intergenics){
	    my $q = "insert into $temp1 values($id,";
	    $q .= join(',',@$aref);
	    $q .= ")";
	    push @sql,$q;

	}
	_run_sql($dbh,@sql) or throw "Error running sql";
	
    }

    # now we create gene-distal regions by shortening the regions in $temp1
    _intergenic_variant($dbh,$temp1,2500);
    _intergenic_variant($dbh,$temp1,5000);
    _intergenic_variant($dbh,$temp1,10000);

    _run_sql($dbh,"drop table $temp1") or throw "cannot drop $temp1 table";

}


sub _intergenic_variant{
    my($dbh,$temp1,$dist) = @_;

    my $final_table = 'intergenic_'.$dist;
    my @sql;
    push @sql,"drop table if exists $final_table";
    my $temp2 = "tmp_".int(rand(1000));
    push @sql,"drop table if exists $temp2";
    push @sql,"create table $temp2 select seq_region_id,feature_start,feature_end from $temp1 where 1=0";
    # only add/subtract $dist if we are not at the chromosome end
    push @sql,"insert into $temp2 select seq_region_id,if(feature_start> 0,feature_start +$dist,feature_start * -1) as feature_start,if(feature_end > 0,feature_end - $dist, feature_end * -1) as feature_end from $temp1";
    
    push @sql,"delete from $temp2 where feature_end <= feature_start";
    push @sql,"create table $final_table select t2.seq_region_id,sr.name as seq_region_name,'$final_table' as feature_type,t2.feature_start,t2.feature_end,'0' as feature_strand from $temp2 t2, seq_region sr where t2.seq_region_id = sr.seq_region_id";
    push @sql,"alter table $final_table add column feature_id int(10) not null auto_increment primary key";

    push @sql,_col_types_and_indices($final_table);

    _run_sql($dbh,@sql) or throw "cannot run sql";
    _run_sql($dbh,"drop table $temp2") or throw "cannot drop $temp2 table";

}


sub _col_types_and_indices{
    my $table = shift;

    my @sql;

    push @sql,"alter table $table modify column feature_start int(10) unsigned";
    push @sql,"alter table $table modify column feature_end int(10) unsigned";
    push @sql,"alter table $table modify column feature_strand tinyint(2)";
    push @sql,"alter table $table modify column feature_type varchar(45)";
    push @sql,"alter table $table add index(seq_region_id)";
    push @sql,"alter table $table add index(seq_region_name)";
    push @sql,"alter table $table add index(feature_strand)";
    push @sql,"alter table $table add index(feature_start)";
    push @sql,"alter table $table add index(feature_end)";
    return @sql;
}

sub _process_regulatory_features {
  my ($dnadb, $efgdb, $workdb, $cell_name) = @_;
      
  _copy_with_rename($efgdb->db_handle,'seq_region',$workdb->db_handle,'func_seq_region'); 

  my @sql = ();

  push @sql, "drop table if exists seq_name_lookup";
  push @sql, "create table seq_name_lookup select distinct seq_region_id,name as seq_region_name from func_seq_region";
  push @sql, "alter table seq_name_lookup add index(seq_region_id)";
  push @sql, "alter table seq_name_lookup add index(seq_region_name)";
  _run_sql($workdb,@sql) or throw "cannot create seq_name_lookup";

  #add seq_region_name to regulatory_features
  @sql = ();
  my $temp1 = "tmp_".int(rand(1000));
  push @sql, "create table $temp1 select r.*,s.seq_region_name from regulatory_feature r,seq_name_lookup s where r.seq_region_id = s.seq_region_id";
  push @sql, "drop table regulatory_feature";
  push @sql, "alter table $temp1 rename as regulatory_feature";
  push @sql, "alter table regulatory_feature add index(feature_set_id)";
  _run_sql($workdb,@sql) or throw "cannot update regulatory_feature table";

  @sql = ();
  # remove all but the current set of features
  my $feature_set_name ='RegulatoryFeatures';
  if($cell_name){
    $feature_set_name .= ':'.$cell_name;
  } 
  
  my $kh = "select feature_set_id from feature_set where name = '$feature_set_name'";
  my $feature_set_id = $workdb->db_handle->selectrow_array($kh);
  if($feature_set_id == 0){ throw "Could not find feature_set_id for $feature_set_name"; }

  push @sql, "delete  from regulatory_feature where feature_set_id != $feature_set_id";
  push @sql, "alter table regulatory_feature add index(regulatory_feature_id)";
  push @sql, "alter table regulatory_feature add index(seq_region_start)";
  push @sql, "alter table regulatory_feature add index(seq_region_end)";
  push @sql, "alter table regulatory_feature add index(seq_region_name)";
  push @sql, "alter table regulatory_feature add index(feature_type_id)";
  _run_sql($workdb,@sql) or throw "could not delete irrelevant regulatory features";
  
  # we remove features which have very long 'whiskers'
  # these tend to contain multiple focus features which results in
  # a group of regulatory features all of which have the same attributes
  @sql = ();
  my $filtered_features_table = 'regulatory_features_filtered';
  push @sql,"drop table if exists $filtered_features_table";
  push @sql,"create table $filtered_features_table select *, 'regulatory_feature' as feature_type from regulatory_feature";
  push @sql, "alter table $filtered_features_table add index(seq_region_name)";
  push @sql, "alter table $filtered_features_table add index(feature_type_id)";
  push @sql, "alter table $filtered_features_table add index(feature_set_id)";
  push @sql, "alter table $filtered_features_table add index(regulatory_feature_id)";
  push @sql, "alter table $filtered_features_table add index(seq_region_start)";
  push @sql, "alter table $filtered_features_table add index(seq_region_end)";
  
  my $temp2 = "tmp_".int(rand(1000));
  my $temp3 = "tmp_".int(rand(1000));
  
  push @sql,"create table $temp3 select ra.regulatory_feature_id , af.* from regulatory_attribute ra, annotated_feature af where ra.attribute_feature_id = af.annotated_feature_id";

  push @sql,"create table $temp2 select j.regulatory_feature_id,count(*) as n_attribs,count(distinct j.feature_set_id) as n_attrib_types,max(j.seq_region_end)-min(j.seq_region_start) as len,min(j.seq_region_start)as attribs_start,max(j.seq_region_end) as attribs_end,rf.binary_string from $temp3 j, $filtered_features_table rf where j.regulatory_feature_id=rf.regulatory_feature_id group by regulatory_feature_id ";
  push @sql,"alter table $temp2 add index(regulatory_feature_id)";
  push @sql,"alter table $temp2 add index(attribs_start)";
  push @sql,"alter table $temp2 add index(attribs_end)";

  push @sql,"delete r from $filtered_features_table r, $temp2 t2 where r.regulatory_feature_id = t2.regulatory_feature_id and t2.len > 5000";

  push @sql, "alter table $filtered_features_table add index(seq_region_name)";
  push @sql, "alter table $filtered_features_table add index(feature_type_id)";
  push @sql, "alter table $filtered_features_table add index(feature_set_id)";
  push @sql, "alter table $filtered_features_table add index(regulatory_feature_id)";
  push @sql, "alter table $filtered_features_table add index(seq_region_start)";
  push @sql, "alter table $filtered_features_table add index(seq_region_end)";
  push @sql, "alter table $filtered_features_table add index(binary_string)";
  # we want to remove features which occur in the same place and have the same
  # attributes as one another, leaving only one representative
  push @sql, "delete b from $filtered_features_table a, $filtered_features_table b, $temp2 t2a, $temp2 t2b where a.seq_region_name=b.seq_region_name and b.binary_string = a.binary_string and a.regulatory_feature_id < b.regulatory_feature_id and a.regulatory_feature_id = t2a.regulatory_feature_id and b.regulatory_feature_id = t2b.regulatory_feature_id and t2a.attribs_start=t2b.attribs_start and t2a.attribs_end = t2b.attribs_end";

  push @sql,"drop table if exists $temp2";

  _run_sql($workdb,@sql) or throw "cannot run sql";

  # reg  features which overlap centromeric repeats tend to contain too many marks# which suggests there is something wrong with the mappings.
  # so we get rid of them

  @sql = ();
  push @sql,"drop table if exists $temp3";
  
  push @sql,"create table $temp3 select f.regulatory_feature_id from $filtered_features_table f, Satellite_centr_repeat s where s.seq_region_name = f.seq_region_name and s.feature_end >= f.seq_region_start and s.feature_start <= f.seq_region_end";
  push @sql,"delete f from $filtered_features_table f,$temp3 t where f.regulatory_feature_id = t.regulatory_feature_id";


  push @sql,"drop table if exists $temp3";
  
  push @sql,"create table $temp3 select f.regulatory_feature_id from $filtered_features_table f, Satellite_repeat s where s.seq_region_name = f.seq_region_name and s.feature_end >= f.seq_region_start and s.feature_start <= f.seq_region_end";
  push @sql,"delete f from $filtered_features_table f,$temp3 t where f.regulatory_feature_id = t.regulatory_feature_id";

  push @sql,"drop table if exists $temp3";
  
  # reg feats on the mitochondrial DNA are unlikely to use the same
  # 'histone code' as the rest of the genome. they might even be artefacts.
  # so we remove them (this should not be necessary as there should not be features on MT..
  push @sql, "delete from $filtered_features_table where seq_region_name = 'MT'";
  _run_sql($workdb,@sql) or throw "cannot run sql";

  @sql=();
  # we make the col names compatible with the genomic features in the analysis
  push @sql,"alter table $filtered_features_table change column seq_region_start feature_start int(10) unsigned";
  push @sql,"alter table $filtered_features_table change column seq_region_end feature_end int(10) unsigned";
  push @sql,"alter table $filtered_features_table change column seq_region_strand feature_strand tinyint(2)";
  push @sql,"alter table $filtered_features_table modify column feature_type varchar(45)";

  _run_sql($workdb,@sql) or throw "cannot run sql";

}

sub _generate_random_frags {
  my ($workdb) = @_;
  
  my $originals_table = "regulatory_features_filtered";
  my $outtable = "mockreg_features_filtered";

  my @sql;

  push @sql,"drop table if exists banned_region_table";
  push @sql,"create table banned_region_table ( seq_region_name varchar(10), feature_start int(11), feature_end int(11))";
  my $dbh = $workdb->db_handle;
  _run_sql($dbh,@sql) or throw "Cannot run sql";

  @sql = ();

  # first we need the lengths of the seq regions so we know the max 
  # tfrag_enc_end value we can generate
  my %enc_len;
  my $q = "select distinct seq_region_name from $originals_table";
  my $col_ref = $dbh->selectcol_arrayref($q) or throw $dbh->errstr;
  my $clause = " where name in('";
  $clause .= join("','",@$col_ref);
  $clause .= "')";

  $q="select name, length from seq_region" .$clause;
  my $aaref = $dbh->selectall_arrayref($q) or throw $dbh->errstr;
  foreach my $aref (@$aaref){
    $enc_len{$aref->[0]} = $aref->[1];
  }

  my $tmp_table = "tmp_".int(rand(1000));
  push @sql, "drop table if exists $tmp_table";
  push @sql, "create table $tmp_table (regulatory_feature_id int(10) unsigned,seq_region_name char(10),feature_start int(10) unsigned,feature_end int(10) unsigned,feature_strand tinyint(2),feature_type varchar(40))"; 

  # next we need the seq_region and length of each original feature
  my %frag_len;
  $q = "select seq_region_name, feature_end - feature_start +1 as length,feature_type,feature_strand,regulatory_feature_id from $originals_table order by (feature_end - feature_start +1) desc ";
  $aaref = $dbh->selectall_arrayref($q) or throw "failed on \n$q\n".$dbh->errstr;
  
  foreach my $aref (@$aaref){
    my $enc = $aref->[0];    
    my $len = $aref->[1];
    my $type = $aref->[2];
    my $ori = $aref->[3];
    my $id =  $aref->[4]; 
    
    my($start,$end) = _get_rand_region($dbh,$enc,$len,%enc_len);

    my $q = "insert into $tmp_table values(";
    $q .= $id.",";
    $q .= "\'".$enc."\'," ;
    $q .= $start.",";
    $q .= $end.",";
    $q .= "\'".$ori."\'," ;
    $q .= "\'".$type."\'" ;
    $q .= ")";

    push @sql,$q;

    #Add frag to the banned list (so it won't be generated again)...
    my $frag_q = "insert into banned_region_table values('$enc',$start,$end) ";
    _run_sql($dbh,$frag_q) or throw("failed on \n$q\n".$dbh->errstr);
    
  }
  
  push @sql, "drop table if exists $outtable";
  push @sql, "create table $outtable select t.*,o.binary_string from $tmp_table t, $originals_table o where t.regulatory_feature_id = o.regulatory_feature_id"; # select f.seq_region_name, f.feature_type,f.feature_start,f.feature_end ,f.feature_strand from temp f";
  
  push @sql, "alter table $outtable add index(regulatory_feature_id)";
  push @sql, "alter table $outtable add index(seq_region_name)";
  push @sql, "alter table $outtable add index(feature_start)";
  push @sql, "alter table $outtable add index(feature_type)";
  push @sql, "alter table $outtable add index(feature_strand)";

  push @sql, "drop table if exists banned_region_table";
  push @sql, "drop table if exists $tmp_table";

  _run_sql($dbh,@sql) or throw "Cannot run sql";


}

sub _calculate_associations{
  my ($dbc, $dumpdir, $gen_feats, $species) = @_;

  #Pass these as parameters?
  my $reg_feat_table = 'regulatory_features_filtered';
  my $mock_reg_table = 'mockreg_features_filtered';
  my $gen_feat_table = 'genomic_features';
  my $overlaps_table = 'reg_feat_gen_feat_overlaps';
  my $mock_olaps_table = 'mock_feat_gen_feat_overlaps';
  my $flags_table = 'regulatory_feature_association_flags';
  my $types_table = 'regulatory_features_classified';

  my $dbh = $dbc->db_handle;

  system("mkdir -p ".$dumpdir) && throw "Error creating $dumpdir folder";

  #@gen_feats is the list of feature classes relevant for classification...
  _create_gen_feats_table($dbh,$gen_feats,$gen_feat_table);
  

  # open a file for writing the results
  my $reg_results_file = $dumpdir.'/reg_results.tab';
  my $mock_results_file = $dumpdir.'/mock_results.tab';
  
  my $reg_ofh;
  my $mock_ofh;
  open($reg_ofh,"> $reg_results_file") or die "couldn't open $reg_results_file";
  open($mock_ofh,"> $mock_results_file") or die "couldn't open $mock_results_file";
  # filename for temp storage of genomic features
  my $gen_file = $dumpdir.'/genomic_features';

  my $q = "select distinct seq_region_name from $gen_feat_table";
  my $chrs_ref = $dbh->selectcol_arrayref($q);
  foreach my $chr (@$chrs_ref){

    my $command = "mysql -h ".$dbc->host." -u ".$dbc->username." -P ".$dbc->port." -p".$dbc->password." ".$dbc->dbname.' -BN -e"'.
		  "select feature_type_id,feature_start,feature_end from  genomic_features where seq_region_name = '$chr' order by feature_start,feature_end".'"'.
		    " > $gen_file";
    system($command) && throw "Error dumping genomics features for $chr";

    _overlap_analysis($dbh,$gen_file,$chr,$reg_feat_table,$reg_ofh,scalar(@{$gen_feats}));
    _overlap_analysis($dbh,$gen_file,$chr,$mock_reg_table,$mock_ofh,scalar(@{$gen_feats}));
  }
  close($reg_ofh);
  close($mock_ofh);

  _load_results($dbh,$reg_results_file,$overlaps_table,$reg_feat_table,$gen_feats);
  _load_results($dbh,$mock_results_file,$mock_olaps_table,$mock_reg_table,$gen_feats);

  my $pat_bits = $dbh->selectrow_array("select length(binary_string) from regulatory_features_filtered limit 1");
  my @pats = _sub_pats($dbh,$reg_feat_table,$pat_bits);
  
  my $patt_count_thresh = 100;

  my $overlap_results_file = $dumpdir.'/overlap_results.tab';
  my $overlap_ofh;
  open($overlap_ofh,"> $overlap_results_file") or throw "couldn't open $overlap_results_file";
  #for each pattern we want the number of reg feats with that pattern  
  # and then for each genomic feature we want the overlap counts for the
  # real and mock reg feats with that pattern
  my($reg_count,$real_count,$mock_count);
  foreach my $pat (@pats){
    $reg_count = $dbh->selectrow_array("select count(1) from  $overlaps_table where pattern like '$pat'");
    if($reg_count >= $patt_count_thresh){
      print $overlap_ofh "$pat\t$reg_count";
      foreach my $gen (@{$gen_feats}){
	$real_count = $dbh->selectrow_array("select count(1) from  $overlaps_table where pattern like '$pat' and $gen");
	$mock_count = $dbh->selectrow_array("select count(1) from  $mock_olaps_table where pattern like '$pat' and $gen");
	print $overlap_ofh "\t$real_count\t$mock_count";
      }
      print $overlap_ofh "\n" or throw "failed to print to file";
    }
  }
  close($overlap_ofh); 

  _create_pattern_results_tables($dbh, $overlap_results_file,$gen_feats);

  my $assoc_thresh = 51; #70;
  my $second_thresh = 50;
  _association_tables($dbh,$gen_feats,'pattern_overlap_summary','pattern_overlap_summary_chi',$assoc_thresh,$second_thresh);

  _create_flags_table($dbh,$flags_table,$gen_feats);

  _create_types_table($dbh,$flags_table,$types_table,$species);

  system("rm -r -f ".$dumpdir) && throw "Error deleting $dumpdir folder";

  
}

sub _create_gen_feats_table{
    my($dbh,$gen_feats_aref,$table)=@_;

    my @sql;
    # create a temporary table of gen_feat,id
    my $temp1 = "tmp_".int(rand(100))."_1";;
    push @sql,"drop table if exists $temp1";
    push @sql,"create table $temp1 (feature_type varchar(60),feature_type_id int(10) unsigned)";
    for(my $i=0;$i < @$gen_feats_aref;$i++){
	push @sql, "insert into $temp1 values( '".$gen_feats_aref->[$i]."','".$i."')";
    }
    push @sql,"alter table $temp1 add index(feature_type)";

    # concatenate all the gen feats into one table
    my $temp2 = "tmp_".int(rand(100))."_2";;
    push @sql,"drop table if exists $temp2";
    push @sql, "create table $temp2 (feature_type varchar(60), seq_region_name varchar(40), feature_strand tinyint(2), feature_start int(10) unsigned, feature_end int(10) unsigned)";
    foreach my $feat (@$gen_feats_aref){
	push @sql, "insert into $temp2 select distinct feature_type, seq_region_name, feature_strand, feature_start, feature_end from $feat";
    }

    push @sql,"alter table $temp2 add index(feature_type)";

    push @sql,"drop table if exists $table";
    push @sql, "create table $table select t2.*,t1.feature_type_id from $temp1 t1, $temp2 t2 where t1.feature_type = t2.feature_type";
    push @sql,"alter table $table add index(seq_region_name)";

    push @sql,"drop table if exists $temp1";
    push @sql,"drop table if exists $temp2";

    _run_sql($dbh,@sql) or throw "Error running sql";


}


sub _overlap_analysis{
    my($dbh,$gen_file,$chr,$reg_feat_table,$ofh,$columns)=@_;

    my $zeros = '0' x $columns;
    my @zero = split('',$zeros);

    my $q = "select regulatory_feature_id,feature_start,feature_end from $reg_feat_table where seq_region_name = '$chr' order by feature_start,feature_end";
    my $reg_aaref = $dbh->selectall_arrayref($q);
    unless(defined $reg_aaref){ throw "failed on \n $q\n".$dbh->errstr; }

    # regulatory features may not be present on all chrs
    if($reg_aaref eq NO_ROWS){return}
    #if($reg_aaref == 0 ){return}

    # create a hash to store the overlap flags
    my %flags;
    foreach my $reg_ref (@$reg_aaref){
        my @zero = split('',$zeros);
	$flags{$reg_ref->[0]} = \@zero;
    }

    open(IN,$gen_file) or throw "failed to open $gen_file";

    while(<IN>){
	chop;
	my($type,$start,$end) = split("\t",$_);

        foreach  my $reg_ref (@$reg_aaref){

            # if the reg and gen feats overlap set the flag for that gen feat
            if($reg_ref->[1] <= $end && $reg_ref->[2] >= $start){
                $flags{$reg_ref->[0]}->[$type] += 1;
	    }

            # trim the regulatory feature array
            # remove any reg feat whose end < this gen feat start
	    if($reg_ref->[2] < $start){shift(@$reg_aaref)}

            # if the start of the reg feat is > end of the gen feat
            # we can move onto the next gen feat
            if($reg_ref->[1] > $end){last}
	}

    }

    while(my($id,$aref)=each(%flags)){
	print $ofh $id."\t".join("\t",@$aref)."\n";
    }
    
}


sub _load_results{
    my($dbh,$results_file,$overlaps_table,$reg_feat_table,$gen_feat_ref)=@_;

    my @sql;

    #create the summary table
    push @sql,"drop table if exists $overlaps_table";
    my $q = "create table $overlaps_table (regulatory_feature_id int(10) unsigned, ";
    $q .= join(' int(1) unsigned, ',@$gen_feat_ref);
    $q .= " int(1) unsigned)";
    push @sql, $q;

    push @sql, "load data local infile '$results_file' into table $overlaps_table";
    push @sql, "alter table $overlaps_table add index(regulatory_feature_id)";

    my $temp1 = "tmp_".int(rand(100));
    push @sql,"drop table if exists $temp1";

    push @sql,"create table $temp1 select o.*,r.binary_string as pattern from $reg_feat_table r, $overlaps_table o where o.regulatory_feature_id = r.regulatory_feature_id";
    push @sql,"drop table $overlaps_table";
    push @sql,"alter table $temp1 rename as $overlaps_table ";
    
    _run_sql($dbh,@sql) or throw "cannot run sql";
}


sub _sub_pats{
    my($dbh,$reg_feat_table,$pat_bits)=@_;

    my $aref = $dbh->selectcol_arrayref("select binary_string from regulatory_features_filtered group by binary_string having count(*)  > 1 ");
    unless(defined $aref && @$aref >0){ die "no binary strings found with query :- \n select binary_string from regulatory_features_filtered group by binary_string having count(*)  > 1"}

    my %pats;
    foreach my $orig (@$aref){
        my $zero_count = $orig =~ tr/0/_/;
        #if($zero_count >= $pat_bits-2){next}#dont want pats with only one mark
        # following the introduction of projected builds and the use of TFBS 
        # as focus features, patterns with only one mark are not necessarily
        # just a focus feature. we therefore want to analyse them
        if($zero_count == $pat_bits){next}
        if($zero_count == $pat_bits - 1){ # just one bit set
            $pats{$orig} = 1;
            next;
        }
 
        my @ch = split('',$orig);
        # we create the set of all patterns which comprise the original but 
        # with one of the set bits made into '_' 
        for(my $i = 0;$i<@ch;$i++){
            if($ch[$i] eq '1'){
                $ch[$i]='_';
                $pats{(join('',@ch))} = 1;
                $ch[$i] = 1;

            }
            $pats{$orig} = 1;
        }
    }

    return keys(%pats);

}

sub _create_pattern_results_tables{
    my($dbh, $overlap_results_file,$gen_feats_ref)=@_;

    my @sql;
    push @sql,"drop table if exists raw_overlap_results";

    my $bitsize=255;
    my $q = "create table raw_overlap_results (pattern varchar(${bitsize}),total_reg_feats int(10)";
    foreach my $gen (@$gen_feats_ref){
	$q .= ",$gen"."_real int(10), $gen"."_mock int(10)";
    }
    $q .= ")";
    push @sql,$q;    
    #Similar as mysqlimport...
    push @sql, "load data local infile '$overlap_results_file' into table raw_overlap_results";

    push @sql,"drop table if exists pattern_overlap_summary";
    push @sql,"drop table if exists pattern_overlap_summary_chi";
    $q = "create table pattern_overlap_summary (pattern varchar(${bitsize})";
    my $q1 = "create table pattern_overlap_summary_chi (pattern varchar(${bitsize})"; 
    my $q2 = "insert into pattern_overlap_summary select pattern";
    my $q3 = "insert into pattern_overlap_summary_chi select pattern";
    foreach my $gen (@$gen_feats_ref){
        # using int(10) to remove non integer part of %age 
	$q .= ",$gen int(10)";
	$q1 .= ",$gen float";
        $q2 .= ",100* $gen".'_real/total_reg_feats '."as $gen";
        #$q3 .= ",($gen"."_real- $gen"."_mock)*($gen"."_real- $gen"."_mock)/ $gen"."_mock as $gen"
        $q3 .= ",($gen"."_real- $gen"."_mock)*($gen"."_real- $gen"."_mock)/ if($gen"."_mock,$gen"."_mock,1) as $gen"
    }
    $q .= ")";
    $q1 .= ")";
    $q2 .= " from raw_overlap_results";
    $q3 .= " from raw_overlap_results";
    push @sql,$q;   
    push @sql,$q1;   
    push @sql,$q2;
    push @sql,$q3;
    push @sql,"alter table pattern_overlap_summary_chi add index(pattern)";
    push @sql,"alter table pattern_overlap_summary add index(pattern)";

    _run_sql($dbh,@sql) or die;

}

sub _association_tables{
    my($dbh, $gen_feat_ref,$summary_table,$chi_table,$assoc_thresh,$second_thresh)=@_;
    
    
    foreach my $feat (@$gen_feat_ref){
	# scan for initial associations
	my $q = "select perc.pattern, perc.$feat from $summary_table perc, $chi_table chi where perc.$feat >= $assoc_thresh and chi.$feat > 8 and perc.pattern = chi.pattern";

	my $aaref = $dbh->selectall_arrayref($q);
	unless(defined $aaref){ throw "failed on \n$q\n".$dbh->errstr}
	foreach my $aref (@$aaref){
	    push @$aref,$feat;
	}

	my @sql;
	my $assoc_table = $feat."_assoc";
	my $not_assoc_table =  $feat."_not_assoc";
	my $zero_one_table = $feat."_0_1";
	push @sql,"drop table if exists $assoc_table";
	push @sql,"create table $assoc_table select * from $summary_table where 1 = 0";

	push @sql,"drop table if exists $not_assoc_table";
	push @sql,"create table $not_assoc_table select * from $summary_table where 1 = 0";
	push @sql,"drop table if exists $zero_one_table";
	push @sql,"create table $zero_one_table select pattern from $summary_table where 1 = 0";
	_run_sql($dbh,@sql) or throw "error running sql";

	foreach my $aref (@$aaref){
	    my $pattern = $aref->[0];

	    _run_sql($dbh, "insert into $assoc_table select * from $summary_table where pattern like '$pattern'") or throw "error running sql";


	    # if none of the patterns have percent < $second_thresh then the original pattern
	    # is OK and it stays in the table.  we remove all the
	    # other patterns which have percent < $assoc_thresh
	    $q = "select pattern from $assoc_table where  $feat < $second_thresh";
	    my $not_ref = $dbh->selectcol_arrayref($q);
	    die $q."\n".$dbh->errstr unless(defined $not_ref);

	    if(@$not_ref > 0){
		foreach my $pat (@$not_ref){
		    _run_sql($dbh,"insert into $not_assoc_table select * from $assoc_table where  pattern = '$pat'") or throw "error running sql";
		    _run_sql($dbh,"delete from $assoc_table where pattern = '$pat'") or throw "error running sql";
		      
		    my $zero_one = _add_nots($pattern,$pat);
		    _run_sql($dbh,"insert into $zero_one_table values('$zero_one')") or throw "error running sql";
		}
	    }


	    @sql=();
	    push @sql,"drop table if exists temp_$$";
	    push @sql,"create table temp_$$ select distinct * from $assoc_table";
	    push @sql,"drop table $assoc_table";
	    push @sql,"alter table temp_$$ rename as $assoc_table";
            push @sql,"alter table $assoc_table add index(pattern)";
            _run_sql($dbh,@sql);
	}

	@sql = ();
        # the assoc table currently contains patterns which lie between
        # the assoc_thresh and the second_thresh
	push @sql,"delete from $assoc_table where $feat < $assoc_thresh";
	# distinct assoc table
	# **** and add back the two gene bits to the patterns ****
	push @sql,"drop table if exists temp_$$";
	push @sql,"create table  temp_$$ select distinct * from $assoc_table";
	push @sql,"drop table $assoc_table";
	push @sql,"alter table temp_$$ rename as $assoc_table";
	#push @sql,"update $assoc_table set pattern = concat(pattern,'__')";
	push @sql,"alter table $assoc_table add index(pattern)";

	# distinct not_assoc table
	# **** and add the two gene bits to the patterns ****
	# distinct not_assoc table
	# **** and add the two gene bits to the patterns ****
	push @sql,"drop table if exists temp_$$";
	push @sql,"create table  temp_$$ select distinct * from $not_assoc_table";
	push @sql,"drop table $not_assoc_table";
	push @sql,"alter table temp_$$ rename as $not_assoc_table";
	#push @sql,"update $not_assoc_table set pattern = concat(pattern,'__')";
	push @sql,"alter table $assoc_table add index(pattern)";

	_run_sql($dbh,@sql) or die;
    }

}

sub _add_nots{
    my($pat,$not) = @_;

    my @pat = split('',$pat);
    my @not = split('',$not);

    for(my $i = 0;$i<@pat;$i++){
        if($pat[$i] ne '1' && $not[$i] eq '1'){
            $pat[$i] = 0;
        }
    }
    return join('',@pat);
}

sub _create_flags_table{
    my($dbh,$flags_table,$gen_feat_ref)=@_;


    my @sql;
    push @sql,"drop table if exists $flags_table";
    push @sql, "create table $flags_table select regulatory_feature_id,binary_string from regulatory_feature";
    push @sql, "alter table $flags_table add index(binary_string)";

    _run_sql($dbh,@sql) or throw "error running sql";
 
   
    @sql = ();
    foreach my $feat (@$gen_feat_ref){

        my $assoc_table = $feat."_assoc";
        my $not_assoc_table =  $feat."_not_assoc";

	#what do we want here? check if the assoc table is empty or see if it exists?
	# hack covering specific cases where there are no associations!!!
	unless( ($feat eq 'intergenic_2500')
		|| ($feat eq 'protein_coding_gene_body')
		|| ($feat eq 'protein_coding_exon1_plus_enhancer')
		|| ($feat eq 'protein_coding_intron1')
		|| ($feat eq 'tss_centred_5000')
		|| ($feat eq 'PolIII_transcribed_gene_plus_enhancer')
		|| ($dbh->selectrow_array("select count(*) from $assoc_table") > 0)){
	    next;
	}

	push @sql, "alter table $flags_table add column $feat int(1) default 0";

	my $q = "select pattern from $assoc_table";
        my $pat_ref = $dbh->selectcol_arrayref($q);
	unless(defined $pat_ref){die "failed on \n$q\n".$dbh->errstr}
        foreach my $pat (@$pat_ref){
	    push @sql, "update $flags_table set $feat = 1 where binary_string like '$pat'";

	}

	$q = "select pattern from $not_assoc_table";
        $pat_ref = $dbh->selectcol_arrayref($q);
	unless(defined $pat_ref){die "failed on \n$q\n".$dbh->errstr}
        foreach my $pat (@$pat_ref){
	    push @sql, "update $flags_table set $feat = 0 where binary_string like '$pat'";
	}
    }

    # hack covering specific cases where there are no associations!!!
#    unless(_column_exists($dbh,$flags_table,'intergenic_2500')){
#    	push @sql, "alter table $flags_table add column intergenic_2500 int(1) default 0";
#    }
#    unless(_column_exists($dbh,$flags_table,'protein_coding_gene_body')){
#	push @sql, "alter table $flags_table add column protein_coding_gene_body int(1) default 0";
#    }
#    unless(_column_exists($dbh,$flags_table,'protein_coding_exon1_plus_enhancer')){
#	push @sql, "alter table $flags_table add column protein_coding_exon1_plus_enhancer int(1) default 0";
#    }
#    unless(_column_exists($dbh,$flags_table,'protein_coding_intron1')){
#	push @sql, "alter table $flags_table add column protein_coding_intron1 int(1) default 0";
#    }
#    unless(_column_exists($dbh,$flags_table,'tss_centred_5000')){
#	push @sql, "alter table $flags_table add column tss_centred_5000  int(1) default 0";
#    }
#    unless(_column_exists($dbh,$flags_table,'PolIII_transcribed_gene_plus_enhancer')){
#      push @sql, "alter table $flags_table add column PolIII_transcribed_gene_plus_enhancer  int(1) default 0";
#    }

    _run_sql($dbh,@sql) or throw "error running sql";

}


sub _column_exists {

  my $dbh = shift;
  my $table = shift;
  my $col = shift;
  
  my $query = "select $col from $table limit 1";
  my $sth;
  
  unless($sth=$dbh->prepare($query)){
    throw("ERROR: preparation of statement $query failed");
  }
  
  unless($sth->execute){
    $sth->finish;
    return(0);
  }

  $sth->finish;
  return(1);
  
}

sub _create_types_table{
    my($dbh,$flags_table,$types_table,$species)=@_;

    my @sql;

    push @sql, "drop table if exists $types_table";

    # The same for all species...
    #if($species eq 'homo_sapiens'){
      push @sql, "create table $types_table select regulatory_feature_id,binary_string,0 as cell_type_specific,tss_centred_5000 as promoter_associated,protein_coding_gene_body as gene_associated, intergenic_2500 as non_gene_associated,0 as unclassified, PolIII_transcribed_gene_plus_enhancer as poliii_transcription_associated from $flags_table";
    #}else{
    #    push @sql, "create table $types_table select regulatory_feature_id,binary_string,0 as cell_type_specific,if(protein_coding_exon1_plus_enhancer+protein_coding_intron1 > 0, 1,0) as promoter_associated,protein_coding_gene_body as gene_associated, intergenic_2500 as non_gene_associated,0 as unclassified from $flags_table";
    #}

    # apply arbitrary rules to resolve conflicts
    push @sql, "update $types_table set promoter_associated = 0 where promoter_associated and gene_associated"; # as used for v58
    #push @sql, "update $types_table set gene_associated = 0 where promoter_associated and gene_associated"; #

    push @sql, "update $types_table set gene_associated = 0 where poliii_transcription_associated and gene_associated";
    push @sql, "update $types_table set non_gene_associated = 0 where poliii_transcription_associated and non_gene_associated";
    push @sql, "update $types_table set poliii_transcription_associated = 0 where poliii_transcription_associated and promoter_associated";

    # use unclassified col to flag conflicts
    push @sql, "update $types_table set unclassified = 1 where promoter_associated and non_gene_associated";
    push @sql, "update $types_table set unclassified = 1 where gene_associated and non_gene_associated";
    # set both of the conflicting cols to 0 where there is a conflict
    push @sql, "update $types_table set promoter_associated = 0 where unclassified";
    push @sql, "update $types_table set gene_associated = 0 where unclassified";    push @sql, "update $types_table set non_gene_associated = 0 where unclassified";

    # set unclassified for rows with no flags
    push @sql, "update $types_table set unclassified = 1 where gene_associated + non_gene_associated + promoter_associated+poliii_transcription_associated = 0";

    _run_sql($dbh,@sql) or throw "error running sql";

    # add the feature_type_id column

    _run_sql($dbh, "alter table $types_table add column feature_type_id int(10) unsigned") or throw "could not add feature_type_id column in $types_table";
    # set the type_id 
    # now that we have moved to single cell line classification there are
    # no Cell type specific classifications
    @sql = ();
    foreach my $ft ('Gene Associated',
                    'Non-Gene Associated',
                    'Promoter Associated',
                    'Unclassified',
                    'PolIII Transcription Associated',
		    ){

      my $ftid = $dbh->selectrow_array("select feature_type_id from feature_type where name = '$ft' and class = 'Regulatory Feature'");
      if(!defined($ftid) || ($ftid <1)){ throw "feature_type_id for $ft does not exist or value is less than 1"}
      my $cts = 0;
      my $col;
      if($ft =~ 'specific'){
	$cts = 1;
	($col) = $ft =~ /(.+ *.+) -.*/;
	$col = lc($col);
	$col =~ tr/ /_/;
	$col =~ tr/-/_/;
      }else{
	$col = lc($ft);
	$col =~ tr/ /_/;
	$col =~ tr/-/_/;
      }
      print "ft $ft col $col\n";
      push @sql, "update $types_table set feature_type_id = $ftid where cell_type_specific = $cts and $col = 1";
      
    }
    _run_sql($dbh,@sql) or throw "cannot run sql";

    # summary report
    # remove and not cell_type_specific... 
    my $res = $dbh->selectrow_array(" select count(*) from  regulatory_features_classified where promoter_associated and not cell_type_specific ");
    _commentary("promoter_associated     $res\n");
    $res = $dbh->selectrow_array(" select count(*) from  regulatory_features_classified where gene_associated and not cell_type_specific");
    _commentary("gene_associated         $res\n");
    $res = $dbh->selectrow_array(" select count(*) from  regulatory_features_classified where non_gene_associated and not cell_type_specific");
    _commentary("non_gene_associated     $res\n");
    $res = $dbh->selectrow_array(" select count(*) from  regulatory_features_classified where poliii_transcription_associated and not cell_type_specific");
    _commentary("PolIII_transcription_associated     $res\n");

    $res = $dbh->selectrow_array(" select count(*) from  regulatory_features_classified where unclassified and not cell_type_specific");
    _commentary("unclassified            $res\n");

    _qc($dbh,);



}


sub _qc{
    my($dbh) = @_;

    my @sql;
    #goto GENE_ASSOC;
    push @sql, "drop table if exists promoter_associated_temp";
    push @sql, "create table promoter_associated_temp select rfc.regulatory_feature_id, seq_region_name,seq_region_start,seq_region_end,rfc.binary_string from regulatory_feature rf, regulatory_features_classified rfc where rfc.promoter_associated and rfc.regulatory_feature_id =rf.regulatory_feature_id order by seq_region_name;";
    push @sql, "alter table promoter_associated_temp add index(seq_region_name)";
    push @sql, "drop table if exists promoter_features";
    push @sql, "create table promoter_features select * from protein_coding_exon1_plus_enhancer ";
    #push @sql, "insert into promoter_features select * from  protein_coding_intron1";
    push @sql, "insert into promoter_features select * from RNA_gene_exon1_plus_enhancer";
    push @sql, " alter table promoter_features add index(seq_region_name)";

    _run_sql($dbh,@sql) or throw "cannot run sql";

    my $res;
    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from promoter_associated_temp pa, promoter_features e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= e.feature_end and pa.seq_region_end >= e.feature_start"); 
    _commentary("promoter_associated features\n");
    _commentary("$res overlap an exon1_plus_2.5kb (both RNA and prot_cod)\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from promoter_associated_temp pa, protein_coding_intron1 e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= e.feature_end and pa.seq_region_end >= e.feature_start");   
    _commentary("$res overlap a protein coding intron1\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from promoter_associated_temp pa, processed_transcript e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= if(e.feature_strand = -1 ,e.feature_end+2500,e.feature_end)  and pa.seq_region_end >= if(e.feature_strand = -1 ,e.feature_start,e.feature_start - 2500)");
    _commentary("$res overlap a 'processed transcript'\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from promoter_associated_temp pa, protein_coding_gene_body e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= e.feature_end and pa.seq_region_end >= e.feature_start");
    _commentary("$res overlap a protein coding gene body\n");
    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from promoter_associated_temp pa, intergenic_2500 e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start > e.feature_start and pa.seq_region_end < e.feature_end");
    _commentary("$res are at least 2500bp from any part of a protein coding gene\n");
    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from promoter_associated_temp pa, pseudogene_exon1_plus_enhancer e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start > e.feature_start and pa.seq_region_end < e.feature_end");
    _commentary("$res overlap a pseudogene exon1_plus_2.5kb \n");

    @sql=();
    push @sql, "drop table if exists gene_associated_temp";
    push @sql, "create table gene_associated_temp select rfc.regulatory_feature_id, seq_region_name,seq_region_start,seq_region_end from regulatory_feature rf, regulatory_features_classified rfc where rfc.gene_associated and rfc.regulatory_feature_id =rf.regulatory_feature_id order by seq_region_name";
    push @sql, "alter table gene_associated_temp add index(seq_region_name)";
    _run_sql($dbh,@sql) or throw "cannot run sql";

    _commentary("gene_associated features\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from gene_associated_temp pa, protein_coding_gene e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= e.feature_end and pa.seq_region_end >= e.feature_start");
    _commentary("$res overlap some part of a protein coding gene\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from gene_associated_temp pa, RNA_gene e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= e.feature_end and pa.seq_region_end >= e.feature_start");
    _commentary("$res overlap some part of an RNA gene\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from gene_associated_temp pa, intergenic_2500 e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start > e.feature_start and pa.seq_region_end < e.feature_end");
    _commentary("$res are at least 2500bp from any part of a protein coding gene\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from gene_associated_temp pa, pseudogene_transcript e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= e.feature_end and pa.seq_region_end >= e.feature_start");
    _commentary("$res overlap some part of a pseudogene\n");

    $res = $dbh->selectrow_array("select count( distinct pa.regulatory_feature_id ) from gene_associated_temp pa, processed_transcript e where pa.seq_region_name=e.seq_region_name and pa.seq_region_start <= e.feature_end and pa.seq_region_end >= e.feature_start");
    _commentary("$res overlap some part of a 'processed_transcript'\n");
}

sub _commentary{
  my ($comm) = @_;
  print $comm;
  print REPORT $comm;
}


# returns start and end 
sub _get_rand_region{
    my ($dbh,$enc,$len,%enc_len) = @_;

    my $bad = 1;
    my $end;
    my $start;
    my $try = 1;
    while($bad){
      $try++;
      if($try > 1000000){
	throw "can't find a random region for seq_region $enc, length $len\n";
      }
      
      # get a random number which is within the encode region
      $end = int(rand($enc_len{$enc}));
      
      # the number must be higher than the length of the frag
      if($end > $len){
	$bad = 0;
      }
      $start = $end - $len +1;
      
      if($bad == 0){
	#Avoid regions that were already taken...
	my $q = "select * from banned_region_table where seq_region_name = '$enc' and $start <= feature_end and $end >= feature_start limit 1 ";
	my $aaref = $dbh->selectall_arrayref($q);
	unless(defined $aaref){die $dbh->errstr} 
	if(scalar(@$aaref) > 0){  $bad = 1; } 
      }
    }
    
    return ($start,$end);

}


#Help function to copy a table between databases with different names...
#Require DBI db handle objects
sub _copy_with_rename{
    my($s_dbh,$source_table,$t_dbh,$targ_table)=@_;

    my $q = "desc $source_table";
    my $desc_aaref=$s_dbh->selectall_arrayref($q);
    unless(defined $desc_aaref){die "failed on:\n $q\n".$s_dbh->errstr}

    _run_sql($t_dbh,"drop table if exists $targ_table");
    $q = " create table $targ_table (";
    my $select;
    foreach my $aref (@$desc_aaref){
	$q  .= $aref->[0].' '.$aref->[1].',';
        $select .=  $aref->[0].',';
    }
    chop $q; #remove last comma
    chop $select;
    $q .= ")";
    print $q."\n";
    _run_sql($t_dbh,$q) or throw "Error creating table";

    $q = "select $select from $source_table";
    my $aaref = $s_dbh->selectall_arrayref($q);
    unless(defined $aaref){die "failed on:\n $q\n".$s_dbh->errstr}
    my @sql;
    foreach my $aref (@$aaref){
	$q= "insert into $targ_table values('".join("','",@$aref);
	$q .= "')";
	push @sql,$q;

    }

    _run_sql($t_dbh,@sql) or throw "Error copying data";

}


# _run_sql
#
#  Arg [1]   : scalar database handle from DBI module
#  Arg [2]   : array of scalars containing lines of text in the format of SQL statements
#  Function  : submits each element of Arg[2] to the $dbh->prepare and 
#              $dbh->execute methods and checks for successful execution. 
#              Returns 0 on failure, 1 on success. Emits error messages.
#  Returntype: int
#  Exceptions: none
#  Example   : _run_sql($dbh, @array);

sub _run_sql{
  my $dbh = shift;
  my (@array)=@_;
  
  my $sth;
  
  foreach my $query(@array){
    
    eval {
      unless($sth=$dbh->prepare($query)){
	warn("Error: preparation of statement failed on: $query\n");
	warn("database error message: ".$dbh->errstr."\n");
	return(0);
      }
      
      unless($sth->execute){ # returns true on success
	warn("Error: statement execution failed on: $query\n");
	warn("statement handle error message:".$sth->errstr."\n");
	return(0);
      }
    };
    if($@){ warn "Error running SQL"; return(0); }

  }
  
  return(1);
}

#Private getter / setter to the Work DB Connection
sub _workdb {
  return $_[0]->_getter_setter('workdb',$_[1]);
}

#Private getter / setter to the cell type
sub _cell_type {
  return $_[0]->_getter_setter('cell_type',$_[1]);
}



1;
