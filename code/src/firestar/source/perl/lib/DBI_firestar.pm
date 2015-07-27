
package DBI_firestar;

use strict;
use DBI;
use FindBin;
use Config::IniFiles;

sub connect {
	my ($class,$config)=@_;
	my $variables=Config::IniFiles->new(-file => $config);
	my($class)=@_;
	my $dbHost=$variables->val('MYSQL','dbHost');
	my $dbName=$variables->val('MYSQL','dbName');
	my $user=$variables->val('MYSQL','user');
	my $pass=$variables->val('MYSQL','pass');
	my $this->{dbi}=DBI->connect("DBI:mysql:database=$dbName;host=$dbHost","$user","$pass",{RaiseError => 1,AutoCommit =>0})
			|| die "Unable to connect to $dbName because $DBI::errstr";
	my $serverInfo = $this->{'mysql_serverinfo'};
	my $serverStat = $this->{'mysql_stat'};
	bless($this);
	return $this;
}

sub fetch{
	my($this,%args)=@_;
	my $query=$args{-query};
	my $sth=$this->{dbi}->prepare($query);
	$sth->execute();
	my @answer=$sth->fetchrow_array();
	$sth->finish();
	return @answer;
}

sub multi_fetch{
	my($this,%args)=@_;
	my $query=$args{-query};
	my $sth=$this->{dbi}->prepare($query);
	$sth->execute();
	my $contador=1;
	my %multi;
	while (my @answer=$sth->fetchrow_array()){
		@{$multi{$contador}}=@answer;
		$contador++;
	}
	$sth->finish();
	return %multi;
}

sub bio_tag{
	my($this,%args)=@_;
	my $query="select COMPID from COMPOUND where TAG=\"$args{-tag}\"";
	return $this->multi_fetch(-query=>$query);
}

sub met_tag{
	my($this,%args)=@_;
	my $query="select COMPID from COMPOUND where METAL_TAG=\"$args{-tag}\"";
	return $this->multi_fetch(-query=>$query);
}

sub extended_web_page_csites_info{
	my ($this,%args)=@_;
	my $query="select CSITEID,NUMCONRES,CONRES,OCCUPANCY,SITETYPE,EVIDENCY,COMPIDS,NUMSEQS,SITEIDS from CSITE35 where CLUSTID=\"$args{-clustid}\" and (SCORE=2 or SCORE=3)";
	return $this->multi_fetch(-query=>$query);
}

sub collapsed_sites{
	my ($this,%args)=@_;
	my $query="select CSITEID,NUMCONRES,CONRES,EVIDENCY,COMPIDS,LIGTYPE from CSITE35 where CLUSTID=\"$args{-clustid}\" and (SCORE=2 or SCORE=3)";
	return $this->multi_fetch(-query=>$query);
}

sub novel_inc_csites{
	my ($this,%args)=@_;
	my $query="select CSITEID,NUMCONRES,CONRES,EVIDENCY,COMPIDS,LIGTYPE,SCORE from CSITE35 where CLUSTID=\"$args{-clustid}\"";
	return $this->multi_fetch(-query=>$query);
}

sub occurrence{
	my ($this,%args)=@_;
	my $query="select OCCURRENCE from CCTEVAL_35 where CSITEID=\"$args{-csiteid}\"";
	my @call=$this->fetch(-query=>$query);
#	$query="select NUMSEQS from CSITE35 where CSITEID=\"$args{-csiteid}\"";
#	my @num=$this->fetch(-query=>$query);
#	print "OCCURRENCE $call[0] numseqs $num[0]\n";
#	return (($call[0]/$num[0])*100);
	return $call[0];
}


sub evo_info{
	my ($this,%args)=@_;
	my $query="select MOTIFNUMS1,MOTIF2,SQSCORES,PCENTID,COMPIDS2,OVERLAP,SIZE2 from COMPARE35 where CLUSTSITEID=\"$args{-csiteid}\"";
	return $this->multi_fetch(-query=>$query);
}

sub csanote{
	my ($this,%args)=@_;
	my $query="select CSANOTE from SITE35 where CADID=\"$args{-cad}\" and SITETYPE=\"CSA\"";
	return $this->fetch(-query=>$query);
}

sub compound_name{
	my ($this,%args)=@_;
	my $query="select NAME from COMPOUND where COMPID=\"$args{-compid}\"";
	my @name=$this->fetch(-query=>$query);
	return $name[0];
}

sub ligtype{
	my ($this,%args)=@_;
	my $query="select LIGTYPE from CSITE35 where CSITEID=\"$args{-csiteid}\" and (SCORE=2 or SCORE=3)";
	my @type=$this->fetch(-query=>$query);
	return $type[0];
}

sub enzyme_codes{
	my ($this,%args)=@_;
	my $query="select EC1,EC2,EC3 from INFOACC where CADID=\"$args{-clustid}\"";
	return $this->fetch(-query=>$query);
}

sub ping{
	my ($this)=@_;
	return($this->{dbi}->ping());
}

sub CSA_LIT{
	my ($this,%args)=@_;
	my $query="select SITEIDS,CLUSTID from CSITE35 where CSITEID=\"$args{-csite}\"";
	my @siteids=$this->fetch(-query=>$query);
	my @parking=split(/ /,$siteids[0]);
	$query="select CADID from CONSENSUS where CLUSTID=\"$siteids[1]\"";
	my %raw=$this->multi_fetch(-query=>$query);
	my %cadids;
	foreach my $i(keys%raw){	$cadids{${$raw{$i}}[0]}=5;}
	my %unicos;
	foreach my $k(@parking){
		$query="select CSANOTE from SITE35 where SITEID=\"$k\" and SITETYPE=\"CSA\"";
        	my @result=$this->fetch(-query=>$query);
		my @answer=grep(/$result[0]/,keys%cadids);
		if (scalar@answer>0){	$unicos{$result[0]}=5;}
	}
	return (keys%unicos);
}

sub evolutive_reliability{
	my ($this,%args)=@_;
	my $query="select CLUSTSITEID from COMPARE35 where CLUSTSITEID=\"$args{-csiteid}\"";
	my %related=$this->multi_fetch(-query=>$query);
	return(scalar(keys%related));
}

sub unique_chain_finder{
	my ($this,%args)=@_;
	my $query="select CADID from SITE35 where SITEID=\"$args{-siteid}\"";
	return $this->fetch(-query=>$query);
}

sub pdb2sequence{
	my ($this,%args)=@_;
	my $query="select PDBSEQ,CLUSTID from CONSENSUS where CADID=\"$args{-cadid}\"";
	return $this->fetch(-query=>$query);
}

sub chain_info{
	my ($this,%args)=@_;
	my $query="select PDBTITLE,UNIACC1,UNITITLE1 from INFOACC where CADID=\"$args{-cadid}\"";
	return $this->fetch(-query=>$query);
}

sub clust_codes{
	my ($this,%args)=@_;
	my $query="select CLUSTCODES from CONSENSUS where CADID=\"$args{-cadid}\"";
	my @codes=$this->fetch(-query=>$query);
	return $codes[0];
}

sub manually_annotated{
	my ($this,$lib)=@_;
	my $query="select ID_COMP,ID from MANUAL";
	my %park=$this->multi_fetch(-query=>$query);
	my %list;
	foreach my $i(keys%park){
		my $id=${$park{$i}}[0];
		$list{$lib->{id2pdbcode}{$id}}=${$park{$i}}[1];
	}
	return %list;
}

sub manual_annotation{
	my ($this,%args)=@_;
	my $query="select DESCRIPTION from MANUAL where ID=\"$args{-id}\"";
	my @annot=$this->fetch(-query=>$query);
	my $annotation="$annot[0] [";
	$query="select REF_ID from M2R where ID=\"$args{-id}\"";
	my %cross=$this->multi_fetch(-query=>$query);
	foreach my $i(keys%cross){
		$query="select INFO from REFER where REF_ID=\"${$cross{$i}}[0]\"";
		my @reference=$this->fetch(-query=>$query);
		if ($i==1){	$annotation.="$reference[0]";}
		else{		$annotation.=",$reference[0]";}
	}
	$annotation.="]";
	return $annotation;
}

sub pharma_annotated{
	my ($this,$lib)=@_;
	my %result;
	my $query="select ID_COMP,EXT_ID from MATCHING";
	my %park=$this->multi_fetch(-query=>$query);
	foreach my $i(keys%park){
		my $id=${$park{$i}}[0];
		if (exists $lib->{cognate}{$lib->{id2pdbcode}{$id}} or exists $lib->{poss_cognate}{$lib->{id2pdbcode}{$id}}){next;}
		$query="select PH_ID from A2M where EXT_ID=\"${$park{$i}}[1]\"";
		my %annot=$this->multi_fetch(-query=>$query);
		if (scalar(keys%annot)>0){push(@{$result{$lib->{id2pdbcode}{$id}}},${$park{$i}}[1]);}
	}
	return %result;
}

sub automatic_annotation{
	my ($this,%args)=@_;
	my @result;
	my $query="select CODE,DB,EXACT_MATCH from MATCHING where EXT_ID=\"$args{-ext_id}\"";
	my @auto=$this->fetch(-query=>$query);
	$query="select NAME from DA_BA where DB=\"$auto[1]\"";
	my @daba=$this->fetch(-query=>$query);
	push(@result,$daba[0],$auto[0]);
	if ($auto[2]==0){	push(@result,"DIFF");}
	else{	push(@result,"SAME");}
	$query="select PH_ID from A2M where EXT_ID=\"$args{-ext_id}\"";
	my %park=$this->multi_fetch(-query=>$query);
	foreach my $i(keys%park){
		my $ph_id=${$park{$i}}[0];
		$query="select INFO from PHARMA_ANNO where PH_ID=\"$ph_id\"";
		my @casual=$this->fetch(-query=>$query);
		push(@result,$casual[0]);
	}
	return @result;
}

sub pharma{
	my ($this,%args)=@_;
	my $flag="NO";
	my $query="select ID_COMP from COMPOUND where COMPID=\"$args{-compid}\"";
	my @park=$this->fetch(-query=>$query);
	my $id_comp=$park[0];
	$query="select EXT_ID from MATCHING where ID_COMP=\"$id_comp\"";
	my %park=$this->multi_fetch(-query=>$query);
	foreach my $i(keys%park){
		my $ext_id=${$park{$i}}[0];
		$query="select PH_ID from A2M where EXT_ID=\"$ext_id\"";
		my %park2=$this->multi_fetch(-query=>$query);
		if ($id_comp == 15788){print "EXT ID $ext_id number ",scalar(keys%park2),"\n";}
		if (scalar(keys%park2)>0){return("YES");}
	}
	$query="select ID from MANUAL where ID_COMP=\"$id_comp\"";
	my %park2=$this->multi_fetch(-query=>$query);
	if (scalar(keys%park2)>0){return("YES");}
	else{return("NO");}
}

sub id2pdb{
	my ($this)=@_;
	my %result;
	my $query="select ID_COMP,COMPID from COMPOUND";
	my %park=$this->multi_fetch(-query=>$query);
	foreach  my $i(keys%park){
		$result{${$park{$i}}[0]}=${$park{$i}}[1];
	}
	return %result;
}


sub close{
	my($this)=@_;
	$this->{dbi}->disconnect();
	return($this);
}

1;
