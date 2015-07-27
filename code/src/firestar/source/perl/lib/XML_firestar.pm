
package XML_firestar;

use strict;
use XML::DOM;
use DBI_firestar;


sub new {
	my($class,%args)=@_;
	my $this={};
	$this->{doc}=XML::DOM::Document->new();
	$this->{root}=$this->{doc}->createElement($args{-root});
	bless($this);
	return $this;
}

sub print {
	my($this)=@_;
	my $out=$this->{root}->toString();
	$out=~s/></>\n</g;
	return ($out);

}

sub app{
	my($this,%args)=@_;
	if (exists $args{-father}){$args{-father}->appendChild($args{-child});}
	else{
		$this->{root}->appendChild($args{-child});
	}
}

sub tex{
	my($this,%args)=@_;
        $args{-tag}->addText($args{-text});
}

sub crea{
	my($this,%args)=@_;
	return $this->{doc}->createElement($args{-tag});

}

sub att{
	my($this,%args)=@_;
        $args{-tag}->setAttribute($args{-att},$args{-val});
}

sub no_CSA{
	my($this)=@_;
	my $fail=$this->crea(-tag=>'NO_catalytic_sites');
	$this->tex(-tag=>$fail,-text=>'There are no CSA predicted for your target');
	$this->app(-child=>$fail);
}

sub no_CCT{
	my($this)=@_;
	my $fail=$this->crea(-tag=>'NO_binding_sites');
	$this->tex(-tag=>$fail,-text=>'There are no binding sites predicted for your target');
	$this->app(-child=>$fail);
}


sub CSA_prediction{
	my ($this,$CSA,$ecs,$pdbs,$orden,$lib,$join,@homologs)=@_;
	my $csa_pred=$this->crea(-tag=>'predicted_catalytic_site');
	$this->app(-father=>$CSA,-child=>$csa_pred);
	my @ecs=split(/ /,$ecs);
	my $ec_list=$this->crea(-tag=>'ec_numbers_list');
	$this->app(-father=>$csa_pred,-child=>$ec_list);
	foreach my $i(@ecs){
		my $ec=$this->crea(-tag=>'EC');
		$this->att(-tag=>$ec,-att=>'number',-val=>$i);
		$this->app(-father=>$ec_list,-child=>$ec);
		if (exists $lib->{ec2go}{$i}){
			my @tmp=@{$lib->{ec2go_full}{$i}};
			foreach my $k(@tmp){
				if($k=~/^EC:(.+)\s>\sGO:(.+); (GO\:\d+)/){
					my $GO=$this->crea(-tag=>'associated_GO');
					$this->att(-tag=>$GO,-att=>'number',-val=>$3);
					$this->tex(-tag=>$GO,-text=>$2);
					$this->app(-father=>$ec,-child=>$GO);
				}
			}
		}

	}
	my $evidence=$this->crea(-tag=>'evidence');
	$this->att(-tag=>$evidence,-att=>'type_of_source',-val=>"literature");
	$this->app(-father=>$csa_pred,-child=>$evidence);
	my $entries=$this->crea(-tag=>'CSA_source_entries');
	$this->app(-father=>$csa_pred,-child=>$entries);
	my @pdbs=split(/ /,$pdbs);
	foreach my $i(@pdbs){
		my $entry=$this->crea(-tag=>'homologous_entry');
		$this->att(-tag=>$entry,-att=>'PDB_id',-val=>$i);
		$this->app(-father=>$entries,-child=>$entry);
	}
	my @residues=split(/ /,$join->{summary_csa}{resnum}[$orden]);
	my @code=split(//,$join->{summary_csa}{resname}[$orden]);
	my $res=$this->crea(-tag=>'predicted_residues_list');
	$this->app(-father=>$csa_pred,-child=>$res);
	for (my $i=0;$i<scalar@residues;$i++){
		my $aa=$this->crea(-tag=>'residue');
		$this->att(-tag=>$aa,-att=>'position',-val=>$residues[$i]);
		$this->att(-tag=>$aa,-att=>'one_letter_code',-val=>$code[$i]);
		$this->app(-father=>$res,-child=>$aa);
	}
	my @score=split(//,$join->{summary_csa}{score}[$orden]);
	my $sum;
	foreach my $i(@score){$sum+=$i;}
	my $sco=sprintf("%.2f",$sum/scalar@score);
	my $score=$this->crea(-tag=>'site_score');
	$this->att(-tag=>$score,-att=>'conservation_score',-val=>$sco);
	$this->att(-tag=>$score,-att=>'max',-val=>"6");
	$this->app(-father=>$csa_pred,-child=>$score);
}


sub CCT_prediction{
	my ($this,$CCT,$GOs,$GOnames,$order,$join,$fire)=@_;
	my $pred=$this->crea(-tag=>'predicted_binding_site');
	$this->app(-father=>$CCT,-child=>$pred);
	my $comp=$this->crea(-tag=>'compound');
	$this->att(-tag=>$comp,-att=>'three_letter_code',-val=>$join->{summary}{compid_code}[$order]);
	$this->tex(-tag=>$comp,-text=>$join->{summary}{compname}[$order]);
	$this->app(-father=>$pred,-child=>$comp);
	my $class=$this->crea(-tag=>'site_classification');
	$this->att(-tag=>$class,-att=>'type',-val=>$join->{summary}{nice_try}[$order]);
	$this->app(-father=>$pred,-child=>$class);
	my $comps_list=$this->crea(-tag=>'templates_compounds_list');
	$this->app(-father=>$pred,-child=>$comps_list);
	my @ligs=split(/ /,$join->{summary}{compid}[$order]);
	$this->{dbi}=DBI_firestar->connect($fire->{config_file});
	foreach my $i(@ligs){
		if ($i=~/(\w+)\((\d+)\)/){
			my $compo=$this->crea(-tag=>'ligand');
			$this->att(-tag=>$compo,-att=>'three_letter_code',-val=>$1);
			$this->att(-tag=>$compo,-att=>'frequency',-val=>$2);
			my $answer=$this->{dbi}->pharma(-compid=>$1);
			if ($answer eq "YES"){$this->att(-tag=>$compo,-att=>'pharmacological_annotation',-val=>'YES');}
			$this->app(-father=>$comps_list,-child=>$compo);
		}
	}
	$this->{dbi}->close();
	my @GOs=split(/ /,$GOs);
	my @GOnames=split(", ",$GOnames);
	my $go_list=$this->crea(-tag=>'associated_GOs_list');
	$this->app(-father=>$pred,-child=>$go_list);
	for (my $i=0; $i<scalar@GOs; $i++){
		my $go=$this->crea(-tag=>'associated_GO');
		$this->att(-tag=>$go,-att=>'term',-val=>$GOs[$i]);
		$this->tex(-tag=>$go,-text=>$GOnames[$i]);
		$this->app(-father=>$go_list,-child=>$go);
	}
	my @residues=split(/;/,$join->{summary}{resfreq}[$order]);
	my $res=$this->crea(-tag=>'predicted_residues_list');
	$this->app(-father=>$pred,-child=>$res);
	foreach my $i(@residues){
		if ($i=~/(\w)\((\d+)\)=(\d\.\d+)/){
			my $aa=$this->crea(-tag=>'residue');
			$this->att(-tag=>$aa,-att=>'position',-val=>$2);
			$this->att(-tag=>$aa,-att=>'one_letter_code',-val=>$1);
			$this->att(-tag=>$aa,-att=>'relative_frequency',-val=>$3);
			$this->app(-father=>$res,-child=>$aa);
		}
	}
	my $sco=$this->crea(-tag=>'site_score');
	$this->app(-father=>$pred,-child=>$sco);
	if ($join->{summary}{score}[$order]=~/(.+%) \[COV: (.+%) SITE: (.+%) Ident: (.+%) - Ali: (.+%)\]/){
		$this->att(-tag=>$sco,-att=>'score',-val=>$1);
		my $par1=$this->crea(-tag=>'coverage');
		$this->att(-tag=>$par1,-att=>'value',-val=>$2);
		$this->app(-father=>$sco,-child=>$par1);
		my $par2=$this->crea(-tag=>'conservation_score');
		$this->att(-tag=>$par2,-att=>'value',-val=>$3);
		$this->app(-father=>$sco,-child=>$par2);
		my $par3=$this->crea(-tag=>'best_template_identities');
		$this->att(-tag=>$par3,-att=>'value',-val=>$4);
		$this->app(-father=>$sco,-child=>$par3);
		my $par4=$this->crea(-tag=>'alignments_merged');
		$this->att(-tag=>$par4,-att=>'number',-val=>$join->{summary}{numOfSites}[$order]);
		$this->att(-tag=>$par4,-att=>'fraction',-val=>$5);
		$this->app(-father=>$sco,-child=>$par4);
	}

	

}


1;
