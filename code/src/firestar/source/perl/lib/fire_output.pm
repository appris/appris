package fire_output;

use strict;
use FindBin;
my $cwd;
BEGIN{
	$cwd=$FindBin::Bin;
}

use lib $cwd;
use File::Copy;
use DBI_firestar;
use fire_summary;
use XML_firestar;

sub new {
        my($class,%args)=@_;
        my $this={};
        bless($this);
        $this->{dbi}=DBI_firestar->connect($args{-config});
        return $this;
}

sub print_results{
	my ($this,$join,$fire,$lib,$pulser)=@_;
	########	ordena los sites según su reliability score (stored in %score_list)
	if (!exists $join->{summary_order} and $fire->{csa_option} eq "NO"){$this->no_results($fire,"YES","NO");return 0;}
	elsif(!exists $join->{csaincr} and $fire->{csa_option} eq "ONLY"){$this->no_results($fire,"NO","YES");return 0;}
	else{
		my @scoreList_unsorted;
		foreach (keys %{$join->{score_list}}){push(@scoreList_unsorted,$join->{score_list}{$_});}
		my @scoreList=sort {$b <=> $a} @scoreList_unsorted;
		foreach my $a (@scoreList){
			foreach my $i(keys %{$join->{score_list}}){
				if ($join->{score_list}{$i}==$a){push(@{$this->{order}},"$i ");delete($join->{score_list}{$i});}
			}
		}
		if ($fire->{output} eq 'appris')								{	$this->APPRIS_printer($join,$fire,$lib);}
		elsif($fire->{output} eq 'text' or $fire->{output} eq 'SIAM' or $fire->{output} eq 'XML')	{	$this->text_printer($join,$fire,$lib,$pulser);}
		elsif($fire->{output} eq 'web')									{	$this->web_printer($fire,$join,$lib);}
#	elsif($fire->{output} eq 'XML'){	$this->xml_printer($join,$fire);}
		return 0;
	}
}

sub no_results{
	my ($this,$fire,$csa,$pred)=@_;
	my $OUTPUT;
	$OUTPUT = \*STDOUT  unless(open($OUTPUT,'>',$fire->{outfile}));
	if ($fire->{output} eq 'appris'){
		print $OUTPUT "######\n";
		#print $OUTPUT ">>>$fire->{queryname}\t0 predicted residues\n";
		print $OUTPUT ">>>Query\t0 predicted residues\n";
		close($OUTPUT);
		return 0;
	}
	if ($fire->{output} eq 'XML'){
		my $xml_fire=XML_firestar->new(-root=>'Prediction');
		print $OUTPUT "Content-type: text/xml \n\n";
		print $OUTPUT '<?xml version="1.0" encoding="UTF-8" ?>'."\n";
		$xml_fire->no_CSA();
		$xml_fire->no_CCT();
		my $out=$xml_fire->print;
		print $OUTPUT $out;
	}
	if ($csa eq "NO" and $pred eq "NO"){
		if($fire->{output} eq 'text'){
			print $OUTPUT "\nfirestar did not detect ligand binding sites or transfer CSA annotation for your query\n";
		}
	}
	elsif($csa eq "YES" and $pred eq "NO"){
		if($fire->{output} eq 'text'){
			print $OUTPUT "\nfirestar did not detect ligand binding sites for your query\n";
		}
	}
	elsif($csa eq "NO" and $pred eq "YES"){
		if($fire->{output} eq 'text'){
			print $OUTPUT "\nfirestar did not detect CSA annotation for your query\n";
		}
	}
	elsif($csa eq "NO" and $pred eq "NOT_PASS"){
		if($fire->{output} eq 'text'){
			print $OUTPUT "\nfirestar detects some binding information, but it doesn't satisfy your filters\n";
		}
	}
	close($OUTPUT);
	return 0;
}



sub APPRIS_printer{
	my ($this,$join,$fire,$lib)=@_;
	my %suma_de_todos;
	my %suma_de_csa;
	my $name;
	my $compounder;
	my $position;
	my $ord=0;
	my $OUTPUT;
	$OUTPUT = \*STDOUT  unless(open($OUTPUT,'>',$fire->{outfile}));
	print $OUTPUT "######\n";
	if ($fire->{csa_option} ne "ONLY"){
		foreach my $incr(@{$this->{order}}){
			my @spl=split(/ /,$incr);
			foreach my $key(@spl){
			#       filtros para APPRIS: 
			#               aqui eliminamos todos los sites de compuestos tageados como NON_COGNATE or POSSIBLE COGNATE
				if ($join->{summary}{nice_try}[$key] eq "NON_COGNATE" or $join->{summary}{nice_try}[$key] eq "POSSIBLE_COGNATE"){next;}
			#               aqui eliminamos todos los sites de union a metales que salgan de alineamientos con un coverage inferior a 65%
				if ($join->{site_score}{$key}[2] < 0.65 && ($join->{summary}{type}[$key] eq "MET" or $join->{summary}{type}[$key] eq "MET_POSS")){next;}
			#               aqui eliminamos todos los sites con score de fiabilidad por debajo del umbral establecido en $cutoff
				if ($join->{summary}{reliability}[$key]<$fire->{cutoff}){next;}
				$ord++;
				my $comptmp=$join->{summary}{compid}[$key];
				my @comps=split(/\s+/,$comptmp);
				foreach my $i (@comps){
					my $candidate=undef;
					if ($i=~/(\w+)\(.+/){$candidate=$1;}
					if (defined $candidate and exists $lib->{cognate}{$candidate}){
						$compounder=$candidate;
						last;
					}
				}
				if (exists $suma_de_todos{$compounder}){
					my $new="other".$ord."-".$compounder;
				#	push(@{$suma_de_todos{$new}},$join->{summary}{resfreq}[$key],$join->{summary}{SQUARE}[$key],$join->{summary}{reliability}[$key]);
					push(@{$suma_de_todos{$new}},$join->{summary}{resfreq}[$key],$join->{site_score}{$key}[1],$join->{summary}{reliability}[$key]);
				}
				else {
				#	push(@{$suma_de_todos{$compounder}},$join->{summary}{resfreq}[$key],$join->{summary}{SQUARE}[$key],$join->{summary}{reliability}[$key]);
					push(@{$suma_de_todos{$compounder}},$join->{summary}{resfreq}[$key],$join->{site_score}{$key}[1],$join->{summary}{reliability}[$key]);
				}
			}
		}
	}
	############# Filtro CSA
	if ($fire->{csa_option} eq "YES" or $fire->{csa_option} eq "ONLY"){
		for(my$i=1;$i<($join->{csaincr}+1);$i++){
			chomp($join->{summary_csa}{resnum}[$i]);
			my @number_composition=split(/ /,$join->{summary_csa}{resnum}[$i]);
			foreach (@number_composition){
				$suma_de_csa{$_}=$join->{all_csa}{$_};
			}
		}
	}
	############# END 
	my @sequence=split(//,$fire->{sequence});
	my %aa_totales;
	if (scalar(keys%suma_de_todos)>0 or scalar(keys%suma_de_csa)>0){
		foreach my $i (keys%suma_de_todos){
			my $flag="NO";
			my $compuesto;
			if ($i=~/other\d+\-(\w+)/){$compuesto=$1;}
			else {$flag="YES";}
			my $site=${$suma_de_todos{$i}}[0];
			my $relia_score=${$suma_de_todos{$i}}[2];
			$site=~s/\s//g;
			chomp($site);
			my @parking=split(/;/,$site);
		#	my @scorez=split(/ /,${$suma_de_todos{$i}}[1]);
		#	my %all_scorez;
		#	foreach (@scorez){
		#		if ($_=~/(\d+)\((\d)\)/){$all_scorez{$1}=$2;}
		#	}
			foreach (@parking){
				$_=~s/\s//g;
				my @info=split(/=/,$_);
				if ($info[0]=~/([A-Z])\((\d+)\)/){$name=$1;$position=$2;}
				if (exists $aa_totales{$position}){
					my $scorium;
				#	if ($flag eq "YES"){$scorium="$i\[$info[1],$all_scorez{$position},$relia_score\]";}
					if ($flag eq "YES"){$scorium="$i\[$info[1],".sprintf("%.1f",${$suma_de_todos{$i}}[1]).",$relia_score\]";}
				#	else {"$compuesto\[$info[1],$all_scorez{$position},$relia_score\]";}
					else {"$compuesto\[$info[1],".sprintf("%.1f",${$suma_de_todos{$i}}[1]).",$relia_score\]";}
					$aa_totales{$position}[1].= "|".$scorium;
				}
				else {
					my $motif=$sequence[($position-7)].$sequence[($position-6)].$sequence[($position-5)].$sequence[($position-4)];
					$motif=$motif.$sequence[($position-3)].$sequence[($position-2)].$sequence[($position-1)].$sequence[$position];
					$motif=$motif.$sequence[($position+1)].$sequence[($position+2)].$sequence[($position+3)];
					$motif=$motif.$sequence[($position+4)].$sequence[($position+5)];
					push(@{$aa_totales{$position}},$motif);
					my $scorium;
				#	if ($flag eq "YES"){$scorium="$i\[$info[1],$all_scorez{$position},$relia_score\]";}
					if ($flag eq "YES"){$scorium="$i\[$info[1],".sprintf("%.1f",${$suma_de_todos{$i}}[1]).",$relia_score\]";}
				#	else {$scorium="$compuesto\[$info[1],$all_scorez{$position},$relia_score\]";}
					else {$scorium="$compuesto\[$info[1],".sprintf("%.1f",${$suma_de_todos{$i}}[1]).",$relia_score\]";}
					push(@{$aa_totales{$position}},$scorium);
				}
			}
		}
		foreach my $posi (keys%suma_de_csa){
			if (exists $aa_totales{$posi}){
				my $scorium="Cat_Site_Atl[1.00,$suma_de_csa{$posi},XXX]";
				$aa_totales{$posi}[1].= "|".$scorium;
			}
			else {
				my $motif=$sequence[($posi-7)].$sequence[($posi-6)].$sequence[($posi-5)].$sequence[($posi-4)];
				$motif=$motif.$sequence[($posi-3)].$sequence[($posi-2)].$sequence[($posi-1)].$sequence[$posi];
				$motif=$motif.$sequence[($posi+1)].$sequence[($posi+2)].$sequence[($posi+3)];
				$motif=$motif.$sequence[($posi+4)].$sequence[($posi+5)];
				push(@{$aa_totales{$posi}},$motif);
				my $scorium="Cat_Site_Atl[1.00,$suma_de_csa{$posi},XXX]";
				push(@{$aa_totales{$posi}},$scorium);
			}
		}
		my @results=sort{$a<=>$b}keys%aa_totales;
		my $resumen=undef;
		foreach (@results){
			if (defined $resumen){$resumen=$resumen.",".$_;}
			else {$resumen=$_;}
			print $OUTPUT "$_\t$aa_totales{$_}[0]\t$aa_totales{$_}[1]\n";
		}
		#print $OUTPUT ">>>$fire->{queryname}\t",scalar(keys%aa_totales),"\t$resumen\n";
		print $OUTPUT ">>>Query\t",scalar(keys%aa_totales),"\t$resumen\n";
	}
	else{
		#print $OUTPUT ">>>$fire->{queryname}\t0 predicted residues\n";
		print $OUTPUT ">>>Query\t0 predicted residues\n";
	}
	
	close($OUTPUT);
	return 0;
}



sub text_printer{
	my ($this,$join,$fire,$lib,$pulser)=@_;
	my $xml_fire;
	if ($fire->{output} eq "XML"){	$xml_fire=XML_firestar->new(-root=>'Prediction');}
	my %resume_compounds;
	my @clustids;
	my @ecgoes;
	if($join->{csaincr}>0){@clustids=keys%{$join->{csa_clust}};}
	my $OUTPUT;
	my $CSA;
	$OUTPUT = \*STDOUT  unless(open($OUTPUT,'>',$fire->{outfile}));
	if (($fire->{csa_option} eq "YES" or $fire->{csa_option} eq "ONLY") and $fire->{output} ne 'SIAM'){
		if ($join->{csaincr} > 0){
			unless ($fire->{output} eq "XML"){	print $OUTPUT "\n";}
			else{	$CSA=$xml_fire->crea(-tag=>'Catalytic_sites_list');}
		}
		elsif($fire->{csa_option} eq "ONLY"){$this->no_results($fire,"NO","YES");}
		elsif ($fire->{output} eq "XML"){	$xml_fire->no_CSA();}
		for(my$i=1;$i<($join->{csaincr}+1);$i++){
		##########################################################
			my $ecs=undef;
			my $key_res=$join->{summary_csa}{resnum}[$i];
			my @control=split(/ /,$key_res);
			$key_res=join(' ',sort@control);
			my @csiteids=split(/ /,$join->{csa_sites}{$key_res});
			my %ECs_output;
			foreach my $csiteid(@csiteids){
				my @source=$this->{dbi}->CSA_LIT(-csite=>$csiteid);
				if (exists $lib->{csa2ec}{$source[0]}){
					my $EC=$lib->{csa2ec}{$source[0]};
					unless (exists $ECs_output{$EC}){
						$ECs_output{$EC}=5;
                                		if (exists $lib->{ec2go}{$EC}){
							my @tmp=@{$lib->{ec2go_full}{$EC}};
							unless (defined $ecs){	$ecs.=join("\t\t\t\t\t",@tmp);}
							else {	$ecs=$ecs."\t\t\t\t\t".join("\t\t\t\t\t",@tmp);}
						}
					}
				}
			}
			my $ec_numbers=join(' ',keys%ECs_output);
			#########################################################
			unless (defined $ecs){$ecs="CAT\t$i\tno EC nor GO information found in our database;";}
			else{chomp($ecs);}
			my %unique;
			foreach my $zz (@clustids){
				$unique{$join->{csa_clust}{$zz}}=5;
			}
			my $tmps=join(" ",keys(%unique));
			if ($fire->{output} ne "XML"){
				print $OUTPUT "CAT\t$i\tEC_number:\t\t$ecs\nCAT\t$i\tEvidence:\t\tLiterature\nCAT\t$i\tHomologs\t\t$tmps\n";
				print $OUTPUT "CAT\t$i\tResidue_positions:\t$join->{summary_csa}{resnum}[$i]\nCAT\t$i\tResidue_composition:\t$join->{summary_csa}{resname}[$i]\n";
				print $OUTPUT "CAT\t$i\tSite_score:\t\t$join->{summary_csa}{score}[$i]\n\n";
			}
			elsif ($fire->{output} eq "XML"){	$xml_fire->CSA_prediction($CSA,$ec_numbers,$tmps,$i,$lib,$join,keys%unique);}
		}
		if ($fire->{output} eq "XML" && $join->{csaincr} > 0){$xml_fire->app(-child=>$CSA);}
	}
	elsif($fire->{output} eq 'SIAM'){
		my %has;
		my %modify;
		my %hasifnot;
		for(my$i=1;$i<($join->{csaincr}+1);$i++){
			my %unique;
			my @total_list;
			my @strange_list;
			foreach my $zz (@clustids){
				$unique{$join->{csa_clust}{$zz}}=5;
			}
			foreach my $zz (keys%unique){
				if (exists $lib->{csa2ec}{$zz}){
					my $EC=$lib->{csa2ec}{$zz};
					my @tmp;
					if (exists $lib->{ec2go}{$EC}){@tmp=@{$lib->{ec2go}{$EC}};}
					foreach my $i (@tmp){
						if (exists $lib->{good_CSA}{$zz}){push(@total_list,$i);}
						else{push(@strange_list,$i);}
					}
				}
			}
			unless (scalar@total_list == 0 and scalar@strange_list == 0){
				my $controller=0;
				my @residuos=split(//,$join->{summary_csa}{score}[$i]);
				foreach my $kkk(@residuos){
					if ($kkk == 6){$controller++;}
				}
				if ($controller == scalar@residuos){
					foreach my $kkk(@total_list){	$has{$kkk}=5;}
					foreach my $kkk(@strange_list){  $has{$kkk}=5;}
				}
				else{
					foreach my $kkk(@total_list){	$modify{$kkk}=5;}
					foreach my $kkk(@strange_list){  $hasifnot{$kkk}=5;}
				}
			}
		}
		my @output;
		$this->GO_printer(\@output,\%has,"HAS");
		$this->GO_printer(\@output,\%modify,"MODIFY");
		$this->GO_printer(\@output,\%hasifnot,"HAS IF NOT SIAM");
		$this->GO_printer(\@output,\%{$fire->{unknown}},"UNKNOWN");
		$this->GO_printer(\@output,\%{$fire->{fail_list}},"NOT");
		foreach my $line (@output){	print $OUTPUT $line;}
		close($OUTPUT);
		return 0;
	}
	if ($fire->{csa_option} ne "ONLY"){
		my $CCT=undef;
		my $ord=0;
		unless ($fire->{output} eq "XML"){print $OUTPUT "\n";}
		foreach my $incr(@{$this->{order}}){
			my @spl=split(/ /,$incr);
			foreach my $key(@spl){
##### aqui estamos excluyendo unos resultados dependiendo de los parametros que el usuario nos ha pasado (CUT-OFF reliability, CSA yes/no/only ...)
				if (exists $join->{shinigami}{$key}){next;}
				if ($fire->{cognate_option} eq "NO" &&  $join->{summary}{nice_try}[$key] eq "NON_COGNATE"){next;}
				if ($join->{summary}{reliability}[$key]<$fire->{cutoff}){next;}
##### end
				$ord++;
				if ($fire->{output} eq "XML"){	unless(defined $CCT){$CCT=$xml_fire->crea(-tag=>'binding_sites_list');}}
				my $comptmp=$join->{summary}{compid}[$key];
				my @comps=split(/\s+/,$comptmp);
				my @compgo;
				my @goname;
				for(@comps){
					my $id=$_;
					$id=~s/\(\d+\)//;
					if(exists$lib->{cif2go}{$id}{go2}){
						push(@compgo,$lib->{cif2go}{$id}{go2});
						push(@goname,$lib->{cif2go}{$id}{name});
					}
				}
				my $compgos=join(" ",@compgo);
				my $goname=join(", ",@goname);
				if ($pulser eq "STRICT"){
					$join->{summary}{compname}[$key]="UNKNOWN; Please check the alignments, this prediction is based on the most conserved residues in the second firestar analysis round";
					if ($join->{summary}{nice_try}[$key] eq "COGNATE"){$join->{summary}{nice_try}[$key]="POSSIBLE_COGNATE";}
					$join->{summary}{compid}[$key]="based on ".$join->{summary}{compid}[$key];
				}
				if ($fire->{output} ne "XML"){
					print $OUTPUT "SITE\t$ord\tCompound_name:\t\t\t$join->{summary}{compname}[$key]\n";
					print $OUTPUT "SITE\t$ord\tSite classification:\t\t$join->{summary}{nice_try}[$key]\n";
					print $OUTPUT "SITE\t$ord\tCompound_id:\t\t\t$join->{summary}{compid}[$key]\n";
					print $OUTPUT "SITE\t$ord\tCompound_GO:\t\t\t$compgos\n";
					print $OUTPUT "SITE\t$ord\tGO_name:\t\t\t$goname\n";
					print $OUTPUT "SITE\t$ord\tResidue_positions:\t\t$join->{summary}{resnum}[$key]\n";
					print $OUTPUT "SITE\t$ord\tResidue_composition:\t\t$join->{summary}{resname}[$key]\n";
					print $OUTPUT "SITE\t$ord\tPer-residue Probab. score:\t$join->{summary}{resfreq}[$key]\n";
					print $OUTPUT "SITE\t$ord\tReliability:\t\t\t$join->{summary}{score}[$key]\n";
					print $OUTPUT "SITE\t$ord\tNumber_of_homologs:\t\t$join->{summary}{numOfSites}[$key]\n\n";
				}
				elsif ($fire->{output} eq "XML"){	$xml_fire->CCT_prediction($CCT,$compgos,$goname,$key,$join,$fire);}
			}
		}
		if ($ord == 0 and ($join->{csaincr}==0 or $fire->{option} eq "NO")){$this->no_results($fire,"YES","NOT_PASS");}
		elsif ($ord == 0 and $fire->{output} eq "XML"){	$xml_fire->no_CCT();}
		elsif ($ord > 0 and $fire->{output} eq "XML"){$xml_fire->app(-child=>$CCT);}
	}
	if ($fire->{output} eq "XML"){
		print $OUTPUT "Content-type: text/xml \n\n";
		print $OUTPUT '<?xml version="1.0" encoding="UTF-8" ?>'."\n";
		my $out=$xml_fire->print;
		print $OUTPUT $out;
	}
	close($OUTPUT);
	return 0;
}


sub web_printer{
	my ($this,$fire,$join,$lib)=@_;
	my $run_type;
	$this->web_headers($fire,$join);
	my $HHS; my $PSI;
	$HHS=\*STDOUT unless(open($HHS,'>',"$fire->{outfile}.hhsall.php"));
	$PSI=\*STDOUT unless(open($PSI,'>',"$fire->{outfile}.psiall.php"));
	print $HHS "$this->{extended_results_header}\n$this->{score_tables}\n";
	print $PSI "$this->{extended_results_header}\n$this->{score_tables}\n";
	my %result=$this->extended_generator($fire);
	my @sorted_keys=sort{$a<=>$b}keys%result;
	if(scalar@sorted_keys>0){
		my $last=scalar@sorted_keys;
		print $PSI "$this->{exResults_table_header}<input type=hidden name=\"query\" value=\"$fire->{queryname}\"></form>\n$this->{download_psi_output}";
		print $HHS "$this->{exResults_table_header}<input type=hidden name=\"query\" value=\"$fire->{queryname}\"></form>\n$this->{download_hhs_output}";
		my $contador_psi=0;
		my $contador_hhs=0;
		for(@sorted_keys){                                                      #####################################
			my $key=$_;                                                     #####################################
			if ($fire->{parameters}{$key}{program} eq "PSI"){print $PSI "<tr><td width=\"99%\">$result{$key}</td></tr>\n";$contador_psi++;}       ##########
			elsif($fire->{parameters}{$key}{program} eq "HHS"){print $HHS "<tr><td width=\"99%\">$result{$key}</td></tr>\n";$contador_hhs++;}       ######
		}                                                               #################
		if ($contador_psi==0){
			print $PSI $this->{no_results_psi_output};
		}
		elsif ($contador_hhs==0){
			print $HHS $this->{no_results_hhs_output};
		}
		else{
			print $PSI "$this->{exResults_table_header}<input type=hidden name=last value=$last>\n$this->{download_psi_output}";
			print $HHS "$this->{exResults_table_header}<input type=hidden name=last value=$last>\n$this->{download_hhs_output}";
		}
	}
	else {
		print $PSI "$this->{exResults_table_header}\n<br>\n$this->{no_results_psi_output}";
		print $HHS "$this->{exResults_table_header}\n<br>\n$this->{no_results_hhs_output}";
	}
	close($HHS);
	close($PSI);
	$this->summary_generator($fire,$join,$lib);
	if (defined $fire->{mail}){	$this->mail_sender($fire);}
	return 0;
}

sub extended_generator{
	my ($this,$fire)=@_;
	my %result;
	my %clustids;
	my %msa_list;
	foreach my $key (keys(%{$fire->{psiout}})){
		if (exists$fire->{csiteid_order}{$key}){
			$fire->{parameters}{$key}{afm_score}=~s/@/C/g;
			my @sfmscore=split(//,$fire->{parameters}{$key}{afm_score});
			my $afm_score;
			my @aln_color=();
			my @position_score=();
			for(@sfmscore){
				$afm_score.=residue_colorer($_,\@aln_color,\@position_score);
			}
			my %siteline=();
			my @tempres=split(//,$fire->{parameters}{$key}{afm_tempt});
			my $position=$fire->{parameters}{$key}{t_start};
			my $tmplate_out="";
			my %sitescore;
			my %sitescore_norm;
			my $global_score=0;
			for(my $i=0;$i<$#tempres+1;$i++){
				my $alnlet=$tempres[$i];
				foreach my$id(keys%{$fire->{num2res}{$key}}){
					if($alnlet =~ /[A-Z]/){
						if(exists$fire->{num2res}{$key}{$id}{$position}){
							my $pos_color;
							if($position_score[$i] > 0.80){ $pos_color=" color=\"white\"";}
							$siteline{$id}.="<td bgcolor=\"$aln_color[$i]\"><font title=\"$position\"$pos_color>$fire->{num2res}{$key}{$id}{$position}</font></td>";
							$sitescore{$id}=$sitescore{$id}+$position_score[$i];
							$sitescore_norm{$id}++;
						}
						else{$siteline{$id}.="<td>-</td>";}
					}
					elsif($alnlet eq "-"){	$siteline{$id}.="<td>-</td>";}
				}
				if($alnlet =~/[A-Z]/){
					$tmplate_out.="<td><font title=\"$position\">$alnlet</font></td>";
					$position++;
					$global_score=$global_score+$position_score[$i];
				}
				else{   
					$tmplate_out.="<td>-</td>";
					$global_score=$global_score+$position_score[$i];
				}
			}
			foreach my $id(keys%{$fire->{num2res}{$key}}){
				if($sitescore_norm{$id}!=0){
					$sitescore{$id}=$sitescore{$id}/$sitescore_norm{$id};
					$sitescore{$id}=sprintf("%.2f",$sitescore{$id});
				}
			}
			my $test=join("",values%siteline);
			if($test =~ /[A-Z]/ and $#tempres >0){
				$global_score=$global_score/($#tempres+1);
				$global_score=sprintf("%.2f",$global_score);
				$clustids{$fire->{parameters}{$key}{template}}=5;
				my @targres=split(//,$fire->{parameters}{$key}{afm_query});
				my $aln_pos=$fire->{parameters}{$key}{q_start};
				my $numbering;
				my $signal;
				my $num_length=0;
				my $query_seq_out;
				for(@targres){
					if($_ ne "-"){
						$query_seq_out.="<td>$_</td>";
						if($aln_pos=~/0$/){
							$num_length=length$aln_pos;
							$numbering.="<td colspan=$num_length>$aln_pos</td>";
							$signal.="<td>|</td>";
						}
						else{
							if($num_length > 1){	$num_length=$num_length-1;}
							else{	$numbering.="<td>.</td>";}
							$signal.="<td> </td>";
						}
						$aln_pos++;
					}
					else{
						if($num_length > 1){$num_length=$num_length-1}
						else{$numbering.="<td> </td>";}
						$signal.="<td> </td>";
						$query_seq_out.="<td>-</td>";
					}
				}
				chomp $fire->{parameters}{$key}{identities};
				my $query_print;
				$msa_list{$key}=$fire->{parameters}{$key}{template};
				$result{$key}.="<table width=\"100%\"border=0.5 margin=\"0px\"><tr>\n";
				$result{$key}.="</td><td bgcolor=pink title=\"Multiple Sequence Alignment\"align=center>ALN&nbsp;$key<br>";
				$result{$key}.="<input type=checkbox name=\"$key\" value=\"$msa_list{$key}\"></td><td bgcolor=\"#AAAAAA\">\n";
				if($fire->{runtype} eq 'str' or $fire->{runtype} eq 'pdb'){
					$query_print=substr($fire->{queryname},0,4);;
					my $straln_file="$fire->{tmpfile}\_$key";
					$result{$key}.="<td bgcolor=\"yellow\" title=\"Structural Alignment Options\"><pre><form method=POST action=\"../straln/$straln_file\_redir.php\" ";
					$result{$key}.="target=_blank name=$fire->{parameters}{$key}{template} style=\"margin:0px;\">$query_print:   from <input type=text size=4 maxlength=4 ";
					$result{$key}.="name=q_start value=$fire->{parameters}{$key}{q_start}> to <input type=text size=4 maxlength=5 name=q_end value=$fire->{parameters}{$key}{q_end}>";
					$result{$key}.="<input type=reset value=\"reset\" name=$fire->{parameters}{$key}{template}><br>$fire->{parameters}{$key}{template}:\tfrom <input type=text ";
					$result{$key}.="size=4 maxlength=4 name=t_start value=$fire->{parameters}{$key}{t_start}> to <input type=text size=4 maxlength=5 name=t_end ";
					$result{$key}.="value=$fire->{parameters}{$key}{t_end}><input type=submit value=\"run LGA\" name=$fire->{parameters}{$key}{template}></form></pre>";
					$this->structure_alignment($straln_file,"fast",$fire->{parameters}{$key}{template},$fire);
				}
				my @info=$this->template_description($fire->{parameters}{$key}{template});
				my $info;
				if ($info[2]=~/\w+/ && $info[1]=~/\w+/){
					$info="<b>$fire->{parameters}{$key}{template}:</b> $info[0]<br><a href=\"http://www.expasy.ch/uniprot/$info[1]\">$info[1]</a>: $info[2]";
				}
				else {	$info="<b>$fire->{parameters}{$key}{template}:</b> $info[0]<br>";}
				$result{$key}.="<td bgcolor=\"#AAAAAA\"></td><td bgcolor=white><table width=\"100%\"><tr><td rowspan=2>Evalue: <font color=blue size=2>$fire->{parameters}{$key}{expect}";
				$result{$key}.="</font></td><td>$fire->{parameters}{$key}{identities}</td></tr>";
				$result{$key}.="<tr><td> Download pairwise alignment: <a href=\"pairwise/$fire->{tmpfile}\_$key.faa\" target=_blank>FASTA</a></td></tr></table></td>\n";
				open(FAA,'>',"$fire->{faatmp}/pairwise/$fire->{tmpfile}\_$key.faa");
				print FAA ">$fire->{queryname}\n$fire->{parameters}{$key}{query_aligned}\n>$fire->{parameters}{$key}{template}\n$fire->{parameters}{$key}{templ_aligned}\n";
				close FAA;
				$result{$key}.="<td bgcolor=\"#AAAAAA\"></td><td bgcolor=\"white\">$info</td></tr></table>\n";
				if($fire->{runtype} eq 'pdb' and $fire->{cluster} eq $fire->{parameters}{$key}{template}){
					$result{$key}.="<div id=\"bigscroll\" style=\"background-color:beige;\"><pre>";
					$result{$key}.="<font color=red>Your query:</font> $fire->{queryname} <font color=red> is inside the cluster:</font> t_$fire->{parameters}{$key}{template}";
					$result{$key}.="<font color=red> in FireDB database</font>";
					$query_print=$fire->{queryname};
				}
				else{	$result{$key}.="<div id=bigscroll><pre>";}
				$result{$key}.="<table border=\"0\" cellspacing=\"0\" celpadding=\"0\">";
				$result{$key}.="<tr><td colspan=2>Query:</td><td>$query_print</td><td>-</td><font color=\"red\">$numbering</font></tr>\n";
				$result{$key}.="<tr><td colspan=3></td><td bgcolor=\"white\">Score</td><font color=\"red\">$signal</font></tr>\n";
				$result{$key}.="<tr><td colspan=2>Query:</td><td>$query_print</td><td>-</td>$query_seq_out<tr>\n";
				$result{$key}.="<tr><td colspan=2>Consensus:</td><td><a href=\"../../Php/FireDB.php?pdbcode=$fire->{parameters}{$key}{template}&cutoff=35\" target=\"_blank\">";
				$result{$key}.="t_$fire->{parameters}{$key}{template}</a></td><td>-</td>$tmplate_out<tr>\n";
				$result{$key}.="<tr><td colspan=2><a href=\"http://firedb.bioinfo.cnio.es/Php/square.php\" target=_blank>SQUARE</a>:</td><td>";
				$result{$key}.="</td><td bgcolor=\"white\"><font color=\"darkblue\">$global_score</font></td><font color=\"darkblue\">$afm_score</font><tr>\n";
				for (my $i=0;$i < scalar(@{$fire->{csiteid_order}{$key}});$i++){
					my $id=$fire->{csiteid_order}{$key}[$i];
					my $value=$fire->{csiteid_values}{$key}[$i];
					my $val_color;
					if($fire->{extended}{$id}{numseqs} < 5){$val_color="black";}
					else{   
						if($value < 20){$val_color="red";}
						if($value > 19){$val_color="yellow";}
						if($value > 49){$val_color="green";}
					}
					if($siteline{$id} =~/[A-Z]/){
						if($fire->{extended}{$id}{sitetype} eq "CCT"){
							my $comp=$this->COMPIDS($fire->{extended}{$id}{compids});
							if($fire->{csite_cases}{$id}>=1 and $value > 2){
								$result{$key}.="<tr><td><a href=\"../../Php/Validate.php?csiteid=$id&cutoff=35\"><font title=\"Evolutively related sites\">";
								$result{$key}.="E=$fire->{csite_cases}{$id}</font></a></td><td bgcolor=\"#AAAAAA\" align=center><font color=\"$val_color\">$value\%</font>";
								$result{$key}.="</td><td>$comp</td><td bgcolor=\"white\"><font color=\"darkblue\">$sitescore{$id}</font></td>$siteline{$id}<tr>\n";
							}
							elsif($value > 2){
								$result{$key}.="<tr><td>&nbsp&nbsp&nbsp</td><td bgcolor=\"#AAAAAA\" align=center><font color=\"$val_color\">$value\%</font></td>";
								$result{$key}.="<td>$comp</td><td bgcolor=\"white\"><font color=\"darkblue\">$sitescore{$id}</font></td>$siteline{$id}<tr>\n";
							}
						}
						elsif($fire->{extended}{$id}{sitetype} eq "CSA"){
							my @pdbid=$this->{dbi}->CSA_LIT(-csite=>$id);
							my $evi;
							if($fire->{extended}{$id}{evidency} eq "lit"){$evi="<font color=\"green\">Literature</font>";}
							elsif($fire->{extended}{$id}{evidency} eq "psi"){$evi="<font color=\"yellow\">PSI-BLAST</font>";}
							$result{$key}.="<tr><td colspan=2 bgcolor=\"#AAAAAA\">$evi</td><td title=\"Catalytic Site Atlas: $pdbid[0]\">";
							$result{$key}.="<a href=\"http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=$pdbid[0]\" target=_blank>";
							$result{$key}.="$fire->{extended}{$id}{sitetype}</a></td><td bgcolor=\"white\"><font color=\"darkblue\">$sitescore{$id}</font></td>$siteline{$id}<tr>\n";
						}
					}
				}
				$result{$key}.="</table></pre></div>\n";
			}
		}
	}
	return %result;
}

sub residue_colorer{
	my ($residue,$color,$pos)=@_;
	my $colorafmscore;
	if($residue eq "1"){
		$colorafmscore="<td  bgcolor=\"#E0FFFF\">$residue</td>";
		push(@{$color},"#E0FFFF");
 		push(@{$pos},"0.45");
 	}
	elsif($residue eq "2"){
		$colorafmscore="<td bgcolor=\"#B0E2FF\">$residue</td>";
		push(@{$color},"#B0E2FF");
		push(@{$pos},"0.60");
	}
	elsif($residue eq "3"){
		$colorafmscore="<td bgcolor=\"#55DDFF\">$residue</td>";
		push(@{$color},"#55DDFF");
		push(@{$pos},"0.75");
	}
	elsif($residue eq "4"){
		$colorafmscore="<td bgcolor=\"#1E90FF\"><font color=white>$residue</font></td>";
		push(@{$color},"#1E90FF");
		push(@{$pos},"0.85");
	}
	elsif($residue eq "5"){
		$colorafmscore="<td bgcolor=\"#1874CD\"><font color=white>$residue</font></td>";
		push(@{$color},"#1874CD");
		push(@{$pos},"0.90");
	}
	elsif($residue eq "C"){
		$colorafmscore="<td bgcolor=\"darkblue\"><font color=white>$residue</font></td>";
		push(@{$color},"darkblue");
		push(@{$pos},"1");
	}
	elsif($residue eq "-"){
		$colorafmscore="<td>$residue</td>";
		push(@{$color},"");
		push(@{$pos},"0");
	}
	else{print STDERR "Wrong scores in alignment. AFM may not be working!!!";}
	return $colorafmscore;
}

sub structure_alignment{
	my ($this,$name,$speed,$template,$fire)=@_;
	open(OUTB,'>',"../straln/$name\_redir.php");
	print OUTB $this->{struct_alig_header1};
	print OUTB "<meta http-equiv=\"refresh\" content=\"2;url=$name.php\" >";
	print OUTB $this->{struct_alig_header2};
	print OUTB "\n\$tmpfile=\"$name\";\n\$querytype=\"$speed\";\n\$filetype=\"$fire->{runtype}\";\n\$template=\"$template\";\n";
	print OUTB "\$target=\"$fire->{infile}\";\n\$target_chain=\"$fire->{chain}\";\n";
	print OUTB $this->{struct_alig_php_form};
	close OUTB;
	open(WAIT,'>',"../straln/$name.php");
	print WAIT $this->{struct_alig_wait_part1};
	print WAIT "<a href=\"$name.php\">this page: $name</a>";
	print WAIT $this->{struct_alig_wait_part2};
	close WAIT;
}

sub template_description{
	my ($this,$clustid)=@_;
	my @info=$this->{dbi}->chain_info(-cadid=>$clustid);
	if (defined $info[1]){return($info[0],$info[1],$info[2]);}
	else{
		my $codes=$this->{dbi}->clust_codes(-cadid=>$clustid);
		my @clustcodes=split(/\s+/,$codes);
		foreach my $cadid(@clustcodes){
			my @info2=$this->{dbi}->chain_info(-cadid=>$cadid);
			if (defined $info2[1]){return($info[0],$info2[1],$info2[2]);}
		}
	}
	return($info[0],"","");
}

sub COMPIDS{
	my ($this,$compiz)=@_;
	my @spl=split(/ /,$compiz);
	my $compids;
	my $first=shift(@spl);
	$first=~s/\(\d+\)//;
	my $name =$this->{dbi}->compound_name(-compid=>$first);
	$compids="<a href=\"../../Php/ligand/index.html?id=$first\" target=_blank><font title=\"$name\">$first</font></a>";		####### CALL 2 LIGAND.php
	if(scalar@spl>0){
		my $counter=0;
		for(my $i=0;$i<scalar(@spl);$i++){
			$counter++;
			$spl[$i]=~s/\(\d+\)//;
			$name =$this->{dbi}->compound_name(-compid=>$spl[$i]);
			$compids.=" <a href=\"../../Php/ligand/index.html?id=$spl[$i]\" target=_blank><font title=\"$name\">$spl[$i]</font></a>";
			if ($counter == 4){last;}
		}
	}
	else {$compids.="\t";}
	$compids.="\t";
	return $compids;
}

sub summary_generator{
	my ($this,$fire,$join,$lib)=@_;
	my $OUTPUT;
	$OUTPUT = \*STDOUT  unless(open($OUTPUT,'>',"$fire->{outfile}.tmp.php"));
	my $variables=Config::IniFiles->new(-file => $fire->{config_file});
        my $blast_bin_path=$variables->val('PROGRAMS','bstbin');
        my $tmpchads=$variables->val('PATHS','square_dir')."/tmpchads";
	open(PN,'>',"$fire->{full}.pn");
	open(SN,'>',"$fire->{full}.sn");
	print PN "$fire->{tmpfile}.chk";
	print SN "$fire->{tmpfile}.faa";
	close PN;close SN;
        `$blast_bin_path/makemat -P $fire->{full}`;
	copy("$fire->{full}.faa","$tmpchads/$fire->{tmpfile}.faa");
	move("$fire->{full}.mtx","$tmpchads/$fire->{tmpfile}.mtx");
	unlink("$fire->{full}.pn","$fire->{full}.sn","$fire->{full}.mn","$fire->{full}.aux");
	print $OUTPUT "$this->{summary_header}$this->{summary_table_header}";
	my $height=$this->iframe_generator($fire,$join);
	print $OUTPUT "<div align=\"center\"><fieldset style=\"width:75%;\"><legend>Prediction summary</legend><iframe width=\"100%\" height=\"$height\" frameborder=\"0\" \n";
	print $OUTPUT "id= \"iframe\" name=\"iframe\" src=\"sites_frame/$fire->{tmpfile}\_ALL.html\">\n$this->{summary_page_iframe_end}";
	my @order=sort{$b<=>$a}(keys%{$join->{summary_order}});
	$this->structural_model($fire,$join,\@order);
	if ($join->{scalar_sites} >0 || $join->{csaincr}>0){
		open(SUM,'>',"$fire->{faatmp}/$fire->{tmpfile}\_sum.html");
		print SUM $this->{supplementary_page_header};
	}
	if($join->{scalar_sites} >0){	print $OUTPUT $this->{beginning_firestar_prediction};}
	if($join->{csaincr}>0){
		if ($join->{csaincr}==1){	print SUM "<h2 class=\"summary\">CATALYTIC SITE PREDICTION</h2>\n";}
		else {	print SUM "<h2 class=\"summary\">CATALYTIC SITE PREDICTIONS</h2>\n";}
		print $OUTPUT "<h3 class=\"summary\">Catalytic site (predicted based on CSA entries):</h3>";
	#	print $OUTPUT "<br><table width=\"600px\" border=\"0\" cellspacing=\"0\">\n\t\t<tr><td colspan=8><b>Catalytic site</b> (predicted based on CSA entries):</td></tr><tr>";
		print SUM "<h4 class=\"summary\">Here cross-references between EC codes and GO terms and TEMPLATES are listed for every CATALYTIC SITE prediction</h4>\n";
		print SUM "<br><table class=\"catalytic\">\n";
	}
	print $OUTPUT "<table width=\"98%\" border=\"0\" cellspacing=\"0\">\n";
	print $OUTPUT "<tr><td colspan=2><br></td></tr>";
	###############################################
	for(my$i=1;$i<$join->{csaincr}+1;$i++){
		my $key_res=$join->{summary_csa}{resnum}[$i];
		my @control=split(/ /,$key_res);
		$key_res=join(' ',sort@control);
		my @csiteids=split(/ /,$join->{csa_sites}{$key_res});
		my %ECs_output;
		my $ecs;
		my $out_controller="NO";
		my %template_list;
		foreach my $csiteid(@csiteids){
			my @source=$this->{dbi}->CSA_LIT(-csite=>$csiteid);
			$template_list{$source[0]}=5;
			if (exists $lib->{csa2ec}{$source[0]}){
				my $EC=$lib->{csa2ec}{$source[0]};
				if (exists $ECs_output{$EC}){next;}
				$ECs_output{$EC}=5;
				if (exists $lib->{ec2go}{$EC}){
					my $html_line="</td></tr><tr><td></td><td style=\"width:4px\"></td><td nowrap colspan=\"6\">";
					my @tmp=@{$lib->{ec2go_full}{$EC}};
					if ($out_controller eq "NO"){
						if (scalar@tmp==1){	$ecs=$tmp[0];}
						else {	$ecs=join($html_line,@tmp);}
						$out_controller="YES";
					}
					else{
						if (scalar@tmp==1){	$ecs.=$html_line.$tmp[0];}
						else{	$ecs.=$html_line.join($html_line,@tmp);}
					}
				}
			}	
		}
		if ($ecs eq ''){print SUM "<tr><td class=\"site\"><b>CAT\t$i</b></td><td>no EC nor GO information found in our database;</td></tr>";}
		else{print SUM "<tr><td class=\"site\"><b>CAT\t$i</b></td><td colspan=\"6\" nowrap>$ecs</td></tr>";}
		print SUM "<tr><td class=\"site\"><b>CAT\t$i</b></td><td style=\"text-align:right\"><b>TEMPLATES FROM</b></td></tr>";
		my $counter_tab=0;
		print SUM "<tr><td class=\"site\"></td><td style=\"text-align:right\"><b>LITERATURE</b></td>";
		foreach my $kk(keys%template_list){
			print SUM "<td class=\"template\"><a href=\"http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=$kk\" target=_blank>$kk</a></td>";
			$counter_tab++;
			my $res=$counter_tab/5;
			if ($res=~/\./ && $counter_tab!=scalar(keys%template_list)){
				print SUM "</tr><tr><td class=\"site\"></td><td></td>";
				$counter_tab=0;
			}
		}
		if ($counter_tab<5 and $counter_tab>0){
			for (my $res=$counter_tab+1; $res<6;$res++){	print SUM "<td class=\"template\"></td>";}
		}
		print SUM "</tr>\n</table>\n<br><br>\n";
		my @rub=split(//,$join->{summary_csa}{score}[$i]);
		my $somma;
		foreach my $z(@rub){$somma=$somma+$z;}
		$somma=sprintf("%.2f",$somma/scalar@rub);
		$somma.=" [max=6]";
		my @res_pak=split(";",$join->{summary_csa}{resresume}[$i]);
		my @spacer;
		foreach my $q (@res_pak){
			if ($q=~/[A-Z]\((\d+)\)/){push(@spacer,$1);}
		}
		my $res_composition = join(" ", @spacer);
		my $ref=$this->{template2model_csa}{$i};
		my $par1="seq_query=".$fire->{psiout}{$ref}[1]."&seq_templ=".$fire->{psiout}{$ref}[3]."&pdbK=".$fire->{psiout}{$ref}[2]."&coor=".$fire->{psiout}{$ref}[4];
		$par1.="&tmp=".$fire->{tmpfile}."&residues=".$res_composition."&psiout_rank=$ref";
		print $OUTPUT "\n<tr><td style=\"vertical-align:top; text-align:left;\" nowrap width=\"24%\">\n";
		print $OUTPUT "CAT     $i      Evidence:</td><td style=\"vertical-align:top; text-align:left;\">Literature</td></tr>\n";
		print $OUTPUT "<tr><td style=\"vertical-align:top; text-align:left;\" nowrap>\n";
		print $OUTPUT "CAT     $i      Residue_resume:</td><td style=\"vertical-align:top; text-align:left;\">$join->{summary_csa}{resresume}[$i]</td></tr>\n";
		print $OUTPUT "<tr><td style=\"vertical-align:top; text-align:left;\" nowrap>\n";
		print $OUTPUT "CAT     $i      AVG_SQUARE_residue_score:</td><td style=\"vertical-align:top; text-align:left;\">$somma</td></tr>\n";
		print $OUTPUT "<tr><td><br></td><td><br></td><td style=\"vertical-align:middle; text-align:center;\">\n";
		print $OUTPUT "<button onClick=\"Model_gen(this)\" value=\"$par1\" style=\"border: 0; background: transparent\">\n";
		print $OUTPUT "<img src=\"../../images/Model_button.png\">\n</button></td><td><a href=\"../../html/firestar_help.html#Model\" target=\"_parent\">\n";
		print $OUTPUT "<img align=\"center\" width=\"18px\" src=\"../../images/Question_mark.png\"></a></td></tr>\n";
	}
	if ($join->{csaincr}>0){
		print $OUTPUT "<tr><td colspan=2>You can find supplementary information about CATALYTIC SITE prediction <a href=\"$fire->{tmpfile}\_sum.html\">HERE</a>.</td></tr>\n</table><br><hr>";
	}
	if($join->{scalar_sites} >0){
		print $OUTPUT "<table width=\"98%\">\n<tr><td colspan=3><b>Ligand Binding sites</b> (predicted based on FireDB entries)<br>\n";
		print $OUTPUT "<b>Please remember that the numbering of sites and pockets is not related !</b><br><br></td></tr>";
		if ($join->{scalar_sites}==1){print SUM "\n<br><br><hr>\n<h2 class=\"summary\">BINDING SITE PREDICTION</h2>\n";}
		else {print SUM "\n<hr>\n<h2 class=\"summary\">BINDING SITE PREDICTIONS</h2>\n";}
	}
	my %pharma_all;
	my $ord=0;
	$lib->pdb2id($fire);
	my $connection_check=$this->{dbi}->ping();
	unless($connection_check == 1){$this->{dbi}=DBI_firestar->connect($fire->{config_file});}
	%{$this->{man_annot}}=$this->{dbi}->manually_annotated($lib);
        %{$this->{pharma}}=$this->{dbi}->pharma_annotated($lib);
	foreach my $incr(@{$this->{order}}){
		my @spl=split(/ /,$incr);
		foreach my $key(@spl){
			if (exists $join->{shinigami}{$key}){next;} ## aqui no imprimimos info acerca de aquellos site que solapan totalmente con otro
			$ord++;
			my $comptmp=$join->{summary}{compid}[$key];
			my @comps=split(/\s+/,$comptmp);
			my @united;
			my @pharma;
			for(@comps){
				my $id=$_;
				$id=~s/\(\d+\)//;
				if(exists$lib->{cif2go}{$id}{go2}){
					push(@united,$lib->{cif2go}{$id}{go2}."(".$lib->{cif2go}{$id}{name}.")");
				}
				if (exists $this->{man_annot}{$id} or exists $this->{pharma}{$id}){
					push(@pharma,$id);
					$pharma_all{$id}=5;
				}
			}
			my $united=join(", ",@united);
			if ($united eq""){$united="-"}
			my @res_pak=split(";",$join->{summary}{resresume}[$key]);
			my @spacer;
			foreach my $q (@res_pak){
				if ($q=~/[A-Z]\((\d+)\)/){push(@spacer,$1);}
			}
			my $res_composition = join(" ", @spacer);
			my $ref=$this->{template2model_site}{$key};
			my $par1="seq_query=".$fire->{psiout}{$ref}[1]."&seq_templ=".$fire->{psiout}{$ref}[3]."&pdbK=".$fire->{psiout}{$ref}[2];
			$par1.="&coor=".$fire->{psiout}{$ref}[4]."&tmp=".$fire->{tmpfile}."&residues=".$res_composition."&psiout_rank=$ref";
			print SUM "<h4 class=\"summary\"><b>SITE        $ord</b></h4>\n";
			print SUM "<table cellspacing=\"5\">\n";
			print $OUTPUT "\n<tr><td style=\"vertical-align:top; text-align:left;\" nowrap width=\"24%\">\n";
			print $OUTPUT "<b style=\"color: #FF0000;\">SITE       $ord    Compound_name:</b></td>\n";
			print $OUTPUT "<td style=\"vertical-align:top; text-align:left;\" width=\"60%\">$join->{summary}{compname}[$key]</td></tr>\n";
			print $OUTPUT "<tr><td style=\"vertical-align:top; text-align:left;\" nowrap>SITE      $ord    Site classification:</td>\n";
			print $OUTPUT "<td style=\"vertical-align:top; text-align:left;\"><b>$join->{summary}{nice_try}[$key]</b>&nbsp;\n";
			print $OUTPUT "<a href=\"../../html/firestar_help.html#Biologic\" target=\"_parent\">\n";
			print $OUTPUT "<img width=\"16px\" src=\"../../images/Question_mark.png\"></a></td></tr>\n";
			print $OUTPUT "<tr><td style=\"vertical-align:top; text-align:left;\" nowrap>SITE      $ord    Reliability:</td>\n";
			print $OUTPUT "<td style=\"vertical-align:top; text-align:left;\">$join->{summary}{score}[$key]\n";
			print $OUTPUT "<a href=\"../../html/firestar_help.html#Reliability\" target=\"_parent\">\n";
			print $OUTPUT "<img width=\"16px\" src=\"../../images/Question_mark.png\"></a></td></tr>\n";
			print $OUTPUT "<tr><td style=\"vertical-align:top; text-align:left;\" nowrap>SITE      $ord    Source alignments:</td>\n";
			print $OUTPUT "<td style=\"vertical-align:top; text-align:left;\">$join->{summary}{numOfSites}[$key]\n";
			print $OUTPUT "--> <a href=\"$fire->{tmpfile}\_sum.html\">here</a> for more information about this site prediction</td></tr>\n";
			print $OUTPUT "<tr><td style=\"vertical-align:top; text-align:left;\" nowrap>SITE      $ord    Residues resume:</td>\n";
			print $OUTPUT "<td style=\"vertical-align:top; text-align:left;\">$join->{summary}{resresume}[$key]</td>\n";
			print $OUTPUT "<td style=\"vertical-align:middle; text-align:center;\"><button onClick=\"Model_gen(this)\" value=\"$par1\" style=\"border: 0; background: transparent\">\n";
			print $OUTPUT "<img src=\"../../images/Model_button.png\">\n</button></td><td><a href=\"../../html/firestar_help.html#Model\" target=\"_parent\">\n";
			print $OUTPUT "<img align=\"center\" width=\"18px\" src=\"../../images/Question_mark.png\"></a></td></tr>\n";
			if (scalar@pharma>0){
				print $OUTPUT "\n<tr><td></td><td style=\"vertical-align:middle; text-align:left;\" width=\"60%\">\n<a href=\"$fire->{tmpfile}\_sum.html\" target=\"_parent\">\n";
				print $OUTPUT "<img align=\"center\" width=\"30px\" src=\"../../images/pharma.jpg\"></a> <-- compounds with pharmacological annotation in this site</td>\n";
				print $OUTPUT "<td style=\"vertical-align:middle; text-align:right;\">\n<a href='#Pred_Summary'><button style=\"border: 0; background: transparent\" ";
				print $OUTPUT "onClick=\"document.getElementById(\'iframe\').src=\'sites_frame/$fire->{tmpfile}\_$key\_pred.html\'\">";
				print $OUTPUT "<img src=\"../../images/Highlight_button.png\"></button></a>\n</td></tr>\n";
			}
			else {
				print $OUTPUT "<tr><td colspan=3 style=\"vertical-align:middle; text-align:right;\">\n";
				print $OUTPUT "<a href='#Pred_Summary'><button style=\"border: 0; background: transparent\" onClick=\"document.getElementById(\'iframe\')";
				print $OUTPUT ".src=\'sites_frame/$fire->{tmpfile}\_$key\_pred.html\'\"><img src=\"../../images/Highlight_button.png\"></button></a>\n</td></tr>\n";
			}
			print SUM "\n<tr><td style=\"vertical-align:top; text-align:left; width: 200px;\" nowrap>\n<b>All ligands_ID and <br>their absolute frequency:</b></td>\n";
			print SUM "<td style=\"vertical-align:middle; text-align:justify;\">$join->{summary}{compid}[$key]</td></tr>\n";
			print SUM "<tr><td style=\"vertical-align:top; text-align:left; width: 200px;\" nowrap><b>Associated GO_terms:</b></td>\n";
			print SUM "<td style=\"vertical-align:middle; text-align:left;\">$united</td></tr>\n";
			if (scalar@pharma>0){
				print SUM "\n<tr><td style=\"vertical-align:top; text-align:left; width: 200px;\" nowrap>\n<b>Compounds with pharmacological<br>annotation</b></td>\n";
				print SUM "<td style=\"vertical-align:middle; text-align:justify;\">";
				foreach my $xxx (@pharma){	print SUM "<a href='#$xxx'>$xxx</a> ";}
				print SUM "\n</td></tr>"
			}
			print $OUTPUT "<tr><td colspan=2><br><br></td></tr>";
			print SUM "</table>\n<hr style=\"border: 1px dashed #000;\" />\n";
		}
	}
	if (scalar(keys%pharma_all)>0){
####################################################################################################################################################################################		
		print SUM "<br>\n<h2 class=\"summary\">COMPOUNDS with PHARMACOLOGICAL ANNOTATION</h2>\n";
		foreach my $xxx(keys%pharma_all){
			$this->COMPOUND_PRINTER($fire,$xxx);
			print SUM "<br><br><hr style=\"border: 1px dashed #FFA345;\" />\n";
		}
	}
	if ($join->{scalar_sites} >0 && $join->{csaincr}>0){	close(SUM);}
	if ($join->{scalar_sites} > 0){
		print $OUTPUT "<tr><td colspan=3 style=\"background-color:black\"></td></tr><tr><td colspan=3 style=\"text-align:right;\">\n";
		print $OUTPUT "<a href='#Pred_Summary'><button onClick=\"document.getElementById(\'iframe\').src=\'sites_frame/$fire->{tmpfile}\_ALL.html\'\">\n";
		print $OUTPUT "REMOVE ALL HIGHLIGHTS</button></a></td></tr>\n";
	}
	if ($join->{scalar_sites} == 0){
		print $OUTPUT "<div class=\"error\"><font color=\"red\"><i>firestar</i> did not detect ligand binding sites for your query <b>$fire->{queryname}</b>.<br>\n";
		print $OUTPUT "Please check the extended <a href=\"$fire->{tmpfile}.psiall.php\" target=_blank>PSI-BLAST output</a> or the \n";
		print $OUTPUT "<a href=\"$fire->{tmpfile}.hhsall.php\" target=_blank>HHsearch output</a> files.<br>\n";
		print $OUTPUT "Although <i>firestar</i> does not have enough information to predict ligand binding sites in this case there may still be clues that\n";
		print $OUTPUT "can be obtained from the alignments between the query sequence and those FireDB templates associated with functional information.</div>\n";
	}
	print $OUTPUT "\t\t</table></pre></fieldset>\n</div>\n</div>\n</body>\n</html>";
	unlink "$fire->{faatmp}/qwe$fire->{tmpfile}.sh";
	return 0;
}

sub structural_model{
	my ($this,$fire,$join,$order)=@_;
	my $numero;
	my @extracted=keys%{$fire->{psiout}};
	for(my$i=1;$i<$join->{csaincr}+1;$i++){
		my @coords=split(";",$join->{summary_csa}{resresume}[$i]);
		my @csa_coords;
		foreach my $q (@coords){
			if ($q=~/[A-Z]\((\d+)\)/){	push(@csa_coords,$1);}
		}
		my @csa_ord_coords=sort { $a <=> $b } @csa_coords;
		my $best_ID=undef;
		for($numero=0;$numero<$#extracted;$numero++){
			if ($fire->{psiout}{$numero}[10] eq 'PSI'){	next;}       ## here we select only HHsearch ali 'cause they're better (CASP9) ##
			elsif ($fire->{psiout}{$numero}[8]=~/Identities=(\d+)%/){
				my $ID=$1;
				if ($ID > 30 && $fire->{psiout}{$numero}[4]<$csa_ord_coords[0] && $fire->{psiout}{$numero}[5]>$csa_ord_coords[$#csa_ord_coords]){
					if ($ID>$best_ID){	$this->{template2model_csa}{$i}=$numero;$best_ID=$ID;}
				}
			}
		}
		unless (defined $best_ID){
			for($numero=0;$numero<$#extracted;$numero++){
				if ($fire->{psiout}{$numero}[10] eq 'PSI'){	next;}
				elsif ($fire->{psiout}{$numero}[4]<$csa_ord_coords[0] && $fire->{psiout}{$numero}[5]>$csa_ord_coords[$#csa_ord_coords]){
					$this->{template2model_csa}{$i}=$numero;
					last;
				}
			}
		}
		unless (exists $this->{template2model_csa}{$i}){	$this->{template2model_csa}{$i}="NO_TEMPLATE";}
	}
	foreach my $incr(@{$order}){
		my @spl=split(/ /,$join->{summary_order}{$incr});
		my @site_coords;
		foreach my $key(@spl){
			my @coords=split(";",$join->{summary}{resresume}[$key]);
			foreach my $q (@coords){
				if ($q=~/[A-Z]\((\d+)\)/){	push(@site_coords,$1);}
			}
			my @site_ord_coords=sort {$a <=> $b} @site_coords;
			my $best_ID=undef;
			for($numero=0;$numero<$#extracted;$numero++){
				if ($fire->{psiout}{$numero}[10] eq 'PSI'){	next;}       ## here we select only HHsearch ali 'cause they're better (CASP9) ##
				elsif ($fire->{psiout}{$numero}[8]=~/Identities=(\d+)%/){
					my $ID=$1;
					if ($ID > 30 && $fire->{psiout}{$numero}[4]<$site_ord_coords[0] && $fire->{psiout}{$numero}[5]>$site_ord_coords[$#site_ord_coords]){
						if ($ID>$best_ID){	$this->{template2model_site}{$key}=$numero;$best_ID=$ID;}
					}
				}
			}
			unless (defined $best_ID){
				for($numero=0;$numero<$#extracted;$numero++){
					if ($fire->{psiout}{$numero}[10] eq 'PSI'){	next;}
					elsif ($fire->{psiout}{$numero}[4]<$site_ord_coords[0] && $fire->{psiout}{$numero}[5]>$site_ord_coords[$#site_ord_coords]){
						$this->{template2model_site}{$key}=$numero;
						last;
					}
				}
			}
			unless (exists $this->{template2model_site}{$key}){	$this->{template2model_site}{$key}="NO_TEMPLATE";}
		}
	}
	return 0;
}

sub iframe_generator{
	my ($this,$fire,$join)=@_;
	my $query_pos;
	my @seq=split(//,$fire->{sequence});
	unshift(@seq,"-");
	my $variables=Config::IniFiles->new(-file => $fire->{config_file});
    my $autoafmpath=$variables->val('PROGRAMS','square_test');
    # BEGIN: jmrc
    # change for APPRIS configuration
	my @autoafm=`$autoafmpath $fire->{config_file} $fire->{sequence} $fire->{sequence} $fire->{tmpfile}`;
	# END: jmrc
	my @seq_conserv=split(//,$autoafm[3]);
	unshift(@seq_conserv,"-");
	my @pockets=pocket_clustering($join);			# uses %all_bind_res
	my @driver=set_driver(length$fire->{sequence});		# los valores de este array definen el numero y longitud de las lineas
	# Display del sequence summary
	my $height;
	if (exists $join->{summary_csa}{resnum}){	$height=(scalar@driver*66)+(scalar@{$join->{summary_csa}{resnum}}*22)+(scalar@pockets*72+50);}
	else {	$height=(scalar@driver*66)+(scalar@pockets*72+50);}
	$height.="px";
	my @global_list=keys(%{$join->{pred_sites_list}});
	if ($join->{csaincr} >0){push (@global_list,$join->{csaincr});}
	## en este bucle que recorre todos los sitios cataliticos predichos se generan automaticamente las paginas que luego se pueden visualizar
	## en la summary page en el primer fieldset (en el cual se carga un iframe)dando al boton de "highlight these residues". A parte de una 
	## mejor visualizacion, se ha añadid una tabla de frecuencias normalizadas que mejor explica al usuario cuales son los residuos mas 
	## conservados en el site predicho y ayuda a explicar los colores que automaticamente se aplican a los residuos !!
	if (scalar@global_list==0 && $join->{csaincr} == 0){
		open (FILE,'>',"$fire->{faatmp}/sites_frame/$fire->{tmpfile}\_ALL.html");
		print FILE "$this->{iframe_header}\t<table class=\"stat\" border=\"0\"  cellspacing=\"0\" cellpadding=\"0\">";
		for(my $j=0;$j<$#driver+1;$j++){
			my $line=$j;
			my $offset=$driver[$j];
			print FILE "<tr><td></td><td></td><td></td>";
			for(my $i=0;$i<$offset;$i++){
				$query_pos=$line*70 +$i+1;
				my $div=$query_pos/10;
				if($div!~/\./){print FILE "<td colspan=10 align=right>$query_pos</td>";}
			}
			my $last=substr($query_pos,-1);
			if($last >5){print FILE "<td colspan=$last align=right>$query_pos</td>";}
			else{print FILE "<td colspan=$last></td>";}
			print FILE "</tr>";
			if ($j==0){print FILE "<tr><td>&nbsp&nbsp</td><td>Query:</td><td>&nbsp&nbsp</td>";}
			else {print FILE "<tr><td>&nbsp&nbsp</td><td>&nbsp&nbsp</td><td>&nbsp&nbsp</td>";}
			for(my $i=0;$i<$offset;$i++){
				$query_pos=$line*70 +$i+1;
				my $bc = escala($seq_conserv[$query_pos]);
				my $fontcolor;
				if($seq_conserv[$query_pos] >4){$fontcolor="white";}
				else{$fontcolor="black";}
				print FILE "<td title=\"Score=$seq_conserv[$query_pos]\" bgcolor=\"$bc\"><font color=\"$fontcolor\">$seq[$query_pos]</font></td>";
			}
			print FILE "</tr>\n<tr><td><br></td><td></td><td></td>";
		}
		print FILE "</table></div></body></html>";
	}
	else{
		for(my$w=0;$w<scalar(@global_list);$w++){
			open (FILE,'>',"$fire->{faatmp}/sites_frame/$fire->{tmpfile}\_$global_list[$w]\_pred.html");
			print FILE "$this->{iframe_header}\t<table><tr><td style=\"vertical-align: top;\"><table class=\"stat\" border=\"0\"  cellspacing=\"0\" cellpadding=\"0\">";
			if ($w==0){
				open (FILE2,'>',"$fire->{faatmp}/sites_frame/$fire->{tmpfile}\_ALL.html");
				print FILE2 "$this->{iframe_header}\t<table class=\"stat\" border=\"0\"  cellspacing=\"0\" cellpadding=\"0\">";
			}
			for(my $j=0;$j<$#driver+1;$j++){
				my $line=$j;
				my $offset=$driver[$j];

				print FILE "<tr><td></td><td></td><td></td>";
				if ($w==0){print FILE2 "<tr><td></td><td></td><td></td>";}
				for(my $i=0;$i<$offset;$i++){
					$query_pos=$line*70 +$i+1;
					my $div=$query_pos/10;
					if($div!~/\./){
						print FILE "<td colspan=10 align=right>$query_pos</td>";
						if ($w==0){print FILE2 "<td colspan=10 align=right>$query_pos</td>";}
					}
				}
				my $last=substr($query_pos,-1);
				if($last >5){
					print FILE "<td colspan=$last align=right>$query_pos</td>";
					if ($w==0){print FILE2 "<td colspan=$last align=right>$query_pos</td>";}
				}
				else{
					print FILE "<td colspan=$last></td>";
					if ($w==0){print FILE2 "<td colspan=$last></td>";}
				}
				print FILE "</tr>";
				if ($w==0){print FILE2 "</tr>";}
				print FILE "<tr><td>&nbsp&nbsp</td><td>Query</td><td>&nbsp&nbsp</td>";
				if ($w==0){print FILE2 "<tr><td>&nbsp&nbsp</td><td>Query</td><td>&nbsp&nbsp</td>";}
				for(my $i=0;$i<$offset;$i++){
					$query_pos=$line*70 +$i+1;
					my $bc = escala($seq_conserv[$query_pos]);
					my $fontcolor;
					if($seq_conserv[$query_pos] >4){	$fontcolor="white";}
					else{$fontcolor="black";}
					print FILE "<td title=\"Score=$seq_conserv[$query_pos]\" bgcolor=\"$bc\"><font color=\"$fontcolor\">$seq[$query_pos]</font></td>";
					if ($w==0){print FILE2 "<td title=\"Score=$seq_conserv[$query_pos]\" bgcolor=\"$bc\"><font color=\"$fontcolor\">$seq[$query_pos]</font></td>";}
				}
				print FILE "</tr>";
				if ($w==0){print FILE2 "</tr>";}
				if(scalar(keys%{$join->{all_csa}})>0){
					print FILE "<tr><td></td><td>CSA</td><td></td>";
					if ($w==0){print FILE2 "<tr><td></td><td>CSA</td><td></td>";}
					$query_pos=0;
					for(my $i=0;$i<$offset;$i++){
						$query_pos=$line*70 +$i+1;
						my $div=($query_pos-1)/5;
						my @div=split(/\./,$div);
						my $div2=$div[0]/2;
						my $bc;
						if($div2=~/\./){	$bc="#D2D2D2";}
						else{	$bc="#E2E2E2";}
						if(exists$join->{all_csa}{$query_pos}){
							print FILE "\n<td title=\"$query_pos\" bgcolor=\"green\"><font color=\"white\">$seq[$query_pos]</font></td>";
							if ($w==0){print FILE2 "\n<td title=\"$query_pos\" bgcolor=\"green\"><font color=\"white\">$seq[$query_pos]</font></td>";}
						}
						else{
							print FILE "<td bgcolor=\"$bc\">.</td>";
							if ($w==0){print FILE2 "<td bgcolor=\"$bc\">.</td>";}
						}
					}
					print FILE "</tr>";
					if ($w==0){print FILE2 "</tr>";}
				}
				my $pocket_num=0;
				for(@pockets){
					$pocket_num++;
					my @spl=split(/ /,$_);
					my %pocket_res=();
					foreach my$id(@spl){
						my @res;
						@res=split(/ /,$join->{all_bind_res}{$id});
						foreach my$res(@res){	$pocket_res{$res}="";}
					}
					my $pocket_letter=$this->{numbers_to_letters}{$pocket_num};
					print FILE "<tr><td></td><td nowrap>Pocket $pocket_letter</td><td></td>";
					if ($w==0){print FILE2 "<tr><td></td><td nowrap>Pocket $pocket_letter</td><td></td>";}
					$query_pos=0;
					for(my $i=0;$i<$offset;$i++){
						$query_pos=$line*70 +$i+1;
						my $div=($query_pos-1)/5;
						my @div=split(/\./,$div);
						my $div2=$div[0]/2;
						my $bc;
						if($div2=~/\./){$bc="#D2D2D2";}
						else{           $bc="#E2E2E2";}
						my $color_mode="OFF";
						if(exists$pocket_res{$query_pos}){
							foreach (@{$join->{pred_sites_list}{$global_list[$w]}}){	if ($_==$query_pos){$color_mode="ON";}}
							if ($color_mode eq "ON"){
								my $style=colouring($global_list[$w],$query_pos,$join);
								print FILE "\n<td $style>$seq[$query_pos]</b></font></td>";
							}
							else {	print FILE "\n<td <td bgcolor=\"$bc\">.</td>";}
							if ($w==0){print FILE2 "\n<td title=\"$query_pos\" bgcolor=\"yellow\">$seq[$query_pos]</td>";}
						}
						else{
							print FILE "<td bgcolor=\"$bc\">.</td>";
							if ($w==0){print FILE2 "<td bgcolor=\"$bc\">.</td>";}
						}
					}
					print FILE "</tr>";
					if ($w==0){print FILE2 "</tr>";}
				}
				print FILE "<tr><td>|</td><td></td><td></td>";
				if ($w==0){print FILE2 "<tr><td>|</td><td></td><td></td>";}
				$query_pos=0;
				for(my $i=0;$i<$offset;$i++){
					if ($w==0){print FILE2 "<td></td>";}
					print FILE "<td></td>";
				}
				print FILE "</tr>";
				if ($w==0){print FILE2 "</tr>";}
			}
			print FILE "$this->{iframe_table_end}\t</table></div></body></html>";
			if ($w==0){print FILE2 "</table></div></body></html>";}
			close(FILE);
			if ($w==0){close(FILE2);}
		}
	}
	return $height;

}

sub pocket_clustering{
	my ($join)=@_;
	my $first_cutoff=0.4;
	my $scnd_cutoff=0.4;
	my %finalPocket=();
	foreach my$id1(keys%{$join->{all_bind_res}}){
		$finalPocket{$id1}=$id1;
		my @res1=split(/ /,$join->{all_bind_res}{$id1});
		my $len1=scalar@res1;
		foreach my$id2(keys%{$join->{all_bind_res}}){
			unless($id1 eq $id2){
				my @res2=split(/ /,$join->{all_bind_res}{$id2});
				my $len2=scalar@res2;
				my $common=0;
				foreach my$uno(@res1){
					foreach my$dos(@res2){
						if($uno eq $dos){	$common++;}
					}
				}
				if($common > $len1 * $first_cutoff or $common > $len2 * $scnd_cutoff){
					if(exists$finalPocket{$id1}){
						$finalPocket{$id1}="$finalPocket{$id1} $id2";
					}
					else{	$finalPocket{$id1}=$id1;}
				}
			}
		}
	}
	my @cpocket=$join->id_clustering(values%finalPocket);
	return @cpocket;
}

sub set_driver{
	my @dri;
	if($_[0]>70){
		my $ln=($_[0]/70);
		my @spl=split(/\./,$ln);
		my $lines=$spl[0];
		my $resto=$_[0]-$lines*70;
		for(my$i=0;$i<$lines;$i++){	push(@dri,"70");}
		if($resto !=0){	push(@dri,$resto);}
	}
	else{	push(@dri,$_[0]);}
	return @dri;
}

sub escala{
	my $bc;
	if($_[0] eq "-"){$bc="#FFFFFF";}
	elsif($_[0] eq "0"){$bc="#EEFFFF";}
	elsif($_[0] eq "1"){$bc="#CCEEFF";}
	elsif($_[0] eq "2"){$bc="#AADDFF";}
	elsif($_[0] eq "3"){$bc="#CCCCFF";}
	elsif($_[0] eq "4"){$bc="#99AAFF";}
	elsif($_[0] eq "5"){$bc="#5599EE";}
	elsif($_[0] eq "6"){$bc="#3377EE";}
	elsif($_[0] eq "7"){$bc="#0055EE";}
	elsif($_[0] eq "8"){$bc="#0000DD";}
	elsif($_[0] eq "9"){$bc="#0000BB";}
	return $bc
}

sub colouring{
	my ($IDENT,$position,$join)=@_;
# Res:$query_pos=<b>$frequenza<\/b>;
	my @residuos=split(/;/,$join->{summary}{resfreq}[$IDENT]);
	foreach my $i (@residuos){
		$i=~s/\s+//g;
		my $pos;my $punt;
		if ($i=~/\w+\((\d+)\)=(\d\.\d+)/){	$pos=$1;$punt=$2;}
#       $pos=~s/Res://;
		if ($pos==$position){
			if ($punt<=0.40){return("title=\"$pos Prob. score=$punt\" bgcolor=\"#FFCC66\"");}
			elsif ($punt<=0.60){return("title=\"$pos Prob. score=$punt\" bgcolor=\"#FF9900\"");}
			elsif ($punt<=0.80){return("title=\"$pos Prob. score=$punt\" bgcolor=\"#FF6600\"");}
			elsif ($punt<=0.99){return("title=\"$pos Prob. score=$punt\" bgcolor=\"#FF3333\"><font color=\"white\"");}
			elsif ($punt==1.00){return("title=\"$pos Prob. score=$punt\" bgcolor=\"red\"><font color=\"white\"><b");}
		}
		else {next;}
	}
}

sub COMPOUND_PRINTER{
	my ($this,$fire,$name)=@_;
	my %pharma_info;
	if (exists $this->{man_annot}{$name}){
		my $annot=$this->{dbi}->manual_annotation(-id=>$this->{man_annot}{$name});
		push(@{$pharma_info{'MANUAL'}},$annot);
	}
	if (exists $this->{pharma}{$name}){
		foreach my $i (@{$this->{pharma}{$name}}){
			my @annot=$this->{dbi}->automatic_annotation(-ext_id=>$i);
			my $DB=shift@annot;
			@{$pharma_info{$DB}}=@annot;
		}
	}
	my $variables=Config::IniFiles->new(-file => $fire->{config_file});
	print SUM "\n\n<br><br><br><a name=\"$name\">\n<table class=\"compounds\">\n<tr><td colspan=4 align=\"center\"><br>PDB COMPOUND <b>$name</b>  ";
	print SUM "<a href=\"",$variables->val('EXT_LINKS','PDB_page'),"$name\">PDB web page</a><br><br></td></tr>\n";
	print SUM "<td class=\"image_compo\" rowspan=\"",scalar(keys%pharma_info),"\"><img width=\"250px\" height=\"250px\" src=\"";
	print SUM $variables->val('EXT_LINKS','PDB_image'),"$name\_620.gif\"/></td>\n";
	foreach my $i(keys%pharma_info){
		if ($i eq "MANUAL"){print SUM "<td class=\"name_daba\">MANUAL</td><td class=\"ext_link\">$name</td><td class=\"description\">${$pharma_info{$i}}[0]";}
		else{
			my $ext_link=shift@{$pharma_info{$i}};
			my $stereo=shift@{$pharma_info{$i}};
			if ($stereo eq "SAME"){$stereo="";}
			else {$stereo="<br>(diff. Stereochem.)";}
			print SUM "<td class=\"name_daba\">$i</td><td class=\"ext_link\"><a href=\"",$variables->val('EXT_LINKS',$i),"$ext_link\">$ext_link</a>$stereo</td>";
			print SUM "<td class=\"description\">",shift@{$pharma_info{$i}};
			foreach my $k(@{$pharma_info{$i}}){
				print SUM "<br>$k";
			}
		}
		print SUM "</td></tr>\n";
	}
	print SUM "</table>\n";
}


sub web_headers{
	my ($this,$fire,$join)=@_;
	$this->{summary_header}="<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n<html>\n<head>\n";
	$this->{summary_header}.="<title>firestar</title>\n<meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\" />\n";
	$this->{summary_header}.="<style type=\"text/css\" media=\"all\">\@import \"../../CSS/style.css\";</style>\n<link rel=\"icon\" type=\"image/x-icon\" href=\"../../images/mini_logo.png\"/>\n";
	$this->{summary_header}.="<link rel=\"stylesheet\" type=\"text/css\" href=\"../../JS/yui/build/fonts/fonts-min.css\" />\n<script type=\"text/javascript\" ";
	$this->{summary_header}.="src=\"../../JS/yui/build/yahoo/yahoo-min.js\"></script>\n<script type=\"text/javascript\" src=\"../../JS/yui/build/event/event-min.js\"></script>\n";
	$this->{summary_header}.="<script src=\"../../JS/yui/build/connection/connection_core-min.js\"></script>\n\n<script type=\"text/javascript\" ";
	$this->{summary_header}.="src=\"../../JS/yui/build/connection/connection-min.js\"></script>\n\n<link rel=\"stylesheet\" type=\"text/css\" ";
	$this->{summary_header}.="href=\"../../JS/yui/build/container/assets/skins/sam/container.css\" />\n<script type=\"text/javascript\" ";
	$this->{summary_header}.="src=\"../../JS/yui/build/yahoo-dom-event/yahoo-dom-event.js\"></script>\n<script type=\"text/javascript\" ";
	$this->{summary_header}.="src=\"../../JS/yui/build/dragdrop/dragdrop-min.js\"></script>\n<script type=\"text/javascript\" src=\"../../JS/yui/build/container/container-min.js\"></script>\n";
	$this->{summary_header}.="<script type=\"text/javascript\" src=\"../../JS/yui/build/connection/connection-min.js\"></script>\n";
	$this->{summary_header}.="<script type=\"text/javascript\" src=\"../../JS/yui/build/animation/animation-min.js\"></script>";
	$this->{summary_header}.="<script type=\"text/javascript\" src=\"../../JS/yui/build/dragdrop/dragdrop-min.js\"></script>";
	$this->{summary_header}.="<script language=\"JavaScript\" src=\"../../Php/Functions.js\"></script>\n</head>\n<body class=\"yui-skin-sam\">\n\t<div id=\"topbar\">\n\t\t<h1>\n";
	$this->{summary_header}.="\t\t\t<table width=\"100%\"><tr><td widht=\"100%\" align=\"left\"><h1><a href=\"../../index.php\">Home</a></h1></td><td widht=\"100%\" align=\"right\">\n";
	$this->{summary_header}.="\t\t\t<a href=\"../../Php/FireStar.php\">Back to <i>firestar</i></a></td></tr></table>\n\t\t</h1>\n\t</div>\n\t";
	$this->{summary_header}.="<div id=\"tabs2\">\t<ul>\t<li><span><i>firestar</i></span></li>\t</ul>\t</div>\n\t<div id=\"main\">\n\t\t<div id=\"subtabs2\">\n\t\t\t&nbsp;\n\t\t</div>\n";
	$this->{summary_header}.="\t\t<div id=\"bodyarea\">\n\n";
	################
	$this->{summary_table_header}="<table><tr><td align=justify>\n<h1><i>firestar</i> server: predicting functional residues from structural templates and alignment reliability</h1>\n";
	$this->{summary_table_header}.="<a name='Pred_Summary'><p><i>firestar</i></a> predicts functionally important residues. Predictions are made by evaluating<br>\n";
	$this->{summary_table_header}.="internally generated alignments using <a href=\"../../Php/square.php\" target=_blank>SQUARE</a></p>\n";
	$this->{summary_table_header}.="Results shown in this page are a consensus summary of the EXTENDED RESULTS (see below).<br><br>\n";
	$this->{summary_table_header}.="<b>REMEMBER: Results will be stored in our server for 3 weeks, than will be deleted !!<br><br>\n";
	$this->{summary_table_header}.="For any doubt about the results format, go to <a href=\"../../html/firestar_help.html\" target=_blank color=red>Help page</a></b></p>\n";
	$this->{summary_table_header}.="</td><td width=\"10%\"></td><td><img src=\"../../images/firestar.png\" width=\"180\"></td></tr></table><br>\n\n";
	################
	$this->{extended_results_header}="<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n<html>\n<head>\n";
	$this->{extended_results_header}.="<title>firestar</title>\n<meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\" />\n<style type=\"text/css\" media=\"all\">";
	$this->{extended_results_header}.="\@import \"../../CSS/style.css\";</style>\n<link rel=\"icon\" type=\"image/x-icon\" href=\"../../images/mini_logo.png\"/>\n";
	$this->{extended_results_header}.="<script language=\"JavaScript\" src=\"../../Php/Functions.js\"></script>\n</head>\n<body>\n\t<div id=\"topbar\">\n\t\t\t<h1>\n\t\t\t";
	$this->{extended_results_header}.="<table width=\"100%\"><tr><td widht=\"100%\" align=\"left\"><h1><a href=\"../../index.php\">Home</a></h1></td><td widht=\"100%\" align=\"right\">\n\n";
	$this->{extended_results_header}.="\t\t\t<a href=\"http://firedb.bioinfo.cnio.es/Php/FireStar.php\">Back to <i>firestar</i></a></td></tr></table>\n\t\t\t</h1>\n\t</div>\n\t";
        $this->{extended_results_header}.="<div id=\"tabs2\">\t<ul>\t<li><span><i>firestar</i></span></li>\t</ul>\t</div>\n\t\t<div id=\"main\">\n\t\t<div id=\"subtabs2\">\n\t\t\t&nbsp;\n\t\t";
	$this->{extended_results_header}.="</div>\n\t\t<div id=\"bodyarea\">\n\n";
	################
	$this->{square_score_text}="Note that sites are usually formed by groups of residues. The user should assess overall conservation, based on all the potential residues important for"; 
	$this->{square_score_text}.="one site. It is common to find binding sites partially conserved.";
	################
	$this->{occurrence_text}="Occurrence is calculated when a given consensus sequence has several representatives in the PDB. It indicates the percentage of PDB representatives that";
	$this->{occurrence_text}.=" bind ligand analogs at same site. Where there are less than 5 representatives the values are marked \"not applicable\".";
	################
	$this->{score_tables}="<h1><i>firestar</i> server: function prediction based on Structural Templates and Alignment Reliability</h1>\n";
	$this->{score_tables}.="<table><tr><td colspan=2 bgcolor=\"#555555\"><font color=white>SQUARE: table of scores</font></td bgcolor=\"#555555\"><td bgcolor=\"#555555\">&nbsp</td>\n";
	$this->{score_tables}.="<td bgcolor=\"#555555\"colspan=2><font color=white>Binding site occurrence and Catalytic Site evidence codes</font></td></tr>";
	$this->{score_tables}.="<tr><td>-&nbsp</td><td bgcolor=white>Gapped or non-conserved position</td><td bgcolor=\"#555555\">&nbsp</td><td bgcolor=\"#AAAAAA\"><font color=green>";
	$this->{score_tables}.="Literature</font></td><td bgcolor=white>Catalytic Site annotated from Literature</td></tr>\n";
	$this->{score_tables}.="<tr><td bgcolor=\"#E0FFFF\">1&nbsp</td><td bgcolor=white>45% reliability</td><td bgcolor=\"#555555\">&nbsp</td><td bgcolor=\"#AAAAAA\"><font color=yellow>";
	$this->{score_tables}.="PSI-blast</font></td><td bgcolor=white>Catalytic Site inferred by similarity</td></tr>\n";
	$this->{score_tables}.="<tr><td bgcolor=\"#B0E2FF\">2&nbsp</td><td bgcolor=white>60% reliability</td><td bgcolor=\"#555555\">&nbsp</td><td bgcolor=\"white\"><td bgcolor=white></td></tr>\n";
	$this->{score_tables}.="<tr><td bgcolor=\"#55DDFF\">3&nbsp</td><td bgcolor=white>75% reliability</td><td bgcolor=\"#555555\">&nbsp</td><td bgcolor=\"#AAAAAA\"><font color=red>";
	$this->{score_tables}.="0-20% occurrence</font></td><td bgcolor=white></td></tr>\n";
	$this->{score_tables}.="<tr><td bgcolor=\"#1E90FF\"><font color=\"white\">4</font>&nbsp</td><td bgcolor=white>85%  reliability</td><td bgcolor=\"#555555\">&nbsp</td><td bgcolor=\"#AAAAAA\">";
	$this->{score_tables}.="<font color=yellow>20-50% occurrence</font></td><td bgcolor=white></td></tr>\n";
	$this->{score_tables}.="<tr><td bgcolor=\"#1874CD\"><font color=\"white\">5</font>&nbsp</td><td bgcolor=white>90%  reliability</td><td bgcolor=\"#555555\">&nbsp</td><td bgcolor=\"#AAAAAA\">";
	$this->{score_tables}.="<font color=green>50-100% occurrence</font></td><td bgcolor=white></td></tr>\n";
	$this->{score_tables}.="<tr><td bgcolor=\"darkblue\"><font color=\"white\">C</font>&nbsp</td><td bgcolor=white>99%  reliability</td><td bgcolor=\"#555555\">&nbsp</td><td bgcolor=\"#AAAAAA\">";
	$this->{score_tables}.="<font color=black>XX% occurrence</font></td><td bgcolor=white>Not applicable</td></tr>\n";
	$this->{score_tables}.="<tr><td  bgcolor=\"#555555\" colspan=2></td><td bgcolor=\"#555555\"></td><td colspan=2  bgcolor=\"#555555\"></td></tr>\n";
	$this->{score_tables}.="<tr><td bgcolor=white colspan=2>$this->{square_score_text}</td><td bgcolor=\"#555555\">&nbsp</td><td colspan=2 bgcolor=white>$this->{occurrence_text}</td></tr>\n";
	$this->{score_tables}.="</table><hr>\n\n";
	################
	$this->{exResults_table_header}="<hr>\n<table>\n<tr>\n<td colspan=2 bgcolor=\"#AAAAAA\"></td>\n<td bgcolor=\"#AAAAAA\"></td>\n<td bgcolor=\"#AAAAAA\"></td>\n</tr>\n<tr>";
	$this->{exResults_table_header}.="<td bgcolor=\"#AAAAAA\"></td>\n<td bgcolor=\"pink\" style=\"text-align:center; width:270px\"><button type=\"button\" onclick=\"MULT_Align()\">";
	$this->{exResults_table_header}.="DISPLAY: Multiple Sequence Alignment</button></td>\n";
	$this->{exResults_table_header}.="<td bgcolor=\"pink\">Here multiple sequence alignments can be generated with <a href=\"https://www.ebi.ac.uk/Tools/msa/muscle/\" target=_blank>MUSCLE</a> ";
	$this->{exResults_table_header}.="selecting your protein(s) of interest marking the checkbox on the left (<b>AT LEAST</b> one).\n"; 
	$this->{exResults_table_header}.="Important residues can be highlighted in the Multiple Alignment. For more information see the <a ";
	$this->{exResults_table_header}.="href=\"http://firedb.bioinfo.cnio.es/html/firestar_help.html#msa\">help page</a>.\n</td>\n<td bgcolor=yellow>Yellow boxes allow ";
	$this->{exResults_table_header}.="generating Structural Alignments with LGA. See the <a href=../../html/firestar_help.html#straln target=_blank>Online Help</a> for more information.</td>\n";
	$this->{exResults_table_header}.="</tr>\n<tr>\n<td colspan=2 bgcolor=\"#AAAAAA\"></td>\n<td bgcolor=\"#AAAAAA\"></td>\n<td bgcolor=\"#AAAAAA\"></td>\n</tr>\n</table><hr>\n\n";
	################
	$this->{struct_alig_header1}="<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n<html>\n<head>\n";
	$this->{struct_alig_header1}.="\n<title>firestar</title>\n";
	################
	$this->{struct_alig_header2}="<meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\" />\n<style type=\"text/css\" media=\"all\">\@import \"../../CSS/style.css\";</style>";
	$this->{struct_alig_header2}.="\n<link rel=\"icon\" type=\"image/x-icon\" href=\"../../images/mini_logo.png\"/>\n</head>\n<h1><table width=\"100%\"><tr><td width=\"100%\" align=\"left\">";
	$this->{struct_alig_header2}.="<h1>Results</h1></td><td width=\"100%\" align=\"right\" nowrap><a href=\"../../Php/FireStar.php\">Back to <i>firestar</i></a></td></tr></table></h1><br>\n";
	$this->{struct_alig_header2}.="<body>\n\t<div id=\"tabs2\">\n\t\t<ul>\t<li><span><i>firestar</i></span></li>\t</ul>\n\t</div>\n\t<div id=\"main\">\n\t\t<div id=\"subtabs2\">\n\t\t\t&nbsp;\n";
	$this->{struct_alig_header2}.="\t\t</div>\n\t\t<div id=\"bodyarea\">\n<?php\n";
	################
	$this->{struct_alig_php_form}="\$q_start=\$_POST[\"q_start\"];\n\$q_end=\$_POST[\"q_end\"];\n\$t_start=\$_POST[\"t_start\"];\n\$t_end=\$_POST[\"t_end\"];\n";
	$this->{struct_alig_php_form}.="\$arguments=\"\$tmpfile \$querytype \$filetype \$template \$target \$target_chain \$q_start \$q_end \$t_start \$t_end\";\n";
	$this->{struct_alig_php_form}.="\$comando=\"../../perl/fireStrAln.pl \$arguments\";\n\$programm=\"\n#!/bin/sh\n../../perl/fireStrAln.pl \$tmpfile \$querytype ";
	$this->{struct_alig_php_form}.="\$filetype \$template \$target \$target_chain \$q_start \$q_end \$t_start \$t_end > ../straln/\$tmpfile.write.php\n";
	$this->{struct_alig_php_form}.="mv -f ../straln/\$tmpfile.write.php ../../tmp/straln/\$tmpfile.php\n\";\n\$open= fopen(\"../straln/\$tmpfile.sh\", \"w\");\n";
	$this->{struct_alig_php_form}.="\$write = fwrite(\$open, \$programm);\nfclose (\$open);\n\$var=shell_exec(\"../straln/\$tmpfile.sh\");\n?>\n\t\t</div>\n\t</div>\n</body>\n</html>\n";
	################
	$this->{struct_alig_wait_part1}="<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n<html>\n<head>\n";
	$this->{struct_alig_wait_part1}.="\n<title>firestar</title>\n<meta http-equiv=\"refresh\" content=\"10\">\n<style type=\"text/css\" media=\"all\">\@import \"../../CSS/style.css\";</style>";
	$this->{struct_alig_wait_part1}.="\n<link rel=\"icon\" type=\"image/x-icon\" href=\"../../images/mini_logo.png\"/>\n</head>\n<body>\n\t<div id=\"tabs2\">\n";
	$this->{struct_alig_wait_part1}.="\t\t<ul>\t<li><span><i>firestar</i></span></li>\t</ul>\n\t</div>\n\t<div id=\"main\">\n\t\t<div id=\"subtabs2\">\n\t\t\t&nbsp;\n\t\t</div>\n\t\t";
	$this->{struct_alig_wait_part1}.="<div id=\"bodyarea\">\n\t\t<div class=note>\n\t\tYour query is being processed. This page will be automatically reloaded <br>\n\t\tResults will be displayed in ";
	################
	$this->{struct_alig_wait_part2}=", copy the link location if you wish to close this window\n\t\t</div>\n\t\t<table align=center>";
	$this->{struct_alig_wait_part2}.="\n\t\t<tr><td bgcolor=white><br>&nbsp;&nbsp;\n\t\t<img width=\"160\" src=\"../../images/searchBar.gif\" align=center>\n\t\t&nbsp;&nbsp;<br></td></tr>\n";
	$this->{struct_alig_wait_part2}.="\t\t</table>\n\t\t</div>\n\t</div>\n</body>\n</html>\n";
	################
	$this->{download_psi_output}="<div class=note><font align=right>Download PSI-BLAST <a href=\"$fire->{tmpfile}\_$fire->{evalue}.psi\" target=_blank>OUTPUT FILE</a></font></div><hr>";
	################
	$this->{download_hhs_output}="<div class=note><font align=right>Download HHsearch <a href=\"$fire->{tmpfile}.hhr\" target=_blank>OUTPUT FILE</a></font></div><hr>";
	################
	$this->{no_results_psi_output}="<div class=\"error2\">There are no hits from PSI-BLAST search. You can check directly \nin this <a href=\"$fire->{tmpfile}\_$fire->{evalue}.psi\""; 
	$this->{no_results_psi_output}.=" target=_blank>page</a> the PSI-BLAST output.";
	################
	$this->{no_results_hhs_output}="<div class=\"error2\">There are no hits from HHsearch analysis. You can check directly \nin this <a href=\"$fire->{tmpfile}.hhr\"";
	$this->{no_results_hhs_output}.=" target=_blank>page</a> the HHsearch output.";
	################
	$this->{iframe_header}="<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n";
	$this->{iframe_header}.="<html style=\"background-color:#F2F2F2;\">\n<head>\n\n<title>FireDB</title>\n<meta http-equiv=\"Content-Type\" ";
	$this->{iframe_header}.="content=\"text/html; charset=ISO-8859-1\" />\n<style type=\"text/css\" media=\"all\">\@import \"../../../CSS/style.css\";</style>\n";
	$this->{iframe_header}.="<link rel=\"icon\" type=\"image/x-icon\" href=\"../../../../images/mini_logo.png\"/>\n</head>\n<body style=\"width:800px; background-color:#F2F2F2;\">\n";
	$this->{iframe_header}.="<h4>Query name&nbsp;&nbsp;&nbsp;&nbsp;$fire->{queryname}</h4>\n";
	################
	$this->{iframe_table_end}="</table></td><td></td><td><a href=\"../../../html/firestar_help.html#Probab\" target=\"_parent\">\n\t\t<img align=\"right\" width=\"18px\" ";
	$this->{iframe_table_end}.="src=\"../../../images/Question_mark.png\"></a><br>\n\t\t<fieldset vertical-align=\"center\"><legend>Legend</legend>\n\t\t<table align=\"center\"><tbody>";
	$this->{iframe_table_end}.="\n\t\t<tr><td bgcolor=\"red\"><font color=\"white\"><b>a.a.</b></td><td width=\"1px\"></td><td>SCORE=1.00</td></tr>\n\t\t";
	$this->{iframe_table_end}.="<tr><td bgcolor=\"#FF3333\"><font color=\"white\">a.a.</td><td width=\"1px\"></td><td>0.81 - 0.99</td></tr>\n\t\t";
	$this->{iframe_table_end}.="<tr><td bgcolor=\"#FF6600\">a.a.</td><td width=\"1px\"></td><td>0.61 - 0.80</td></tr>\n\t\t<tr><td bgcolor=\"#FF9900\">a.a.</td><td width=\"1px\"></td>";
	$this->{iframe_table_end}.="<td>0.41 - 0.60</td></tr>\n\t\t<tr><td bgcolor=\"#FFCC66\">a.a.</td><td width=\"1px\"></td><td nowrap> SCORE < 0.40</td></tr>\n";
        $this->{iframe_table_end}.="\n\t\t</tbody></table>\n\t\t</fieldset></td></tr>\n\t\t<tr><td></td></tr>\n";
	################
	$this->{summary_page_iframe_end}="</iframe></fieldset></div>\n<br>\n<fieldset style=\"width:97%;\"><legend>Text summary</legend>\n<div id=\"lolapalooza\">\n<br><i>firestar</i> results<br>\n";
	$this->{summary_page_iframe_end}.="================<br>\n";
	################
	$this->{supplementary_page_header}="<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n";
	$this->{supplementary_page_header}.="<html style=\"background-color:#F2F2F2;\">\n<head>\n<title>firestar suppl. info</title>\n";
	$this->{supplementary_page_header}.="<meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\" />\n<style type=\"text/css\" ";
	$this->{supplementary_page_header}.="media=\"all\">\@import \"../../CSS/summary.css\";</style>\n<link rel=\"icon\" type=\"image/x-icon\" href=\"../../../images/mini_logo.png\"/>\n";
	$this->{supplementary_page_header}.="</head>\n<body class=\"summary\">\n\t<h1 class=\"summary\"><i>firestar</i> SUPPLEMENTARY info page</h1><img id=\"logo\" src=\"../../images/firestar.png\">";
        $this->{supplementary_page_header}.="\n\t<br><br><h2 class=\"summary\">QUERY: $fire->{queryname}<h2>\n\t<p id=\"sequence\">$fire->{sequence}</p><br><hr>\n\t";
	################
	$this->{beginning_firestar_prediction}="<table width=100%><tbody>\n\t<tr><td style=\"text-align:justify;\"><i>firestar</i> predicts a total of ";
	$this->{beginning_firestar_prediction}.="<font color=\"blue\"><b>$join->{scalar_sites}</b></font> \n\tbinding sites ordered by number of homologs found. Sites are related to the pockets in ";
	$this->{beginning_firestar_prediction}.="the summary figure.\n\tSome pockets may be made up of more than one site and sites themselves may overlap.<br><br>\n\t";
	$this->{beginning_firestar_prediction}.="<font color=\"red\">Note that small molecule ligands [in particular Na, Cl-, K and SO4--] might not be cognate ligands.</font></td>\n\t";
	$this->{beginning_firestar_prediction}.="<td style=\"width: 100px;\"></td>\n\t<td bgcolor=\"#aaddff\" style=\"text-align: justify; width: 450px; vertical-align: top; padding: 15px;\">\n\t";
	$this->{beginning_firestar_prediction}.="<b>See detailed information on the alignments and templates <i>firestar</i> used to obtain predictions in the:<ul><li>EXTENDED\n\t";
	$this->{beginning_firestar_prediction}.="<a href=\"$fire->{tmpfile}.psiall.php\" target=_blank>PSI-BLAST RESULTS</a> page</li><li>EXTENDED\n\t";
	$this->{beginning_firestar_prediction}.="<a href=\"$fire->{tmpfile}.hhsall.php\" target=_blank>HHsearch RESULTS</a> page</li></ul></b></td></tr>\n\t</tbody></table><br><hr>\n";
	################
	%{$this->{numbers_to_letters}}=(1=>"A",2=>"B",3=>"C",4=>"D",5=>"E",6=>"F",7=>"G",8=>"H",9=>"I",10=>"J",11=>"K",12=>"L",13=>"M",14=>"N",15=>"O",16=>"P",17=>"Q",18=>"R",19=>"S",
					20=>"T",21=>"U",22=>"V",23=>"W",24=>"X",25=>"Y",26=>"Z",27=>"AA",28=>"AB",29=>"AC",30=>"AD",31=>"AE",32=>"AF",33=>"AG",34=>"AH",35=>"AI",36=>"AJ",
					37=>"AK",38=>"AL",39=>"AM",40=>"AN",41=>"AO",42=>"AP",43=>"AQ",44=>"AR",45=>"AS",46=>"AT",47=>"AU",48=>"AV",49=>"AW",50=>"AX",51=>"AY",52=>"AZ",
					53=>"BA",54=>"BB",55=>"BC",56=>"BD",57=>"BE",58=>"BF",59=>"BG",60=>"BH");
	return 0;
}

sub GO_printer{
	my ($this,$content,$list,$tag)=@_;
	foreach my $kkk (keys%{$list}){
		unless (exists $this->{already_seen_GOs}{$kkk}){
			$this->{already_seen_GOs}{$kkk}=5;
			if ($tag ne "UNKNOWN"){	push(@{$content},"$kkk\t$tag\n");}
		}
	}
}

sub mail_sender{
	my ($this,$fire)=@_;
	require Mail::Sender;
	my $remit="firestar\@bioinfo.cnio.es";
	my $asunto = "firestar results are now available";
	my $cuerpo = " 'Your results for your query $fire->{queryname} are now available at http://firedb.bioinfo.cnio.es/tmp/faatmp/$fire->{tmpfile}.php \n";
	$cuerpo.="\n\nPlease remember that results will be stored in our server only for 3 weeks, than will be deleted !!\n";
	$cuerpo.="\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	$cuerpo.="FireDB: a compendium of biological and pharmacologically relevant ligands. Maietta P, Lopez G, Carro A, Pingilley BJ, Leon LG, Valencia A, Tress ML\n";
	$cuerpo.="\thttp://nar.oxfordjournals.org/content/42/D1/D267.abstract\n\n";
	$cuerpo.="FireDB-a database of functionally important residues from proteins of known structure. Lopez G, Valencia A, Tress M;\n";
	$cuerpo.="\thttp://nar.oxfordjournals.org/content/35/suppl_1/D219.abstract\n\n";
	$cuerpo.="firestar-prediction of functionally important residues using structural templates and alignment reliability.";
	$cuerpo.=" Lopez G, Valencia A, Tress M;\n";
	$cuerpo.="\thttp://nar.oxfordjournals.org/content/35/suppl_2/W573.full?keytype=ref&ijkey=vMKUHaDnlVj9v0G\n\n";
	$cuerpo.="firestar-advances in the prediction of functionally important residues. Lopez G, Maietta P, Rodriguez JM, Valencia A and Tress ML;\n";
	$cuerpo.="\thttp://nar.oxfordjournals.org/content/39/suppl_2/W235.full\n'";
	my $sender=new Mail::Sender {smtp=>'flash.cnio.es',from => $remit,to => $fire->{mail},replyto => 'pmaietta@cnio.es'};
	$sender->MailMsg({subject=>$asunto,msg=>$cuerpo});
}




sub DESTROY {
        my ($this)=@_;
        $this->{dbi}->close();
}



1;
