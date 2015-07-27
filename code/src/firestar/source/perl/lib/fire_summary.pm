package fire_summary;

use strict;
use FindBin;
my $cwd;
BEGIN{
	$cwd=$FindBin::Bin;
}

use lib $cwd;
use fire_libraries;
use DBI_firestar;
use File::Copy;
use POSIX;

sub new {
	my($class,%args)=@_;
	my $this={};
	$this->{dbi}=DBI_firestar->connect($args{-config});
	bless($this);
	return $this;
}

sub CSA_summary{
	my ($this,$fire)=@_;
#~~~~~~~~~~	Here we collapse all the overlapping CSA annotations
	my ($all_csa_res,$all_csa_score,$csa_res_freq)=$this->collapse_overlapping_csa($fire);
	foreach my $key_res(keys%{$this->{csa_sites}}){
		my %csa;
		my @spl=split(/ /,$key_res);
		my @csiteids=split(/ /,$this->{csa_sites}{$key_res});
		foreach my $csiteid(@csiteids){
			my @new_spl=split(/ /,$all_csa_res->{$csiteid});
			my @scores=split(//,$all_csa_score->{$csiteid});
			for(my$i=0;$i<$#new_spl+1;$i++){
				my $num=$new_spl[$i];
				if ((exists$csa{$num} && $scores[$i]>$csa{$num}) || !exists$csa{$num}){$csa{$num}=$scores[$i];}
			}
			my $clustid=$fire->{all_csite_info}{$csiteid};  ## solo para CSA $all_csite_info{$csiteid} contiene el clustid; normalmente hay compuestos 
			my @csa_note=$this->{dbi}->CSA_LIT(-csite=>$csiteid);
			$this->{csa_clust}{$clustid}=$csa_note[0];
		}
		my $query_pos=0;
		$this->{csaincr}++;
		my @seq=split(//,$fire->{sequence});
		for(@seq){                      ## en este bucle se pasa por toda la secuencia y se guarda en el hash %summary_csa toda la informacion de
			if($_ ne "-"){          ## posicion, nombre, etc etc. 
				$query_pos++;
				if(exists$csa{$query_pos}){
					if (exists $this->{summary_csa}{resnum}[$this->{csaincr}]){
						$this->{summary_csa}{resnum}[$this->{csaincr}].=" $query_pos";
						$this->{summary_csa}{resname}[$this->{csaincr}].=$_;
						$this->{summary_csa}{resresume}[$this->{csaincr}].=" $_($query_pos);";
						$this->{summary_csa}{score}[$this->{csaincr}].=$csa{$query_pos};
						$this->{all_csa}{$query_pos}=5;
					}
					else {
						$this->{summary_csa}{resnum}[$this->{csaincr}]=$query_pos;
						$this->{summary_csa}{resname}[$this->{csaincr}]=$_;
						$this->{summary_csa}{score}[$this->{csaincr}]=$csa{$query_pos};
						$this->{all_csa}{$query_pos}=5;
						$this->{summary_csa}{resresume}[$this->{csaincr}]="$_($query_pos);";
					}
				}
			}
		}
		if ($key_res ne $this->{summary_csa}{resnum}[$this->{csaincr}]){
			$this->{csa_sites}{$this->{summary_csa}{resnum}[$this->{csaincr}]}=$this->{csa_sites}{$key_res};
			delete$this->{csa_sites}{$key_res};
		}
        }
}


sub binding_prediction_summary{
	my ($this,$fire,$lib,$pulser)=@_;
#~~~~~~~~~~	The residue_clustering method collapses all the single predictions that overlap. 
#~~~~~~~~~~	Please notice that we are collapsing independently metal - nonmetal - nonmetalnoncognate sites
	my @clustered_csites=$this->residue_clustering(%{$fire->{all_met_res}});
	push(@clustered_csites,$this->residue_clustering(%{$fire->{all_nom_res}}));
	push(@clustered_csites,$this->residue_clustering(%{$fire->{all_nom_nocog_res}}));
	my $hitNumOfSites=0;
	for(@clustered_csites){
		my @spl=split(/ /,$_);
		my $tmpcases=scalar@spl;
		if($hitNumOfSites < $tmpcases){$hitNumOfSites = $tmpcases;}
	}
	my $numbering_filtered_residues;
	$fire->{longest_ali}=ceil($fire->{longest_ali}/2);
	my $autoincr=0;
	for(@clustered_csites){
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# cut-offs manually established
		my $res_hit_gold60=1;
		my $res_hit_gold35=1;
		my $res_cut_off=0.25;
		my $res_hit=1;
		my $hit_mean_cutoff=1;                          # pej 1
		my $hit_score_cutoff=3;
		my $relative_num_of_cases=0;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		my $ids=$_;
		my $hit_score=0;
		my $hit_mean=0;
		my @spl=split(/ /,$_);
		my $numberOFsites=scalar@spl;
		my %all_aa=();
		my %compids=();
		my %res;
		my %all_scorez=();
		my $posit_score;
		my @position_array;
		my @gold_duo;
		my @biggest_length_ever;			###     TRASH
		my $site_length;				###     TRASH
		my %global_SQUARE_best_score;			###     TRASH
		my $FINAL_SQUARE_SCORE_COMBO;			###     TRASH
		my $global_clus_pred_tag;			###	TRASH
		my $coverage_tag;
		my %well60;
		my %well35;
		my @score_identit;
		foreach my$id(@spl){
			my @alignment=split (/_/,$id);
			push (@position_array,$alignment[1]);
			my @res;
			if(exists$fire->{all_met_res}{$id}){@res=split(/ /,$fire->{all_met_res}{$id});$coverage_tag="MET";}
			elsif(exists$fire->{all_nom_res}{$id}){@res=split(/ /,$fire->{all_nom_res}{$id});$coverage_tag="NO_MET";}
			elsif(exists$fire->{all_nom_nocog_res}{$id}){@res=split(/ /,$fire->{all_nom_nocog_res}{$id});$coverage_tag="NO_MET";}
			my @tempora=split(/ /,$fire->{score_app}{$id});
			my %scorez;
			foreach (@tempora){
				if ($_=~/(\d+)\((\d)\)/){$scorez{$1}=$2;}
			}
			my @temporal=split(/_/,$id);
			my $alig=$temporal[1];
			my $length_ali=($fire->{psiout}{$alig}[5]-$fire->{psiout}{$alig}[4]+1);
		#~~~~~~~~~~	This method selects the residues that appear in alignments with high %ID
		#~~~~~~~~~~	if present. In this case it also raises the residue frequency cut_off, $res_hit_gold60 and $res_hit_gold35
			
			weighter($fire,\$res_cut_off,$alig,$length_ali,\@res,\$res_hit_gold60,\$res_hit_gold35,\%well60,\%well35,\@score_identit);
			foreach my$res(@res){
				if(exists$res{$res}){	### generamos un hash %res con CLAVE: pos absoluta del residuo y VALOR: el numero de veces que aparece
					if ($length_ali>=$fire->{longest_ali}){$res{$res}++;}  ## en los clusters 
					if($res{$res}>$res_hit){$res_hit=$res{$res};}	## $res_hit es un umbral del numero de veces que aparece un residuo
				}
				elsif ($length_ali>= $fire->{longest_ali}){$res{$res}=1;}
				if(exists$res{$res} && $all_scorez{$res}<$scorez{$res}){$all_scorez{$res}=$scorez{$res};}
				else{$all_scorez{$res}=$scorez{$res};}
				if (exists $global_SQUARE_best_score{$res}){
					if ($global_SQUARE_best_score{$res}<$fire->{global_SQUARE_MEAN}{$id}{$res}){
						$global_SQUARE_best_score{$res}=$fire->{global_SQUARE_MEAN}{$id}{$res};
					}
				}
				else {
					$global_SQUARE_best_score{$res}=$fire->{global_SQUARE_MEAN}{$id}{$res};
				}
			}
			my @comps=split(/ /,$fire->{all_csite_info}{$id});	### aqui recuperamos todos los compuestos que hacen binding con el site
			if ($global_clus_pred_tag ne "YES"){$global_clus_pred_tag=$fire->{all_csite_TAG}{$id};}
			foreach my$compid(@comps){
				my @spl=split(/\(/,$compid);
				$spl[1]=~s/\)//;
				$compid=$spl[0];
				$compids{$compid}++;	### here we are storing the compound frequency among all the sites
			}
			if($hit_score < $fire->{all_csite_mean}{$id}){		## aqui guardamos la media del score de SQUARE de los aa del site mas alta entre
				$hit_score = $fire->{all_csite_mean}{$id};	## todos los sitios clusterizados juntos !!
				$gold_duo[0]=$hit_score;
			}
			push(@biggest_length_ever,$fire->{all_csite_coverage}{$id});
			if($hit_mean < $fire->{all_csite_score}{$id}){		## aqui guardamos el score mas alto entre los sites, calculado antes como la media de 
				$hit_mean = $fire->{all_csite_score}{$id};	## los scores de los aa del site/media de los scores de todos los aa de la sec
			}						## en el alineamiento de donde se saco el site
		}
		if (scalar(keys%res)==0){next;}
		my $corrector_HHS;
		foreach my $k (sort { $a <=> $b} keys %{$fire->{psiout}}){
			if ($fire->{psiout}{$k}[10] eq "HHS"){
				$corrector_HHS=$k-1;
				last;
			}
		}
		$numbering_filtered_residues=0;
		############### CHANGE in the CALCULATION OF THE COVERAGE OF THE SITE ##############
		## Aquí estamos calculando el promedio de la longitud de los sites que se han usado#
		## para obtener la prediccion de estos residuos 
		my $mean_site_length;
		foreach my $x (@biggest_length_ever){
			$mean_site_length=$mean_site_length+$x;
		}
		if ($coverage_tag ne "MET"){
			$mean_site_length=$mean_site_length/scalar@biggest_length_ever;
		}
		else {
			my @ord=sort{$b<=>$a}(@biggest_length_ever);
	#		$mean_site_length=$ord[0];
			################ CHANGE ##########################
			my %median;
			foreach my $sizes(@biggest_length_ever){$median{$sizes}++;}
			my @oreo=sort {$median{$b} <=> $median{$a}} keys %median;
			$mean_site_length=$oreo[0];
		}
		$gold_duo[1]=$mean_site_length;
		###################################################################################
		my @ordered_identit=sort { $b <=> $a} @score_identit;
		my $best_alignment_position=3000;
		foreach my $z (@position_array){
			if ($z>$corrector_HHS){$z=$z-$corrector_HHS;}
			if ($z<$best_alignment_position){$best_alignment_position=$z;}
		}
		if($numberOFsites/$hitNumOfSites > $relative_num_of_cases and $hit_mean >= $hit_mean_cutoff and $hit_score >= $hit_score_cutoff){
	#	aqui tenemos 3 filtros: el primero es una comparacion entre el tamano del cluster actual y el tamano del cluster mas grande (calculado antes). 
	#	Este tiene que ser mayor de un umbral pre-establecido (ahora 0).
	#	El segundo establece que la mejor media de los Scores de SQUARE de los sites del cluster tiene que superar un umbral pre-establecido (ahora 1).
	#	El tercero compara el mejor score de todos los sites con un umbral pre-establecido (ahora 3).
			$autoincr++;
			push(@{$this->{site_score}{$autoincr}},$ordered_identit[0],$gold_duo[0],$gold_duo[1]);
			$this->{summary}{numOfSites}[$autoincr]=$numberOFsites;
			$this->{summary}{score}[$autoincr]=$hit_mean;
			if(exists$this->{summary_order}{$numberOFsites}){$this->{summary_order}{$numberOFsites}.=" $autoincr";}
			else{$this->{summary_order}{$numberOFsites}=$autoincr;}
		#~~~~~~~~~	Here we are calcultaing the relative frequency of every residue in the prediction
		#~~~~~~~~~	The residues that appear in close homologues have an increased frequency (more weight)
			my %all_res_freq;
			my %well;
			my $res_hit_gold;
			if (scalar(keys %well60) != 0){
				%well=%well60;
				$res_hit_gold=$res_hit_gold60;
			}
			elsif(scalar(keys %well35) != 0){
				%well=%well35;
				$res_hit_gold=$res_hit_gold35;
			}
			foreach my $num(keys%res){
				if (exists $well{$num}){
					my $good_freq=$well{$num}/$res_hit_gold;
					my $norm_freq=$res{$num}/$res_hit;
					$all_res_freq{$num}=($good_freq+$norm_freq)/2;
				}
				else{$all_res_freq{$num}=$res{$num}/($res_hit+$res_hit_gold);} 
									### aqui normalizamos la frecuencia de aparicion de todos los residuos:
			}						### dividimos el numero de veces que aparece por el numero max de veces que 
			
		#~~~~~~~~~	this method select the possible ligand. The choice is based on the frequency 
		#~~~~~~~~~	but also taking in account 
			$this->compound_selection_site_tagging($lib,\%compids,$global_clus_pred_tag,$autoincr,$pulser);
			
		#~~~~~~~~~	 
			if ($this->{summary}{type}[$autoincr] eq "MET" or $this->{summary}{type}[$autoincr] eq "POSS_MET"){
				if ($res_cut_off==0.25){$res_cut_off=0.35;}
				else {$res_cut_off=0.25;}
			}
			my $query_pos=0;
			my @seq=split(//,$fire->{sequence});
			for(@seq){
				if($_ ne "-"){
					$query_pos++;
					if(exists $all_res_freq{$query_pos} and $all_res_freq{$query_pos} ne "none"){
						if($all_res_freq{$query_pos} > $res_cut_off){	
							$numbering_filtered_residues++;
							$FINAL_SQUARE_SCORE_COMBO+=$global_SQUARE_best_score{$query_pos};
							if(exists $this->{summary}{resnum}[$autoincr]){
								$this->{summary}{resnum}[$autoincr].=" $query_pos";
								$this->{summary}{resname}[$autoincr].=" $_";
								my $frequenza=sprintf "%.2f", $all_res_freq{$query_pos};
								$this->{summary}{resresume}[$autoincr].="$_($query_pos); ";
								$this->{summary}{resfreq}[$autoincr].="$_($query_pos)=$frequenza; ";
								$this->{summary}{SQUARE}[$autoincr].=" $query_pos($all_scorez{$query_pos})";
								push (@{$this->{pred_sites_list}{$autoincr}},$query_pos);
								$this->{all_bind_res}{$autoincr}.=" $query_pos";
							}
							else{	$this->{summary}{resnum}[$autoincr]=$query_pos;
								$this->{summary}{resname}[$autoincr]=$_;
								my $frequenza=sprintf "%.2f", $all_res_freq{$query_pos};
								$this->{summary}{resresume}[$autoincr]="$_($query_pos); ";
								$this->{summary}{resfreq}[$autoincr]="$_($query_pos)=$frequenza; ";
								$this->{summary}{SQUARE}[$autoincr]="$query_pos($all_scorez{$query_pos})";
								push (@{$this->{pred_sites_list}{$autoincr}},$query_pos);
								$this->{all_bind_res}{$autoincr}=$query_pos;
							}
						}
					}
				}
			}
		}
		else{next;}
		if ($numbering_filtered_residues==0){$FINAL_SQUARE_SCORE_COMBO=0;}
		else{
			$FINAL_SQUARE_SCORE_COMBO=$FINAL_SQUARE_SCORE_COMBO/$numbering_filtered_residues;	###     TRASH
		}												###     TRASH
		push(@{$this->{site_score}{$autoincr}},$numbering_filtered_residues);					###     TRASH
		${$this->{site_score}{$autoincr}}[1]=$FINAL_SQUARE_SCORE_COMBO;						###     TRASH
	}
}

sub collapse_overlapping_csa{
	my ($this,$fire)=@_;
	my %all_csa_res; my %all_csa_score; my %csa_res_freq;
	if (scalar(keys%{$fire->{all_csa_res_60}})>5){
		%all_csa_res=%{$fire->{all_csa_res_60}};
		%all_csa_score=%{$fire->{all_csa_score_60}};
		%csa_res_freq=%{$fire->{csa_res_freq_60}};
	}
	elsif(scalar(keys%{$fire->{all_csa_res_30}})>15){
		%all_csa_res=%{$fire->{all_csa_res_30}};
		%all_csa_score=%{$fire->{all_csa_score_30}};
		%csa_res_freq=%{$fire->{csa_res_freq_30}};
	}
	else {
		%all_csa_res=%{$fire->{all_csa_res}};
		%all_csa_score=%{$fire->{all_csa_score}};
		%csa_res_freq=%{$fire->{csa_res_freq}};
	}
	my %csa_res;
	my %csa_res_aa;
	foreach my $csiteid(keys%all_csa_res){          ## este es un bucle para clusterizar los CSA, ya que no se ha hecho antes con la subrutina
		my $key_res=$all_csa_res{$csiteid};     ## residue_clustering !!!!
		my @csa_res=split(/ /,$key_res);        ## en este array estan todas las posiciones de los aa en la sec query
		my @csa_res_ord=sort@csa_res;
		my $all_aa;
		my $new_string;
		foreach my$num(@csa_res_ord){
			if ($num ne '' && $new_string eq ''){$new_string="$num";}
			elsif($num ne ''){$new_string.=" $num";}
			$all_aa=$all_aa.$fire->{all_csa_aa}{$num};      ## en $all_aa guardamos todos los aa, codigo de una letra
		}
		if(exists$this->{csa_sites}{$new_string}){              ## generamos un hash con CLAVE: "23 45 75 .." y VALOR: el Cat. Site ID
			$this->{csa_sites}{$new_string}.=" $csiteid";
		}
		else{$this->{csa_sites}{$new_string}=$csiteid;}
	}       
	my $flag="RED";
	while ($flag eq "RED"){
		my $change="ON";
		foreach my $csa_block (keys%{$this->{csa_sites}}){
			my $csa_block2;
			my @key_res2;
			my @key_res=split(/ /,$csa_block);
			foreach $csa_block2 (keys%{$this->{csa_sites}}){
				my $iguales=0;
				if ($csa_block2 eq $csa_block or $csa_block2 eq ''){next;}
				@key_res2=split(/ /,$csa_block2);
				foreach my $iii(@key_res){
					foreach my $aaa(@key_res2){
						if ($iii==$aaa){$iguales++;}
					}
				}
				if ($iguales==scalar@key_res && defined $csa_block){
					$this->{csa_sites}{$csa_block2}.=" $this->{csa_sites}{$csa_block}";
					delete$this->{csa_sites}{$csa_block};
					$change="OFF";
					last;
				}
				if (scalar@key_res < scalar@key_res2 && $iguales/scalar@key_res >= 0.33){
					my %common;
					foreach my $iii(@key_res){$common{$iii}=5;}
					foreach my $iii(@key_res2){$common{$iii}=5;}
					my $new_block=join(' ',sort(keys%common));
					$this->{csa_sites}{$new_block}=$this->{csa_sites}{$csa_block2}." ".$this->{csa_sites}{$csa_block};
					delete$this->{csa_sites}{$csa_block};
					delete$this->{csa_sites}{$csa_block2};
					$change="OFF";
					last;
				}
			}
			if ($change eq "OFF"){last;}
		}
		if ($change eq "ON"){$flag="GREEN";}
	}
	return(\%all_csa_res,\%all_csa_score,\%csa_res_freq)
}


sub residue_clustering{
	my ($this,%all_csite_res)=@_;
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # cut-offs manually established
	my $first_cutoff=0.6;
	my $scnd_cutoff=0.6;
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	my %all_csite_res=%{$refe};
	my %finalClstr=();
	foreach my $id1(keys%all_csite_res){            ### a este hash se asigna un valor distinto cada llamada; CLAVE: CSITE_ID; VALOR: pos absol AA;
		$finalClstr{$id1}=$id1;                         ### %finalClstr es el hash final donde se guardan por cada clust-id los otros
		my @res1=split(/ /,$all_csite_res{$id1});       ### IDs de los clusters con los que comparten el 60% de los aa !!
		my $len1=scalar@res1;
		my @aparca=split(/_/,$id1);
		my $CSITE1=$aparca[0];
		my $ligand_type1=$this->{dbi}->ligtype(-csiteid=>$CSITE1);
		foreach my $id2(keys%all_csite_res){    ### aqui empieza una comparacion uno contra todos de los sites
			unless($id1 eq $id2){
				my @res2=split(/ /,$all_csite_res{$id2});
				my $len2=scalar@res2;
				my $common=0;
				foreach my $uno(@res1){
					foreach my$dos(@res2){
						if($uno eq $dos){ # si la posicion de un aminoacido del primer site es igual a la de otro se aumenta
							$common++;      ### esta variable.
							last;
						}
					}
				}
				@aparca=();
				@aparca=split(/_/,$id2);
				my $CSITE2=$aparca[0];
				my $ligand_type2=$this->{dbi}->ligtype(-csiteid=>$CSITE2);
				if ($ligand_type1 eq $ligand_type2 && $ligand_type1 eq "NOM"){
					if($common >= $len1 * $first_cutoff or $common >= $len2 * $scnd_cutoff){
													## si el numero de los aa en comun es superior
						if(exists$finalClstr{$id1}){                            ## al 60% de la longitud de ambos sites, se
							$finalClstr{$id1}="$finalClstr{$id1} $id2";     ## ponen juntos
						}
					else{   $finalClstr{$id1}=$id2;}
					}
				}
				else {
					if($common >= $len1 * $first_cutoff and $common >= $len2 * $scnd_cutoff){
													## si el numero de los aa en comun es superior
						if(exists$finalClstr{$id1}){                            ## al 60% de la longitud de ambos sites, se
							$finalClstr{$id1}="$finalClstr{$id1} $id2";     ## ponen juntos
						}
						else{   $finalClstr{$id1}=$id2;}
					}
				}
			}
		}
	}
	my @cbind=$this->id_clustering(values%finalClstr);
	return @cbind;
}


sub id_clustering{
	my ($this,@global)=@_;
	my $start;
	my $end;
	do{
		my %uniq_id=();
		$start=scalar@global;
		for(@global){
			my @spl=split(/ /,$_);
			my $join=join(" ",@spl);
			for my$id(@spl){
				if(exists$uniq_id{$id}){$uniq_id{$id}="$uniq_id{$id} $join";}
				else{$uniq_id{$id}=$join;}
			}
		}
		@global=();
		for(values%uniq_id){
			my @spl = split(/ /,$_);
			my @sort = uniqsort(@spl);
			my $sort = join(" ",@sort);
			push(@global,$sort);
		}
		@global=uniqsort(@global);
		$end=scalar@global;
	}until ($start==$end);
	return @global;
}

sub uniqsort{
	my %uniq=();
	foreach my$key(@_){
		chomp $key;
		$uniq{$key}="";
	}
	my @out=sort{$a cmp $b}(keys(%uniq));
	return @out;
}


sub weighter{
	#########################################################################################################################
	# Aquí damos un peso diferente a los residuos, que depende de la frecuencia relativa de aparicion en los alineamientos	#
	# y del tipo de alineamiento. Si unos residuos resultan aparecer en alineamientos con e-value bueno y % de identidad	#
	# alto (> 60 o > 35) tendran un peso mayor en la prediccion respecto a los demas.     					#
	#########################################################################################################################
	my($fire,$cut_off,$alig,$len_ali,$residues,$hit_60,$hit_35,$well60,$well35,$score_identit)=@_;
	my $probo;
	if ($fire->{psiout}{$alig}[10] eq "HHS" && $fire->{psiout}{$alig}[8]=~/Ide.+Probab=(.+)/){$probo=$1;}
	if ($fire->{psiout}{$alig}[10] eq "PSI" && ($fire->{psiout}{$alig}[9] eq "0.0" || $fire->{psiout}{$alig}[9]=~/e\-(\d+)/)){
		my $evalor;
		my $identit=0;
		if ($fire->{psiout}{$alig}[9]=~/e\-(\d+)/){$evalor=sprintf("%.10f", $fire->{psiout}{$alig}[9]);}
		if ($fire->{psiout}{$alig}[8]=~/\d+\/\d+\s+\((\d+)%\)/){$identit=$1;}
		if ($fire->{psiout}{$alig}[9]!=0 && ($evalor>0.005 or $len_ali<=$fire->{longest_ali})){$identit=20;}
		push (@{$score_identit},$identit);
		if ($identit>59 && $len_ali>= $fire->{longest_ali}){
			$$cut_off=0.45;
			foreach my $res(@{$residues}){
				if(exists $$well60{$res}){
					$$well60{$res}++;
					if($$well60{$res}>$$hit_60){$$hit_60=$$well60{$res};}
				}
				else{$$well60{$res}=1;}
			}
		}
		elsif (scalar(keys%{$well60})==0 && $identit>=30 && $len_ali>= $fire->{longest_ali}){
			$$cut_off=0.45;
			foreach my $res(@{$residues}){
				if(exists $$well35{$res}){
					$$well35{$res}++;
					if($$well35{$res}>$$hit_35){$$hit_35=$$well35{$res};}
				}
				else{$$well35{$res}=1;}
			}
		}
	}
	elsif ($fire->{psiout}{$alig}[10] eq "PSI" && $fire->{psiout}{$alig}[8]=~/\d+\/\d+\s+\((\d+)%\)/){
		my $identit=$1;
		if ($fire->{psiout}{$alig}[9]>0.005){$identit=20;}
		push (@{$score_identit},$identit);
	}
	elsif ($fire->{psiout}{$alig}[10] eq "HHS" && $probo>=95){
		my $identit=0;
		if ($fire->{psiout}{$alig}[8]=~/Identities=(\d+)%,/){$identit=$1;}
		if ($len_ali<= $fire->{longest_ali}){$identit=20;}
		push (@{$score_identit},$identit);
		if ($identit>59 && $len_ali>= $fire->{longest_ali}){
			$$cut_off=0.45;
			foreach my $res(@{$residues}){
				if(exists $$well60{$res}){
					$$well60{$res}++;
					if($$well60{$res}>$$hit_60){$$hit_60=$$well60{$res};}
				}
				else{$$well60{$res}=1;}
			}
		}
		elsif (scalar(keys%{$well60})==0 && $identit>=30 && $len_ali>= $fire->{longest_ali}){
			$$cut_off=0.45;
			foreach my $res(@{$residues}){
				if(exists $$well35{$res}){
					$$well35{$res}++;
					if($$well35{$res}>$$hit_35){$$hit_35=$$well35{$res};}
				}
				else{$$well35{$res}=1;}
			}
		}
	}
	elsif ($fire->{psiout}{$alig}[10] eq "HHS" && $fire->{psiout}{$alig}[8]=~/Identities=(\d+)%,/){
		my $identit=$1;
		if ($probo<95){push (@{$score_identit},20);}
		else {push (@{$score_identit},$identit);}
	}
	return 0;
}


sub compound_selection_site_tagging{
	my ($this,$lib,$compids,$TAG,$autoinc,$pulser)=@_;
	my %comp_freq;
	foreach my$compid(keys%{$compids}){
		my $freq=$compids->{$compid};
		if(exists$comp_freq{$freq}){$comp_freq{$freq}.=" $compid";}
		else{$comp_freq{$freq}=$compid;}
	}
	#	en este hash %comp_freq la CLAVE es el numero de veces que aparece un compuesto en contacto en los sites y el VALOR es el nombre de los
	#	compuestos con esa frecuencia
	my @all_comp;
	my @freqs=sort{$b<=>$a}(keys%comp_freq);
	foreach my$freq(@freqs){
		my @spl=split(/ /,$comp_freq{$freq});
		for(@spl){
			push(@all_comp,"$_($freq)");	
		}
	}
	my $all_comp=join(" ",@all_comp);
	#	en la variable $all_comp se guardan todos los compuestos del cluster ordenados de mayor a menor frecuencia en el formato nombre(frec)
	#################	SELECTION OF THE COMPOUND BASED ON THE "GENERAL TYPE" OF THE SITE	###############
	my @compid;
	my $first_poss;
	my $flag_poss_cog="DOWN";
	my $flag_cog="DOWN";
	for (my $xy=0;$xy<scalar@all_comp;$xy++){
		@compid=();
		@compid=split(/\(/,$all_comp[$xy]);     ## aqui elegimos el compuesto que esta haciendo binding (el que aparezca
		$compid[1]=~s/\)//;                     ## mas veces haciendo binding en los CSITES
		if (scalar@compid==0){next;}
		if (exists $lib->{metals}{$compid[0]}){last;}
		if ($TAG eq "YES" && exists $lib->{cognate}{$compid[0]}){
			$flag_cog="UP";
			last;
		}
		elsif ($TAG eq "YES" && exists $lib->{poss_cognate}{$compid[0]} && $flag_poss_cog eq "DOWN"){
			$first_poss=$compid[0];
			$flag_poss_cog="UP";
		}
		elsif ($TAG eq "NON" && !exists $lib->{cognate}{$compid[0]}){last;}
	}
#	if ($global_clus_pred_tag eq "YES" && $flag_cog eq "DOWN"){$compid[0]=$first_poss;}
	if ($TAG eq "YES" && $flag_cog eq "DOWN" && $compid[0] ne "FEO"){$compid[0]=$first_poss;}		## delete this line after FireDB update
	my $compname=$this->{dbi}->compound_name(-compid=>$compid[0]);
	$this->{summary}{compid_code}[$autoinc]=$compid[0];
	if(exists $lib->{metals}{$compid[0]} && exists $lib->{poss_cognate}{$compid[0]}){
		$this->{summary}{nice_try}[$autoinc]="POSSIBLE_COGNATE";
		$this->{summary}{type}[$autoinc]="MET_POSS";
	}
	elsif (exists $lib->{metals}{$compid[0]} && $pulser eq "STRICT"){
		$this->{summary}{nice_try}[$autoinc]="POSSIBLE_COGNATE";
		$this->{summary}{type}[$autoinc]="MET_POSS";
	}
	elsif (exists $lib->{poss_cognate}{$compid[0]}){
		$this->{summary}{nice_try}[$autoinc]="POSSIBLE_COGNATE";
		$this->{summary}{type}[$autoinc]="NOM";
	}
	elsif (exists $lib->{metals}{$compid[0]}){
		$this->{summary}{nice_try}[$autoinc]="COGNATE";
		$this->{summary}{type}[$autoinc]="MET";
	}
	elsif (exists $lib->{cognate}{$compid[0]} and $pulser eq "STRICT"){
		$this->{summary}{nice_try}[$autoinc]="POSSIBLE_COGNATE";
		$this->{summary}{type}[$autoinc]="NOM";
	}
	elsif (exists $lib->{cognate}{$compid[0]}){
		$this->{summary}{nice_try}[$autoinc]="COGNATE";
		$this->{summary}{type}[$autoinc]="NOM_COG";
	}
	else {
		$this->{summary}{nice_try}[$autoinc]="NON_COGNATE";
		$this->{summary}{type}[$autoinc]="NOM";
	}
	if ($compid[0] eq "FEO"){
		$this->{summary}{nice_try}[$autoinc]="NON_COGNATE";
		$this->{summary}{type}[$autoinc]="NOM";
	}
	$this->{summary}{compid}[$autoinc]=$all_comp;
	$this->{summary}{compname}[$autoinc]=$compname;
}

sub scoring_and_final_list{
	my ($this)=@_;
	my @order=sort{$b<=>$a}(keys%{$this->{summary_order}});
	my $joinSites=join(" ",values%{$this->{summary_order}});
	my @joinSites=split(/ /,$joinSites);
	$this->{scalar_sites}=scalar@joinSites;
###################################### NEW SITE' SCORE #############################################

## Hemos guardado informacion sobre cada cluster of sites en el hash %site_score con clave el/los ID(s) del clustered site.
## En concreto se trata del mejor Score de SQUARE entre todos los sites clusterizados y su respectivo coverage, mas la
## mejor posicion en el hash de alineamientos psiout (no necesariamente esta informacion tiene que llegar del mismo alineamiento
## de donde se han sacado los otros dos parametros. Se normalizan (a parte del coverage) y se hace una media aritmetica, luego convertida en 
## porcentaje. Este score teoricamente va de 0 a 3, pero si un site llega hasta aqui, no puede llegar a ser 0.
	foreach my $god_hand(@order){
		my @spl=split(/ /,$this->{summary_order}{$god_hand});
		foreach my $z (@spl){
			my $new_coverage;
			if ($this->{site_score}{$z}[3]>$this->{site_score}{$z}[2]){$new_coverage=1;}
			else {$new_coverage=$this->{site_score}{$z}[3]/$this->{site_score}{$z}[2];}
			my $R_Ali=rali($this->{summary}{numOfSites}[$z],$order[0]);
			#my $park=$new_coverage+2*((${$site_score{$z}}[1])/6)+((${$site_score{$z}}[0])/100)+$R_Ali;
			$this->{score_list}{$z}=sprintf("%.3f",($new_coverage+2*(($this->{site_score}{$z}[1])/6)+(($this->{site_score}{$z}[0])/100)+$R_Ali));
			#$score_list{$z}=sprintf("%.3f",$park);
			my $output=sprintf("%.1f",($this->{score_list}{$z}/5)*100);
			$this->{summary}{reliability}[$z]=$output;                                      # aqui guardamos el reliability score para el filtro !!
			$output.="% [COV: ".ceil($new_coverage*100)."% SITE: ";
			$output.=floor(($this->{site_score}{$z}[1])/6*100)."% Ident: ".$this->{site_score}{$z}[0]."% - Ali: ";
			$output.=ceil(sprintf("%.1f",($R_Ali*100)))."%]";
			$this->{summary}{score}[$z]=$output;
		}
	}
######################## ELIMINATING REDUNDANT BINDING SITE #########################################
## este es un bucle para eliminar clustered sites que solapan perfectamente con otro site mas grande
## esto puede pasar si no se han clusterizado al principio y luego se han quitado algunos residuos por
## los filtros. No se eliminan tal cual los sites, se han puesto reglas para dar la prioridad a los
## COGNATE !! (No deberia verificarse esta situacion, pero en el caso de que haya misma longitud entre
## un site y otro, nos quedamos con el site con el score mas alto)
	foreach my $id1(keys%{$this->{all_bind_res}}){
		my @res1=split(/ /,$this->{all_bind_res}{$id1});
		my $len1=scalar@res1;
		my $type1=$this->{summary}{type}[$id1];
		my $score1= $this->{summary}{score}[$id1];
		foreach my$id2(keys%{$this->{all_bind_res}}){
			unless($id1 eq $id2){
				my @res2=split(/ /,$this->{all_bind_res}{$id2});
				my $len2=scalar@res2;
				my $type2=$this->{summary}{type}[$id2];
				my $score2= $this->{summary}{score}[$id2];
				my $common=0;
				foreach my$uno(@res1){
					foreach my$dos(@res2){
						if($uno eq $dos){$common++;}
					}
				}
				if ($common == $len1 && $common == $len2){
					if ($type1 eq $type2 && $score1 > $score2){$this->{shinigami}{$id2}=4;}
					elsif ($type1 eq $type2 && $score1 < $score2){$this->{shinigami}{$id1}=4;}
				}
				else {
					if($common == $len1 && $type1 eq "NOM"){$this->{shinigami}{$id1}=4;}
					elsif($common == $len2 && $type2 eq "NOM"){$this->{shinigami}{$id2}=4;}
					elsif($common == $len1 && $type1 eq "MET" && $type1 eq $type2){$this->{shinigami}{$id1}=4;}
					elsif($common == $len2 && $type2 eq "MET" && $type1 eq $type2){$this->{shinigami}{$id2}=4;}
					elsif($common == $len1 && $type1 eq "MET_POSS"){$this->{shinigami}{$id1}=4;}
					elsif($common == $len2 && $type2 eq "MET_POSS"){$this->{shinigami}{$id2}=4;}
					elsif($common == $len1 && $type1 eq "NOM_COG" && $type1 eq $type2){$this->{shinigami}{$id1}=4;}
					elsif($common == $len2 && $type2 eq "NOM_COG" && $type1 eq $type2){$this->{shinigami}{$id2}=4;}
				}
			}
		}
	}
}


sub rali{
	my $number=$_[0];
	my $MAX=$_[1];
	my $b;
	if ($MAX>100){
		my $multiplier=100/$MAX;
		$b=($multiplier*$number)/100;
		#     if ($number>100){$b=1;}
		#     else{$b=$number/100;}
	}
	elsif ($MAX>10){
		$b=$number/$MAX;
		$b=sprintf("%.2f",$b);
	}
	else {
		$b=$number/10;
	}
	return($b);
}





sub DESTROY {
        my ($this)=@_;
        $this->{dbi}->close();
}

1;
