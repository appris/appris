package fire_CSA_study;

use strict;
use FindBin;
my $cwd;
BEGIN{
	$cwd=$FindBin::Bin;
}

use lib $cwd;
use fire_libraries;
use DBI_firestar;
use Config::IniFiles;
use File::Copy;
use Getopt::Long qw(GetOptionsFromArray);

sub new {
	my($class,%args)=@_;
	my $this={};
	bless($this);
	return $this;
}

sub input_parameters_eva{
	my ($this,$entry)=@_;
	my ($queryname)=undef;
	my ($sequence)=undef;
	my ($evalue)=10;
	my ($outfile)=undef;
	my ($cutoff)=0;
	my ($csa_option)='YES';
	my ($cognate_option)='YES';
	my ($output)='text';
	my ($infile)=undef;
	my ($runtype)=undef;
	my ($config)="$this->{cwd}/../CONFIG_fire_var.ini";
	my ($chain)="A";
	my ($pdbcode)=undef;
	my ($tmpfile)=undef;
	my ($mail)=undef;
	my ($unicode)=undef;

	&GetOptionsFromArray(
		\@{$entry},
		'q=s'		=> \$queryname,
		's=s'		=> \$sequence,
		'e=f'		=> \$evalue,
		'o=s'		=> \$outfile,
		'cut=i'		=> \$cutoff,
		'csa=s'		=> \$csa_option,
		'cog=s'		=> \$cognate_option,
		'opt=s'		=> \$output,
		'f=s'		=> \$infile,
		'type=s'	=> \$runtype,
		'chain=s'	=> \$chain,
		'conf=s'	=> \$config,
		'pdb=s'		=> \$pdbcode,
		'm=s'		=> \$mail,
		'tmpfile=s'	=> \$tmpfile,
		'uni=s'		=> \$unicode
	);
	
	
	unless ((defined $queryname and defined $sequence) or defined $infile or defined $pdbcode or defined $unicode){
		print "\n\n\n\n Some required input parameters are missing\n\n\n\n\n";
		return "FAIL";
	}

	if (defined $queryname and defined $sequence and defined $infile){
		print "\n\n\n\n input file and input sequence can't be introduced at the same time\n\n\n\n\n";
		return "FAIL";
	}

	if (defined $queryname and defined $sequence and defined $pdbcode){
		print "\n\n\n\n input sequence and PDB code can't be introduced at the same time\n\n\n\n\n";
		return "FAIL";
	}
	
	if (defined $queryname and defined $sequence and defined $unicode){
		print "\n\n\n\n Uniprot accession code and input sequence can't be introduced at the same time\n\n\n\n\n";
		return "FAIL";
	}

	if (defined $infile and defined $pdbcode){
		print "\n\n\n\n input file and PDB code can't be introduced at the same time\n\n\n\n\n";
		return "FAIL";
	}

	if (defined $infile and defined $unicode){
		print "\n\n\n\n input file and Uniprot accession code can't be introduced at the same time\n\n\n\n\n";
		return "FAIL";
	}
	
	if (defined $pdbcode and defined $unicode){
		print "\n\n\n\n Uniprot accession code and PDB code can't be introduced at the same time\n\n\n\n\n";
		return "FAIL";
	}
	
	if (defined $infile){
		unless (defined $runtype and ($runtype=~/str/i or $runtype=~/fas/i)){
			print "\n\n\n\n please specify if the file you're introducing is a pdb (-type=str) or a fasta (-type=fas) type\n\n\n\n\n";
			return "FAIL";
		}
	}
	
	#~~~~~~~~~~~	CONNECTION 2 FireDB	~~~~~~~~~~~~~~~
	
	$this->{config_file}=$config;
	$this->{dbi}=DBI_firestar->connect($this->{config_file});
	
	#~~~~~~~~~~     IMPORTANT PATHS     ~~~~~~~~~~~~~~~~~

	my $variables=Config::IniFiles->new(-file => $this->{config_file});

	$this->{home}=$variables->val('PATHS','home');
	$this->{faatmp}=$variables->val('PATHS','faatmp');

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if (defined $infile and $runtype eq "fas"){	$this->fasta_reader($infile);}
	elsif (defined $infile and $runtype eq "str"){
		my $answer=$this->pdb_reader($infile,$chain);
		if ($answer eq "NO"){
			print "\n\n\n\n the PDB file you introduced doesn't seem to be a valid one. Please check it or check the chain you selected\n\n\n\n\n";
			return "FAIL";
		}
		
	}
	elsif (defined $pdbcode){
		my $answer=$this->pdbcode_retriever($pdbcode,$chain);
		$runtype='pdb';
		if ($answer eq "NO"){
			print "\n\n\n\n firestar can't retrieve any protein sequence from FireDB using the code you introduced. Please check if the release date\n";
			print "is later than the current release date of FireDB (",$variables->val('DATABASES','release'),") or if your code is wrong\n\n\n\n\n"; 
			return "FAIL";
		}
	}
	
	elsif (defined $queryname and defined $sequence){
		$this->{queryname}=$queryname;
		$this->{sequence}=$sequence;
	}

	elsif (defined $unicode){
		my $answer=$this->uniprot_retriever($unicode,$variables->val('EXT_LINKS','DB_fetch'));
		if ($answer eq "NO"){
			print "\n\n\n\n the Uniprot accession code you introduced doesn't seem to be a valid one. Please check it\n\n\n\n\n";
			return "FAIL";
		}
	}
	
	else {
		print "\n\n\n\n something is wrong with your input parameters; please check the available documentation and run again firestar\n\n\n\n\n";
		return "FAIL";
	}

	$this->{evalue}=$evalue;
	$this->{cutoff}=$cutoff;
	$this->{csa_option}=$csa_option;
	$this->{cognate_option}=$cognate_option;
	$this->{output}=$output;
	$this->{infile}=$infile;
	$this->{chain}=$chain;
	$this->{runtype}=$runtype;
	$this->{pdbcode}=$pdbcode;
	$this->{tmpfile}=$tmpfile;
	$this->{mail}=$mail;
	if ($output eq "web"){$this->{outfile}=$this->{tmpfile};}
	else{$this->{outfile}=$outfile;}
	if ($this->{output} eq 'web'){
		open (PHP,"$this->{faatmp}/$this->{tmpfile}.php");			## aqui estamos guardando el fichero de espera $tmpfile.php para poder imprimir poco a poco
		while(<PHP>){								## que es lo que esta haciendo en ese momento firestar por detras
			if ($_!~/class="statuslog"/){push(@{$this->{header_php}},$_);}
			else {push(@{$this->{header_php}},$_);last;}
		}
		close(PHP);
	}
	return "OK";
}

sub homologs_searching{
	my ($this)=@_;
	if ($this->{output} eq 'appris'){$this->{tmpfile}=$this->{queryname};}
	unless (defined $this->{tmpfile}){$this->{tmpfile}=$this->RANDNAME;}
	$this->{full}="$this->{faatmp}/$this->{tmpfile}";
	open(TMP,'>',"$this->{full}.faa") or die "Error $this->{full} no se abre";
	print TMP ">Query\n$this->{sequence}";close TMP;
	my $program_choice='BOTH';
	if ($this->{output} ne "appris"){
		open (PSIDONE,"$this->{faatmp}/FAA_LOG.txt");
		while (<PSIDONE>){
			chomp($_);
			my @check=split(/\t/,$_);
			my $cache=$check[0];
			if ($this->{sequence} eq $check[1] && -e "$this->{faatmp}/$cache.hhr"){
				my $parking=`ls $this->{faatmp}/$cache\_*.psi`;chomp($parking);
				my @split=split(/_/,$parking);
				my $old_evalue=$split[$#split];$old_evalue=~s/\.psi$//;
				if ($old_evalue == $this->{evalue}){
					copy "$this->{faatmp}/$cache.hhr","$this->{full}.hhr";
					copy "$parking","$this->{full}\_$this->{evalue}.psi";
					if ($this->{output} eq 'web'){
						$this->user_info_printer("  Same sequence's analysis found .... information recovered  ");
						$this->user_info_printer("gif");
						copy "$this->{faatmp}/$cache.chk","$this->{full}.chk";
					}
					return 0;
				}
				else{
					copy "$this->{faatmp}/$cache.hhr","$this->{full}.hhr";
					if ($this->{output} eq 'web'){
						$this->user_info_printer("  Same sequence's analysis found with different parameters, HHsearch output recovered  ");
						$this->user_info_printer("gif");
					}
					$program_choice='PSI';
				}
			}
			elsif ($this->{sequence} eq $check[1] && length($this->{sequence})>2000 && -e -s "$this->{faatmp}/$cache\_$this->{evalue}.psi"){
				copy "$this->{faatmp}/$cache\_$this->{evalue}.psi","$this->{full}\_$this->{evalue}.psi";
				if ($this->{output} eq 'web'){
					$this->user_info_printer("  Same sequence's analysis found .... information recovered  ");
					$this->user_info_printer("gif");
					copy "$this->{faatmp}/$cache.chk","$this->{full}.chk";
				}
				return 0;
                        }
                }
                close(PSIDONE);
        }
        else{   #esta parte se ha modificado para APPRIS !!! corre HHblits solo si la prot tiene menos de 251 aa y corre en AHSOKA !!
		open (PSIDONE,"$this->{faatmp}/FAA_LOG.txt");
		while (<PSIDONE>){
			chomp($_);
			my @check=split(/\t/,$_);
			my $already_run=$check[0];
			if ($this->{sequence} eq $check[1] && length($this->{sequence})<=250 && -e "$this->{faatmp}/$already_run.hhr" && -e "$this->{faatmp}/$already_run.psi"){
				$this->{tmpfile}=$already_run;  #### OJO !!! esta parte es superimportante porque cambiamos el nombre temporal !!!!!!!
				$this->{full}="$this->{faatmp}/$this->{tmpfile}";
                        	return 0;
			}
			elsif($this->{sequence} eq $check[1] && length($this->{sequence})>250 && -e "$this->{faatmp}/$already_run.psi"){
				$this->{tmpfile}=$already_run;  #### OJO !!! esta parte es superimportante porque cambiamos el nombre temporal !!!!!!!
				$this->{full}="$this->{faatmp}/$this->{tmpfile}";
                        	return 0;
			}
		}
		
		close(PSIDONE);
	}
	open (PSIDONE,'>>',"$this->{faatmp}/FAA_LOG.txt");
	print PSIDONE "$this->{tmpfile}\t$this->{sequence}\n";
	close(PSIDONE);
	unless ($program_choice eq 'PSI'){
		if (length($this->{sequence})>2000 or (length($this->{sequence})>250 and $this->{output} eq "appris")){	$program_choice='PSI';}
		else{	$program_choice='BOTH';}
	}
	if ($this->{output} eq 'web' and $program_choice eq 'PSI'){
		$this->user_info_printer("  Starting PSI-BLAST analysis only. Sequence is too long (",length($this->{sequence})," a.a.) for HHsearch analysis ...  ");
		$this->user_info_printer("wait_for");
	}
	elsif ($this->{output} eq 'web'){
		$this->user_info_printer("  Starting PSI-BLAST and HHsearch analysis  ");
		$this->user_info_printer("wait_for");
	}
	$this->fork_analysis($program_choice);
	if ($this->{output} eq 'web'){	$this->user_info_printer("gif");}
	return 0;
}


sub fork_analysis{
	my ($this,$choice)=@_;
	if ($choice eq 'BOTH'){
		my $control="FALSE";
		while ($control eq "FALSE"){
			my $pid1 = fork();
			if($pid1){
				$this->HHsearcher();
				open (R,'>',"$this->{full}.hhr_tmpa");
				close(R);
				waitpid($pid1,0);
			}
			elsif($pid1 == 0){
				$this->psiBlaster();
				open (T,'>',"$this->{full}.psi_tmpa");
				close(T);
				exit 0;
			}
			if (-e "$this->{full}.hhr"){
				if (-e -s "$this->{full}\_$this->{evalue}.psi"){$control="TRUE";}
                        }
			unlink "$this->{full}.psi_tmpa";
			unlink "$this->{full}.hhr_tmpa";
			if (-e "$this->{full}.chk"){
				unlink ("$this->{full}.hhm");
				unlink ("$this->{full}.a3m");
				unless ($this->{output} eq 'web'){unlink ("$this->{full}.chk");}
			}
			return 0;
		}
	}
	else{
		$this->psiBlaster();
		if (-e "$this->{full}.chk"){
			unlink ("$this->{full}.chk");
		}
		return 0;
	}
}


sub psiBlaster{
	my ($this)=@_;
	my $variables=Config::IniFiles->new(-file => $this->{config_file});
	my $nrdb=$variables->val('DATABASES','nrdb');
	my $release_date=$variables->val('DATABASES','release');
	my $path2db=$variables->val('DATABASES','blast_path');
	my $blast_bin_path=$variables->val('PROGRAMS','bstbin');
	my $name=$this->{full};
	if ($this->{output} ne 'appris'){	$name="$this->{full}\_$this->{evalue}";}
	`$blast_bin_path/blastpgp -a4 -C $this->{full}.chk -d $path2db/$nrdb -e0.01 -F F -h0.01 -j3 -b0 -v50 -i $this->{full}.faa 2> /dev/null`;
	`$blast_bin_path/blastpgp -a2 -R $this->{full}.chk -d $path2db/fdbTptDB_$release_date -F F -e$this->{evalue} -v0 -i $this->{full}.faa -o $name.psi 2> /dev/null`;
	return 0;
}


sub HHsearcher{
	my ($this)=@_;
	my $variables=Config::IniFiles->new(-file => $this->{config_file});
	my $release_date=$variables->val('DATABASES','release');
	my $path2db=$variables->val('DATABASES','hhdb_path');
	my $hhbdb=$variables->val('DATABASES','hhbdb');
	my $hh_nr_db=$variables->val('DATABASES','hhprof');
	`hhblits -cpu 4 -i $this->{full}.faa -d $path2db/$hh_nr_db -oa3m $this->{full}.a3m -o /dev/null`;
	`perl $ENV{"HHLIB"}/scripts/addss.pl -i $this->{full}.a3m`;
	`hhmake -i $this->{full}.a3m`;
	`hhsearch -cpu 4 -d $path2db/$hhbdb$release_date\_hhm_db -i $this->{full}.hhm`;
	return 0;
}


sub parse_psiblast{
	my ($this)=@_;
	if ($this->{output} eq 'web'){	$this->user_info_printer("  Extracting information from PSI-BLAST search  ");}
	if ($this->{output} ne "appris"){
		open(PSI,"$this->{full}\_$this->{evalue}.psi") or print "no pude abrir $this->{full}\_$this->{evalue}.psi";
	}
	else{open(PSI,"$this->{full}.psi") or print "no pude abrir $this->{full}.psi";}
	my @psifile=<PSI>;close PSI;
	my @evalue=grep{$_=~/^ Score/}@psifile;
	my @psiout=grep{$_=~/^>/ or $_=~/^ Identities/ or $_=~/^Query:/ or $_=~/^Sbjct:/}@psifile;
	for(@psiout){
		$_=~s/^>/>\n/;
	}
	my $allpsi=join("",@psiout);
	@psiout=split(/>/,$allpsi);
	my $keyId=0;
	my %psiout;
	for(my$i=0;$i<100;$i++){
		my $ln=$psiout[$i];
		$ln=~s/^\n//;   $ln=~s/\n$//;
		my @ev=split(/,/,$evalue[$i-1]);
		my $expect=$ev[1];      $expect=~s/\s+//g;      $expect=~s/Expect=//;
		my @ALN=split(/\n/,$ln);
		my $template=shift@ALN;
		if($#ALN >1){
			my @ID=();
			for(my $i=0;$i<$#ALN+1;$i++){
				if($ALN[$i]=~/Identities/){push(@ID,$i);}
			}
			my @SINGLE=();
			my $queryStart;
			my $tempStart;
			my $queryEnd;
			my $tempEnd;
			for(my $j=0;$j<$#ID+1;$j++){
	                        if(exists$ID[$j+1]){@SINGLE=@ALN[$ID[$j]..$ID[$j+1]-1];}
				else{@SINGLE=@ALN[$ID[$j]..$#ALN];}
				my $identities=shift@SINGLE;
				my @temp=grep{$_=~/^Sbjct/}@SINGLE;
				my @query=grep{$_=~/^Query/}@SINGLE;
				my @SEtemp=split(/\s+/,$temp[0]);
				my @SEquery=split(/\s+/,$query[0]);
				$queryStart=$SEquery[1];                # Residuo de comienzo queryet
				$tempStart=$SEtemp[1];                  # Residuo de comienzo template
				@SEquery=split(/\s+/,$query[$#query]);
				@SEtemp=split(/\s+/,$temp[$#temp]);
				$queryEnd=$SEquery[3];                  # Residuo final queryet
				$tempEnd=$SEtemp[3];                    # Residuo final template
				my $tempSeq="";
				my $querySeq="";
				for(@query){
					my @spl=split(/\s+/,$_);
					$querySeq=$querySeq.$spl[2];
				}
				for(@temp){
					my @spl=split(/\s+/,$_);
					$tempSeq=$tempSeq.$spl[2];
				}
				$keyId++;
				$psiout{$keyId}[0]=$this->{tmpfile};
				$psiout{$keyId}[1]=$querySeq;
				$psiout{$keyId}[2]=$template;
				$psiout{$keyId}[3]=$tempSeq;
				$psiout{$keyId}[4]=$queryStart;
				$psiout{$keyId}[5]=$queryEnd;
				$psiout{$keyId}[6]=$tempStart;
				$psiout{$keyId}[7]=$tempEnd;
				$psiout{$keyId}[8]=$identities;
				$psiout{$keyId}[9]=$expect;
				$psiout{$keyId}[10]='PSI';
			}
		}
	}
	%{$this->{psiout}}=%psiout;
	if ($this->{output} eq 'web'){	$this->user_info_printer("gif");}
	return 0;
}


sub parse_hhsearch{
	my ($this)=@_;
	if ($this->{output} eq 'web'){	$this->user_info_printer("  Extracting information from HHsearch result  ");}
	my %psiout=%{$this->{psiout}};
	my $querySeq;
	my $template;
	my $tempSeq;
	my $identities;
	my $expect;
	my $flag1="OFF";
	my @coords_query;
	my @coords_templ;
	my $keyId=scalar(keys%psiout)+1;
	open(HHR,"$this->{full}.hhr") or print "no pude abrir $this->{full}.hhr\n\n";
	while (defined(my $line=<HHR>)){
		if ($line=~/>.+/){
			if (defined $template){
				$psiout{$keyId}[0]=$this->{queryname};
				$psiout{$keyId}[1]=$querySeq;
				$psiout{$keyId}[2]=$template;
				$psiout{$keyId}[3]=$tempSeq;
				$psiout{$keyId}[4]=$coords_query[0];
				$psiout{$keyId}[5]=$coords_query[$#coords_query];
				$psiout{$keyId}[6]=$coords_templ[0];
				$psiout{$keyId}[7]=$coords_templ[$#coords_templ];
				$psiout{$keyId}[8]=$identities;
				$psiout{$keyId}[9]=$expect;
				$psiout{$keyId}[10]='HHS';
				$keyId++;
				$flag1="OFF";
				$querySeq="";$tempSeq="";@coords_templ=();@coords_query=();
				$template=$line;$template=~s/>//;chomp($template);
			}
			else {$template=$line;$template=~s/>//;chomp($template);}
		}
		elsif($line=~/^Probab/){
			if ($flag1 eq "ON"){
				$psiout{$keyId}[0]=$this->{queryname};          #V
				$psiout{$keyId}[1]=$querySeq;
				$psiout{$keyId}[2]=$template;
				$psiout{$keyId}[3]=$tempSeq;
				$psiout{$keyId}[4]=$coords_query[0];
				$psiout{$keyId}[5]=$coords_query[$#coords_query];
				$psiout{$keyId}[6]=$coords_templ[0];
				$psiout{$keyId}[7]=$coords_templ[$#coords_templ];
				$psiout{$keyId}[8]=$identities;         #V
				$psiout{$keyId}[9]=$expect;             #V
				$psiout{$keyId}[10]='HHS';
				$keyId++;
			}
			$flag1="ON";
			my @parking=split(/\s+/,$line);
			$expect=$parking[1];$expect=~s/E-value=//;
			my $Probab=$parking[0];
			$Probab=~s/Probab=//;
			if ($Probab<5){close(HHR);last;}
			$identities=join(", ",$parking[4],$parking[0]);
		}
		elsif($line=~/^Q $this->{queryname}/ || $line=~/^Q Query/ || ($line=~/^Q/ && $line!~/ss_pred|ss_conf|Consensus|Query/)){
			my @parking=split(/\s+/,$line);
			push(@coords_query,$parking[2],$parking[4]);
			$querySeq.=$parking[3];
		}
		elsif(defined $template && $line=~/^T $template/){
			my @parking=split(/\s+/,$line);
			push(@coords_templ,$parking[2],$parking[4]);
			$tempSeq.=$parking[3];
		}
		else{next;}
	}
	close(HHR);
	%{$this->{psiout}}=%psiout;
	if ($this->{output} eq 'web'){	$this->user_info_printer("gif");}
	return 0;
}

sub site_info_extractor{
	my($this,$lib,$pulser)=@_;
	my %psiout=%{$this->{psiout}};
	$this->{longest_ali}=0;
	$this->{already_seen_ECs}=();
	my $variables=Config::IniFiles->new(-file => $this->{config_file});
	my $afmpath=$variables->val('PROGRAMS','square');
	my @spit;
	push(@spit,"PASS-THROUGH\tPROGRAM\tTEMPLATE\tIDENTITIES (%)\tEXPECTED VALUE\tMEAN SQUARE SCORE\tCSA SOURCE\tRESIDUES SCORE SUM\tRESIDUES SCORE MEAN\tRESIDUES(SQUARE SCORE)");
	foreach my $key (keys%psiout){
		my %parameters;
		$parameters{query}		=$psiout{$key}[0];
		$parameters{query_aligned}	=$psiout{$key}[1];
		$parameters{template}		=$psiout{$key}[2];		# nombre del template encontrado
		$parameters{templ_aligned}	=$psiout{$key}[3];
		$parameters{q_start}		=$psiout{$key}[4];		# coordenada absoluta del primer aminoacido de la query en este alineamiento
		$parameters{q_end}		=$psiout{$key}[5];		# coordenada absoluta del ultimo aminoacido de la query en este alineamiento
		$parameters{t_start}		=$psiout{$key}[6];		# coordenada absoluta del primer aminoacido del TEMPLATE en este alineamiento
		$parameters{t_end}		=$psiout{$key}[7];		# coordenada absoluta del ultimo aminoacido del TEMPLATE en este alineamiento
		$parameters{len_ali}		=$psiout{$key}[5]-$psiout{$key}[4]+1;
		$parameters{identities}		=$psiout{$key}[8];		# identities del alineamiento (y probabilidad en el caso de HHsearch)
		$parameters{expect}		=$psiout{$key}[9];		# score of the alignments extracted from the program output
		$parameters{program}		=$psiout{$key}[10];		# program from where the alignment has been extracted
		if ($parameters{identities}=~/\d+\/\d+\s+\((\d+)%\)/ && $psiout{$key}[10] eq "PSI"){$parameters{identit}=$1;}
		elsif ($parameters{identities}=~/Identities=(\d+)%,/ && $psiout{$key}[10] eq "HHS"){$parameters{identit}=$1;}
	#	if ($parameters{identit}>45){next;}			# line used for firestar evaluation
		my @afmout=`$afmpath $psiout{$key}[1] $psiout{$key}[3] $parameters{template}`;		# aqui se lanza el programa another_fire_mess_web; 
		for(@afmout){chomp $_;$_=~s/\s+//g;}
		$parameters{afm_query}=$afmout[1];
		$parameters{afm_tempt}=$afmout[0];
		$parameters{afm_score}=$afmout[2];
		if ($this->{output} eq 'web'){
			%{$this->{parameters}{$key}}=%parameters;
			$this->extended_info($key);
		}
		$parameters{afm_score}=~s/@/6/g;$parameters{afm_score}=~s/-/0/g;
		my @afm_score=split(//,$parameters{afm_score});	# este array contiene los score de SQUARE por cada posicion 	e.g.	221----1-------12343
		my $afm_mean=0;
		for(my $i=0;$i<$#afm_score+1;$i++){
			if($afm_score[$i] eq "-"){$afm_score[$i]=0;}		# e.g.	221----1 => 22100001
			elsif($afm_score[$i] eq "@"){$afm_score[$i]=6;}		# e.g.	@34----2 => 63400002
			$afm_mean=$afm_mean+$afm_score[$i];	## al final del bucle en $afm_mean no tenemos una media sino una suma de todos los scores.
		}
		#	aqui estamos evaluando el alineamiento: si hay por lo menos un a.a. con un score de SQUARE de 2, entra
		if($afm_mean >1){
			$parameters{afm_mean} = $afm_mean/length($parameters{afm_score}); # aqui si que hacemos la media aritmetica
			my $connection_check=$this->{dbi}->ping();
			unless($connection_check == 1){$this->{dbi}=DBI_firestar->connect($this->{config_file});}
			my %binding_sites=$this->{dbi}->collapsed_sites(-clustid=>$parameters{template});
			foreach my $i(keys%binding_sites){
				my @csite=@{$binding_sites{$i}};
				if($csite[5] ne 'CSA' or $csite[3] ne 'lit'){next;}
				push(@csite,"${$binding_sites{$i}}[0]\_$key");
				my $line=$this->CSA_filtering(\%parameters,\@csite,$lib);
				if ($line ne ""){push(@spit,$line);}
			}
		}
	}
	return(@spit);
}

sub afm_parser {
	my ($pm)=@_;
	my @afm_query=split(//,$pm->{afm_query});     # este array contiene los a.a. de la query, alineados           e.g.    KLIEQAKKWGHPAIAVTDHA
	my @afm_tempt=split(//,$pm->{afm_tempt});     # este array contiene los a.a. del template, alineados          e.g.    EMVLKAIELDFDEYSIVEHA
	my @afm_score=split(//,$pm->{afm_score});
	my %afm_temp_score;
	my %afm_temp_aa;
	my %afm_targ_aa;
	my %afm_targ_num;
	my $t_count=0;
	my $q_count=0;
	for(my $i=0;$i<$#afm_tempt+1;$i++){ # en este bucle genera 4 hashes en los cuales se guardan informacion relativa:
		if($afm_tempt[$i] =~ /[A-Z]/){
			my $t_pos=$pm->{t_start}+$t_count;
			$t_count++;
			if($afm_query[$i] =~ /[A-Z]/){
				my $q_pos=$pm->{q_start}+$q_count;
				$q_count++;
				$afm_targ_num{$t_pos}=$q_pos;		# CLAVE: pos. absoluta templ VALOR: pos. absoluta query
				$afm_temp_aa{$t_pos}=$afm_tempt[$i];	# CLAVE: pos. absoluta templ VALOR: a.a. correspondiente template
				$afm_targ_aa{$t_pos}=$afm_query[$i];	# CLAVE: pos. absoluta templ VALOR: a.a. correspondiente query
				$afm_temp_score{$t_pos}=$afm_score[$i];	# CLAVE: pos. absoluta templ VALOR: SQUARE score de la posicion
			}
		}
		elsif($afm_tempt[$i] !~ /[A-Z]/ and $afm_query[$i] =~ /[A-Z]/){
			$q_count++;
		}
	}
	return(\%afm_targ_num,\%afm_temp_aa,\%afm_targ_aa,\%afm_temp_score);
}


sub CSA_filtering{
	my ($this,$para,$site,$lib)=@_;
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# cut-offs manually established
	my $loose_cutoff=3;
	my $tight_cutoff=4;
	my $cutPercentage_csa=65.99;
	my $mean_score_cutoff=4;
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	my @numconres=split(/ /,$site->[1]);		# absolute position in the consensus sequence of the binding residues
	my @conres=split(//,$site->[2]);		# one letter code aminoacid in the consensus sequence of the binding residues
	my $csiteid=$site->[6];
	my $ali_cov=$para->{len_ali}/length $this->{sequence};
	my ($afm_targ_num,$afm_temp_aa,$afm_targ_aa,$afm_temp_score)=afm_parser($para);
	my @source=$this->{dbi}->CSA_LIT(-csite=>$site->[0]);
	my $line;
	if (scalar@source ==0){return($line);}
	my $mean_score;
	my @error;
	my %good_res;
	my %bad_res;
	for(my $i=0;$i<$#numconres+1;$i++){
                my $t_pos=$numconres[$i];
                if(exists$afm_targ_num->{$t_pos} and $afm_temp_score->{$t_pos} >= $tight_cutoff){
			$good_res{$t_pos}=5;
                }
                elsif(exists $afm_targ_num->{$t_pos} && $afm_temp_score->{$t_pos} >= $loose_cutoff && $afm_temp_aa->{$t_pos} eq $afm_targ_aa->{$t_pos}){
			$good_res{$t_pos}=5;
                }
                if(exists $afm_targ_num->{$t_pos} and $afm_temp_score->{$t_pos} < $loose_cutoff){
			$bad_res{$t_pos}=5
		}
                $mean_score=$mean_score+$afm_temp_score->{$t_pos};
        }
        if(scalar(keys%good_res)>0){$mean_score=$mean_score/($#numconres+1);}
        else{$mean_score=0};
        my $pass="NO";
        if (length$this->{sequence} <= 120 && $ali_cov >= 0.5){$pass="YES";}
        elsif (length$this->{sequence} > 120 && $para->{len_ali}>70){$pass="YES";}
	my $mean_SQ_score=sprintf("%.2f", $para->{afm_mean});
#       if(($#good_res+1)==($#numconres+1) and $mean_score>4)
        if((scalar(keys%good_res)*100/($#numconres+1)) > $cutPercentage_csa && $mean_score>$mean_score_cutoff){
		$line="YES\t$para->{program}\t$para->{template}\t$para->{identit}\t$para->{expect}\t$mean_SQ_score\t$source[0]\t";
		my $parking;
		my $sum;
		for(my $i=0;$i<$#numconres+1;$i++){
			my $t_pos=$numconres[$i];
			my $resi=$conres[$i];
			if(exists$afm_temp_score->{$t_pos}){	$parking.="$resi($afm_temp_score->{$t_pos}) ";$sum+=$afm_temp_score->{$t_pos};}
			else {	$parking.="$resi(NA) ";}
		}
		$line.="$sum\t".sprintf("%.2f",$sum/scalar@numconres)."\t$parking";
	}
	elsif (scalar(keys%good_res)>0){
		$line="NO\t$para->{program}\t$para->{template}\t$para->{identit}\t$para->{expect}\t$mean_SQ_score\t$source[0]\t";
		my $parking;
		my $sum;
		for(my $i=0;$i<$#numconres+1;$i++){
			my $t_pos=$numconres[$i];
			my $resi=$conres[$i];
			if(exists$afm_temp_score->{$t_pos}){	$parking.="$resi($afm_temp_score->{$t_pos}) ";$sum+=$afm_temp_score->{$t_pos};}
			else {	$parking.="$resi(NA) ";}
		}
		$line.="$sum\t".sprintf("%.2f",$sum/scalar@numconres)."\t$parking";
	}
	elsif (scalar(keys%bad_res)>0){
		$line="NO\t$para->{program}\t$para->{template}\t$para->{identit}\t$para->{expect}\t$mean_SQ_score\t$source[0]\t";
		my $parking;
		my $sum;
		for(my $i=0;$i<$#numconres+1;$i++){
			my $t_pos=$numconres[$i];
			my $resi=$conres[$i];
			if(exists$afm_temp_score->{$t_pos}){	$parking.="$resi($afm_temp_score->{$t_pos}) ";$sum+=$afm_temp_score->{$t_pos};}
			else {	$parking.="$resi(NA) ";}
		}
		$line.="$sum\t".sprintf("%.2f",$sum/scalar@numconres)."\t$parking";
	}
	return($line);
}

sub RANDNAME{
	my ($this)=@_;
	my $flag="CLOSE";
	my $chunk;
	while ($flag eq "CLOSE"){
		my $randnum=rand(1);
		$randnum=~s/^0\.//;
		$chunk=substr ($randnum,0,12);
		unless (-e "$this->{faatmp}/$chunk.faa"){$flag="OPEN";}
	}
	return($chunk);
}

sub fasta_reader{
	my ($this,$file)=@_;
	open(F,$file);
	my $seq;
	while(<F>){
		if ($_=~/^\n/){next;}
		chomp($_);
		if($_=~/^>.+/){
			$_=~s/\|/ /g;
			$_=~s/>//g;
			my @parking=split(/\s+/,$_);
			if (length$parking[0]<18){$this->{queryname}=$parking[0];}
			else{$this->{queryname}=substr($parking[0],0,16);}
		}
		else{	$seq.=$_;}
	}
	close(F);
	$this->{sequence}=$seq;
	return 0;
}

sub pdb_reader{
	my ($this,$file,$chain)=@_;
	my @pdbseq;
	if (defined $chain){
		$chain=uc$chain;
		@pdbseq=`$this->{cwd}/getpdbseq.pl -f $file -s -n -x -c $chain`;
	}
	else {  @pdbseq=`$this->{cwd}/getpdbseq.pl -f $file -s -n -x`;}
	if(scalar@pdbseq==0){return "NO";}
	chomp $pdbseq[0];
	$pdbseq[0]=~s/>//;$pdbseq[0]=~s/_//;
	$this->{queryname}=$pdbseq[0];
	chomp $pdbseq[1];
	$this->{sequence}=$pdbseq[1];
	return "YES";
}

sub pdbcode_retriever{
	my ($this,$code,$chain)=@_;
	$code=lc$code;
	$chain=uc$chain;
	my @info=$this->{dbi}->pdb2sequence(-cadid=>"$code$chain");
	if (scalar@info==0){return "NO";}
	$this->{queryname}="$code$chain";
	$this->{sequence}=$info[0];
	$this->{sequence}=~s/\-//g;
	$this->{cluster}=$info[1];
	return "YES";
}

sub uniprot_retriever{
	my ($this,$code,$link)=@_;
	$link.=$code;
	`wget -O /tmp/$code.tmp '$link' 2> /dev/null`;
	open (F,"/tmp/$code.tmp");
	my @fetch=<F>;
	close(F);
	my @error=grep(/No result found/,@fetch);
	if (scalar@error > 0){	return "NO";}
	my $flag="OFF";
	my $seq;
	foreach (@fetch){
		if ($_=~/^SQ\s+SEQUENCE/){	$flag="ON";next;}
		elsif ($_=~/^\/\/\n/){	$flag="OFF";last;}
		if ($flag eq "OFF"){	next;}
		else{
			chomp($_);
			$_=~s/\s//g;
			$seq.=$_;
		}
	}
	$this->{queryname}=$code;
	$this->{sequence}=$seq;
	return "YES";
}

sub user_info_printer {
	my ($this,$string)=@_;
        if ($string =~/wait_for$/){
                open(TEMPO,'>',"$this->{full}.update");
                foreach (@{$this->{header_php}}){print TEMPO $_;}
                print TEMPO "<img width=\"14\" src=\"../../images/loading2.gif\" </img>";
        }
        elsif ($string =~/gif$/){
                push (@{$this->{header_php}},"<img width=\"12\" src=\"../../images/chopped.png\" </img> \n");
                open(TEMPO,'>',"$this->{full}.update");
                foreach (@{$this->{header_php}}){print TEMPO $_;}
        }
        else { 
                push (@{$this->{header_php}},$string);
                open(TEMPO,'>',"$this->{full}.update");
                foreach (@{$this->{header_php}}){print TEMPO $_;}
        }
        print TEMPO "\n\t\t</div>\n\t</div>\n</body>\n</html>";;
        close(TEMPO);
        rename "$this->{full}.update","$this->{full}.php";
}





sub DESTROY {
	my ($this)=@_;
	$this->{dbi}->close();
}





1;
