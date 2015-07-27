package fire_libraries;

use strict;
use FindBin;
my $cwd;
BEGIN{
	$cwd=$FindBin::Bin;
}

use lib $cwd;
use DBI_firestar;
use FindBin;
use Config::IniFiles;

sub new {
	my($class,%args)=@_;
	my $this={};
	bless($this);
	my $path=$this->lib($args{-file});
	$this->{path}=$path;
	$this->{dbi}=DBI_firestar->connect($args{-file});
	return $this;
}

sub ec2go{
#EC:1.1.1.102 > GO:3-dehydrosphinganine reductase activity ; GO:0047560
	my ($this)=@_;
	open(EC2GO,"$this->{path}/ec2go");
        do{
		my $line=readline(EC2GO);
		if($line=~/^EC:(.+)\s>\sGO:.+; (GO\:\d+)/){
			push(@{$this->{ec2go}{$1}},$2);
			push(@{$this->{ec2go_full}{$1}},$line);
		}
	}until eof;
	close(EC2GO);
	#################### TEMPORAL CORRECTION OF CSA ANNOTATION
	open(CSA,"$this->{path}/good_lit_CSA");
	do{
		my $line=readline(CSA);
		chomp($line);
		$this->{good_CSA}{$line}=5;
		}until eof;
	close CSA;
	####################
	return 0;
}

sub cif2go{
	my ($this)=@_;
	my %list;
	open(CIF2GO,"$this->{path}/cif2go");
	do{
		my $line=readline(CIF2GO);
		if($line!~/^#/){
			chomp $line;
			my @spl=split(/\t/,$line);
			my $cif=$spl[0];
			$list{$cif}{go1}=$spl[1];
			$list{$cif}{go2}=$spl[2];
			$list{$cif}{name}=$spl[3];
		}
	}until eof;
	close(CIF2GO);
	%{$this->{cif2go}}=%list;
	return 0;
}

sub tag {
	my ($this,%args)=@_;
	my %compos;
	my %list=$this->{dbi}->bio_tag(-tag=>$args{-type});
	foreach my $i(keys%list){
		$compos{${$list{$i}}[0]}=5;
	}
	if ($args{-type} eq "COGNATE"){
		%{$this->{cognate}}=%compos;
	}
	else{%{$this->{poss_cognate}}=%compos;}
	return 0;
}

sub metal{
	my ($this,%args)=@_;
	my %compos;
	my %list=$this->{dbi}->met_tag(-tag=>$args{-type});
	foreach my $i(keys%list){
		$compos{${$list{$i}}[0]}=5;
	}
	%{$this->{metals}}=%compos;
	return 0;	
}

sub lib{
	my ($this,$config_file)=@_;
	my $cwd=$FindBin::Bin;
	my $variables=Config::IniFiles->new(-file => $config_file);
	my $home=$variables->val('PATHS','home');
	return("$home/lib");
}

sub csa2ec{
	my ($this)=@_;
	#12as    6.3.1.1
	open(CSA2EC,"$this->{path}/csa2ec");
	do{
		my $line=readline(CSA2EC);
		chomp $line;
		my @spl=split(/\t/,$line);
		unless ($spl[1] eq "-.-.-.-"){
			$spl[1]=~s/\.\-//g;
			$this->{csa2ec}{$spl[0]}=$spl[1];
		}
	}until eof;
	close(CSA2EC);
	return 0;
}

sub pdb2id{
	my ($this,$fire)=@_;
	my $connection_check=$this->{dbi}->ping();
	unless($connection_check == 1){$this->{dbi}=DBI_firestar->connect($fire->{config_file});}
	%{$this->{id2pdbcode}}=$this->{dbi}->id2pdb();
	return 0;
}


sub DESTROY {
	my ($this)=@_;
	$this->{dbi}->close();
}




1;
