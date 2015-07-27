
package results_parser;
use strict;


sub new {
	my($class,%args)=@_;
	my $this={};
	bless($this);
	return $this;
}

sub parse {
	my ($this,$target)=@_;
	open(F,"$this->{dir_target}/$target.res");
	my $flag_CAT="DOWN";
	my $number;
	while(<F>){
		if($_=~/^\n/){$number=undef;next;}
		if($_=~/firestar did not detect ligand binding sites for your query/){last;}
		chomp($_);
		if ($_=~/^CAT.+/){
			my @line=split(/\t+/,$_);
			$number=$line[1];
			if($line[2]=~/EC_number/){
				unless ($line[3] eq "CAT" or $line[3]=~/no EC nor GO information/){
					$flag_CAT="UP";
					$this->{$target}{CAT}{$number}{ec_number}{$line[3]};
				}
			}
			elsif($line[2]=~/Evidence/){$flag_CAT="DOWN";$this->{$target}{CAT}{$number}{evidence}=$line[3];}
			elsif($line[2]=~/Homologs/){
				$this->{$target}{CAT}{$number}{homologs}=$line[3];
			}
			elsif($line[2]=~/Residue_positions/){$this->{$target}{CAT}{$number}{res_posi}=$line[3];}
			elsif($line[2]=~/Residue_composition/){$this->{$target}{CAT}{$number}{residues}=$line[3];}
			elsif($line[2]=~/Site_score/){$this->{$target}{CAT}{$number}{score}=$line[3];}
		}
		if ($_=~/^\t+.+/ and $flag_CAT eq "UP" and defined $number){
			$_=s/\t//g;
			$this->{$target}{CAT}{$number}{ec_number}.="\t$_";
		}
		if ($_=~/^SITE.+/){
			my @line=split(/\t+/,$_);
			$number=$line[1];
			if($line[2]=~/Compound_name/){
				if ($line[3]=~/UNKNOWN/){
					$this->{$target}{SITE}{$number}{comp_name}="unknown";
				}
			}
			elsif($line[2]=~/Compound_type/){$this->{$target}{SITE}{$number}{tag}=$line[3];}
			elsif($line[2]=~/Compound_id/){$this->{$target}{SITE}{$number}{all_comp_ids}=$line[3];}
			elsif($line[2]=~/Compound_GO/){
				if ($line[3] ne '' and $this->{$target}{SITE}{$number}{comp_name} ne 'unknown'){
					$this->{$target}{SITE}{$number}{GO_terms}=$line[3];
				}
			}
			elsif($line[2]=~/GO_name/){
				if ($line[3] ne '' and $this->{$target}{SITE}{$number}{comp_name} ne 'unknown'){
					$this->{$target}{SITE}{$number}{GO_name}=$line[3];
				}
			}
			elsif($line[2]=~/Residue_positions/){$this->{$target}{SITE}{$number}{res_posi}=$line[3];}
			elsif($line[2]=~/Residue_composition/){$this->{$target}{SITE}{$number}{residues}=$line[3];}
			elsif($line[2]=~/Probab\. score/){$this->{$target}{SITE}{$number}{res_freqs}=$line[3];}
			elsif($line[2]=~/Reliability/){
				if ($line[3]=~/(.+)% \[COV: (.+)% SITE: (.+)% Ident: (.+)% - Ali: (.+)%\]/){
					$this->{$target}{SITE}{$number}{reli}=$1;
					$this->{$target}{SITE}{$number}{cov}=$2;
					$this->{$target}{SITE}{$number}{square}=$3;
					$this->{$target}{SITE}{$number}{ident}=$4;
					$this->{$target}{SITE}{$number}{alignments}=$5;
				}
			}
			elsif($line[2]=~/Number_of_homologs/){$this->{$target}{SITE}{$number}{homologs}=$line[3];}
		}
					
	}
	close(F);
	return 0;
}



1;
