=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::RunCentipede

=head1 DESCRIPTION

=cut


package Bio::EnsEMBL::Funcgen::RunnableDB::RunCentipede;

use warnings;
use strict;
use Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor;
use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;

sub fetch_input {   # nothing to fetch... just the parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();
  
  my $work_dir = $self->param('work_dir');
  my $matrix = $self->param('matrix');
  my $dnase = $self->param('dnase');

  my $file = $work_dir."/output/".$matrix."_".$dnase.".counts";
  throw "expected file $file does not exist" if(! -e $file);
  
  return 1;
}

sub run {   # Create Groups and Analysis and stuff...
  my $self = shift @_;

  my $hlen = 100;
  my $work_dir = $self->param('work_dir');
  my $matrix = $self->param('matrix');
  my $dnase = $self->param('dnase');
  my $is_male = $self->param('is_male');

  my $file = $work_dir."/output/".$matrix."_".$dnase.".counts";
  open(FO,">".$file.".r");
  print FO "library(CENTIPEDE);\n";
  #print FO "pdf(file='".$file.".pdf')\n";
  print FO 'data<-read.table("'.$file.'");'."\n";
  print FO 'fit <- fitCentipede(Xlist = list(DNase=as.matrix(data[,7:dim(data)[2]])), Y=cbind(rep(1, dim(data)[1]), data[,5]));'."\n";
  print FO 'write.table(data[which(fit$PostPr>0.99),1:6],file="'.$file.'.sites",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);'."\n";
  #print FO 'plotProfile(fit$LambdaParList[[1]])'."\n";
  #print FO "dev.off();\n";
  close FO;

  system("/software/bin/R-2.11.1 CMD BATCH --slave ${file}.r ${file}.Rout");

  system("gzip $file");
  #system("rm -f $file");

  return 1;
}


sub write_output {  # Nothing to do here
  return 1;
}

1;
