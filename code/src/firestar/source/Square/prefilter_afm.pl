#!/usr/bin/perl -w

use DBI;
use FindBin;
my $cwd=$FindBin::Bin;
use Config::IniFiles;
my $variables=Config::IniFiles->new(-file => "$cwd/../CONFIG_fire_var.ini");

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Prueba de conexion con base de datos
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my $dbHost=$variables->val('MYSQL','dbHost');
my $dbName=$variables->val('MYSQL','dbName');
my $user=$variables->val('MYSQL','user');
my $pass=$variables->val('MYSQL','pass');
my $sth;                #objeto DBI
my $dataHandle = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost","$user","$pass",{RaiseError => 1,AutoCommit =>0})
                                || die "Unable to connect to $dbName because $DBI::errstr";
my $serverInfo = $dataHandle->{'mysql_serverinfo'};
my $serverStat = $dataHandle->{'mysql_stat'};
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $ssearch=$variables->val('PROGRAMS','ssearch');

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $queryAln=$ARGV[0];	# Query aligned
my $temptAln=$ARGV[1];	# Template aligned
my $template=$ARGV[2];	# Template code

my $clustid;
my @identity;
my $temptSeq;
my %sites;

$sth = $dataHandle -> prepare("select CLUSTID,CONSEQ from CONSENSUS where CADID=\"$template\"");
$sth -> execute();
my @tmp = $sth -> fetchrow_array();
$clustid=$tmp[0];
$temptSeq=$tmp[1];
$sth -> finish();

my $template_raw_seq=$temptAln;	$template_raw_seq=~s/-//g;
open(USER,">user_seq_tmp_file.tmp");
print USER ">user\n$template_raw_seq";
close USER;

open(DBSEQ,">db_seq_tmp_file.tmp");
print DBSEQ ">dbseq\n$temptSeq";
close DBSEQ;
my @ssearch =`$ssearch -d1 user_seq_tmp_file.tmp db_seq_tmp_file.tmp -a -h -m=3 -Q`;

`rm user_seq_tmp_file.tmp db_seq_tmp_file.tmp`;
@identity=grep{$_=~ /identity/}@ssearch;
@ssearch=@ssearch[25..$#ssearch-7];
@ssearch=grep{$_=~ /^dbseq/ or $_=~/^user/}@ssearch;

#print @identity;
#print @ssearch;		# stores ssearch output
my $line;
my $dbseq;
my $userseq;
while(scalar@ssearch>0){
	chomp $ssearch[0];
	if($ssearch[0]=~/^dbseq/){
		$line=substr($ssearch[0],7,60);
		$dbseq=$dbseq.$line;
		$line=~s/./-/g;
		$userseq=$userseq.$line;
		shift@ssearch;
		}
	elsif($ssearch[0]=~/^user/){
		$line=substr($ssearch[0],7,60);
		$line=~s/ /-/g;
		$userseq=$userseq.$line;
		shift@ssearch;
		if($ssearch[0]=~/^dbseq/){
			$line=substr($ssearch[0],7,60);
			$line=~s/ /-/g;
			$dbseq=$dbseq.$line;
			shift@ssearch;
			}
		elsif($ssearch[0]=~/^user/){
			$line=~s/./-/g;
			$dbseq=$dbseq.$line;
			}
		}
	};

# All printed stuff is here
chomp $dbseq;
print "$queryAln\n$temptAln\n";
print "$userseq\n$dbseq\n$clustid\n";
print @identity;

$dataHandle -> disconnect();

END;
