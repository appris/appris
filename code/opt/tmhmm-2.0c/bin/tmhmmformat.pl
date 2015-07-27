#!/usr/local/bin/perl -w

# This is version 2.0c of tmhmmformat.pl

# Take output from decodeanhmm -N1 -PrintNumbers -PrintScore -PrintStat [-plp]
# Make GFF output, save files, etc

# Flush immediately
$| = 1;

$gnuplot = "/usr/local/bin/gnuplot-3.7 ";
$epsmove = "/sbin/grep -v \'translate\$\' ";
$ps2gif = "/usr/local/bin/gs -dNOPAUSE -q -sDEVICE=ppm -g710x350 -r72 -sOutputFile=- - -c quit | /usr/freeware/bin/ppmtogif - ";


# OPTION PARSING ##########################################
use Getopt::Long;

$opt_html = 0;       # Make HTML output
$opt_short = 0;      # Make short output format
$opt_plot = 1;       # Make plots
$opt_v1 = 0;         # Use old model (version 1)
$opt_signal = 10;    # Cut-off used for producing a warning of possible
                     # Signal sequence
$opt_Nterm = 60;     # Number of bases to consider for signal peptides
                     # in the Nterm


$opt_workdir = ".";  # Working dir.
$opt_wwwdir = ".";   # The place where the www server looks for files
                     # (The www name for the working dir)
$opt_serverhome = "."; # This is the location of various stationary files
                       # (Help files etc)

# Process options
$result = GetOptions ('workdir=s','wwwdir=s','html!','short!','plot!','v1!','signal=i','serverhome=s','Nterm=i');
die ("Error on command line") unless $result;

###########################################

%gffwords = ( "i", "inside", "M", "TMhelix", "o", "outside", "O", "outside");

$opt_workdir .= "/" if ($opt_workdir =~ /[^\/]$/);
$opt_wwwdir .= "/" if ($opt_wwwdir =~ /[^\/]$/);
$opt_serverhome =~ s/\/$//;

if ($opt_v1) { $version = "TMHMM1.0"; }
else {         $version = "TMHMM2.0"; }


# HTML header
if ($opt_html) {
    print <<END;
<HTML>
<HEAD>
<TITLE>TMHMM result</TITLE>
</HEAD>
<BODY BGCOLOR="#FFFFFF">
<h2>TMHMM result</h2>
<A HREF="$opt_serverhome/TMHMM2.0.guide.html#output">HELP</A>
with output formats
<P>
END
    print "<hr><pre>\n" if ($opt_short );
}


# Read input
#print STDERR "tmhmmformat: Reading input\n";


while (<>) {
    # Read id
    if ( /^>/ ) {
	chop;
	$id = $_;
	$id =~ s/^>\s*//;
	$id =~ s/\s.*$//;
	$fn = id2filename($_);
	$plp = $plpread = 0;
    }
    # Read plp output if present
    elsif (/^\#/) {
	if ( !$plp && $opt_plot ) {
	    open(PLPFILE,">$opt_workdir$fn.plp");
	    print PLPFILE "$_";
	}
	if ( $plp ) { # second ^#
	    @keys = split(/\s+/,$_);
	    shift @keys;
	    $N = $#keys+1;
	    $i = 0;
	    print PLPFILE "# AA\tinside\tmembr\toutside\n" if ($opt_plot);
	    for $k (@keys) {
		$col{$k} = $i++;
	    }
	    $Mstart = 0.;
	    $aaexp = 0.;
	    $pos = 0;
	}
	$plp=1;
    }
    elsif ( /^%/ ) {
	# End of plp output
	if ($plp) {
	    $plpread = 1;
	    $plp = 0;
	    close(PLPFILE) if ($opt_plot);
	}
	# Read other constants
	if (/^\%len\s+(\S+)/) { $len = $1; }
	elsif (/^%score (\S+) ([0-9\.]+) \(([0-9\.]+)/) {
	    $score{$1}=$2;
	    $normscore{$1}=$3;
	}
	# LAST line
	elsif (/^%pred NB\(0\): /) {
	    @pred = split(/,\s+/,$');
	    &print_stuff;
	}
    }
    # Is it plp output still?
    elsif ($plp) {
	++$pos;
	@x = split(/\s+/,$_);
	$lett = shift @x;

	# Check if labels are present in plp output
	$lab = shift @x if ( $#x==$N );

	# Now @x contains the probabilities
	$ins = $x[$col{"i"}];
	$mem = $x[$col{"M"}];
	$out = $x[$col{"o"}]+$x[$col{"O"}];

	# Stat
	$aaexp += $mem;
	$Mstart += $mem if ( $pos<=$opt_Nterm );
	$istart = $ins if ($pos==1);

	# Print for plot:
	print PLPFILE "$pos $lett\t$ins\t$mem\t$out\n" if ($opt_plot);
    }
}

# HTML footer
if ($opt_html) {
    print "</pre><P>" if ($opt_short );
    print "<hr>\n</BODY>\n</HTML>\n"
}




sub print_stuff {

    # Count number of predicted helices by N-best
    $ntmh=0;
    for $p (@pred) { ++$ntmh if ($p=~/^M/); }

    if ($opt_short) {
	printf "%s\tlen=%d",$id,$len;
	if ($plpread) {
	    printf "\tExpAA=%.2f",$aaexp;
	    printf "\tFirst%d=%.2f",$opt_Nterm,$Mstart;
	}
	printf "\tPredHel=%d",$ntmh;
	print "\tTopology=";
	for $p (@pred) {
	    ( $type, $b, $e ) = split(/\s+/,$p);
	    $type="o" if ($type eq "O");
	    if ($type eq "M") { print "$b-${e}"; }
	    else { print $type; }
	}
	print "\n"
    }
    else {
	print "<hr>\n<pre>\n" if ($opt_html);
	print "# $id Length: $len\n";
	print "# $id Number of predicted TMHs:  $ntmh\n";
	if ($plpread) {
	    print "# $id Exp number of AAs in TMHs: $aaexp\n";
	    print "# $id Exp number, first $opt_Nterm AAs:  $Mstart\n";
	    print "# $id Total prob of N-in:        $istart\n";
	    if ( $Mstart>$opt_signal ) {
		print "# $id POSSIBLE N-term signal sequence\n";
	    }
	}
	for $p (@pred) {
	    ( $type, $b, $e ) = split(/\s+/,$p);
	    printf "$id\t$version\t$gffwords{$type}\t%6d%6d\n",$b,$e;
	}
	print "</pre><P>\n" if ($opt_html);

	# Make a plot
	if ( $opt_plot && $plpread ) {
	    make_plot($opt_workdir.$fn,$id,$len,\@pred);
	    # Make the eps plot
	    system("$gnuplot $opt_workdir$fn.gnuplot");
	    if ($opt_html) {
		# Make gif file
		system("$epsmove $opt_workdir$fn.eps | $ps2gif > $opt_workdir$fn.gif");
		$file = "$opt_wwwdir$fn";
		print "<IMG SRC=\"$file.gif\"><P>\n";
		print "# <A HREF=\"$file.eps\">plot</A> in postscript, \n";
		print "<A HREF=\"$file.gnuplot\">script</A> for making the plot in gnuplot, \n";
		print "<A HREF=\"$file.plp\">data</A> for plot<br>\n";
	    }
	}
    }
}


# Take a sequence ID and return a safe file name
sub id2filename {
    $fn = $_[0];
    $fn =~ s/^>\s*//;
    $fn =~ s/\s.*\n*$//;
    $fn =~ s/[^0-9a-zA-Z_%-\.\~]/_/g;
    $fn = "x" if ( length($fn)==0 );
    $newfn = $fn;

    # Check if filename is already used
    if ( defined($filenames{$fn}) ) {
	$newfn .= $filenames{$fn};
	$filenames{$fn} += 1;
    }
    else {
	$filenames{$fn} = 1;
    }
    return $newfn;
}





##############################
# Make plot                  #
##############################

sub make_plot {
    $filename = $_[0];
    $id = $_[1];
    $len = $_[2];
    @features = @{$_[3]};

    $y{"i"} = 1.07;
    $y{"M"} = 1.09;
    $y{"o"} = 1.11;

    $lt{"i"} = "lt 3";
    $lt{"M"} = "lt 1";
    $lt{"o"} = "lt 4";

    $lw{"i"} = "lw 10";
    $lw{"M"} = "lw 40";
    $lw{"o"} = "lw 10";


    # Read N-best topology
    if ( @features ) {
	# Make 'arrows' for the plot
	if ( $#features>=0 ) {
	    $arrows="";
	    foreach $f (@features) {
		( $type, $b, $e ) = split(/\s+/,$f);
		$type = "o" if ( $type eq "O" );
		$yy = $y{$type};
		$ll = $lt{$type};
		$ww = $lw{$type};
		$arrows .= "set arrow from $b,$yy to $e,$yy nohead $ll $ww\n";
	    }
	}
    }

    open(PLOTFILE,">$filename.gnuplot");
    print PLOTFILE $arrows if (defined($arrows));

    print PLOTFILE <<END;
set key below
set title \"TMHMM posterior probabilities for $id\"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel \"position\"
set ylabel \"probability\"
set xrange [1:$len]
# Make the ps plot
set term postscript eps color solid \"Helvetica\" 30
set output \"$filename.eps\"
plot \"$filename.plp\" using 1:4 title \"transmembrane\" with impulses $lt{"M"} lw 2, \\
\"\" using 1:3 title \"inside\" with line $lt{"i"} lw 2, \\
\"\" using 1:5 title \"outside\" with line $lt{"o"} lw 2
exit
END
    close(PLOTFILE);
}
