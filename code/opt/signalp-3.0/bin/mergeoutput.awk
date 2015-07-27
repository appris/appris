#! /usr/bin/nawk -f

BEGIN {
	OFS="\n";

	if (www=="") { # WWW header provided by WebFace
		printf "*********************** SignalP 3.0 predictions "
		print "***********************"
	}	

	methodstr["nn"] = "neural networks (NN)"
	methodstr["hmm"] = "hidden Markov models (HMM)"
	methodstr["nn+hmm"] = methodstr["nn"] " and " methodstr["hmm"]
	typestr["euk"] = "eukaryotes"
	typestr["gram-"] = "Gram-negative bacteria"
	typestr["gram+"] = "Gram-positive bacteria"
	printf "Using %s trained on %s\n\n", 
		methodstr[method], typestr[type]

	if (graphics) 
		initgraphics() 
	
	# Read up to the first '>' in the output files:
	if (method ~ "nn") 
		print_entry(output_nn_file)
	if (method ~ "hmm") 
		print_entry(output_hmm_file)
}

/^>/{	#read from FASTAFILE
	header=$0
	id = substr($1,2)
	N++
	
	if (www) 
# udkommenteret den 6 januar 2004 JDB
#		print "<hr>"
		print "<br>"
	else {
		for (i=1;i<=70;i++)
			printf "-"
		print ""
	}		
	if (www) print "<PRE><B>"
	print header
	if (www) print "</B></PRE>"
	print ""
	
	if (method ~ "nn") {
		if (www) print "<b>"
		print "SignalP-NN result:"
		if (www) print "</b><br>"
		if (format=="full")
			writedata(data_nn_prefix "." N)
		if (graphics)	
			plotsequence(data_nn_prefix "." N, "nn", N, id)
		print_entry(output_nn_file)
	}
	
	if (method ~ "hmm") {
		if (www) print "<b>"
		print "SignalP-HMM result:"
		if (www) print "</b><br>"
		if (format=="full")
			writedata(data_hmm_prefix "." N)
		if (graphics)	
			plotsequence(data_hmm_prefix "." N, "hmm", N, id)
		print_entry(output_hmm_file)
	}
	print ""
}

END {
	if (www && graphics) {
		if (graphics == "ps") {
			print "# <A HREF=\"" linkrel psfile "\">Postscript",
				"graphics</A><br>"
		}		
		# Find relative name of gnu script file:
		gfname = gf
		sub(/^.*\//,"",gfname)		
		print "# <A HREF=\"" linkrel gfname "\">gnuplot script</A>",
			"for making the plot(s) <br>"
	}
}


function initgraphics() {
	print   "#Gnuplot commands for SignalP job",
        	"set data style lines",
        	"set yrange [-0.2:1.1]",
        	"set ytics 0,0.2,1",
		"set size 1,0.7",
        	"set format y \"%.1f\"" >gf

	linkrel = www "/" tmpdirname "/" 
	
	if (graphics == "ps") {
		psfile = graph_prefix ".ps"
		print 	"set term post landscape color",
        		"set output \"" psfile "\"" >gf
	}

}	

function writedata(data) {
	if (www) print "<PRE>"
	while ((getline < data) > 0) print
	close(data)
	if (www) print "</PRE>"
}

function plotsequence(data, thismethod, N, id) {
	if (thismethod == "nn") 
		title = "SignalP-NN prediction (" type " networks): " id 
	if (thismethod == "hmm") 
		title = "SignalP-HMM prediction (" type " models): " id 

	print   "set title \"" title "\"",
                "set xlabel \"Position\"",
                "set ylabel \"Score\"",
                "set nolabel" >gf

	while ((getline < data) > 0) {
		if ($0 !~ /^[>#]/)
			printf "set label \"%s\" at %d,-0.1 center\n", 
				$2, $1 >gf
	}
	close(data)

	giffile = graph_prefix "." thismethod "." N ".gif"
	epsfile = graph_prefix "." thismethod "." N ".eps"
	
	if (graphics ~ "gif") {			# gif or gif+eps
		print 	"set output",		# make sure pipe is flushed
			"set term pbm small color",
			"set output \"|" ppmtogif ">" giffile" \"" >gf
	}	
	if (graphics == "eps") {		# eps only
		print   "set term postscript eps color",
			"set output \"" epsfile "\"" >gf
	}		
	
	if (thismethod == "nn") {
	        print   "plot   \"" data "\" u 1:3 t \"C score\" w imp, \\",
                        "       \"" data "\" u 1:4 t \"S score\", \\",
                        "       \"" data "\" u 1:5 t \"Y score\", \\",
                        "       0.5 t \"\" with dots" >gf
	}
	if (thismethod == "hmm") {
#	        print   "plot   \"" data "\" u 1:3 t \"C prob.\" w imp, \\",
#                        "       \"" data "\" u 1:4 t \"S prob.\", \\",
#                        "       0.5 t \"\" with dots" >gf
	        print   "plot   \"" data "\" u 1:3 t \"Cleavage prob.\" w imp, \\",
                        "       \"" data "\" u 1:5 t \"n-region prob.\", \\",
                        "       \"" data "\" u 1:6 t \"h-region prob.\", \\",
                        "       \"" data "\" u 1:7 t \"c-region prob.\" 5, \\",
                        "       0.5 t \"\" with dots" >gf
	}	
	

	if (graphics == "gif+eps") {
		print   "set term postscript eps color",
			"set output \"" epsfile "\"",
			"replot" >gf
	}		

	print "set output" >gf	# HN bugfix 05-12-15, wait for pipe

	if (www) { 	#include graphics	
		print "<br>"
		if (graphics ~ "gif")
			print "<IMG SRC=\"" linkrel giffile "\"><P>";
		print "# <A HREF=\"" linkrel data "\">data</A>" 
		if (graphics ~ "eps")
			print "# <A HREF=\"" linkrel epsfile "\">plot</A>",
				"in EPS format";
		print "<P>"	
	}

}

function print_entry(file) {	
	OK=0
	# Print previous header:
	if (buffer[file]) {
		if (www) print "<PRE>"
		print buffer[file]
	}	
	while ((getline line < file) > 0) {
		OK=1
		if (line ~ /^>/) { 
			if (www && buffer[file]) print "</PRE>"
			buffer[file] = line 
			return 
		}
		else print line
	}
	buffer[file] = ""
	if (www) print "</PRE>"
	if (!OK) {print "No more entries in " file; exit 1}	
}
