#! /usr/bin/nawk -f

BEGIN {
	nparts = 5	# No. of partitions to average over
	nf = 9		# No. of fields in each partition file
	longpred["S"] = "Signal peptide"
	longpred["A"] = "Signal anchor"
	longpred["Q"] = "Non-secretory protein"
	Cmaxcut = 0.5  # <  these two are arbitrary
	Stotcut = 0.5  # <
	
        if (format=="short") {
                printf "                # SignalP-HMM %s predictions\n", type
                print "# name      !  Cmax  pos ?  "\
                        "Sprob ?"
        }
	
}	

/^#/ { 
	if (N>0) seqout()
	N++
	name = $2
	pos = 0		# Initialize position
	Cmax = Cmaxpos = 0
	datafile = datafileprefix "." N
	if (format != "short") {
		print ">" name > datafile
		print "# pos  aa    C       S      n-reg   h-reg   c-reg" > datafile
	}	
	getline		# Read the label header line
	if (!columns_set) {
		N_columns = 9
		for (i=2; i<=N_columns; i++) 
			col[$i] = i
		columns_set = 1
	}
	next
}

{
	pos++
	C = S = A = Q = nreg = hreg = creg = 0
	for (i=0; i<nparts; i++) {
		C += $(nf*i+col["C"])		
		# Signal peptide score is S (ini) +n+H+c
		S += $(nf*i+col["S"]) + $(nf*i+col["n"]) + $(nf*i+col["H"]) + $(nf*i+col["c"])
		A += $(nf*i+col["A"])
		Q += $(nf*i+col["Q"])
		nreg += $(nf*i+col["S"]) + $(nf*i+col["n"])
		hreg += $(nf*i+col["H"])
		creg += $(nf*i+col["c"])
	}
	C = C/nparts
	if (C>Cmax) { Cmax=C; Cmaxpos=pos }
	S = S/nparts
	A = A/nparts
	Q = Q/nparts
	nreg = nreg/nparts; hreg = hreg/nparts; creg = creg/nparts
	if (pos == 1) { 
		Stot = S
		Atot = A
		Qtot = Q
	}	
	if (format != "short") {
		printf "%5d   %1s   %5.3f   %5.3f   %5.3f   %5.3f   %5.3f\n", 
			pos, $1, C, S, nreg, hreg, creg > datafile
	}
}

function seqout() {
	close (datafile)

	if (Stot > Atot && Stot >  Qtot) pred = "S"
	else if (Atot > Stot && Atot >  Qtot) pred = "A"
	else pred = "Q"
	
	if (format=="short") {
		printf "%-10s  ", name
		printf "%1s  ", pred
		printf "%5.3f %3d %1s  ", 
			Cmax, Cmaxpos, (Cmax>Cmaxcut) ? "Y" : "N"
		printf "%5.3f %1s  \n", 
			Stot, (Stot>Stotcut) ? "Y" : "N"
		return
	}
	
	print ">" name
	print "Prediction:", longpred[pred] 
	printf "Signal peptide probability: %5.3f\n", Stot
	if (type == "euk") 
		printf "Signal anchor probability: %5.3f\n", Atot
	printf "Max cleavage site probability: %5.3f between pos. %2d and %2d\n", 
		Cmax, Cmaxpos-1, Cmaxpos
}

END {
	seqout()
}
