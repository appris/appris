#!/usr/bin/R

# libraries
library('getopt')

# Input parameters
spec = matrix(c(
	'inpath',	'i',	1,	"character",
	'specie',	's',	1,	"character",
	'verbose', 	'v',	2,	"integer",
	'help',		'h',	0,	"logical"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Required parametes
if ( is.null(opt$inpath) || is.null(opt$specie) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# if help was asked for print a friendly message and exit with a non-zero error code
if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}
if ( is.null(opt$verbose) ) {
	opt$verbose = FALSE
}
#print some progress messages to stderr, if requested.
if ( opt$verbose ) {
	write("writing...", stderr())
}

###################
# local functions #
###################
fun1 <- function(x,col,fn) {
  m <- fn(x[,col])
  t <- x[x[,col] == m,]
  return (t)
}

fun_equal <- function(x,col,list) {
  out <- NULL
  for ( equ in list ) {
    out <- rbind(out,x[x[,col] == equ,])
  }
  return (out)
}

fun_retrieve_vals <- function(x) {
  if ( !is.na(x) && x != "-" ) {
    n <- as.numeric(sub(".([0-9]+)$","", x))
    d <- as.numeric(sub("^([0-9]+).","", x))
  }
  return (list(n,d))
}

number_ticks <- function(n) {
  function(limits) pretty(limits, n)
}

########
# main #
########

# set workspace

# Lynx
#setwd("/home/jmrodriguez/projects/Lynx/pardinus23A/data/p23A.v3.14Oct2013")
#table_score = read.table("appris_data.appris_score.lynx_pardinus.txt", fill=TRUE, col.names=c("gene_id","transc_id","status","biotype","codons","ccds_id","firestar","matador3d","corsair","spade","inertia","cexonic","thump","crash_sp","crash_tp","appris","aa_len"), colClasses=c("firestar"="numeric","matador3d"="numeric","corsair"="numeric","spade"="character","thump"="character") )
# Human
#setwd("/home/jmrodriguez/projects/Encode/gencode15/data/g15.v4.16Oct2013")
#table_prin = read.table("../g15.v3.15Jul2013/appris_data.principal.homo_sapiens.tsv", col.names=c("chr","gene_id","transc_id","transc_name","ccds_id","appris_label"), stringsAsFactors=FALSE )
#table_score = read.table("appris_data.appris_score.homo_sapiens.txt", fill=TRUE, col.names=c("gene_id","transc_id","transl","status","biotype","codons","ccds_id","firestar","matador3d","corsair","spade","inertia","cexonic","thump","crash_sp","crash_tp","appris","aa_len"), colClasses=c("firestar"="numeric","matador3d"="numeric","corsair"="numeric","spade"="character","thump"="character"),  stringsAsFactors=FALSE )
# Mouse
#setwd("/home/jmrodriguez/projects/Encode/musculus70/data/e70.v4.29Oct2013/")
#table_prin = read.table("appris_data.principal.mus_musculus.tsv", col.names=c("chr","gene_id","transc_id","transc_name","ccds_id","appris_label"), stringsAsFactors=FALSE )
#table_score = read.table("appris_data.appris_score.mus_musculus.txt", fill=TRUE, col.names=c("gene_id","transc_id","transl","status","biotype","codons","ccds_id","firestar","matador3d","corsair","spade","inertia","cexonic","thump","crash_sp","crash_tp","appris","aa_len"), colClasses=c("firestar"="numeric","matador3d"="numeric","corsair"="numeric","spade"="character","thump"="character"),  stringsAsFactors=FALSE )
# Rat
#setwd("/home/jmrodriguez/projects/Encode/norvegicus70/data/e70.v3.10Jul2013/")
#table_prin = read.table("appris_data.principal.rattus_norvegicus.tsv", fill=TRUE, col.names=c("chr","gene_id","transc_id","transc_name","ccds_id","appris_label"), stringsAsFactors=FALSE )
#table_score = read.table("appris_data.appris_score.rattus_norvegicus.txt", fill=TRUE, col.names=c("gene_id","transc_id","transl","status","biotype","codons","ccds_id","firestar","matador3d","corsair","spade","inertia","cexonic","thump","crash_sp","crash_tp","appris","aa_len"), colClasses=c("firestar"="numeric","matador3d"="numeric","corsair"="numeric","spade"="character","thump"="character"),  stringsAsFactors=FALSE )

setwd(opt$inpath)
file_pri <- paste("appris_data.principal.",opt$specie,".tsv", sep="")
file_score <- paste("appris_data.appris_score.",opt$specie,".txt", sep="")
table_prin <- read.table(file_pri, col.names=c("chr","gene_id","transc_id","transc_name","ccds_id","appris_label"), stringsAsFactors=FALSE )
table_score <- read.table(file_score, fill=TRUE, col.names=c("gene_id","transc_id","transl","status","biotype","codons","ccds_id","firestar","matador3d","corsair","spade","inertia","cexonic","thump","crash_sp","crash_tp","appris","aa_len"), colClasses=c("firestar"="numeric","matador3d"="numeric","corsair"="numeric","spade"="character","thump"="character"),  stringsAsFactors=FALSE)

# scan data per gene
genelist <- unique(table_score[,c('gene_id')])
principalstatspergene <- NULL
totalstatspertransc <- NULL
disagreecorsairmatador3d <- NULL
disagreematador3dcorsair <- NULL
for ( g_id in genelist ) {
  
  # init values
  func_residues <- 0
  struc_score <- 0
  conser_score <- 0
  num_domains <- 0
  num_dam_domains <- 0
  num_tmh <- 0
  num_dam_tmh <- 0
  
  # get principal isoform
  principal_id <- table_prin[table_prin$gene_id == g_id,c('transc_id')]
	  
  # extract method values of principal isoform
  firestar <- table_score[table_score$transc_id == principal_id,c('firestar')]
  matador3d <- table_score[table_score$transc_id == principal_id,c('matador3d')]
  corsair <- table_score[table_score$transc_id == principal_id,c('corsair')]
  spade <- table_score[table_score$transc_id == principal_id,c('spade')]
  thump <- table_score[table_score$transc_id == principal_id,c('thump')]
  if ( !is.na(firestar) && firestar != "-" ) {
    func_residues <- firestar
  }
  if ( !is.na(matador3d) && matador3d != "-" ) {
    struc_score <- matador3d
  }
  if ( !is.na(corsair) && corsair != "-" ) {
    conser_score <- corsair
  }
  num_domains <- fun_retrieve_vals(spade)
  num_tmh <- fun_retrieve_vals(thump)
  
  # join row
	#statspergene <- data.frame('gene_id'=g_id,'transc_id'=principal_id,'func_residues'=func_residues,'num_domains'=num_domains,'num_dam_domains'=num_dam_domains,'num_tmh'=num_tmh, 'num_dam_tmh'=num_dam_tmh)
  statspergene <- data.frame('gene_id'=g_id,'transc_id'=principal_id,'func_residues'=func_residues,'num_domains'=as.integer(num_domains[1]),'num_dam_domains'=as.integer(num_domains[2]),'num_tmh'=as.integer(num_tmh[1]), 'num_dam_tmh'=as.integer(num_tmh[2]) )
	principalstatspergene <- rbind(principalstatspergene, statspergene)
  
  # disagrements between Matador3D and CORSAIR
  if ( struc_score == 0 && (conser_score - struc_score > 10) ) {
    disagreepergene <- data.frame('gene_id'=g_id,'transc_id'=principal_id,'struc_score'=struc_score,'conser_score'=conser_score)
    disagreecorsairmatador3d <- rbind(disagreecorsairmatador3d, disagreepergene)
  } else if ( conser_score == 0 && (struc_score - conser_score > 10) ) {
    disagreepergene <- data.frame('gene_id'=g_id,'transc_id'=principal_id,'struc_score'=struc_score,'conser_score'=conser_score)
    disagreematador3dcorsair <- rbind(disagreematador3dcorsair, disagreepergene)
  }
  
  # init values
  transc_ids <- NULL
  transc_funcres <- NULL
  transc_struc <- NULL
  transc_conser <- NULL
  transc_domain <- NULL
  transc_dam_domain <- NULL
  transc_tmh <- NULL
  transc_dam_tmh <- NULL
  
  # extract values per transcript
  transc_rep <- table_score[table_score$gene_id == g_id,c('transc_id','firestar','spade','thump')]
  transc_ids <- transc_rep$transc_id
  transc_funcres <- transc_rep$firestar
  for ( t_spade in transc_rep$spade ) {
    t_domains <- fun_retrieve_vals(as.character(t_spade))
    transc_domain <- rbind(transc_domain,t_domains[[1]])
    transc_dam_domain <- rbind(transc_dam_domain,t_domains[[2]])
  }
  for ( t_thump in transc_rep$thump ) {
    t_tmh <- fun_retrieve_vals(as.character(t_thump))
    transc_tmh <- rbind(transc_tmh,t_tmh[[1]])
    transc_dam_tmh <- rbind(transc_dam_tmh,t_tmh[[2]])
  }
  
  # join row
  statspertransc <- data.frame('gene_id'=g_id, 'principal_id'=principal_id, 'transc_id'=transc_ids, 'func_residues'=transc_funcres, 'num_domains'=transc_domain, 'num_dam_domains'=transc_dam_domain, 'num_tmh'=transc_tmh, 'num_dam_tmh'=transc_dam_tmh )
  totalstatspertransc <- rbind(totalstatspertransc, statspertransc)
  
}

######################################
# FREQ. DISTRIBUTION PER TRANSCRIPTS #
######################################


# Frecuency distribution of Maximum functional residues - firestar (save into pdf)
totalfuncresiduespertransc <- totalstatspertransc$func_residues
bmin <- 1
bmax <- 433  # there are more values: sort(unique(totalfuncresiduespergene)) -> [0 .. 432]
#bmax <- as.integer(range(totalfuncresiduespergene)[2])
bfuncres = seq(bmin,bmax,by=1) # personalized
totalfuncresiduespertransc.cut = cut(totalfuncresiduespertransc, bfuncres, right=FALSE)
totalfuncresiduespertransc.freq = table(totalfuncresiduespertransc.cut)
cbind(totalfuncresiduespertransc.freq)
#pdf(file='freq_funcresiduespertransc.pdf')
yprettyaxis = round(c(0, pretty(range(totalfuncresiduespertransc.freq)[2])[2],10))
bp <- barplot(totalfuncresiduespertransc.freq, col="red", main="Frequency of num. functional residues per transcript", xlab="Num. functional residues", ylab="Num. transcripts", yaxp=yprettyaxis)
text(bp, bmin, round(totalfuncresiduespertransc.freq), cex=0.8, pos=3)
#dev.off()



# Frecuency distribution of pfam domains - spade (save into pdf)
domainspertransc <- totalstatspertransc$num_domains
damdomainspertransc <- totalstatspertransc$num_dam_domains
uniondomain <- c(unique(domainspertransc),unique(damdomainspertransc))
bmin <- 1
bmax <- 50 # there are more values: sort(unique(uniondomain)) -> [0 .. 80, 265]
#bmax <- as.integer(range(uniondomain)[2])
#buniondomain = c(seq(bmin,5,by=1),pretty(bmax+1,20))
buniondomain = seq(bmin,bmax,by=1) # personalized
domainspertransc.cut = cut(domainspertransc, buniondomain, right=FALSE)
domainspertransc.freq = table(domainspertransc.cut)
cbind(domainspertransc.freq)
damdomainspertransc.cut = cut(damdomainspertransc, buniondomain, right=FALSE)
damdomainspertransc.freq = table(damdomainspertransc.cut)
cbind(damdomainspertransc.freq)
totaldompergene <- rbind(domainspertransc.freq, damdomainspertransc.freq)
#pdf(file='freq_domainspertransc.pdf')
yprettyaxis = c(bmin, pretty(range(totaldompergene)[2])[2],10)
bp <- barplot(totaldompergene, col=c("red","lightblue"), main="Frequency of num. domains per transcript", xlab="Num. domains", ylab="Num. transcript", legend=c("whole domains", "damaged domains"), yaxp=yprettyaxis, beside=TRUE)
text(bp, bmin, round(totaldompergene), cex=0.8, pos=3)
#dev.off()



# Frecuency distribution of trans-membrane helices - thump (save into pdf)
tmhpertransc <- totalstatspertransc$num_tmh
damtmhpertransc <- totalstatspertransc$num_dam_tmh
uniontmh <- c(unique(tmhpertransc),unique(damtmhpertransc))
bmin <- 1
bmax <- 21
#bmax <- as.integer(range(uniontmh)[2])
buniontmh = seq(bmin,bmax,by=1) # personalized
tmhpertransc.cut = cut(tmhpertransc, buniontmh, right=FALSE)
tmhpertransc.freq = table(tmhpertransc.cut)
cbind(tmhpertransc.freq)
damtmhpertransc.cut = cut(damtmhpertransc, buniontmh, right=FALSE)
damtmhpertransc.freq = table(damtmhpertransc.cut)
cbind(damtmhpertransc.freq)
totaltmhpergene <- rbind(tmhpertransc.freq, damtmhpertransc.freq)
#pdf(file='freq_tmhpertransc.pdf')
yprettyaxis = round(c(0, pretty(range(totaltmhpergene)[2])[2],10))
bp <- barplot(totaltmhpergene, col=c("red","lightblue"), main="Frequency of num. trans-membrane helices per transcript", xlab="Num. trans-membrane helices", ylab="Num. genes", legend=c("trans-membrane helices", "damaged trans-membrane helices"), yaxp=yprettyaxis, beside=TRUE)
text(bp, bmin, round(totaltmhpergene), cex=0.8, pos=3)
#dev.off()



###############################
# FREQ. DISTRIBUTION PER GENE #
###############################

# Frecuency distribution of Maximum functional residues - firestar (save into pdf)
totalfuncresiduespergene <- principalstatspergene$func_residues
bmin <- 1
bmax <- 432  # there are more values: sort(unique(totalfuncresiduespergene)) -> [0 .. 432]
#bmax <- as.integer(range(totalfuncresiduespergene)[2])
bfuncres = seq(bmin,bmax,by=1) # personalized
totalfuncresiduespergene.cut = cut(totalfuncresiduespergene, bfuncres, right=FALSE)
totalfuncresiduespergene.freq = table(totalfuncresiduespergene.cut)
cbind(totalfuncresiduespergene.freq)
#pdf(file='freq_funcresiduespergene.pdf')
yprettyaxis = round(c(0, pretty(range(totalfuncresiduespergene.freq)[2])[2],10))
bp <- barplot(totalfuncresiduespergene.freq, col="red", main="Frequency of num. functional residues per gene", xlab="Num. functional residues", ylab="Num. genes", yaxp=yprettyaxis)
text(bp, bmin, round(totalfuncresiduespergene.freq), cex=0.8, pos=3)
#dev.off()


# Frecuency distribution of pfam domains - spade (save into pdf)
domainspergene <- principalstatspergene$num_domains
damdomainspergene <- principalstatspergene$num_dam_domains
uniondomain <- c(unique(domainspergene),unique(damdomainspergene))
bmin <- 1
bmax <- 50 # there are more values: sort(unique(uniondomain)) -> [0 .. 80, 265]
#bmax <- as.integer(range(uniondomain)[2])
#buniondomain = c(seq(bmin,5,by=1),pretty(bmax+1,20))
buniondomain = seq(bmin,bmax,by=1) # personalized
domainspergene.cut = cut(domainspergene, buniondomain, right=FALSE)
domainspergene.freq = table(domainspergene.cut)
cbind(domainspergene.freq)
damdomainspergene.cut = cut(damdomainspergene, buniondomain, right=FALSE)
damdomainspergene.freq = table(damdomainspergene.cut)
cbind(damdomainspergene.freq)
totaldompergene <- rbind(domainspergene.freq, damdomainspergene.freq)
#pdf(file='freq_domainspergene.pdf')
yprettyaxis = round(c(0, pretty(range(totaldompergene)[2])[2],10))
bp <- barplot(totaldompergene, col=c("red","lightblue"), main="Frequency of num. domains per gene", xlab="Num. domains", ylab="Num. genes", legend=c("whole domains", "damaged domains"), yaxp=yprettyaxis, beside=TRUE)
text(bp, bmin, round(totaldompergene), cex=0.8, pos=3)
#dev.off()



# Frecuency distribution of trans-membrane helices - thump (save into pdf)
tmhpergene <- principalstatspergene$num_tmh
damtmhpergene <- principalstatspergene$num_dam_tmh
uniontmh <- c(unique(tmhpergene),unique(damtmhpergene))
bmin <- 1
bmax <- 21
#bmax <- as.integer(range(uniontmh)[2])
buniontmh = seq(bmin,bmax,by=1) # personalized
tmhpergene.cut = cut(tmhpergene, buniontmh, right=FALSE)
tmhpergene.freq = table(tmhpergene.cut)
cbind(tmhpergene.freq)
damtmhpergene.cut = cut(damtmhpergene, buniontmh, right=FALSE)
damtmhpergene.freq = table(damtmhpergene.cut)
cbind(damtmhpergene.freq)
totaltmhpergene <- rbind(tmhpergene.freq, damtmhpergene.freq)
#pdf(file='freq_tmhpergene.pdf')
yprettyaxis = round(c(0, pretty(range(totaltmhpergene)[2])[2],10))
bp <- barplot(totaltmhpergene, col=c("red","lightblue"), main="Frequency of num. trans-membrane helices per gene", xlab="Num. trans-membrane helices", ylab="Num. genes", legend=c("trans-membrane helices", "damaged trans-membrane helices"), yaxp=yprettyaxis, beside=TRUE)
text(bp, bmin, round(totaltmhpergene), cex=0.8, pos=3)
#dev.off()




# Disagrement between corsair and matador3d
slices <- c(length(genelist),length(unique(disagreecorsairmatador3d$gene_id)))
lbls <- c("Agree", "Disagree")
pct <- round(slices/sum(slices)*100, 2)
lbls <- paste(lbls, slices) # add value to labels
lbls <- paste(lbls, " (",pct,sep="") # add percents to labels
lbls <- paste(lbls,"%) ",sep="") # ad % to labels
#pdf(file='dis_corsairmatador3d.pdf')
pie(slices,labels = lbls, col=rainbow(length(lbls)),main="Disagrement between CORSAIR and Matador3d")
#dev.off()





# Disagrement between matador3d and corsair
slices <- c(length(genelist),length(unique(disagreematador3dcorsair$gene_id)))
lbls <- c("Agree", "Disagree")
pct <- round(slices/sum(slices)*100, 2)
lbls <- paste(lbls, slices) # add value to labels
lbls <- paste(lbls, " (",pct,sep="") # add percents to labels
lbls <- paste(lbls,"%) ",sep="") # ad % to labels
#pdf(file='dis_matador3dcorsair.pdf')
pie(slices,labels = lbls, col=rainbow(length(lbls)),main="Disagrement between Matador3d and CORSAIR")
#dev.off()

#signal success and exit.
q(status=0);
