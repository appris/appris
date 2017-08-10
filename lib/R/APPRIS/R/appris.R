################
# Dependences #
###############
library(R6)
library(plyr)

#' APPRIS class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords principal isoforms
#' @format \code{\link{R6Class}} object.
#' @examples
#' appris = APPRIS$new("/Users/jmrodriguez/projects/APPRIS/data/homo_sapiens/ens81.v15/")
#' appris_m28a3$gene_decision()

APPRIS <- R6Class( "APPRIS",

# PUBLIC members

public = list(

  # Attributes
  dir = NA,
  df_scores = NA,
  df_labels = NA,
  genes = NA,
  num_transl = NA,
  num_ccds = NA,
  # methods = c("firestar","matador3d","spade","corsair","thump","crash","appris"),
  methods = c("firestar","matador3d","spade","corsair","appris"),

# Create the class
  initialize = function(path) {
    setwd(path)
    # create dataframes
    df_scores = read.table('appris_data.appris.txt',
                               sep="\t",
                               fill = TRUE,
                               col.names=c(
                                 "gene_id",
                                 "gene_name",
                                 "transc_id",
                                 "transl_id",
                                 "translation",
                                 "flags",
                                 "no_codons",
                                 "ccds_id",
                                 "tsl",
                                 "transl_len",
                                 "firestar",
                                 "matador3d",
                                 "corsair",
                                 "spade",
                                 "thump",
                                 "crash",
                                 "inertia",
                                 "proteo",
                                 "score",
                                 "appris"),
                               stringsAsFactors=TRUE )
    scores = df_scores[df_scores$translation == "TRANSLATION",]

    df_labels = read.table('appris_data.appris_label.txt',
                       sep="\t",
                       fill = TRUE,
                       col.names=c(
                         "gene_id",
                         "gene_name",
                         "transc_id",
                         "transl_id",
                         "translation",
                         "flags",
                         "no_codons",
                         "ccds_id",
                         "tsl",
                         "transl_len",
                         "firestar",
                         "matador3d",
                         "corsair",
                         "spade",
                         "thump",
                         "crash",
                         "inertia",
                         "proteo",
                         "score",
                         "appris"),
                       stringsAsFactors=TRUE )
    labels = df_labels[df_labels$translation == "TRANSLATION",]

    # get the list of all translated genes with gene name, variants and ccds
    genes =  labels[labels$translation == "TRANSLATION",c("gene_id","gene_name","transc_id","ccds_id")]

    # get the genes with the num. transcript
    num_transl =  count( unique(genes[,c("gene_id","transc_id") ])$gene_id )
    colnames(num_transl)[1] = "gene_id"
    colnames(num_transl)[2] = "num_transl"

    # get the genes with the num. ccds
    c =  count( unique(genes[genes$ccds_id != "-",c("gene_id","ccds_id")])$gene_id )
    colnames(c)[1] = "gene_id"
    colnames(c)[2] = "num_ccds"
    g = data.frame( gene_id=unique(genes$gene_id) )
    num_ccds = merge.data.frame( g, c, by = "gene_id", all = TRUE)
    num_ccds[ is.na(num_ccds$num_ccds), c("num_ccds")] = 0

    self$dir = path
    self$df_labels = labels
    self$df_scores = scores
    self$genes = genes
    self$num_transl = num_transl
    self$num_ccds = num_ccds
    self$greet()
  },
  greet = function() {
    cat(paste0("APPRIS class has been initialized from ", self$dir, " files.\n"))
  },
# gene_decision()
# Assign for each gene if there is a decision in the method.
  gene_decision = function() {
    labels = self$df_labels
    genes = self$genes
    num_transl = self$num_transl
    num_ccds = self$num_ccds

    # init df with the list of all translated genes, and gene_names, and methods
    # merge with the num. translated transcripts
    # merge with the num. ccds
    decision = unique( data.frame(gene_id=genes$gene_id, gene_name=genes$gene_name) )
    decision = merge.data.frame( decision, num_transl, by = "gene_id", all = TRUE)
    decision = merge.data.frame( decision, num_ccds, by = "gene_id", all = TRUE)
    decision$firestar = "NO"
    decision$matador3d = "NO"
    decision$spade = "NO"
    decision$corsair = "NO"
    decision$appris = NA

    # get the list of genes where a method gives a decision (at least say one NO)
    gd_f = labels[labels$firestar == "NO",c("gene_id")]
    gd_m = labels[labels$matador3d == "NO",c("gene_id")]
    gd_s = labels[labels$spade == "NO",c("gene_id")]
    gd_c = labels[labels$corsair == "NO",c("gene_id")]
    # fill with "YES" for the genes where there was a decision
    decision$firestar[ decision$gene_id %in% gd_f ] = "YES"
    decision$matador3d[ decision$gene_id %in% gd_m ] = "YES"
    decision$spade[ decision$gene_id %in% gd_s ] = "YES"
    decision$corsair[ decision$gene_id %in% gd_c ] = "YES"

    # get the list of genes where a method gives a decision (at least say one NO)
    gd_a1 = unique( labels[labels$appris == "PRINCIPAL:1",c("gene_id")] )
    gd_a2 = unique( labels[labels$appris == "PRINCIPAL:2",c("gene_id")] )
    gd_a3 = unique( labels[labels$appris == "PRINCIPAL:3",c("gene_id")] )
    gd_a4 = unique( labels[labels$appris == "PRINCIPAL:4",c("gene_id")] )
    gd_a5 = unique( labels[labels$appris == "PRINCIPAL:5",c("gene_id")] )
    # fill with "YES" for the genes where there was a decision
    decision$appris[ decision$gene_id %in% gd_a1 ] = "PRINCIPAL:1"
    decision$appris[ decision$gene_id %in% gd_a2 ] = "PRINCIPAL:2"
    decision$appris[ decision$gene_id %in% gd_a3 ] = "PRINCIPAL:3"
    decision$appris[ decision$gene_id %in% gd_a4 ] = "PRINCIPAL:4"
    decision$appris[ decision$gene_id %in% gd_a5 ] = "PRINCIPAL:5"

    # convert method labels to factor
    decision$firestar = as.factor(decision$firestar)
    decision$matador3d = as.factor(decision$matador3d)
    decision$spade = as.factor(decision$spade)
    decision$corsair = as.factor(decision$corsair)
    decision$appris = as.factor(decision$appris)

    return (decision)
  },
# gene_decision()
# Assign for each gene if there is a decision in the method.
  ccds_reject = function() {
    labels = self$df_labels
    genes = self$genes
    num_transl = self$num_transl
    num_ccds = self$num_ccds
    methods = self$methods

    # get the list of genes with unique CCDS
    a = self$gene_decision()
    g = a[a$num_ccds == 1,]

    # genes/transcripts/methos for unique CCDS
    b = labels[ labels$gene_id %in% g$gene_id, ]

    # get the list of genes where a method rejexts unique CCDS
    decision = b[ b$ccds_id != '-' & (b$firestar == "NO" | b$matador3d == "NO" | b$spade == "NO" | b$corsair == "NO" | b$appris == "MINOR"),]
    decision[, methods] = sapply(decision[, methods], as.character)
    decision[, methods][decision[, methods] == "YES"] = "PRIN"
    decision[, methods][decision[, methods] == "MINOR"] = "YES"
    decision[, methods][decision[, methods] == "NO"] = "YES"
    decision[, methods][decision[, methods] != "YES"] = "NO"
    decision[,methods] = lapply(decision[,methods], as.factor)

    return (decision)
  },
# appris_decision()
# Get the Num. genes with PRINCIPAL:1,PRINCIPAL:2, etc..
  appris_decision = function() {
    gd = self$gene_decision()

    # get the Num. genes with PRINCIPAL:1,PRINCIPAL:2, etc..
    g_single = nrow( gd[gd$num_transl == 1 & gd$appris=="PRINCIPAL:1", ])
    g_p1 = nrow( gd[gd$num_transl > 1 & gd$appris=="PRINCIPAL:1", ])
    g_p2 = nrow( gd[gd$num_transl > 1 & gd$appris=="PRINCIPAL:2", ])
    g_p3 = nrow( gd[gd$num_transl > 1 & gd$appris=="PRINCIPAL:3", ])
    g_p4 = nrow( gd[gd$num_transl > 1 & gd$appris=="PRINCIPAL:4", ])
    g_p5 = nrow( gd[gd$num_transl > 1 & gd$appris=="PRINCIPAL:5", ])
    g_sm = nrow( gd )

    df = data.frame(
      label = c('PRINCIPAL:5', 'PRINCIPAL:4', 'PRINCIPAL:3', 'PRINCIPAL:2', 'PRINCIPAL:1', 'single genes'),
      value = c(g_p5/g_sm, g_p4/g_sm, g_p3/g_sm, g_p2/g_sm, g_p1/g_sm, g_single/g_sm)
    )
    df$label = factor(df$label, levels = df$label)

    return(df)
  },
# appris_decision()
# Get the Num. genes with PRINCIPAL:1,PRINCIPAL:2, etc..
  graph_appris_decision = function(dsname) {
    df = self$appris_decision()
    df$pos = rev( cumsum(rev(df$value)) -0.02)

    bp = ggplot(df, aes(x=dsname, y=value, fill=label, width = 0.5))+
      geom_bar(width=1, colour="black", stat="identity", position="fill") +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) +
      scale_fill_brewer(palette="PuBu") +
      labs(fill="APPRIS annotation") +
      scale_y_continuous(labels = percent_format()) +
      geom_text(aes(x = c(1.3, 1.3, 1, 1, 1, 1), y = pos, label = percent(value)), size=4)
    bp

    return(bp)
  }

),
###################
# PRIVATE members #
###################
private = list(
  queue = list(),
  length = function() base::length(private$queue)
)
)
