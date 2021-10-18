### known issues:
# may not handle variants with multiple annotations


##################################
# Load packages and parse inputs #
##################################

library("optparse")

# define argument list and data types
option_list <- list(
  make_option(c("--mutations"), type="character", default = NULL, 
              help="Mutation table", metavar="character"),
  make_option(c("--tumor_ID"), type="character", default = NULL, 
               help="Used as prefix", metavar="character"),
  make_option(c("--hotspots"), type="character", default = NULL, 
              help="Hotspot table (single amino acid hotspots)", metavar="character"),
  make_option(c("--noncoding_hotspots"), type="character", default = NULL, 
              help="Non-coding hotspot table", metavar="character"),
  make_option(c("--splicesite_hotspots"), type="character", default = NULL, 
              help="Splice hotspot table", metavar="character"),
  make_option(c("--indel_hotspots"), type="character", default = NULL,
              help="Indel hotspot table", metavar="character"),
  make_option(c("--skip_3D"), type="character", default = "yes",
              help="[yes|no] Omit 3D hotspots. Default yes.", metavar="character")
)

# parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

# check if mandatory arguments were provided, show help menu otherwise
if ( is.null(opt[["mutations"]]) | is.null(opt[["tumor_ID"]]) | is.null(opt[["hotspots"]]) | is.null(opt[["splicesite_hotspots"]]) | is.null(opt[["indel_hotspots"]]) ) {
  print_help(opt_parser)
  print(opt)
  stop("One or more mandatory arguments missing.", call=FALSE)
}


muts <- read.delim(as.character(opt[["mutations"]]), as.is = T)
colnames(muts) <- sub("ANN....", "", colnames(muts))

# split AA col
muts$AA.ref <- gsub("[0-9].+", "", sub("p.", "", muts$AA))
muts$AA.pos <- gsub("p.[^0-9]+([0-9]+).+", "\\1", muts$AA)

muts$AA.ref.end <- gsub("[0-9].+", "", sub("^[^_]+_", "", sub("p.", "", muts$AA)))
muts$AA.pos.end <- gsub("[^0-9]+([0-9]+).+", "\\1", sub("^[^_]+_", "", sub("p.", "", muts$AA)))

# for splice sites, use the codon position. Extract the position after "c." -> can be followed by "-" or "*". Omit "n."
muts$c.pos <- "."
muts$c.pos[grepl("^c\\.", muts$HGVS_C)] <- gsub("^(c\\.[0-9\\*-][0-9]*).+$", "\\1", muts$HGVS_C[grepl("^c\\.", muts$HGVS_C)])



########################
### exonic hotspots  ###
########################
hotspots <- read.delim(as.character(opt[["hotspots"]]), as.is = T)

# if the mutation input has no 'chr' prefix, remove it from the hotspot lists
if(!any(grepl("chr", muts$CHROM))){
  hotspots$CHROM[which(grepl("chr", hotspots$CHROM))] <- gsub("chr", "",hotspots$CHROM[which(grepl("chr", hotspots$CHROM))])
}

# remove hotspot3D
if ( opt[["skip_3D"]] == "yes" ) {
  hotspots <- hotspots[which(hotspots$hotspot_type != "HOTSPOT3D"),]
}

# go over each mut row and fetch same chrom/gene/aa_ref/aa_pos in hotspots. 
# Consider only missense_variant/stop_gained/initiator_codon_variant/start_lost because other effects are not real single aa hotspot changes, and we don't want to be tagging those. Omit indels.
muts$hotspot <- apply(muts, 1, function(x){
  if (grepl("frameshift|inframe", x["EFFECT"])) {
    "."
  } else if (!grepl("missense_variant|stop_gained|initiator_codon_variant|start_lost", x["EFFECT"])) {
      "."
  } else if (length(as.character(unlist(hotspots[which(hotspots$CHROM == x["CHROM"] & hotspots$GENE == x["GENE"] & hotspots$aa_ref == x["AA.ref"] & hotspots$aa_pos == x["AA.pos"]),which(colnames(hotspots) == "hotspot_type")]))) == 0) {
    "."
  } else {
    as.character(unlist(hotspots[which(hotspots$CHROM == x["CHROM"] & hotspots$GENE == x["GENE"] & hotspots$aa_ref == x["AA.ref"] & hotspots$aa_pos == x["AA.pos"]), which(colnames(hotspots) == "hotspot_type")]))
  }
})

# add "aa" sufix to muts$hot so we can distinguish it downstream
muts$hotspot <- paste0(muts$hotspot, "aa")

# save VCF where muts$hot had a hit, split HOTSPOT and HOTSPOT3D
# for <SnpSift annotate> we need VCF cols #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
# FILTER must be all PASS, everything that is in INFO will be passed to the annotated vcf.
# if there are no hotspot hits, save empty header so that we have an output for every sample
if (length(rownames(muts[grepl("HOTSPOT", muts$hotspot),])) > 0) {
  vcf_hotspot <- muts[grepl("HOTSPOT", muts$hotspot),]
  vcf_hotspot$QUAL <- "."
  vcf_hotspot$FILTER <- "PASS"
  vcf_hotspot$INFO <- paste0(vcf_hotspot$hotspot, ";", vcf_hotspot$hotspot, "_GENE=", vcf_hotspot$GENE, ";", vcf_hotspot$hotspot, "_AA=", vcf_hotspot$AA.ref, vcf_hotspot$AA.pos)
  vcf_hotspot <- vcf_hotspot[c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]
  colnames(vcf_hotspot) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  write.table(vcf_hotspot, file = paste0(opt[["tumor_ID"]],"_hotspot.vcf"), row.names = F, quote = F, sep = '\t')
} else {
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = paste0(opt[["tumor_ID"]],"_hotspot.vcf"))
}


#############################
### non-coding hotspots   ###
#############################
noncoding <- read.delim(as.character(opt[["noncoding_hotspots"]]), as.is = T)

# if the mutation input has no 'chr' prefix, remove it from the hotspot lists
if(!any(grepl("chr", muts$CHROM))){
  noncoding$CHROM[which(grepl("chr", noncoding$CHROM))] <- gsub("chr", "",noncoding$CHROM[which(grepl("chr", noncoding$CHROM))])
}

# noncoding hotspots need to be an exact match of chr, pos and alt base. Just save the 'noncoding hotspots' rows with an exact match.
noncoding$matchID <- paste(noncoding$CHROM,noncoding$POS,noncoding$REF,noncoding$ALT,sep=";")
muts$matchID <- paste(muts$CHROM, muts$POS, muts$REF, muts$ALT,sep=";")

if(any(muts$matchID %in% noncoding$matchID)){
  vcf_hotspotNC <- noncoding[which(noncoding$matchID %in% muts$matchID),]
  vcf_hotspotNC$FILTER <- "PASS"
  vcf_hotspotNC <- vcf_hotspotNC[c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]
  colnames(vcf_hotspotNC) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  write.table(vcf_hotspotNC, file = paste0(opt[["tumor_ID"]],"_hotspot_noncoding.vcf"), row.names = F, quote = F, sep = '\t')
} else {
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = paste0(opt[["tumor_ID"]],"_hotspot_noncoding.vcf"))
}


############################
### splicesite hotspots  ###
############################
splice <- read.delim(as.character(opt[["splicesite_hotspots"]]), as.is = T)

# if the mutation input has no 'chr' prefix, remove it from the hotspot lists
if(!any(grepl("chr", muts$CHROM))){
  splice$CHROM[which(grepl("chr", splice$CHROM))] <- gsub("chr", "",splice$CHROM[which(grepl("chr", splice$CHROM))])
}

# go over each mut row and fetch same chrom/gene/pos/ in splice. Consider only splice_acceptor_variant/splice_donor_variant because other effects are not real splicing effects, and we don't want to be tagging those.

muts$hotsplice <- apply(muts, 1, function(x){
  if (!grepl("splice_acceptor_variant|splice_donor_variant", x["EFFECT"])) {
    "."
  } else if (length(as.character(unlist(splice[which(splice$CHROM == x["CHROM"] & splice$GENE == x["GENE"] & splice$c_pos == x["c.pos"]), which(colnames(splice) == "hotspot_type")]))) == 0) {
    "."
  } else {
    as.character(unlist(splice[which(splice$CHROM == x["CHROM"] & splice$GENE == x["GENE"] & splice$c_pos == x["c.pos"]), which(colnames(splice) == "hotspot_type")]))
  }
})

# add "sp" sufix to muts$hot so we can distinguish it downstream
muts$hotsplice <- paste0(muts$hotsplice, "sp")

# save VCF where muts$hotsplice had a hit
# for <SnpSift annotate> we need VCF cols #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
# FILTER must be all PASS, everything that is in INFO will be passed to the annotated vcf.
# if ther are no hits, save empty header so that we have an output for every sample
if (length(rownames(muts[grepl("HOTSPOT", muts$hotsplice),])) > 0) {
  vcf_hotsplice <- muts[grepl("HOTSPOT", muts$hotsplice),]
  vcf_hotsplice$QUAL <- "."
  vcf_hotsplice$FILTER <- "PASS"
  vcf_hotsplice$INFO <- paste0(vcf_hotsplice$hotsplice, ";", vcf_hotsplice$hotsplice, "_GENE=", vcf_hotsplice$GENE, ";", vcf_hotsplice$hotsplice, "_c.pos=", vcf_hotsplice$c.pos)
  vcf_hotsplice <- vcf_hotsplice[c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]
  colnames(vcf_hotsplice) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  write.table(vcf_hotsplice, file = paste0(opt[["tumor_ID"]],"_hotspot_sp.vcf"), row.names = F, quote = F, sep = '\t')
} else {
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = paste0(opt[["tumor_ID"]],"_hotspot_sp.vcf"))
}




###############################
### in-frame indel hotspots ###
###############################

indel <- read.delim(as.character(opt[["indel_hotspots"]]), as.is = T)

# if the mutation input has no 'chr' prefix, remove it from the hotspot lists
if(!any(grepl("chr", muts$CHROM))){
  indel$CHROM[which(grepl("chr", indel$CHROM))] <- gsub("chr", "",indel$CHROM[which(grepl("chr", indel$CHROM))])
}

# go over each mut row and fetch same chrom/gene/aa_start/aa_end/ in indel. Consider only inframe indels, we don't want to be tagging other frameshifts.

muts$hotindel <- apply(muts, 1, function(x){
  if (!grepl("inframe", x["EFFECT"])) {
    "."
  } else if (length(as.character(unlist(indel[which(indel$CHROM == x["CHROM"] & indel$GENE == x["GENE"] & indel$aa_start <= x["AA.pos"] & indel$aa_end >= x["AA.pos.end"]), which(colnames(indel) == "hotspot_type")]))) == 1) {
    "HOTSPOT"
  } else {
  "."
  }
})

# add "indel" sufix to muts$hot so we can distinguish it downstream
muts$hotindel <- paste0(muts$hotindel, "indel")

# save VCF where muts$hotindel had a hit
# for <SnpSift annotate> we need VCF cols #CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO
# FILTER must be all PASS, everything that is in INFO will be passed to the annotated vcf.
# if not hits, save empty header so that we have an output for every sample
if (length(rownames(muts[grepl("HOTSPOT", muts$hotindel),])) > 0) {
  vcf_hotindel <- muts[grepl("HOTSPOT", muts$hotindel),]
  vcf_hotindel$QUAL <- "."
  vcf_hotindel$FILTER <- "PASS"
  vcf_hotindel$INFO <- paste0(vcf_hotindel$hotindel, ";", vcf_hotindel$hotindel, "_GENE=", vcf_hotindel$GENE)
  vcf_hotindel <- vcf_hotindel[c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]
  colnames(vcf_hotindel) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  write.table(vcf_hotindel, file = paste0(opt[["tumor_ID"]],"_hotspot_indel.vcf"), row.names = F, quote = F, sep = '\t')
} else {
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = paste0(opt[["tumor_ID"]],"_hotspot_indel.vcf"))
}
