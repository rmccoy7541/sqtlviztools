library(dplyr)
library(purrr)
library(data.table)
library(leafcutter)
library(stringr)
library(readr)
library(sqtlviztools)
library(qvalue)

setwd("~/Desktop/sqtlviztools-master/")

permutation_sqtl_file <- "WHLBLD.permutations_full.txt.gz"
vcf_file <- "sqtl.vcf"
cluster_counts_file <- "Ne-sQTL_perind_numers.counts.gz"
srr2subject_file <- "srr2subject.txt"
sprime_file <- "sprime_calls.txt.gz"

## reformat VCF
input_vcf <- fread(vcf_file) %>%
  setnames(., "#CHROM", "CHROM")
input_vcf[, 10:ncol(input_vcf)] <- as.data.frame(apply(input_vcf[, 10:ncol(input_vcf)], MAR = c(1, 2), FUN = function(x) if (x == "0/0") {return(0)} else if (x == "0/1") {return(1)} else if (x == "1/1") {return(2)} else if (x == "./.") {return("NA")} ))

variant_list <- input_vcf$ID

input_vcf[, ID := paste(ID, CHROM, POS, sep = ".")]

input_vcf[, CHROM := as.character(CHROM)]
input_vcf[!grepl("chr", CHROM), CHROM := paste0("chr", CHROM)]
setnames(input_vcf, gsub("_.*", "", colnames(input_vcf)))

fwrite(input_vcf, "example/example.vcf", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## reformat cluster counts
input_cc <- fread(cluster_counts_file)
input_cc[!(grepl("^chr", V1)), V1 := paste0("chr", V1)]

srr2subject <- fread(srr2subject_file)
new_colnames <- srr2subject$submitted_subject_id[match(colnames(input_cc), srr2subject$Run)]
setnames(input_cc, new_colnames)

fwrite(input_cc, "example/cluster_counts_temp.txt", sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
system("head -1 example/cluster_counts_temp.txt | cut -f2- > example/cluster_counts.txt")
system("sed -e1,1d example/cluster_counts_temp.txt >> example/cluster_counts.txt")
system("rm example/cluster_counts_temp.txt")

## reformat permutation-pass sQTL results
gtp <- fread(permutation_sqtl_file) %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                "p", "beta", "emp_p", "adj_p")) %>%
  setorder(., adj_p)

gtp[, qval := qvalue(gtp$adj_p)$qvalues]

neand <- fread(sprime_file)[vindija_match == "match" | altai_match == "match"] %>%
  mutate(., var_id = paste(CHROM, POS, REF, ALT, "b37", sep = "_")) %>%
  as.data.table()

gtp[, logP := -log10(adj_p)]
setorder(gtp, logP)
gtp[, expectedP := rev(-log10(ppoints(n = length(gtp$adj_p))))]
gtp[, is_neand := variant_id %in% neand$var_id]
gtp[, dummy2 := paste(variant_id, variant_chrom, var_start, sep = ".")]
gtp[, dummy3 := paste(variant_chrom, var_start, sep = ":")]

gtp <- gtp[variant_id %in% variant_list]

pp_all <- gtp[, c("intron_cluster", "dummy2", "adj_p")] %>%
  setnames(., c("pid", "dummy2", "bpval"))
fwrite(pp_all, "example/example_permutations.all.0.05.bh.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# https://groups.google.com/forum/#!topic/leafcutter-users/Hb0S2aZBWtw
pp_full <- gtp[, c("intron_cluster", "dummy3", "distance", "adj_p", "beta")]
fwrite(pp_full, "example/example_permutations.full.0.05.bh.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

input_ids <- gtp[, c("dummy3", "variant_id")] %>%
  setnames(., c("dummy2", "RS_id"))
fwrite(input_ids, "example/snp_ids_rs_ids.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## reformat bed files
exons <- fread("gencode_hg19_all_exons.bed.gz")
fwrite(exons, "example/gencode_hg19_all_exons.bed", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

introns <- fread("gencode_hg19_all_introns.bed.gz")[, c(1:7, 9, 8, 10)]
fwrite(introns, "example/gencode_hg19_all_introns.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

five_prime <- fread("gencode_hg19_fiveprime.bed.gz")
fwrite(five_prime, "example/gencode_hg19_fiveprime.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

three_prime <- fread("gencode_hg19_threeprime.bed.gz")
fwrite(three_prime, "example/gencode_hg19_threeprime.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# clear workspace
rm(list = ls())

### original file starts here ###

options(echo = TRUE)

# example input files
VCF <- "example/example.vcf"
clusters_table <- "example/cluster_counts.txt"
permutation_res <- "example/example_permutations.all.0.05.bh.txt"
# p values and betas
permutation_full_res <- "example/example_permutations.full.0.05.bh.txt"

# annotation - created by Leafcutter
annotation_code <- "example/gencode_hg19"
exon_file <- paste0(annotation_code, "_all_exons.bed")
all_introns <- paste0(annotation_code, "_all_introns.bed")
threeprime_file <- paste0( annotation_code, "_threeprime.bed")
fiveprime_file <- paste0( annotation_code, "_fiveprime.bed")

exons_table <- if (!is.null( exon_file )) {
  cat("Loading exons from", exon_file, "\n")
  as.data.frame(fread(exon_file))
} else {
  cat("No exon_file provided.\n")
  NULL
}

# read in cluster counts
clusters <- read.table(clusters_table, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# read in junction x snp results
print("reading in results")

res <- as.data.frame(fread(permutation_res, header = TRUE), stringsAsFactors = FALSE)
genotypes <- unique(res$dummy2)

######
# PREPARE VCF GENOTYPES
######

vcf <- read.table(file = VCF, header = TRUE, check.names = FALSE)

vcf[vcf$REF == TRUE]$REF <- "T"
vcf[vcf$ALT == TRUE]$ALT <- "T"

vcf_meta <- sqtlviztools::get_vcf_meta(vcf)

##################
# PREPARE CLUSTERS
##################

# from significant associations
sigClusters <- str_split_fixed(res[, 1], ":", 4)[, 4]

introns <- leafcutter::get_intron_meta(row.names(clusters))
keepClusters <- match(introns$clu, sigClusters)

# remove non-significant (or non-GWAS SNP-associated) clusters
introns <- introns[!is.na(keepClusters),]
clusters <- clusters[!is.na(keepClusters),]

# rearrange sample columns in clusters so they match the VCF
vcf <- cbind(vcf[, 1:9], vcf[, 10:ncol(vcf)][, colnames(vcf[, 10:ncol(vcf)]) %in% colnames(clusters)])
samples <- names(vcf)[10:ncol(vcf)]
clusters <- clusters[, samples]

introns_to_plot <- get_intron_meta(row.names(clusters))

# thin down junctions in clusters - remove low contributing junctions
juncProp <- function(cluster){
  cluster$prop <- cluster$meanCount / sum(cluster$meanCount)
  return(cluster)
}

splitClusters <- introns_to_plot %>%
  mutate(
    clu = factor(.$clu, levels = unique(.$clu)),
    meanCount = rowMeans(clusters)) %>%
  split(.$clu) %>%
  purrr::map_df(juncProp) %>%
  mutate(clu = as.character(.$clu))

# thinning out clusters - turn off for now - this needs to be worked on
thinClusters <- FALSE

if (thinClusters == TRUE) {
  introns_to_plot <- introns_to_plot[splitClusters$prop >= 0.01,]
  clusters <- clusters[splitClusters$prop >= 0.01,]
  introns <- introns[splitClusters$prop >= 0.01,]
}

####################
# ANNOTATE JUNCTIONS
####################

# functions are now in sqtlviztools
intersects <- sqtlviztools::intersect_introns(introns)

threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]

print("Annotating junctions")

uniqueClusters <- unique(introns$clu)

# for debugging
save.image("debug.RData")

annotatedClusters <- purrr::map_df(seq_along(uniqueClusters),
                                   ~annotate_single_cluster(introns,
                                                            clu = uniqueClusters[.],
                                                            cluIndex = .,
                                                            fiveprime = fiveprime_intersect,
                                                            threeprime = threeprime_intersect,
                                                            bothSS = all.introns_intersect
                                   )
)

annotatedClusters[is.na(annotatedClusters$gene),]$gene <- "."
annotatedClusters[is.na(annotatedClusters$ensemblID),]$ensemblID <- "."

#################
# PREPARE RESULTS - MOST SIGNIFICANT SNP x JUNCTION
#################

# Bind together metadata with original results
sigJunctions <- cbind(get_intron_meta(res[, 1]),
                      res[, c(1, 2, 3)],
                      get_snp_meta(res[, 2]))

# sometimes there will be duplicates - remove!
sigJunctions <- dplyr::distinct(sigJunctions)

#names(sigJunctions)[8] <- "bpval"
# present most significant junction for each SNP?
# or most significant SNP for each junction?
resultsByCluster <- dplyr::group_by(sigJunctions[order(sigJunctions$bpval),], clu) %>%
  dplyr::summarise(chr = first(chr),
                   start = min(start),
                   end = max(end),
                   snp = first(snp_ID),
                   snp_chr = first(snp_chr),
                   pos = first(snp_pos),
                   FDR = first(bpval)) %>%
  dplyr::arrange(FDR)

####
## PREPARE FOR SHINY
####

code <- "example"
annotation_code <- "gencode_hg19"

resultsByCluster$gene <- annotatedClusters$gene[match(resultsByCluster$clu, annotatedClusters$clusterID)]
resultsByCluster$SNP_pos <- paste0(resultsByCluster$snp_chr, ":", resultsByCluster$pos)

resultsByCluster$cluster_pos = paste0(resultsByCluster$chr, ":", resultsByCluster$start, "-", resultsByCluster$end)

resultsToPlot <- as.data.frame(select(resultsByCluster,
                                      SNP = snp,
                                      SNP_pos,
                                      gene = gene,
                                      cluster_pos,
                                      q = FDR
))

row.names(resultsToPlot) <- resultsByCluster$clu

resultsToPlot$q <- signif(resultsToPlot$q, digits = 3)

#########
## GET BETAS FOR EACH JUNCTION
#########

perm_full <- fread(permutation_full_res, data.table = FALSE)
names(perm_full) <- c("clusterID", "SNP","X","FDR", "Beta")

junctionsNeeded <- introns_to_plot %>%
  mutate(chr = gsub("chr", "", chr)) %>%
  mutate(cluster = paste(chr, start, end, clu, sep = ":")) %>%
  pull(cluster)

snps <- read.table(file = "example/snp_ids_rs_ids.txt", stringsAsFactors = FALSE, header = TRUE)

perm_full$RS_id <- snps$RS_id[ match(perm_full$SNP, snps$dummy2)]

perm_clean <- perm_full %>%
  filter(clusterID %in% junctionsNeeded) %>%
  filter(!is.na(RS_id)) %>%
  filter(clusterID %in% junctionsNeeded) %>%
  mutate(SNP_ID = gsub(":", ".", paste0(RS_id, ".", SNP)))

# sort out SNP IDs
perm_clean <- select(perm_clean, clusterID, SNP = SNP_ID, Beta, FDR)

perm_clean <- cbind(perm_clean,
                    get_intron_meta(perm_clean$clusterID),
                    get_snp_meta(perm_clean$SNP))

perm_clean$snp_chr <- paste0("chr", perm_clean$snp_chr)

# junction table - each junction with Beta, P value and annotation
junctionTable <- resultsToPlot %>%
  mutate(clu = row.names(resultsToPlot)) %>%
  left_join(introns_to_plot, by = "clu") %>%
  rename(snp_ID = SNP) %>%
  left_join(perm_clean,
            by = c("chr" = "snp_chr", "start", "end", "snp_ID", "clu", "middle")
  ) %>%
  mutate(coord = paste0(chr, ":", start, "-", end)) %>%
  left_join(annotatedClusters,
            by = c("clu" = "clusterID", "coord" )
  ) %>%
  select(clu, coord, verdict, Beta, q = FDR) %>%
  mutate(Beta = signif(Beta, digits = 3),
         q = signif(q, digits = 3)) %>%
  mutate(Beta = ifelse(is.na(Beta), ".", Beta),
         q = ifelse(is.na(q), ".", q))

#########
## SAVE OBJECTS
#########
rm(perm_full)

print("saving objects")
save(annotatedClusters, # every junction needed
     sigJunctions, # every junction x SNP interaction
     resultsToPlot, #significant clusters and the most significant SNP
     clusters, # junction counts for each sample
     vcf,# the genotypes of each sample
     vcf_meta, # the vcf metadata
     introns_to_plot, # all the intron positions
     exons_table, # the annotation
     junctionTable, # the junctions to display for each cluster
     annotation_code,
     code,
     file = "shiny/sQTL_results.Rdata"
)


