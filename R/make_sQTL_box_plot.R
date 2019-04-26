#' Make genotype x junction count box plots
#'
#' @import ggplot2
#' @import ggbeeswarm
#' @import dplyr
#' @export

saveSwarmPlot <- function(myPlot) {
  pdf("~/Desktop/swarm_plot.pdf", height = 5, width = 3)
  print(myPlot)
  dev.off()
}

make_sQTL_box_plot <- function(
  cluster_to_plot,
  all_junctions = NA,
  junction_to_plot = NA,
  main_title = NA,
  exons_table = NULL,
  vcf = NULL,
  vcf_meta = NULL,
  cluster_ids = NULL,
  counts = NULL,
  introns = NULL,
  snp_pos=NA,
  junctionTable = NA,
  snp = snp ){

  stopifnot(snp %in% vcf_meta$SNP)

  if( !grepl("^chr", junction_to_plot)){
    junction_to_plot <- paste0("chr", junction_to_plot)
  }

  # subset VCF and get genotype groups
  vcfIndex <- which( vcf_meta$SNP == snp)
  VCF <- vcf[vcfIndex,10:ncol(vcf)]

  VCF_meta <- vcf_meta[vcfIndex,]
  meta <- as.data.frame(t(VCF))

  print(VCF_meta)
  meta$group=as.factor(meta[,1])
  group_names <- c(0,1,2)
  names(group_names) <- c(
    paste0( VCF_meta$REF, "/", VCF_meta$REF),
    paste0( VCF_meta$REF, "/", VCF_meta$ALT),
    paste0( VCF_meta$ALT, "/", VCF_meta$ALT)
  )

  y <- t(counts[ cluster_ids==cluster_to_plot, ])
  # for each sample divide each junction count by the total for that sample
  normalisedCounts <- as.data.frame(sweep(y, 1, rowSums(y), "/"))
  genotypes <- as.data.frame(t(VCF))
  names(genotypes)[1] <- "geno"
  normalisedCounts$genotypeCode <- genotypes$geno[ match( row.names(normalisedCounts), row.names(genotypes))]
  normalisedCounts <- normalisedCounts[ complete.cases(normalisedCounts),]

  normalisedCounts$genotype <- names(group_names)[ match(normalisedCounts$genotypeCode, group_names)]

  print(head(normalisedCounts))

  toPlot <- dplyr::select( normalisedCounts,
                    junction = junction_to_plot,
                    geno =  "genotype")

  toPlot$geno <- factor(toPlot$geno, levels = (names(group_names))) # this was reversed

  # get junction information for title
  junc <- dplyr::mutate(junctionTable,
                 j = paste0( gsub("-",":", coord),":", clu)
                 ) %>%
    dplyr::filter( j == junction_to_plot )

  values <- paste( signif(as.numeric(junc$Beta),3)," p =", signif(as.numeric(junc$q), 3) )

  # Function to produce summary statistics (mean and +/- sd)
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m - sd(x)
    ymax <- m + sd(x)
    return(c(y = m, ymin = ymin, ymax = ymax))
  }
  
  # toPlot$geno <- factor(toPlot$geno, rev(levels(toPlot$geno)))
  
  plot <- ggplot(data = toPlot, aes(x = geno, y = junction, group = geno, color = geno)) +
    geom_quasirandom(size = 1) +
    stat_summary(fun.data = data_summary, color = "black") +
    theme_classic() +
    theme(axis.title.y = element_text(angle = 90), legend.position = "none", axis.text = element_text(color = "black")) +
    ylab("relative intron representation") +
    xlab("") +
    labs(title = junc$coord,
         subtitle = bquote(beta == .(values))) + 
    scale_color_brewer(palette = "Dark2") +
    # coord_flip() +
    NULL
  
  saveSwarmPlot(plot)
  
  return(plot)
}
