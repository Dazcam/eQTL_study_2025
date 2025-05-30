#--------------------------------------------------------------------------------------
#
#    Plot TSS densities, and eQTL / snATAC peak overlaps
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  PMID: Ziffra 34616060; O'Brien 30419947; Bryois 35915177
#  Running locally for now

## Consideration for eQTL overlap

# Method 1: SNP-eGene Pair Match: stringent definition of overlap
# Method 2: eGene Match with Linkage Disequilibrium (LD) Consideration

##  -----------------------------------------------------------------------------------

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(stringr)
library(readr)
library(cowplot)
library(readxl)
library(LDlinkR)

root_dir <- '~/Desktop/eQTL_study_2025/'
tensorqtl_dir <- paste0(root_dir, 'results/04TENSORQTL/')
fig_dir <- paste0(root_dir, 'results/figs/')
table_dir <- paste0(root_dir, 'results/tables/')
ziffra_dir <- paste0(root_dir, '/resources/public_datasets/ziffra_2021/')
obrien_dir <- paste0(root_dir, '/resources/public_datasets/obrien_2018/')
bryois_dir <- paste0(root_dir, '/resources/public_datasets/bryois_2022/')
ziffra_url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03209-8/MediaObjects/41586_2021_3209_MOESM4_ESM.xlsx'
obrien_url <- 'https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-018-1567-1/MediaObjects/13059_2018_1567_MOESM1_ESM.xlsx'
bryois_url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-022-01128-z/MediaObjects/41593_2022_1128_MOESM3_ESM.xlsx'
ziffra_dest_file <- paste0(ziffra_dir, "Ziffra_2021_supp_tables_2_13.xlsx")
obrien_dest_file <- paste0(obrien_dir, "Obrien_2018_supp_tables_1_10.xlsx")
bryois_dest_file <- paste0(obrien_dir, "Bryois_2022_supp_tables_1_7.xlsx")
cell_types <- c('RG', 'ExN-1', 'InN-1', 'ExN-2', 'InN-2', 'Endo-Peri', 'MG', 'Mig-N')
dir.create(ziffra_dir, recursive = TRUE)
dir.create(obrien_dir, recursive = TRUE)
dir.create(bryois_dir, recursive = TRUE)

## Load public data files
for (dest_file in c(ziffra_dest_file, obrien_dest_file, bryois_dest_file)) {
  
  url <- case_when(
    str_detect(dest_file, 'Ziffra') ~ ziffra_url,
    str_detect(dest_file, 'Obrien') ~ obrien_url,
    str_detect(dest_file, 'Bryois') ~ bryois_url)
  
  if (!file.exists(dest_file)) {
    download.file(url, destfile = dest_file, mode = "wb")
    message("File downloaded to: ", dest_file)
    } else {
    message("File already exists: ", dest_file)
    }
}


##  Load Ziffra peaks  -----------------------------------------------------------------
# peak list = all peaks key; mac2 all peaks per cell type; specific; specific peaks per cell type
cat('\nLoading Ziffra 2021 snATAC-seq fetal brain data ... \n')
options(scipen = 999) # required to prevent peak coords. being abbr. in sci' notation
peak_list <- read_excel(paste0(ziffra_dir, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST2 AllPrimaryPeaks') %>%
  dplyr::select(seqnames, start, end, peak_name) 
enhancer_peaks <- read_excel(paste0(ziffra_dir, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST5 EnhancerPeaks') %>%
  dplyr::select(peak_name) 
lookup_peaks <- read_excel(paste0(ziffra_dir, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST3 MACSpeaks_byCelltype')

# Convert peak_list to GRanges
peaks_gr <- makeGRangesFromDataFrame(
  peak_list,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

# cat('\n\nAnnotating peaks with SNP ... \n\n\n')
# snps_withPeak_ann <- ChIPseeker::annotatePeak(snps_withPeak_gr, tssRegion = c(-1000, 100),
#                                               TxDb = txdb, annoDb = "org.Hs.eg.db", level = 'gene')

cell_types <- c('Endo-Peri', 'MG', 'Mig-N')
overlap_list <- list()
overlap_report_tbl <- tibble()
for (cell_type in cell_types) {
  message("Processing eQTL data for: ", cell_type)
  
  # Load eQTL data
  eqtl_tbl <- read_tsv(paste0(tensorqtl_dir, cell_type, "/", cell_type, "_perm.cis_qtl.txt.gz")) %>%
    filter(!is.na(start_distance), qval < 0.05)
  
  if (nrow(eqtl_tbl) == 0) {
    message("No significant eQTLs found for: ", cell_type)
    next
  }
  
  # Get SNP coordinates from Ensembl
  ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  snp_lookup_hg38 <- getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
    filters = "snp_filter",
    values = eqtl_tbl$variant_id,
    mart = ensembl,
    uniqueRows = TRUE
  )
  
  snp_lookup_hg38_filt <- snp_lookup_hg38 |>
    filter(chr_name %in% c(1:22, "X", "Y")) |>
    mutate(chr_name = paste0("chr", chr_name)) |>
    distinct()
  
  message(nrow(snp_lookup_hg38) - nrow(snp_lookup_hg38_filt), 
          ' eQTL removed due to non-standard chr annotation.')

  # Convert eQTL SNPs to GRanges
  snp_gr <- makeGRangesFromDataFrame(
    snp_lookup_hg38_filt,
    seqnames.field = "chr_name",
    start.field = "chrom_start",
    end.field = "chrom_end",
    keep.extra.columns = TRUE
  )
  
  # Find overlaps
  overlaps <- findOverlaps(snp_gr, peaks_gr)
  
  # Extract overlapping SNPs and peaks
  if (length(overlaps) > 0) {
    
    overlap_tbl <- data.frame(
      
      cell_type = cell_type,
      variant_id = snp_lookup_hg38$refsnp_id[queryHits(overlaps)],
      chr = seqnames(snp_gr)[queryHits(overlaps)],
      snp_start = start(snp_gr)[queryHits(overlaps)],
      peak_name = peak_list$peak_name[subjectHits(overlaps)],
      peak_chr = seqnames(peaks_gr)[subjectHits(overlaps)],
      peak_start = start(peaks_gr)[subjectHits(overlaps)],
      peak_end = end(peaks_gr)[subjectHits(overlaps)]
      ) |>
      as_tibble() |>
      mutate(enhancer_peak = ifelse(peak_name %in% enhancer_peaks$peak_name, TRUE, FALSE))
    
    # Merge with eQTL data to include additional info
    overlap_tbl <- overlap_tbl %>%
      left_join(eqtl_tbl %>% dplyr::select(variant_id, phenotype_id, qval, start_distance), 
                by = "variant_id")
    
    overlap_list[[cell_type]] <- overlap_tbl
    
    overlap_report_tbl <- bind_rows(overlap_report_tbl, tibble(
      cell_type = cell_type,
      sc_eqtl = nrow(eqtl_tbl),
      sc_eqtl_filt = nrow(snp_lookup_hg38_filt),
      ziffra_peaks = nrow(peak_list),
      num_overlaps = nrow(overlap_tbl),
      num_enhancers = sum(overlap_tbl$enhancer_peak)
    ))
    
    message("Found ", nrow(overlap_tbl), " eQTL-peak overlaps for: ", cell_type)
  } else {
    message("No overlaps found for: ", cell_type)
  }
  
}
overlap_list

## Overlapping with O-Brien Bulk eQTL  ------------------------------------------------
message("\nLoading O'Brien (2018) bulk fetal brain eQTL data ... \n")
options(scipen = 999) # required to prevent peak coords. being abbr. in sci' notation
obrien_eqtl_tbl <- read_excel(obrien_dest_file, sheet = 'TblS2', skip = 2) 

## Method 1
message("\nRunning overlap method 1 on O'Brien (2018) bulk eQTL data ... \n")
bulk_overlap1_tbl <- tibble()
obrien_list <- list()
slope_tbl <- tibble()
for (cell_type in cell_types) {
  
  message('Reading eQTL data for: ', cell_type)
  eqtl_tbl <- read_tsv(paste0(tensorqtl_dir, cell_type, '/', cell_type, '_perm.cis_qtl.txt.gz')) |>
    filter(qval < 0.05) |>
    dplyr::select(rsID = variant_id, ensembl_id_sc = phenotype_id, slope_sc = slope, qval_sc = qval)
  
  message('eQTL detected at FDR < 0.05: ', nrow(eqtl_tbl))
  
  bulk_overlap <- obrien_eqtl_tbl |>
    inner_join(eqtl_tbl, by = join_by(rsID)) |>
    mutate(gene_match = ifelse(ENSEMBL_geneID == ensembl_id_sc, T, F)) |>
    mutate(direction_match = sign(slope) == sign(slope_sc)) |>
    dplyr::select(Chromosome, Position, rsID, ENSEMBL_geneID, ensembl_id_sc, 
                  SYMBOL, gene_match, direction_match, slope, slope_sc, qval, qval_sc)
    message('Overlapping snp: ', nrow(bulk_overlap))
    message('Overlapping snp and gene: ', sum(bulk_overlap$gene_match))
    message('Direction of effect match: ', sum(bulk_overlap$direction_match))

    bulk_overlap1_tbl <- bind_rows(bulk_overlap1_tbl, tibble(
      cell_type = cell_type,
      sc_eqtl = nrow(eqtl_tbl),
      bulk_eqtl = nrow(obrien_eqtl_tbl),
      overlapping_snp = nrow(bulk_overlap),
      overlapping_snp_and_gene = sum(bulk_overlap$gene_match),
      direction_effect_match = sum(bulk_overlap$direction_match)
    ))
    
    # Check slope
    slope_tbl <- bind_rows(slope_tbl, bulk_overlap |>
      summarise(slope_min = min(slope_sc),
                slope_max = max(slope_sc),
                slope_range = max(slope_sc, na.rm = TRUE) - min(slope_sc, na.rm = TRUE)) |>
      mutate(cell_type = cell_type) |>
      relocate(cell_type))
    
    obrien_list[[cell_type]] <- bulk_overlap
    
  
}
bulk_overlap1_tbl
slope_tbl

## Method 2
message("\nRunning overlap method 2 on O'Brien (2018) bulk eQTL data ... \n")
# Set up LDlink API token (get from LDlink website)
token <- "your_ldlink_token" # Replace with your actual token



## Create TSS density plot  -----------------------------------------------------------
plt_list <- list()
for (cell_type in cell_types) {
  
  message('Reading eQTL data for: ', cell_type)
  eqtl_tbl <- read_tsv(paste0(tensorqtl_dir, cell_type, '/', cell_type, '_perm.cis_qtl.txt.gz')) %>%
    filter(!is.na(start_distance)) %>%
    mutate(distance_kb = start_distance / 1000) %>%
    filter(qval < 0.05)
  
  message('eQTL detected at FDR < 0.05: ', nrow(eqtl_tbl))

  if (nrow(eqtl_tbl) > 0) {
    
    message('Plotting.')
    
    eqtl_count <- nrow(eqtl_tbl)
    
    eqtl_density_plt <- eqtl_tbl  %>%
      ggplot(aes(x = distance_kb)) +
      geom_density() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      annotate("text", x = 100, y = Inf, 
               label = paste("eQTLs:", eqtl_count), 
               hjust = 0, vjust = 1.5, size = 4, fontface = "bold") +
      scale_x_continuous(limits = c(-600, 600), 
                         breaks = seq(-600, 600, by = 200),
                         labels = seq(-600, 600, by = 200)) +
      labs(x = "Distance to TSS (kb)",
           y = "Density",
           title = cell_type) +
      theme_minimal() +
      theme(
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold")
      )
    
    plt_list[[cell_type]] <- eqtl_density_plt
    
  } else {
    
    message('Skipped plotting.')
    
  }

}

## Overlapping with Bryois adult sn-eQTL  ------------------------------------------------
message("\nLoading Bryois (2022) adult brain sn-eQTL data ... \n")
options(scipen = 999) # required to prevent peak coords. being abbr. in sci' notation
bryois_eqtl_tbl <- read_excel(bryois_dest_file, sheet = 'Table S2', skip = 3) |>
  dplyr::rename(rsID = SNP,
                adult_cell_type = cell_type) |>
  filter(adj_p < 0.05)

bryois_list <- list()
bryois_report1_tbl <- tibble()
for (cell_type in cell_types) {
  
  message('Reading eQTL data for: ', cell_type)
  eqtl_tbl <- read_tsv(paste0(tensorqtl_dir, cell_type, '/', cell_type, '_perm.cis_qtl.txt.gz')) |>
    filter(qval < 0.05) |>
    dplyr::select(rsID = variant_id, ensembl_id_sc = phenotype_id, slope_sc = slope, qval_sc = qval)
  
  message('eQTL detected at FDR < 0.05: ', nrow(eqtl_tbl))
  
  bryois_overlap <- bryois_eqtl_tbl |>
    inner_join(eqtl_tbl, by = join_by(rsID)) |>
    mutate(gene_match = ifelse(ensembl == ensembl_id_sc, T, F)) |>
    mutate(direction_match = sign(beta) == sign(slope_sc)) |>
    dplyr::select(adult_cell_type, rsID, symbol, ensembl, ensembl_id_sc, 
                  gene_match, direction_match, beta, slope_sc, adj_p, qval_sc)
  message('Overlapping snp: ', nrow(bryois_overlap))
  message('Overlapping snp and gene: ', sum(bryois_overlap$gene_match))
  message('Direction of effect match: ', sum(bryois_overlap$direction_match))
  
  bryois_report1_tbl <- bind_rows(bryois_report1_tbl, tibble(
    cell_type = cell_type,
    sc_eqtl = nrow(eqtl_tbl),
    bryois_eqtl = nrow(bryois_eqtl_tbl),
    overlapping_snp = nrow(bryois_overlap),
    overlapping_snp_and_gene = sum(bryois_overlap$gene_match),
    direction_effect_match = sum(bryois_overlap$direction_match)
  ))

  bryois_list[[cell_type]] <- bryois_overlap

}


openxlsx::write.xlsx(overlap_list, paste0(table_dir, 'eQTL_ziffra_peak_overlaps.xlsx'))
openxlsx::write.xlsx(bryois_list, paste0(table_dir, 'eQTL_bryois_overlaps.xlsx'))
openxlsx::write.xlsx(obrien_list, paste0(table_dir, 'eQTL_obrien_overlaps.xlsx'))
eqtl_plts <- cowplot::plot_grid(plotlist = plt_list)
ggsave(paste0(fig_dir, 'eqtl_plts.png'), eqtl_plts, width = 12, height = 12, 
       dpi = 300, units = "in", bg = "white")




#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
