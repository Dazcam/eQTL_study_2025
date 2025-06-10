suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("susieR"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("arrow"))

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--phenotype_meta"), type="character", default=NULL,
                        help="Phenotype metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--sample_meta"), type="character", default=NULL,
                        help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--expression_matrix"), type="character", default=NULL,
                        help="Expression matrix file path with gene phenotype-id in rownames and sample-is in columnnames", metavar = "type"),
  optparse::make_option(c("--phenotype_list"), type="character", default=NULL,
                        help="Path to the phenotype list file.", metavar = "type"),
  optparse::make_option(c("--genotype_matrix"), type="character", default=NULL,
                        help="Genotype dosage matrix extracted from VCF.", metavar = "type"),
  optparse::make_option(c("--covariates"), type="character", default=NULL,
                        help="Path to covariates file in QTLtools format.", metavar = "type"),
  optparse::make_option(c("--out_prefix"), type="character", default="./finemapping_output",
                        help="Prefix of the output files.", metavar = "type"),
  optparse::make_option(c("--qtl_group"), type="character", default=NULL,
                        help="Value of the current qtl_group.", metavar = "type"),
  optparse::make_option(c("--cisdistance"), type="integer", default=1000000, 
                        help="Cis distance in bases from center of gene. [default \"%default\"]", metavar = "number"),
  optparse::make_option(c("--chunk"), type="character", default="1 1", 
                        help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type"),
#  optparse::make_option(c("--eqtlutils"), type="character", default=NULL,
#              help="Optional path to the eQTLUtils R package location. If not specified then eQTLUtils is assumed to be installed in the container. [default \"%default\"]", metavar = "type"),
  optparse::make_option(c("--write_full_susie"), type="character", default="true",
                        help="If 'true' then full SuSiE output will not be written to disk. Setting this to 'false' will apply credible set connected components based filtering to SuSiE results. [default \"%default\"]", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Print all options
print(opt)

#### Local testing  -----
# test_data <- FALSE  # Set to FALSE for real data
# 
# if (test_data) {
#   root <- "~/Desktop/test/susie_test"
#   exp_mat     <- file.path(root, "GEUVADIS_cqn.tsv")
#   gene_meta   <- file.path(root, "GEUVADIS_phenotype_metadata.tsv")
#   smpl_lst    <- file.path(root, "GEUVADIS_sample_metadata.tsv")
#   eqtl_to_test <- file.path(root, "GEUVADIS_test_ge.permuted.tsv.gz")
#   covar       <- file.path(root, "GEUVADIS_test_ge.covariates.txt")
#   geno_mat    <- file.path(root, "LCL.dose.tsv.gz")
# } else {
#   root <- "~/Desktop/test/susie_test"
#   exp_mat     <- file.path(root, "ExN-1_tmm_susie_ready.bed")
#   gene_meta   <- file.path(root, "ExN-1_gene_meta.tsv")
#   smpl_lst    <- file.path(root, "ExN-1_samp_lst_susie_ready.txt")
#   eqtl_to_test <- file.path(root, "ExN-1_perm.cis_qtl.txt.gz")
#   covar       <- file.path(root, "ExN-1_covariates_susie_ready.txt")
#   geno_mat    <- file.path(root, "chrALL_final.filt.dose.tsv.gz")
# }

#### Helper Functions  ----------
importQtlmapCovariates <- function(covariates_path){
  pc_matrix = read.table(covariates_path, check.names = F, header = T, stringsAsFactors = F)
  
  # Added this. Testdata does not have chracter qualitative covariates
  pc_matrix[pc_matrix == "M"] <- 0
  pc_matrix[pc_matrix == "F"] <- 1
  # Error of colinear cov matrix for male/female eQTL. All sex is 0/1
  pc_matrix <- pc_matrix[apply(pc_matrix[-1], 1, var) != 0, ]
  pc_transpose = t(pc_matrix[,-1])
  colnames(pc_transpose) = pc_matrix$id
  pc_df = dplyr::mutate(as.data.frame(pc_transpose), genotype_id = rownames(pc_transpose)) %>%
    dplyr::as_tibble() %>% 
    dplyr::select(genotype_id, dplyr::everything())
  
  #Make PCA matrix
  pc_matrix = as.matrix(dplyr::select(pc_df,-genotype_id))
  rownames(pc_matrix) = pc_df$genotype_id
  return(pc_matrix)
}


# importQtlmapPermutedPvalues <- function(perm_path){
#   tbl = read.table(perm_path, check.names = F, header = T, stringsAsFactors = F) %>%
#     dplyr::as_tibble() %>%
#     dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
#     dplyr::mutate(group_id = molecular_trait_object_id)
#   return(tbl)
# }


splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}


splitIntoChunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = floor(n_total/(n_chunks))
  batches = splitIntoBatches(n_total,chunk_size)
  batches[batches > n_chunks] = n_chunks
  selected_batch = batches == chunk_number
  return(selected_batch)
}

finemapPhenotype <- function(phenotype_id, se, genotype_file, covariates, cis_distance){
  message("Processing phenotype: ", phenotype_id)
  
  #Extract phenotype from SE
  gene_vector = eQTLUtils::extractPhentypeFromSE(phenotype_id, se, "counts") %>%
    dplyr::mutate(phenotype_value_std = qnorm((rank(phenotype_value, na.last = "keep") - 0.5) / sum(!is.na(phenotype_value))))
  selected_phenotype = phenotype_id
  gene_meta = dplyr::filter(SummarizedExperiment::rowData(se) %>% as.data.frame(), phenotype_id == selected_phenotype)
  
  #Rearrange samples in the covariates matrix
  covariates_matrix = cbind(covariates[gene_vector$genotype_id,], 1)
  
  #Import genotype matrix
  genotype_matrix = eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = gene_meta$chromosome, 
    start = gene_meta$phenotype_pos - cis_distance, 
    end = gene_meta$phenotype_pos + cis_distance, 
    dosage_file = genotype_file)
  
  #Residualise gene expression and genotype matrix
  hat = diag(nrow(covariates_matrix)) - covariates_matrix %*% solve(crossprod(covariates_matrix)) %*% t(covariates_matrix)
  expression_vector = hat %*% gene_vector$phenotype_value_std
  names(expression_vector) = gene_vector$genotype_id
  
  gt_matrix = genotype_matrix[,names(expression_vector)]
  
  #Exclude variants with no alternative alleles
  gt_matrix = gt_matrix[rowSums(round(gt_matrix,0)) != 0,]
  
  #Standardise genotypes
  gt_std = t(gt_matrix - apply(gt_matrix, 1, mean))
  gt_hat = hat %*% gt_std
  
  # Fit finemapping model
  fitted <- susieR::susie(gt_hat, expression_vector,
                          L = 10,
                          estimate_residual_variance = TRUE, 
                          estimate_prior_variance = TRUE,
                          scaled_prior_variance = 0.1,
                          verbose = TRUE,
                          compute_univariate_zscore = TRUE,
                          min_abs_corr = 0)
  fitted$variant_id = rownames(gt_matrix)
  return(fitted)
}

# Reported issue for susieR:susie_get_posterior_mean/sd  when include_idx is a single number, i.e. only one causal variant detected
susie_get_posterior_mean_custom = function (res, prior_tol = 1e-9) {
  
  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V))
    include_idx = which(res$V > prior_tol)
  else
    include_idx = 1:nrow(res$alpha)
  
  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0)
    return(colSums((res$alpha*res$mu)[include_idx,,drop=FALSE])/
             res$X_column_scale_factors)
  else
    return(numeric(ncol(res$mu)))
}
susie_get_posterior_sd_custom = function (res, prior_tol = 1e-9) {
  
  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V))
    include_idx = which(res$V > prior_tol)
  else
    include_idx = 1:nrow(res$alpha)
  
  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0)
    return(sqrt(colSums((res$alpha * res$mu2 -
                           (res$alpha*res$mu)^2)[include_idx,,drop=FALSE]))/
             (res$X_column_scale_factors))
  else
    return(numeric(ncol(res$mu)))
}

extractResults <- function(susie_object){
  credible_sets = susie_object$sets$cs
  cs_list = list()
  susie_object$sets$purity = dplyr::as_tibble(susie_object$sets$purity) %>%
    dplyr::mutate(
      cs_id = rownames(susie_object$sets$purity),
      cs_size = NA,
      cs_log10bf = NA,
      overlapped = NA
    )
  added_variants = c()
  # overlap checks if CS variants already exist in previous CS
  for (index in seq_along(credible_sets)){
    cs_variants = credible_sets[[index]]
    cs_id = susie_object$sets$cs_index[[index]]
    
    is_overlapped = any(cs_variants %in% added_variants)
    susie_object$sets$purity$overlapped[index] = is_overlapped
    susie_object$sets$purity$cs_size[index] = length(cs_variants)
    susie_object$sets$purity$cs_log10bf[index] = log10(exp(susie_object$lbf[cs_id]))
    if (!is_overlapped) {
      cs_list[[index]] = dplyr::tibble(cs_id = paste0("L", cs_id),
                                       variant_id = susie_object$variant_id[cs_variants])
      added_variants = append(added_variants, cs_variants)
    }
  }
  df = purrr::map_df(cs_list, identity)
  
  #Extract purity values for all sets
  purity_res = susie_object$sets$purity
  
  #Sometimes all the PIP values are 0 and there are no purity values, then skip this step
  if(nrow(purity_res) > 0){
    purity_df = dplyr::as_tibble(purity_res) %>%
      dplyr::filter(!overlapped) %>%
      dplyr::mutate(
        cs_avg_r2 = mean.abs.corr^2,
        cs_min_r2 = min.abs.corr^2,
        low_purity = min.abs.corr < 0.5
      )  %>%
      dplyr::select(cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity) 
  } else{
    purity_df = dplyr::tibble()
  }
  
  #Extract betas and standard errors
  # mean_vec = susieR::susie_get_posterior_mean(susie_object)
  mean_vec = susie_get_posterior_mean_custom(susie_object)
  # sd_vec = susieR::susie_get_posterior_sd(susie_object)
  sd_vec = susie_get_posterior_sd_custom(susie_object)
  alpha_mat = t(susie_object$alpha)
  colnames(alpha_mat) = paste0("alpha", seq(ncol(alpha_mat)))
  mean_mat = t(susie_object$alpha * susie_object$mu) / susie_object$X_column_scale_factors
  colnames(mean_mat) = paste0("mean", seq(ncol(mean_mat)))
  sd_mat = sqrt(t(susie_object$alpha * susie_object$mu2 - (susie_object$alpha * susie_object$mu)^2)) / (susie_object$X_column_scale_factors)
  colnames(sd_mat) = paste0("sd", seq(ncol(sd_mat)))
  posterior_df = dplyr::tibble(variant_id = names(mean_vec), 
                               pip = susie_object$pip,
                               z = susie_object$z,
                               posterior_mean = mean_vec, 
                               posterior_sd = sd_vec) %>%
    dplyr::bind_cols(purrr::map(list(alpha_mat, mean_mat, sd_mat), dplyr::as_tibble))
  
  if(nrow(df) > 0 & nrow(purity_df) > 0){
    cs_df = purity_df
    variant_df = dplyr::left_join(posterior_df, df, by = "variant_id") %>%
      dplyr::left_join(cs_df, by = "cs_id")
  } else{
    cs_df = NULL
    variant_df = NULL
  }
  
  return(list(cs_df = cs_df, variant_df = variant_df))
}


#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = utils::read.csv(opt$phenotype_meta, sep = "\t", stringsAsFactors = F)
covariates_matrix = importQtlmapCovariates(opt$covariates)

#Exclude covariates with zero variance
#exclude_cov = apply(covariates_matrix, 2, sd) != 0
#covariates_matrix = covariates_matrix[,exclude_cov]

#Import list of phenotypes for finemapping
#phenotype_table = importQtlmapPermutedPvalues(opt$phenotype_list)
phenotype_table = read_tsv(opt$phenotype_list) # TensorQTL automatically calcs qvals
#filtered_list = dplyr::filter(phenotype_table, p_fdr < 0.05) # Changed from 0.01 to 0.05 (should probs add to shell script as param)
filtered_list = dplyr::filter(phenotype_table, qval < 0.05)

#phenotype_list = dplyr::semi_join(phenotype_meta, filtered_list, by = "group_id")
#phenotype_list = dplyr::semi_join(phenotype_meta, filtered_list, by = join_by("group_id" == "phenotype_id")) # TensorQTL output different
phenotype_list = dplyr::semi_join(phenotype_meta, filtered_list, by = c("group_id" = "phenotype_id")) # No join_by() in dplyr 1.0.8 

message("Number of phenotypes included for analysis: ", nrow(phenotype_list))

#Keep only those phenotypes that are present in the expression matrix
phenotype_list = dplyr::filter(phenotype_list, phenotype_id %in% expression_matrix$phenotype_id)


#### I added this #### 
# Track cell type rather than qtl_group or study_id
cell_type <- str_extract(basename(opt$expression_matrix), "^[^_]+")
message("cell_type is set to: ", cell_type)

#eQTLUtils::makeSummarizedExperimentFromCountMatrix doesn't work unless the additional cols are added
# NAs for now may add metadata later
# What is genotype_id? Set it to sample ID as genotype is used later
# Error in covariates[gene_vector$genotype_id, ] : subscript out of bounds
# Calls: <Anonymous> ... finemapPhenotype -> cbind -> standardGeneric -> eval -> eval -> eval
if (!is.data.frame(sample_metadata)) {
  sample_metadata <- data.frame(sample_id = sample_metadata, stringsAsFactors = FALSE)
}
sample_metadata$genotype_id <- sample_metadata$sample_id  
sample_metadata$qtl_group <- cell_type  # Match test data
sample_metadata$sex <- NA_character_
str(sample_metadata)
####


#Set parameters
cis_distance = opt$cisdistance
genotype_file = opt$genotype_matrix
#study_id = sample_metadata$study[1] # Track cell type as study_id
study_id = cell_type

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                         row_data = phenotype_meta, 
                                                         col_data = sample_metadata, 
                                                         quant_method = "gene_counts",
                                                         reformat = FALSE)

#If qtl_group is not specified, then use the first value in the qtl_group column of the sample metadata
if(is.null(opt$qtl_group)){
  opt$qtl_group = se$qtl_group[1]
}

#Split phenotype list into chunks
chunk_vector = strsplit(opt$chunk, split = " ") %>% unlist() %>% as.numeric()
chunk_id = chunk_vector[1]
n_chunks = chunk_vector[2]
selected_chunk = splitIntoChunks(chunk_id, n_chunks, length(phenotype_list$phenotype_id))
selected_phenotypes = phenotype_list$phenotype_id[selected_chunk] %>% setNames(as.list(.), .)

#Only proceed if the there are more than 0 phenotypes
message("Number of overall unique group_ids: ", length(unique(phenotype_list$group_id)))
message("Number of phenotypes in the batch: ", length(selected_phenotypes))

#Check that the qtl_group is valid and subset
assertthat::assert_that(opt$qtl_group %in% unique(se$qtl_group))
selected_qtl_group = eQTLUtils::subsetSEByColumnValue(se, "qtl_group", 'ExN-1')

#Apply finemapping to all genes
results = purrr::map(selected_phenotypes, ~finemapPhenotype(., selected_qtl_group, 
                                                            genotype_file, covariates_matrix, cis_distance))

#Define fine-mapped regions
region_df = dplyr::transmute(phenotype_list, phenotype_id, finemapped_region = paste0("chr", chromosome, ":", phenotype_pos - cis_distance, "-", phenotype_pos + cis_distance))
# message("region df OK!")

#Extract credible sets from finemapping results
res = purrr::map(results, extractResults) %>%
  purrr::transpose()
# message("res OK!")

#Extract information about all variants
variant_df <- purrr::map_df(res$variant_df, identity, .id = "phenotype_id")
if(nrow(variant_df) > 0){
  variant_df <- variant_df %>%
    dplyr::left_join(region_df, by = "phenotype_id") %>%
    tidyr::separate(variant_id, c("chr", "pos", "ref", "alt"),sep = "_", remove = FALSE) %>%
    dplyr::mutate(chr = stringr::str_remove_all(chr, "chr")) %>%
    dplyr::mutate(cs_index = cs_id) %>%
    dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_"))
}
# message("variant df OK!")

#Extraxt information about credible sets
cs_df <- purrr::map_df(res$cs_df, identity, .id = "phenotype_id")
# message("cs df OK!")

if(nrow(cs_df) > 0){
  cs_df = dplyr::left_join(cs_df, region_df, by = "phenotype_id") %>%
    dplyr::mutate(cs_index = cs_id) %>%
    dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_")) %>%
    dplyr::select(phenotype_id, cs_id, cs_index, finemapped_region, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity)
  
  #Extract information about variants that belong to a credible set
  in_cs_variant_df <- dplyr::filter(variant_df, !is.na(cs_index) & !low_purity) %>%
    dplyr::select(phenotype_id, variant_id, chr, pos, ref, alt, cs_id, cs_index, finemapped_region, pip, z, cs_min_r2, cs_avg_r2, cs_size, posterior_mean, posterior_sd, cs_log10bf)
} else{
  #Initialize empty tibbles with correct column names
  in_cs_variant_df = dplyr::tibble(
    phenotype_id = character(),
    variant_id = character(),
    chr = numeric(),
    pos = numeric(),
    ref = numeric(),
    alt = numeric(),
    cs_id = numeric(),
    cs_index = numeric(),
    finemapped_region = numeric(),
    pip = numeric(),
    z = numeric(),
    cs_min_r2 = numeric(),
    cs_avg_r2 = numeric(),
    cs_size = numeric(),
    posterior_mean = numeric(),
    posterior_sd = numeric(),
    cs_log10bf = numeric()
  )
  
  cs_df = dplyr::tibble(
    phenotype_id = numeric(),
    cs_id = numeric(),
    cs_index = numeric(),
    finemapped_region = numeric(),
    cs_log10bf = numeric(),
    cs_avg_r2 = numeric(),
    cs_min_r2 = numeric(),
    cs_size = numeric(),
    low_purity = numeric()
  )
}

#Extract information about all variants
if(nrow(variant_df) > 0){
  variant_df <- dplyr::select(variant_df, phenotype_id, variant_id, chr, pos, ref, alt, cs_id, cs_index, low_purity, finemapped_region, pip, z, posterior_mean, posterior_sd, alpha1:sd10)
}

#Export high purity credible set results only
write.table(in_cs_variant_df, paste0(opt$out_prefix, ".cred.hp.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Export all other results
write.table(cs_df, paste0(opt$out_prefix, ".cred.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(variant_df, paste0(opt$out_prefix, ".snp.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



