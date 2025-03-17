#!/usr/bin/env /nfs/sw/R/R-4.2.3/lib64/R/bin/Rscript

require(hdf5r)
require(dplyr)
require(purrr)
require(inspre)

ncores=16
# k562_norm_sc_fn = '/gpfs/commons/groups/knowles_lab/gwps/data/K562_essential_normalized_singlecell_01.h5ad'
# save_dag_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_res_dag.Rdata'
# save_nodag_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_res_nodag.Rdata'
# save_guide_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_guide_data.Rdata'
# k562_norm_sc_fn = '/gpfs/commons/groups/knowles_lab/gwps/data/K562_gwps_normalized_singlecell_01.h5ad'
# save_dag_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_all_res_dag.Rdata'
# save_nodag_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_all_res_nodag.Rdata'
# save_guide_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_all_guide_data.Rdata'
k562_norm_sc_fn = '/gpfs/commons/groups/knowles_lab/gwps/data/rpe1_normalized_singlecell_01.h5ad'
save_dag_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_rpe1_res_dag.Rdata'
save_nodag_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_rpe1_res_nodag.Rdata'
save_guide_res = '/gpfs/commons/groups/knowles_lab/gwps/saved_rdata/gwps_rpe1_guide_data.Rdata'


# At some point this gets auto-converted to 'non.targeting' by R...
ntc = 'non-targeting'

# Load files.
process_start <- Sys.time()
cat("Loading files.\n")
hfile_sc <- H5File$new(k562_norm_sc_fn, "r")
obs <- inspre::parse_hdf5_df(hfile_sc, 'obs')
var <- inspre::parse_hdf5_df(hfile_sc, 'var')

# Load NTC into memory.
cat("Finding and loading NTC guides.\n")
genes_guides <- filter(obs, gene != ntc) %>% distinct(gene_id, gene_transcript) %>%
  filter(gene_id %in% var$gene_id)
cells_ntc <- obs$gene == ntc
n_ntc <- sum(cells_ntc)
X_ntc <- hfile_sc[['X']][, cells_ntc]
rownames(X_ntc) <- var$gene_id
colnames(X_ntc) <- obs$gene_transcript[cells_ntc]

# Calculate guide effects.
cat("Calculating guide effect sizes and selecting effective guides.\n")
guide_effects <- map2(genes_guides$gene_transcript, genes_guides$gene_id,
                      ~calc_inst_effect_h5X(.x, .y, hfile_sc[['X']], X_ntc,
                                            obs$gene_transcript, var$gene_id)) %>% list_rbind()

guide_effects <- guide_effects %>% mutate(
  Z = inst_cor/cor_se, p = pt(Z, df = n-2), p_adj = p.adjust(p, method='fdr')) %>%
  arrange(target, p_adj) %>% filter(!duplicated(target)) %>%
  left_join(select(var, gene_id, cv), by = c("target" = "gene_id"))
save(guide_effects, file=save_guide_res)

# Filter on beta < -0.75 and at least 50 targeting guides.
keep_guides <- filter(guide_effects, inst_beta < -0.75, n > n_ntc+50)
# keep_guides <- filter(guide_effects, inst_beta < -2, n > n_ntc+50)
cat("Keeping ", nrow(keep_guides), "genes/guides.\n")
# Reorder to match the order the variables appear in the hfile,
# prevents us from having to do this every time later.
keep_guides <- keep_guides[match(var$gene_id[var$gene_id %in% keep_guides$target],
                                 keep_guides$target), ]
targets <- as.list(keep_guides$target)
names(targets) <- keep_guides$inst_id
gc_res <- gc(reset=TRUE)
cat("Prepocessing time: ", difftime(Sys.time(), process_start, units="mins"),
    "(m). Memory used: ", gc_res[11] + gc_res[12], "MB\n")

# Fit the model with DAG=T and DAG=F.
dag_start <- Sys.time()
cat("Fitting inspre with DAG=TRUE.\n")
res_dag <- fit_inspre_from_h5X(hfile_sc[['X']], X_ntc, obs$gene_transcript, var$gene_id, targets,
                               max_med_ratio = 50, cv_folds = 5, ncores = ncores,
                               DAG = TRUE, min_nz = 0.03/2, lambda_min_ratio = 0.1, nlambda = 10)
save(res_dag, file=save_dag_res)
gc_res <- gc(reset=TRUE)
cat("DAG time: ", difftime(Sys.time(), dag_start, units="mins"),
    "(m). Memory used: ", gc_res[11] + gc_res[12], "MB\n")

# nodag_start <- Sys.time()
# cat("Fitting inspre with DAG=FALSE.\n")
# res_nodag <- fit_inspre_from_h5X(hfile_sc[['X']], X_ntc, obs$gene_transcript, var$gene_id, targets,
#                                   max_med_ratio = 50, lambda_min_ratio = 0.1, nlambda = 5,
#                                  cv_fold = 5, ncores = ncores, DAG = FALSE)
# save(res_nodag, file=save_nodag_res)
# gc_res <- gc(reset=TRUE)
# cat("No DAG time: ", difftime(Sys.time(), nodag_start, units="mins"),
#     "(m). Memory used: ", gc_res[11] + gc_res[12], "MB\n")
