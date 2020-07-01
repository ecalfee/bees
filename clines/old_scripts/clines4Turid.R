# make snp cline file for Webster group to share:
# other data comes from plot_clines.R
load("../local_ancestry/results/sites_rpos.RData")
sites_clines <- sites_rpos %>%
  mutate(snp_index = 1:nrow(.)) %>%
  left_join(., clines$params, by = "snp_index") %>%
  dplyr::select(snp_index, scaffold, chr, chr_n, pos, major, minor, 
                cM_Mb, pos_cM, term, estimate, std.error, statistic, p.value, conf.low, conf.high)

sites_clines %>%
  pivot_wider(., id_cols = c("snp_index", "scaffold", "chr", "chr_n", "pos", "major", "minor", 
                             "cM_Mb", "pos_cM", "term"),
              names_from = "term", names_sep = ".",
              values_from = c("estimate", "std.error", "p.value", "conf.low", "conf.high")) %>%
  write.table("results/ind_snp_nls_multistart_clines/clines_4_Turid.txt",
              col.names = T, row.names = F,
              sep = "\t", quote = F)