for (a in 1:2){
  hets_small_sample[[a]] %>%
    cbind(hybrid_freqs[[a]][ , c("group_snp", "round_pos")], .) %>%
    group_by(round_pos) %>%
    summarise(A = mean(SA.A, na.rm = T), M = mean(SA.M, na.rm = T),
              C = mean(SA.C, na.rm = T)) %>%
    tidyr::gather(., "ancestry", "het", ACM) %>%
    ggplot(., aes(x = round_pos, y = het, color = ancestry)) +
    geom_point(size = .5) +
    xlab("Position (bp)") +
    ylab("pi") +
    scale_color_manual(values = col_ACM) +
    geom_smooth(span = smoother, method = "loess", fill = "black") +
    theme_classic() +
    ggtitle(paste0(outlier_regions$outlier_names[a], " - S. America pi within ancestries"))
  ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[a], "_SA_meanOver", k, "bp.png"),
         height = 4, width = 8, device = "png")
  hets_small_sample[[a]] %>%
    cbind(hybrid_freqs[[a]][ , c("group_snp", "round_pos")], .) %>%
    group_by(round_pos) %>%
    summarise(A = mean(NA.A, na.rm = T), M = mean(NA.M, na.rm = T),
              C = mean(NA.C, na.rm = T)) %>%
    tidyr::gather(., "ancestry", "het", ACM) %>%
    ggplot(., aes(x = round_pos, y = het, color = ancestry)) +
    geom_point(size = .5) +
    xlab("Position (bp)") +
    ylab("pi") +
    scale_color_manual(values = col_ACM) +
    geom_smooth(span = smoother, method = "loess", fill = "black") +
    theme_classic() +
    ggtitle(paste0(outlier_regions$outlier_names[a], " - N. America pi within ancestries"))
  ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[a], "_NA_meanOver", k, "bp.png"),
         height = 4, width = 8, device = "png")
  hets_small_sample[[a]] %>%
    cbind(hybrid_freqs[[a]][ , c("group_snp", "round_pos")], .) %>%
    group_by(round_pos) %>%
    summarise(A = mean(A, na.rm = T), M = mean(M, na.rm = T),
              C = mean(C, na.rm = T)) %>%
    tidyr::gather(., "ancestry", "het", ACM) %>%
    ggplot(., aes(x = round_pos, y = het, color = ancestry)) +
    geom_point(size = .5) +
    xlab("Position (bp)") +
    ylab("pi") +
    scale_color_manual(values = col_ACM) +
    geom_smooth(span = smoother, method = "loess", fill = "black") +
    theme_classic() +
    ggtitle(paste0(outlier_regions$outlier_names[a], " - pi within ancestries ref. panel"))
  ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[a], "_ACM_meanOver", k, "bp.png"),
         height = 4, width = 8, device = "png")
}

hets_small_sample[[1]] %>%
  cbind(hybrid_freqs[[1]][ , c("group_snp", "round_pos")], .) %>%
  group_by(round_pos) %>%
  summarise(A = mean(A, na.rm = T), SA.A = mean(SA.A, na.rm = T),
            NA.A = mean(NA.A, na.rm = T)) %>%
  tidyr::gather(., "group", "het", c("A", "SA.A", "NA.A")) %>%
  mutate(group = ifelse(group=="SA.A", "S. America", ifelse(group=="NA.A", "N. America", group))) %>%
  ggplot(., aes(x = round_pos, y = het, color = group)) +
  geom_point(size = .5) +
  xlab("Position (bp)") +
  ylab("pi") +
  scale_color_manual(values = c(col_NA_SA_both, col_ACM), name = "Population") +
  geom_smooth(span = smoother, method = "loess", fill = "black") +
  theme_classic()
ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[1], "_Api_meanOver", k, "bp.png"),
       height = 4, width = 8, device = "png")
hets_small_sample[[2]] %>%
  cbind(hybrid_freqs[[2]][ , c("group_snp", "round_pos")], .) %>%
  group_by(round_pos) %>%
  summarise(M = mean(M, na.rm = T), SA.M = mean(SA.M, na.rm = T),
            NA.M = mean(NA.M, na.rm = T)) %>%
  tidyr::gather(., "group", "het", c("M", "SA.M", "NA.M")) %>%
  mutate(group = ifelse(group=="SA.M", "S. America", ifelse(group=="NA.M", "N. America", group))) %>%
  ggplot(., aes(x = round_pos, y = het, color = group)) +
  geom_point(size = .5) +
  xlab("Position (bp)") +
  ylab("pi") +
  scale_color_manual(values = c(col_NA_SA_both, col_ACM), name = "Population") +
  geom_smooth(span = smoother, method = "loess", fill = "black") +
  theme_classic()
ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[2], "_Mpi_meanOver", k, "bp.png"),
       height = 4, width = 8, device = "png")


