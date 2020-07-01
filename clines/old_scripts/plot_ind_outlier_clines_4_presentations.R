########## ---------------------- MAKE SOME PLOTS OF CLINES ---------------- ##############

# load plot_clines.R for functions/data
# plot a set of clines from simulations:
mean_cline_mvn <- group_by(mvn$params, term) %>%
  summarise(mean = mean(estimate))
plot(meta.AR.order.by.lat$lat, 
     rep(0, length(meta.AR.order.by.lat$lat)),
     col = NULL,
     ylim = c(0, 1),
     xlim = range(d_A$lat[d_A$group == "AR_2018"]),
     ylab = "A ancestry",
     xlab = "latitude")
curve(logistic3(mu = mean_cline_mvn$mean[mean_cline_mvn$term == "mu"], 
                b = mean_cline_mvn$mean[mean_cline_mvn$term == "b"], 
                x), 
      from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
      to = range(d_A$lat[d_A$group == "AR_2018"])[2],
      n = 100,
      add = T,
      lwd = 2, 
      lty = 2)


# plot 1 steep cline from my data vs. 1 steep cline in the simulated data and median snp cline
outlier_clines <- c(which.max(clines$params$estimate[clines$params$term == "b"]),
                    which.min(clines$params$estimate[clines$params$term == "b"]),
                    which.max(clines$params$estimate[clines$params$term == "mu"]),
                    which.min(clines$params$estimate[clines$params$term == "mu"]))

outlier_mvn_clines <- c(which.max(mvn$params$estimate[mvn$params$term == "b"]),
                        which.min(mvn$params$estimate[mvn$params$term == "b"]),
                        which.max(mvn$params$estimate[mvn$params$term == "mu"]),
                        which.min(mvn$params$estimate[mvn$params$term == "mu"]))

plot_random_clines <- function(d, n, seed = 500, color = "darkgrey"){
  # set frame
  plot(abs(meta.AR.order.by.lat$lat), 
       rep(0, length(meta.AR.order.by.lat$lat)),
       col = NULL,
       ylim = c(0, 1),
       xlim = range(d_A$abs_lat),#[d_A$group == "AR_2018"]),
       ylab = "African ancestry frequency",
       xlab = "Degrees latitude from equator")
  # set seed
  set.seed(seed)
  # sample data
  sample_data <- sample(1:(nrow(d$params)/2), n, replace = F)
  # plot
  sapply(sample_data, function(i)
    curve(logistic3(mu = d$params$estimate[d$params$term == "mu" & d$params$snp_index == i], 
                    b = d$params$estimate[d$params$term == "b" & d$params$snp_index == i], 
                    -x), 
          from = range(d_A$abs_lat)[1],
          to = range(d_A$abs_lat)[2],
          #from = range(d_A$abs_lat[d_A$group == "AR_2018"])[1], 
          #to = range(d_A$abs_lat[d_A$group == "AR_2018"])[2],
          n = 100,
          col = "darkgrey",
          add = T,
          lwd = 1, 
          lty = 1))
}
# just the background random clines
png("plots/random_clines_100.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
dev.off()

png("plots/random_clines_100_plus_steep.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1, function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],
        
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1))
dev.off()

png("plots/random_clines_100_plus_introgress.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(2:3, function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1))
dev.off()

# all 3 outliers
png("plots/random_clines_100_plus_3outliers.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = -clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = -clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
        to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
}
)
dev.off()

# ADD Points to all of these plots!
png("plots/random_clines_100_plus_steep_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
})
dev.off()

png("plots/random_clines_100_plus_introgress_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(2:3, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
})
dev.off()

# all 3 outliers
png("plots/random_clines_100_plus_3outliers_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
}
)
dev.off()


# what does the N. American cline look like for the 'steep' SNP? Also looks steep.
# should exclude Avalon -- such a consistent outlier point can really affect cline shape
png("plots/random_clines_100_plus_3outliers_points_and_NA_at_steep_cline_with_and_without_Avalon.png", units = "in", res = 300, height = 6, width = 6)

plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = -clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = -clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
        to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
}
)


meta.CA.order.by.lat <- meta.pop %>%
  arrange(lat) %>%
  filter(zone == "N. America")

fit_pink_with_Avalon = nls_multstart(A ~ logistic3(x = lat, 
                                                   b = b, mu = mu),
                                     start_lower = list(b = -1, mu = min(d_A$abs_lat)),
                                     start_upper = list(b = 1, mu = max(d_A$abs_lat)),
                                     supp_errors = 'Y',
                                     iter = 250,
                                     convergence_count = 100,
                                     data = data.frame(A = unname(t(A[outlier_clines[1], 
                                                                      meta.CA.order.by.lat$population])),
                                                       lat = meta.CA.order.by.lat$lat, 
                                                       stringsAsFactors = F))
fit_pink = nls_multstart(A ~ logistic3(x = lat, 
                                       b = b, mu = mu),
                         start_lower = list(b = -1, mu = min(d_A$abs_lat)),
                         start_upper = list(b = 1, mu = max(d_A$abs_lat)),
                         supp_errors = 'Y',
                         iter = 250,
                         convergence_count = 100,
                         data = data.frame(A = unname(t(A[outlier_clines[1], 
                                                          meta.CA.order.by.lat$population[meta.CA.order.by.lat$population != "Avalon_2014"]])), 
                                           lat = meta.CA.order.by.lat$lat[meta.CA.order.by.lat$population != "Avalon_2014"], 
                                           stringsAsFactors = F))

curve(logistic3(mu = unlist(tidy(fit_pink)[tidy(fit_pink)$term == "mu", "estimate"]), 
                b = unlist(tidy(fit_pink)[tidy(fit_pink)$term == "b", "estimate"]), 
                x), 
      from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
      to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
      n = 100,
      col = "deeppink",
      add = T,
      lwd = 3, 
      lty = 1)
curve(logistic3(mu = unlist(tidy(fit_pink_with_Avalon)[tidy(fit_pink_with_Avalon)$term == "mu", "estimate"]), 
                b = unlist(tidy(fit_pink_with_Avalon)[tidy(fit_pink_with_Avalon)$term == "b", "estimate"]), 
                x), 
      from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
      to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
      n = 100,
      col = "deeppink",
      add = T,
      lwd = 3, 
      lty = 2)
dev.off()

# I can alternatively plot absolute latitude:
curve(logistic3(mu = -mean_cline_mvn$mean[mean_cline_mvn$term == "mu"], 
                b = -mean_cline_mvn$mean[mean_cline_mvn$term == "b"], 
                x), 
      from = abs(range(d_A$lat[d_A$group == "AR_2018"])[2]), 
      to = abs(range(d_A$lat[d_A$group == "AR_2018"])[1]),
      n = 100,
      col = "blue",
      add = F,
      lwd = 2, 
      lty = 2)

# the top clines with data from N. America too:
meta.pop.A <- meta.pop %>%
  mutate(meanA = apply(A[ , meta.pop$population], 2, mean)) %>%
  mutate(abs_lat = abs(lat))
png("plots/high_A_SA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$lat[meta.pop.A$zone == "S. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "S. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "S. America"], 
       A[outlier_clines[2], meta.pop$zone == "S. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["S. America"])
dev.off()

png("plots/high_A_NA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$abs_lat[meta.pop.A$zone == "N. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "N. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[2], meta.pop$zone == "N. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["N. America"])
dev.off()

# low SA outlier:
png("plots/low_A_SA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$lat[meta.pop.A$zone == "S. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "S. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "S. America"], 
       A[outlier_clines[3], meta.pop$zone == "S. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["S. America"])
dev.off()

png("plots/low_A_NA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$abs_lat[meta.pop.A$zone == "N. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "N. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[3], meta.pop$zone == "N. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["N. America"])
dev.off()
