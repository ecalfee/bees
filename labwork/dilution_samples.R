# dilution concentrations
# goal = concentrations of 0.25 to 0.5 ng/uL
# original qubit concentration reads:
x1 <- c(2.94, 4.84, 2.90, 4.25)
x2 <- c(3.01, 1.19, 1.98, 3.03, 1.94, 6.10, 3.48, 4.84, 8.25, 4.69, 5.65, 3.46)
x3 <- c(5.50, 7.20, 5.00, 15.20, 4.93, 7.90, 4.42, 3.84, 11.20, 5.00, 5.95, 6.25, 
        10.10, 4.95, 6.75, 4.35, 4.57, 4.52, 2.90, 5.70, 5.25, 7.80, 12.10, 2.61)
x4 <- c(3.36, 2.73, 5.15, 14.10, 3.14, 5.35, 5.75, 4.15, 2.40, 3.13, 2.82, 5.50, 
        3.69, 2.59, 2.81, 3.95, 5.85, 3.67, 4.29, 4.32, 4.78, 3.61, 1.64)
x_all <- list(x1=x1, x2=x2, x3=x3, x4=x4)
# new concentration = (x * volume_x) / (volume_x + volume_buffer)
buffer2add = function(x, vol_x, target) {
  vol_b = vol_x * (x/target - 1)
    return(vol_b)
}
bufferAndDNA2add = function(x, vol_tot, target){
  # at least 3 microliters, otherwise round value to end up with around vol_tot total
  vol_x = ceiling(vol_tot * target / x) 
  vol_x = ifelse(vol_x<3, 3, vol_x)
  vol_b = buffer2add(x = x, vol_x = vol_x, target = target) # find volume of buffer to add
  return(c(DNA = vol_x, buffer=round(vol_b, digits = 3))) # round to 3 digits
}
# how much buffer should I add to 5 nanoliters of each library to get 0.35 ng/uL concentration?
#v5 = lapply(x_all, function(x) buffer2add(x = x, vol_x = 5, target = 0.35))
#v3 = lapply(x_all, function(x) buffer2add(x = x, vol_x = 3, target = 0.35))
#v10 = lapply(x_all, function(x) buffer2add(x = x, vol_x = 10, target = 0.35))

# how much buffer and DNA should I add to get ~40 uL at 0.5 ng/uL concentration?
# without having to pippet less than 3uL at any point
v = lapply(x_all, function(x) t(sapply(x, function(i) 
  bufferAndDNA2add(x = i, vol_tot = 40, target = 0.5))))
ids = lapply(1:4, function(i) read.table(paste0("../labwork/ids2extract_", i)))
extract = lapply(1:4, function(i) rep(i, nrow(ids[[i]])))
extractN = lapply(1:4, function(i) 1:nrow(ids[[i]])) 
d = data.frame(extract = unlist(extract), N = unlist(extractN),
               ids = do.call(rbind, ids)$x, ng_uL = unlist(x_all),
               do.call(rbind, v))
rownames(d) <- NULL
d$tot = round(d$DNA + d$buffer)
write.table(d, "../labwork/dilution0.5_numbers_extract1-4.txt", 
            row.names = F, col.names = T, sep = "\t", quote = F)

# first measurements were off, I did an initial dilution to get samples in the ~ 5-10ng/uL range
# and then re-measured DNA concentrations to do a final dilution down to 0.5 concentration
new <- read.table("../labwork/dilution_new_QUBIT_first_dilution_extract_1-4.txt",
                  header = T, stringsAsFactors = F)
new$DNA_2_add = 5
new$buffer_2_add = buffer2add(x = new$Qubit_HS_ng_uL, vol_x = new$DNA_2_add, target = 0.5)
new$"DNA:buffer_D1" = paste(new$DNA_added, new$buffer_added, sep = ":")
#new$"DNA:buffer_D0.5" = paste(new$DNA_2_add, new$buffer_2_add, sep = ":")
write.table(new[, c("extract", "N", "ids", "DNA:buffer_D1", "Qubit_HS_ng_uL", "DNA_2_add", "buffer_2_add")], "../labwork/dilution0.5_numbers_extract1-4_NEW.txt", 
            row.names = F, col.names = T, sep = "\t", quote = F)

# one at a time
lane5 <- c(5.76, 8.86, 5.68, 4.8, 4.52, 5.16, 5.8, 6.56, 4.48, 5.32, 7.32, 7.84, 5.84, 7.88, 5.80, 16.6, 6.48, 8.2, 6.72, 7.88, 6.4, 5.84, 7.76, 7.48, 5.3, 7.2, 5.72, 8.04, 9.28, 5.68, 5.48, 9.36, 5.52, 4.96, 4.48, 6.16, 6.0, 6.44, 7.08, 7.84, 6.28, 7.28, 6.12, 5.64, 6.8, 7.48, 2.68, 7.48, 6.64, 6.52, 6.36, 7.08, 6.16, 2.78, 3.75, 2.63)
data.frame(undil_conc = lane5, t(sapply(lane5, function(x) bufferAndDNA2add(x = x, vol_tot = 100, target = 0.375))))
redo5 <- data.frame(ID = c("12-9", "12-10", "12-11", "13-7", "13-10", "13-15", "13-17", "13-24"), 
                    dilA = c(0.584, 9.76, 0.696, 0.576, 0.592, 0.572, 0.724, 0.588))
data.frame(redo5, t(sapply(redo5$dilA, function(x) 
  bufferAndDNA2add(x = x, vol_tot = 100, target = 0.5)))) %>%
  dplyr::arrange(desc(dilA))


lane2 <- c(3.20, 3.31, 2.72, 2.89, 3.4, 2.68, 2.87, 3.06, 2.11, 2.72, 2.11, 3.51, 3.64, 2.18, 2.97, 2.14, 3.58, 2.69, 3.22, 4.40, 3.30, 2.63, 2.31, 2.94, 2.7, 2.34, 2.15, 3.12, 3.58, 2.67, 1.99, 2.82, 3.18, 2.44, 1.97, 3.11, 2.13, 2.18, 2.46, 1.71, 2.6, 2.27, 2.09, 4.08, 2.72, 2.51, 4.52, 3.39, 2.72, 2.1, 2.44, 3.19, 1.8, 2.0, 1.39, 1.59)
data.frame(undil_conc = lane2, t(sapply(lane2, function(x) bufferAndDNA2add(x = x, vol_tot = 100, target = 0.375))))
lane3 <- c(3.65, 1.64, 4.46, 3.9, 2.63, 3.98, 3.1, 3.22, 2.66, 5.36, 3.67, 3.7, 3.48, 3.52, 3.81, 1.54, 2.63, 3.95, 3.36, 4.28, 3.03, 2.24, 4.04, 4.2, 3.2, 2.68, 3.21, 2.81, 3.16, 4.72, 3.12, 3.38, 2.93, 3.03, 3.91, 3.46, 2.55, 5.48, 3.68, 2.98, 3.82, 3.77, 1.42, 3.46, 4.28, 3.18, 2.93, 3.17, 2.68, 3.68, 3.17, 2.59, 4.64, 4.32, 1.83, 2.51)
data.frame(undil_conc = lane3, t(sapply(lane3, function(x) bufferAndDNA2add(x = x, vol_tot = 100, target = 0.375))))
lane4 <- c(3.75, 2.16, 3.75, 3.5, 2.86, 2.46, 2.57, 3.4, 2.97, 3.01, 2.64, 1.96, 2.47, 3.69, 2.94, 5.04, 3.34, 2.96, 2.79, 2.82, 3.0, 2.96, 3.35, 3.06, 2.96, 3.16, 4.04, 3.06, 3.26, 3.47, 4.44, 2.7, 3.57, 2.54, 2.06, 3.03, 2.93, 2.83, 2.2, 1.32, 2.02, 2.58, 2.25, 2.04, 3.25, 3.22, 2.78, 2.98, 2.53, 2.45, 1.94, 2.05, 2.63, 3.58, 2.33, 1.55)
data.frame(undil_conc = lane4, t(sapply(lane4, function(x) bufferAndDNA2add(x = x, vol_tot = 100, target = 0.375))))


