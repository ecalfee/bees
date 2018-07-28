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

# how much buffer and DNA should I add to get ~40 uL at 0.35 ng/uL concentration?
# without having to pippet less than 3uL at any point
v = lapply(x_all, function(x) t(sapply(x, function(i) 
  bufferAndDNA2add(x = i, vol_tot = 40, target = 0.35))))
ids = lapply(1:4, function(i) read.table(paste0("../labwork/ids2extract_", i)))
extract = lapply(1:4, function(i) rep(i, nrow(ids[[i]])))
extractN = lapply(1:4, function(i) 1:nrow(ids[[i]])) 
d = data.frame(extract = unlist(extract), N = unlist(extractN),
               ids = do.call(rbind, ids)$x, ng_uL = unlist(x_all),
               do.call(rbind, v))
rownames(d) <- NULL
d$tot = round(d$DNA + d$buffer)
write.table(d, "../labwork/dilution_numbers_extract1-4.txt", 
            row.names = F, col.names = T, sep = "\t", quote = F)
