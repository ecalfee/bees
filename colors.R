library(RColorBrewer)
library(scales)
# color palettes to use in R plots:
# IBM color-blind palette + green, grey and black
# don't mix the green and grey
# also not clear if faded versions of the same colors are distinguishable
# https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000-%2322542F-%23464343-%23111111
#col_ibm <- c("#23648FFF", "#23785EF0", "#23DC267F", "#23FE6100", 
#  "#23FFB000", "#2322542F", "#23464343", "#23111111")
col_ibm <- c("#648FFF", "#785EF0", "#DC267F",
             "#FE6100", "#FFB000", "#22542F",
             "#464343", "#111111")
col_blind <- cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#display.brewer.all(colorblindFriendly = TRUE)
dark2 <- brewer.pal(n = 8, "Dark2")
paired <- brewer.pal(n = 12, "Paired")
#col_ACM <- col_ibm[c(3,5,4)]
#col_NA_SA_both <- col_ibm[c(1, 2, 7)]
ACM = c("A", "C", "M")
col_ACM = dark2[c(4,3,1)]
names(col_ACM) = c("A", "C", "M")
col_ACM_all = dark2[c(2,4,3,1)] 
names(col_ACM_all) <- c("Combined", "A", "C", "M")
col_NA_SA_both <- dark2[c(8, 6, 5)]
names(col_NA_SA_both) <- c("Combined", "N. America", "S. America")
#viridisLite::viridis(3)
# alternatively I could use plasma colors for ancestry
brewer.pal(n = 5, "Paired") # paired is color-blind friendly
dark_brew <- brewer.pal(n = 8, "Dark2")
#col_ACM <- dark_brew[c(2,1,3)]
#col_ACM <- dark_brew[c(4,1,3)]
#col_NA_SA_both <- dark_brew[c(5,1)]

#plot(1:3, 1:3, cex = 5, pch = 20, col = viridis(3, begin = 0.3, end = 1))
col_FDR <- c("grey", "lightgrey", "#440154FF", "#FDE725FF", "#2FB47CFF")
names(col_FDR) <- c("n.s. - even chr", "n.s. - odd chr", "0.1", "0.05", "0.01")
#plot(1:5, 1:5, cex = 5, pch = 20, col = magma(5, begin = 0, end = 1))

#https://venngage.com/blog/color-blind-friendly-palette/
# retro color palette:
retro <- c("#ee442f", "#601a4a", "#63abce")
#plot(1:3, 1:3, cex = 5, pch = 20, col = retro)
#col_ACM <- retro

# larger set of colors for Fst contrasts:
col_fst <- c(hue_pal()(6), hue_pal()(6)[4:6])
#col_fst <- viridis(9)
#col_fst <- plasma(9)
fst_names <- paste("fst", c("M_A", "C_A", "C_M", "A.NA_A", "A.SA_A", "A.NA_SA", "M.NA_M", "M.SA_M", "M.NA_SA"), sep = ".")
names(col_fst) <- fst_names
labels_col_fst <- c("A vs. M", "A vs. C", "M vs. C",
                    "A vs. N. America (within A)", "A vs. S. America (within A)", "N. America vs. S. America (within A)",
                    "M vs. N. America (within M)", "M vs. S. America (within M)", "N. America vs. S. America (within M)")
names(labels_col_fst) <- fst_names

col_pi_predictions <- dark2[c(2,8)]
names(col_pi_predictions) <- c("observed", "predicted_ref")
col_pi_predictions3 <- dark2[c(2,7,8)] 
names(col_pi_predictions3) <- c("observed", "predicted_admix", "predicted_ref")


shapes_sig <- c(17,16,15,4)
names(shapes_sig) <- c("0.01", "0.05", "0.1", "n.s.")
shapes_ACM_NA_SA <- c(0,1,2,17,19)
names(shapes_ACM_NA_SA) <- c("A", "C", "M", "N. America", "S. America") #c(col_ACM, col_NA_SA_both)
col_geo_label <- c(dark2[c(6, 5)], dark2[4], dark2[2], dark2[3], paired[2], dark2[8], dark2[1], dark2[7])
names(col_geo_label) <- c("N. America", "S. America", "Kenya", "S. Africa", "Croatia", "Germany", "Slovenia", "Poland", "Spain") 

col_low_high_A <- c("purple", "red")
names(col_low_high_A) <- c("Low A", "High A")


col_qtl <- c("deeppink", "orange", "blue", "purple")
names(col_qtl) <- c("Defense response", "Grooming", "Hygienic behavior", "Varroa Sensitive Hygiene")

col_top2_outliers <- c("orange", "blue", "black")
shapes_top2_outliers <- c(15, 17, 1)
names(col_top2_outliers) <- c("Chr1:10921910", "Chr11:14490150", "Genomewide mean")
names(shapes_top2_outliers) <- c("Chr1:10921910", "Chr11:14490150", "Genomewide mean")
