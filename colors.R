library(RColorBrewer)
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
col_ACM <- col_ibm[c(3,5,4)]
col_NA_SA_both <- col_ibm[c(1, 2, 7)]
#display.brewer.all(colorblindFriendly = TRUE)
#viridisLite::viridis(3)
# alternatively I could use plasma colors for ancestry
brewer.pal(n = 5, "Paired") # paired is color-blind friendly
dark_brew <- brewer.pal(n = 8, "Dark2")
#col_ACM <- dark_brew[c(2,1,3)]
#col_ACM <- dark_brew[c(4,1,3)]
#col_NA_SA_both <- dark_brew[c(5,1)]