library(dplyr)
#d <- read.table('filter_by_group.txt', sep=" ", stringsAsFactors = F) %>%
#pull out from the ALL.fam file for CA (2014) and Harpur SNPs 
#all the bees from CA (2014 Santiago) and C, M, A (Harpur) ancestries
d <- read.table('ALL.fam', sep=" ", stringsAsFactors = F) %>%
  select(V1, V2) %>%
  mutate(V3 = substr(V1, 1, 3)) %>%
  filter(V1!="Y") %>%
  filter(V1!="Cerana")
write.table(d, 'CA2014_CMA_bees.txt', sep = ' ', col.names = F, row.names = F, quote = F)

#creating a new position file with 1 morgan = 100,000,000 assumption roughly true by
#multiplying map distance in morgans by 100,000,000
#p <- read.table('CALI_Harpur_SNPs.bim', sep='\t', stringsAsFactors = F) 
p <- read.table('CA2014_CMA.bim', sep='\t', stringsAsFactors = F) # for full dataset
p1 <- p %>%
  filter(V1==1) %>%
  mutate(V7=100000000*V3) %>%
  select(V2, V7)
p1.text <- p %>%
  filter(V1==1)

write.table(p1, 'morgan_pos_chr_1.txt', sep = ' ', col.names = F, row.names = F, quote = F)

p12 <- p %>%
  filter(V1==12) %>%
  mutate(V7=100000000*V3) %>%
  select(V2, V7)

write.table(p12, 'morgan_pos_chr_12.txt', sep = ' ', col.names = F, row.names = F, quote = F)

