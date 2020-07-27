library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
source("../colors.R") # get color palette

# load data
load(file = "results/wings_bees_all.RData") # loads bees_all; equivalent to reading in S1 as below
#bees_all <- read.table("../../bee_manuscript/files_supp/S1_Table_Sample_information.txt",
#                       header = T, stringsAsFactors = F, sep = "\t")

wings_all = bees_all %>%
  mutate(ruttner = (publication == "Published in Ruttner 1988, Oberursel Collection, Germany"))


# summarise mean wing length by ancestry for reference bees
acm.summary <- wings_all %>%
  filter(ruttner) %>% # only reference bees
  group_by(population) %>% # A/C/M
  summarise(mean_cm = mean(wing_length_cm),
            n = n())

# wing lengths from hybrid zone
wings <- read.table("results/bee_wing_lengths_cm_1.28.20.txt",
                    header = T, stringsAsFactors = F)

# # combine data
wings.meta <- wings %>%
  left_join(., dplyr::select(meta.ind, -c(date, time)), by = "Bee_ID") %>%
  left_join(., d_admix_ACM_labelled[ , c("Bee_ID", ACM)], by = "Bee_ID")


# wing length predicted by A ancestry proportion genomewide
m_wing <- with(filter(wings_all, !ruttner), lm(wing_length_cm ~ A))
summary(m_wing)
glance(m_wing)
tidy(m_wing)
# test for different effect by continent
m_wing_SA <- with(filter(wings_all, !ruttner) %>%
                    mutate(SA = as.numeric(location == "Argentina")), 
                  lm(wing_length_cm ~ A*SA))
summary(m_wing_SA)
glance(m_wing_SA)
tidy(m_wing_SA)

# plot
# wing length predicted by mean genomewide ancestry
p_wing_by_A <- filter(wings_all, !ruttner & !is.na(wing_length_cm)) %>% # all wing lengths from the hybrid zone
  mutate(zone = ifelse(location == "Argentina", "S. America", "N. America")) %>%
  ggplot(., aes(x = A, y = wing_length_cm*10)) + # plot in units mm
  annotate("rect", xmin = 0, xmax = 1, # linear fit on range 0,1
                ymin = 7.5, ymax = 9.5, fill = "grey", alpha = .2) + # grey background
  geom_point(size = 1, alpha = .75, aes(color = zone)) + # hybrid zone bees
  theme_classic() +
  geom_line(data = data.frame(A = c(0, 1), # linear fit on range 0,1
                   wing_length_cm = c(coefficients(m_wing)["(Intercept)"],
                   coefficients(m_wing)["(Intercept)"] + 
                                coefficients(m_wing)["A"])),
               color = "black") +
  xlab("A ancestry proportion") +
  ylab("Wing length (mm)") +
  geom_jitter(data = filter(wings_all, ruttner) %>% # reference bees
                dplyr::select(wing_length_cm, population) %>%
               left_join(., data.frame(population = c("A", "C", "M"),
                                       A = c(1.03, -0.03, -.03), # A ancestry for ref pops plot slightly outside of range
                                       stringsAsFactors = F),
                         by = "population"), width = 0.015, size = 1, alpha = 0.75, 
              aes(color = population)) +
  scale_color_manual(values = c(col_NA_SA_both, col_ACM), name = "")
p_wing_by_A
ggsave("plots/wing_length_by_A_ancestry.png", p_wing_by_A, height = 4, width = 5.2)
ggsave("../../bee_manuscript/figures/wing_length_by_A_ancestry.png", p_wing_by_A, height = 4, width = 5.2)
ggsave("../../bee_manuscript/figures_supp/wing_length_by_A_ancestry.tiff", 
       p_wing_by_A, compression = "lzw", type = "cairo",
       height = 4, width = 5.2)


# add model predictions and residuals to wings dataframe
d.wings <- wings_all %>%
  filter(!ruttner, !is.na(wing_length_cm)) %>%
  mutate(model_predicted_cm = stats::predict(m_wing),
         model_residual_cm = stats::residuals(m_wing)) %>%
  arrange(abs(latitude)) %>%
  dplyr::select(Bee_ID, wing_length_cm, model_predicted_cm, model_residual_cm)

save(list = c("d.wings", "m_wing", "m_wing_SA"), file = "results/wing_fits.RData")