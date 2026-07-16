library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("~/Desktop/photobiont turnover")

#load photobiont community data
OTU_matrix <- read.csv("matrix.csv")
#load metadata
phys_land_data <- read_csv("physio_land_all_correct.csv")

#compute proportion of each genotype in each thallus 

#remove sample ID column for now
sample_ID <- OTU_matrix$Sample_ID
rownames(OTU_matrix) <- OTU_matrix$Sample_ID
OTU_matrix <- OTU_matrix[,2:23]
col_names <- colnames(OTU_matrix)

#### Calculating proportion ####
#for each x, I want it to divide it by the sum of the entire row

read_sums <- rowSums(OTU_matrix)
proportion_matrix <- OTU_matrix/read_sums

#reload to get sample ID column
OTU_matrix <- read.csv("matrix.csv")

#merge to metadata
flacap_treb <- merge(OTU_matrix, phys_land_data, by = "Sample_ID")

#OTU as key (column names), sequences as value (row values)
flacap_treb_long <- flacap_treb %>%
  gather(OTU, sequences, 2:23)

flacap_treb_long$Sample_ID <- as.character(flacap_treb_long$Sample_ID)

glimpse(flacap_treb_long)

#order OTUs by hand
flacap_treb_long$OTU <- fct_relevel(flacap_treb_long$OTU, 
                                    "A13", "A19", "A33", "A40", "A44", "A46",
                                    "C03",
                                    "S01","S02","S05","S06","S10",
                                    "I01","I02","I05","I08")

#picking colors for OTUs
OTU_colors <- hcl.colors(n = 22, palette = "Roma")
OTU_colors <- rev(OTU_colors)

flacap_treb_long <- flacap_treb_long%>%
  arrange(desc((ndvi10)))

flacap_treb_long <- flacap_treb_long%>%
  mutate(Sample_ID = fct_reorder(Sample_ID, ndvi10, .desc = TRUE))

flacap_treb_long_filtered <- flacap_treb_long%>%
  filter(Project == "Urban_LCCMR" | Project == "Urban_Lindsey")

pop_key = c("Source_2023" = "Source_2023", 
            "Suburban transplant_6m" = "Transplant_6m",
            "Urban transplant_6m" = "Transplant_6m",
            "Source_2022" = "Source_2022", 
            "Suburban transplant_10m" = "Transplant_10m",
            "Suburban transplant_5m" = "Transplant_5m",
            "Urban transplant_10m" = "Transplant_10m",
            "Urban transplant_5m" = "Transplant_5m")

flacap_treb_long_filtered$Population <- fct_relevel(flacap_treb_long_filtered$Population,
                                                    "Source_2022", "Source_2023", "Transplant_5m", 
                                                    "Transplant_6m", "Transplant_10m")
#2023 Transplants
flacap_treb_long_filtered%>%
#  filter(Project == "Urban_Lindsey")%>%
    filter(Project == "Urban_LCCMR")%>%
  ggplot(aes(x = Sample_ID, y = sequences, fill = OTU))+
  facet_grid(Project~Population, scales = "free", space = "free")+
  geom_col(position = "fill", show.legend = TRUE)+
  #  facet_grid(~Population, scales = "free_x")+
  scale_fill_manual(values = OTU_colors)+
  theme(axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(),
        strip.placement = "outside")+
  labs(x = "Thalli from high to low NDVI 10m",
       y = "proportion")

ggsave("2023_facet.svg",
       plot = last_plot(), dpi = 400, width = 17, height = 10, units = "cm")


#2022 transplants
flacap_treb_long_filtered%>%
  filter(Project == "Urban_Lindsey")%>%
  #  filter(Project == "Urban_LCCMR")%>%
  ggplot(aes(x = Sample_ID, y = sequences, fill = OTU))+
  facet_grid(Project~Population, scales = "free", space = "free")+
  geom_col(position = "fill", show.legend = FALSE)+
  #  facet_grid(~Population, scales = "free_x")+
  scale_fill_manual(values = OTU_colors)+
  theme(axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(),
        strip.placement = "outside")+
  labs(x = "Thalli from high to low NDVI 10m",
       y = "proportion")

ggsave("2022_facet.svg",
       plot = last_plot(), dpi = 400, width = 17, height = 10, units = "cm")


