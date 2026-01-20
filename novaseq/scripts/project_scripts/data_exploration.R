#!/usr/bin/env Rscript

#### load libraries, directories and files ####

# load necessary libraries
library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(tidyr)


# specify path to BOLDigger output directories
miseq_boldigger_dir <- "miseq/COI/BOLDigger"
novaseq_boldigger_dir <- "novaseq/COI/BOLDigger"

# specify path to generated meta directory
gen_meta <- "metadata/generated_meta"

# specify path to figures output
figures_dir <- "novaseq/figures"
if(!dir.exists(figures_dir)) dir.create(figures_dir)


# read the generated meta file COI_demultiplexed_summary.csv from demultiplexed_summary.R
metadata_sum <- read.csv(file = file.path(gen_meta, "COI_demultiplexed_summary.csv"),
                         sep = ",", header = TRUE)

# read the combined read and taxonomy table from BOLDigger_tax_table.R
miseq_taxtab <- read.csv(file = file.path(miseq_boldigger_dir, "COI_motu_table_tax_counts_species_fullname.txt"),
                         sep = "\t", header = TRUE)

novaseq_taxtab <- read.csv(file = file.path(novaseq_boldigger_dir, "COI_motu_table_tax_counts_species_fullname.txt"),
                         sep = "\t", header = TRUE)




# remove all rows where genus and species are NA
# miseq_taxtab <- miseq_taxtab[!(is.na(miseq_taxtab$genus) & is.na(miseq_taxtab$species)), ]
miseq_taxtab <- miseq_taxtab[!is.na(miseq_taxtab$species), ]

# remove all rows where genus and species are NA
# novaseq_taxtab <- novaseq_taxtab[!(is.na(novaseq_taxtab$genus) & is.na(novaseq_taxtab$species)), ]
novaseq_taxtab <- novaseq_taxtab[!is.na(novaseq_taxtab$species), ]

# create subset of taxonomic  information to extract
miseq_taxinfo <- data.frame(miseq_taxtab[c("MOTU", "phylum", "family", "genus", "species")])
novaseq_taxinfo <- data.frame(novaseq_taxtab[c("MOTU", "phylum", "family", "genus", "species")])

# merge novaseq and miseq phylum category from species subset
comb_species_taxtab <- bind_rows(novaseq_taxtab %>% mutate(dataset = "2022"),
                                 miseq_taxtab %>% mutate(dataset = "2018-2020"))

# fill missing phyla from miseq
comb_species_taxtab <- comb_species_taxtab %>% mutate(phylum = factor(phylum, levels = sort(unique(phylum))))
count_species_taxtab <- comb_species_taxtab %>% count(phylum, dataset) %>% complete(phylum, dataset, fill = list(n = 0))

#### extract pyhylum information from the total MOTU taxonomy table ####
colors <- c("#CAFFBF", "#fdffb6", "#ffd6a5", "#ffadad", "#ffc6ff", "#bdb2ff", "#a0c4ff")
gradient_colors <- colorRampPalette(colors)(21)

novaseq_colors <- c("#CAFFBF", "#fdffb6", "#ffd6a5", "#ffadad", "#ffc6ff", "#bdb2ff", "#a0c4ff")
novaseq_gradient <- colorRampPalette(colors)(26)

comb_phylum_barplot <- ggplot(count_species_taxtab, aes(x = phylum, y = n, fill = dataset)) +
  geom_col(color = "#0F30b4", width = 0.9, position = position_dodge2(width = 0.9, preserve = "single")) +
  
  geom_text(aes(label = ifelse(n == 0, "", n)), vjust = -0.7, hjust = 0.5, color = "#0F30b4", size = 5,  position = position_dodge2(width = 0.9, preserve = "single")) +

  ylab("No. of identified species") +
  
  theme_minimal() + 

  scale_fill_manual(values = c("2022" = "#ffd6a5", "2018-2020" = "#a0c4ff")) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        
        axis.line.x = element_line(color = "#0F30b4", size = 0.6),
        axis.line.y = element_line(color = "#0F30b4", size = 0.6),
        
        axis.text = element_text(color = "#0F30b4"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "#0F30b4", size = 16), 
        axis.title = element_text(color = "#0F30b4", size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 19, b = 0, l = 19), size = 20),
        axis.ticks = element_line(color = "#0F30b4"),
        axis.ticks.length = unit(0.2, "cm"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        plot.background = element_rect(color = "#e8e8e8", fill = "#e8e8e8"),
        
        legend.position = c(0.92, 0.92),
        legend.background = element_rect(colour = "#0F30b4", fill = "#e8e8e8", linetype='solid'),
        legend.title = element_blank(),
        legend.text = element_text(color = "#0F30b4", size = 18),
        legend.key.size = unit(0.8, "cm"))

ggsave("combined_phylum_species.png", plot = comb_phylum_barplot, 
       dpi = 400, path = figures_dir, width = 18, height = 12) 

novaseq_phylum_barplot <- ggplot(novaseq_taxtab, aes(x = phylum, fill = phylum)) +
  geom_bar(color = "#0F30b4") + 
  
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, color = "#0F30b4", size = 5) +
  
  scale_fill_manual(values = novaseq_gradient) +
  
  ylab("No. of identified species") +
  
  theme_minimal() + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        
        axis.line.x = element_line(color = "#0F30b4", size = 0.6),
        axis.line.y = element_line(color = "#0F30b4", size = 0.6),
        
        axis.text = element_text(color = "#0F30b4"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "#0F30b4", size = 16), 
        axis.title = element_text(color = "#0F30b4", size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 19, b = 0, l = 19), size = 20),
        axis.ticks = element_line(color = "#0F30b4"),
        axis.ticks.length = unit(0.2, "cm"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        plot.background = element_rect(fill = "#e8e8e8", color = "#e8e8e8"),
        legend.position = "none")

ggsave("phylum_species_2022.png", plot = novaseq_phylum_barplot, 
       dpi = 400, path = figures_dir, width = 12, height = 12)  


#### refine taxonomic tables and assign to respective year ####
# extract the sample names from the samples retrieved at specific dates
# this way, the ARMS units retrieved earlier are excluded
samples_2018 <- metadata_sum$Gene_COI[metadata_sum$DateCollected == "2019-05-27"]
samples_2019 <- metadata_sum$Gene_COI[metadata_sum$DateCollected == "2020-07-16"]
samples_2020 <- metadata_sum$Gene_COI[metadata_sum$DateCollected == "2021-06-08"]

# create file of reads from each year
taxtab_2018 <- data.frame(miseq_taxinfo, miseq_taxtab[colnames(miseq_taxtab) %in% samples_2018])
taxtab_2019 <- data.frame(miseq_taxinfo, miseq_taxtab[colnames(miseq_taxtab) %in% samples_2019])
taxtab_2020 <- data.frame(miseq_taxinfo, miseq_taxtab[colnames(miseq_taxtab) %in% samples_2020])

# remove zero reads and singleton rows
taxtab_2018 <- taxtab_2018[rowSums(taxtab_2018[sapply(taxtab_2018, is.numeric)], na.rm = TRUE) > 1, ]
taxtab_2019 <- taxtab_2019[rowSums(taxtab_2019[sapply(taxtab_2019, is.numeric)], na.rm = TRUE) > 1, ]
taxtab_2020 <- taxtab_2020[rowSums(taxtab_2020[sapply(taxtab_2020, is.numeric)], na.rm = TRUE) > 1, ]
taxtab_2022 <- novaseq_taxinfo[rowSums(novaseq_taxtab[sapply(novaseq_taxtab, is.numeric)], na.rm = TRUE) > 1, ]

#### create bar plot of species found 2018, 2019 and 2020 ####

species_df <- data.frame(
  species = c(sum(taxtab_2018$species != "NA", na.rm = TRUE),
    sum(taxtab_2019$species != "NA", na.rm = TRUE),
    sum(taxtab_2020$species != "NA", na.rm = TRUE),
    sum(taxtab_2022$species != "NA", na.rm = TRUE)),
  year = c("2018", "2019", "2020", "2022"))


species_sum_barplot <- ggplot(species_df, aes(x = year, y = species, fill = year)) +
  
  geom_col(color = "#0F30b4",
           width = 0.55) +
  
  geom_text(aes(label = species), hjust = -0.5, color = "#0F30b4", size = 5) + 
  
  scale_fill_manual(values = c("#CAFFBF", "#ffd6a5", "#ffc6ff", "#a0c4ff")) +

  ylab("No. of identified species") +

  coord_flip() +
  
theme(axis.text.x = element_text(size = 18),
        
        axis.line.x = element_line(color = "#0F30b4", size = 0.6),
        axis.line.y = element_line(color = "#0F30b4", size = 0.6),
        
        axis.text = element_text(color = "#0F30b4"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "#0F30b4", size = 16, margin = margin(t = 0, r = 19, b = 0, l = 19)), 
        axis.title = element_text(color = "#0F30b4", size = 16),
        axis.title.x = element_text(margin = margin(t = 19, r = 0, b = 19, l = 0), size = 20),
        axis.ticks = element_line(color = "#0F30b4"),
        axis.ticks.length = unit(0.2, "cm"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        plot.background = element_rect(fill = "#e8e8e8", color = "#e8e8e8"),
        legend.position = "none")

ggsave("species_sum_barplot_2018_2022.png", plot = species_sum_barplot, 
       dpi = 400, path = figures_dir, width = 18, height = 8)                           
  
#### create venn diagrams of genera and species for 2018, 2019 and 2020 ####


# create list of species from each year with NAs removed
species_list <- list(miseq_taxtab$species[!is.na(miseq_taxtab$species)],
                   novaseq_taxtab$species[!is.na(novaseq_taxtab$species)])

# create venn diagram from list
species_venn <- ggVennDiagram(list(miseq_taxtab$species, novaseq_taxtab$species),
              category.names = c("Species 2018-2020", "Species 2022"),
              label_alpha = 0, set_color = "#0F30b4", label_color = "#0F30b4",
              label_txtWidth = 0.5, edge_size = 0.5, set_size = 8, label_size = 8) + 
  scale_fill_gradient(low = "#CAFFBF", high = "#ffadad") +
  scale_x_continuous(expand = expansion(mult = 0.4)) +
  theme(legend.position = c(0.88, 0.5),
  plot.background = element_rect(fill = "#e8e8e8", color = "#e8e8e8"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text = element_text(color = "#0F30b4", size = 16),
  legend.title = element_blank(),
  legend.key.size = unit(1.5, "cm"))

# save venn diagram to file in figures_dir
ggsave("venn_species_miseq_vs_novaseq.png", plot = species_venn, 
       dpi = 400, path = figures_dir, width = 12, height = 12)


