#!/usr/bin/env Rscript

#### load libraries, directories and files ####

# load necessary libraries
library(ggVennDiagram)
library(ggplot2)



# specify path to BOLDigger output directory
boldigger_dir <- "miseq/COI/BOLDigger"

# specify path to generated meta directory
gen_meta <- "metadata/generated_meta"

# specify path to figures output
figures_dir <- "miseq/figures"
if(!dir.exists(figures_dir)) dir.create(figures_dir)



# read the generated meta file COI_demultiplexed_summary.csv from demultiplexed_summary.R
metadata_sum <- read.csv(file = file.path(gen_meta, "COI_demultiplexed_summary.csv"),
                         sep = ",", header = TRUE)

# read the combined read and taxonomy table from BOLDigger_tax_table.R
MOTU_tax_tab <- read.csv(file = file.path(boldigger_dir, "COI_motu_table_tax_counts_species_fullname.txt"),
                         sep = "\t", header = TRUE)




# remove all rows where genus and species are NA
MOTU_tax_tab <- MOTU_tax_tab[!(is.na(MOTU_tax_tab$genus) & is.na(MOTU_tax_tab$species)), ]
species_tax_tab <- MOTU_tax_tab[!is.na(MOTU_tax_tab$species), ]

# create subset of taxonomic  information to extract
taxonomy_info <- data.frame(MOTU_tax_tab[c("MOTU", "phylum", "family", "genus", "species")])

#### extract pyhylum information from the total MOTU taxonomy table ####
colors <- c("#CAFFBF", "#fdffb6", "#ffd6a5", "#ffadad", "#ffc6ff", "#bdb2ff", "#a0c4ff")
gradient_colors <- colorRampPalette(colors)(21)

phylum_sum_barplot <- ggplot(species_tax_tab, aes(x = phylum, fill = phylum)) +
  geom_bar(color = "#0F30b4") + 
  
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.4, color = "#0F30b4", size = 3.3) +
  
  scale_fill_manual(values = gradient_colors) +
  
  ylab("No. of identified species") +
  
  theme_minimal() + 
  
  theme(plot.background = element_rect(fill = "#e8e8e8"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
        
        legend.position = "none",
        axis.line.x = element_line(color = "#0F30b4", size = 0.5),
        axis.line.y = element_line(color = "#0F30b4", size = 0.5),
        
        axis.text = element_text(color = "#0F30b4"),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(color = "#0F30b4", size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.ticks = element_line(color = "#0F30b4"),
        axis.ticks.length = unit(0.2, "cm"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

phylum_sum_barplot

ggsave("phylum_species_2018_2020.png", plot = phylum_sum_barplot, 
       dpi = 400, path = figures_dir)  

#### refine taxonomic tables and assign to respective year ####
# extract the sample names from the samples retrieved at specific dates
# this way, the ARMS units retrieved earlier are excluded
samples_2018 <- metadata_sum$Gene_COI[metadata_sum$DateCollected == "2019-05-27"]
samples_2019 <- metadata_sum$Gene_COI[metadata_sum$DateCollected == "2020-07-16"]
samples_2020 <- metadata_sum$Gene_COI[metadata_sum$DateCollected == "2021-06-08"]

# create file of reads from each year
taxtab_2018 <- data.frame(taxonomy_info, MOTU_tax_tab[colnames(MOTU_tax_tab) %in% samples_2018])
taxtab_2019 <- data.frame(taxonomy_info, MOTU_tax_tab[colnames(MOTU_tax_tab) %in% samples_2019])
taxtab_2020 <- data.frame(taxonomy_info, MOTU_tax_tab[colnames(MOTU_tax_tab) %in% samples_2020])

# remove zero reads and singleton rows
taxtab_2018 <- taxtab_2018[rowSums(taxtab_2018[sapply(taxtab_2018, is.numeric)], na.rm = TRUE) > 1, ]
taxtab_2019 <- taxtab_2019[rowSums(taxtab_2019[sapply(taxtab_2019, is.numeric)], na.rm = TRUE) > 1, ]
taxtab_2020 <- taxtab_2020[rowSums(taxtab_2020[sapply(taxtab_2020, is.numeric)], na.rm = TRUE) > 1, ]


#### create bar plot of species found 2018, 2019 and 2020 ####

species_df <- data.frame(
  species = c(sum(taxtab_2018$species != "NA", na.rm = TRUE),
    sum(taxtab_2019$species != "NA", na.rm = TRUE),
    sum(taxtab_2020$species != "NA", na.rm = TRUE)),
  year = c("2018", "2019", "2020"))

colors <- c("#CAFFBF", "#fdffb6", "#ffd6a5", "#ffadad", "#ffc6ff", "#bdb2ff", "#a0c4ff")
gradient_colors <- colorRampPalette(colors)(21)

# species_sum_barplot <- 
ggplot(species_df, aes(x = year, y = species, fill = year)) +
  
  geom_col(color = "#0F30b4", 
           fill = c("#fdffb6", "#ffadad", "#bdb2ff"),
           width = 0.55) +
  
  geom_text(aes(label = species), vjust = -0.5, color = "#0F30b4", size = 6) + 
  
  ylab("No. of identified species") +
  
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        
        legend.position = "none",
        axis.line.x = element_line(color = "#0F30b4", size = 0.5),
        axis.line.y = element_line(color = "#0F30b4", size = 0.5),
        
        axis.text = element_text(color = "#0F30b4", size = 18),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(color = "#0F30b4", size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.ticks = element_line(color = "#0F30b4"),
        axis.ticks.length = unit(0.2, "cm"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        plot.margin = unit(c(2, 2, 2, 2), "cm"))

ggsave("species_sum_barplot_2018_2020.png", plot = species_sum_barplot, 
       dpi = 400, path = figures_dir)                           
  
#### create venn diagrams of genera and species for 2018, 2019 and 2020 ####

# create list of species from each year with NAs removed
family_list <- list(taxtab_2018$family[!is.na(MOTU_tax_tab$family)],
                     taxtab_2019$family[!is.na(MOTU_tax_tab$family)],
                     taxtab_2020$family[!is.na(MOTU_tax_tab$family)])

# create venn diagram from list
family_venn <- ggVennDiagram(list(taxtab_2018$family, taxtab_2019$family, taxtab_2020$family),
                              category.names = c("Family 2018", "Family 2019", "Family 2020"),
                              label_alpha = 0) + 
  scale_fill_gradient(low = "#CAFFBF", high = "#FFC6FF") +
  scale_x_continuous(expand = expansion(mult = .3)) +
  labs(caption = Sys.Date())

# save venn diagram to file in figures_dir
ggsave("venn_family_2018_to_2020.png", plot = family_venn, 
       dpi = 300, path = figures_dir)



# create list of genera from each year with NAs removed
genera_list <- list(taxtab_2018$genus[!is.na(MOTU_tax_tab$genus)],
                   taxtab_2019$genus[!is.na(MOTU_tax_tab$genus)],
                   taxtab_2020$genus[!is.na(MOTU_tax_tab$genus)])

# create venn diagram from list
genera_venn <- ggVennDiagram(genera_list,
                             category.names = c("Genera 2018", "Genera 2019", "Genera 2020"),
                             label_alpha = 0) + 
  scale_fill_gradient(low = "#FDFFB6", high = "#BDB2FF") +
  scale_x_continuous(expand = expansion(mult = .3)) +
  labs(caption = Sys.Date())

# save venn diagram to file in figures_dir
ggsave("venn_genera_2018_to_2020.png", plot = genera_venn, 
       dpi = 300, path = figures_dir)



# create list of species from each year with NAs removed
species_list <- list(taxtab_2018$species[!is.na(MOTU_tax_tab$species)],
                   taxtab_2019$species[!is.na(MOTU_tax_tab$species)],
                   taxtab_2020$species[!is.na(MOTU_tax_tab$species)])

# create venn diagram from list
species_venn <- ggVennDiagram(list(taxtab_2018$species, taxtab_2019$species, taxtab_2020$species),
              category.names = c("Species 2018", "Species 2019", "Species 2020"),
              label_alpha = 0) + 
  scale_fill_gradient(low = "#FFD6A5", high = "#A0C4FF") +
  scale_x_continuous(expand = expansion(mult = .3)) +
  labs(caption = Sys.Date())

# save venn diagram to file in figures_dir
ggsave("venn_species_2018_to_2020.png", plot = species_venn, 
       dpi = 300, path = figures_dir)

