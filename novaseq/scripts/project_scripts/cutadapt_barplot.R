#!/usr/bin/env Rscript

library(ggplot2)

cutadapt_subset <- "novaseq/COI/cutadapt_subset"

filenames <- list.files(cutadapt_subset, pattern = ".txt", full.names = TRUE)

fwd_adapters <- c()
rev_adapters <- c()
var_names <- c()

for (file in filenames){
    table <- read.table(file, header = TRUE, sep = "\t")
    name <- sub("_ABBPOSTA_cutadapt_summary.txt", "", file)
    name <- sub("novaseq/COI/cutadapt_subset/", "", name)
    assign(paste0(name), table)

    fwd_adapters <- append(fwd_adapters, table$w.adapters)
    rev_adapters <- append(rev_adapters, table$w.adapters2)
    var_names <- append(var_names, name)
}

nsamples = 2191581
nbars = 63
nsubsamples = 34787


adapters <- rbind(fwd_adapters, rev_adapters)
colnames(adapters) <- seq(nsubsamples, nsamples, nsubsamples)

tot_reads <- table$in_reads



fwd_data <- data.frame(
    name = colnames(adapters),
    value = adapters[1,]/tot_reads
    )


# fwd adapter distribution
fwd_plot <- ggplot(fwd_data, aes(x = factor(name, levels = name), y = value)) +

  geom_col(color = "#0F30b4", 
           fill = "#e6e6e6",
           width = 0.7,
           linewidth = 0.3) +

  ylab("Proportion of primer sequence present") +
  xlab("Sequence order") +
  ggtitle("Primer content in VH1_2022_SF400 forward sequence") + 

  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     expand = c(0,0)) +

  scale_x_discrete(breaks = function(x) x[seq(0, length(x), by = 9)]) +

  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        
        plot.title = element_text(color = "#0F30b4", 
                                  size = 12,
                                  margin = margin(t = 0, r = 0, b = 0, l = 25)),

        legend.position = "none",
        axis.line.x = element_line(color = "#0F30b4", size = 0.4, linetype = 1),
        axis.line.y = element_line(color = "#0F30b4", size = 0.4),
        axis.text.y = element_text(color = "#0F30b4", size = 10),
        axis.text.x = element_text(color = "#0F30b4", size = 10),
        axis.title = element_text(color = "#0F30b4", size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.ticks = element_line(color = "#0F30b4"),
        axis.ticks.length = unit(0.2, "cm"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        plot.margin = unit(c(2, 2, 2, 2), "cm"))


rev_data <- data.frame(
    name = colnames(adapters),
    value = adapters[2,]/tot_reads
    )


# fwd adapter distribution
rev_plot <- ggplot(rev_data, aes(x = factor(name, levels = name), y = value)) +

  geom_col(color = "#0F30b4", 
           fill = "#e6e6e6",
           width = 0.7,
           linewidth = 0.3) +

  ylab("Proportion of primer sequence present") +
  xlab("Sequence order") +
  ggtitle("Primer content in VH1_2022_SF400 reverse sequence") + 

  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     expand = c(0,0)) +

  scale_x_discrete(breaks = function(x) x[seq(0, length(x), by = 9)]) +

  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        
        plot.title = element_text(color = "#0F30b4", 
                                  size = 12,
                                  margin = margin(t = 0, r = 0, b = 0, l = 25)),

        legend.position = "none",
        axis.line.x = element_line(color = "#0F30b4", size = 0.4, linetype = 1),
        axis.line.y = element_line(color = "#0F30b4", size = 0.4),
        axis.text.y = element_text(color = "#0F30b4", size = 10),
        axis.text.x = element_text(color = "#0F30b4", size = 10),
        axis.title = element_text(color = "#0F30b4", size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.ticks = element_line(color = "#0F30b4"),
        axis.ticks.length = unit(0.2, "cm"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        plot.margin = unit(c(2, 2, 2, 2), "cm"))


ggsave("fwd_VH1_2022_SF400_primers.png", plot = fwd_plot, 
       dpi = 300, path = cutadapt_subset)

ggsave("rev_VH1_2022_SF400_primers.png", plot = rev_plot, 
       dpi = 300, path = cutadapt_subset)
