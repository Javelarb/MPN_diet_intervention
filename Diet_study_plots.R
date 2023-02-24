library(tidyverse)
library(ggpubr)
library(vegan)

options(scipen=10000, digits = 12)
set.seed(seed = 999)

metadata = read.delim("/media/julio/Storage/MPN/MPN_metadata.txt", row.names=1, check.names = F)
metadata$Subject = as.factor(metadata$Subject)
OTU_table = read.delim("/media/julio/Storage/MPN/metaphlan3/OTU_table_relab.txt", header = T, check.names = F, comment.char = "#", row.names = 1)
read_counts_NH = read.table("/media/julio/Storage/MPN/NH_read_counts.txt", check.names = F)

OTU_table_counts = as.data.frame(t(OTU_table[,-1])) %>% rownames_to_column() %>% mutate(., rowname = str_replace(rowname, "_relab", "")) %>% mutate(., rowname = str_replace(rowname, "_", "\\.")) %>% merge(., read_counts_NH, by.x = "rowname", by.y = "V1") %>% column_to_rownames("rowname")
OTU_table_counts = ((OTU_table_counts[,-ncol(OTU_table_counts)])*OTU_table_counts$V2)/100

species_OTU_table = OTU_table_counts[,grep("s__", colnames(OTU_table_counts))] 

#Check to see if there are any contaminants. Only 10 taxa should be there. See:
# https://files.zymoresearch.com/protocols/_d6305_d6306_zymobiomics_microbial_community_dna_standard.pdf
pos_ctrl =species_OTU_table[grep("Comm.std", rownames(species_OTU_table)),] %>% .[,!(colSums(.) == 0)] #Comm.std-001-3 seems like it got cross contaminated. Lots of phages in zymo's standard, maybe we should tell them.
neg_ctrl =species_OTU_table[grep("Neg.ctrl", rownames(species_OTU_table)),] %>% .[,!(colSums(.) == 0)]

#Filter out contaminants.
contaminants = c("k__Viruses|p__Viruses_unclassified|c__Viruses_unclassified|o__Herpesvirales|f__Alloherpesviridae|g__Cyprinivirus|s__Cyprinid_herpesvirus_3", colnames(neg_ctrl))
species_OTU_table = species_OTU_table[, !(colnames(species_OTU_table) %in% contaminants)]

species_OTU_table_noctrls = species_OTU_table[-grep("Neg.ctrl", rownames(species_OTU_table)),] %>% 
  .[-grep("Comm.std", rownames(.)),] %>% merge(metadata, ., by = "row.names") %>% column_to_rownames("Row.names")
species_OTU_table_noctrls = species_OTU_table_noctrls[!(species_OTU_table_noctrls$Subject == "13"),]
species_OTU_table_noctrls = species_OTU_table_noctrls[!(row.names(species_OTU_table_noctrls) %in% c("FEA03-9", "FEA24-6.2", "FEA24-9.2", "FEA25-1.2", "FEA26-1.2", "FEA26-15.2", "FEA26-6.2")),]

species_OTU_table_noctrls[,(ncol(metadata)+1):ncol(species_OTU_table_noctrls)] = (species_OTU_table_noctrls[,(ncol(metadata)+1):ncol(species_OTU_table_noctrls)]/rowSums(species_OTU_table_noctrls[,(ncol(metadata)+1):ncol(species_OTU_table_noctrls)]))*100

# Alpha diversity  
## Shannon
alpha_div = as.data.frame(diversity(species_OTU_table_noctrls[,(ncol(metadata)+1):ncol(species_OTU_table_noctrls)], index = "shannon"))

alpha_div = merge(alpha_div, metadata, by = "row.names")
names(alpha_div)[[2]] = "Shannon"

alpha_stats = alpha_div %>% group_by(Diet, week) %>% summarise(Avg_shan = mean(Shannon), sd_shan = sd(Shannon), n = n(), stderr_shan = sd_shan/sqrt(n))

plot1 = ggplot(data = alpha_div) +
  aes(x = week, y = Shannon) +
  annotate("rect", xmin=3, xmax= 13, ymin = -Inf, ymax = Inf, fill = "#FFEFE0") +
  geom_path(data = alpha_stats, inherit.aes = F, aes(x = week, y = Avg_shan, color = Diet), linewidth = 1) + 
  geom_errorbar(data = alpha_stats, inherit.aes = F, aes(x = week, y = Avg_shan, color = Diet, ymin = Avg_shan - stderr_shan, ymax = Avg_shan + stderr_shan), width = .5, size = 1) +
  geom_point(data = alpha_stats, inherit.aes = F, aes(x = week, y = Avg_shan, color = Diet), size = 2, shape = 15) +
  theme_classic() +
  ylim(2,3.3) +
  labs(y = "Shannon Index", x = "Week") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  scale_color_manual(values = c("deeppink1", "black"))

# Beta Diversity
MDS_out = metaMDS(species_OTU_table_noctrls[,(ncol(metadata)+1):ncol(species_OTU_table_noctrls)], trymax = 999)
MDS_points = as.data.frame(MDS_out$points) %>% merge(., metadata, by = "row.names") %>% arrange(., Subject, week)
MDS_points = MDS_points %>% arrange(Subject, week)

beta_stats = MDS_points %>% group_by(Diet, week) %>% summarise(Avg = mean(MDS1), sd = sd(MDS1), n = n(), stderr = sd/sqrt(n))

plot2 = ggplot(data = MDS_points) +
  aes(x = week, y = MDS1) +
  annotate("rect", xmin=3, xmax= 13, ymin = -Inf, ymax = Inf, fill = "#FFEFE0") +
  geom_path(data = beta_stats, inherit.aes = F, aes(x = week, y = Avg, group = Diet, color = Diet), linewidth = 1) + 
  geom_errorbar(data = beta_stats, inherit.aes = F, aes(x = week, y = Avg, color = Diet, ymin = Avg - stderr, ymax = Avg + stderr), width = .5, linewidth = 1) +
  geom_point(data = beta_stats, inherit.aes = F, aes(x = week, y = Avg, color = Diet), shape = 15, size = 2) +
  theme_classic() +
  labs(y = "MDS1", x = "Week") +
  ylim(-1,1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  scale_color_manual(values = c("deeppink1", "black"))

ggsave(filename = "/media/julio/Storage/MPN/Diet_study_Figure_5.png", plot = ggarrange(plot1, plot2, nrow = 1, labels = c("A", "B"), common.legend = T, legend = "right"), 
       device = "png", width = 9, height = 4, dpi = 600, bg = "white")
