# Figure Eggnog mapper

library(ggplot2)
library(dplyr)
library(tidyr)
library(KEGGREST)

dir.create("../figures/Fig6/", showWarnings = FALSE, recursive = TRUE)


# Read in eggnog mapper data
eggnog = read.csv('../data/zenodo/eggnog_results/eggnog_results_all.emapper.annotations',
                  header=T, sep='\t')
eggnog$id = eggnog$query
sort(table(eggnog$COG_category))

# Read in AMR data
amr = read.csv('../data/zenodo/amrfinder-output.tsv', header=T, sep='\t')
amr$id = amr$Protein.id
amr_proteins = amr$id

# Defense/antidefense
defensefinder = read.csv('../data/zenodo/defensefinder_results/all_genes.tsv',
                         header=T,
                         sep='\t')
antidefense = defensefinder[which(defensefinder$activity=="Antidefense"),]
defense = defensefinder[which(defensefinder$activity=="Defense"),]
antidefense_proteins = antidefense$hit_id
defense_proteins = defense$hit_id

# Combine these with all functional annotations
amr_df <- data.frame(id = amr_proteins) %>%
  mutate(is_amr = TRUE)
defense_df <- data.frame(id = defense_proteins) %>%
  mutate(is_defense = TRUE)
antidefense_df <- data.frame(id = antidefense_proteins) %>%
  mutate(is_antidefense = TRUE)

# Left join both to your eggnog dataframe
eggnog_annotated <- eggnog %>%
  left_join(amr_df, by = "id") %>%
  left_join(defense_df, by = "id") %>%
  left_join(antidefense_df, by = "id") %>%
  mutate(
    is_amr = replace_na(is_amr, FALSE),
    is_defense = replace_na(is_defense, FALSE),
    is_antidefense = replace_na(is_antidefense, FALSE)
  )

# Read in CAI data
cai = read.csv('../data/zenodo/cai_results_prodigal_meta.tsv',
               header=T,
               sep='\t')
cai$id = gsub(".*\\|", "", cai$gene_id)
cai$plasmid = sapply(cai$id, 
                     function(x) paste(head(strsplit(x, "_")[[1]], 2), collapse = "_"))
cai$position = as.numeric(gsub(".*_", "", cai$gene_id))
plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T,
                         row.names = 1)
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]
cai$PTU = plasmid_table[cai$plasmid, "PTU"]
cai$HostGenus = plasmid_table[cai$plasmid, "Genus"]
cai$HostSpecies = plasmid_table[cai$plasmid, "Species"]

# Mob types (not yet used)
mob.types = read.csv('../data/PTU-mob-types.csv', header=T)


cai.eggnog = merge(cai, eggnog_annotated, by="id", all.x=TRUE)




# Plot AMR / defense / anti-defense
cai.eggnog$beta.lactamase = grepl("Beta-lactamase",cai.eggnog$Description)
ggplot(cai.eggnog[cai.eggnog$is_amr==TRUE,], aes(beta.lactamase, cai))+
  geom_boxplot()

# Plot COG one-letter categories by CAI
# Only top categories
top.categories = head(names(sort(table(cai.eggnog$COG_category), decreasing=TRUE)), n=10)
ggplot(cai.eggnog[which(cai.eggnog$COG_category %in% top.categories),], aes(COG_category, cai))+
  geom_boxplot()

# Seems far too general
# Use KEGG KO categories instead
cai.eggnog.by.KEGG.KO = cai.eggnog %>% group_by(KEGG_ko) %>%
  summarise(mean_cai=mean(cai),
            n=length(cai),
            COG=unique(COG_category))
cai.eggnog.by.KEGG.KO = cai.eggnog.by.KEGG.KO[order(cai.eggnog.by.KEGG.KO$mean_cai, decreasing = TRUE),]
# Keep only those with >=10 examples
cai.eggnog.by.KEGG.KO.10 = cai.eggnog.by.KEGG.KO %>% filter(n>9)

# Histogram of these categories
p.histogram = ggplot(cai.eggnog.by.KEGG.KO.10, aes(mean_cai))+
  geom_histogram(binwidth = 0.01, fill="lightgrey", colour="black")+
  theme_bw()+
  theme(axis.text=element_text(colour="black"))+
  theme(panel.grid = element_blank())+
  xlab("Mean codon adaptation index")+
  ylab("Count")
ggsave("../figures/Fig8/fig8-histogram.pdf", width=6, height=4)

plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T,
                         row.names = 1)
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]




