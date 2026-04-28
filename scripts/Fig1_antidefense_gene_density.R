# Plot antidefense gene density
library(ggplot2)
library(dplyr)

# Plasmid table
plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T, 
                         row.names = 1)
# Subset to those we have reliable leading region for
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]

# Read in defensefinder results
defensefinder = read.csv('../data/TableSX_defensefinder_results_prodigal_meta.tsv',
                         header=T,
                         sep='\t')
defensefinder$plasmid = sub("_[0-9]+$", "", defensefinder$sys_beg)
defensefinder$position = as.numeric(gsub(".*_", "", defensefinder$sys_beg))
# Add PTU
defensefinder$plasmid = sub("_[0-9]+$", "", defensefinder$sys_beg)
defensefinder$PTU = plasmid_table[defensefinder$plasmid, "PTU"]



# Values for plotting
colour.values = read.csv('colour-values-antidefense.csv',
                         header=T)
colour.values = colour.values %>% arrange(Subtype)
antidefense = defensefinder[which(defensefinder$activity=="Antidefense"),]
# Change "Other" to "PsiAB" - only works for this dataset!
antidefense$type[antidefense$type=="Other"] = "psiAB"

# Exclude one Anti_Pycsar for purpose of plotting
antidefense.plot = antidefense[which(antidefense$type!="Anti_Pycsar"),]
  p.all = ggplot(antidefense.plot, aes(position, fill=type))+
    geom_histogram()+theme_bw()+
    xlab("Start of system (ORF position)")+
    ylab("Frequency")+
    theme(legend.position = c(0.8, 0.8),
          legend.background = element_rect(colour="black"))+
    scale_fill_manual(values=c("#fecc5c","#fd8d3c","#e31a1c"))

# Normalise facet by total number of plasmids within each PTU?
  # Add majority MOB type for PTUs
  mob.types = read.csv('../data/PTU-mob-types.csv',
                       header=T)
  antidefense.plot = antidefense.plot %>% left_join(mob.types, by="PTU")

  facetPlot = function(df){
    return(ggplot(df, aes(position, fill=type))+
      geom_histogram()+theme_bw()+
      xlab("Start of system (ORF position)")+
      facet_wrap(~PTU, scales="free_y")+
      theme(legend.position = "none")+
      theme(panel.grid=element_blank())+
      scale_fill_manual(values=c("#fecc5c","#fd8d3c","#e31a1c"))+
      ylab("Frequency"))
  }
  
p.facet.MOBF = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBF"))
p.facet.MOBH = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBH"))
p.facet.MOBFandH = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBH" | Majority.MOB=="MOBF"))

p.facet.MOBP = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBP"))
p.facet.all = facetPlot(antidefense.plot)
ggsave(p.all, file="../figures/Fig1/fig-facet-antidefense-all.pdf", width=6, height=4.5)

ggsave(p.facet.MOBF, file="../figures/Fig1/fig-facet-antidefense-MOBF.pdf", width=6, height=4)
ggsave(p.facet.MOBH, file="../figures/Fig1/fig-facet-antidefense-MOBH.pdf", width=6, height=4)
ggsave(p.facet.MOBFandH, file="../figures/Fig1/fig-facet-antidefense-MOBFandH.pdf", width=6, height=4) # to get ratio right

ggsave(p.facet.MOBP, file="../figures/Fig1/fig-facet-antidefense-MOBP.pdf", width=6, height=4)
# These are then combined and edited to make figure





# For each PTU, calculate the proportion of its plasmids 
# which carry an anti-defense system
plasmid_table$antidefense.detected = ifelse(rownames(plasmid_table) %in% unique(antidefense$plasmid),
                                                                                "yes", "no")
antidefense.PTU.count = plasmid_table %>% group_by(PTU) %>%
  summarise(n.antidefense=length(antidefense.detected[antidefense.detected=="yes"]),
            n.total=length(antidefense.detected)) %>%
  mutate(prop=n.antidefense/n.total * 100)

# Give genus counts for each plasmid
plasmid_table$Genus.simple = sapply(plasmid_table$Genus, 
                                    function(x) ifelse(x %in% head(names(sort(table(plasmid_table$Genus), 
                                                                             decreasing = TRUE)), n=6), x, "Other"))
table(plasmid_table$Genus.simple)
PTU.by.genus = plasmid_table %>% group_by(PTU, Genus.simple) %>%
  summarise(n=length(Genus.simple)) %>% mutate(prop=n/sum(n))
PTU.by.genus$Genus.simple = ordered(PTU.by.genus$Genus.simple,
                                    levels=c("Escherichia", "Shigella", "Salmonella",
                                             "Citrobacter", "Enterobacter", "Klebsiella", "Other"))
p.genus.props = ggplot(PTU.by.genus, aes(PTU, prop, fill=Genus.simple))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=c("#2B8CBE", "#A6BDDB", "#A6D854",  "#66C2A5","#FC8D62", "#E78AC3", "grey"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()
ggsave(p.genus.props, file='../figures/Fig1/fig-prop-genus-counts.pdf', width=4, height=4)


# Redo for all HMM hits evalue <1e-10 - just to check results are similar
# Read in defensefinder results
defensefinder = read.csv('../data/defensefinder_results/all_hmmer.tsv',
                         header=T,
                         sep='\t')
defensefinder = defensefinder[which(defensefinder$i_eval<1e-10),]
defensefinder$plasmid = sub("_[0-9]+$", "", defensefinder$hit_id)
defensefinder$position = as.numeric(gsub(".*_", "", defensefinder$hit_id))
# Add PTU
defensefinder$PTU = plasmid_table[defensefinder$plasmid, "PTU"]

# Read in genes, store which are defense/antidefense
df_genes= read.csv('../data/defensefinder_results/all_genes.tsv',
                   header=T,
                   sep='\t')
gene_names_antidefense = unique(df_genes$hit_gene_ref[which(df_genes$activity=="Antidefense")])
gene_names_defense = unique(df_genes$hit_gene_ref[which(df_genes$activity=="Defense")])

antidefense = defensefinder[which(defensefinder$gene_name %in% gene_names_antidefense),]
gene_names_antiRM = unique(df_genes$hit_gene_ref[which(df_genes$type=="Anti_RM")])
gene_names_antiCRISPR = unique(df_genes$hit_gene_ref[which(df_genes$type=="Anti_CRISPR")])
gene_names_psiAB = unique(df_genes$hit_gene_ref[which(df_genes$type=="Other")])

antidefense$type = sapply(antidefense$gene_name, 
                          function(x) ifelse(x %in% gene_names_antiRM, "Anti_RM",
                                             ifelse(x %in% gene_names_antiCRISPR, "Anti_CRISPR", 
                                                    ifelse(x %in% gene_names_psiAB, "psiAB", "Anti_Pycsar"))))
# Exclude anti-pycsar for plotting
antidefense.plot = antidefense[which(antidefense$type!="Anti_Pycsar"),]

p.all = ggplot(antidefense.plot, aes(position, fill=type))+
  geom_histogram()+theme_bw()+
  xlab("Start of system (ORF position)")+
  ylab("Frequency")+
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(colour="black"))+
  scale_fill_manual(values=c("#fecc5c","#fd8d3c","#e31a1c"))
