library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Create output directory for figures
dir.create("../figures/Fig5/", showWarnings = FALSE, recursive = TRUE)
dir.create("../figures/Fig7/", showWarnings = FALSE, recursive = TRUE)
dir.create("../figures/Supp/", showWarnings = FALSE, recursive = TRUE)

# Read in CAI data
cai = read.csv('../data/cai_prodigal_meta.tsv',
               header=T,
               sep='\t')
defensefinder = read.csv('../data/TableSX_defensefinder_results_prodigal_meta.tsv',
                         header=T,
                         sep='\t')
plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T,
                         row.names = 1)
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]

cai$id = gsub(".*\\|", "", cai$gene_id)
cai$plasmid = gsub("\\.fna.*", "", gsub(".*\\/", "", cai$gene_id))
cai$position = as.numeric(gsub(".*_", "", cai$gene_id))
cai$PTU = plasmid_table[cai$plasmid, "PTU"]
cai$HostGenus = plasmid_table[cai$plasmid, "Genus"]
cai$HostSpecies = plasmid_table[cai$plasmid, "Species"]

# Mob types
mob.types = read.csv('../data/PTU-mob-types.csv', header=T)



cai.ptu = cai %>% group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai))
# Add a length threshold for each PTU
cai.ptu.thresholds = cai %>%
  group_by(PTU, plasmid) %>%
  summarise(max=max(position)) %>%
  group_by(PTU) %>%
  summarise(percentile=quantile(max, 0.25)) 
cai.ptu = cai.ptu %>% left_join(cai.ptu.thresholds, by="PTU")

cai.ptu = cai.ptu %>% left_join(mob.types, by="PTU")



pdf('../figures/Fig5/fig5-cai_results.pdf')
ggplot(cai.ptu %>% filter(position <= percentile), aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU, scales="free_x")+
  scale_colour_continuous(low="white", high="red")+
  xlab("ORF position")+
  ylab("Mean CAI")+
  theme_bw()+
  theme(panel.grid=element_blank())
dev.off()

plotCAI = function(df){
  return(ggplot(df, aes(position, mean.cai))+
           geom_point(size=1, aes(colour=mean.cai))+
           facet_wrap(~PTU, scales="free_x")+
           scale_colour_continuous(low="white", high="red", limits=c(0.15, 0.425))+
           xlab("ORF position")+
           ylab("Mean CAI")+
           theme_bw()+
           theme(panel.grid=element_blank()))
}

empty_plot <- ggplot() +
  geom_blank() +
  theme_minimal()

# Plot by MOB type
cai.ptu.plot = cai.ptu[which(cai.ptu$position<100),]
p.F.1 = plotCAI(cai.ptu.plot[which(cai.ptu.plot$PTU %in% c("PTU-E81", "PTU-E84", "PTU-FE")),]) +
  theme(legend.position = "none")
p.F.2 = plotCAI(cai.ptu.plot[which(cai.ptu.plot$PTU %in% c("PTU-FE", "PTU-N1")),]) +
  theme(legend.position = "none")
p.P.1 = plotCAI(cai.ptu.plot[which(cai.ptu.plot$PTU %in% c("PTU-BOKZ", "PTU-I1", "PTU-I2")),]) +
  theme(legend.position = "none")+
  ylab("")+
  xlab("")
p.P.2 = plotCAI(cai.ptu.plot[which(cai.ptu.plot$PTU %in% c("PTU-LM", "PTU-X1", "PTU-X3")),]) +
  theme(legend.position = "none")+
  ylab("")+
  xlab("")
p.H = plotCAI(cai.ptu.plot[which(cai.ptu.plot$PTU %in% c("PTU-C", "PTU-HI1A")),]) +
  theme(legend.position = "none")

p.F = plotCAI(cai.ptu.plot[which(cai.ptu.plot$Majority.MOB=="MOBF"),]) +
  theme(legend.position = "none")+
  ylim(c(0.15, 0.45))
p.P = plotCAI(cai.ptu.plot[which(cai.ptu.plot$Majority.MOB=="MOBP"),]) +
  theme(legend.position = "none")+
  ylab("")+
  ylim(c(0.15, 0.45))
# Missing data for MOBH to preserve facet structure
df.H <- cai.ptu.plot[cai.ptu.plot$PTU %in% c("PTU-C", "PTU-HI1A"), ]
# Make sure the faceting variable includes a third level for the missing panel
df.H$PTU <- factor(df.H$PTU, levels = c("PTU-C", "PTU-HI1A", "MISSING", "MISSING2", "MISSING3", "MISSING4"))

p.H <- plotCAI(df.H) +
  theme(legend.position = c(0.8, 0.2)) +
  ylab("") +
  ylim(c(0.15, 0.45)) +
  facet_wrap(~PTU, drop = FALSE, scales="free_x")
p.CAI.combined = plot_grid(p.F, p.P, p.H, axis='bt',align='h', nrow=1,
          rel_widths = c(3,3.2,3.2))
ggsave(file="../figures/Fig5/fig5-CAI-combined.pdf", p.CAI.combined, width=12, height=4)

# And average over all PTUs at each position
p.average.over.all = cai.ptu %>% group_by(position) %>%
  summarise(mean.cai=mean(mean.cai)) %>%
  ggplot(aes(position, mean.cai))+
  geom_point(size=1, aes(colour=mean.cai))+
  scale_colour_continuous(low="white", high="red", limits=c(0.15, 0.425))+
  xlab("ORF position")+
  ylab("Mean CAI")+
  theme_bw()+
  xlim(c(0,100))+
  theme(panel.grid=element_blank())+
  theme(legend.position = "none")
ggsave(file="../figures/Fig5/fig5-CAI-average-over-PTUs.pdf", width=3, height=1.8)

# Flip them as well
cai.flipped <- cai %>%
  group_by(plasmid) %>%
  mutate(position = -(max(position) + 1 - position))%>% # 1 -> -N, N -> -1
  ungroup()
cai.both = rbind(cai, cai.flipped)
cai.both.max.positions = cai.both %>%
  group_by(PTU, plasmid) %>%
  summarise(max=max(position)) %>%
  group_by(PTU) %>%
  summarise(percentile=quantile(max, 0.25)) # Choose a percentile threshold e.g. 0.2 means at this length 75% of plasmids in the PTU are at least this length
cai.ptu.both = cai.both %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai)) 
cai.ptu.both = cai.ptu.both %>%
  left_join(cai.both.max.positions, by = "PTU")

cai.ptu.both.filter = cai.ptu.both %>%
  filter(position <= percentile/2 & position >= -percentile/2)



# Make plots for 
MOBF.PLASMIDS = c("PTU-N1", "PTU-FE", "PTU-FS")
MOBP.PLASMIDS = c("PTU-I1", "PTU-I2", "PTU-BOKZ", "PTU-LM", 
                  "PTU-X1", "PTU-X3", "PTU-E81", "PTU-E84")
MOBH.PLASMIDS = c("PTU-HI1A", "PTU-C")



p.mobF = ggplot(cai.ptu.both.filter %>% filter(PTU %in% MOBF.PLASMIDS), aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU, scales="free_x", nrow=1)+
  scale_colour_continuous(low="white", high="red")
p.mobP = ggplot(cai.ptu.both.filter %>% filter(PTU %in% MOBP.PLASMIDS), aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU, nrow=1,scales="free_x")+
  scale_colour_continuous(low="white", high="red")
p.mobH = ggplot(cai.ptu.both.filter %>% filter(PTU %in% MOBH.PLASMIDS), aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU,nrow=1, scales="free_x")+
  scale_colour_continuous(low="white", high="red")


##############################
# Just Escherichia coli plasmids
##############################
cai.ecoli = cai %>% filter(HostSpecies=="Escherichia coli")
cai.ptu.ecoli = cai.ecoli %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai))



# Try inverting so that we have the lagging region going from -1 to -N
cai.ecoli.flipped <- cai.ecoli %>%
  group_by(plasmid) %>%
  mutate(position = -(max(position) + 1 - position))%>%
  ungroup()
cai.ptu.ecoli.flipped = cai.ecoli.flipped %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai))

cai.ecoli.both = rbind(cai.ecoli.flipped, cai.ecoli)
cai.ptu.ecoli.both = cai.ecoli.both %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai)) 



# Check maximum genes for each PTU and divide by two
cai.ecoli.both.max.positions = cai.ecoli.both %>%
  group_by(PTU, plasmid) %>%
  summarise(max=max(position)) %>%
  group_by(PTU) %>%
  summarise(percentile=quantile(max, 0.25)) # Choose a percentile threshold e.g. 0.2 means at this length 75% of plasmids in the PTU are at least this length
cai.ptu.ecoli.both = cai.ptu.ecoli.both %>%
  left_join(cai.ecoli.both.max.positions, by = "PTU")
# Filter out 
cai.ptu.ecoli.both.constrained = cai.ptu.ecoli.both %>%
  filter(position <= percentile/2 & position >= -percentile/2)

ggplot(cai.ptu.ecoli.both.constrained, aes(position, mean.cai))+
  geom_point(size=0.5)+
  facet_wrap(~PTU, scales="free_x")


## FIGURE 7
# Analysing defensefinder results
antidefense = defensefinder[which(defensefinder$activity=="Antidefense"),]
defense = defensefinder[which(defensefinder$activity=="Defense"),]

antidefense_proteins <- antidefense %>%
  separate_rows(protein_in_syst, sep = ",") %>%
  pull(protein_in_syst) %>%
  unique()
defense_proteins <- defense %>%
  separate_rows(protein_in_syst, sep = ",") %>%
  pull(protein_in_syst) %>%
  unique()
cai <- cai %>%
  mutate(protein_type = if_else(id %in% antidefense_proteins, "antidefense", 
                               ifelse(id %in% defense_proteins, "defense", "other")),
         leading_region = if_else(position < 20, " In first 25 ORFs", "Beyond first 25 ORFs"))


# Plot of antidefense genes by position, using sliding window
cai.rolling.antidefense = cai %>%filter(protein_type=="antidefense") %>%
  arrange(position) %>% 
  mutate(CAI_roll = zoo::rollmean(cai, k = 15, fill = NA, align = "center")) %>%
  group_by(position) %>% 
  summarise(mean_cai=mean(CAI_roll))
cai.rolling.other = cai %>%filter(protein_type=="other") %>%
  arrange(position) %>% 
  mutate(CAI_roll = zoo::rollmean(cai, k = 15, fill = NA, align = "center")) %>%
  group_by(position) %>% 
  summarise(mean_cai=mean(CAI_roll))
p.cai.rolling = ggplot(cai.rolling.antidefense, aes(position, mean_cai))+
  geom_line(colour="red")+
  geom_point(shape=21, fill="red")+
  xlim(c(0,100))+
  ylim(c(0.25,0.47))+
  theme_bw()+
  theme(panel.grid = element_blank())
p.cai.rolling = p.cai.rolling + 
  geom_point(data=cai.rolling.other)+
  geom_line(data=cai.rolling.other)+
  xlim(c(0,100))+
  ylim(c(0.25,0.47))+
  xlab("ORF position")+
  ylab("Mean CAI")

# Add the types of gene
defensefinder.expanded = defensefinder %>%
  separate_rows(protein_in_syst, sep = ",") 
defensefinder.expanded$id = defensefinder.expanded$protein_in_syst # map the IDs
# Add in tra genes (from hmmsearch against CONJScan hmms from MacSyfinder)
conj.hits = read.csv('/Users/Liam/Downloads/new_1751/prodigal-output-meta/all_results.tsv', sep='\t', header=T)
# keep those with e-value <1e-10
conj.hits.e10 = conj.hits[which(conj.hits$evalue<1e-10),]
conj.hits.e10$id = conj.hits.e10$target
conj.hits.e10$type = "CONJScan"
conj.hits.e10$subtype = conj.hits.e10$query
conj.hits.e10$activity = "CONJScan"

# Combine defensefinder and Conj hits
info.to.add = rbind(defensefinder.expanded[,c("id", "type", "subtype", "activity")],
                    conj.hits.e10[,c("id", "type", "subtype", "activity")])

cai.combined = cai %>% left_join(info.to.add, by="id") # used to be defensefiner.expanded
# Assign non defense/antidefense to "Other" category
# Look at the tables
table(cai.combined$activity)
sort(table(cai.combined$type[which(cai.combined$activity=="Defense")]))

# Add some more detailed subtypes
cai.combined$type.for.plot = cai.combined$type
cai.combined$type.for.plot[which(cai.combined$activity=="Defense")] = sapply(cai.combined$type.for.plot[which(cai.combined$activity=="Defense")],
                                           function(x) 
                                             ifelse(grepl("CBASS", x), "CBASS", ifelse(grep("RM", x), "RM", "Other defense")))
cai.combined$type.for.plot[which(cai.combined$activity=="Defense" & is.na(cai.combined$type.for.plot))] = "Other defense" # hacky...
# Change "Other" anti-defense category  to psiAB (after inspection, makes up all the cases)
cai.combined$type.for.plot[which(cai.combined$activity=="Antidefense" & cai.combined$type.for.plot=="Other")] = "psiAB"
# Add Other
# Add Unknown category
cai.combined$activity[is.na(cai.combined$activity)] = "Other"
cai.combined$type.for.plot[which(cai.combined$activity=="Other")] = "Unknown"

# Order
cai.combined$type.for.plot = ordered(cai.combined$type.for.plot, 
                                     levels=c("Anti_CRISPR", "Anti_RM", "psiAB", "RM", "CBASS", "Other defense", "Unknown"))
# Drop the anti-Pycsar
cai.combined = cai.combined[which(cai.combined$id!="NZ_CP023144.1_94"),]
# Compute medians
medians <- cai.combined %>%
  group_by(activity, type.for.plot) %>%
  summarise(median_cai = median(cai, na.rm = TRUE), .groups = "drop")
pd <- position_dodge2(width = 0.75, preserve = "single")

# Ecoli K-12 optimal gene set CAIs
ecoli.optimal.set.of.genes = read.csv('../data/Ecoli-K12-optimal-weights.txt', header=F)
colnames(ecoli.optimal.set.of.genes) = "cai"
ecoli.optimal.set.of.genes$activity = " Highly expressed\nchromosomal genes"
ecoli.optimal.set.of.genes$type.for.plot = "Unknown"
p.fig.7.antidefense.types = ggplot(cai.combined, aes(activity, group=type.for.plot, cai, fill=type.for.plot))+
  geom_boxplot(position=pd)+
  scale_fill_manual(values=c("#fecc5c","#fd8d3c","#e31a1c", "#ece2f0", "#1c9099", "#a6bddb",  "grey"))+
  theme_bw()+
  geom_text(data = medians,
            aes( y=median_cai+0.013, label = round(median_cai, 2)),
            size = 3,
            position=pd,
            inherit.aes = TRUE) +
  geom_boxplot(data=ecoli.optimal.set.of.genes)+
  ylab("CAI")+
  xlab("")+
  labs(fill="Category")
ggsave(p.fig.7.antidefense.types, file='../figures/Fig7/figure-7-2026-gene-types.pdf', 
       width=6, height=4)
