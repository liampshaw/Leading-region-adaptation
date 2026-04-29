library(data.table)
library(zoo)
library(ggplot2)
library(dplyr)
library(ggpubfigs)

# Create output directory for figures
dir.create("../figures/Fig2/", showWarnings = FALSE, recursive = TRUE)

# Plasmid table
plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T, 
                         row.names = 1)
# Subset 
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]

# GC content
gc = fread("../data/zenodo/seq_stats/avg_gc.tsv",
           col.names = c("plasmid_long", "position", "gc"))
gc$plasmid = gsub("_leading.fa", "", gc$plasmid_long)
gc$PTU = plasmid_table[gc$plasmid, "PTU"]

mob.types = read.csv('../data/PTU-mob-types.csv',
                     header=T)
gc = gc %>% left_join(mob.types, by="PTU")

# Mean GC content at each position for PTU
gc.by.PTU = gc %>% group_by(position, PTU, Majority.MOB) %>%
  summarise(
    mean.gc = mean(gc),
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(PTU) %>%
  mutate(prop = n / max(n)) %>%
  ungroup()

gc.by.PTU$position.kb = gc.by.PTU$position/1000
gcFacetPlot = function(df){
 return(ggplot(df, aes(position.kb, mean.gc))+
          geom_line(aes(alpha=prop))+
          facet_wrap(~PTU)+
          xlim(c(0,100))+
          xlab("Position (kb)")+
          ylab("Mean GC content")+
          theme_bw()+
          theme(panel.grid=element_blank())+
          theme(legend.position = "none")+
          ylim(c(0.35, 0.6))) 
}
p.gc.facet.MOBF = gcFacetPlot(gc.by.PTU[which(gc.by.PTU$Majority.MOB=="MOBF"),])
p.gc.facet.MOBP = gcFacetPlot(gc.by.PTU[which(gc.by.PTU$Majority.MOB=="MOBP"),])
p.gc.facet.MOBHandP = gcFacetPlot(gc.by.PTU[which(gc.by.PTU$Majority.MOB=="MOBH" | gc.by.PTU$Majority.MOB=="MOBP"),])
ggsave(p.gc.facet.MOBF, file="../figures/Fig2/fig-facet-GC-MOBF.pdf", width=6, height=4)
ggsave(p.gc.facet.MOBP, file="../figures/Fig2/fig-facet-GC-MOBP.pdf", width=6, height=4)
ggsave(p.gc.facet.MOBHandP, file="../figures/Fig2/fig-facet-GC-MOBHandP.pdf", width=6, height=6)
ggsave(p.gc.facet.MOBF+theme(legend.position = c(0.8, 0.8)), file="../figures/Fig2/fig-facet-GC-MOBF-with-legend.pdf", width=6, height=4)

# Average per MOB type
gc.by.MOB = gc.by.PTU %>% group_by(position.kb, Majority.MOB) %>%
  summarise(mean.gc=mean(mean.gc))
ggplot(gc.by.MOB, aes(position.kb, mean.gc))+
  geom_line(aes(colour=Majority.MOB))+  
  xlim(c(0,100))+
  xlab("Position (kb)")+
  ylab("Mean GC content")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position = "none")+
  ylim(c(0.35, 0.6))

# Take average GC for first 25kb vs. rest of plasmid
gc.by.PTU$leading.region = ifelse(gc.by.PTU$position.kb<26, "Leading 25kb",
                                  "Rest of plasmid")
gc.by.MOB = gc.by.PTU %>% group_by(PTU, leading.region, Majority.MOB) %>%
  summarise(mean.gc=mean(mean.gc)) %>%
  group_by(PTU, Majority.MOB) %>%
  mutate(diff=mean.gc[leading.region=="Leading 25kb"]-mean.gc[leading.region=="Rest of plasmid"]) %>%
  filter(leading.region=="Leading 25kb")
gc.by.MOB$Majority.MOB = ordered(gc.by.MOB$Majority.MOB,
                                 levels=c("MOBF", "MOBP", "MOBH"))
p.gc.by.MOB = ggplot(gc.by.MOB, aes(Majority.MOB, group=PTU, fill=Majority.MOB, y=diff))+
  geom_bar(stat="identity", width=1, position = position_dodge2(width = 0.8, preserve = "single"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=friendly_pal("contrast_three"))+
  xlab("")+
  ggrepel::geom_text_repel(aes(label=PTU), size=2)+
  ylab("Mean GC")+
  theme(legend.position = "none")
ggsave(p.gc.by.MOB, file="../figures/Fig2/figure-GC-by-MOB.pdf", width=4, height=2)

