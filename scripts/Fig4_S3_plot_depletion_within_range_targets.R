library(ggplot2)

# Create output directory for figures
dir.create("../figures/Fig4/", showWarnings = FALSE, recursive = TRUE)
dir.create("../figures/Supp/", showWarnings = FALSE, recursive = TRUE)

# Plasmid table
plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T,
                         row.names = 1)
# Subset to leading region plasimds
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]

# Mob types
mob.types = read.csv('../data/PTU-mob-types.csv',
                     header=T)

scores_average_all = NULL
for (PTU in unique(plasmid_table$PTU)){
  for (k in c(5,6)){
    scores_df = read.csv(paste0('../data/zenodo/targets_leading_lagging/', PTU, '_k', k, '_scores.csv'), header=T)
    n_targets = nrow(scores_df)
    scores_average <- scores_df %>%
      rowwise() %>%
      mutate(score_avg = mean(c_across(starts_with("score"))),
             sd=sd(c_across(starts_with("score")))) %>%
      ungroup() %>%
      select(kmer, score_avg, sd)
    counts_df = read.csv(paste0('../data/zenodo/targets_leading_lagging//', PTU, '_k', k, '_counts.csv'), header=T)
    counts_average <- counts_df %>%
      rowwise() %>%
      mutate(count_avg = mean(c_across(starts_with("count")))) %>%
      ungroup() %>%
      select(kmer, count_avg)
    scores_average = scores_average %>% left_join(counts_average, by="kmer")
    scores_average$PTU = PTU
    scores_average$k = nchar(scores_average$kmer)
    scores_average_all = rbind(scores_average_all, scores_average)
      }
}

scores_average_all = data.frame(scores_average_all)
scores_average_all = scores_average_all %>% left_join(mob.types, by="PTU")
scores_average_all$Majority.MOB = ordered(scores_average_all$Majority.MOB,
                                          levels=c("MOBF", "MOBP", "MOBH"))
p.k.5 = ggplot(scores_average_all %>% filter(k==5), aes(Majority.MOB, group=PTU, score_avg, fill=Majority.MOB))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_boxplot(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=friendly_pal("contrast_three"))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid = element_blank())+
  ylab("Relative representation score\n(leading 10kb vs. lagging 10kb)")+
  xlab("")
p.k.6 = ggplot(scores_average_all %>% filter(k==6), aes(Majority.MOB, group=PTU, score_avg, fill=Majority.MOB))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_boxplot(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=friendly_pal("contrast_three"))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid = element_blank())+
  ylab("Relative representation score\n(leading 10kb vs. lagging 10kb)")+
  xlab("")

p.targets.combined = plot_grid(p.k.5+ggtitle("(a) 5-bp within-range targets"), 
          p.k.6+ggtitle("(b) 6-bp within-range targets"), nrow=1)
ggsave(file="../figures/Fig4/fig4-within-range-targets.pdf", width=10, height=4, p.targets.combined)



#colnames(scores_average_all) = c("PTU", "k", "n.targets", "median.score", "median.count")
scores_average_all$median.score= as.numeric(scores_average_all$median.score)
#scores_average_all$median.count = as.numeric(results$median.count)

# Look at mean score for 6-bp motifs across PTUs
scores_average_k6 = scores_average_all %>% filter(k==6)
scores_average_k6_means = scores_average_k6 %>% group_by(kmer) %>%
  summarise(mean=mean(score_avg),
            median=median(score_avg))

scores_average_k6_means$kmer = ordered(scores_average_k6_means$kmer,
                                       levels=scores_average_k6_means$kmer[order(scores_average_k6_means$mean)])


p.average.score.per.target = ggplot(scores_average_k6_means, aes(kmer, mean))+
  geom_bar(stat="identity")+
  geom_hline(yintercept = mean(scores_average_k6_means$mean), colour='red', linetype='dashed')+
  theme_bw()+
  xlab("6-bp RM target")+
  ylab("Mean difference in representation score across all plasmids\n (leading 10kb vs. lagging 10kb)")+
  coord_flip()
ggsave(file="../figures/Supp/FigS4-within-range-targets.pdf", width=8, height=8, p.average.score.per.target)
