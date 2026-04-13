# Plotting palindrome depletion
library(data.table)
library(zoo)
library(ggplot2)
library(dplyr)
library(ggpubfigs)

# Create output directory for figures
dir.create("../figures/Fig3/", showWarnings = FALSE, recursive = TRUE)

# Plasmid table
plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T,
                         row.names = 1)
# Subset to leading region plasimds
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]

theoreticalPalindromeExpectation = function(p, k){# p is GC content, k is k-mer length
  return ((p**2 - p + 1/2)**(k/2))
}

# GC content
gc = fread("../data/seq_stats/avg_gc.tsv",
           col.names = c("plasmid_long", "position", "gc"))
gc$plasmid = gsub("_leading.fa", "", gc$plasmid_long)
gc$PTU = plasmid_table[gc$plasmid, "PTU"]

mob.types = read.csv('../data/PTU-mob-types.csv',
                     header=T)
# y-limits for plots
# 4: -30, 20
# 6: -15, 5
# 8: -4, 3
Y.LIMITS = list("4"=c(-30, 20), 
                "6"=c(-15, 5),
                "8"=c(-4,3))


palindromeFacetPlot = function(df){
  return(ggplot(df, aes(position.kb, mean.depletion.per.kb))+
           geom_hline(linetype='dashed', yintercept = 0)+
           facet_wrap(~PTU)+
           geom_line(aes(alpha=prop))+
           xlim(c(0,100))+
           theme_bw()+
           theme(panel.grid = element_blank())+
           ylab("Palindrome depletion")+
           theme(legend.position = "none"))
}

makeDepletionPlot = function(input_file, output_plot_filename, k){
  pal.dt = fread(input_file,
                 col.names=c("plasmid_long", "position", "density"))
  pal.dt$plasmid = gsub("_leading.fa", "", gsub(".*\\/", "", pal.dt$plasmid_long))
  pal.dt$PTU = plasmid_table[pal.dt$plasmid, "PTU"]
  
  # Add GC content
  pal.dt[gc, on = .(plasmid, position), gc := i.gc]
  # Calculate theoretical density
  pal.dt$theoretical.density = theoreticalPalindromeExpectation(pal.dt$gc, k) 
  pal.dt$depletion = pal.dt$density-pal.dt$theoretical.density
  
  dt.by.PTU = pal.dt %>% group_by(position, PTU) %>%
    summarise(mean.depletion.per.kb=mean(depletion) * 1000,    n = n(),
              .groups = "drop"
    ) %>%
    group_by(PTU) %>%
    mutate(prop = n / max(n)) %>%
    ungroup()
  
  dt.by.PTU = dt.by.PTU %>% left_join(mob.types, by="PTU")
  dt.by.PTU$Majority.MOB = ordered(dt.by.PTU$Majority.MOB, 
                                   levels=c("MOBF", "MOBP", "MOBH"))
  # Sort out for plotting
  dt.by.PTU$position.kb = dt.by.PTU$position/1000
  
  p.F.1 = palindromeFacetPlot(dt.by.PTU[which(dt.by.PTU$PTU %in% c("PTU-E81", "PTU-E84", "PTU-FE")),])+
    ylim(Y.LIMITS[[as.character(k)]])+
    facet_wrap(~PTU, nrow=1)
  p.F.2 = palindromeFacetPlot(dt.by.PTU[which(dt.by.PTU$PTU %in% c("PTU-FS", "PTU-N1")),])+
    ylim(Y.LIMITS[[as.character(k)]])+
    facet_wrap(~PTU, nrow=1)
  p.P.1 = palindromeFacetPlot(dt.by.PTU[which(dt.by.PTU$PTU %in% c("PTU-BOKZ", "PTU-I1", "PTU-I2")),])+
    ylim(Y.LIMITS[[as.character(k)]])+
    facet_wrap(~PTU, nrow=1)
  p.P.2 = palindromeFacetPlot(dt.by.PTU[which(dt.by.PTU$PTU %in% c("PTU-LM", "PTU-X1", "PTU-X3")),])+
    ylim(Y.LIMITS[[as.character(k)]])+
    facet_wrap(~PTU, nrow=1)
  
  p.H = palindromeFacetPlot(dt.by.PTU[which(dt.by.PTU$Majority.MOB=="MOBH"),])+
    ylim(Y.LIMITS[[as.character(k)]])+
    facet_wrap(~PTU, nrow=1)
  empty_plot <- ggplot() + 
    theme_void()
  
  p.top = cowplot::plot_grid(p.F.1, p.P.1, p.H, nrow=1, rel_widths=c(3,3,2))
  p.bottom = cowplot::plot_grid(cowplot::plot_grid(p.F.2, empty_plot, rel_widths=c(2.5,1)),
                                p.P.2, empty_plot, nrow=1, rel_widths=c(3,3,2))
  p.combined = cowplot::plot_grid(p.top, p.bottom, nrow=2)
  ggsave(file=output_plot_filename, width=11, height=4, p.combined)
}

# Make figures for  k=6 (used in Fig 3) and also k=4 (Fig S1) and k=8 (Fig S2)
makeDepletionPlot("../data/seq_stats/avg_4_step1000_window1000.tsv", 
                  "../figures/Fig3/figS1-palindrome-depletion-k4.pdf",
                  4)
makeDepletionPlot("../data/seq_stats/avg_6_step1000_window1000.tsv", 
                  "../figures/Fig3/fig-palindrome-depletion-k6.pdf",
                  6)
makeDepletionPlot("../data/seq_stats/avg_8_step1000_window1000.tsv", 
                  "../figures/Fig3/figS2-palindrome-depletion-k8.pdf",
                  8)