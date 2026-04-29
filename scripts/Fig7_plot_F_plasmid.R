# Plot gggenes
library(stringr)
library(dplyr)
library(zoo)
library(ggplot2)
library(gggenes)
library(cowplot)
library(circlize)
library(ggsignif)


# Leading region of NC_002483.1 is from 66118 in upstream direction.

# Create output directory for figures
dir.create("../figures/Fig7/", showWarnings = FALSE, recursive = TRUE)

f_plasmid_genes =  read.csv('../data/zenodo/F-plasmid/F-plasmid-genes-NCBI-spreadsheet.csv', 
                          header=T)
# Add features
f_plasmid_features = read.csv('../data/zenodo/F-plasmid/F-plasmid-features.csv',
                              header=T)

# Plot the whole thig
p.cai = ggplot(f_plasmid_genes, aes((start/1000), cai, fill=cai))+
  scale_fill_continuous(low="white", high="red")+
  ylim(c(0.1,0.5))+
  geom_point(shape=21, size=3)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  xlab("Position (kb)")


# Plot the leading region
p.leading.region.closeup = ggplot(f_plasmid_genes %>% filter(pseudo==FALSE & strand=="+"), aes(fill=cai, xmin=start, xmax=end, ymin=0.95, ymax=1))+
  geom_rect()+
 scale_fill_continuous(low="white", high="red")+
  ylim(c(0,1.5))+
  geom_rect(data=f_plasmid_features, aes(xmin=start, xmax=end, ymin=0.95, ymax=1), fill="#0072B2")+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x=element_line()
  )+
  ylim(c(0.5,1.5))+
  xlim(c(52500,66410))



# We need to reorientate starting from CDS 68 and then going down 
# so we go in direction of transfer strand
f_plasmid_genes$new.position = (68-f_plasmid_genes$idx) %% 106 # is this right?
f_plasmid_genes$new.start = 66118 - f_plasmid_genes$start
f_plasmid_genes$new.end = 66118 - f_plasmid_genes$end

# Reorient so that if negative, they go positive (using length of F plasmid)
f_plasmid_genes$new.start[which(f_plasmid_genes$new.start<0)] = f_plasmid_genes$new.start[which(f_plasmid_genes$new.start<0)]+99159
f_plasmid_genes$new.end[which(f_plasmid_genes$new.end<0)] = f_plasmid_genes$new.end[which(f_plasmid_genes$new.end<0)]+99159

# get rid of pseudogenes
f_plasmid_genes = f_plasmid_genes %>% filter(pseudo==FALSE)
# add rolling average
cai.df <- f_plasmid_genes %>%
  arrange(new.start) %>%  
  mutate(cai_roll = zoo::rollmean(cai, k = 15, fill = NA, align = "center"))
cai.df = cai.df[!is.na(cai.df$cai_roll),] # remove missing values (start and end because of averaging)

p.cai.positive.strand = ggplot(f_plasmid_genes %>% filter(strand=="+"), aes((new.start/1000), cai, fill=cai))+
  scale_fill_continuous(low="white", high="red")
p.cai.positive.strand = p.cai.positive.strand + 
  ylim(c(0.1,0.5))+
  geom_point(shape=21, size=3)+
  #geom_line(data=cai.df, aes(new.start/1000, cai_roll), size=0.5, colour="red", size=2, alpha=0.5)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  xlab("Position (kb)")
p.cai.negative.strand = ggplot(f_plasmid_genes %>% filter(strand=="-"), aes((new.start/1000), cai, fill=cai))+
  scale_fill_continuous(low="white", high="red")
p.cai.negative.strand = p.cai.negative.strand + 
  ylim(c(0.1,0.5))+
  geom_point(shape=21, size=3)+
  #geom_line(data=cai.df, aes(new.start/1000, cai_roll), size=0.5, colour="red", size=2, alpha=0.5)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  xlab("Position (kb)")



# Add transfer region, oriT, leading region
df_rect_tra <- data.frame(
  xmin = 66568/1000,
  xmax = 98798/1000,
  ymin = -Inf,
  ymax = 0.1   # adjust as needed
)
df_rect_leading <- data.frame(
  xmin = 859/1000,
  xmax = 12334/1000,
  ymin = -Inf,
  ymax = 0.1   # adjust as needed
)
df_rect_oriT <- data.frame(
  xmin = -400/1000,
  xmax = 400/1000,
  ymin = -Inf,
  ymax = 0.12  # adjust as needed
)

p.cai.positive.strand = p.cai.positive.strand + geom_rect(data=df_rect_tra, 
                  aes(xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, 
                  fill="darkgrey")+
  geom_rect(data=df_rect_leading,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill="#009E73")+
  geom_rect(data=df_rect_oriT,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill="black")
# Compare tra genes to leading region
tra_genes = f_plasmid_genes[which(f_plasmid_genes$idx>69 & f_plasmid_genes$strand=="+"),]
leading_region_genes = f_plasmid_genes[which(f_plasmid_genes$idx<69 & f_plasmid_genes$idx>48 & f_plasmid_genes$strand=="+"),]
plot.type.df= data.frame(cai=c(tra_genes$cai, leading_region_genes$cai),
                         type=c(rep("Transfer region", length(tra_genes$cai)),
                                rep("Leading region", length(leading_region_genes$cai))))
plot.type.df$type = ordered(plot.type.df$type,
                            levels=c("Transfer region", "Leading region"))
p.boxplot = ggplot(plot.type.df, aes(type, cai, fill=cai))+
  geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA)+
  ggbeeswarm::geom_quasirandom(shape=21, size=3, alpha=0.8)+
  scale_fill_gradient(low = "white", high = "red", limits=c(min(f_plasmid_genes$cai), 
                                                            max(f_plasmid_genes$cai))) +
  theme_bw()+
  theme(panel.grid=element_blank())+
  xlab("")+
  theme(legend.position = "none")+
  ylim(c(0.1, 0.5))+
  ggsignif::geom_signif(comparisons=list(c("Leading region", "Transfer region")),
                          test="wilcox.test")
  

p.cai.combine = plot_grid(p.cai.positive.strand+
                            theme(legend.position = "none"), 
                          p.boxplot, rel_widths=c(3,1.7), nrow=1)



p.genes = ggplot(f_plasmid_genes %>% filter(pseudo==FALSE),
       aes(xmin = new.start/1000,
           xmax = new.end/1000,
           y = "",
           fill = cai,
           forward = strand == "+")) +
  geom_gene_arrow(arrowhead_height = unit(4, "mm"),
                  arrowhead_width = unit(1, "mm")) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_genes()+
  xlim(c(0,20))+
  ggrepel::geom_text_repel(aes(label=gene, x=new.start))+
  xlab("Position (kb)")+
  theme(legend.position = "none")




p.F.plasmid = plot_grid(p.cai.combine, p.leading.region.closeup, nrow=2)
ggsave(file='../figures/Fig7/fig7-F-plasmid-panel2.pdf', width=7, height=5, p.F.plasmid)

ggsave(file='../figures/Fig7/fig7-F-plasmid-with-legend.pdf', p.cai)



# Use circlize to plot genes in ring
tick_by = 10000
plasmid_length = 99159
cai_min <- min(f_plasmid_genes$cai)
cai_max <- max(f_plasmid_genes$cai)
col_fun <- colorRamp2(
c(cai_min, cai_max),
c( "white", "red")
)

pdf('../figures/Fig7/fig7-F-plasmid-circular-panel-1.pdf', width=5, height=5)
circos.clear()
plasmid_length = 99159
circos.par(start.degree = 90, gap.degree = 0,
           cell.padding = c(0,0,0,0),
           track.margin = c(0.01, 0.01))
circos.initialize(
  factors = "plasmid",
  xlim = c(0, plasmid_length)
)
circos.trackPlotRegion(
  ylim = c(0, 2),
  track.height = 0.3,
  bg.border = NA
)

circos.axis(
  h = "bottom",
  major.at = seq(0, 100000, by = tick_by),
  labels = paste0(seq(0, 100000, by = tick_by) / 1000, " kb"),
  labels.cex = 0.6,
  major.tick.length = -0.2,
  lwd = 1
)

for(i in seq_len(nrow(f_plasmid_genes))) {
  ybottom <- ifelse(f_plasmid_genes$strand[i] == "+", 0.55, 0.05)
  ytop    <- ifelse(f_plasmid_genes$strand[i] == "+", 0.95, 0.45)
  circos.rect(
    f_plasmid_genes$start[i],
    ybottom,
    f_plasmid_genes$end[i],
    ytop,
    col =col_fun(f_plasmid_genes$cai[i]) ,
    border = NA
  )
}
# Transfer region annotaitno
annot_start <- 66500
annot_end   <- 98709
annot_label <- "Transfer region"
circos.rect(
  xleft = annot_start,
  ybottom = 1.5,
  xright = annot_end,
  ytop = 1.7,
  col = "darkgrey",
  border = NA
)
# oriT
annot_start <- 66407	
annot_end   <- 66118
annot_label <- "oriT"
circos.rect(
  xleft = annot_start,
  ybottom = 1.2,
  xright = annot_end,
  ytop = 2,
  col = "black",
  border = NA
)
# Leading region annotaion
annot_start <- 53784
annot_end   <- 66000
annot_label <- "Leading region"
circos.rect(
  xleft = annot_start,
  ybottom = 1.5,
  xright = annot_end,
  ytop = 1.7,
  col = "#009E73",
  border = NA
)

dev.off()

