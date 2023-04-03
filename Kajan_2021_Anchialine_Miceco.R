#R codes used to generate plots in Kajan et al. 2021

#Insert data----



#Fig.2. Hydrographical profile in a depth profile of anchialine cave-----
library(tidyverse)
library(ggplot2)
library(ggpubr)


dataVert_final_n <- dataVert_final %>% 
  select(Site, Depth, Salinity, Temperature, pH, DO)

env_df<-dataVert_final_n %>% gather("index","value",
                                    -Site, -Depth,factor_key=T)


VGvertical <- subset(env_df, Site=="VG")
BPvertical <- subset(env_df, Site=="BP")
ZVPvertical <- subset(env_df, Site=="ZVP")
GKVvertical <- subset(env_df, Site=="GKV")

zone_data <- tibble(ymin = 4, ymax = 6.5, xmin = -Inf, xmax = Inf)

all_plot_VG <- ggplot(VGvertical, aes(x = value, y = Depth)) +
  geom_lineh(size=0.4, color="gray") +
  geom_point(size=1, alpha = 0.7) +
  scale_y_reverse(breaks = seq(0, 25, by =5)) +
  facet_geochem_gridh(vars(index), scales = "free") +
  labs(x = NULL, y = "Depth (m)") +
  theme(axis.text.x = element_text(angle=90,vjust=0.3),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        panel.spacing = unit(1, "lines"),
        axis.line = element_line(color="black"),
        axis.title = element_text(size=8),
        panel.grid.major = element_line(size = 0.3, colour="gray",linetype="dashed"),
        plot.title = element_text(size = 8, hjust = 0.5, face="bold"),
        strip.background = element_rect( fill="white"))+
  ggtitle("Cave VG")


all_plot_VG_2 <- all_plot_VG+geom_rect(
  mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  data = zone_data, 
  alpha = 0.2,
  fill = "steelblue4",
  inherit.aes = FALSE
)



grid.arrange(all_plot_VG_2 , all_plot_BP_2, all_plot_GKV_2, all_plot_ZVP_2,nrow=4)

Ver_plot_grid <- ggarrange(all_plot_VG_2 , all_plot_BP_2, all_plot_GKV_2, all_plot_ZVP_2,
                           labels = c("a","b", "c","d"), font.label = list(size = 12, color = "black", face = "plain", family = NULL),
                           hjust = -0.5,vjust= 1.45,
                           nrow = 4)

ggsave(filename = "Vertical_parameters_AC.tiff", plot= Ver_plot_grid  ,  width = 10, height = 25, units = "cm")








#Fig.5.a Taxonomic composition-----
phylumRA_16S<-phylumCounts_16S/(apply(otu_table_16S,1,sum))*100
phylumRA_AB_16S<-round(as.data.frame(phylumRA_16S[,apply(phylumRA_16S,2,max)>=1]), digits=1)

phylumRA_AB_16S_other<-phylumRA_B[,apply(phylumRA_16S,2,max)<1]
Others <- round((as.data.frame(apply(phylumRA_AB_16S_other,1,sum))), digits=1)
colnames(Others) <- paste("Others")

phylumRA_AB_16S$Others <- Others$Others
phylum_all <- as.data.frame(t(phylumRA_AB_16S))

phylum_all <- phylum_all[c("Acidobacteria", "Actinobacteria","Bacteroidetes", "Calditrichaeota",
                           "Chloroflexi","Cloacimonetes","Epsilonbacteraeota",
                           "Firmicutes","Gemmatimonadetes","Latescibacteria","Magnetococcia",
                           "Marinimicrobia (SAR406 clade)",
                           "Nitrospirae","Omnitrophicaeota","Patescibacteria","Planctomycetes",
                           "Alphaproteobacteria","Deltaproteobacteria","Gammaproteobacteria",
                           "Rokubacteria","Verrucomicrobia","Others"),]


phylum_all$ID <- rownames(phylum_all)
phylum_all.melt <- melt(phylum_all, id.var = 'ID')
phylum_all.melt$ID <- factor(phylum_all$ID, levels = phylum_all$ID)



palette1=c("coral","coral4","darkgoldenrod1","cadetblue","cadetblue1","gray39")



png("Relative_abundance_phylum_16S.png",1400,900)
ggplot(phylum_all.melt, aes(x = variable, y = value, fill = ID)) +
  geom_bar(aes(), stat = 'identity', position = "stack") +
  theme(legend.title=element_blank(),
        axis.line = element_line(size=0.4),
        axis.text = element_text(size=25),
        axis.title = element_text(size=25),
        plot.title = element_text(size = 35),
        legend.text = element_text(color = "black", size = 25,margin = margin(t = 5)),
        legend.position = 'bottom', 
        legend.spacing.x = unit(1, 'cm'),
        legend.key.size = unit(1, "cm"),
        legend.key.height=unit(1, "cm")) + 
  guides(fill=guide_legend(ncol=6,nrow=2, byrow=TRUE)) +
  scale_fill_manual(values=palette1)+
  #scale_fill_hue(palette1) +
  scale_y_continuous(limits=c(0,100.2), expand=c(0,0)) + 
  ylab("Percentage of classified OTUs (%)") +
  xlab("Samples")

dev.off()




#Fig.5.b PCoA--------
otu_table_16S.pc<-pcoa(vegdist(decostand(otu_table_16St,"hellinger"),"bray"))

plot(hclust(vegdist(decostand(otu_table_16St,"hellinger"),"bray")))

otu_table_16Spc<-data.frame(otu_table_16S.pc$vectors)

pc1_16S=round(otu_table_16S.pc$values$Relative_eig[1], 3)*100
pc2_16S=round(otu_table_16S.pc$values$Relative_eig[2], 3)*100

otu_table_16Spc$cave<-as.factor(rep(c("VG","BP","GKV","ZVP"), c(3,3,3,3)))
otu_table_16Spc$halo<-as.factor(c("AH","H","BH","AH","H","BH",
                                "AH","H","BH","AH","H","BH"))


otu_table_16S_df <- otu_table_16Spc[,c("Axis.1", "Axis.2", "cave")]
otu_table_16S_find_hull <- function(otu_table_16Spc) otu_table_16Spc[chull(otu_table_16Spc$Axis.1, otu_table_16Spc$Axis.2), ]
otu_table_16S_hulls <- ddply(otu_table_16S_df, "cave", otu_table_16S_find_hull)




p1_16S <- ggplot(data=otu_table_16Spc, aes(x= Axis.1 * -1 , y= Axis.2 * -1))+ 
  geom_point(aes(shape=cave, colour=halo), cex=1, stroke=1.2,alpha = 0.8) + 
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  #geom_polygon(data=otu_table_16S_hulls, alpha=.1, aes(group=cave), fill="royalblue4")+
  labs(x = sprintf("PCo1 (%.1f%%)", pc1_16S), y = sprintf("PCo2 (%.1f%%)", pc2_16S))+
  scale_colour_manual(values=c("cadetblue3", "gray", "blue"), 
                      breaks=c("AH","H","BH"))+
  scale_shape_manual(values=c(0:3),
                     breaks=c("VG","BP","GKV","ZVP"))+
  geom_mark_ellipse(aes(label = cave, group=cave), linetype = 1, fill=c("#9AA490", "#9AA490","#9AA490",
                                                                                 '#AA643E','#AA643E','#AA643E',
                                                                                 '#B9A984','#B9A984','#B9A984',
                                                                                 '#205858', '#205858','#205858'),
                                                                                 label.fontsize = 6, expand = .025, label.fill=NA,
                    label.buffer = unit(5, 'mm'))+
  theme_bw()+
  guides(color=guide_legend(title= "Salinity gradient",order=1), 
         shape=guide_legend(title= "Cave", order=2))+
  theme(legend.title = element_text(colour="black", size=7, face="bold"),
        legend.text = element_text(colour="black", size=7, face="bold"),
        axis.text=element_text(size=7,face="bold",colour="black"),
        axis.title=element_text(size=7,face="bold", vjust=1, colour="black"),
        legend.position = "bottom")+
  scale_y_continuous(breaks=c(-.3,-.2,-.1,0,.1,.2,.3,.4), limits = c(-.35, .35))+
  scale_x_continuous(breaks=c(-.3,-.2,-.1,0,.1,.2,.3), limits = c(-.31, .31))+
  guides(col = guide_legend(label.position="right",title= "Sampling point",
                            title.position = "top",override.aes = list(size=4)),  shape=FALSE)+
  labs(tag = "b")








#Fig.5.c Venn diagram------
library(VennDiagram)

group.OTUsZVP <- rownames(otu_table_16S[which(rowSums(otu_table_16S[,1:3])>0),])
group.OTUsVG <- rownames(otu_table_16S[which(rowSums(otu_table_16S[,4:6])>0),])
group.OTUsGKV <- rownames(otu_table_16S[which(rowSums(otu_table_16S[,7:9])>0),])
group.OTUsBP <- rownames(otu_table_16S[which(rowSums(otu_table_16S[,10:12])>0),])

venn.diagram(list(ZVP = group.OTUsZVP,  GKV = group.OTUsGKV, VG = group.OTUsVG, BP = group.OTUsBP  ),
             fill = c("coral","coral4","cadetblue","lightblue"),
             cex = 1.5,filename="AC_Bacteria_Venn_all.png")









#Fig.6. Bubble plot------
abun_16S<-read.table("16S_abundant_AC_final_R.txt", header=T, sep="\t")

abund_16S_df <- abun_16S %>% 
  gather(sample, value, VG1:ZVP3)

abund_16S_dfF <- data.frame(abund_16S_df,
                           cave = rep(c("VG","BP","GKV","ZVP"), each=117),
                           halo = rep(c("AH","H","BH"), each=39))

abund_16S_dfF$cave<-factor(abund_16S_dfF$cave, c("VG","BP","GKV","ZVP"))
abund_16S_dfF$halo<-factor(abund_16S_dfF$halo, c("AH","H","BH"))

Halocline_rep <- factor(rep(c("AH","H","BH"), 4))

geomPOINT <- ggplot(abund_16S_dfF, aes(x=factor(sample, level=c("VG1",  "VG2",  "VG3","BP1",  "BP2",  "BP3", "GKV1", "GKV2", "GKV3",
                                                               "ZVP1", "ZVP2", "ZVP3")), X, size = value)) + 
  geom_point(alpha=0.8, aes(color=halo))+
  facet_grid(.~cave, scales = "free_x")+
  scale_y_discrete(limits = rev)+
  theme(axis.title = element_blank(),
        legend.position="bottom", legend.box = "horizontal",
        panel.background = element_blank(),
        panel.grid.major.x = element_line(size = 0.1, colour="gray"),
        axis.ticks = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.text.y =element_text(size=8),
        #axis.text.x = element_text(angle = 90, vjust=.5, hjust=.5),
        axis.text.x =element_blank(),
        panel.grid.minor.x = element_line(colour="black", size = (1.5)),
        axis.ticks.x = element_line(colour = "gray",
                                    size = 10),
        axis.ticks.y = element_line(colour = "gray",
                                    size = 4.5),
        strip.background.x = element_blank(),
        strip.text.x = element_blank())+
  labs(size="Relative abundance (%)",
       col="Sampling point")+
  scale_size_continuous(range = c(1, 8))+
  scale_color_manual(values=c("AH"="cadetblue3", "H"='gray', "BH"='blue'))
#scale_x_discrete(labels= Halocline_rep)

h2 <- ggplot(abund_16S_dfF)+
  geom_bar(mapping = aes(x = sample, y = 1, fill = cave), 
           stat = "identity", 
           width = 1)+
  theme_void()+
  theme(panel.spacing.x = unit(1, "mm"),
        plot.margin = unit(c(0.1,0.1,.1,.1), "cm"))+
  facet_grid(.~cave, scales = "free_x")+
  scale_fill_manual(values=c("VG"="#9AA490", "BP"='#AA643E', "GKV"='#B9A984',"ZVP"='#205858'))



legend <- plot_grid(get_legend(h2), get_legend(geomPOINT_F), ncol = 1)
h1 <- geomPOINT_F  + theme(legend.position="bottom", legend.box = "horizontal")+
  guides(size = guide_legend(order = 1, label.position="bottom", title.position = "top"),col = guide_legend(order = 2, label.position="bottom",
                                                                                                            title.position = "top",override.aes = list(size=8)), ncol=2, alpha=FALSE)
h2 <- h2 + theme(legend.position = "none")
plot <- plot_grid(h2, h1 , align = "v", ncol = 1, axis = "tb", rel_heights = c(0.5, 15))
plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 1.5))



zone_Alpha <- tibble(ymin = 25.5, ymax = 32.5, xmin = -Inf, xmax = Inf)
zone_Gamma <- tibble(ymin = 5.5, ymax = 14.5, xmin = -Inf, xmax = Inf)
zone_Delta <- tibble(ymin = 19.5, ymax = 20.5, xmin = -Inf, xmax = Inf)

geomPOINT_F <- geomPOINT +geom_rect(
  mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  data = zone_Alpha , 
  alpha = 0.2,
  fill = "gray",
  inherit.aes = FALSE
)+
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_Gamma, 
    alpha = 0.2,
    fill = "gray",
    inherit.aes = FALSE
  )+
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_Delta, 
    alpha = 0.2,
    fill = "gray",
    inherit.aes = FALSE
  )



ggsave(filename = "Relative_abundance_bubble_plot_16S.tiff", plot= plot  ,  width = 20, height = 20, units = "cm")










#Fig.7. Co-inertia analysis-----
otu_table_18St_n <- otu_table_18St %>% select("VG1",  "VG2",  "VG3","BP1",  "BP2",  "BP3", "GKV1", "GKV2", "GKV3",
                                                    "ZVP1", "ZVP2", "ZVP3")


otu_table_16St_n <- otu_table_16St %>% select("VG1",  "VG2",  "VG3","BP1",  "BP2",  "BP3", "GKV1", "GKV2", "GKV3",
                                                          "ZVP1", "ZVP2", "ZVP3")

rownames(otu_table_18St_n)==rownames(otu_table_16St_n)

table_16S.pco<-dudi.pco(vegdist(decostand(otu_table_16St_n,"hellinger"), "bray"),scannf = FALSE, nf = 2)
table_18S.pco<-dudi.pco(vegdist(decostand(otu_table_18St_n,"hellinger"), "bray"),scannf = FALSE, nf = 2)


coia.16S.18S <- coinertia(table_16S.pco, table_18S.pco,
                        scannf = FALSE,
                        nf = 2)
summary(coia.16S.18S)


RV.rtest(table_16S.pco$tab, table_18S.pco$tab, 999)

# Relative variation on first eigenvalue
coia.16S.18S$eig[2] / sum(coia.16S.18S$eig)


coia.16S.18S_table <- cbind(pro=data.frame(coia.16S.18S$mY),
                   euk=data.frame(coia.16S.18S$mX))

coia.16S.18S_table_rn <- coia.16S.18S_table %>% rename(x1 = pro.NorS1,
                                     y1 = pro.NorS2,
                                     x2 = euk.NorS1,
                                     y2 = euk.NorS2)


groups <- as.factor(rep(c("VG", "BP", "GKV", "ZVP"), c(3,3,3,3)))


text_CIA <- grobTree(textGrob("RV = 0.8369 p = 0.001", x=0.8,  y=0.98, hjust=0,
                                             gp=gpar(col="black", fontsize=13, fontface="bold")))

CIA_plot <- plot.arr(coia.16S.18S_table_rn)+
  theme_bw()+
  theme(plot.title = element_text(face = "bold", size = 12,hjust = 0.5),
        legend.position = "none")+
  xlab("Axis 1 (38.71% explained)") +
  ylab("Axis 2 (24.71% explained)") +
  annotation_custom(text_CIA )+
  scale_color_manual(values=c("VG"="#9AA490", "BP"='#AA643E', "GKV"='#B9A984',"ZVP"='#205858'))

ggsave(filename = "CIA_16S_18S_AC.tiff", plot=CIA_plot ,  width = 25, height = 20, units = "cm")








#plot.arr function---
#Source: https://stackoverflow.com/questions/39004617/connecting-large-points-with-lines-and-arrows-in-ggplot2
plot.arr <- function(df, pointsiz=2, pointsiz.scale.factor=100){
  
  ## calculate weights to adjust for aspect ratio
  norm2 <- function(v) sqrt(sum(v^2))
  w <- c(diff(range(df[ ,c("x1","x2")])), diff(range(df[ ,c("y1","y2")])))
  w <- w/norm2(w)
  
  ## use "elliptical" norm to account for different scales on x vs. y axes
  norm2w <- function(v) sqrt(sum((v/w)^2))
  
  ## compute normalized direction vectors, using "elliptical" norm
  direc <- do.call("rbind",lapply(1:nrow(df), function(i) {
    vec <- with(df[i, ], c(dx=x2-x1, dy=y2-y1))
    data.frame(as.list(vec/norm2w(vec)))
  }))
  
  ## "shift back" endpoints:
  ## translate endpoints towards startpoints by a fixed length;
  ## translation direction is given by the normalized vectors;
  ## translation length is proportional to the overall size of the plot
  ## along both x and y directions
  ## pointsiz.scale.factor can be decreased/increased for larger/smaller pointsizes
  
  epsil <- direc * diff(range(df)) / pointsiz.scale.factor
  df$xend2 <- df$x2 - epsil$dx
  df$yend2 <- df$y2 - epsil$dy
  
  g <- ggplot(df) + 
    geom_hline(yintercept=0,color='darkgray', alpha=0.8)+
    geom_vline(xintercept=0,color='darkgray', alpha=0.8)+
    #geom_point(aes(x=x1* -1, y=y1* -1)) + 
    #geom_point(aes(x=x2* -1, y=y2* -1)) +
    geom_segment(aes(x=x1* -1, y=y1* -1, xend=x2* -1, yend=y2* -1, colour=groups), 
                 arrow = arrow(length=unit(0.30,"cm")), size=1.4)+
    geom_point(aes(x=mid_X* -1, y=mid_Y* -1))+
    #geom_text(aes(x=mid_X* -1, y=mid_Y* -1),label=row.names(df),hjust=0, vjust=0) +
    coord_cartesian()+
    theme(legend.position = "none",
          plot.margin=unit(c(1,1,1,1), "cm"))+
    geom_label_repel(aes(label = row.names(df),x=mid_X* -1, y=mid_Y* -1),
                     box.padding   = 0.25, 
                     point.padding = 0.3,
                     segment.color = 'grey50',
                     direction="both",size=2,alpha=.7)
  
  print(g)
}
