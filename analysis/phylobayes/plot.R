library(ggplot2)
library(cowplot)

black = "black"
blue = "#439FC3"
green = "#76A88D"
orange = "#F9842A"
red = "#B52A24"
pink = "#CC79A7"
purple = "#8E7790"

setwd("/Users/sarah/Google Drive/GS/lab/divergence_timing_manuscript/analysis/phylobayes")
doud = read.csv("outputs/model_by_model_norm.csv")
m <- ceiling(max(max(doud$distance1), max(doud$distance2)))
doud_branch_lengths = ggplot(doud, aes(distance1, distance2)) + geom_abline(intercept=0, slope=1, color = "grey", linetype="dashed") + geom_point(size=1.0, color="grey") 
doud_branch_lengths = doud_branch_lengths + facet_grid(divergence_level~ model_pair)
doud_branch_lengths = doud_branch_lengths + theme(strip.background = element_blank(), strip.text = element_text(size=8), legend.position="bottom", panel.spacing = unit(2, "lines"), panel.border = element_rect(colour = "black"))
doud_branch_lengths = doud_branch_lengths + xlim(c(0,m)) + ylim(c(0,m)) + coord_fixed() 
doud_branch_lengths = doud_branch_lengths + xlab("model 1") + ylab("model 2") 
doud_branch_lengths = doud_branch_lengths + geom_point(size = 0.5, data = subset(doud, from_ref == 'WSN'), color=blue) + geom_point(size = 0.5, data = subset(doud, from_ref == 'Perth'), color=orange)
doud_branch_lengths
ggsave("outputs/branch_lengths.png", width =12)
