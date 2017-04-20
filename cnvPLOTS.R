####### GRAPH GENERATION #######
rm(list = ls())
###load data from cnvANNO.R script
load(file="cnvANNO.RData")

###load required libraries
require(methods)

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}
if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}
###histgram of number of exons per CNV
y <- as.data.frame(table(x$CNV_ID)) 
bin <- 1
scale.y <- bin * (length(y$Freq)) ##for density plot scaling
q <- quantile(y$Freq)

###histogram of exons per cnv - all types
p1 <- ggplot(y, aes(x = Freq)) +
        geom_histogram(aes(y = ..count..), binwidth = bin, fill = "steelblue3") +
        geom_vline(aes(xintercept= q[3], linetype = "Median"), size = 1, alpha = 0.9, colour = "darkorange2", show.legend = TRUE) +
        scale_linetype_manual(name = "", values = c(2, 2), guide = guide_legend(override.aes = list(color = "darkorange2"))) +
        labs(list(title = "Number of Exons per CNV Call (Frequency Plot)", y = "Frequency/Count", x = "Exons per CNV", fill = "CNV Type")) +
        theme(panel.border = element_blank(), axis.line = element_line(colour="black"), panel.grid.major = element_line(colour = "gray90")) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 14,colour = "gray10",margin = margin(0,0,10,0))) +
        theme(legend.key = element_blank(), axis.title = element_text(size = 13)) +
        scale_x_continuous(expand = c(0,0), limits = c(0.1,round(max(y$Freq),-2)), breaks = seq(0,round(max(y$Freq),-2), by = 50)) +
        scale_y_continuous(expand = c(0,0), limits = c(NA,150), breaks = seq(0,150, by = 30))
      

###DEL/DUP Count Comparisons
p2 <- ggplot(x, aes(x = x$CNV)) +
        geom_bar(fill = c("steelblue3", "goldenrod2"), colour = c("gray12","gray12")) +
        labs(list(title = "Total Number Targets - CNV Calls", x = "CNV Type", y = "", fill = "CNV Type")) +
        theme(panel.border = element_blank(), axis.line = element_line(colour="black"), panel.grid.major = element_line(colour = "gray90")) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 14,colour = "gray10",margin = margin(0,0,10,0))) +
        theme(axis.title = element_text(size = 13), plot.margin = margin(10,25,0,0), axis.text.y = element_text(size=8)) +
        scale_y_continuous(breaks = pretty(1:max(c(length(x$CNV[x$CNV == "DUP"]), length(x$CNV[x$CNV == "DEL"]))), n = 10)) +
        coord_flip()

p3 <- ggplot(x,aes(x = SAMPLE, fill = CNV)) +
        geom_bar(colour = "gray12", width = 0.7) +
        #geom_text(aes(x= SAMPLE, y=value, label=value)) +
        labs(list(title = "CNV Calls Per Sample", y = "Count", x = "Sample", fill = "CNV Type")) +
        theme(panel.border = element_blank(), axis.line = element_line(colour="black"), panel.grid.major = element_line(colour = "gray90")) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 14,colour = "gray10",margin = margin(0,0,10,0))) +
        theme(axis.title = element_text(size = 13), plot.margin = margin(10,25,0,0), axis.text.x = element_blank()) +
        theme(panel.grid.major.x = element_blank(), axis.line.x = element_blank(), legend.title = element_text(size = 10)) +
        theme(axis.ticks.x = element_blank(), axis.title.x = element_text(margin = margin(b = 10))) +
        scale_fill_manual(values = c("steelblue3", "goldenrod2"))+
        scale_y_continuous(expand = c(0,0), limits = c(NA,800), breaks = pretty(1:800, n = 11))

z.DUP <- as.data.frame(table(x$CNV_ID[x$CNV == "DUP"]))
z.DEL <- as.data.frame(table(x$CNV_ID[x$CNV == "DEL"]))

p4 <- ggplot() +
        geom_density(aes(x = z.DUP$Freq, fill = "DUP"), alpha = 0.5, size=0.1) +
        geom_density(aes(x = z.DEL$Freq, fill = "DEL"), alpha = 0.5, size=0.1) +
        labs(list(title = "Density Plot - Deletions & Duplications", x = "Exons per CNV", y = "Density", fill = "CNV Type")) +
        scale_fill_manual(values = c("steelblue3", "goldenrod2")) +
        theme(panel.border = element_blank(), axis.line = element_line(colour="black"), panel.grid.major = element_line(colour = "gray90")) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 14,colour = "gray10",margin = margin(0,0,10,0))) +
        theme(legend.key = element_blank(), axis.title = element_text(size = 13),legend.title = element_text(size = 10)) +
        scale_x_continuous(expand = c(0.005,0), limits = c(NA,500)) +
        scale_y_continuous(expand = c(0,0), limits = c(NA,0.07))


p5 <- ggplot() +
        geom_density(aes(x = z.DUP$Freq, fill = "DUP"), alpha = 0.5, size=0.1) +
        #geom_density(aes(x = z.DEL$Freq, fill = "DEL"), alpha = 0.5, size=0.1) +
        labs(list(title = "Density Plot - Deletions & Duplications", x = "Exons per CNV", y = "Density", fill = "CNV Type")) +
        scale_fill_manual(values = c("goldenrod2")) +
        theme(panel.border = element_blank(), axis.line = element_line(colour="black"), panel.grid.major = element_line(colour = "gray90")) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 14,colour = "gray10",margin = margin(0,0,10,0))) +
        theme(legend.key = element_blank(), axis.title = element_text(size = 13),legend.title = element_text(size = 10)) +
        scale_x_continuous(expand = c(0.005,0), limits = c(NA,500)) +
        scale_y_continuous(expand = c(0,0), limits = c(NA,0.07))

p6 <- ggplot() +
        #geom_density(aes(x = z.DUP$Freq, fill = "DUP"), alpha = 0.5, size=0.1) +
        geom_density(aes(x = z.DEL$Freq, fill = "DEL"), alpha = 0.5, size=0.1) +
        labs(list(title = "Density Plot - Deletions & Duplications", x = "Exons per CNV", y = "Density", fill = "CNV Type")) +
        scale_fill_manual(values = c("steelblue3")) +
        theme(panel.border = element_blank(), axis.line = element_line(colour="black"), panel.grid.major = element_line(colour = "gray90")) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 14,colour = "gray10",margin = margin(0,0,10,0))) +
        theme(legend.key = element_blank(), axis.title = element_text(size = 13),legend.title = element_text(size = 10)) +
        scale_x_continuous(expand = c(0.005,0), limits = c(NA,500)) +
        scale_y_continuous(expand = c(0,0), limits = c(NA,0.07))

mid <- plot_grid(p5,p6,p4, nrow=1, labels = c("B","C","D"))
bot <- plot_grid(p1,p2, nrow=1, labels = c("E","F"))
plot1 <- plot_grid(p3,mid,bot, ncol=1, labels = c("A","",""))

png(file = "exons_per_cnv_freqplot.png", width = 16, height = 9, units = "in", res = 300)
print(p1)
dev.off()

png(file = "total_cnv_calls_bytype.png", width = 16, height = 9, units = "in", res = 300)
print(p2)
dev.off()

png(file = "cnv_calls_persample_bytype.png", width = 16, height = 9, units = "in", res = 300)
print(p3)
dev.off()

png(file = "exons_per_cnv_densityplot.png", width = 16, height = 9, units = "in", res = 300)
print(p4)
dev.off()

png(file = "combined_plot.png", width = 20, height = 12, units = "in", res = 300)
print(plot1)
dev.off()
