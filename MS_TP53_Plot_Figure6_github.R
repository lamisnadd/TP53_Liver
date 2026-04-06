#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.6
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
FigDir<-"./Figure6/"


####
# rm(list=ls())
myinf1 = "./data/TCGA_LIHC_combined_info.txt"
myinf2 = "./data/TCGA_Firehose_RNASeqV2_ImmGen_representative_TILs_CLS.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
cnum = ncol(info)/4
info = info[, 1:cnum] - info[,(cnum+1):(2*cnum)]
colnames(info) = gsub("_up.ES","", colnames(info))
se = grep("LIHC__", row.names(info))
info = info[se,]
row.names(info) = gsub("LIHC__", "", row.names(info))

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]


#--------------------------------------
mytag = rep("", nrow(info))
mys1 = data$p53.score
mys2 = data$tis.score
mytag[mys1>=median(mys1, na.rm=T) & mys2>=median(mys2, na.rm=T)] = "HH"
mytag[mys1>=median(mys1, na.rm=T) & mys2<median(mys2, na.rm=T)] = "HL"
mytag[mys1<median(mys1, na.rm=T) & mys2>=median(mys2, na.rm=T)] = "LH"
mytag[mys1<median(mys1, na.rm=T) & mys2<median(mys2, na.rm=T)] = "LL"
# mytag[mys1<median(mys1, na.rm=T)] = "P53L"
table(mytag)
info$mytag = mytag

xx = apply(!is.na(info), 2, sum)
se = which(xx>100)
info = info[,se]

#--------------------------------------

par(mfrow=c(2,3))
for(k in 1:6)
{
  cat("\r", k)
  myf = info[,k]
  myList = list(NULL)
  myList[[1]] = myf[info$mytag=="HH"]
  myList[[2]] = myf[info$mytag=="HL"]
  myList[[3]] = myf[info$mytag=="LH"]
  myList[[4]] = myf[info$mytag=="LL"]
  boxplot(myList, outlier=F, main=colnames(info)[k])
}




P_score<-p_score<-P_tis<-p_tis<-MeasureName<-NULL
for (i in 1 :(length(info[1,])-1)){
  MeasureName[i]<-colnames(info)[i]
  P_tis[i]<-cor.test(info[,i],data$tis.score)$estimate
  P_score[i]<- cor.test(info[,i],data$p53.score)$estimate
  p_tis[i]<-cor.test(info[,i],data$tis.score)$p.value
  p_score[i]<- cor.test(info[,i],data$p53.score)$p.value
}

plot_cor <- data.frame(
  Measure = MeasureName,
  Correlation = c(P_score, P_tis),
  P_value = c(p_score, p_tis),
  Score_Type = rep(c("TP53", "TIS"), each = length(MeasureName))
)
plot_cor <- plot_cor %>%
  filter(!is.na(Correlation), !is.na(P_value))
library(ggplot2)

plot_cor$Measure <- factor(
  plot_cor$Measure,
  levels = plot_cor %>%
    filter(Score_Type == "TP53") %>%
    arrange(Correlation) %>%
    pull(Measure)
)

P<-ggplot(plot_cor,
          aes(y =Measure, 
              x = Score_Type,
              size = -log10(P_value),
              color = Correlation)) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(
    colors = c("blue4", "blue", "grey80", "red", "red4"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1)
  )+
  scale_size(range = c(2, 7)) +
  theme_classic() +
  labs(
    x = "Score",
    y = NULL,
    color = "Correlation",
    size = "-log10(p-value)",
    title = ""
  ) +
  
  theme(
    # axis.text.x = element_text(size = 12,angle = 45, hjust = 1),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16)  ,
    # ↓ Reduce legend size
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing = unit(0.2, "cm")
  )



graphics.off()
pdf(paste0(FigDir,"dot_immunecells.pdf"), width = 3.2, height = 2.5) 
print(P)
dev.off()

####################################
myinf1 = "./data/TCGA_LIHC_Freq_GenomicEvents.txt"
myinf2 = "./data/Thorsson_2018_TCGA_immunelandscape.csv"

data = read.table(myinf1, sep="\t", header=T, row.names=1)
row.names(data) = substr(row.names(data), 1, 12)
#------------------------
xx = apply(is.na(data), 1, sum)
se = which(xx==0)
mydat = data[se,]
for(k in 1:ncol(mydat))
{
  mydat[,k] = ifelse(mydat[,k]>0, 1, 0)
}
data = mydat

#------------------------
 info = read.table(myinf2, sep=",", header=T, row.names=1, quote="")
se = which(info$TCGA.Study=="LIHC")
info = info[se,]
subset = c(4:31,36:63)
subset = c(4:31)
info = info[, subset]
# info$
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

info$Group = ifelse(data$TP53__MUT==1, "Mut", "Wt")
raw.data = info

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
myFig = "/Fig0_Boxplot_TIL_TP53Mut_vs_WT_TCGALIHC.pdf"

data = raw.data
mycat = data$Group
mys = data$Lymphocyte.Infiltration.Signature.Score
data = data.frame(mys, mycat)

xx1 = data$mys[data$mycat=="Mut"]
xx2 = data$mys[data$mycat=="Wt"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))

myp = wilcox.test(xx1, xx2)$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = ifelse(myp<0.05, paste("P=", myp, sep=""), "P > 0.1")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

data$mycat <- factor(data$mycat, levels = c("Mut", "Wt") )

mycol = c('darkgreen', 'coral1')[1:2]
mycol =c("red3", "grey50")
mygg1 <- ggplot(data, aes(x=mycat, y=mys,  color=mycat)) + 
  geom_boxplot(size=1, outlier.shape = NA) +
  scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  labs(title='',x= 'TP53 mutation', y = 'TIL level')+
  # scale_y_continuous(limits = c(0, 0.6)) +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )+  annotate("text", x = 1.5, y = 2.8, label = txt.myp, color = "black", size = 4)

dev.off()
pdf(file = paste0(FigDir,myFig), width = 2.2, height = 3) 
mygg1
dev.off()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
myFig = "/Fig0_Boxplot_Leukocyte_TP53Mut_vs_WT_TCGALIHC.pdf"

data = raw.data
mycat = data$Group
mys = data$Leukocyte.Fraction
data = data.frame(mys, mycat)

xx1 = data$mys[data$mycat=="Mut"]
xx2 = data$mys[data$mycat=="Wt"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))

myp = wilcox.test(xx1, xx2)$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = ifelse(myp<0.05, paste("P=", myp, sep=""), "P > 0.1")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)
data$mycat <- factor(data$mycat, levels = c("Mut", "Wt") )

mycol =c("sienna2", "palegreen4")
mycol =c("red3", "grey50")
mygg2 <- ggplot(data, aes(x=mycat, y=mys,  color=mycat)) + 
  geom_boxplot(size=1, outlier.shape = NA) +
  scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  labs(title='',x= 'TP53 mutation', y = 'Leukocyte fraction')+
  scale_y_continuous(limits = c(0, 0.6)) +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )+  annotate("text", x = 1.5, y = 0.55, label = txt.myp, color = "black", size = 4)

mygg2

dev.off()
pdf(file = paste0(FigDir,myFig), width = 2.2, height = 3) 
mygg2
dev.off()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
myFig = "/Fig0_Boxplot_IFNg_TP53Mut_vs_WT_TCGALIHC.pdf"

data = raw.data
mycat = data$Group
mys = data$IFN.gamma.Response
data = data.frame(mys, mycat)

xx1 = data$mys[data$mycat=="Mut"]
xx2 = data$mys[data$mycat=="Wt"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))

myp = wilcox.test(xx1, xx2)$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = ifelse(myp<0.05, paste("P=", myp, sep=""), "P > 0.1")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)


data$mycat <- factor(data$mycat, levels = c("Mut", "Wt") )

mygg <- ggplot(data, aes(x=mycat, y=mys,  fill=mycat)) + 
  geom_boxplot(size=1, outlier.shape = NA) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= 'TP53 mutation', y = 'IFN-gamma response')+
  # scale_y_continuous(limits = c(0, 0.6)) +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )

mycol =c("sienna2", "palegreen4")
mycol =c("red3", "grey50")
mygg3 <- ggplot(data, aes(x=mycat, y=mys,  color=mycat)) + 
  geom_boxplot(size=1, outlier.shape = NA) +
  scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  labs(title='',x= 'TP53 mutation', y = 'IFN-gamma response')+
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )+  annotate("text", x = 1.5, y = 1.3, label = txt.myp, color = "black", size = 4)

mygg3

dev.off()
pdf(file =paste0(FigDir,myFig), width = 2.2, height = 3) 
mygg3
dev.off()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##genomic and imunological features 
# rm(list=ls()


#########
myinf1 = "./data/TCGA_LIHC_combined_info.txt"
myinf2 = "./data/Thorsson_2018_TCGA_immunelandscape.csv"


data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info = read.table(myinf2, sep=",", header=T, row.names=1, quote="")
se = which(info$TCGA.Study=="LIHC")
info = info[se,]


comxx = intersect(row.names(data), row.names(info))


data = data[comxx,]
info = info[comxx,]

subset = c(4:31, 36:63)
subset = c(4:28)
info = info[, subset]

mytag = rep("", nrow(info))
mys1 = data$p53.score
mys2 = data$tis.score
mytag[mys1>=median(mys1, na.rm=T) & mys2>=median(mys2, na.rm=T)] = "HH"
mytag[mys1>=median(mys1, na.rm=T) & mys2<median(mys2, na.rm=T)] = "HL"
mytag[mys1<median(mys1, na.rm=T) & mys2>=median(mys2, na.rm=T)] = "LH"
mytag[mys1<median(mys1, na.rm=T) & mys2<median(mys2, na.rm=T)] = "LL"

# mytag[mys1<median(mys1, na.rm=T)] = "P53L"
table(mytag)

info$mytag = mytag

p53_group = rep("", nrow(info))
mys1 = data$p53.score
p53_group[mys1>=median(mys1, na.rm=T) ] = "High"
p53_group[mys1<median(mys1, na.rm=T) ] = "Low"
table(p53_group)
info$p53_group = p53_group


xx = apply(!is.na(info), 2, sum)
se = which(xx>100)
info = info[,se]


xx = apply(!is.na(info), 2, sum)
se = which(xx>100)
info = info[,se]

#--------------------------------------
par(mfrow=c(4,6))
for(k in 1:22)
{
  cat("\r", k)
  myf = info[,k]
  myList = list(NULL)
  myList[[1]] = myf[info$mytag=="HH"]
  myList[[2]] = myf[info$mytag=="HL"]
  myList[[3]] = myf[info$mytag=="LH"]
  myList[[4]] = myf[info$mytag=="LL"]
  boxplot(myList, outlier=F, main=colnames(info)[k])
}


P_score<-p_score<-P_tis<-p_tis<-MeasureName<-NULL
for (i in 1 :(length(info[1,])-2)){
  # for (i in 1 :22){
  MeasureName[i]<-colnames(info)[i]
  P_tis[i]<-cor.test(info[,i],data$tis.score)$estimate
  P_score[i]<- cor.test(info[,i],data$p53.score)$estimate
  p_tis[i]<-cor.test(info[,i],data$tis.score)$p.value
  p_score[i]<- cor.test(info[,i],data$p53.score)$p.value
}

plot_cor <- data.frame(
  Measure = MeasureName,
  Correlation = c(P_score, P_tis),
  P_value = c(p_score, p_tis),
  Score_Type = rep(c("TP53", "TIS"), each = length(MeasureName))
)
plot_cor <- plot_cor %>%
  filter(!is.na(Correlation), !is.na(P_value))
library(ggplot2)

plot_cor$Measure <- factor(
  plot_cor$Measure,
  levels = plot_cor %>%
    filter(Score_Type == "TP53") %>%
    arrange(Correlation) %>%
    pull(Measure)
)

plot_cor$Measure <- factor(
  plot_cor$Measure,
  levels = plot_cor %>%
    filter(Score_Type == "TIS") %>%
    arrange(Correlation) %>%
    pull(Measure)
)

P<-ggplot(plot_cor,
       aes(y = Measure,
           x = Score_Type,
           size = -log10(P_value),
           color = Correlation)) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(
    colors = c("blue4", "blue", "white", "red", "red4"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1)
  )+
  scale_size(range = c(2, 7)) +
  theme_classic() +
  labs(
    x = "Score",
    y = NULL,
    color = "Correlation",
    size = "-log10(p-value)",
    title = ""
  ) +
  # theme(
  #   axis.text.x = element_text(angle = 45, hjust = 1)
  # )+
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 13),
    plot.title   = element_text(size = 17, hjust = 0.5)  # center title
  )

graphics.off()
pdf(paste0(FigDir,"dot_immunemeasures4.pdf"), width = 5.8, height = 13) 
print(P)
dev.off()

graphics.off()
pdf(paste0(FigDir,"dot_immunemeasures3.pdf"), width = 5.8, height = 7) 
print(P)
dev.off()


#+++++++++++++++++++++++++++++++++++




#--------------------------------------
######################
myinf1 = "./data/TCGA_LIHC_combined_info.txt"
myinf2 = "./data/Thorsson_2018_TCGA_immunelandscape.csv"

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info = read.table(myinf2, sep=",", header=T, row.names=1, quote="")
se = which(info$TCGA.Study=="LIHC")
info = info[se,]


comxx = intersect(row.names(data), row.names(info))


data = data[comxx,]
info = info[comxx,]

subset = c(4:28)
info = info[, subset]


# Categorize p53 and TIS into low/high based on their medians
median_p53 <- median(mys1, na.rm = TRUE)
median_tis <- median(mys2, na.rm = TRUE)
plot(mys1,mys2)


library(ggpubr)

df <- data.frame(P53 = mys1,
                 TIS = mys2)

ggscatter(df,
          x = "P53",
          y = "TIS",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "P53",
          ylab = "TIS")


# Categorize p53 (low/high) and TIS (low/high)
p53_group <- ifelse(mys1 >= median_p53, "High", "Low")
tis_group <- ifelse(mys2 >= median_tis, "High", "Low")

# Create a data frame for plotting
plot_data <- data.frame(p53_group, tis_group)


# Create a contingency table for p53 and TIS groups
contingency_table <- table(p53_group, tis_group)

# Perform the Chi-squared test
chisq_test <- chisq.test(contingency_table)

# Check the p-value
chisq_test$p.value
# Remove samples with NA values in either p53 or TIS
plot_data <- data.frame(p53_group, tis_group)
plot_data <- plot_data[!is.na(plot_data$p53_group) & !is.na(plot_data$tis_group), ]

# Create the stacked bar plot
P<-ggplot(plot_data, aes(x = p53_group, fill = tis_group)) +
  geom_bar(position = "stack") + 
  labs(title = paste("p-value = ", 
                     round(chisq_test$p.value, 4), sep = ""), 
       x = "TP53 score", 
       y = "Number",
       fill = "TIS group" ) +
  scale_fill_manual(values = c("salmon4", "skyblue3")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(size = 16, hjust = 0.5)  # center title
  )


graphics.off()
pdf(paste0(FigDir,"bar_tis_score.pdf"), width = 3, height = 2.2) 
print(P)
dev.off()

mytag = rep("", nrow(info))
mys1 = data$p53.score
mys2 = data$tis.score
mytag[mys1>=median(mys1, na.rm=T) & mys2>=median(mys2, na.rm=T)] = "HH"
mytag[mys1>=median(mys1, na.rm=T) & mys2<median(mys2, na.rm=T)] = "HL"
mytag[mys1<median(mys1, na.rm=T) & mys2>=median(mys2, na.rm=T)] = "LH"
mytag[mys1<median(mys1, na.rm=T) & mys2<median(mys2, na.rm=T)] = "LL"
# mytag[mys1<median(mys1, na.rm=T)] = "P53L"
table(mytag)

info$mytag = mytag

p53_group = rep("", nrow(info))
mys1 = data$p53.score
p53_group[mys1>=median(mys1, na.rm=T) ] = "High"
p53_group[mys1<median(mys1, na.rm=T) ] = "Low"
table(p53_group)
info$p53_group = p53_group


xx = apply(!is.na(info), 2, sum)
se = which(xx>100)
info = info[,se]


xx = apply(!is.na(info), 2, sum)
se = which(xx>100)
info = info[,se]


####
p53_features <- c(
  "Proliferation", ###
  "Aneuploidy.Score", ###
  "Fraction.Altered",
  "Number.of.Segments", ###
  "Nonsilent.Mutation.Rate",
  "Silent.Mutation.Rate",
  "SNV.Neoantigens",
  "Indel.Neoantigens",
  "Homologous.Recombination.Defects", ###
  "Intratumor.Heterogeneity", ###
  "CTA.Score", ###
  "Wound.Healing" ###
)

p53_features <- c(
  "Proliferation", ###
  "Aneuploidy.Score", ###
  "Number.of.Segments", ###
  "Homologous.Recombination.Defects", ###
  "Intratumor.Heterogeneity", ###
  "CTA.Score", ###
  "Wound.Healing" ###
)

t <- c(
  "Proliferation", ###
  "Aneuploidy", ###
  "Number of Segments", ###
  "HRD",
  "Intratumor Heterogeneity", ###
  "CTA Score", ###
  "Wound Healing" )


## score
p53_features <- c(
"Proliferation",
"Homologous.Recombination.Defects",
"Wound.Healing")

t <- c(
  "Proliferation",
  "HRD",
  "Wound Healing")

my_comparisons <- list(
  c("HH", "HL"),
  c("LH", "LL"),
  c("HH", "LH"),
  c("HL", "LL"))




library(ggplot2)
library(dplyr)
library(tidyr)

f="Aneuploidy.Score"
i=2
for (f in p53_features) {
  i=i+1

  df_long <- info %>%
    dplyr::filter(mytag %in% c("HH", "HL", "LH", "LL")) %>%
    dplyr::select(mytag,  all_of(f)) %>%
    dplyr::rename(value = all_of(f))
    
  p <- ggplot(df_long, aes(x = mytag, y = value, color = mytag)) +
    # geom_boxplot(size=1,outlier.size = 0.3,linewidth=0.5) +
    geom_boxplot(size=1, outlier.shape = NA) +
    # geom_jitter(width = 0.15, size = 0.3, alpha = 0.7) +
    scale_fill_manual(values="white") +
    scale_color_manual(values = c(
      "LL" = "#0D3B66",    # deeper, slightly muted blue for low/low
      "LH" = "#2E8B90",    # strong forest green for low/high
      "HL" = "palegreen3",    # bright yellow-green for high/low
      "HH" = "lightgoldenrod2"    # rich golden yellow for high/high
    ))+
    labs(
      title = t[i],
      y = t[i],
      x = "Patient group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12),       
      axis.text.y = element_text(size = 12),       
      axis.title.x = element_text(size = 14),      
      axis.title.y = element_text(size = 14),   
      legend.position = "none",
      plot.title = element_text(hjust = 0.5))+
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",      
      label = "p.signif",
      size = 3.5   # Increase size of p-values
    ) +
   
    coord_cartesian(
      ylim = c(
        min(df_long$value, na.rm = TRUE),
        max(df_long$value, na.rm = TRUE) * 1.5
      )
    )
  
  assign(paste0("fig4G",f),p)
  
  graphics.off()
  pdf(paste0(FigDir,"fig4G",f,".pdf"), width = 2.4, height = 3) 
  print(p)
  dev.off()
}





immune_features <- c(
  # Immune activation / inflammation
  "Lymphocyte.Infiltration.Signature.Score", ##
  "IFN.gamma.Response",##
  "TGF.beta.Response",##
  
  # T cell repertoire
  "TCR.Shannon", ##
  "TCR.Richness", ##
  "TCR.Evenness",  ## higher in HL
  
  # Immune cell fractions
  "T.Cells.CD8", ##
  # "T.Cells.CD4.Memory.Activated",
  "T.Cells.Regulatory.Tregs",##
  # "NK.Cells.Activated",
  "Macrophages.M1", ##
  "Macrophages.M2", ## lower in HH
  
  # Overall immune content
  "Leukocyte.Fraction", ## 
  "Lymphocytes"##
)

key_immune <- c(
  "T.Cells.CD8",
  "IFN.gamma.Response",
  "TCR.Shannon",
  "Lymphocyte.Infiltration.Signature.Score",
  "NK.Cells.Activated",
  "T.Cells.Regulatory.Tregs"
)
key_immune <- c(
  "T.Cells.CD8",
  "IFN.gamma.Response",
  "TCR.Shannon",
  "Lymphocyte.Infiltration.Signature.Score",
  "NK.Cells.Activated",
  "T.Cells.Regulatory.Tregs"
)
##tis
key_immune <- c(
  # "TCR.Richness",
# "Leukocyte.Fraction",
"Lymphocyte.Infiltration.Signature.Score",
"Macrophage.Regulation",
# "TCR.Shannon",
"IFN.gamma.Response"
)



t <- c(
  # "TCR Richness",
       # "Leukocyte Fraction",
       "TIL score",
       "Macrophage Regulation",
       # "TCR Shannon",
       "IFN-gamma response")


my_comparisons <- list(
  c("HH", "HL"),
  c("LH", "LL"),
  c("HH", "LH"),
  c("HL", "LL"))




library(ggplot2)
library(dplyr)
library(tidyr)

f="Macrophage.Regulation"
i=2
for (f in key_immune) {
  i=i+1
  df_long <- info %>%
    dplyr::filter(mytag %in% c("HH", "HL", "LH","LL")) %>%
    dplyr::select(mytag, all_of(f)) %>%
    dplyr::rename(value = all_of(f))
  
  
  p <- ggplot(df_long, aes(x = mytag, y = value, color = mytag)) +
    geom_boxplot(size=1, outlier.shape = NA) +
    # geom_jitter(width = 0.15, size = 0.3, alpha = 0.7) +
    scale_fill_manual(values="white") +
    scale_color_manual(values = c(
      "LL" = "#0D3B66",    # deeper, slightly muted blue for low/low
      "LH" = "#2E8B90",    # strong forest green for low/high
      "HL" = "palegreen3",    # bright yellow-green for high/low
      "HH" = "lightgoldenrod2"    # rich golden yellow for high/high
    ))+
    labs(
      title = t[i],
      y = t[i],
      x ="Patient group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12),       
      axis.text.y = element_text(size = 12),       
      axis.title.x = element_text(size = 14),      
      axis.title.y = element_text(size = 14),   
      legend.position = "none",
      plot.title = element_text(hjust = 0.5))+
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",      
      label = "p.signif",
      size = 3.5   # Increase size of p-values
    )  +
    coord_cartesian(
      ylim = c(
        min(df_long$value, na.rm = TRUE),
        max(df_long$value, na.rm = TRUE) * 2.2
      )
    )
  
  assign(paste0("fig4H",f),p)
  
  graphics.off()
  pdf(paste0(FigDir,"fig4H_",f,".pdf"), width = 2.4, height = 3) 
  print(p)
  dev.off()
}






library(ggplot2)
library(dplyr)
library(ggpubr)


info$p53_group <- factor(info$p53_group, levels = c("Low", "High"))
my_comparisons <- list(
  c("Low", "High"))
for (f in key_immune ) {
  
  df_long <- info %>%
    dplyr::filter(!is.na(p53_group)) %>%
    dplyr::select(p53_group, all_of(f)) %>%
    dplyr::rename(value = all_of(f))
  
  max_y <- max(df_long$value, na.rm = TRUE)
  min_y <- min(df_long$value, na.rm = TRUE)
  
  p <- ggplot(df_long, aes(x = p53_group, y = value, color = p53_group)) +
    geom_boxplot(size = 0.6, outlier.size = 0.3) +
    scale_color_manual(
      values = c(
        "Low"  = "grey40",
        "High" = "lightsalmon2"
      )
    ) +
    labs(
      title = f,
      y = f,
      x = "TP53 score"
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",      
      label = "p.signif",
      size = 3.5   # Increase size of p-values
    )  +
    coord_cartesian(ylim = c(min_y, max_y * 1.25))
  
  pdf(paste0(FigDir, "fig4I_", f, ".pdf"), width = 2.5, height = 3)
  print(p)
  dev.off()
}


