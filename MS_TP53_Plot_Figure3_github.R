#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.3
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

figDir<-"./Figure3/"
myinf1 = "./data/Hoshida_GSE15654__GenomicEvent_iRAS.txt"
myinf2 = "./data/Hoshida_GSE15654/Clinical_info.txt"
min(info$prediction.confidence.p.value[which(info$prediction=="Poor prognosis")])
data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(2*cnum)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp
se = grep("mul.adj__", colnames(data))
data = data[, se]
for(i in 1:ncol(data))
{
  data[,i] = data[,i]/sd(abs(data[,i]))
}
colnames(data) = gsub("mul.adj__", "", colnames(data))

##
info = read.table(myinf2, sep="\t", header=T, row.names=1)

##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(data)		## n=216

info$score = data$TP53__MUT


idx = rank(info$score)
mycat = rep("Int", nrow(info))
mycat[idx>161] = "Poor"
mycat[idx<=60] = "Good"
info$P53_Score_Prediction = mycat
table(info$P53_Score_Prediction)




library(survival)
# PV1 = HR1 = HR.lo = HR.hi = rep(0,6)

info$logp <- log(info$score + 1e-6)  # avoid log(0)

info_no_varices<-info[which(info$presence.of.varices==0),]
info_varices<-info[which(info$presence.of.varices==1),]

# Named list of Cox models
cox_models <- list(
  cox_score = coxph(Surv(days.to.hcc, hcc) ~ score, data = info),
  cox_score_OS = coxph(Surv(days.to.death, death) ~ score, data = info),
  
  cox_score_cat = coxph(Surv(days.to.hcc, hcc) ~ P53_Score_Prediction, data = info),
  cox_score_cat_OS = coxph(Surv(days.to.death, death) ~ P53_Score_Prediction, data = info),
  
  # Prediction only
  cox_pred = coxph(Surv(days.to.hcc, hcc) ~ prediction..p.0.05., data = info),
  cox_pred_OS = coxph(Surv(days.to.death, death) ~ prediction..p.0.05. ,data = info),
  
  # Combined prediction + score
  cox_combined = coxph(Surv(days.to.hcc, hcc) ~  P53_Score_Prediction+prediction..p.0.05. , data = info),
  cox_combined_OS = coxph(Surv(days.to.death, death) ~ prediction..p.0.05. + P53_Score_Prediction, data = info),
  
  # Bilirubin
  cox_bilirubin = coxph(Surv(days.to.hcc, hcc) ~ bilirubin...1.0mg.dl, data = info),
  cox_bilirubin_combined = coxph(Surv(days.to.hcc, hcc) ~ bilirubin...1.0mg.dl + score, data = info),
  
  # Varices
  cox_varices = coxph(Surv(days.to.hcc, hcc) ~ presence.of.varices, data = info),
  cox_varices_combined = coxph(Surv(days.to.hcc, hcc) ~ presence.of.varices +  score, data = info),
  
  # Platelet
  cox_platelet = coxph(Surv(days.to.hcc, hcc) ~ platelet...100.000.mm3, data = info),
  cox_platelet_combined = coxph(Surv(days.to.hcc, hcc) ~ platelet...100.000.mm3 + score, data = info),
  
  # Bilirubin
  cox_bilirubin_OS = coxph(Surv(days.to.death, death) ~ bilirubin...1.0mg.dl, data = info),
  cox_bilirubin_combined_OS = coxph(Surv(days.to.death, death) ~ bilirubin...1.0mg.dl + score, data = info),
  
  # Varices
  cox_varices_OS = coxph(Surv(days.to.death, death) ~ presence.of.varices, data = info),
  cox_varices_combined_OS = coxph(Surv(days.to.death, death) ~ presence.of.varices + score, data = info),
  
  # Platelet
  cox_platelet_OS = coxph(Surv(days.to.death, death) ~ platelet...100.000.mm3, data = info),
  cox_platelet_combined_OS = coxph(Surv(days.to.death, death) ~ platelet...100.000.mm3 + score, data = info)

)


cox_subgroup_models <- list(
  
  # All samples
  all = coxph(Surv(days.to.hcc, hcc) ~ score, data = info),
  
  # Bilirubin subgroups
  bilirubin_high = coxph(
    Surv(days.to.hcc, hcc) ~ score,
    data = subset(info, bilirubin...1.0mg.dl == 1)
  ),
  bilirubin_low = coxph(
    Surv(days.to.hcc, hcc) ~ score,
    data = subset(info, bilirubin...1.0mg.dl == 0)
  ),
  
  # Platelet subgroups
  platelet_high = coxph(
    Surv(days.to.hcc, hcc) ~ score,
    data = subset(info, platelet...100.000.mm3 == 0)
  ),
  platelet_low = coxph(
    Surv(days.to.hcc, hcc) ~ score,
    data = subset(info, platelet...100.000.mm3 == 1)
  ),
  
  # Varices subgroups
  varices_yes = coxph(
    Surv(days.to.hcc, hcc) ~ score,
    data = subset(info, presence.of.varices == 1)
  ),
  varices_no = coxph(
    Surv(days.to.hcc, hcc) ~ score,
    data = subset(info, presence.of.varices == 0)
  )
)

cox_subgroup_models_OS <- list(
  
  # All samples
  all = coxph(Surv(days.to.death, death) ~ score, data = info),
  
  # Bilirubin subgroups
  bilirubin_high = coxph(
    Surv(days.to.death, death) ~ score,
    data = subset(info, bilirubin...1.0mg.dl == 1)
  ),
  bilirubin_low = coxph(
    Surv(days.to.death, death) ~ score,
    data = subset(info, bilirubin...1.0mg.dl == 0)
  ),
  
  # Platelet subgroups
  platelet_high = coxph(
    Surv(days.to.death, death) ~ score,
    data = subset(info, platelet...100.000.mm3 == 0)
  ),
  platelet_low = coxph(
    Surv(days.to.death, death) ~ score,
    data = subset(info, platelet...100.000.mm3 == 1)
  ),
  
  # Varices subgroups
  varices_yes = coxph(
    Surv(days.to.death, death) ~ score,
    data = subset(info, presence.of.varices == 1)
  ),
  varices_no = coxph(
    Surv(days.to.death, death) ~ score,
    data = subset(info, presence.of.varices == 0)
  )
)

# Helper function for one Cox model
extract_cox_stats <- function(cox_model) {
  s <- summary(cox_model)
  data.frame(
    Variable = rownames(s$coefficients),
    HR = s$conf.int[,"exp(coef)"],
    HR.lo = s$conf.int[,"lower .95"],
    HR.hi = s$conf.int[,"upper .95"],
    PV1 = s$coefficients[,"Pr(>|z|)"],
    N = cox_model$n,           # sample size
    row.names = NULL
  )
}

# Apply to all models and combine
raw.data_subgroups <- do.call(rbind, lapply(names(cox_subgroup_models), function(nm) {
  df <- extract_cox_stats(cox_subgroup_models[[nm]])
  df$Model <- nm
  df
}))

raw.data_subgroups_OS <- do.call(rbind, lapply(names(cox_subgroup_models_OS), function(nm) {
  df <- extract_cox_stats(cox_subgroup_models_OS[[nm]])
  df$Model <- nm
  df
}))


# Apply to all models and combine
raw.data <- do.call(rbind, lapply(names(cox_models), function(nm) {
  df <- extract_cox_stats(cox_models[[nm]])
  df$Model <- nm
  df
}))

# Reorder columns nicely
raw.data <- raw.data[, c("Model", "Variable", "HR", "HR.lo", "HR.hi", "PV1")]
raw.data_subgroups<-raw.data_subgroups[, c("Model", "Variable","N", "HR", "HR.lo", "HR.hi", "PV1")]
raw.data_subgroups_OS<-raw.data_subgroups_OS[, c("Model", "Variable", "N","HR", "HR.lo", "HR.hi", "PV1")]

library(ggplot2)
library(dplyr)


raw.data <- raw.data %>%
  mutate(Variable = paste(Model, Variable, sep = ": "),
         Signif = ifelse(PV1 < 0.05, "*", ""))

raw.data_subgroups <- raw.data_subgroups %>%
  mutate(Variable = paste(Model, Variable, sep = ": "),
         Signif = ifelse(PV1 < 0.05, "*", ""))

raw.data_subgroups_OS <- raw.data_subgroups_OS %>%
  mutate(Variable = paste(Model, Variable, sep = ": "),
         Signif = ifelse(PV1 < 0.05, "*", ""))



data<-raw.data[c(7,8,3,4),]
ForestPlot3 <- function(data, ids2, Title){
  
  library(dplyr); library(tibble); library(forestplot); library(grid)
  
  res <- data[,3:ncol(data)]; colnames(res) <- c("HR1","HR.lo","HR.hi","PV1")
  mhr=res$HR1; mlw=res$HR.lo; mup=res$HR.hi; myp=res$PV1
  myp <- ifelse(myp<0.001,"<0.001",ifelse(myp>0.1,">0.1",formatC(signif(myp,1),format="g")))
  hr <- paste0(formatC(as.numeric(mhr), digits=2, format="f"),
               " (", formatC(as.numeric(mlw), digits=2, format="f"),
               " to ", formatC(as.numeric(mup), digits=2, format="f"), ")")
  
  
  data <- tibble(Feature=ids2,mean=mhr,lower=mlw,upper=mup,hr=hr,coxphP=myp)
  data <- bind_rows(
    tibble(Feature="Feature",mean=NA,lower=NA,upper=NA,hr="HR (95% CI)",coxphP="P-value"),
    tibble(Feature="Hoshida",mean=NA,lower=NA,upper=NA,hr="",coxphP=""),
    data[1:2,],
    tibble(Feature="TP53",mean=NA,lower=NA,upper=NA,hr="",coxphP=""),
    data[3:4,]
  )
  
  forestplot(
    data,
    labeltext = data[, c("Feature","hr","coxphP")],
    # mean = data$mean, lower = data$lower, upper = data$upper,
    graph.pos = 2,
    zero = 1,
    is.summary = data$Feature %in% c("Feature","Hoshida","TP53"),
    xlog = FALSE,
    xlab = paste0("Hazard ratio: ", Title),
    ci.vertices = TRUE,
    hrzl_lines = TRUE,
    boxsize = 0.2,
    line.margin = 0.1,
    mar = unit(rep(0.1, 4), "mm"),
    graphwidth = "auto",
    graphhight = "auto",
    align = rep("l", 3),
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.9),
      title = gpar(cex = 0.9),
      ticks = gpar(cex = 0.9),
      xlab = gpar(cex = 0.9)
    ),
    col = fpColors(box = "black", line = "black")
  ) |> fp_set_zebra_style("#EFEFEF")
}


FP_PR_Hoshida_P53_independent<-ForestPlot3(raw.data[c(7,8,3,4),],c("Hoshida (Int vs Good)","Hoshida (Int vs Good)","TP53 (Int vs Good)","TP53 (Poor vs Good)"),"PFS")
FP_OS_Hoshida_P53_independent<-ForestPlot3(raw.data[c(9,10,5,6),],c("Hoshida (Int vs Good)","Hoshida (Int vs Good)","TP53 (Int vs Good)","TP53 (Poor vs Good)"),"OS")

graphics.off()
# pdf(paste0(figDir,"FP_Progression_Hoshida_P53_cat.pdf"),width=5,height=2)
cairo_pdf(paste0(figDir,"FP_PR_Hoshida_P53_independent.pdf"), width = 6, height = 2.2) 
print(FP_PR_Hoshida_P53_independent)
dev.off()

graphics.off()
# pdf(paste0(figDir,"FP_OS_subgroups.pdf"),width=5,height=2)
cairo_pdf(paste0(figDir,"FP_OS_Hoshida_P53_independent.pdf"), width = 6, height = 2.2) 
print(FP_OS_Hoshida_P53_independent)
dev.off()



ForestPlot4 <- function(data, ids2, Title){
  
  library(dplyr); library(tibble); library(forestplot); library(grid)
  
  # Extract results
  res <- data[, 3:ncol(data)]
  colnames(res) <- c("N","HR1","HR.lo","HR.hi","PV1")
  
  nsample <- data$N
  mhr <- res$HR1; mlw <- res$HR.lo; mup <- res$HR.hi; myp <- res$PV1
  myp <- dplyr::case_when(
    myp < 0.001 ~ "<0.001",
    myp > 0.1   ~ ">0.1",
    TRUE        ~ formatC(signif(myp,1),format="g")
  )
  hr <- paste0(formatC(as.numeric(mhr), digits=2, format="f"),
               " (", formatC(as.numeric(mlw), digits=2, format="f"),
               " to ", formatC(as.numeric(mup), digits=2, format="f"), ")")
  
  # Base data
  data <- tibble(Feature=ids2, N=as.character(nsample), mean=mhr, lower=mlw, upper=mup, hr=hr, coxphP=myp)
  
    # Insert group headers
  data <- bind_rows(
    tibble(Feature="Feature", N="N", mean=NA, lower=NA, upper=NA, hr="HR (95%CI)", coxphP="P-value"),
    tibble(Feature="All", N="", mean=NA, lower=NA, upper=NA, hr="", coxphP=""),
    data[1,],                     # row 1 under All
    tibble(Feature="Bilirubin", N="", mean=NA, lower=NA, upper=NA, hr="", coxphP=""),
    data[2:3,],                   # rows 2-3 under Bilirubin
    tibble(Feature="Platelet", N="", mean=NA, lower=NA, upper=NA, hr="", coxphP=""),
    data[4:5,],                   # rows 4-5 under Platelet
    tibble(Feature="Varices", N="", mean=NA, lower=NA, upper=NA, hr="", coxphP=""),
    data[6:7,]                    # rows 6-7 under Varices
  )
  
  # Define summary rows for bold
  is_summary <- data$Feature %in% c("Feature","All","Bilirubin","Platelet","Varices")
  
  # Plot
  forestplot(
    data,
    labeltext = data[, c("Feature","N","hr","coxphP")],
    # mean = data$mean, lower = data$lower, upper = data$upper,
    graph.pos = 3,
    zero = 1,
    is.summary = is_summary,
    xlog = FALSE,
    xlab = paste0("Hazard ratio: ", Title),
    ci.vertices = TRUE,
    hrzl_lines = TRUE,
    boxsize = 0.2,
    line.margin = 0.1,
    mar = unit(rep(0.1,4),"mm"),
    graphwidth = "auto",
    graphhight = "auto",
    align = c("l","c","l","l"),
    txt_gp = fpTxtGp(label=gpar(cex=0.9), title=gpar(cex=0.9),
                     ticks=gpar(cex=0.9), xlab=gpar(cex=0.9)),
    col = fpColors(box="black", line="black")
  ) |> fp_set_zebra_style("#EFEFEF") -> P2
  
  return(P2)
}


FP_PR_subgroups<-ForestPlot4(raw.data_subgroups,c("All Samples (PFS)","Bilirubin high","Bilirubin normal","Platelet normal","Platelet low","Varices exist","No varices"),"PFS")
FP_OS_subgroups<-ForestPlot4(raw.data_subgroups_OS,c("All Samples (OS)","Bilirubin high","Bilirubin normal","Platelet normal","Platelet low","Varices exist","No varices"),"OS")

graphics.off()
# pdf(paste0(figDir,"FP_Progression_Hoshida_P53_cat.pdf"),width=5,height=2)
cairo_pdf(paste0(figDir,"FP_PR_subgroups.pdf"), width = 6.1, height = 3.2) 
print(FP_PR_subgroups)
dev.off()

graphics.off()
 # pdf(paste0(figDir,"FP_OS_subgroups.pdf"),width=5,height=2)
cairo_pdf(paste0(figDir,"FP_OS_subgroups.pdf"), width = 6.1, height = 3.2) 
print(FP_OS_subgroups)
dev.off()


#############

##### fig2b
#--------------------------
library(ggplot2)
my_comparisons <- list(c( "Good prognosis","Poor prognosis"))
median_score <- median(info$score, na.rm = TRUE)

offset <- 0.2 * diff(range(info$score))   # space between line and text


library(ggplot2)
library(survival)
library(survminer)



plot_score_group2 <- function(
    data,
    group_var,
    group_labels = NULL,
    colors = c("aquamarine4", "grey40","orchid4"),
    xlab = "",
    ylab = "TP53 score",
    title = "",
    add_median_line = TRUE,
    add_labels = FALSE,
    quadrant_labels 
) {
  
  data[[group_var]] <- factor(data[[group_var]])
  
  if (!is.null(group_labels)) {
    levels(data[[group_var]]) <- group_labels
  }
  
  # order scores
  score_sorted <- sort(data$score, na.last = NA)
  
  cut_good  <- score_sorted[60]
  cut_poor  <- score_sorted[161]
  

  median_score <- median(data$score, na.rm = TRUE)
  my_comparisons <- list(c("Good", "Intermediate"),c("Intermediate","Poor"))
  
  p <- ggplot(data, aes_string(x = group_var, y = "score", color = group_var)) +
    # geom_violin(fill = "grey93", color = NA, trim = FALSE, width = 0.7) +
    geom_boxplot(width = 0.6,  alpha = 0.8,outlier.size = 0.5, fill = NA,size=0.8) +
    geom_jitter(width = 0.15,size = 0.3, alpha = 0.7) +
    scale_color_manual(values = colors) +
    # scale_x_discrete(expand = expansion(add = 1))+
     scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    # --- optional median line ---
    { if (add_median_line)
      geom_hline(
        # yintercept = c(cut_good,median_score, cut_poor),  
        yintercept = c(cut_good, cut_poor),
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.8
      )
      
    } +labs(x = xlab, y = ylab, title = title) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 10),       
          axis.text.y = element_text(size = 10),       
          axis.title.x = element_text(size = 12),      
          axis.title.y = element_text(size = 12),      
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"))+ 
    stat_compare_means(
        comparisons = my_comparisons,
        method = "wilcox.test",
        label = "p.signif",
        size = 6, 
        tip.length = 0.03
      )
      
  return(p)
}


fig2bb <- plot_score_group2(
  data = info,
  group_var = "prediction..p.0.05.",
  group_labels = c("Good", "Intermediate","Poor"),
  # colors = c("grey30", "aquamarine4","orchid4"),
  xlab = "Prediction (Hoshida signature)",
  title = "",
  quadrant_labels=c("High TP53\nGood prognosis", "Low TP53\nGood prognosis", "High TP53\nPoor prognosis", "Low TP53\nPoor prognosis")
)

graphics.off()
pdf(paste0(figDir,"box_score_hoshida3cat.pdf"),width=3,height=3)
print(fig2bb)
dev.off()



plot_score_group <- function(
    data,
    group_var,
    group_labels = NULL,
    colors = c("aquamarine4", "orchid4"),
    xlab = "",
    ylab = "TP53 score",
    title = "",
    add_median_line = TRUE,
    add_labels = FALSE,
    quadrant_labels 
) {
  
  data[[group_var]] <- factor(data[[group_var]])
  
  if (!is.null(group_labels)) {
    levels(data[[group_var]]) <- group_labels
  }
  
  median_score <- median(data$score, na.rm = TRUE)
  offset <- 0.18 * diff(range(data$score))   # spacing above/below median
  
  my_comparisons <- list(levels(data[[group_var]]))
  
  p <- ggplot(data, aes_string(x = group_var, y = "score", color = group_var)) +
    
    geom_boxplot(width = 0.6, outlier.size = 0.5, alpha = 0.7, fill = NA) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
    scale_color_manual(values = colors) +
    # scale_x_discrete(expand = expansion(add = 1))+
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    # --- optional median line ---
    { if (add_median_line)
      geom_hline(yintercept = median_score,
                 linetype = "dashed",
                 color = "grey40",
                 linewidth = 0.8)
    } +
    
    # --- quadrant annotations ---
    # annotate("text", x = 0.1, y = 4, 
    #          label = quadrant_labels[1], size = 3.5,  hjust = 0,fontface = "bold") +
    # annotate("text", x = 0.1, y = -2, 
    #          label = quadrant_labels[2], size = 3.5,  hjust = 0,fontface = "bold") +
    # annotate("text", x = 2.9, y = 4, 
    #          label = quadrant_labels[3], size = 3.5,  hjust = 1,fontface = "bold") +
    # annotate("text", x = 2.9, y = -2, 
    #          label = quadrant_labels[4], size = 3.5,  hjust = 1,fontface = "bold") +
    # 
  # --- Wilcoxon p-value ---
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    size = 6, 
    tip.length = 0.03
  ) +
    
    labs(x = xlab, y = ylab, title = title) +
    
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 12),       
          axis.text.y = element_text(size = 12),       
          axis.title.x = element_text(size = 14),      
          axis.title.y = element_text(size = 14), 
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")
          # ,
          # plot.margin = margin(-7, 5, 0, 5)
    )
  
  return(p)
}



info$bilirubin_group <- factor(info$bilirubin...1.0mg.dl,levels = c(0,1),labels = c("<1.0 mg/dL", ">1.0 mg/dL"))

fig_bilirubin <- plot_score_group(data = info,
  group_var = "bilirubin_group",
  colors = c("lightskyblue3", "lightsalmon4"),
  xlab = "Bilirubin level",
  quadrant_labels=c("High TP53\nNormal Bilirubin", "Low TP53\nNormal Bilirubin", "High TP53\nHigh Bilirubin", "Low TP53\nHigh Bilirubin")
)


info$platelet_group <- factor(
  ifelse(info$platelet...100.000.mm3 > 0, "<100k/mm³", ">100k/mm³")
)

fig_platelets <- plot_score_group(
  data = info,
  group_var = "platelet_group",
  colors = c("lightskyblue3", "lightsalmon4"),
  xlab = "Platelet count",
  quadrant_labels=c( "High TP53\nLow Platelet", "Low TP53\nLow Platelet","High TP53\nNormal Platelet", "Low TP53\nNormal Platelet")
)

info$varices_group <- factor(
  info$presence.of.varices,
  levels = c(0, 1),
  labels = c("No varices", "Varices present")
)

fig_varices <- plot_score_group(
  data = info[!is.na(info$presence.of.varices), ],
  group_var = "presence.of.varices",
  group_labels = c("No varices", "Varices present"),
  colors = c("lightskyblue3", "lightsalmon4"),
  xlab = "Presence of varices",
  quadrant_labels=c( "High TP53\nNo varices", "Low TP53\nNo varices","High TP53\nVarices present", "Low TP53\nVarices present")
)


graphics.off()
pdf(paste0(figDir,"box_score_varices.pdf"),width = 2.8,height = 3)
fig_varices
dev.off()

graphics.off()
pdf(paste0(figDir,"box_score_bilirubin.pdf"),width = 2.8,height = 3)
fig_bilirubin
dev.off()

graphics.off()
pdf(paste0(figDir,"box_score_platelets.pdf"),width = 2.8,height = 3)
fig_platelets
dev.off()
##############

################
library(survival)
library(survminer)

info$prediction <- factor(info$prediction, levels = c("Good prognosis", "Poor prognosis"))
info$score_group <- ifelse(info$score > median(info$score, na.rm=TRUE), "High", "Low")
info$score_group <- factor(info$score_group, levels = c("Low", "High"))
info$varices <- factor(info$presence.of.varices,
                       levels = c(0, 1),
                       labels = c("No varices", "Varices"))


info$P53_Score_Prediction <- factor(info$P53_Score_Prediction,
                         levels = c("Good", "Int","Poor"))

info$months.to.hcc<-(info$days.to.hcc)/30.0
fit_var_hcc <- survfit(Surv(months.to.hcc, hcc) ~ P53_Score_Prediction, data = info)



library(survival)
library(survminer)

# KM fit for progression to HCC
fit_var_hcc <- survfit(
  Surv(months.to.hcc, hcc) ~ P53_Score_Prediction,
  data = info
)

# Pairwise log-rank tests
pw <- pairwise_survdiff(
  Surv(months.to.hcc, hcc) ~ P53_Score_Prediction,
  data = info,
  p.adjust.method = "BH"
)

fmt_p <- function(p) {ifelse(p < 0.001,"< 0.001", formatC(signif(p, 1), format = "f", drop0trailing = TRUE))}

p_int_vs_low  <- fmt_p(pw$p.value["Int", "Good"])
p_high_vs_low <- fmt_p(pw$p.value["Poor", "Good"])



library(survival)


cox <- coxph(Surv(months.to.hcc, hcc) ~ P53_Score_Prediction,data = info)
hr <- exp(coef(cox))
# HRs vs reference ("Good")
hr_int  <- sprintf("%.1f", hr["P53_Score_PredictionInt"])
hr_high <- sprintf("%.1f", hr["P53_Score_PredictionPoor"])

# Combine p-values + HRs
label_int_vs_low  <- paste0("Int vs Good: p = ", p_int_vs_low, ", HR = ", hr_int)
label_high_vs_low <- paste0("Poor vs Good: p = ", p_high_vs_low, ", HR = ", hr_high)



# KM plot with pairwise p-values
KM_p53_3cat_progression <- ggsurvplot(
  fit_var_hcc,
  data = info,
  risk.table = FALSE,
  palette = c("grey50", "sandybrown", "salmon3"),
  xlab = "Time (Months)",
  ylab = "Progression free survival",
  legend.title = "",
  legend.labs = c(
    "Low TP53\nscore",
    "Intermediate\nTP53 score",
    "High TP53\nscore"
  ),
  ggtheme = theme_classic() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
)

KM_p53_3cat_progression$plot <- KM_p53_3cat_progression$plot +
  annotate(
    "text",
    x = 0,
    y = 0.1,
    hjust = 0,
    size = 4,
    label = paste0( "n=", length(info$prediction),"\n",
                    label_int_vs_low, "\n",
                    label_high_vs_low)
  )


graphics.off()
pdf(paste0(figDir,"KM_p53_3cat_progression.pdf"),width =4,height = 3.8)
print(KM_p53_3cat_progression$plot)
dev.off()





library(survival)
library(survminer)

info$months.to.death<-(info$days.to.death)/30.0
# Fit KM
fit_var_hcc <- survfit(
  Surv(months.to.death, death) ~ P53_Score_Prediction,
  data = info
)

# Pairwise log-rank tests
pw <- pairwise_survdiff(
  Surv(months.to.death, death) ~ P53_Score_Prediction,
  data = info,
  p.adjust.method = "BH"
)

# Format p-values
fmt_p <- function(p) {ifelse(p < 0.001,"< 0.001", formatC(signif(p, 1), format = "f", drop0trailing = TRUE))}


p_int_vs_low  <- fmt_p(pw$p.value["Int", "Good"])
p_high_vs_low <- fmt_p(pw$p.value["Poor", "Good"])


cox <- coxph(Surv(months.to.death, death) ~ P53_Score_Prediction,data = info)
hr <- exp(coef(cox))
# HRs vs reference ("Good")
hr_int  <- sprintf("%.1f", hr["P53_Score_PredictionInt"])
hr_high <- sprintf("%.1f", hr["P53_Score_PredictionPoor"])

# Combine p-values + HRs
label_int_vs_low  <- paste0("Int vs Good: p = ", p_int_vs_low, ", HR = ", hr_int)
label_high_vs_low <- paste0("Poor vs Good: p = ", p_high_vs_low, ", HR = ", hr_high)


# KM plot with pairwise p-values
KM_p53_3cat_OS <- ggsurvplot(
  fit_var_hcc,
  data = info,
  risk.table = FALSE,
  palette = c("grey50", "sandybrown", "salmon3"),
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  legend.title = "",
  legend.labs = c("Low TP53\nscore", "Intermediate\nTP53 score", "High TP53\nscore"),
  ggtheme = theme_classic()+
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
)

KM_p53_3cat_OS$plot <- KM_p53_3cat_OS$plot +
  annotate(
    "text",
    x = 0,
    y = 0.1,
    hjust = 0,
    size = 4,
    label = paste0( "n=", length(info$prediction),"\n",
                    label_int_vs_low, "\n",
                    label_high_vs_low)
  )

# Save
graphics.off()
pdf(paste0(figDir, "KM_p53_3cat_OS.pdf"), width = 4, height = 3.8)
print(KM_p53_3cat_OS$plot)
dev.off()





info1<-info[which(info$prediction..p.0.05.=="Intermediate prognosis"),]
info1<-info1[which(info1$P53_Score_Prediction!="Int"),]
fit_var_hcc <- survfit(Surv(months.to.death, death) ~ P53_Score_Prediction, data = info1)

cox <- coxph(Surv(months.to.death, death) ~ P53_Score_Prediction, data = info1)
hr <- exp(coef(cox))
# HRs vs reference ("Good")
hr_high <- sprintf("%.1f", hr["P53_Score_PredictionPoor"])

# Combine p-values + HRs
# label_int_vs_low  <- paste0("Int vs Good: p = ", p_int_vs_low, ", HR = ", hr_int)
label_high_vs_low <- paste0("HR = ", hr_high)


KM_p53_2cat_OS_in_intermediate<-ggsurvplot(
  fit_var_hcc,
  data = info1,
  pval = TRUE,
  risk.table = F,
  xlab = "Time (Months)",
  palette = c("grey50", "salmon3"),
  title = "OS within intermediate\n Hoshida prognosis",
  legend.title = "",
  legend.labs = c("Low TP53 score",  "High TP53 score"),
  ggtheme = theme_classic() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      plot.title   = element_text(size = 14, hjust = 0.5) # center title
    )
)

KM_p53_2cat_OS_in_intermediate$plot <- KM_p53_2cat_OS_in_intermediate$plot +
  annotate(
    "text",
    x = 5,
    y = 0.35,
    hjust = 0,
    size = 4.5,
    label = paste0( "n=", length(info1$prediction),"\n",
                    label_high_vs_low)
  )

graphics.off()
pdf(paste0(figDir,"KM_p53_2cat_OS_in_intermediate.pdf"),width=4,height = 4)
print(KM_p53_2cat_OS_in_intermediate$plot)
dev.off()
 


#################
##Now combine the Tp53 score with other measures
library(survival)
library(survminer)

# Dichotomize TP53
median_score <- median(info$score, na.rm = TRUE)
info$TP53_group <- factor(ifelse(info$score > median_score, "High TP53", "Low TP53"))
info$prediction <- factor(info$prediction, levels = c("Good prognosis", "Poor prognosis"))
# Clinical variables
info$bilirubin_group <- factor(info$bilirubin...1.0mg.dl,
                               levels = c(0,1),
                               labels = c("≤1.0 mg/dL", ">1.0 mg/dL"))
info$platelet_group <- factor(ifelse(info$platelet...100.000.mm3 > 0,  "≤100k/mm³",">100k/mm³"))

info$varices_group <- factor(info$presence.of.varices,
                             levels = c(0,1),
                             labels = c("No varices", "Varices"))

info$bilirubin_TP53 <- interaction(info$bilirubin_group, info$TP53_group, sep="_")
info$platelet_TP53 <- interaction(info$platelet_group, info$TP53_group, sep="_")
info$varices_TP53  <- interaction(info$varices_group, info$TP53_group, sep="_")
info$Prediction_TP53  <- interaction(info$prediction, info$TP53_group, sep="_")



###################
