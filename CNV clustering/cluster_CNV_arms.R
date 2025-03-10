library(gplots)
library(RColorBrewer)
library(ggplot2)

data <- data.frame(read.table("chr_arms_data.txt", header=T, sep="\t", row.names=1))
CNV_matrix <- data.frame(data)

rownames(CNV_matrix)[rownames(CNV_matrix)=="10110"]="25P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="1985"]="2M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="2294"]="6P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="3105"]="27M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="3243"]="30P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="3245"]="30M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="3425"]="28P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="3992"]="17P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="4017"]="17M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="4059"]="18M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="4407"]="14M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="4410"]="14P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="4783"]="3P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="5028"]="3M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="5184"]="26P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="5318"]="26M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="5662"]="10P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="5962"]="5P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="5963"]="5M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6003"]="12P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6464"]="23M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6483"]="16P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6489"]="7P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6493"]="7M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6619"]="24P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6683"]="24M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6733"]="21P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6782"]="15M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6852"]="20P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="6884"]="13P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="7189"]="8P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="7446"]="9M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="8161"]="22P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="8162"]="22M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="8344"]="1M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="849"]="29M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="9146"]="11P"
rownames(CNV_matrix)[rownames(CNV_matrix)=="9147"]="11M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="9397"]="4M"
rownames(CNV_matrix)[rownames(CNV_matrix)=="9733"]="19M"
colors <- as.character(c("red", "white", "blue"))

rsc_CDKN1B <- rownames(CNV_matrix)
rsc_CDKN1B[rsc_CDKN1B=="26M"] <- "Mutated"
rsc_CDKN1B[rsc_CDKN1B=="20P"] <- "Mutated"
rsc_CDKN1B[rsc_CDKN1B=="9M"] <- "Mutated"
rsc_CDKN1B[rsc_CDKN1B!="Mutated"] <- "WT"

rsc_CDKN1B_cols <- rsc_CDKN1B

rsc_CDKN1B_cols[rsc_CDKN1B_cols=="Mutated"] <- "#00A676"
rsc_CDKN1B_cols[rsc_CDKN1B_cols=="WT"] <- "#FFFFFF"

manhattan <- function(x) {return(dist(x, method="manhattan"))}

heatmap.2(as.matrix(CNV_matrix), Colv=NA, col=colors, scale="none", 
          distfun=manhattan, trace="none",
          RowSideColors = rsc_CDKN1B_cols)

rsc_cluster <- rownames(CNV_matrix)
rsc_cluster[rsc_cluster=="26M"] <- "Loss18"
rsc_cluster[rsc_cluster=="12P"] <- "Loss18"
rsc_cluster[rsc_cluster=="11P"] <- "Loss18"
rsc_cluster[rsc_cluster=="17P"] <- "Loss18"
rsc_cluster[rsc_cluster=="17M"] <- "Loss18"
rsc_cluster[rsc_cluster=="26P"] <- "Loss18"
rsc_cluster[rsc_cluster=="5M"] <- "Loss18"
rsc_cluster[rsc_cluster=="10P"] <- "Loss18"
rsc_cluster[rsc_cluster=="7P"] <- "Loss18"
rsc_cluster[rsc_cluster=="16P"] <- "Loss18"
rsc_cluster[rsc_cluster=="21P"] <- "Loss18"
rsc_cluster[rsc_cluster=="3P"] <- "Loss18"
rsc_cluster[rsc_cluster=="13P"] <- "Loss18"
rsc_cluster[rsc_cluster=="5P"] <- "Loss18"

rsc_cluster[rsc_cluster=="19M"] <- "Silent"
rsc_cluster[rsc_cluster=="2M"] <- "Silent"
rsc_cluster[rsc_cluster=="11M"] <- "Silent"
rsc_cluster[rsc_cluster=="27M"] <- "Silent"
rsc_cluster[rsc_cluster=="29M"] <- "Silent"
rsc_cluster[rsc_cluster=="14P"] <- "Silent"
rsc_cluster[rsc_cluster=="28P"] <- "Silent"
rsc_cluster[rsc_cluster=="8P"] <- "Silent"
rsc_cluster[rsc_cluster=="30P"] <- "Silent"
rsc_cluster[rsc_cluster=="25P"] <- "Silent"
rsc_cluster[rsc_cluster=="30M"] <- "Silent"
rsc_cluster[rsc_cluster=="24P"] <- "Silent"
rsc_cluster[rsc_cluster=="24M"] <- "Silent"

rsc_cluster[rsc_cluster=="3M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="1M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="22P"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="9M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="15M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="22M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="20P"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="6P"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="14M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="18M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="23M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="4M"] <- "MultiCNV"
rsc_cluster[rsc_cluster=="7M"] <- "MultiCNV"

rsc_cluster_cols <- rsc_cluster
rsc_cluster_cols[rsc_cluster_cols=="MultiCNV"] <- "#E0D0C1"
rsc_cluster_cols[rsc_cluster_cols=="Silent"] <- "#A76D60"
rsc_cluster_cols[rsc_cluster_cols=="Loss18"] <- "#601700"

heatmap.2(as.matrix(CNV_matrix), Colv=NA, col=colors, scale="none", 
          distfun=manhattan, trace="none",
          RowSideColors = rsc_cluster_cols)

#
# Just a sanity check to ensure that we do not get wildly different results 
# if we use euclidean distance
#
heatmap.2(as.matrix(CNV_matrix), Colv=NA, col=colors, scale="none", 
          trace="none",
          RowSideColors = rsc_cluster_cols)


rsc_survival <- rownames(CNV_matrix)
rsc_survival[rsc_survival=="1M"] <- "Short"
rsc_survival[rsc_survival=="2M"] <- "Short"
rsc_survival[rsc_survival=="3P"] <- "Short"
rsc_survival[rsc_survival=="3M"] <- "Short"
rsc_survival[rsc_survival=="4M"] <- "Short"
rsc_survival[rsc_survival=="5P"] <- "Short"
rsc_survival[rsc_survival=="5M"] <- "Short"
rsc_survival[rsc_survival=="6P"] <- "Short"
rsc_survival[rsc_survival=="7P"] <- "Short"
rsc_survival[rsc_survival=="7M"] <- "Short"
rsc_survival[rsc_survival=="8P"] <- "Short"
rsc_survival[rsc_survival=="9M"] <- "Short"
rsc_survival[rsc_survival=="10P"] <- "Short"
rsc_survival[rsc_survival=="11P"] <- "Short"
rsc_survival[rsc_survival=="11M"] <- "Short"
rsc_survival[rsc_survival=="12P"] <- "Short"
rsc_survival[rsc_survival=="13P"] <- "Short"
rsc_survival[rsc_survival=="14P"] <- "Short"
rsc_survival[rsc_survival=="14M"] <- "Short"
rsc_survival[rsc_survival!="Short"] <- "Long"

rsc_survival_cols <- rsc_survival
rsc_survival_cols[rsc_survival_cols=="Long"] <- "#000000"
rsc_survival_cols[rsc_survival_cols=="Short"] <- "#999999"

heatmap.2(as.matrix(CNV_matrix), Colv=NA, col=colors, scale="none", 
          distfun=manhattan, trace="none",
          RowSideColors = rsc_survival_cols)

affected_arms <- rowSums(abs(CNV_matrix))
affected_arms_data <- data.frame(affected_arms, rsc_CDKN1B, rsc_cluster, rsc_survival)

ggplot(affected_arms_data, aes(x=rsc_cluster, y=affected_arms)) + geom_boxplot() + 
  theme_classic() + xlab("CNV cluster") + ylab("Number of affected chromosome arms")

ggplot(affected_arms_data, aes(x=rsc_survival, y=affected_arms)) + geom_boxplot() + 
  theme_classic() + xlab("Survival group") + ylab("Number of affected chromosome arms")
