library(DESeq2)
library(tximport)
library(gplots)
library(ggplot2)

# Read in data
samples <- read.table("../samples.txt", sep="\t", header=T)
files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample
tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ 1)

#Retain features with at least 5 counts in at least 5 (i.e. 25%) of the samples
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]

dds <- DESeq(ddsTxi)

rld <- vst(dds)

topgenes_n <- dim(rld)[1]*0.10 #10% most variable
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),topgenes_n)

data <- assay(rld[topVarGenes])


#
# Replace the sample names with the names used for publication
#
colnames_new <- c(1:20)
for (i in c(1:20)) {
  colnames_new[i] <- strsplit(colnames(data)[i], "-")[[1]][3]
}
colnames(data) <- colnames_new
colnames(data)[colnames(data)=="10110"]="25P"
colnames(data)[colnames(data)=="1985"]="2M"
colnames(data)[colnames(data)=="2294"]="6P"
colnames(data)[colnames(data)=="3105"]="27M"
colnames(data)[colnames(data)=="3243"]="30P"
colnames(data)[colnames(data)=="3245"]="30M"
colnames(data)[colnames(data)=="3425"]="28P"
colnames(data)[colnames(data)=="3992"]="17P"
colnames(data)[colnames(data)=="4017"]="17M"
colnames(data)[colnames(data)=="4059"]="18M"
colnames(data)[colnames(data)=="4407"]="14M"
colnames(data)[colnames(data)=="4410"]="14P"
colnames(data)[colnames(data)=="4783"]="3P"
colnames(data)[colnames(data)=="5028"]="3M"
colnames(data)[colnames(data)=="5184"]="26P"
colnames(data)[colnames(data)=="5318"]="26M"
colnames(data)[colnames(data)=="5662"]="10P"
colnames(data)[colnames(data)=="5962"]="5P"
colnames(data)[colnames(data)=="5963"]="5M"
colnames(data)[colnames(data)=="6003"]="12P"
colnames(data)[colnames(data)=="6464"]="23M"
colnames(data)[colnames(data)=="6483"]="16P"
colnames(data)[colnames(data)=="6489"]="7P"
colnames(data)[colnames(data)=="6493"]="7M"
colnames(data)[colnames(data)=="6619"]="24P"
colnames(data)[colnames(data)=="6683"]="24M"
colnames(data)[colnames(data)=="6733"]="21P"
colnames(data)[colnames(data)=="6782"]="15M"
colnames(data)[colnames(data)=="6852"]="20P"
colnames(data)[colnames(data)=="6884"]="13P"
colnames(data)[colnames(data)=="7189"]="8P"
colnames(data)[colnames(data)=="7446"]="9M"
colnames(data)[colnames(data)=="8161"]="22P"
colnames(data)[colnames(data)=="8162"]="22M"
colnames(data)[colnames(data)=="8344"]="1M"
colnames(data)[colnames(data)=="849"]="29M"
colnames(data)[colnames(data)=="9146"]="11P"
colnames(data)[colnames(data)=="9147"]="11M"
colnames(data)[colnames(data)=="9397"]="4M"
colnames(data)[colnames(data)=="9733"]="19M"

#
# Generate 4 heatmaps including different information
# These are subsequently merged manually for the manuscript figure
#
rsc_group <- samples$group
rsc_group[rsc_group=="multi"] <- "#E0D0C1"
rsc_group[rsc_group=="loss18"] <- "#601700"
rsc_group[rsc_group=="silent"] <- "#A76D60"
rsc_survival <- samples$survival
rsc_survival[rsc_survival=="long"] <- "#000000"
rsc_survival[rsc_survival=="short"] <- "#999999"
rsc_site <- samples$site
rsc_site[rsc_site=="primary"] <- "#0000AA"
rsc_site[rsc_site=="metastasis"] <- "#00AA00"
rsc_cdkn1b <- samples$CDKN1B
rsc_cdkn1b[rsc_cdkn1b=="wt"] <- "#AAAAAA"
rsc_cdkn1b[rsc_cdkn1b=="mutated"] <- "#AA0000"

data <- t(scale(t(data)))

heatmap.2(t(data), trace="none", 
          density.info ="none", RowSideColors=rsc_group, col=colorRampPalette(c('green',"black",'red'))(32),
          labCol=NA, dendrogram="row")

heatmap.2(t(data), trace="none", 
          density.info ="none", RowSideColors=rsc_survival, col=colorRampPalette(c('green',"black",'red'))(32),
          labCol=NA, dendrogram="row")

heatmap.2(t(data), trace="none", 
          density.info ="none", RowSideColors=rsc_site, col=colorRampPalette(c('green',"black",'red'))(32),
          labCol=NA, dendrogram="row")

heatmap.2(t(data), trace="none", 
          density.info ="none", RowSideColors=rsc_cdkn1b, col=colorRampPalette(c('green',"black",'red'))(32),
          labCol=NA, dendrogram="row", cexRow = 2)

#
# Consensus clustering
#
library(ConsensusClusterPlus)
consensusresults = ConsensusClusterPlus(data,maxK=6,reps=10000,pItem=0.8,pFeature=0.8,title="consensusclustering",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

heatmap.2(t(data), trace="none", 
          density.info ="none", RowSideColors=rsc_group, col=colorRampPalette(c('green',"black",'red'))(32),
          labCol=NA, dendrogram="row", Rowv=as.dendrogram(consensusresults[[4]]$consensusTree))

#
# Principal components analysis
#x
dds <- DESeq(ddsTxi)
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("group"))
plotPCA(vsd, intgroup = c("site"))

