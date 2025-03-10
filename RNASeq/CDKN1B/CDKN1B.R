library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db", character.only=TRUE)

#
# CDKN1B
#
samples <- read.table("../samples.txt", sep="\t", header=T)
samples$CDKN1B <- relevel(factor(samples$CDKN1B, levels=c("mutated", "wt")), ref="wt")

files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample
tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ CDKN1B)
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

resShrink <- lfcShrink(dds, res=res, type="apeglm", coef="CDKN1B_mutated_vs_wt")
resOrderedShrink <- resShrink[order(resShrink$pvalue),]
write.table(resOrderedShrink, "CDKN1B.txt")

symbolmap=mapIds(org.Hs.eg.db, rownames(res), keytype="ENSEMBL", column="SYMBOL")

EnhancedVolcano(resShrink, rownames(res), x="log2FoldChange", y="padj", pCutoff=0.05)
EnhancedVolcano(resShrink, symbolmap, x="log2FoldChange", y="padj", pCutoff=0.05)
