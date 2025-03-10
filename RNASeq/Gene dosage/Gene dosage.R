library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db", character.only=TRUE)
library(limma)

#
# CHR 18
#
samples <- read.table("../samples.txt", sep="\t", header=T)
samples$chr18_loss <- relevel(factor(samples$chr18_loss, levels=c("yes", "no")), ref="no")

files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample
tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ chr18_loss)
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

resShrink <- lfcShrink(dds, res=res, type="apeglm", coef="chr18_loss_yes_vs_no")
resOrderedShrink <- resShrink[order(resShrink$pvalue),]
write.table(resOrderedShrink, "chr18.txt")

symbolmap=mapIds(org.Hs.eg.db, rownames(res), keytype="ENSEMBL", column="SYMBOL")

EnhancedVolcano(resShrink, rownames(res), x="log2FoldChange", y="padj", pCutoff=0.05)

#
#
# Figure 4A
#
#
EnhancedVolcano(resShrink, symbolmap, x="log2FoldChange", y="padj", pCutoff=0.05)

#
# Figure 4B
#
genes_18q <- c("ENSG00000122490","ENSG00000264247","ENSG00000101546","ENSG00000197971","ENSG00000141759","ENSG00000060069","ENSG00000170677","ENSG00000166573","ENSG00000166377","ENSG00000287693","ENSG00000075336","ENSG00000101544","ENSG00000215421","ENSG00000133313","ENSG00000166479","ENSG00000176225","ENSG00000226742")
name_map <- data.frame(read.table("ensembl112mart.txt", sep="\t", header=T))
res18q <- res[rownames(res) %in% genes_18q,]
res18q <- res18q[order(res18q$log2FoldChange),]
res18q$x <- c(1:17)
res18q$names <- rownames(res18q)
res18q$symbol <- mapIds(org.Hs.eg.db, res18q$names, keytype="ENSEMBL", column="SYMBOL")

ggplot(data=res18q, aes(y=log2FoldChange, x=x, label=symbol, color=padj)) + geom_point() + 
  geom_text_repel() + xlab("") + theme_classic()
#
# CHR 11
#
samples <- read.table("../samples.txt", sep="\t", header=T)
samples$chr11_loss <- relevel(factor(samples$chr11_loss, levels=c("yes", "no")), ref="no")

files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample
tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ chr11_loss)
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

resShrink <- lfcShrink(dds, res=res, type="apeglm", coef="chr11_loss_yes_vs_no")
resOrderedShrink <- resShrink[order(resShrink$pvalue),]
write.table(resOrderedShrink, "chr11.txt")

symbolmap=mapIds(org.Hs.eg.db, rownames(res), keytype="ENSEMBL", column="SYMBOL")

EnhancedVolcano(resShrink, rownames(res), x="log2FoldChange", y="padj", pCutoff=0.05)

#
# Figure 4C
#
EnhancedVolcano(resShrink, symbolmap, x="log2FoldChange", y="padj", pCutoff=0.05, selectLab = c("ATM"))


