library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db", character.only=TRUE)

#
# Loss18 group
#
samples <- read.table("../samples.txt", sep="\t", header=T)
samples$loss18 <- relevel(factor(samples$loss18, levels=c("1", "0")), ref="0")

files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ loss18)
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

resShrink <- lfcShrink(dds, res=res, type="apeglm", coef="loss18_1_vs_0")
resOrderedShrink <- resShrink[order(resShrink$pvalue),]

write.table(resOrderedShrink, "loss18_group.txt")

symbolmap=mapIds(org.Hs.eg.db, rownames(res), keytype="ENSEMBL", column="SYMBOL")

EnhancedVolcano(resShrink, rownames(res), x="log2FoldChange", y="padj", pCutoff=0.05)

#
# 3A
#

EnhancedVolcano(resShrink, symbolmap, x="log2FoldChange", y="padj", pCutoff=0.05)

resnona <- na.omit(resShrink)
enrich <- enrichGO(rownames(resnona[resnona$padj<0.05,]), org.Hs.eg.db, keyType="ENSEMBL", ont="ALL", universe=rownames(resnona))

enrichup <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange>0.5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))
enrichdown <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange< -0.5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))

enrichdown <- mutate(enrichdown, q=-log(p.adjust, base=10))
enrichup <- mutate(enrichup, q=-log(p.adjust, base=10))
enrich <- mutate(enrich, q=-log(p.adjust, base=10))
write.table(enrichdown, "loss18_group_GO_down.txt")
write.table(enrichup, "loss18_group_GO_up.txt")
write.table(enrich, "loss18_group_GO.txt")

#
#3D
#
barplot(simplify(enrichdown), x="q")

#
#3E
#
barplot(simplify(enrichup), x="q")
#barplot(simplify(enrich), x="q")

#
# MultiCNV group
#
library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db", character.only=TRUE)

samples <- read.table("../samples.txt", sep="\t", header=T)
samples$multi <- relevel(factor(samples$multi, levels=c("1", "0")), ref="0")

files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ multi)
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

resShrink <- lfcShrink(dds, res=res, type="apeglm", coef="multi_1_vs_0")
resOrderedShrink <- resShrink[order(resShrink$pvalue),]
write.table(resOrderedShrink, "multi_group.txt")

symbolmap=mapIds(org.Hs.eg.db, rownames(resShrink), keytype="ENSEMBL", column="SYMBOL")

EnhancedVolcano(resShrink, rownames(resShrink), x="log2FoldChange", y="padj", pCutoff=0.05)

#
#3B
#
EnhancedVolcano(resShrink, symbolmap, x="log2FoldChange", y="padj", pCutoff=0.05)

resnona <- na.omit(resShrink)
enrich <- enrichGO(rownames(resnona[resnona$padj<0.05,]), org.Hs.eg.db, keyType="ENSEMBL", ont="ALL", universe=rownames(resnona))
enrichup <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange>0-5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))
enrichdown <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange< -0.5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))

enrichdown <- mutate(enrichdown, q=-log(p.adjust, base=10))
enrichup <- mutate(enrichup, q=-log(p.adjust, base=10))
enrich <- mutate(enrich, q=-log(p.adjust, base=10))
write.table(enrichdown, "multicnv_group_GO_down.txt")
write.table(enrichup, "multicnv_group_GO_up.txt")
#barplot(simplify(enrichdown), x="q") #none found

#barplot(simplify(enrichup), x="q") #none found
#barplot(simplify(enrich), x="q")

#
# Silent group
#
library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db", character.only=TRUE)

# All_samples
samples <- read.table("../samples.txt", sep="\t", header=T)
samples$silent <- relevel(factor(samples$silent, levels=c("1", "0")), ref="0")

files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ silent)
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]


resShrink <- lfcShrink(dds, res=res, type="apeglm", coef="silent_1_vs_0")
resOrderedShrink <- resShrink[order(resShrink$pvalue),]

write.table(resOrderedShrink, "silent_group.txt")
symbolmap=mapIds(org.Hs.eg.db, rownames(resShrink), keytype="ENSEMBL", column="SYMBOL")

EnhancedVolcano(resShrink, rownames(resShrink), x="log2FoldChange", y="padj", pCutoff=0.05)

#
#3C
#
EnhancedVolcano(resShrink, symbolmap, x="log2FoldChange", y="padj", pCutoff=0.05)

resnona <- na.omit(resShrink)

enrich <- enrichGO(rownames(resnona[resnona$padj<0.05,]), org.Hs.eg.db, keyType="ENSEMBL", ont="ALL", universe=rownames(resnona))
enrichup <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange>0.5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))
enrichdown <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange< -0.5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))

enrichdown <- mutate(enrichdown, q=-log(p.adjust, base=10))
enrichup <- mutate(enrichup, q=-log(p.adjust, base=10))
enrich <- mutate(enrich, q=-log(p.adjust, base=10))
write.table(enrichdown, "silent_group_GO_down.txt")
write.table(enrichup, "silent_group_GO_up.txt")
#barplot(simplify(enrichdown), x="q") #none found
#
#3F
#
barplot(simplify(enrichup), x="q")
#barplot(simplify(enrich), x="q")