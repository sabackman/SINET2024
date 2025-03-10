#
# Site
#
library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db", character.only=TRUE)

# All_samples
samples <- read.table("../samples.txt", sep="\t", header=T)
samples$site <- relevel(factor(samples$site, levels=c("metastasis", "primary")), ref="primary")

files <- file.path("../salmon_output", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.table(file.path("../salmon_output", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ site)
keep <- rowSums(counts(ddsTxi) >= 5) >= 5
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

resShrink <- lfcShrink(dds, res=res, type="apeglm", coef="site_metastasis_vs_primary")
resOrderedShrink <- resShrink[order(resShrink$pvalue),]
write.table(resOrderedShrink, "site.txt")

symbolmap=mapIds(org.Hs.eg.db, rownames(resShrink), keytype="ENSEMBL", column="SYMBOL")

EnhancedVolcano(resShrink, rownames(resShrink), x="log2FoldChange", y="padj", pCutoff=0.05)
EnhancedVolcano(resShrink, symbolmap, x="log2FoldChange", y="padj", pCutoff=0.05)

resnona <- na.omit(resShrink)

enrichup <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange>0.5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))
enrichdown <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange< -0.5,]), org.Hs.eg.db, keyType="ENSEMBL", ont="BP", universe=rownames(resnona))

enrichdown <- mutate(enrichdown, q=-log(p.adjust, base=10))
enrichup <- mutate(enrichup, q=-log(p.adjust, base=10))
write.table(enrichdown, "site_GO_down.txt")
write.table(enrichup, "site_GO_up.txt")
#barplot(simplify(enrichdown), x="q") #none found
barplot(simplify(enrichup), x="q")
