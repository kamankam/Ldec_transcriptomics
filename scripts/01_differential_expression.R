## read DGE
## read gene expression raw counts
counts = read.table('./ldec_dge_counts.tsv', header = T, row.names = 1)
## read the samples data
phen   = read.table('./Design_table.tsv', header = T, row.names = 1)
## read the genes annotation data
gff_anno <- read.table('./ldec_gff_annotated.tsv', header = T, sep = '\t',
                       quote = '',
                       stringsAsFactors = F)
## extract unique annotations
gen_anno <- unique.data.frame( gff_anno[,c('gene_id', 'product')] )
row.names( gen_anno ) <- gen_anno$gene_id

## order expression according to phen 
counts <- counts[,row.names(phen)]
#View(counts[,row.names(phen)])

## Analyzing RNA-seq data with DESeq2
## according to https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
require(DESeq2)

## make the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = phen,
                              design= ~ Tissue + Treat) ## the design formula
## make the DEG test
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
## extract the results for a particular comparison
res <- results(dds, name="Treat_B.bassiana_vs_0")
# or to shrink log fold changes association with condition:
## shrinking log fold changes helps to lower FPs
resbb  <- lfcShrink(dds, coef="Treat_B.bassiana_vs_0", type="apeglm")
resmr  <- lfcShrink(dds, coef="Treat_M.robertsii_vs_0", type="apeglm")
resfat <- lfcShrink(dds, coef="Tissue_Haemocytes_vs_FatBody", type="apeglm")

res_bb_df <- as.data.frame(resbb)
res_bb_df <- cbind(res_bb_df, gen_anno[row.names(res_bb_df), ])
#View(res_df)

res_mr_df <- as.data.frame(resmr)
res_mr_df <- cbind(res_mr_df, gen_anno[row.names(res_mr_df), ])
#View(res_df2)

res_ft_df <- as.data.frame(resfat)
res_ft_df <- cbind(res_ft_df, gen_anno[row.names(res_ft_df), ])

## MAplot
plotMA(resfat, main = "Tissue_Haemocytes_vs_FatBody")
plotMA(resbb, main = "B.bassiana vs. control" )
plotMA(resmr, main = "M.robertsii vs. control")

## extract the transformed counts
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
#head(assay(vsd), 3)

## make PCA plot of the samples
#DESeq2::plotPCA(vsd, intgroup=c("Tissue", "Treat")) + theme(aspect.ratio=1/1)
## plot the dispersion estimates
#plotDispEsts(dds)

## customize the PCA plot of the samples
pcaData <- plotPCA(vsd, intgroup=c("Tissue", "Treat"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Treat, shape=Tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme(aspect.ratio=1/1)
## haemocytes:

## order expression according to phen 
counts_hae <- counts[,phen$Tissue == 'Haemocytes']

## make the DESeqDataSet object
dds_hae <- DESeqDataSetFromMatrix(countData = counts_hae,
                                  colData = phen[phen$Tissue == 'Haemocytes',],
                                  design= ~Treat) ## the design formula

## make the DEG test
dds_hae <- DESeq(dds_hae)
resultsNames(dds_hae) # lists the coefficients

## extract the results for a particular comparison
#res_hae_bb <- results(dds_hae, name="Treat_B.bassiana_vs_0")

# or to shrink log fold changes association with condition:
## shrinking log fold changes helps to lower FPs
res_hae_bb <- lfcShrink(dds_hae, coef="Treat_B.bassiana_vs_0", type="apeglm")
res_hae_mr <- lfcShrink(dds_hae, coef="Treat_M.robertsii_vs_0", type="apeglm")

res_hae_bb_df <- as.data.frame(res_hae_bb)
res_hae_bb_df <- cbind(res_hae_bb_df, gen_anno[row.names(res_hae_bb_df), ])
#View(res_hae_bb_df)

res_hae_mr_df <- as.data.frame(res_hae_mr)
res_hae_mr_df <- cbind(res_hae_mr_df, gen_anno[row.names(res_hae_mr_df), ])
#View(res_hae_mr_df)

## MAplot
plotMA(res_hae_mr, main = "Haemocytes. M.robertsii vs. control")
plotMA(res_hae_bb, main = "Haemocytes. B.bassiana vs. control" )

## extract the transformed counts
vsd_hae <- vst(dds_hae, blind=FALSE)
rld_hae <- rlog(dds_hae, blind=FALSE)

## make PCA plot of the samples
DESeq2::plotPCA(vsd_hae, intgroup=c("Treat")) + theme(aspect.ratio=1/1)

## plot the dispersion estimates
#plotDispEsts(dds_hae)

## fat body without C17:

## order expression according to phen 
counts_fat <- counts[,phen$Tissue != 'Haemocytes']
counts_fat <- counts_fat[, names(counts_fat) != 'C17']

## make the DESeqDataSet object
phen2 <- phen[phen$Tissue != 'Haemocytes',]
phen2 <- phen2[row.names(phen2) != 'C17',]

dds_fat <- DESeqDataSetFromMatrix(countData = counts_fat,
                                  colData = phen2,
                                  design= ~Treat) ## the design formula

## make the DEG test
dds_fat <- DESeq(dds_fat)
resultsNames(dds_fat) # lists the coefficients

## extract the results for a particular comparison
#res_fat_bb <- results(dds_fat, name="Treat_B.bassiana_vs_0")

# or to shrink log fold changes association with condition:
## shrinking log fold changes helps to lower FPs
res_fat2_bb <- lfcShrink(dds_fat, coef="Treat_B.bassiana_vs_0", type="apeglm")
res_fat2_mr <- lfcShrink(dds_fat, coef="Treat_M.robertsii_vs_0", type="apeglm")

res_fat2_bb_df <- as.data.frame(res_fat2_bb)
res_fat2_bb_df <- cbind(res_fat2_bb_df, gen_anno[row.names(res_fat2_bb_df), ])
#View(res_fat2_bb_df)

res_fat2_mr_df <- as.data.frame(res_fat2_mr)
res_fat2_mr_df <- cbind(res_fat2_mr_df, gen_anno[row.names(res_fat2_mr_df), ])
#View(res_fat2_mr_df)

## MAplot
plotMA(res_fat2_mr, main = "Fat body (without C17). M.robertsii vs. control")
plotMA(res_fat2_bb, main = "Fat body (without C17). B.bassiana vs. control" )

## extract the transformed counts
vsd_fat <- vst(dds_fat, blind=FALSE)
rld_fat <- rlog(dds_fat, blind=FALSE)

## make PCA plot of the samples
DESeq2::plotPCA(vsd_fat, intgroup=c("Treat")) + theme(aspect.ratio=1/1)

## plot the dispersion estimates
##plotDispEsts(dds_fat)

#B. bassiana vs M. robertsii

bbas_vs_mrob_all <- results(
  dds,
  contrast = c('Treat', 'B.bassiana', 'M.robertsii')
)

bbas_vs_mrob_hae <- results(
  dds_hae,
  contrast = c('Treat', 'B.bassiana', 'M.robertsii')
)

bbas_vs_mrob_fat <- results(
  dds_fat,
  contrast = c('Treat', 'B.bassiana', 'M.robertsii')
)

bbas_vs_mrob_all2 <- lfcShrink(
  dds,
  contrast = c('Treat', 'B.bassiana', 'M.robertsii'),
  type = 'ashr'
)

bbas_vs_mrob_hae2 <- lfcShrink(
  dds_hae,
  contrast = c('Treat', 'B.bassiana', 'M.robertsii'),
  type = 'ashr'
)

bbas_vs_mrob_fat2 <- lfcShrink(
  dds_fat,
  contrast = c('Treat', 'B.bassiana', 'M.robertsii'),
  type = 'ashr'
)

bbas_vs_mrob_all2_df <- as.data.frame(bbas_vs_mrob_all2)
bbas_vs_mrob_all2_df <- cbind(
  bbas_vs_mrob_all2_df,
  gen_anno[row.names(bbas_vs_mrob_all2_df), ]
)

bbas_vs_mrob_hae2_df <- as.data.frame(bbas_vs_mrob_hae2)
bbas_vs_mrob_hae2_df <- cbind(
  bbas_vs_mrob_hae2_df,
  gen_anno[row.names(bbas_vs_mrob_hae2_df), ]
)

bbas_vs_mrob_fat2_df <- as.data.frame(bbas_vs_mrob_fat2)
bbas_vs_mrob_fat2_df <- cbind(
  bbas_vs_mrob_fat2_df,
  gen_anno[row.names(bbas_vs_mrob_fat2_df), ]
)
