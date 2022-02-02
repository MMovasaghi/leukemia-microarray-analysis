################## In the name of Allah #######################
#                                                             #
#          Introduction to Bioinformatics, bahman 1400        #
#                    Sharif University                        #
#                                                             #
#  Name:    Mohammad Hosein Movasaghinia                      #
#  Uni-num: 400200919                                         #
#                                                             #
##############################################################*

#### Project ####

library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(gplots)
library(ggplot2)

#### Set working directory ####
curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(sub(paste0("/",sub("(.+)/","",curD)),"",curD))

#### Set global variables ####
series <- "GSE48558"
platform <- "GPL6244"

#### Load data ####
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = 'data/')

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#### Set group ####
gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               "00000000000000000000")

sml <- strsplit(gsms, split="")[[1]]

#### filter by X and add Group to data ####
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
gs <- factor(sml)
groups <- make.names(c("normal", "test"))
levels(gs) <- groups
gset$group <- gs

sname <- gset$source_name_ch1
sname <- sub(" ", "", sname)
sname <- gsub("[+]", "p", sname)
sname.gs <- factor(sname)

#### Expression Boxplot ####
ex <- exprs(gset)

#pdf("result/boxplot2.pdf")
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)
palette(c("blue", "red"))
par(mar=c(7,4,2,1))
boxplot(ex[,ord], boxwex=0.6, notch=T, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
#dev.off()


#### Normalize, if required ####
ex <- normalizeQuantiles(ex)
exprs(gset) <- ex

#### Expression Boxplot after Normalize ####
ex <- exprs(gset)
#pdf("result/boxplot2.pdf")
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)
palette(c("blue", "red"))
par(mar=c(7,4,2,1))
boxplot(ex[,ord], boxwex=0.6, notch=T, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
#dev.off()

#### Dim. reduction and Scale expressions ####
pc <- prcomp(ex)
plot(pc)
plot(pc$x[,1:2])
pcr <- data.frame(pc$r[,1:3], Group = gset$group)
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3)

ex.scale <- t(scale(t(ex), scale = F))
pc_scale <- prcomp(ex.scale)
plot(pc_scale)
plot(pc_scale$x[,1:2])
pcr_scale <- data.frame(pc_scale$r[,1:3], Group = gset$group)
ggplot(pcr_scale, aes(PC1, PC2, color=Group)) + geom_point(size=3)


#### Heatmap and Correlation of samples ####
ex.scale.cor <- cor(ex.scale)
#pdf("result/pheatmap.pdf", width = 10, height = 10)
pheatmap(ex.scale.cor, 
         labels_row = gset$source_name_ch1, 
         labels_col = gset$source_name_ch1, 
         color = bluered(255), border_color = NA)
#dev.off()


#### Finding sorted correlation with AML Patient ####
df <- data.frame(ex.scale.cor)
a <- t(df)
colnames(a) <- gset$source_name_ch1
a <- t(a)
colnames(a) <- gset$source_name_ch1
cols <- gset$source_name_ch1
cols <- cols[cols != "AML Patient"]
b <- subset(a, select=cols)
b <- t(b)
cols <- gset$source_name_ch1
cols <- cols[cols == "AML Patient"]
b <- subset(b, select=cols)
b <- t(b)
bB <- b
maxCor <- data.frame(Cells = unique(colnames(bB)[max.col(bB,ties.method="first")]), 
                    CorWithAML = unique(rowMax(bB)))
for (x in 1:4) {
  cols <- colnames(bB)
  cols <- cols[cols != maxCor[[1]][x]]
  bB <- subset(bB, select=cols)
  maxCor <- rbind(maxCor, 
                  c(gene = unique(colnames(bB)[max.col(bB,ties.method="first")]), 
                            cor = unique(rowMax(bB))))
}


#### Diff Expression Analysis ####
##### based on all test-normal #####
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)
## fitting linear model to data
fit <- lmFit(gset, design)

cont.matrix <- makeContrasts(test-normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
## adjust by: fault-discovery-rate or Benjamini-hochberg
## sort by adj.P.Val
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID","adj.P.Val","logFC", "B"))
write.table(tT, "result/dea/dea_test-normal_B.txt", row.names=F, sep="\t", quote=F)
### Top Gene Expression mining 
aml.up <- subset(tT, logFC > 1 & adj.P.Val > 0.05)
aml.up.genes <- unique(as.character(strsplit2(unique(aml.up$Gene.symbol), "///")))
write.table(aml.up.genes, "result/dea/dea_test-normal_Up.txt", 
            quote=F, row.names=F, col.names=F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val > 0.05)
aml.down.genes <- unique(as.character(strsplit2(unique(aml.down$Gene.symbol), "///")))
write.table(aml.down.genes, "result/dea/dea_test-normal_Down.txt", 
            quote=F, row.names=F, col.names=F)

##### based on top correlated cells #####
design <- model.matrix(~source_name_ch1 + 0, gset)
sfl <- factor(sname.gs)
colnames(design) <- levels(sfl)
## fitting linear model to data
fit <- lmFit(gset, design)
###### AMLPatient-Monocytes ######
cont.matrix <- makeContrasts(AMLPatient-Monocytes, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
## adjust by: fault-discovery-rate or Benjamini-hochberg
## sort by adj.P.Val
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID","adj.P.Val","logFC", "B"))
write.table(tT, "result/dea/dea_AMLPatient-Monocytes_B.txt", 
            row.names=F, sep="\t", quote=F)
### Top Gene Expression mining 
aml.up <- subset(tT, logFC > 1 & adj.P.Val > 0.05)
aml.up.genes <- unique(as.character(strsplit2(unique(aml.up$Gene.symbol), "///")))
write.table(aml.up.genes, "result/dea/dea_AMLPatient-Monocytes_Up.txt", 
            quote=F, row.names=F, col.names=F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val > 0.05)
aml.down.genes <- unique(as.character(strsplit2(unique(aml.down$Gene.symbol), "///")))
write.table(aml.down.genes, "result/dea/dea_AMLPatient-Monocytes_Down.txt", 
            quote=F, row.names=F, col.names=F)
###### AMLPatient-CD34+HSPC ######
cont.matrix <- makeContrasts(AMLPatient-CD34pHSPC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
## adjust by: fault-discovery-rate or Benjamini-hochberg
## sort by adj.P.Val
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID","adj.P.Val","logFC", "B"))
write.table(tT, "result/dea/dea_AMLPatient-CD34pHSPC_B.txt", 
            row.names=F, sep="\t", quote=F)
### Top Gene Expression mining 
aml.up <- subset(tT, logFC > 1 & adj.P.Val > 0.05)
aml.up.genes <- unique(as.character(strsplit2(unique(aml.up$Gene.symbol), "///")))
write.table(aml.up.genes, "result/dea/dea_AMLPatient-CD34pHSPC_Up.txt", 
            quote=F, row.names=F, col.names=F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val > 0.05)
aml.down.genes <- unique(as.character(strsplit2(unique(aml.down$Gene.symbol), "///")))
write.table(aml.down.genes, "result/dea/dea_AMLPatient-CD34pHSPC_Down.txt", 
            quote=F, row.names=F, col.names=F)





