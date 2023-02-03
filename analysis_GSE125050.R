library(edgeR)
library(limma)
library(DESeq)
library(DESeq2)
library(vidger)
library(ggplot2)
library(MASS)

#3d pca kodu

plotPCA3D <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  message("Generating plotly plot")
  p <- plotly::plot_ly(data = d,
                       x = ~PC1,
                       y = ~PC2,
                       z = ~PC3,
                       color = group,
                       mode = "markers",
                       type = "scatter3d")
  return(p)
}


#veriyi indirme ve çekme


url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE125050&format=file"
utils::download.file(url, destfile="GSE125050_RAW.tar", mode="wb") #for download data from GEO
utils::untar("GSE125050_RAW.tar", exdir = ".")#for untar
files <- c("GSM3561843_sample1.tsv","GSM3561844_sample2.tsv","GSM3561845_sample3.tsv","GSM3561846_sample4.tsv","GSM3561847_sample5.tsv","GSM3561848_sample6.tsv","GSM3561849_sample7.tsv","GSM3561850_sample8.tsv","GSM3561851_sample9.tsv","GSM3561852_sample10.tsv","GSM3561853_sample11.tsv","GSM3561854_sample12.tsv","GSM3561855_sample13.tsv","GSM3561856_sample14.tsv","GSM3561857_sample15.tsv","GSM3561858_sample16.tsv","GSM3561859_sample17.tsv","GSM3561860_sample18.tsv","GSM3561861_sample19.tsv","GSM3561862_sample20.tsv","GSM3561863_sample21.tsv","GSM3561864_sample22.tsv","GSM3561865_sample23.tsv","GSM3561866_sample24.tsv","GSM3561867_sample25.tsv","GSM3561868_sample26.tsv","GSM3561869_sample27.tsv","GSM3561870_sample28.tsv","GSM3561871_sample29.tsv","GSM3561872_sample30.tsv","GSM3561873_sample31.tsv","GSM3561874_sample32.tsv","GSM3561875_sample33.tsv","GSM3561876_sample34.tsv","GSM3561877_sample35.tsv","GSM3561878_sample36.tsv","GSM3561879_sample37.tsv","GSM3561880_sample38.tsv","GSM3561881_sample39.tsv","GSM3561882_sample40.tsv","GSM3561883_sample41.tsv","GSM3561884_sample42.tsv","GSM3561885_sample43.tsv","GSM3561886_sample44.tsv","GSM3561887_sample45.tsv","GSM3561888_sample46.tsv","GSM3561889_sample47.tsv","GSM3561890_sample48.tsv","GSM3561891_sample49.tsv","GSM3561892_sample50.tsv","GSM3561893_sample51.tsv","GSM3561894_sample52.tsv","GSM3561895_sample53.tsv","GSM3561896_sample54.tsv","GSM3561897_sample55.tsv","GSM3561898_sample56.tsv","GSM3561899_sample57.tsv","GSM3561900_sample58.tsv","GSM3561901_sample59.tsv","GSM3561902_sample60.tsv","GSM3561903_sample61.tsv","GSM3561904_sample62.tsv","GSM3561905_sample63.tsv","GSM3561906_sample64.tsv","GSM3561907_sample65.tsv","GSM3561908_sample66.tsv","GSM3561909_sample67.tsv","GSM3561910_sample68.tsv","GSM3561911_sample69.tsv","GSM3561912_sample70.tsv","GSM3561913_sample71.tsv","GSM3561914_sample72.tsv","GSM3561915_sample73.tsv","GSM3561916_sample74.tsv","GSM3561917_sample75.tsv","GSM3561918_sample76.tsv","GSM3561919_sample77.tsv","GSM3561920_sample78.tsv","GSM3561921_sample79.tsv","GSM3561922_sample80.tsv","GSM3561923_sample81.tsv","GSM3561924_sample82.tsv","GSM3561925_sample83.tsv","GSM3561926_sample84.tsv","GSM3561927_sample85.tsv","GSM3561928_sample86.tsv","GSM3561929_sample87.tsv","GSM3561930_sample88.tsv","GSM3561931_sample89.tsv","GSM3561932_sample90.tsv","GSM3561933_sample91.tsv","GSM3561934_sample92.tsv","GSM3561935_sample93.tsv","GSM3561936_sample94.tsv","GSM3561937_sample95.tsv","GSM3561938_sample96.tsv","GSM3561939_sample97.tsv","GSM3561940_sample98.tsv","GSM3561941_sample99.tsv","GSM3561942_sample100.tsv","GSM3561943_sample101.tsv","GSM3561944_sample102.tsv","GSM3561945_sample103.tsv","GSM3561946_sample104.tsv","GSM3561947_sample105.tsv","GSM3561948_sample106.tsv","GSM3561949_sample107.tsv","GSM3561950_sample108.tsv","GSM3561951_sample109.tsv","GSM3561952_sample110.tsv","GSM3561953_sample111.tsv","GSM3561954_sample112.tsv","GSM3561955_sample113.tsv")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

x <- readDGE(files, columns=c(1,3),comment.char ="#") #extract data

samplenames <- substring(colnames(x), 12, nchar(colnames(x)))

colnames(x) <- samplenames

group_all <- as.factor(c("Control myeloid","Control endothelial","Control neuron", "AD neuron", "Control myeloid", "AD myeloid","Control myeloid","Control endothelial", "AD endothelial","Control neuron","AD neuron","AD neuron","Control endothelial","Control astrocyte","Control astrocyte","Control neuron","Control neuron","Control endothelial","AD neuron","Control neuron","AD myeloid","AD neuron","AD endothelia","Control neuron","AD neuron","Control endothelial","Control neuron","AD astrocyte","Control neuron","AD astrocyte","Control neuron","Control endothelial","Control endothelial","AD endothelial","AD myeloid","Control endothelial","Control astrocyte","AD astrocyte","AD neuron","Control myeloid","AD neuron","AD myeloid","Control astrocyte","AD myeloid","Control neuron","AD myeloid","AD neuron","AD endothelial","AD endothelial","AD neuron","Control neuron","Control astrocyte","Control endothelial","Control neuron","Control endothelial","Control astrocyte","AD astrocyte","AD neuron","AD myeloid","Control neuron","AD neuron","Control neuron","Control astrocyte","Control neuron","AD neuron","AD neuron","AD endothelial","Control neuron","Control myeloid","Control endothelial","Control myeloid","AD myeloid","AD neuron","AD neuron","Control myeloid","Control myeloid","Control astrocyte","AD endothelial","AD neuron","Control neuron","AD astrocyte","Control myeloid","Control myeloid","Control astrocyte","Control myeloid","Control endothelial","Control myeloid","Control endothelial","Control myeloid","Control endothelial","AD neuron","Control neuron","AD neuron","Control myeloid","AD neuron","AD myeloid","Control myeloid","Control astrocyte","Control endothelial","AD myeloid","AD endothelial","Control astrocyte","Control neuron","AD endothelial","Control neuron","Control endothelial","Control astrocyte","Control neuron","AD astrocyte","AD neuron","AD endothelial","Control endothelial","AD astrocyte"))

x$samples$group_all <- group_all

dds_counts_all_data <- DESeqDataSetFromMatrix(countData = x$counts, colData = x$samples, design = ~group_all) #transform to deseq

dds_all_data <- DESeq(dds_counts_all_data)

norm_all_data<-vst(dds_all_data) #normalization for data

DESeq2::plotPCA(norm_all_data, intgroup = c("group"))#plot pca graph

#neuron and astrocyte together

neu_ast<- c("GSM3561845_sample3.tsv","GSM3561846_sample4.tsv","GSM3561852_sample10.tsv","GSM3561853_sample11.tsv","GSM3561854_sample12.tsv","GSM3561858_sample16.tsv","GSM3561859_sample17.tsv","GSM3561861_sample19.tsv","GSM3561862_sample20.tsv","GSM3561864_sample22.tsv","GSM3561866_sample24.tsv","GSM3561867_sample25.tsv","GSM3561869_sample27.tsv","GSM3561871_sample29.tsv","GSM3561873_sample31.tsv","GSM3561881_sample39.tsv","GSM3561883_sample41.tsv","GSM3561887_sample45.tsv","GSM3561889_sample47.tsv","GSM3561892_sample50.tsv","GSM3561893_sample51.tsv","GSM3561896_sample54.tsv","GSM3561900_sample58.tsv","GSM3561902_sample60.tsv","GSM3561903_sample61.tsv","GSM3561904_sample62.tsv","GSM3561906_sample64.tsv","GSM3561907_sample65.tsv","GSM3561908_sample66.tsv","GSM3561910_sample68.tsv","GSM3561915_sample73.tsv","GSM3561916_sample74.tsv","GSM3561921_sample79.tsv","GSM3561922_sample80.tsv","GSM3561933_sample91.tsv","GSM3561934_sample92.tsv","GSM3561935_sample93.tsv","GSM3561937_sample95.tsv","GSM3561945_sample103.tsv","GSM3561947_sample105.tsv","GSM3561950_sample108.tsv","GSM3561952_sample110.tsv","GSM3561856_sample14.tsv","GSM3561857_sample15.tsv","GSM3561870_sample28.tsv","GSM3561872_sample30.tsv","GSM3561879_sample37.tsv","GSM3561880_sample38.tsv","GSM3561885_sample43.tsv","GSM3561894_sample52.tsv","GSM3561898_sample56.tsv","GSM3561899_sample57.tsv","GSM3561905_sample63.tsv","GSM3561919_sample77.tsv","GSM3561923_sample81.tsv","GSM3561926_sample84.tsv","GSM3561940_sample98.tsv","GSM3561944_sample102.tsv","GSM3561949_sample107.tsv","GSM3561951_sample109.tsv","GSM3561955_sample113.tsv")

x_neu_ast<- readDGE(neu_ast, columns=c(1,3),comment.char ="#")

samplenames <- substring(colnames(x_neu_ast), 12, nchar(colnames(x_neu_ast)))

colnames(x_neu_ast) <- samplenames

group_neu_ast <- as.factor(c("Controlneuron","ADneuron","Controlneuron","ADneuron","ADneuron","Controlneuron","Controlneuron","ADneuron","Controlneuron","ADneuron","Controlneuron","ADneuron","Controlneuron","Controlneuron","Controlneuron","ADneuron","ADneuron","Controlneuron","ADneuron","ADneuron","Controlneuron","Controlneuron","ADneuron","Controlneuron","ADneuron","Controlneuron","Controlneuron","ADneuron","ADneuron","Controlneuron","ADneuron","ADneuron","ADneuron","Controlneuron","ADneuron","Controlneuron","ADneuron","ADneuron","Controlneuron","Controlneuron","Controlneuron","ADneuron","Controlastrocyte","Controlastrocyte","ADastrocyte","ADastrocyte","Controlastrocyte","ADastrocyte","Controlastrocyte","Controlastrocyte","Controlastrocyte","ADastrocyte","Controlastrocyte","Controlastrocyte","ADastrocyte","Controlastrocyte","Controlastrocyte","Controlastrocyte","Controlastrocyte","ADastrocyte","ADastrocyte"))

x_neu_ast$samples$group_neu_ast<-group_neu_ast

dds_counts_neu_ast <- DESeqDataSetFromMatrix(countData = x_neu_ast$counts, colData = x_neu_ast$samples, design = ~group_neu_ast) #transform to deseq

dds_neu_ast <- DESeq(dds_counts_neu_ast)

norm_neu_ast<-vst(dds_neu_ast) #normalization for data

DESeq2::plotPCA(norm_neu_ast, intgroup = c("group_neu_ast"))#plot pca graph

#neuron cells

neuron<-c("GSM3561845_sample3.tsv","GSM3561852_sample10.tsv","GSM3561858_sample16.tsv","GSM3561859_sample17.tsv","GSM3561862_sample20.tsv","GSM3561866_sample24.tsv","GSM3561869_sample27.tsv","GSM3561871_sample29.tsv","GSM3561873_sample31.tsv","GSM3561887_sample45.tsv","GSM3561893_sample51.tsv","GSM3561896_sample54.tsv","GSM3561902_sample60.tsv","GSM3561904_sample62.tsv","GSM3561906_sample64.tsv","GSM3561910_sample68.tsv","GSM3561922_sample80.tsv","GSM3561934_sample92.tsv","GSM3561945_sample103.tsv","GSM3561947_sample105.tsv","GSM3561950_sample108.tsv","GSM3561846_sample4.tsv","GSM3561853_sample11.tsv","GSM3561854_sample12.tsv","GSM3561861_sample19.tsv","GSM3561864_sample22.tsv","GSM3561867_sample25.tsv","GSM3561881_sample39.tsv","GSM3561883_sample41.tsv","GSM3561889_sample47.tsv","GSM3561892_sample50.tsv","GSM3561900_sample58.tsv","GSM3561903_sample61.tsv","GSM3561907_sample65.tsv","GSM3561908_sample66.tsv","GSM3561915_sample73.tsv","GSM3561916_sample74.tsv","GSM3561921_sample79.tsv","GSM3561933_sample91.tsv","GSM3561935_sample93.tsv","GSM3561937_sample95.tsv","GSM3561952_sample110.tsv")

group_neuron<-as.factor(c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD"))

x_neuron<-readDGE(neuron, columns=c(1,3),comment.char ="#")

x_neuron$samples$group_neuron <- as.factor(group_neuron)

colnames(x_neuron) <- substring(colnames(x_neuron), 12, nchar(colnames(x_neuron)))

dds_counts_neuron <- DESeqDataSetFromMatrix(countData = x_neuron$counts, colData = x_neuron$samples, design = ~group_neuron) #transform to deseq

dds_neuron <- DESeq(dds_counts_neuron)

norm_neuron<-vst(dds_neuron) #normalization for data

DESeq2::plotPCA(norm_neuron, intgroup = c("group_neuron"))#plot pca graph

plotPCA3D(norm_neuron,intgroup = c("group_neuron"))

vsBoxPlot(
  data = dds_neuron, d.factor = 'group_neuron', type = 'deseq', 
  title = TRUE, legend = TRUE, grid = TRUE)

#astrocyte

astrocyte<-c("GSM3561856_sample14.tsv","GSM3561857_sample15.tsv","GSM3561879_sample37.tsv","GSM3561885_sample43.tsv","GSM3561894_sample52.tsv","GSM3561898_sample56.tsv","GSM3561905_sample63.tsv","GSM3561919_sample77.tsv","GSM3561926_sample84.tsv","GSM3561940_sample98.tsv","GSM3561944_sample102.tsv","GSM3561949_sample107.tsv","GSM3561870_sample28.tsv","GSM3561872_sample30.tsv","GSM3561880_sample38.tsv","GSM3561899_sample57.tsv","GSM3561923_sample81.tsv","GSM3561951_sample109.tsv","GSM3561955_sample113.tsv")

group_astrocyte<-as.factor(c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","AD","AD","AD","AD","AD","AD","AD"))

x_astrocyte<-readDGE(astrocyte, columns=c(1,3),comment.char ="#")

x_astrocyte$samples$group_astrocyte <- as.factor(group_astrocyte)

colnames(x_astrocyte) <- substring(colnames(x_astrocyte), 12, nchar(colnames(x_astrocyte)))

dds_counts_astrocyte <- DESeqDataSetFromMatrix(countData = x_astrocyte$counts, colData = x_astrocyte$samples, design = ~group_astrocyte)

dds_astrocyte <- DESeq(dds_counts_astrocyte)

norm_astrocyte<-vst(dds_astrocyte) #normalization for data

#plotPCA(norm_astrocyte, intgroup = c("group_astrocyte"))#plot pca graph

DESeq2::plotPCA(norm_astrocyte, intgroup = c("group_astrocyte"))+ theme_bw() + geom_text(aes(label = colnames(norm_astrocyte)))

vsBoxPlot(
  data = dds_astrocyte, d.factor = 'group_astrocyte', type = 'deseq', 
  title = TRUE, legend = TRUE, grid = TRUE)




#astrocyte

###outliersız pca ,outliier olan örnekleri çıkartıp analize öyle devam ettim.

astrocyte<-c("GSM3561856_sample14.tsv","GSM3561857_sample15.tsv","GSM3561879_sample37.tsv","GSM3561885_sample43.tsv","GSM3561894_sample52.tsv","GSM3561898_sample56.tsv","GSM3561905_sample63.tsv","GSM3561919_sample77.tsv","GSM3561926_sample84.tsv","GSM3561944_sample102.tsv","GSM3561949_sample107.tsv","GSM3561870_sample28.tsv","GSM3561872_sample30.tsv","GSM3561880_sample38.tsv","GSM3561899_sample57.tsv","GSM3561923_sample81.tsv","GSM3561955_sample113.tsv")

group_astrocyte<-as.factor(c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","AD","AD","AD","AD","AD","AD"))

x_astrocyte<-readDGE(astrocyte, columns=c(1,3),comment.char ="#")

x_astrocyte$samples$group_astrocyte <- as.factor(group_astrocyte)

colnames(x_astrocyte) <- substring(colnames(x_astrocyte), 12, nchar(colnames(x_astrocyte)))

dds_counts_astrocyte <- DESeqDataSetFromMatrix(countData = x_astrocyte$counts, colData = x_astrocyte$samples, design = ~group_astrocyte)

dds_astrocyte <- DESeq(dds_counts_astrocyte)



#for neuron results

#forpvalue

res_neuron01<-results(dds_neuron)

ressig_neuronpvalue01 <- subset(res_neuron01, res_neuron01$pvalue < 0.01 )

ressig_neuronpvalue01<-subset(ressig_neuronpvalue01,ressig_neuronpvalue01$log2FoldChange< -1|ressig_neuronpvalue01$log2FoldChange>1)

write.csv(ressig_neuronpvalue01,"ressig_neuronpvalue01.txt")

#forpvalue

res_neuron05<-results(dds_neuron)

ressig_neuronpvalue05 <- subset(res_neuron05, res_neuron05$pvalue < 0.05 )

ressig_neuronpvalue05<-subset(ressig_neuronpvalue05,ressig_neuronpvalue05$log2FoldChange< -1|ressig_neuronpvalue05$log2FoldChange>1)

write.csv(ressig_neuronpvalue05,"ressig_neuronpvalue05.txt")

#forpadjBH

res_neuron<-results(dds_neuron)

ressig_neuronpadj <- subset(res_neuron, res_neuron$padj < 0.01 )

ressig_neuronpadj<-subset(ressig_neuronpadj,ressig_neuronpadj$log2FoldChange< -1 | ressig_neuronpadj$log2FoldChange>1 )

write.csv(ressig_neuronpadj,"ressig_neuronpadjBH.txt")

#forpadjBonferroni

res_neuron<-results(dds_neuron)

newcutoff<-0.01/length(rownames(x_neuron$counts))

ressig_neuronpadjBonf <- subset(res_neuron, res_neuron$pvalue < newcutoff )

ressig_neuronpadjBonf<-subset(ressig_neuronpadjBonf,ressig_neuronpadjBonf$log2FoldChange< -1 | ressig_neuronpadjBonf$log2FoldChange>1 )

write.csv(ressig_neuronpadjBonf,"ressig_neuronpadjBonf.txt")


#for astrocyte results

#forpval <0.01

res_astrocyte01<-results(dds_astrocyte)

ressig_astrocytepval01 <- subset(res_astrocyte01, res_astrocyte01$pvalue < 0.01 )

ressig_astrocytepval01<-subset(ressig_astrocytepval01,ressig_astrocytepval01$log2FoldChange< -1|ressig_astrocytepval01$log2FoldChange>1)

write.csv(ressig_astrocytepval01,"ressig_astrocytepval01.txt")

#forpval <0.05

res_astrocyte05<-results(dds_astrocyte)

ressig_astrocytepval05 <- subset(res_astrocyte05, res_astrocyte05$pvalue < 0.05 )

ressig_astrocytepval05<-subset(ressig_astrocytepval05,ressig_astrocytepval05$log2FoldChange < -1|ressig_astrocytepval05$log2FoldChange >1)

write.csv(ressig_astrocytepval05,"ressig_astrocytepval05.txt")

#forpadjBH

res_astrocyte<-results(dds_astrocyte)

ressig_astrocytepadj <- subset(res_astrocyte, res_astrocyte$pvadj < 0.01 )

ressig_astrocytepadj<-subset(ressig_astrocytepadj,ressig_astrocytepadj$log2FoldChange< -1|ressig_astrocytepadj$log2FoldChange>1)

write.csv(ressig_astrocytepadj,"ressig_astrocytepadjBH.txt")

#forpadjBonferroni

res_astrocyte<-results(dds_astrocyte)

newcutoff<-0.01/length(rownames(x_astrocyte$counts))

ressig_astrocytepadjBonf <- subset(res_astrocyte, res_astrocyte$pvalue < newcutoff )

ressig_astrocytepadjBonf<-subset(ressig_astrocytepadjBonf,ressig_astrocytepadjBonf$log2FoldChange< -1|ressig_astrocytepadjBonf$log2FoldChange>1)

write.csv(ressig_astrocytepadjBonf,"ressig_astrocytepadjBonf.txt")

write.csv(log2(x_astrocyte$counts),"deneme.csv")


#lda neuron

neuron<-c("GSM3561845_sample3.tsv","GSM3561852_sample10.tsv","GSM3561858_sample16.tsv","GSM3561859_sample17.tsv","GSM3561862_sample20.tsv","GSM3561866_sample24.tsv","GSM3561869_sample27.tsv","GSM3561871_sample29.tsv","GSM3561873_sample31.tsv","GSM3561887_sample45.tsv","GSM3561893_sample51.tsv","GSM3561896_sample54.tsv","GSM3561902_sample60.tsv","GSM3561904_sample62.tsv","GSM3561906_sample64.tsv","GSM3561910_sample68.tsv","GSM3561922_sample80.tsv","GSM3561934_sample92.tsv","GSM3561945_sample103.tsv","GSM3561947_sample105.tsv","GSM3561950_sample108.tsv","GSM3561846_sample4.tsv","GSM3561853_sample11.tsv","GSM3561854_sample12.tsv","GSM3561861_sample19.tsv","GSM3561864_sample22.tsv","GSM3561867_sample25.tsv","GSM3561881_sample39.tsv","GSM3561883_sample41.tsv","GSM3561889_sample47.tsv","GSM3561892_sample50.tsv","GSM3561900_sample58.tsv","GSM3561903_sample61.tsv","GSM3561907_sample65.tsv","GSM3561908_sample66.tsv","GSM3561915_sample73.tsv","GSM3561916_sample74.tsv","GSM3561921_sample79.tsv","GSM3561933_sample91.tsv","GSM3561935_sample93.tsv","GSM3561937_sample95.tsv","GSM3561952_sample110.tsv")

group_neuron<-as.factor(c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD"))

x_neuron<-readDGE(neuron, columns=c(1,3),comment.char ="#")

data_new1<-x_neuron$counts


data_new1 = data_new1[as.logical(rowSums(data_new1 != 0)), ]

data_new1<-matrix(data_new1,ncol = 42)

t<-t(data_new1[1:29668,])

pe1<-as.data.frame(t(data_new1))

v1<-lda(pe1[,1:29668],group_neuron)

plot(v1)



#lda astrocyte


astrocyte<-c("GSM3561856_sample14.tsv","GSM3561857_sample15.tsv","GSM3561879_sample37.tsv","GSM3561885_sample43.tsv","GSM3561894_sample52.tsv","GSM3561898_sample56.tsv","GSM3561905_sample63.tsv","GSM3561919_sample77.tsv","GSM3561926_sample84.tsv","GSM3561944_sample102.tsv","GSM3561949_sample107.tsv","GSM3561870_sample28.tsv","GSM3561872_sample30.tsv","GSM3561880_sample38.tsv","GSM3561899_sample57.tsv","GSM3561923_sample81.tsv","GSM3561955_sample113.tsv")

group_astrocyte<-as.factor(c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","AD","AD","AD","AD","AD","AD"))

x_astrocyte<-readDGE(astrocyte, columns=c(1,3),comment.char ="#")

data_new<-x_astrocyte$counts

data_new = data_new[as.logical(rowSums(data_new != 0)), ]

data_new<-matrix(data_new,ncol = 17)

t<-t(data_new[1:150,])

pe<-as.data.frame(t(data_new))

v1<-lda(pe[,1:150],group_astrocyte)

plot(v1)


#Keypathwayminer için çekme
#neuron
res<-ressig_neuronpvalue05

x<-rownames(res)
x<-matrix(x)

y<-rownames(x_neuron$counts)
y<-matrix(y)
bse<-y
r<-matrix(,ncol = 2,nrow = 30727)
r[,1]=y[,1]
r[,2]=rep(c(0))

for (i in 1:2079){
  for (j in 1:30727) {
    if(x[i]==y[j,1]){
      r[j,2]<-1
      
    }
    
  }
}

rownames(r)<-r[,1]
r<-r[,-1]

write.csv(r,"neuronke.txt")

#astrocyte

res<-ressig_astrocytepval05

x<-rownames(res)
x<-matrix(x)

y<-rownames(x_astrocyte$counts)
y<-matrix(y)
bse<-y
r<-matrix(,ncol = 2,nrow = 30727)
r[,1]=y[,1]
r[,2]=rep(c(0))

for (i in 1:338){
  for (j in 1:30727) {
    if(x[i]==y[j,1]){
      r[j,2]<-1
      
    }
    
  }
}

rownames(r)<-r[,1]
r<-r[,-1]

write.csv(r,"astrocyte_kpm.txt")

#entrezID to genesymbol

library(org.Hs.eg.db)

t<-c(mapIds(org.Hs.eg.db, y, 'SYMBOL','ENTREZID'))

write.csv(t,"genesymbol.csv")

#astrocyte genesymbol

rnameastrocyte<-rownames(ressig_astrocytepval05)

t<-c(mapIds(org.Hs.eg.db, rnameastrocyte, 'SYMBOL','ENTREZID'))

write.csv(t,"astrocytegenesymbol.csv")
