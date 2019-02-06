library(DESeq2)

a = read.table (file = "superEnhancer_counts.txt",header = T,sep = "\t",row.names = 1)

coldata = read.table(file = "superEnhancer_names.txt",header = T,sep = "\t")[,]
dso = DESeqDataSetFromMatrix(countData = a,colData = coldata,design =  ~ Sample)
dsoA<-DESeq(dso )
write.table((results(dsoA)),file="results.txt",sep = "\t",quote = F)

pdf(file = "PCA.pdf")
plotPCA(rlog(dsoA),intgroup=c("pair"),ntop = 1000)

dev.off()

