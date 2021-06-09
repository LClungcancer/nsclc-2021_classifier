library(ConsensusClusterPlus)
library(pheatmap)
library(survival)
library(survminer)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

##illumina 
indata<-read.table("illumina_cibersort.txt",sep = "\t",header = T,row.names = 1)
dim(indata)

subtype <- ConsensusClusterPlus(d = as.matrix(indata),
                                maxK = 6, # 原文参数
                                pItem = 0.8, # 默认参数
                                pFeature = 1, # 默认参数
                                reps = 1000, # 原文参数
                                clusterAlg = "km", # 原文参数
                                innerLinkage = "ward.D", # 原文参数
                                finalLinkage = "ward.D", # 原文参数
                                distance = "euclidean", # 原文参数
                                seed = 123456,
                                plot = "png", #或png
                                writeTable = TRUE,
                                title = "km_ConsensusCluster") 
geneclust <- subtype[[3]]
table(geneclust$consensusClass)
rt<-as.data.frame(geneclust$consensusClass)
colnames(rt)[1]<-"cluster"

##for Submap input
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

skcm.immunotherapy.logNC <- read.table("nsclc.txt",sep = "\t",
                                       row.names = 1,header = T,check.names = F,stringsAsFactors = F)
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.logNC=skcm.immunotherapy.logNC[which(apply(skcm.immunotherapy.logNC,1,function(x){return(sum(x>1))})>ncol(skcm.immunotherapy.logNC)*0.9),]

skcm.immunotherapy.info <- read.table("nsclc_cluster.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
skcm.immunotherapy.info$sample<-rownames(skcm.immunotherapy.info)
skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$cluster),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3),times=as.character(table(skcm.immunotherapy.info$cluster))) 

tmp <-read.csv("illumina.csv",header = T,row.names = 1)
tmp=tmp[which(apply(tmp,1,function(x){return(sum(x>1))})>ncol(tmp)*0.9),]
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC))

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

samples.C1 <- rownames(rt[which(rt$cluster == "1"),])
samples.C2 <- rownames(rt[which(rt$cluster == "2"),])
samples.C3 <- rownames(rt[which(rt$cluster == "3"),])

sam_info <- data.frame("cluster"=c(samples.C1,samples.C2,samples.C3),row.names = c(samples.C1,samples.C2,samples.C3))
sam_info$rank <- rep(c(1,2,3),times=c(length(samples.C1),length(samples.C2),length(samples.C3))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

gct_file <- "illumina.for.SubMap.gct"
cls_file <- "illumina.for.SubMap.cls"

in_gct <- tmp[GENELIST,rownames(sam_info)] # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

###gct cls files in GenePattern submap module
###for the output



###Affy
indata<-read.table("affy_cibersort.txt",sep = "\t",header = T,row.names = 1)

subtype <- ConsensusClusterPlus(d = as.matrix(indata),
                                maxK = 6, # 原文参数
                                pItem = 0.8, # 默认参数
                                pFeature = 1, # 默认参数
                                reps = 1000, # 原文参数
                                clusterAlg = "km", # 原文参数
                                innerLinkage = "ward.D", # 原文参数
                                finalLinkage = "ward.D", # 原文参数
                                distance = "euclidean", # 原文参数
                                seed = 123456,
                                plot = "png", #或png
                                writeTable = TRUE,
                                title = "km_ConsensusCluster") 

geneclust <- subtype[[3]]
table(geneclust$consensusClass)

rt<-as.data.frame(geneclust$consensusClass)
colnames(rt)[1]<-"cluster"
rt$sample<-rownames(rt)

##
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

skcm.immunotherapy.logNC <- read.table("nsclc.txt",sep = "\t",
                                       row.names = 1,header = T,check.names = F,stringsAsFactors = F)
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.logNC=skcm.immunotherapy.logNC[which(apply(skcm.immunotherapy.logNC,1,function(x){return(sum(x>1))})>ncol(skcm.immunotherapy.logNC)*0.9),]

skcm.immunotherapy.info <- read.table("nsclc_cluster.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
skcm.immunotherapy.info$sample<-rownames(skcm.immunotherapy.info)
skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$cluster),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3),times=as.character(table(skcm.immunotherapy.info$cluster))) 


tmp <-read.csv("affy_output_combined_expr.csv",header = T,row.names = 1)
tmp=tmp[which(apply(tmp,1,function(x){return(sum(x>1))})>ncol(tmp)*0.9),]
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

samples.C1 <- rownames(rt[which(rt$cluster == "1"),])
samples.C2 <- rownames(rt[which(rt$cluster == "2"),])
samples.C3 <- rownames(rt[which(rt$cluster == "3"),])

sam_info <- data.frame("cluster"=c(samples.C1,samples.C2,samples.C3),row.names = c(samples.C1,samples.C2,samples.C3))
sam_info$rank <- rep(c(1,2,3),times=c(length(samples.C1),length(samples.C2),length(samples.C3))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

gct_file <- "affy.for.SubMap.gct"
cls_file <- "affy.for.SubMap.cls"

in_gct <- tmp[GENELIST,rownames(sam_info)] # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

###Result for GenePattern submap





###Calculation for immune-parameters
##Illumina
library(GSVA)
library(estimate)
library(ggplot2)
library(grid)
library(gridExtra)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 

immunity <- read.csv("41.csv", header = T)

immunity <- immunity %>% 
  split(., .$CellType) %>% 
  lapply(., function(x)(x$genesymbol))
immunity <- lapply(immunity, unique)

expr<-read.csv("illumina.csv",header = T,row.names = 1)
tcga_gsva <- as.data.frame(t(gsva(as.matrix(expr), immunity, method = "ssgsea")))
write.csv(tcga_gsva,"illumina_41gsva.csv",quote = F)


filterCommonGenes(input.f="illumina_uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")


scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="illumina_scores.txt",sep="\t",quote=F,col.names=F)


ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}

gene_expression<-read.table("illumina_uniq.symbol.txt",row.names=1,header=TRUE, sep="\t", dec = ".",check.names=FALSE)
sample_names<-names(gene_expression)

## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
# For different 
IPSG<-read.table("IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes<-as.vector(unique(IPSG$NAME))

IPS<-NULL
MHC<-NULL
CP<-NULL
EC<-NULL
SC<-NULL
AZ<-NULL

GVEC<-row.names(gene_expression)

VEC<-as.vector(IPSG$GENE)

ind<-which(is.na(match(VEC,GVEC)))

MISSING_GENES<-VEC[ind]
dat<-IPSG[ind,]
if (length(MISSING_GENES)>0) {
  cat("differently named or missing genes: ",MISSING_GENES,"\n")
}
for (x in 1:length(ind)) {
  print(IPSG[ind,])
}

for (i in 1:length(sample_names)) {	
  GE<-gene_expression[[i]]
  mGE<-mean(GE)
  sGE<-sd(GE)
  Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1<-IPSG$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k<-k+1
  }
  WG<-MIG*WEIGHT
  WG[is.na(WG)]<-0   
  MHC[i]<-mean(WG[1:10])
  CP[i]<-mean(WG[11:20])
  EC[i]<-mean(WG[21:24])
  SC[i]<-mean(WG[25:26])
  AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i]<-ipsmap(AZ[i])
}

DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
write.table(DF,file="illumina_IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")

###Combined the above results as illumina_train.txt

##Affy
immunity <- read.csv("41.csv", header = T)

immunity <- immunity %>% 
  split(., .$CellType) %>% 
  lapply(., function(x)(x$genesymbol))
immunity <- lapply(immunity, unique)

expr<-read.csv("affy_output_combined_expr.csv",header = T,row.names = 1)
tcga_gsva <- as.data.frame(t(gsva(as.matrix(expr), immunity, method = "ssgsea")))
write.csv(tcga_gsva,"illumina_41gsva.csv",quote = F)


filterCommonGenes(input.f="affy_uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")


scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="affy_scores.txt",sep="\t",quote=F,col.names=F)


ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}

gene_expression<-read.table("affy_uniq.symbol.txt",row.names=1,header=TRUE, sep="\t", dec = ".",check.names=FALSE)
sample_names<-names(gene_expression)

## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
# For different 
IPSG<-read.table("IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes<-as.vector(unique(IPSG$NAME))

IPS<-NULL
MHC<-NULL
CP<-NULL
EC<-NULL
SC<-NULL
AZ<-NULL

GVEC<-row.names(gene_expression)

VEC<-as.vector(IPSG$GENE)

ind<-which(is.na(match(VEC,GVEC)))

MISSING_GENES<-VEC[ind]
dat<-IPSG[ind,]
if (length(MISSING_GENES)>0) {
  cat("differently named or missing genes: ",MISSING_GENES,"\n")
}
for (x in 1:length(ind)) {
  print(IPSG[ind,])
}

for (i in 1:length(sample_names)) {	
  GE<-gene_expression[[i]]
  mGE<-mean(GE)
  sGE<-sd(GE)
  Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1<-IPSG$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k<-k+1
  }
  WG<-MIG*WEIGHT
  WG[is.na(WG)]<-0   
  MHC[i]<-mean(WG[1:10])
  CP[i]<-mean(WG[11:20])
  EC[i]<-mean(WG[21:24])
  SC[i]<-mean(WG[25:26])
  AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i]<-ipsmap(AZ[i])
}

DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
write.table(DF,file="affy_IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")

##Combined the above results as affy_traindata.txt


### Verify the classifier

