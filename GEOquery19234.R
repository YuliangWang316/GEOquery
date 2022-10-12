library(GEOquery)
gds <- getGEO("GSE19234")
exp<-cbind(gds$GSE19234_series_matrix.txt.gz@assayData$exprs,as.data.frame(gds$GSE19234_series_matrix.txt.gz@featureData@data$`Gene Symbol`))
exp<-na.omit(exp)
colnames(exp)[length(colnames(exp))]<-"Symbol"
gene_symbol<-unique(exp$Symbol)
library(dplyr)
a<-filter(exp,Symbol == gene_symbol[1])
a<-a[,-length(colnames(a))]
b<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
rownames(b)<-gene_symbol[1]
remove(a)
for (i in 2:length(gene_symbol)) {
  a<-filter(exp,Symbol == gene_symbol[i])
  a<-a[,-length(colnames(a))]
  c<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
  remove(a)
  rownames(c)<-gene_symbol[i]
  b<-rbind(b,c)
  remove(c)
}
library(GSVA)
genset<-read.table("c:/Users/xjmik/Downloads/Treg_1C_signature.txt",sep = "\t",header = TRUE)
geneset_list<-as.list(genset)
remove(exp,genset,gene_symbol)


remove(i)


gsva<-gsva(expr = as.matrix(b),gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20)
gsva<-as.data.frame(t(gsva))
remove(b,geneset_list)
clinicaldata<-as.data.frame(cbind(gds$GSE19234_series_matrix.txt.gz@phenoData@data$characteristics_ch1.1,gds$GSE19234_series_matrix.txt.gz@phenoData@data$characteristics_ch1.5))
rownames(clinicaldata)<-rownames(gds$GSE19234_series_matrix.txt.gz@phenoData@data)
a<-as.data.frame(t(as.data.frame(strsplit(clinicaldata[,1],":"))))
a<-as.data.frame(a[,-1])
rownames(a)<-rownames(clinicaldata)
colnames(a)<-"OS"
for (i in 2:length(colnames(clinicaldata))) {
  b<-as.data.frame(t(as.data.frame(strsplit(clinicaldata[,i],":"))))
  b<-as.data.frame(b[,-1])
  rownames(b)<-rownames(clinicaldata)
  a<-cbind(a,b)
}
remove(b,i)
colnames(a)<-c("OS","events")
remove(clinicaldata)
remove(gds)
df<-cbind(a,gsva)
remove(a,gsva)
library(survival)
library(survminer)
library(dplyr)

df$OS<-as.numeric(df$OS)
df$events<-as.numeric(df$events)
res.cut<-surv_cutpoint(df,time = "OS",event = "events",variables = "Treg_1C" )
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(OS,events)~ Treg_1C,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
res.cox<-coxph(Surv(OS,events) ~ Treg_1C,data=res.cat)
summary(res.cox)
test.ph<-cox.zph(res.cox)
ggcoxzph(test.ph)
