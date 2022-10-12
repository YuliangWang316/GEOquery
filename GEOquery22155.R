library(GEOquery)
gds <- getGEO("GSE22155")
exp<-cbind(gds$`GSE22155-GPL6102_series_matrix.txt.gz`@assayData$exprs,as.data.frame(gds$`GSE22155-GPL6102_series_matrix.txt.gz`@featureData@data$Symbol))
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

exp<-cbind(gds$`GSE22155-GPL6947_series_matrix.txt.gz`@assayData$exprs,as.data.frame(gds$`GSE22155-GPL6947_series_matrix.txt.gz`@featureData@data$Symbol))
colnames(exp)[length(colnames(exp))]<-"Symbol"
gene_symbol<-unique(exp$Symbol)
library(dplyr)
a<-filter(exp,Symbol == gene_symbol[1])
a<-a[,-length(colnames(a))]
d<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
rownames(d)<-gene_symbol[1]
remove(a)
for (i in 2:length(gene_symbol)) {
  a<-filter(exp,Symbol == gene_symbol[i])
  a<-a[,-length(colnames(a))]
  c<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
  remove(a)
  rownames(c)<-gene_symbol[i]
  d<-rbind(d,c)
  remove(c)
}
remove(gene_symbol,exp)
remove(i)
e<-intersect(rownames(b),rownames(d))
b_new<-b[e,]
d_new<-d[e,]
f<-cbind(b_new,d_new)
remove(b,b_new,d,d_new,e)

gsva<-gsva(expr = as.matrix(f),gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20)
gsva<-as.data.frame(t(gsva))
remove(f,geneset_list)
clinicaldata1<-as.data.frame(cbind(gds$`GSE22155-GPL6102_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.1,gds$`GSE22155-GPL6102_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.2,gds$`GSE22155-GPL6102_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.3,gds$`GSE22155-GPL6102_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.4,gds$`GSE22155-GPL6102_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.11))
rownames(clinicaldata1)<-rownames(gds$`GSE22155-GPL6102_series_matrix.txt.gz`@phenoData@data)
clinicaldata2<-as.data.frame(cbind(gds$`GSE22155-GPL6947_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.1,gds$`GSE22155-GPL6947_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.2,gds$`GSE22155-GPL6947_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.3,gds$`GSE22155-GPL6947_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.4,gds$`GSE22155-GPL6947_series_matrix.txt.gz`@phenoData@data$characteristics_ch1.11))
rownames(clinicaldata2)<-rownames(gds$`GSE22155-GPL6947_series_matrix.txt.gz`@phenoData@data)
clinicaldata<-rbind(clinicaldata1,clinicaldata2)
remove(clinicaldata1,clinicaldata2)
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
colnames(a)<-c("OS","events","gender","age","stage")
remove(clinicaldata)
remove(gds)
df<-cbind(a,gsva)
remove(a,gsva)
library(survival)
library(survminer)
library(dplyr)
df<-dplyr::filter(df,OS != " -")
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
