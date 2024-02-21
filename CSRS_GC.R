library(iClusterPlus)
library(GSVA)
library(survival)
library(survminer)
library(parallel)
library(dplyr)
library(readxl)
library(maftools)
library(psych)
library(NMF)
cop_gene = c("ATP7A","ATP7B", "CDKN2A","DBT","DLAT","DLD","DLST","FDX1","GCSH","GLS", "LIAS","LIPT1" ,"MTF1", "PDHA1","PDHB","SLC31A1")

######################### Cuproptosis signature subtype ########################
setwd("./data_set")
acrg_cli <- read.csv("acrg_cli.csv")   # ACRG cohort curated from GSE62254
exp62 <- readRDS("exp62.Rds")  # ACRG cohort curated from GSE62254
exp62 <- exp62[,acrg_cli$Tumor_Sample_Barcode]
exp62t <- as.data.frame(t(exp62))
result = iCluster2(x=list( t(exp62[cop_gene,acrg_cli$Tumor_Sample_Barcode])), K=4,lambda = 0.1,method = "lasso") 
plotiCluster(result)
acrg_clis = cbind(acrg_cli,cluster = result$clusters)
acrg_clis$CSC_cluster[acrg_clis$cluster == "1"] = "CSC1"   # The subtype name may be need adjust
acrg_clis$CSC_cluster[acrg_clis$cluster == "2"] = "CSC2"
acrg_clis$CSC_cluster[acrg_clis$cluster == "3"] = "CSC3"
acrg_clis$CSC_cluster[acrg_clis$cluster == "4"] = "CSC4"
ggsurvplot(survfit(Surv(Surv_time,Surv_status)~cluster,acrg_clis),pval = TRUE,risk.table = TRUE,pval.method = T,palette =  ggsci::pal_nejm()(6),ggtheme = theme_bw(), title = "survival")

acrg_clis$Mol_Subtype <- acrg_clis$Mol..Subtype..0.MSS.TP53...1.MSS.TP53...2...MSI..3..EMT
annotation_col = acrg_clis[,c("Tumor_Sample_Barcode","MSI", "EBV_group","OS_status", "Mol_Subtype","pStage", "Gender","AGE_group","CSC_cluster")]
dplyr::arrange(annotation_col,CSC_cluster)->annotation_col
rownames(annotation_col) <- annotation_col$Tumor_Sample_Barcode
annotation_col[-1]->annotation_col
ann_colors = list(
  CSC_cluster = c("CSC1"="#BC3C29FF","CSC2"="#0072B5FF","CSC3"="#E18727FF", "CSC4"="#20854EFF"),
  MSI = c("negative" = "#00D8EF", "positive" = "#DDF0E4"),
  EBV_group = c("0" = "#DDF0E4", "1" = "#E199FF"),
  OS_status = c("0" = "#DDF0E4", "1" = "#3484BE"),
  Mol_Subtype = c("0"="#DDF0E4","1" ="#2CA25F","2" = "#EAB736" ,"3" = "#73A7F9"),
  pStage= c("I"="grey70","II"="lightblue3","III"="dodgerblue2","IV"="blue4"),
  Gender=c(F="#DDF0E4",M="#EAB736"), AGE_group=c(xiao60="#DDF0E4",dadeng60="#62BCFF")
)
norm_4heatmap = norm_4heatmap <- function(exp) {
  scale(exp)->exp1
  exp1[exp1< -2] <- -2
  exp1[exp1> 2] <- 2
  normalization<-function(x){
    return(-1+2*(x-min(x))/(max(x)-min(x)))}
  exp1 <- normalization(exp1)
  return(exp1)
}
aa = t(norm_4heatmap(exp62t[rownames(annotation_col),cop_gene]))
pheatmap(aa,
         show_colnames = F,cluster_cols = F,border=FALSE,
         cluster_rows = T,clustering_method = "ward.D2",treeheight_row = 0,
         annotation_col = annotation_col,annotation_colors = ann_colors)

ggplot(acrg_clis,aes(x=CSC_cluster,fill=Mol_Subtype))+geom_bar(stat ="count",position = "fill")+
  scale_fill_manual(values = c("#DDF0E4","#2CA25F", "#EAB736" , "#73A7F9"))+labs(x = '', y = 'Fraction of Patients')+coord_flip()
ggplot(acrg_clis,aes(x=CSC_cluster,fill=LAUREN.1.intestinal..2.diffuse..3.mixed))+geom_bar(stat ="count",position = "fill")+
  scale_fill_manual(values = c("#DDF0E4","#2CA25F", "#EAB736" , "#73A7F9"))+labs(x = '', y = 'Fraction of Patients')+coord_flip()


############################## KEGG pathway Heatmap in 4 subtype #######################
mean_data = read.csv("Fig3a_pathway.csv",row.names = 1)
annotation_col = data.frame(dataset=c(replicate(4,c("ACRG","TCGA","Yonsei"))),CSC_cluster=c(replicate(3,"CSC1"),replicate(3,"CSC2") ,replicate(3,"CSC3"),replicate(3,"CSC4")))
rownames(annotation_col)=colnames(mean_data)
ann_colors = list(
  CSC_cluster = c("CSC1"="#BC3C29FF","CSC2"="#0072B5FF","CSC3"="#E18727FF", "CSC4"="#20854EFF")
)
pheatmap::pheatmap(mean_data,cluster_cols = F,annotation_col = annotation_col, annotation_colors =  ann_colors,  clustering_method = "ward.D",show_rownames = T,
                   show_colnames = F,treeheight_row=20, gaps_col = c(3, 6, 9) )


############################## Immune Cell Infiltration ########################
gene_set <- as.data.frame(read_excel("mmc3.xlsx",sheet=1))
gene_set$`Cell type` = gsub(" ","_", gene_set$`Cell type`)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])   
acrg_im28 <- GSVA::gsva(as.matrix(exp62), list, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE) 
acrg_im28 = as.data.frame(t(acrg_im28 ))
aa <- sapply(colnames(acrg_im28)[1:28],function(i)kruskal.test(acrg_im28[,i],acrg_clisc$CSC_cluster)$p.value)
exp_cell = cbind( acrg_im28[,names(aa[aa<0.01])],"group"=acrg_clisc$CSC_cluster)
melt_cell <- reshape2::melt(exp_cell,id.vars="group")
ggplot(data = melt_cell,aes(x=variable,y=value,fill=group))+geom_boxplot(width=0.7,outlier.size = 0.4)+labs(title="",x="",y="Estimated AUC Value")+ #scale_y_continuous(limits = c(0,0.6))+
  theme(axis.title=element_text(face="bold",size=11),axis.text.x =element_text(face="bold",size=9,angle=45,hjust = 0.9,vjust=0.95))+scale_fill_manual(values = ggsci::pal_nejm()(6))+stat_compare_means(aes(group = group), label = "p.signif")


######### Construction of the cuproptosis signature risk score(CSRS) ########### 
dt<-as.matrix(scale(exp62t[acrg_clis$Tumor_Sample_Barcode,cop_gene]))  # 原始PCA方法
y=eigen(cor(dt))  # 用标准化后的dt，求出cor(dt)的特征值和特征向量
y$values 
y$values/sum(y$values)
y$vectors[,1:10] 
U<-as.matrix(y$vectors)
PC <- dt %*% U    
CSRS=0.0
for (i in 1:5)
  CSRS=(y$values[i]*PC[,i])/(sum(y$values[1:5]))+CSRS

acrg_clisc = cbind(acrg_clis,CSRS =  CSRS)
cutpoint = surv_cutpoint(acrg_clisc,time="Surv_time",event="Surv_status",variables="CSRS")
cutpoint
bb = surv_categorize(cutpoint, variables = "CSRS", labels = c("alow", "high"))
acrg_clisc$CSRS_group = bb$CSRS
survfit <- survfit(Surv(Surv_time,Surv_status)~CSRS_group, acrg_clisc)
pdf("surv_CSRS_acrg.pdf",height = 5,width = 5,onefile = F)
ggsurvplot(survfit,pval = TRUE,risk.table = TRUE,legend.title = "clin681",
           palette =  c("#B24745FF","#374E55FF"),ggtheme = theme_bw(), title = "survival")
dev.off()
pdf("CSRS_distribution.pdf",height = 4,width = 5,onefile = F)
plot(cutpoint)
dev.off()

########################## ggalluvial   #######################################
library(ggalluvial)
group39 <- as.data.frame( group_by(acrg_clisc, CSC_cluster, Mol_Subtype, CSRS_subgroup) %>% summarise(., count = n()) )
pdf("ggalluvial1.pdf",width = 6,height = 5)
ggplot(group39,aes(axis1 = CSC_cluster, axis2 = Mol_Subtype, axis3 = CSRS_subgroup,y= count,fill=CSC_cluster)) +
  scale_x_discrete(limits = c("MS cluster", "Mol subtype", "MS score"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = as.factor(CSC_cluster))) +
  geom_stratum() + geom_text(stat = "stratum", infer.label = TRUE) +
  theme_minimal() + scale_fill_manual(values = ggsci::pal_nejm()(6) )+
  ggtitle("Patients in the ACRG cohort")
dev.off()

############################ Tumor Mutation in CSRS_subgroup #######################
tcga_clisc = read.csv("tcga_clisc.csv")
tcga_maf = read.maf(maf="tcga_stad_hg19.txt",clinicalData =tcga_clisc)
tcga_maf_alow <-subsetMaf(maf=tcga_maf, tsb=tcga_clisc$Tumor_Sample_Barcode[tcga_clisc$CSRS_subgroup== "alow"], isTCGA=T)
tcga_maf_high <-subsetMaf(maf=tcga_maf, tsb=tcga_clisc$Tumor_Sample_Barcode[tcga_clisc$CSRS_subgroup== "high"], isTCGA=T)
genes <- c("TP53","CDH1","SMAD4","PIK3CA","RHOA","ARID1A","KRAS","MUC6","APC","COL11A1", "BCOR","EYA4","BNC2","RNF43","ABCA10","CTNNB1","MACF1","SMAD2","SOHLH2","RASA1","FAM46D","PLB1","CNGA4","EIF2C4","ERBB2","PTPRC","PTEN","ERBB3","HLA-B","PTPN23")     
genes <- intersect(getGeneSummary(tcga_maf)$Hugo_Symbol,genes) #调整瀑布图的基因排列按照突变由多到少排序
pdf("waterfall_tcga.pdf",height = 10,width = 10)
coOncoplot(m1 = tcga_maf_alow, m2 = tcga_maf_high, m1Name = 'alow', m2Name = 'high', genes = genes, removeNonMutated = F,borderCol= NULL,
           clinicalFeatures1 = c("AGE_group","SEX","STAGE", "SUBTYPE","imm_Immune.Subtype","MFP_an","CSC_cluster","CSRS_subgroup"),
           annotationColor1 = list(
             CSRS_subgroup = c("alow" = "#B24745FF", "high" = "#374E55FF"),
             CSC_cluster = c("CSC1"="#BC3C29FF","CSC2"="#0072B5FF","CSC3"="#E18727FF", "CSC4"="#20854EFF"),
             MSI = c("MSS" = "#DDF0E4","MSI" = "#00D8EF"),
             MFP_an = c("D" = "#DDF0E4", "F" = "#2CA25F","IE" = "#CD534CFF","IE/F" = "#7AA6DCFF"),
             OS_status = c("0" = "#DDF0E4", "1" = "#3484BE"),
             SUBTYPE = c("STAD_CIN"="#DDF0E4","STAD_EBV" ="#2CA25F","STAD_GS" = "#EAB736" ,"STAD_MSI" = "#73A7F9"),
             imm_Immune.Subtype = c("C1"="#7570B3","C2" ="#E7298A","C3" = "#66A61E" ,"C4" = "#00A1D5","C6" = "#6A6599FF"),
             STAGE= c("I"="grey70","II"="lightblue3","III"="dodgerblue2","IV"="blue4"),
             SEX=c(Female="#DDF0E4",Male="#EAB736"), AGE_group=c(xiao60="#DDF0E4",dadeng60="#62BCFF")),
           clinicalFeatures2 = c("AGE_group","SEX","STAGE", "SUBTYPE","imm_Immune.Subtype","MFP_an","CSC_cluster","CSRS_subgroup"),
           annotationColor2 = list(
             CSRS_subgroup = c("alow" = "#B24745FF", "high" = "#374E55FF"),
             CSC_cluster = c("CSC1"="#BC3C29FF","CSC2"="#0072B5FF","CSC3"="#E18727FF", "CSC4"="#20854EFF"),
             MSI = c("MSS" = "#DDF0E4","MSI" = "#00D8EF"),
             MFP_an = c("D" = "#DDF0E4", "F" = "#2CA25F","IE" = "#CD534CFF","IE/F" = "#7AA6DCFF"),
             OS_status = c("0" = "#DDF0E4", "1" = "#3484BE"),
             SUBTYPE = c("STAD_CIN"="#DDF0E4","STAD_EBV" ="#2CA25F","STAD_GS" = "#EAB736" ,"STAD_MSI" = "#73A7F9"),
             imm_Immune.Subtype = c("C1"="#7570B3","C2" ="#E7298A","C3" = "#66A61E" ,"C4" = "#00A1D5","C5" = "#6A6599FF"),
             STAGE= c("I"="grey70","II"="lightblue3","III"="dodgerblue2","IV"="blue4"),
             SEX=c(Female="#DDF0E4",Male="#EAB736"), AGE_group=c(xiao60="#DDF0E4",dadeng60="#62BCFF")),
)
dev.off()

maf_cell_tnm <- trinucleotideMatrix(maf=tcga_maf,prefix = 'chr', add = TRUE, ref_genome="BSgenome.Hsapiens.UCSC.hg19")
maf_cell_sig <- extractSignatures(mat=maf_cell_tnm,n=5, parallel = 14, plotBestFitRes=FALSE)
laml.og30.cosm = compareSignatures(nmfRes = maf_cell_sig, sig_db = "legacy")
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
maftools::plotSignatures(nmfRes = maf_cell_sig, title_size = 1.2, sig_db = "legacy",contributions = FALSE)

