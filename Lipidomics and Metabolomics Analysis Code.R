install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
rm(list = ls())
185.199.111.133 raw.githubusercontent.com
conda activate qiime2-2020.8
library(tidyverse)
library(data.table)
library(pmp)
library(attempt)
#install.packages("hacksaw")
library(randomForest)
library("RColorBrewer")
library(ggpubr)
library(psych)
library(ropls)
library(ComplexHeatmap)
library(readxl)
QCdo <- function(x2,delRSD=0.5,deloutlier="none",
                 order=1:ncol(x2)) {
  
  
  
  x2<-as.data.frame(x2)
  batch=rep(1,nrow(x2))
  t(x2[,-c(1:2)])->l
  l
  x2[,-c(1:2)] <-  t(QCRSC(df =t(x2[,-c(1:2)])%>%as.data.frame(), 
                           order = order, 
                           batch =batch,
                           qc_label='Q',
                           classes = x2$Groups,
                           spar = 1,minQC = 4))
  
  ##QC CV变异系数 cv=0.5
  cv=delRSD
  if(cv!=1){apply(x2[x2$Groups=='Q',-c(1:2)],2,function(x){
    x<-x[!is.na(x)&x!=0]
    if(length(x)==0)return(1)
    return(sd(x)/mean(x))
  })->y3
    x2[,c(1,2,which(y3<cv)+2)]->x2
  }
  # IQR 异常值
  if(deloutlier=="IQR"){apply(x2[,-c(1:2)],2,function(x){
    quantile(x, c(0.25, 0.75),na.rm = T)->x4
    
    x[!is.na(x)&x>(3*(x4[2]-x4[1])+x4[2])]<-NA
    x[!is.na(x)&x<(x4[1])-3*(x4[2]-x4[1])]<-NA
    return(x)
  })->x2[,-c(1:2)]
    x2[,-c(1:2)][is.na(x2[,-c(1:2)])]<-0
    narm=0.8
    x2[,which((colSums(x2!=0))>(narm*nrow(x2)))]->x2
  } 
  return(x2)
}

  
  
  
delRSD<- function(x2,delRSD=0.5) {
  cv=delRSD
  if(cv!=1){apply(x2[x2$Groups=='Q',-c(1:2)],2,function(x){
    x<-x[!is.na(x)&x!=0]
    if(length(x)==0)return(1)
    return(sd(x)/mean(x))
  })->y3
    x2[,c(1,2,which(y3<cv)+2)]->x2
  }
  return(x2)
  
}

deloutlier <- function(x2,deloutlier="IQR") {
  # IQR 异常值
  if(deloutlier=="IQR"){apply(x2[,-c(1:2)],2,function(x){
    quantile(x, c(0.25, 0.75),na.rm = T)->x4
    
    x[!is.na(x)&x>(3*(x4[2]-x4[1])+x4[2])]<-NA
    x[!is.na(x)&x<(x4[1])-3*(x4[2]-x4[1])]<-NA
    return(x)
  })->x2[,-c(1:2)]
    x2[,-c(1:2)][is.na(x2[,-c(1:2)])]<-0
    narm=0.8
    x2[,which((colSums(x2!=0))>(narm*nrow(x2)))]->x2
  } 
  return(x2)
}

imputevery <- function(i=X,minto=1,scale=0) {
  if(minto){
    3:ncol(x)%>% map(function(i){ x %>% split(.$Group) %>%
        purrr::transpose()%>% .[[i]] %>%
        (function(x) combn(names(x), m=1, function(y) { 
          if(!is.null(attempt::if_all(flatten_dbl(x[y[1]]),~ .x!=0, ~ return(1)))) {return(3)}
          if(!is.null(attempt::if_all(flatten_dbl(x[y[1]]),~ .x==0, ~ return(2)))) {return(2)}
          return(0)
        }
        )
        )})%>%data.frame%>%t->o
    colnames(o)<-names(table(x$Groups))
    sum(rowSums(o==2)>0)
    rowSums(o==3)>0
    x<-x[order(x$Groups),]
    x%>%as.data.frame%>%.[,-c(1:2)] %>%.[.>0.01] %>%.[.>0.01]%>%min%>%{./2}->min
    rownames(o)<-colnames(x)[-c(1:2)]
    apply(o, 1, function(o) {
      if(sum(o==3)){if(sum(o==2))return(1); return(0)};return(0)
    })->d
    d[d==1]
    (which(d==1)+2)->mi
    x<-as.data.frame(x)
    x[,c(mi,2)]->y
    i=1
    #最小值的1/2 替代
    o<-as.data.frame(o)
    o[(which((rowSums(o==3)>0)&((rowSums(o==2)>0)))),]->o2
    colnames(o)->ocolname
    #for (i in seq_along(mi)) {
    # ocolname[o2[i,]==1]->nam
    #y[y$Groups%in%nam,i]<-min
    #}
    if((nrow(o2)!=0)){
      x[,c(2,(which((rowSums(o==3)>0)&((rowSums(o==2)>0)))+2))]%>%
        as.data.frame()->data
      if(ncol(data)==2){
        data[data$Groups== 
               colnames(o)[o[(which((rowSums(o==3)>0)&((rowSums(o==2)>0)))),]%>%
                             as.data.frame()%>%{.==2}]
             ,2]<-min
      }
      if(ncol(data)!=2){
        for (i in 2:ncol(data)) {
          data[data$Groups%in% 
                 colnames(o)[o[(which((rowSums(o==3)>0)&((rowSums(o==2)>0)))),]%>%
                               as.data.frame()%>%.[i-1,]%>%{.==2}]
               ,i]<-min
        }
      }
      cbind(x[,-c((which((rowSums(o==3)>0)&((rowSums(o==2)>0)))+2))],data[,-1])->data
      return(data)
      
    }
    if((nrow(o2)==0)){
      return(x)
    }
  }
}

scaleevery <- function(x2,scale="Internalstandard",Internalstandard.lable="Docosane") {
  if(scale=="pqn") t(pqn_normalisation(t(x2[,-c(1:2)]),x2$Groups,qc_label='Q'))->x2[,-c(1:2)]
  if(scale=="median"){
    da<-"median"
    eval(parse(text=da))(as.numeric(x2[3,-c(1:2)]))->meding
    for (i in 1:nrow(x)) {
      x2[i,-c(1:2)]<-x2[i,-c(1:2)]/eval(parse(text=da))(as.numeric(x2[i,-c(1:2)]))
    }
    meding*x2[,-c(1:2)]->x2[,-c(1:2)]
    
  }
  if(scale=="sum"){
    da<-"sum"
    eval(parse(text=da))(as.numeric(x2[3,-c(1:2)]))->meding
    for (i in 1:nrow(x)) {
      x2[i,-c(1:2)]<-x2[i,-c(1:2)]/eval(parse(text=da))(as.numeric(x2[i,-c(1:2)]))
    }
    meding*x2[,-c(1:2)]->x2[,-c(1:2)]
    
  }
  
  if(scale=="Internalstandard"){
    #内标 加 归一化
    da<-"median"
    x2
    which(colnames(x2)%in%Internalstandard.lable)->k
    x2<-as.data.frame(x2)
    x2[,which(colnames(x2)==Internalstandard.lable)]->y
    x2[which(y!=0),-c(1:2)]->Scale 
    x2[which(y==0),-c(1:2)]->noScale 
    median(as.numeric(x2[,which(colnames(x2)==Internalstandard.lable)]))->meding
    # T
    for (i in which(y!=0)) {
      x2[i,-c(1:2)]<-x2[i,-c(1:2)]/median(x2[i,which(colnames(x2)==Internalstandard.lable)][x2[i,which(colnames(x2)==Internalstandard.lable)]!=0])
    }
    meding*x2[which(y!=0),-c(1:2)]->x2[which(y!=0),-c(1:2)]
    # F
    eval(parse(text=da))(apply(x2[which(y!=0),-c(1:2)],1,eval(parse(text=da))))->meding
    apply(x2[which(y==0),-c(1:2)],1,function(x){x/median(x)})*meding->x2[which(y==0),-c(1:2)]
    
  }
  return(x2)
}

op<-function(x,y,x2, orthoI =1,scaleC ="pareto",predI = 2,crossvalI=NA){
  if((x2$Groups%>%table()%>%{.==1}%>%sum())>0){cat("using plsda ,because your group have only one");orthoI=NA}
  rownames(x2)<-x2$name
  opls(x2[x2$Groups%in%c(x,y),-c(1:2)],
       as.factor(x2$Groups[x2$Groups%in%c(x,y)]),
       orthoI =orthoI,scaleC =scaleC,
       crossvalI = ifelse(is.na(crossvalI),length(x2$Groups[x2$Groups%in%c(x,y)])-1,crossvalI) 
       ,plotL=FALSE,printL=FALSE)->l
  dir.create(paste(x,"vs",y))
  pdf(paste(x," vs ",y,"/",x,"vs",y," opls",".pdf",sep = ""))
  plot(l)
  dev.off()
  l%>%return
}

test_data<-function(x2,modoel="t.test",rev=F) {
  if(x2$Groups%>%table()%>%{.<2}%>%sum()%>%{.>0}){cat("you group is only one");stop()}
  if(x2$Groups%>%table()%>%length()%>%{.<2}){cat("you group is only one");stop()}
  if(modoel=="t.test"){
    return(2:(ncol(x2)-1)%>%map(function(i){(x2[,-c(1)] %>% split(.$Groups) %>%
                                               purrr::transpose() %>%
                                               .[[i]] %>%
                                               (function(x) combn(names(x), m=2, function(y){
                                                 if(mean(purrr::flatten_dbl(x[y[1]])%>%
                                                         as.numeric())!= mean(purrr::flatten_dbl(x[y[2]])%>%as.numeric())){
                                                   return( t.test(purrr::flatten_dbl(x[y[1]])%>%as.numeric(),
                                                                  purrr::flatten_dbl(x[y[2]])%>%as.numeric())[["p.value"]])
                                                   
                                                 }
                                                 
                                                 return(1.1)} ))%>%
                                               unlist)})%>%
             data.frame%>%t%>%as.data.frame())
  }
  if(modoel=="wilcox.test"){
    return(2:(ncol(x2)-1)%>%map(function(i){(x2[,-c(1)] %>% split(.$Groups) %>%
                                               purrr::transpose() %>%
                                               .[[i]] %>%
                                               (function(x) combn(names(x), m=2, function(y){
                                                 if(mean(purrr::flatten_dbl(x[y[1]])%>%
                                                         as.numeric())!= mean(purrr::flatten_dbl(x[y[2]])%>%as.numeric())){
                                                   return( wilcox.test(purrr::flatten_dbl(x[y[1]])%>%as.numeric(),
                                                                       purrr::flatten_dbl(x[y[2]])%>%as.numeric())[["p.value"]])
                                                   
                                                 }
                                                 
                                                 return(1.1)} ))%>%
                                               unlist)})%>%
             data.frame%>%t%>%as.data.frame())
  }
  
  
  
  if(modoel=="kruskal.test")return(x2[,-c(1:2)]%>% map(~ kruskal.test(.x~ x2$Groups ))%>%
                                     lapply(function(x){x[["p.value"]]})%>%as.numeric())#{return(}
  if(modoel=="aov.test")return(x2[,-c(1:2)]%>%
                                 map(~ aov(.x~ x2$Groups ))%>%map(~summary(.x)[[1]][["Pr(>F)"]][1])%>%as.numeric())#return()
  if(modoel=="FC"){return(2:(ncol(x2)-1)%>%map(function(i){(x2[,-c(1)] %>% split(.$Groups) %>%
                                                              purrr::transpose() %>%
                                                              .[[i]] %>%
                                                              (function(x) combn(names(x), m=2, function(y) 
                                                                ifelse(rev,c((2^mean(c(flatten_dbl(x[y[1]])))-1)/
                                                                               (2^mean(c(flatten_dbl(x[y[2]])))-1)),
                                                                       c((2^mean(c(flatten_dbl(x[y[2]])))-1)/
                                                                           (2^mean(c(flatten_dbl(x[y[1]])))-1))
                                                                )
                                                              ))%>%
                                                              unlist)})%>%
                            data.frame%>%t)}
  cat("t.test  wilcox.test kruskal.test aov.test")
  
  
}



#数据读取
"C:/Users/靓仔/Desktop/1"%>%setwd()
data<-read.csv("C:/Users/靓仔/Documents/更正14.csv",header=F)
#data<-read_excel("C:/Users/联想/Desktop/师兄 陈序/血清.xlsx",sheet =1,range="E1:BO462",col_names = TRUE)
#data<-cbind(k,data)
data[1,1]="name"
colnames(data)<-data[1,]
data<-data[-1,]
#data<-data[order(-data$"M2-XQ-QC1"),]
#data[!duplicated(data$`Metabolite name`),]->data
dir.create("C:/Users/靓仔/Desktop/1/data",recursive = T)
setwd("C:/Users/靓仔/Desktop/1/data")
data%>%writexl::write_xlsx(.,"date.xlsx")

#开始

setwd("C:/Users/靓仔/Desktop/1/data")
data<-readxl::read_xlsx("date.xlsx")
dir.create("C:/Users/靓仔/Desktop/1/resultM2",recursive = T)
setwd("C:/Users/靓仔/Desktop/1/resultM2")
#k<-data[1]
#ee1<-which(names(data)=='M14-XQ-QC1')
#ee2<-which(names(data)=='M14-XQ-D97')
#data<-data[,((ee1):(ee2))]
#data<-cbind(k,data)

dir.create("data")
#data%>%writexl::write_xlsx(.,"data/data.xlsx")
#data<-readxl::read_xlsx("data/data.xlsx")
#data<-data[,-(11:12)]
#data.table(name=colnames(data)[-1],
#           Groups=c(rep("Q",1),rep("A",10),rep("Q",1)
#                    ,rep("B",10),rep("Q",1),rep("D",10)
#                    ),#,rep("Q",4)
#           t(data[,-1]))->date

#colnames(date)[-c(1:2)]<-data$`Name`
data%>%writexl::write_xlsx(.,"data/date.xlsx")


data[,1:3]->group
#group[1:(ee2+1-ee1),3]=1:(ee2+1-ee1)
colnames(group)[3]="order"
group%>%writexl::write_xlsx(.,"data/group.xlsx")
group<-readxl::read_xlsx("data/group.xlsx")
data%>%as.data.frame()->x
x[is.na(x)]<-0
x[x==""]<-0
dim(x)
#找出全部缺失值和非全部缺失 2全部缺失，3全部有，0两个中间
x[,3:dim(x)[1]]<-as.data.frame(lapply(x[,3:dim(x)[1]],as.numeric))
colnames(x)[2]<-"Groups"
cx<-colnames(x)
rx<-rownames(x)
x[,-(1:2)]<-as.data.frame(lapply(x[,-(1:2)],as.numeric))
rownames(x)<-rx
colnames(x)<-cx
imputevery(x,T)->x
x->x2
#缺失值删除 输入narm
{ 
  narm=0.6
  x[,which((colSums(x!=0))>(narm*nrow(x)))]->x2
  nrow(x)->qso
  x2[1:qso,]->x2
}
#x2[,-c(1:2)]<-x2[,-c(1:2)]*( 10^6)
scaleevery(x2,"median")->x2#标准化 pqn median sum 没有QC样本选 median
#"Internalstandard","Docosane"
x2[x2=="NaN"]=0
x2[,-c(1:2)][x2[,-c(1:2)]==0]<-NA
x2[,-c(1:2)] <- t(mv_imputation(df=t(x2[,-c(1:2)]), 
                                method="rf", rowmax = 1, colmax = 1))#缺失值插补

group[order(group$name),]->group
group[(group$name%in%x2$name),]->group
# QC 矫正 QC CV变异系数 cv=0.5  IQR 异常值
QCdo(x2,1,"none",order=group[order(group$name),]$order)->x2 #如果没有QC 可以不走
cv <- function(.x){sd(.x)/mean(.x)}
x2[x2$Groups=="Q",-c(1:2)]%>%apply(2,cv)->k
x2[,c(1,2,which(k<0.5)+2)]->x2

# 插补 QC矫正后有缺失值
x2[,-c(1:2)][x2[,-c(1:2)]==0]<-NA
x2[,-c(1:2)] <- t(mv_imputation(df=as.matrix(t(x2[,-c(1:2)])), method="rf",
                                rowmax = 1, colmax = 1))


#箱线图
rownames(x2)<-x2$name
x2[order(x2$Groups),]->smaple
smaple[,-c(1:2)]%>%t->data2
colnames(data2)<-smaple$name
boxplot(log2(data2+1), las = 2,
        col=rep(brewer.pal(n = length(table(smaple$Groups)), name = "RdBu"),
                table(smaple$Groups))
)
#log 化
#x2<-x2[-3,]
x2[,-c(1:2)]<-log2(x2[,-c(1:2)]+1)
colnames(x2)->name



pca_result <- opls(x = x2[,-c(1:2)],predI = 2,scaleC ="pareto")


scoreMN <- pca_result@scoreMN
scoreMN<-as.data.frame(scoreMN)
scoreMN <- data.frame(scoreMN,x2$Groups)
scoreMN$samples <- (x2$name)
pca<-ggscatter(scoreMN, x = "p1", y = "p2",
               color = "black", 
               fill = "x2.Groups",
               shape = 21,
               ellipse = T, 
               size = 5,
               ellipse.level = 0.95,
               mean.point = F,
               alpha = 0.8,
               font.label = 10, repel = TRUE,
               xlab = paste0("PC1 (",pca_result@modelDF$R2X[1]*100,"%)"),
               ylab = paste0("PC2 (",pca_result@modelDF$R2X[2]*100,"%)"),
               font.main = c(16, "plain", "black"),            
               font.x = c(16, "plain", "black"),                  
               font.y = c(16, "plain", "black"),                
               legend = "top",  
               legend.title = "Treatment",  
               font.legend = c(16, "plain", "black"),               
               rotate = F,                                         
               ticks = T,   
               tickslab = T, 
               label="samples",
               palette = "aaas") +
  theme(axis.text = element_text(size = 16))
dir.create("output")
pdf("output/pca.pdf")
pca
dev.off()
x2->x

x->x2
#标准化
colnames(x2)->name
x2[x2$Groups!="Q",]->x2
x2[x2$Groups!="Q",-c(1:2)]%>%t%>%scale()%>%t->x2[x2$Groups!="Q",-c(1:2)]
colnames(x2)[-c(1:2)]<-paste("A",1:length(colnames(x2)[-c(1:2)]),sep = "")
randomForest(as.factor(Groups) ~ ., data=x2[,-1])->rf_classifier
importance(rf_classifier)->rf
rf<-as.data.frame(rf)
rownames(rf)<-name[-c(1:2)]
colnames(x2)<-name
pdf("output/rf.pdf",width = 10,height = 8)
dotchart(rf[order(rf$MeanDecreaseGini),][(length(rf$MeanDecreaseGini)- 
                                            min(30,length(rf$MeanDecreaseGini))):length(rf$MeanDecreaseGini)],
         labels = name[-c(1:2)][order(rf$MeanDecreaseGini)][
           (length(rf$MeanDecreaseGini)- min(30,length(rf$MeanDecreaseGini))):length(rf$MeanDecreaseGini)])


dev.off()
#opls
x->x2
x2<-x2[x2$Groups!="Q",]
combn(names(table(x2$Groups)),2)->x3
#对照组作分母
#ee<-x3[,1:3]
#x3<-x3[,-(1:3)]
#ee[1:2,]<-ee[2:1,]
#ee<-as.data.frame(ee)


#x3<-cbind(x3,ee)
#标准化
op<-function(x,y){
  opls(x2[x2$Groups%in%c(x,y),-c(1:2)],
       as.factor(x2$Groups[x2$Groups%in%c(x,y)]),predI = 2,
       orthoI = 1,scaleC ="pareto",
       crossvalI = length(x2$Groups[x2$Groups%in%c(x,y)])-1)->l
  dir.create(paste(x,"vs",y))
  pdf(paste(x," vs ",y,"/",x,"vs",y," opls",".pdf",sep = ""))
  plot(l)
  dev.off()
  l%>%return
}
map2(t(combn(names(table(x2$Groups)),2)[1,]),
     t(combn(names(table(x2$Groups)),2)[2,]),
     ~op(.x ,.y))->oplsda
data <- as.data.frame(oplsda[[1]]@scoreMN)
o1 <- oplsda[[1]]@orthoScoreMN[,1]
data$o1 <- o1
sink("opls data.txt") 
oplsda[[1]]
"vip"
oplsda[[1]]@vipVn
oplsda[[1]]@scoreMN
oplsda[[1]]@orthoScoreMN
sink()

# oplsda plsda
combn(names(table(x2$Groups)),2)->x3
# 导出vip值
map(oplsda,getVipVn)%>% reduce(cbind)->vip
vip<-data.frame(vip,row.names = NULL)
# 各种检验

#去QC
x2<-x2[x2$Groups!="Q",]
x2[,-c(1:2)]->dataa
#t.test
#ee<-which(names(x2)=='Docosane')
#x2<-x2[,-ee]
2:(ncol(x2)-1)%>%map(function(i){(x2[,-c(1)] %>% split(.$Groups) %>%
                                    purrr::transpose() %>%
                                    .[[i]] %>%
                                    (function(x) combn(names(x), m=2, function(y) 
                                      t.test(flatten_dbl(x[y[1]]), flatten_dbl(x[y[2]]))[["p.value"]]))%>%
                                    unlist)})%>%
  data.frame%>%t->o.t
#wilcox.test
2:(ncol(x2)-1)%>%map(function(i){(x2[,-c(1)] %>% split(.$Groups) %>%
                                    purrr::transpose() %>%
                                    .[[i]] %>%
                                    (function(x) combn(names(x), m=2, function(y) 
                                      wilcox.test(flatten_dbl(x[y[1]]), flatten_dbl(x[y[2]]))[["p.value"]]))%>%
                                    unlist)})%>%
  data.frame%>%t->w.t
# kw
x=3
kw.models<-x2[,-c(1,2)]%>%
  map(~ kruskal.test(.x~ x2$Groups ))%>%
  lapply(function(x){x[["p.value"]]})
kw.models<-as.numeric(kw.models)
#aov
aov.models =x2[,-c(1,2)]%>%
  map(~ aov(.x~ x2$Groups ))%>%map(~summary(.x)[[1]][["Pr(>F)"]][1])
aov.models<-as.numeric(aov.models)

#fc
2:(ncol(x2)-1)%>%map(function(i){(x2[,-c(1)] %>% split(.$Groups) %>%
                                    purrr::transpose() %>%
                                    .[[i]] %>%
                                    (function(x) combn(names(x), m=2, function(y) 
                                      (2^mean(c(flatten_dbl(x[y[2]])))-1)/
                                        (2^mean(c(flatten_dbl(x[y[1]])))-1)))%>%
                                    unlist)})%>%
  data.frame%>%t->fc
#2:(ncol(x2)-1)%>%map(function(i){(x2[,-c(1)] %>% split(.$Groups) %>%
#                                    purrr::transpose() %>%
#                                    .[[i]] %>%
#                                    (function(x) combn(names(x), m=2, function(y) 
#                                      (2^mean(c(flatten_dbl(x[y[1]])))-1)/
#                                        (2^mean(c(flatten_dbl(x[y[2]])))-1)))%>%
#                                    unlist)})%>%
#  data.frame%>%t->fc1
#fc<-fc[,-(4:6)]
#fc1<-fc1[,-(1:3)]
#fc<-cbind(fc,fc1)

# ok


# 合并 给名字

# cbind everyone
end<-"FC"
colnames(fc)<- paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],end)
end<-"t.test"
colnames(o.t)<- paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],end)
end<-"wilcox.test"
colnames(w.t)<- paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],end)
end<-"oplsds vip"
colnames(vip)<- paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],end)
end<-"t.test.pavl"
apply(o.t, 2,function(p) {
  p.adjust(p,method = "fdr")#fdr <-优先 or BH
  
})->o.t.pavl

colnames(o.t.pavl)<- paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],end)
end<-"wilcox.test.pavl"
apply(w.t, 2,function(p) {
  p.adjust(p,method =  "fdr")
})->w.t.pavl
colnames(w.t.pavl)<- paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],end)




end<-"oplsds vip"

rf$MeanDecreaseGini->rf

data.table(colnames(x2)[-c(1:2)],t(dataa),fc,vip,o.t,o.t.pavl,w.t,w.t.pavl,kw.models,aov.models,rf)%>%as.data.table->data
colnames(data)[1]<-"name"
colnames(data)[2:(nrow(x2)+1)]<-x2$name


name[-c(1:2)]->rownames(data)
rbind(paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"FC")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"t.test")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"w.test")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"t.test.pavl")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"w.test.pavl")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"oplsds vip")
)->names
DT::datatable(data)




writexl::write_xlsx(data , "statedata.xlsx")








#```{r }
case_character_type <- function(p, fc, vip, robots = TRUE){
  case_when(
    vip > 1 & p <0.05 &fc>1    ~ "up",
    vip > 1 & p <0.05 &fc<1    ~ "down",
    
    TRUE  ~ "other"
  )
}
i=1#第几个
#y2<-1
y2=1#FC阈值
vip2=1
#vip2<-1
p2<-0.05
# p2<-0.05
case_character_type <- function(p, fc, vip, robots = TRUE){
  case_when(
    vip > vip2 & p <p2 &fc>y2 ~ "up",
    vip > vip2 & p <p2 &fc<y2   ~ "down",
    
    TRUE  ~ "other"
  )
}



library(ggrepel)
volcano <- function(df=data, x=names[1,i],y=names[2,i],vip=names[6,i],f=y2) {
  ggplot(df,aes(x = log2(.data[[x]]),-log10(.data[[y]]) ,size=.data[[vip]],col=col,label =name )) + 
    geom_point() + 
    labs(x="log2(FC)",y="-log10(p)")+
    ggrepel::geom_text_repel(data=subset(df,col!="other"),
                             segment.size  = 0.2,
                             segment.color = "grey50",
                             show.legend = FALSE 
                             
    )  +
    theme(plot.background = element_rect(fill="white"),
          panel.background = element_rect(fill='white', colour='gray'),          
          strip.text.x=element_text(size=rel(1.2), family="serif", angle=-90),
          strip.text.y=element_text(size=rel(1.2),  family="serif") ,
          axis.ticks.length.y = unit(.25, "cm"),
          axis.ticks.length.x = unit(.25, "cm"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(face="bold", color = "black", size=20, hjust=0.5, vjust=1),
          axis.text.y = element_text(face="bold", color = "black", size=20, hjust=0, vjust=0.5),
          axis.title=element_text(size=14,face="bold"),
          axis.line.x=element_line( color='black',size = 0.75, lineend = 'square'),
          axis.line.y=element_line( color='black',size = 0.75),
          axis.ticks.length=unit(.4,"lines"),#                
          axis.ticks = element_line( color='black',size = 1.5),
          legend.title = element_text(colour="black", size=16, face="bold"),
          legend.text = element_text(color="black", size = 16),
          panel.border = element_rect(colour = "black",size = 5,fill=NA)
    )+
    guides(col = guide_legend(override.aes = list(size=5)),title.theme = element_text(
      size = 15 ))+ 
    geom_vline(xintercept=c(log2(f),-log2(f)),lty=1,col="black",lwd=1) +#添加横线
    geom_hline(yintercept =-log10(0.05),lty=1,col="black",lwd=1)+#添加竖线
    scale_color_manual(values =c("#619cff","grey","#f8766d")) 
}
i=2
for (i in 1:ncol(names)) {
  data$col<-"no"
  
  data %>%
    mutate(col = case_character_type(.data[[names[2,i]]] ,  .data[[names[1,i]]] 
                                     ,.data[[names[6,i]]] 
                                     , robots = FALSE)) ->data
  
  volcano(data, names[1,i],names[2,i],names[6,i],y2)->volcanoplot
  paste(str_remove(names[1,i]," FC"),"/volcanoplot.pdf"
        ,sep = "")->save
  pdf(save,width = 10,height = 8)
  volcanoplot%>%plot
  dev.off()  
  #
}
#corplot
#heatmap
rbind(paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"FC")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"t.test")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"w.test")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"t.test.pavl")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"w.test.pavl")
      ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"oplsds vip")
)->names

i=1#第几个
#y2<-1

case_character_type <- function(p, fc, vip, robots = TRUE){
  case_when(
    vip > vip2 & p <p2 &fc>y2 ~ "up",
    vip > vip2 & p <p2 &fc<y2   ~ "down",
    
    TRUE  ~ "other"
  )
}

library(dplyr)

#ee<-which(names(x2)=='L-Serine')
#x2<-x2[,-ee]

for (i in 1:ncol(names)) {
  data$col<-"other"
  print(i)
  data %>%
    mutate(col = case_character_type(.data[[names[2,i]]] ,  .data[[names[1,i]]] ,
                                     .data[[names[6,i]]] , robots = FALSE)) ->data
  if({{data$col!="other"}%>%sum}!=0){
    data[data$col!="other",]%>%as.data.frame()->data2
    data2$name->name
    x2$name[c(x2$Groups%in%
                c(t(combn(names(table(x2$Groups)),2))[i,1],
                  t(combn(names(table(x2$Groups)),2))[i,2]))]->group
    data2[,c(1,which(colnames(data2)%in%group))]->x3
    data2[,c(1,which(colnames(data2)%in%group),
             which(colnames(data2)%in%names[,i])
    )]->x4
    x3<-x3 %>%
      group_by(name) %>%
      filter(row_number() == 1) %>%
      ungroup()
    x4<-x4 %>%
      group_by(name) %>%
      filter(row_number() == 1) %>%
      ungroup()
    paste(str_remove(names[1,i]," FC"),"/data.xlsx"
          ,sep = "")->save
    x4%>%writexl::write_xlsx(save)
    x3<-as.data.frame(x3)
    rownames(x3)<-x3$name               
    corr.test(t(x3[,-c(1)]), method = "spearman",adjust= "none")->cor
    cor$r->corr
    cor$p->corp
    corp->corp2
    corp[corp2>0.05]<-""
    corp[corp2<0.05]<-"*"
    corp[corp2<0.01]<-"***"
    corp<-as.matrix(corp)
    col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    library(ComplexHeatmap)
    
    ComplexHeatmap::Heatmap(as.matrix(corr), name = "correlation", 
                            col = col_fun, rect_gp = gpar(type = "none"),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              grid.rect(x = x, y = y, width = width, height = height, 
                                        gp = gpar(col = "grey", fill = NA))
                              grid.circle(x = x, y = y, r = abs(as.matrix(corr)[i, j])*
                                            min(unit.c(width, height))/2, 
                                          gp = gpar(fill = col_fun(as.matrix(corr)[i, j]), col = NA))
                              grid.text(corp[i, j], x, y-unit(.004,"npc"), gp = gpar(fontsize = 10),just="center")
                              
                            }
                            
    )->het
    paste(str_remove(names[1,i]," FC"),"/corrplot.pdf"
          ,sep = "")->save
    pdf(save,width = 10,height = 8)
    draw(het)
    dev.off()
    #heatmap
    rbind(paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"FC")
          ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"t.test")
          ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"w.test")
          ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"t.test.pavl")
          ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"w.test.pavl")
          ,paste(t(combn(names(table(x2$Groups)),2))[,1],"vs",t(combn(names(table(x2$Groups)),2))[,2],"oplsds vip")
    )->names
    #x3<-x3[-1,]
    x3[,-c(1)]%>%t%>%scale()%>%t->x3[,-c(1)]
    
    rownames(x3)<-x3$name               
    col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    Heatmap(((x3[,-c(1)])%>%as.matrix()), name = "scale data", col = col_fun,
            column_split = x2$Groups[x2$Groups%in%
                                       c(t(combn(names(table(x2$Groups)),2))[i,1],
                                         t(combn(names(table(x2$Groups)),2))[i,2])],
            top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:3))),
            
            
    )->het
    paste(str_remove(names[1,i]," FC"),"/corr heatmap.pdf"
          ,sep = "")->save
    pdf(save,width = 10,height = 8)
    draw(het)
    dev.off()
  }
}



i=12

dat <- function(i){data$col<-"other"
data %>%
  mutate(col = case_character_type(.data[[names[2,i]]] ,  .data[[names[1,i]]] ,
                                   .data[[names[6,i]]] , robots = FALSE)) ->data
if({{data$col!="other"}%>%sum}!=0){
  data[data$col!="other",]%>%as.data.frame()->data2
  data2$name->name
  x2$name[c(x2$Groups%in%
              c(t(combn(names(table(x2$Groups)),2))[1,1],
                t(combn(names(table(x2$Groups)),2))[1,2]))]->group
  data2[,c(1,which(colnames(data2)%in%group))]->x3
  data2[,c(1,
           which(colnames(data2)%in%names[1,i])
  )]->x4
  return(x4)
  print(i)
}
print(i)

return((data[,  which(colnames(data2)%in%names[1,i])
]))

}

map(1:ncol(names),~dat(.x))->l

l%>% reduce(full_join, by = "name")->y
y%>%writexl::write_xlsx("diff.xlsx")

















dat <- function(i){data$col<-"other"
data %>%
  mutate(col = case_character_type(.data[[names[2,i]]] ,  .data[[names[1,i]]] ,
                                   .data[[names[6,i]]] , robots = FALSE)) ->data
if({{data$col!="other"}%>%sum}!=0){
  data[data$col!="other",]%>%as.data.frame()->data2
  data2$name->name
  x2$name[c(x2$Groups%in%
              c(t(combn(names(table(x2$Groups)),2))[1,1],
                t(combn(names(table(x2$Groups)),2))[1,2]))]->group
  data2[,c(1,which(colnames(data2)%in%group))]->x3
  data2[,c(1,1,
           which(colnames(data2)%in%names[1,i])
  )]->x4
  colnames(x4)[2]<-names[1,i]%>%str_remove_all(" FC")
  return(x4)
  print(i)
}
print(i)

return(0)

}

map(1:ncol(names),~dat(.x))->l

l%>% reduce(full_join, by = "name")->y
y%>%writexl::write_xlsx("diffname.xlsx")
x2[,-c(1:3)]=2^x2[,-c(1:3)]-1
write.csv(x2,"x2.csv")
