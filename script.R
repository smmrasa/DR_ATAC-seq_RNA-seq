#Figure 1A (PCA plot and the calculations for Figure 1B)
options(stringsAsFactors = FALSE)
setwd("./RNA_seq")
ss<-read.csv("SampleSheet.csv")
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(factoextra)
library(grid)
library(gridExtra)
#drawing PCA plot and pearson correlation for each tissue
pdf("./Plot/PCAplot_table_final2_NormalizedFilteredCounts.pdf")
for(t in unique(ss$Tissue)){
  if(t=="SI"){next}
  f<-ss[ss$Tissue==t & ss$FileCategory=="Count",]
  s<-ss[ss$Tissue==t & ss$FileCategory=="SampleSheet",]
  s<-read.csv(s$FileName)
  count<-read.csv(f$FileName,row.names = 1)
  count<-count[apply(count > 10, 1, sum)>=3,]
  #PCA
  iqr<-summary(apply(log2(count+0.1), 1, IQR))
  iq<-iqr[4]
  countPCA <- count[apply(log2(count+0.1), 1, IQR) >= iq, ]
  t.log<-log(countPCA+0.1)
  df<-t.log[ , ! apply( t.log , 2 , function(x) all(is.na(x)) ) ]
  pca<-prcomp(t(df),center = TRUE,scale. = TRUE)
  
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- c(paste( "PC1(", as.character(percentage)[1] ,"%)", sep=""),paste( "PC2(", as.character(percentage)[2] ,"%)", sep="")) 
  p1<-cbind(as.data.frame(pca$x[,1:2]),Group=s$Group)
  rownames(p1)<-s$SampleName
  
  #table for comparing group based on the PCA
  young<-grep("Young",p1$Group)
  LTAL<-grep("LTAL",p1$Group)
  LTDR<-grep("LTDR",p1$Group)
  STDR<-grep("STDR",p1$Group)
  young.pc1<-mean(p1$PC1[young])
  young.pc2<-mean(p1$PC2[young])
  
  LTAL.pc1<-mean(p1$PC1[LTAL])
  LTAL.pc2<-mean(p1$PC2[LTAL])
  
  LTDR.pc1<-mean(p1$PC1[LTDR])
  LTDR.pc2<-mean(p1$PC2[LTDR])
  
  STDR.pc1<-mean(p1$PC1[STDR])
  STDR.pc2<-mean(p1$PC2[STDR])
  
  #young as point 0
  pc1.Aging<-LTAL.pc1-young.pc1
  pc2.Aging<-LTAL.pc2-young.pc2
  
  pc1.DRl<-LTDR.pc1-young.pc1
  pc2.DRl<-LTDR.pc2-young.pc2
  
  pc1.DRs<-STDR.pc1-young.pc1
  pc2.DRs<-STDR.pc2-young.pc2
  
  #caculating the distance from young
  
  l.Age.y=sqrt((pc1.Aging^2)+(pc2.Aging^2))
  l.DRs.y<-sqrt((pc1.DRs^2)+(pc2.DRs^2))
  l.DRl.y<-sqrt((pc1.DRl^2)+(pc2.DRl^2))
  
  a.age<-pc2.Aging/pc1.Aging
  at.age<-atan(a.age)
  
  a.DRs<-pc2.DRs/pc1.DRs
  at.DRs<-atan(a.DRs)
  
  a.DRl<-pc2.DRl/pc1.DRl
  at.DRl<-atan(a.DRl)
  #  
  alpha.DRl<-at.DRl-at.age
  alpha.DRs<-at.DRs-at.age
  
  #The height of the perpendicular line
  x.s<-l.DRs.y*sin(alpha.DRs)
  x.l<-l.DRl.y*sin(alpha.DRl)
  
  x.s.norm<-x.s/l.Age.y
  x.l.norm<-x.l/l.Age.y
  #the distance between young and the prependicular line
  ytoP.DRs<-sqrt(l.DRs.y^2-x.s^2)
  ytoP.DRl<-sqrt(l.DRl.y^2-x.l^2)
  
  r.DRs<-ytoP.DRs/l.Age.y
  r.DRl<-ytoP.DRl/l.Age.y
  
  #distance between DRLT and DRST
  LT.ST<-sqrt((LTDR.pc1-STDR.pc1)^2+(LTDR.pc2-STDR.pc2)^2)
  
  #plot
  #adding sudo points to have circular ellipse
  s<-rbind(s[,c("FileName","SampleName","Group","GroupPriority")],data.frame(FileName=NA,SampleName=c("Young_sudo1","Young_sudo2","LTAL_sudo1","LTAL_sudo2",
                                                                                                      "LTDR_sudo1","LTDR_sudo2","STDR_sudo1","STDR_sudo2"),
                                                                             Group=c(paste(t,c("Young","Young","LTAL","LTAL","LTDR","LTDR","STDR","STDR"),sep="_")),GroupPriority=NA))
  
  pca<-prcomp(t(df),center = TRUE,scale. = TRUE)
  
  sudo=5
  pca$x<-rbind(pca$x,c(young.pc1+sudo,young.pc2-sudo,rep(0,ncol(pca$x)-2)))
  pca$x<-rbind(pca$x,c(young.pc1-sudo,young.pc2+sudo,rep(0,ncol(pca$x)-2)))
  pca$x<-rbind(pca$x,c(LTAL.pc1+sudo,LTAL.pc2-sudo,rep(0,ncol(pca$x)-2)))
  pca$x<-rbind(pca$x,c(LTAL.pc1-sudo,LTAL.pc2+sudo,rep(0,ncol(pca$x)-2)))
  pca$x<-rbind(pca$x,c(LTDR.pc1+sudo,LTDR.pc2-sudo,rep(0,ncol(pca$x)-2)))
  pca$x<-rbind(pca$x,c(LTDR.pc1-sudo,LTDR.pc2+sudo,rep(0,ncol(pca$x)-2)))
  pca$x<-rbind(pca$x,c(STDR.pc1+sudo,STDR.pc2-sudo,rep(0,ncol(pca$x)-2)))
  pca$x<-rbind(pca$x,c(STDR.pc1-sudo,STDR.pc2+sudo,rep(0,ncol(pca$x)-2)))
  rownames(pca$x)[(nrow(pca$x)-7):nrow(pca$x)]<-c("Young_sudo1","Young_sudo2","LTAL_sudo1","LTAL_sudo2",
                                                  "LTDR_sudo1","LTDR_sudo2","STDR_sudo1","STDR_sudo2")
  
  p<-fviz(pca, element="ind", axes=c(1,2), geom=c("point"),
          habillage = s$Group, repel = TRUE, palette = "Dark2",
          addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.999,mean.point=FALSE,alpha=0) + ggtitle(t)
  ggp<-ggplot_build(p)  
  xr<-max(ggp$layout$panel_params[[1]]$x.range)-min(ggp$layout$panel_params[[1]]$x.range)
  yr<-max(ggp$layout$panel_params[[1]]$y.range)-min(ggp$layout$panel_params[[1]]$y.range)
  del<-abs(xr-yr)
  p<-p+
    coord_fixed(ratio = 1, xlim = NULL, ylim = c(min(ggp$layout$panel_params[[1]]$y.range)-(del/2),max(ggp$layout$panel_params[[1]]$y.range)+(del/2)), expand = TRUE, clip = "on")
  
  p<-p+
    geom_text(aes(x=young.pc1,y=young.pc2,label="Young",col=unique(p1$Group[young])))+
    geom_text(aes(x=LTAL.pc1,y=LTAL.pc2,label="LTAL",col=unique(p1$Group[LTAL])))+
    geom_text(aes(x=LTDR.pc1,y=LTDR.pc2,label="LTDR",col=unique(p1$Group[LTDR])))+
    geom_text(aes(x=STDR.pc1,y=STDR.pc2,label="STDR",col=unique(p1$Group[STDR])))
  
  p<-p+geom_point(aes(x=pca$x[,1],y=pca$x[,2],col=s$Group ) ,alpha=c(rep(1,nrow(p1)),rep(0,8)))
  print(p)
  
  tab<-data.frame(Young_Old=round(l.Age.y,2),
                  DRLT.AgingRelated=paste0(round(ytoP.DRl),"(",round((1-r.DRl)*100),"%)"),DRST.AgingRelated=paste0(round(ytoP.DRs),"(",round((1-r.DRs)*100),"%)"),
                  LTDR.NotAgingRelated=round(abs(x.l)), STDR.NotAgingRelated=round(abs(x.s)))
  grid.newpage()
  grid.table(tab,rows=NULL,theme=ttheme_default(base_size = 5))
  
  
  
}
dev.off()







#################################################################################
##################################################################################
#figure 3C S2A and S2B
#plotting the final network (with all of the tissues)

library("igraph")
#function to calculate the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

options(stringsAsFactors = FALSE)
dir<-read.table("Directed connection_DR_inflamation_pathway.txt",skip=1,sep="\t",header=TRUE)
res<-read.csv("GSEAResults2_56UR.csv")
res$UpstreamRegulator<-gsub("LPS","lipopolysaccharide",res$UpstreamRegulator)
res$UpstreamRegulator<-gsub("IFH1","IFIH1",res$UpstreamRegulator)
res$UpstreamRegulator<-gsub("IFN beta","IFN Beta",res$UpstreamRegulator)

res<-res[res$UpstreamRegulator %in% c(dir$From.Molecule.s.,dir$To.Molecule.s.), ]
res$Tissue<-factor(res$Tissue)
net.LT<-data.frame()
net.ST<-data.frame()
net.Aging<-data.frame()

for(u in unique(res$UpstreamRegulator)){
  tmp<-res[res$UpstreamRegulator==u,] 
  tmp<-tmp[tmp$pval<0.05,]  
  #LTDR
  tmp2<-tmp[tmp$Group=="LTAL_vs_LTDR",]
  d<-getmode(sign(tmp2$NES))
  tmp2<-tmp2[sign(tmp2$NES)==d,]
  tmp2$pval<-(-log10(tmp2$pval))
  s<-sum(tmp2$pval)
  m<-matrix(rep(0,length(levels(res$Tissue))),nrow = 1)
  colnames(m)<-levels(res$Tissue)
  p<-match(tmp2$Tissue,colnames(m))
  m[1,p]<-tmp2$pval
  net.LT<-rbind(net.LT,cbind(UpstreamRegulator=u,p.value.sum=s,Direction=d,as.data.frame(m)))
  
  #STDR
  tmp2<-tmp[tmp$Group=="LTAL_vs_STDR",]
  d<-getmode(sign(tmp2$NES))
  tmp2<-tmp2[sign(tmp2$NES)==d,]
  tmp2$pval<-(-log10(tmp2$pval))
  s<-sum(tmp2$pval)
  m<-matrix(rep(0,length(levels(res$Tissue))),nrow = 1)
  colnames(m)<-levels(res$Tissue)
  p<-match(tmp2$Tissue,colnames(m))
  m[1,p]<-tmp2$pval
  net.ST<-rbind(net.ST,cbind(UpstreamRegulator=u,p.value.sum=s,Direction=d,as.data.frame(m)))
  
  #Aging
  tmp2<-tmp[tmp$Group=="Young_vs_LTAL",]
  d<-getmode(sign(tmp2$NES))
  tmp2<-tmp2[sign(tmp2$NES)==d,]
  tmp2$pval<-(-log10(tmp2$pval))
  s<-sum(tmp2$pval)
  m<-matrix(rep(0,length(levels(res$Tissue))),nrow = 1)
  colnames(m)<-levels(res$Tissue)
  p<-match(tmp2$Tissue,colnames(m))
  m[1,p]<-tmp2$pval
  net.Aging<-rbind(net.Aging,cbind(UpstreamRegulator=u,p.value.sum=s,Direction=d,as.data.frame(m)))
  
}


#Drawing network
#LTDR
coor<-read.csv("coordinates.csv")
coor<-t(coor)
coor<-coor[nrow(coor):1,]
co<-t(sapply(net.LT$UpstreamRegulator,function(x) which(coor==x,arr.ind = TRUE)))
colnames(co)<-c("y","x")
co<-co[,c(2,1)]

net.LT$Direction<-factor(net.LT$Direction)
levels(net.LT$Direction)<-c("Down","Up")
nodes<-data.frame(id=net.LT$UpstreamRegulator,name=paste0(net.LT$UpstreamRegulator,"(",net.LT$Direction,")"),size=net.LT$p.value.sum,co)
tmp<-dir[dir$From.Molecule.s. %in% nodes$id & dir$To.Molecule.s. %in% nodes$id,]
tmp<-tmp[tmp$From.Molecule.s.!=tmp$To.Molecule.s.,]
tmp<-tmp[!duplicated(paste(tmp$From.Molecule.s.,tmp$To.Molecule.s.)),]
links<-data.frame(from=tmp$From.Molecule.s.,to=tmp$To.Molecule.s.,type=tmp$Relationship.Type)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 

V(net)$size<-V(net)$size
E(net)$arrow.size<-0.5
E(net)$arrow.width<-0.5
E(net)$width<-1
#E(net)$color<-"black"
V(net)$label <- V(net)$name
V(net)$label.cex<-0.7
V(net)$label.dist<-2
V(net)$label.degree<-(-pi/2)

#l <- layout_with_fr(net)
#plot(net, layout=l,main="")
library(RColorBrewer)

#layout.reingold.tilford(net)
pie.val=sapply(1:nrow(net.LT),function(x) list(as.numeric(net.LT[x,4:11])))

pdf("LTDR_network_rotate.pdf")
g<-plot(net, vertex.shape="pie", vertex.pie=pie.val,
        vertex.pie.color=list(brewer.pal(8, "Set1")),main="DRLT")

legend("topleft",legend = levels(res$Tissue),fill=brewer.pal(8, "Set1"),cex=0.8)



dev.off()

write.csv(net.LT,"Net.LT.csv")
#scales

nodes<-data.frame(id=seq(5,30,by = 5),name=seq(5,30,by = 5),size=seq(5,30,by = 5))
net2 <- graph_from_data_frame(d=links[0,],vertices=nodes, directed=T) 

V(net2)$size<-V(net2)$size
#E(net)$color<-"black"
V(net2)$label <- V(net2)$name
V(net2)$label.cex<-0.7
V(net2)$label.dist<-2
V(net2)$label.degree<-(-pi/2)

#l <- layout_with_fr(net)
pdf("ScaleForNetwork.pdf")
plot(net2)
dev.off()
#layout.reingold.tilford(net)
pie.val=sapply(1:nrow(net.LT),function(x) list(as.numeric(net.LT[x,4:11])))

pdf("LTDR_network_rotate.pdf")
g<-plot(net, vertex.shape="pie", vertex.pie=pie.val,
        vertex.pie.color=list(brewer.pal(8, "Set1")),main="DRLT")

legend("topleft",legend = levels(res$Tissue),fill=brewer.pal(8, "Set1"),cex=0.8)



dev.off()


#STDR
# coor<-read.csv("coordinates.csv")
# co<-t(sapply(net.ST$UpstreamRegulator,function(x) which(coor==x,arr.ind = TRUE)))
# colnames(co)<-c("y","x")
# co<-co[,c(2,1)]

net.ST$Direction<-factor(net.ST$Direction)
levels(net.ST$Direction)<-c("Down","Up")
nodes<-data.frame(id=net.ST$UpstreamRegulator,name=paste0(net.ST$UpstreamRegulator,"(",net.ST$Direction,")"),size=net.ST$p.value.sum,co)
tmp<-dir[dir$From.Molecule.s. %in% nodes$id & dir$To.Molecule.s. %in% nodes$id,]
tmp<-tmp[tmp$From.Molecule.s.!=tmp$To.Molecule.s.,]
tmp<-tmp[!duplicated(paste(tmp$From.Molecule.s.,tmp$To.Molecule.s.)),]
links<-data.frame(from=tmp$From.Molecule.s.,to=tmp$To.Molecule.s.,type=tmp$Relationship.Type)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 

V(net)$size<-V(net)$size
E(net)$arrow.size<-0.5
E(net)$arrow.width<-0.5
E(net)$width<-1
#E(net)$color<-"black"
V(net)$label <- V(net)$name
V(net)$label.cex<-0.7
V(net)$label.dist<-2
V(net)$label.degree<-(-pi/2)

#l <- layout_with_fr(net)
#plot(net, layout=l,main="")
library(RColorBrewer)

#layout.reingold.tilford(net)
pie.val=sapply(1:nrow(net.ST),function(x) list(as.numeric(net.ST[x,4:11])))

pdf("STDR_network_rotate.pdf")
plot(net, vertex.shape="pie", vertex.pie=pie.val,
     vertex.pie.color=list(brewer.pal(8, "Set1")),main="DRST")

legend("topleft",legend = levels(res$Tissue),fill=brewer.pal(8, "Set1"),cex=0.8)
dev.off()

write.csv(net.ST,"Net.ST.csv")


#Aging
# coor<-read.csv("coordinates.csv")
# co<-t(sapply(net.Aging$UpstreamRegulator,function(x) which(coor==x,arr.ind = TRUE)))
# colnames(co)<-c("y","x")
# co<-co[,c(2,1)]

net.Aging$Direction<-factor(net.Aging$Direction)
levels(net.Aging$Direction)<-c("Down","Up")
nodes<-data.frame(id=net.Aging$UpstreamRegulator,name=paste0(net.Aging$UpstreamRegulator,"(",net.Aging$Direction,")"),size=net.Aging$p.value.sum,co)
tmp<-dir[dir$From.Molecule.s. %in% nodes$id & dir$To.Molecule.s. %in% nodes$id,]
tmp<-tmp[tmp$From.Molecule.s.!=tmp$To.Molecule.s.,]
tmp<-tmp[!duplicated(paste(tmp$From.Molecule.s.,tmp$To.Molecule.s.)),]
links<-data.frame(from=tmp$From.Molecule.s.,to=tmp$To.Molecule.s.,type=tmp$Relationship.Type)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 

V(net)$size<-V(net)$size
E(net)$arrow.size<-0.5
E(net)$arrow.width<-0.5
E(net)$width<-1
#E(net)$color<-"black"
V(net)$label <- V(net)$name
V(net)$label.cex<-0.7
V(net)$label.dist<-2
V(net)$label.degree<-(-pi/2)

#l <- layout_with_fr(net)
#plot(net, layout=l,main="")
library(RColorBrewer)

#layout.reingold.tilford(net)
pie.val=sapply(1:nrow(net.Aging),function(x) list(as.numeric(net.Aging[x,4:11])))

pdf("Aging_network_rotate.pdf")
plot(net, vertex.shape="pie", vertex.pie=pie.val,
     vertex.pie.color=list(brewer.pal(8, "Set1")),main="Aging")

legend("topleft",legend = levels(res$Tissue),fill=brewer.pal(8, "Set1"),cex=0.8)
dev.off()

write.csv(net.Aging,"Net.Aging.csv")


##################################################################################
###############################################
#Figure 5A, PCA plot ATAC-seq
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(factoextra)
library(grid)
library(gridExtra)
setwd("ATAC_seq")
pdf("PCA_PerasonCor_allTissue.pdf")
#kidney  
count<-read.csv("KidneyAllSamples_RPM.csv",row.names = 1)

t<-"Kidney"
rownames(count)<-count$position
count<-count[,-1]
colnames(count)<-gsub("DRLT","LTDR",colnames(count))
sam<-gsub("_RPM","" ,colnames(count))
gro<-sapply(sam, function(x) strsplit(x,"_")[[1]][2])
s<-data.frame(Sample=sam,Group=gro)
s$Group<-factor(s$Group,levels=c("Young","LTAL","LTDR"))
iqr<-summary(apply(log2(count+0.1), 1, IQR))
iq<-iqr[4]
countPCA <- count[apply(log2(count+0.1), 1, IQR) >= iq, ]
t.log<-log(countPCA+0.1)
df<-t.log[ , ! apply( t.log , 2 , function(x) all(is.na(x)) ) ]
pca<-prcomp(t(df),center = TRUE,scale. = TRUE)

percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- c(paste( "PC1(", as.character(percentage)[1] ,"%)", sep=""),paste( "PC2(", as.character(percentage)[2] ,"%)", sep="")) 
p1<-cbind(as.data.frame(pca$x[,1:2]),Group=s$Group)
rownames(p1)<-s$SampleName

#plot
#adding sudo points to have circular ellipse
s<-rbind(s,data.frame(Sample=c("Young_sudo1","Young_sudo2","LTAL_sudo1","LTAL_sudo2",
                               "LTDR_sudo1","LTDR_sudo2"),
                      Group=c(c("Young","Young","LTAL","LTAL","LTDR","LTDR"))))

pca<-prcomp(t(df),center = TRUE,scale. = TRUE)

sudo=5
pca$x<-rbind(pca$x,c(young.pc1+sudo,young.pc2-sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(young.pc1-sudo,young.pc2+sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTAL.pc1+sudo,LTAL.pc2-sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTAL.pc1-sudo,LTAL.pc2+sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTDR.pc1+sudo,LTDR.pc2-sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTDR.pc1-sudo,LTDR.pc2+sudo,rep(0,ncol(pca$x)-2)))
rownames(pca$x)[(nrow(pca$x)-5):nrow(pca$x)]<-c("Young_sudo1","Young_sudo2","LTAL_sudo1","LTAL_sudo2",
                                                "LTDR_sudo1","LTDR_sudo2")

p<-fviz(pca, element="ind", axes=c(1,2), geom=c("point"),
        habillage = s$Group, repel = TRUE, palette = "Dark2",
        addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.999,mean.point=FALSE,alpha=0) + ggtitle(t)
ggp<-ggplot_build(p)  
xr<-max(ggp$layout$panel_params[[1]]$x.range)-min(ggp$layout$panel_params[[1]]$x.range)
yr<-max(ggp$layout$panel_params[[1]]$y.range)-min(ggp$layout$panel_params[[1]]$y.range)
del<-abs(xr-yr)
p<-p+
  coord_fixed(ratio = 1, xlim = NULL, ylim = c(min(ggp$layout$panel_params[[1]]$y.range)-(del/2),max(ggp$layout$panel_params[[1]]$y.range)+(del/2)), expand = TRUE, clip = "on")

p<-p+
  geom_text(aes(x=young.pc1,y=young.pc2,label="Young",col=unique(p1$Group[young])))+
  geom_text(aes(x=LTAL.pc1,y=LTAL.pc2,label="LTAL",col=unique(p1$Group[LTAL])))+
  geom_text(aes(x=LTDR.pc1,y=LTDR.pc2,label="LTDR",col=unique(p1$Group[LTDR])))
p<-p+geom_point(aes(x=pca$x[,1],y=pca$x[,2],col=s$Group ) ,alpha=c(rep(1,nrow(p1)),rep(0,6)))
print(p)


#Liver  
count<-read.csv("LiverAllSamples_withrepRemove_RPM.csv",row.names = 1)

t<-"Liver"
rownames(count)<-count$position
count<-count[,-1]
colnames(count)<-gsub("DRLT","LTDR",colnames(count))
sam<-gsub("_RPM","" ,colnames(count))
gro<-sapply(sam, function(x) strsplit(x,"_")[[1]][2])
s<-data.frame(Sample=sam,Group=gro)
s$Group<-factor(s$Group,levels=c("Young","LTAL","LTDR"))
iqr<-summary(apply(log2(count+0.1), 1, IQR))
iq<-iqr[4]
countPCA <- count[apply(log2(count+0.1), 1, IQR) >= iq, ]
t.log<-log(countPCA+0.1)
df<-t.log[ , ! apply( t.log , 2 , function(x) all(is.na(x)) ) ]
pca<-prcomp(t(df),center = TRUE,scale. = TRUE)

percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- c(paste( "PC1(", as.character(percentage)[1] ,"%)", sep=""),paste( "PC2(", as.character(percentage)[2] ,"%)", sep="")) 
p1<-cbind(as.data.frame(pca$x[,1:2]),Group=s$Group)
rownames(p1)<-s$SampleName



#plot
#adding sudo points to have circular ellipse
s<-rbind(s,data.frame(Sample=c("Young_sudo1","Young_sudo2","LTAL_sudo1","LTAL_sudo2",
                               "LTDR_sudo1","LTDR_sudo2"),
                      Group=c(c("Young","Young","LTAL","LTAL","LTDR","LTDR"))))

pca<-prcomp(t(df),center = TRUE,scale. = TRUE)

sudo=5
pca$x<-rbind(pca$x,c(young.pc1+sudo,young.pc2-sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(young.pc1-sudo,young.pc2+sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTAL.pc1+sudo,LTAL.pc2-sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTAL.pc1-sudo,LTAL.pc2+sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTDR.pc1+sudo,LTDR.pc2-sudo,rep(0,ncol(pca$x)-2)))
pca$x<-rbind(pca$x,c(LTDR.pc1-sudo,LTDR.pc2+sudo,rep(0,ncol(pca$x)-2)))
rownames(pca$x)[(nrow(pca$x)-5):nrow(pca$x)]<-c("Young_sudo1","Young_sudo2","LTAL_sudo1","LTAL_sudo2",
                                                "LTDR_sudo1","LTDR_sudo2")

p<-fviz(pca, element="ind", axes=c(1,2), geom=c("point"),
        habillage = s$Group, repel = TRUE, palette = "Dark2",
        addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.999,mean.point=FALSE,alpha=0) + ggtitle(t)
ggp<-ggplot_build(p)  
xr<-max(ggp$layout$panel_params[[1]]$x.range)-min(ggp$layout$panel_params[[1]]$x.range)
yr<-max(ggp$layout$panel_params[[1]]$y.range)-min(ggp$layout$panel_params[[1]]$y.range)
del<-abs(xr-yr)
p<-p+
  coord_fixed(ratio = 1, xlim = NULL, ylim = c(min(ggp$layout$panel_params[[1]]$y.range)-(del/2),max(ggp$layout$panel_params[[1]]$y.range)+(del/2)), expand = TRUE, clip = "on")

p<-p+
  geom_text(aes(x=young.pc1,y=young.pc2,label="Young",col=unique(p1$Group[young])))+
  geom_text(aes(x=LTAL.pc1,y=LTAL.pc2,label="LTAL",col=unique(p1$Group[LTAL])))+
  geom_text(aes(x=LTDR.pc1,y=LTDR.pc2,label="LTDR",col=unique(p1$Group[LTDR])))
p<-p+geom_point(aes(x=pca$x[,1],y=pca$x[,2],col=s$Group ) ,alpha=c(rep(1,nrow(p1)),rep(0,6)))
print(p)


dev.off()

######################################################################################
####################################################################################
#Figure 5D, 5E, 5F, and 5G

library(methylKit)
library(genomation)
library(ggplot2)
options(stringsAsFactors = FALSE)
t<-"Kidney"


#promoters (longest isoforms)
res<-read.csv(paste0(t,"_promoter_refflat_longestisoforman_rpkm.csv"))
res<-res[apply(res[,7:ncol(res)],1,function(x) sum(x>0))>2 & apply(res[,7:ncol(res)],1,IQR)>0,]

colnames(res)<-gsub("LTDR","DRLT",colnames(res))
young<-grep("Young",colnames(res))
old<-grep("LTAL",colnames(res))
DR<-grep("DRLT",colnames(res))
res2<-data.frame(gene=res$gene,Group="Young",mean=apply(res[,young],1,mean))
res2<-rbind(res2,data.frame(gene=res$gene,Group="ALLT",mean=apply(res[,old],1,mean)))
res2<-rbind(res2,data.frame(gene=res$gene,Group="DRLT",mean=apply(res[,DR],1,mean)))
res2$Group<-factor(res2$Group,levels=c("Young","ALLT","DRLT"))


res3<-res2
res3$mean<-res3$mean+1
res3$mean<-log10(res3$mean)

library(mixtools)

norm<-read.csv("../RNA seq/With_removing_replicates/FPKM_calculatedfromHtseqCount.csv")
rownames(norm)<-norm$gene

pdf(paste0("Promoter_hist_modeling_3",t,".pdf"))
for (gr in c("Young", "ALLT",  "DRLT" )){
  res4<-res3[res3$Group==gr,]
  
  wait1 <- normalmixEM(res4$mean)
  par(xpd = FALSE)
  plot(wait1, breaks = 100,density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8, main2=gr, xlab2="Log10 FPKM",xlim=c(0,1))
  plot(wait1, breaks = 100,density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8, main2=gr, xlab2="Log10 FPKM",xlim=c(0,1))
  abline(v = wait1$mu[1]-wait1$sigma[1],lwd=2, lty=2)
  abline(v = wait1$mu[2]-wait1$sigma[2],lwd=2, lty=2)
  
  abline(v = wait1$mu[1]+wait1$sigma[1],lwd=2, lty=2)
  abline(v = wait1$mu[2]+wait1$sigma[2],lwd=2, lty=2)
  
  #ploting the expression of each normal distribution
  
  d1<-res4[res4$mean>=(wait1$mu[1]-wait1$sigma[1]) & res4$mean<=(wait1$mu[1]+wait1$sigma[1]),]
  d2<-res4[res4$mean>=(wait1$mu[2]-wait1$sigma[2]) & res4$mean<=(wait1$mu[2]+wait1$sigma[2]),]
  
  
  exp<-norm[,grep(t,colnames(norm))]
  
  colnames(exp)<-gsub("LTDR","DRLT",colnames(exp))
  
  young<-grep("Young",colnames(exp))
  old<-grep("LTAL",colnames(exp))
  DR<-grep("DRLT",colnames(exp))
  exp.mean<-data.frame(Young=apply(exp[,young],1,function(x) mean(as.numeric(x))),ALLT=apply(exp[,old],1,function(x) mean(as.numeric(x))),DRLT=apply(exp[,DR],1,function(x) mean(as.numeric(x))))
  
  exp.1<-exp.mean[rownames(exp.mean) %in% d1$gene,]
  exp.2<-exp.mean[rownames(exp.mean) %in% d2$gene,]
  
  pl<-data.frame(Group=gr,FPKM=exp.1[,gr],GeneGroup="Peak1")
  pl<-rbind(pl,data.frame(Group=gr,FPKM=exp.2[,gr],GeneGroup="Peak2"))
  pval=wilcox.test(as.numeric(pl[pl$GeneGroup=="Peak1","FPKM"]),as.numeric(pl[pl$GeneGroup=="Peak2","FPKM"]))$p.value
  #pl$FPKM<-log10(pl$FPKM+1)
  pl$GeneGroup<-factor(pl$GeneGroup,levels=c("Peak1","Peak2"))
  p<-ggplot(pl,aes(x=GeneGroup,y=FPKM))+
    geom_boxplot(outlier.shape = NA)+
    coord_cartesian(ylim =c(0,55))+
    ggtitle(paste("pvalue=",pval))
  print(p)
  
  
  
}
dev.off()
###########################################################################################
########################################################################################
#figure S4B S4C S4D S4F
library(methylKit)
library(genomation)
library(ggplot2)
options(stringsAsFactors = FALSE)
t<-"Liver"


#promoters (longest isoforms)


res<-read.csv(paste0(t,"_promoter_refflat_longestisoforman_rpkm.csv"))
res<-res[,! colnames(res) %in% c( "Liver_Young_rep3","Liver_LTDR_rep2" )]
res<-res[apply(res[,7:ncol(res)],1,function(x) sum(x>0))>2 & apply(res[,7:ncol(res)],1,IQR)>0,]

colnames(res)<-gsub("LTDR","DRLT",colnames(res))
young<-grep("Young",colnames(res))
old<-grep("LTAL",colnames(res))
DR<-grep("DRLT",colnames(res))
res2<-data.frame(gene=res$gene,Group="Young",mean=apply(res[,young],1,mean))
res2<-rbind(res2,data.frame(gene=res$gene,Group="ALLT",mean=apply(res[,old],1,mean)))
res2<-rbind(res2,data.frame(gene=res$gene,Group="DRLT",mean=apply(res[,DR],1,mean)))
res2$Group<-factor(res2$Group,levels=c("Young","ALLT","DRLT"))


res3<-res2
res3$mean<-res3$mean+1
res3$mean<-log10(res3$mean)

library(mixtools)

norm<-read.csv("../RNA_seq/FPKM_calculatedfromHtseqCount.csv")
rownames(norm)<-norm$gene

pdf(paste0("Promoter_hist_modeling_3",t,".pdf"))
for (gr in c("Young", "ALLT",  "DRLT" )){
  res4<-res3[res3$Group==gr,]
  
  wait1 <- normalmixEM(res4$mean)
  par(xpd = FALSE)
  plot(wait1, breaks = 100,density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8, main2=gr, xlab2="Log10 FPKM",xlim=c(0,1))
  plot(wait1, breaks = 100,density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8, main2=gr, xlab2="Log10 FPKM",xlim=c(0,1))
  abline(v = wait1$mu[1]-wait1$sigma[1],lwd=2, lty=2)
  abline(v = wait1$mu[2]-wait1$sigma[2],lwd=2, lty=2)
  
  abline(v = wait1$mu[1]+wait1$sigma[1],lwd=2, lty=2)
  abline(v = wait1$mu[2]+wait1$sigma[2],lwd=2, lty=2)
  
  #ploting the expression of each normal distribution
  
  d1<-res4[res4$mean>=(wait1$mu[1]-wait1$sigma[1]) & res4$mean<=(wait1$mu[1]+wait1$sigma[1]),]
  d2<-res4[res4$mean>=(wait1$mu[2]-wait1$sigma[2]) & res4$mean<=(wait1$mu[2]+wait1$sigma[2]),]
  
  
  exp<-norm[,grep(t,colnames(norm))]
  
  colnames(exp)<-gsub("LTDR","DRLT",colnames(exp))
  
  young<-grep("Young",colnames(exp))
  old<-grep("LTAL",colnames(exp))
  DR<-grep("DRLT",colnames(exp))
  exp.mean<-data.frame(Young=apply(exp[,young],1,function(x) mean(as.numeric(x))),ALLT=apply(exp[,old],1,function(x) mean(as.numeric(x))),DRLT=apply(exp[,DR],1,function(x) mean(as.numeric(x))))
  
  exp.1<-exp.mean[rownames(exp.mean) %in% d1$gene,]
  exp.2<-exp.mean[rownames(exp.mean) %in% d2$gene,]
  
  pl<-data.frame(Group=gr,FPKM=exp.1[,gr],GeneGroup="Peak1")
  pl<-rbind(pl,data.frame(Group=gr,FPKM=exp.2[,gr],GeneGroup="Peak2"))
  pval=wilcox.test(as.numeric(pl[pl$GeneGroup=="Peak1","FPKM"]),as.numeric(pl[pl$GeneGroup=="Peak2","FPKM"]))$p.value
  #pl$FPKM<-log10(pl$FPKM+1)
  pl$GeneGroup<-factor(pl$GeneGroup,levels=c("Peak1","Peak2"))
  p<-ggplot(pl,aes(x=GeneGroup,y=FPKM))+
    geom_boxplot(outlier.shape = NA)+
    coord_cartesian(ylim =c(0,55))+
    ggtitle(paste("pvalue=",pval))
  print(p)
  
  
  
}
dev.off()










