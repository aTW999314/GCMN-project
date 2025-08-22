
#differential
library(limma)

setwd("/home/Tanwei")
rt=read.table('rna-seq.txt',header=T,sep="\t",check.names=F)
row.names(rt)=rt[,1]
rt=as.matrix(rt)
row.names(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames = dimnames)
rt=avereps(rt)
str(rt)
rt1=as.data.frame(rt)
rt2=log2(rt1+1)
write.table(rt1,file="symbol.txt",sep='\t',row.names = TRUE,quote = F)


class<-c(rep("NS",8),rep("GCMN",18))
design<-model.matrix(~factor(class))
colnames(design)<-c("NS","GCMN")
str(rt1)
fit<-lmFit(rt1,design)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.table(allDiff,file="limmaTab-all.txt",sep="\t",quote=F)

#write table
diffLab<-allDiff[with(allDiff, P.Value<0.05),]
diffLab <- diffLab[abs(diffLab$logFC) > 2, ]
write.table(diffLab,file="diffEXp.txt",sep="\t",quote=F)
###??????
diffExpLevel<-rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)
write.table(diffExpQvalue,file="diffExpLevel.txt",sep="\t",quote=F)

######volcano####
rm(list = ls())
library(RColorBrewer)
library(ggplot2)
library(ggrepel) 

df=diffLab
df=read.table('diffEXp1.txt',header=T,sep="\t",check.names=F,row.names = 1)
df$threshold = factor(ifelse(df$P.Value  < 0.05 & abs(df$logFC) >= 2, ifelse(df$logFC >= 2 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) 
fix(df)
pdf(file="volplot-bulk-rna.pdf",height=12,width=9)
ggplot(df,aes(x=logFC,y= -log10(adj.P.Val),color=threshold))+
  geom_point(data = df[df$adj.P.Val<0.05&abs(df$logFC)>2,],size = 5)+ 
  geom_point(data = df[df$adj.P.Val>0.05|abs(df$logFC)<2,],size = 4)+
  scale_color_manual(values=c("#FC4E2A","#4393C3","#00000033"))+#ȷ????????ɫ
  geom_text_repel(
    data = df[df$adj.P.Val<0.05&abs(df$logFC)>2,],
    aes(label = gene),
    size = 4.5,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#???ӹ?ע?ĵ??Ļ?????
  ylab('-log10 (P.value)')+#?޸?y??????
  xlab('log2 (FoldChange)')+#?޸?x??????
  geom_vline(xintercept=c(-2,2),lty=3,col="black",lwd=0.5) +#|logFoldChange|>0.25
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5) +#padj<0.05
  theme_classic(  
    base_line_size = 1 
  )+
  theme(axis.title.x = element_text(size = 15, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 15,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black",
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13, 
                                   color = "black", 
                                   face = "bold", 
                                   vjust = 0.5,
                                   hjust = 0.5, 
                                   angle = 0),
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )

dev.off()
