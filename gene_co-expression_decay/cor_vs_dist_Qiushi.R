##main script examine the gene co-expression decays with distance among genes, with pathway genes color coded##

#####START HERE#####

options(stringsAsFactors = FALSE)
setwd("/Users/Qiushi/Dropbox/Mac/Documents/croseus/plot/genome_expression/gene_co-expression_decay_croseus")
library(plyr)
library(ggplot2)


################# generate plotting data for distances 10 mb or less ##############

#read in pathway allocations, make sure header is path1 path2 dist cor compare
path10less <- read.table("path_10mbless_convert.txt", sep="\t", header=T)

#read in cor and distance data
decay10less <- read.table("mergelist_nodup_noself_sortdist_10mbless.txt", sep="\t", header=T)

#prep background data
rng1 <- c(0,10,100,1000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000)
rng2 <- seq(600000,10000000, by=200000)
rng <- append(rng1, rng2)
labels <- rng[2:66]

newdata<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q05=quantile(value.y, 0.05))
newdata1<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q25=quantile(value.y, 0.25))
newdata2<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q5=quantile(value.y, 0.50))
newdata3<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q75=quantile(value.y, 0.75))
newdata4<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q95=quantile(value.y, 0.95))
newdata5<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,mean=mean(value.y))

a <- Reduce(function(x, y) merge(x, y, all=TRUE), list(newdata, newdata1, newdata2, newdata3, newdata4, newdata5))
colnames(a)[1] <- "bins"
###sort a
a <- a[order(a$bins),]

#prep pathway data
same <- path10less[path10less$compare=="same",]
diff <- path10less[path10less$compare=="different",]

df <- data.frame(x = cut(same$dist, breaks=rng, labels=labels,include.lowest = TRUE), y= same$cor, z=same$path1)
df1 <- data.frame(x = cut(diff$dist, breaks=rng, labels=labels,include.lowest = TRUE), y= diff$cor)


################# generate plotting data for distances 10 mb or more ##############

#read in pathway allocations
path10more <- read.table("path_10mbmore_convert.txt", sep="\t", header=T)

#read in cor and distance data
decay10more <- read.table("mergelist_nodup_noself_sortdist_10mbmore.txt", sep="\t", header=T)

#prep background data
rng0 <- seq(10000000,78000000, by=1400000)
labels0 <- rng0[2:49]

newdata0<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q05=quantile(value.y, 0.05))
newdata10<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q25=quantile(value.y, 0.25))
newdata20<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q5=quantile(value.y, 0.50))
newdata30<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q75=quantile(value.y, 0.75))
newdata40<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q95=quantile(value.y, 0.95))
newdata50<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,mean=mean(value.y))
a0 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(newdata0, newdata10, newdata20, newdata30, newdata40, newdata50))
colnames(a0)[1] <- "bins"
a0 <- a0[order(a0$bins),] 

#prep pathway data
same0 <- path10more[path10more$compare=="same",]
diff0 <- path10more[path10more$compare=="different",]

df0 <- data.frame(x = cut(same0$dist, breaks=rng0, labels=labels0,include.lowest = TRUE), y= same0$cor,z=same0$path1)
df10 <- data.frame(x = cut(diff0$dist, breaks=rng0, labels=labels0,include.lowest = TRUE), y= diff0$cor)


########################################### merge 10mb or less with 10mb or more #########

aboth <- rbind(a,a0)
sameboth <- rbind(same,same0)
diffboth <- rbind(diff,diff0)
labelsboth <- append(labels,labels0)

#################################    plot it all ####



#pdf("10mb.col.rev.final1.pdf",h=7,w=10)
ggplot(aboth, aes(x = as.numeric(log10(labelsboth)))) +
### plot from 100bp dist
#ggplot(aboth[-1,], aes(x = as.numeric(log10(labelsboth[-1])))) +
  geom_point(data=diffboth, aes(x=as.numeric(log10(diffboth$dist)), y=diffboth$cor, color=factor(diffboth$compare), fill= factor(diffboth$compare)), size=2, shape=21, alpha = 0.75)+
  geom_point(data=sameboth, aes(x=as.numeric(log10(sameboth$dist)), y=sameboth$cor, color=factor(sameboth$path1), fill= factor(sameboth$path1)), shape=21, size=2, alpha = 0.75)+ 
  scale_color_manual(values = c('Aerial alkaloid'= "black",'Iridoid'=rgb(0.5,0,0), 'lipid'="red", 'MEP'="purple", 'Root alkaloid'="brown", 'terpene'=rgb(0,0.5,1), 'different'="black")) +
  scale_fill_manual(values = c('Aerial alkaloid'= "black",'Iridoid'=rgb(0.5,0,0), 'lipid'="red", 'MEP'="purple", 'Root alkaloid'="brown", 'terpene'=rgb(0,0.5,1),'different'="white"))+
  geom_ribbon(aes(x = as.numeric(log10(labelsboth)),ymax=q95, ymin=q75), alpha  = .3)+
  geom_ribbon(aes(x = as.numeric(log10(labelsboth)),ymax=q75, ymin=q25), alpha  = .2)+
  geom_ribbon(aes(x = as.numeric(log10(labelsboth)),ymax=q25, ymin=q05), alpha  = .3)+
  geom_line(aes(x = as.numeric(log10(labelsboth)),y=aboth$mean),size=0.25)+
  scale_x_continuous(name="Intergenic distance log10[D]", expand = c(0.005, 0))+
  scale_y_continuous(name=expression("Pairwise gene expression correlation" (rho^2)), expand = c(0, 0.07), breaks=c(0,0.25,0.5,0.75,1))+
  geom_vline(xintercept=log10(100000), linetype="dashed", color = "black", size=0.25)+
  geom_vline(xintercept=log10(200000), linetype="dashed", color = "black", size=0.25)+
  geom_vline(xintercept=log10(500000), linetype="dashed", color = "black", size=0.25)+
  geom_vline(xintercept=log10(1000000), linetype="dashed", color = "black", size=0.25)+
  geom_vline(xintercept=log10(10000000), linetype="dashed", color = "black", size=0.25)+
  annotate("text", label="0.1Mb", x=log10(80000), y=1.1, angle=90, size=3)+
  annotate("text", label="0.2Mb", x=log10(160000), y=1.1, angle=90, size=3)+
  annotate("text", label="0.5Mb", x=log10(400000), y=1.1, angle=90, size=3)+
  annotate("text", label="1.0Mb", x=log10(800000), y=1.1, angle=90, size=3)+
  annotate("text", label="10.0Mb", x=log10(8000000), y=1.1, angle=90, size=3)+
  #annotate("text", label="<10Mb", x=log10(4500000), y=0.95, size=3)+
  #annotate("text", label=">10Mb", x=log10(130000000), y=0.95, size=3)+
  annotate("text", label="", x=log10(130000000), y=0.9, size=3)+
  #geom_rect(xmin = log10(10000000), xmax = log10(297000000), ymin = -0.85, ymax = 1.09, fill=NA, color="grey50", size=0.25)+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), panel.background = element_blank(), text = element_text(size= 22))
