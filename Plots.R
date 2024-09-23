######################################################
########## Plotting script for OR paper   ############
########## Kimia Vahdat: kvahdat@ncsu.edu ############
##########         Date: 9/2/2022        ############
######################################################

# Load ggplot2
library(ggplot2)
library(MASS)
library(dplyr)
library(reshape2)
library(reshape)


### Color Palette

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### reading the data ###
setwd("/Users/kimia/Library/CloudStorage/OneDrive-NorthCarolinaStateUniversity/Google Drive/Kimia/writing presentation/Dissertation")

### features boxplots ####
test2=read.csv("./MLCI2.csv")
sapply(test2, class)

factor(test2$Method)
# 
l=factor(c("Opt HOIF FIB CV","HOIF FIB CV","HOIF FIB" ,"Barton","LamQian", "LOOBoot", "RepCV"))
l
l=relevel(l,ref = "Opt HOIF FIB CV")
#l=relevel(l, ref = "Linear")
l
# test2$DGF=factor(test2$DGF,levels = l)

lx=factor(c("Linear","Poly","Complex"))
lx=relevel(lx, ref = "Poly")
lx=relevel(lx, ref = "Linear")
lx

factor(test2$Method)

neworder <- c("Linear","Poly","Complex")
library(plyr)  ## or dplyr (transform -> mutate)
test22 <- arrange(transform(test2,
                           DGF=factor(DGF,levels=neworder)),DGF)

test23 <- arrange(transform(test22,
                            Noise=factor(Noise,levels=c("low","high"))),Noise)

png("./MLCIPlus.png", units="cm",width=90,height=60,res=350 )
p<-ggplot(test23 , aes(x=DGF, y=Method, group=Noise)) +
  scale_x_discrete(breaks=as.character(lx))+
  geom_errorbarh(aes(xmin=lci,xmax=uci), height=0.4, size=3)+
  
  scale_y_discrete(breaks=as.character(l),limits=rev(as.character(l)))+
  
  theme_light()+ #theme( axis.text.x = element_blank())+#,
                   #     axis.text.y = element_blank())+
  theme( panel.background = element_rect(fill = "transparent", color=NA),
         panel.border = element_blank(),
         panel.grid.minor = element_blank(),
         # Change axis line
         axis.line = element_line(colour = "black",size = 3))
#print(p)
p=p+facet_wrap(Noise~DGF,nrow = 2, scales="free",
               strip.position = "top"#,
               #labeller = as_labeller(c(EC="EC",ER="ER",SC="SC",SR="SR", LM="", SVM="")) 
               )+
  theme(legend.text = element_blank())+
  theme(strip.text = element_blank(),strip.background = element_blank(),strip.placement = "outside", 
        legend.text = element_blank(), legend.title = element_blank())+
  
  
  theme( axis.title.x = element_blank() ,  axis.title.y = element_blank()) 

p + theme(panel.spacing = unit(1, "lines"), legend.title = element_blank(), legend.text = element_text())
#print(p)

dev.off()


#############################

library(forcats)
# setwd("C:/Users/Administrator/OneDrive - North Carolina State University/Google Drive/Kimia/results/OPTE results/plots")
test=read.csv("./SimCI.csv")
trcost=data.frame(matrix(c("Sc1","Sc2","Sc3","Sc4",175.01,188.5,191.4,193.7), nrow = 4, ncol = 2))
colnames(trcost)=c("Sc","Exp")
trcost$Exp=as.numeric(trcost$Exp)

png("./SimApp10.png", units="cm",width=40,height=20,res=350 )
p<-ggplot(test 
          , aes(x=Sc, y=avg,group=Method))+ 
  geom_point(aes(color=Method)) + 
  #geom_line(aes(y=NLCI, color=Type), linetype="dashed")+
  #geom_line( aes(y=NUCI, color=Type), linetype="dashed")+
  geom_smooth(aes(ymin=lci,ymax=uci,color=Method,fill =Method),
              stat = "identity")+
  geom_line( aes(y=avg,color=Method), size=1.5)+
  #geom_abline(slope=0,intercept = 0,aes(color="red"), size=1)+ ylim(c(-1,1))+
  #scale_color_manual(breaks=c("SOFS-GAFS","SOFS-RFE"),values=c(cbPalette[c(2,3,2,3)])) + 
  theme( panel.background = element_rect(fill = "transparent",colour = NA),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         # Change axis line
         axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c(cbPalette[c(2,3,8)])) +
  scale_colour_manual(values=c(cbPalette[c(2,3,8)]))+
  #labs(x="Sample Size",y= "Normalized Performance") 
  theme(legend.text = element_text(size=30),
        legend.title = element_text(size=30,face="bold"),
        axis.text = element_text(size=20,face="bold"))
p= p +
geom_line(aes(y=c(trcost$Exp,trcost$Exp,trcost$Exp)), size=1.5)
#print(p)
# p+facet_wrap(~Dataset, scales="free",strip.position = "top" )+theme(strip.background = element_rect(),
#                                                                     strip.placement = "outside", strip.text = element_text(size = 30,face="bold"))
p
dev.off()












#############################

library(forcats)
# setwd("C:/Users/Administrator/OneDrive - North Carolina State University/Google Drive/Kimia/results/OPTE results/plots")
test=read.csv("./SimCI_bias.csv")
trcost=data.frame(matrix(c("Sc1","Sc2","Sc3","Sc4",182.2,178.7,181.95,199.26), nrow = 4, ncol = 2))
colnames(trcost)=c("Sc","Exp")

trcost$Exp=as.numeric(trcost$Exp)
#trcost$n= rep("10",4)
#trcost=rbind(trcost,trcost)
#trcost$n[5:8]=rep("50",4)
test$n=as.character(test$n)
test$Method=as.factor(test$Method)
lx=factor(c("Crude","Bias-corrected"))
# lx=relevel(lx, ref = "IU-inflated")
lx=relevel(lx, ref = "Crude")
lx
colnames(test)
png("./SimCI_1123.png", units="cm",width=60,height=20,res=350 )
p<-ggplot(test 
          , aes(x=Sc, y=mean,group=Method))+ 
  geom_point(aes(color=Method)) + 
  #geom_line(aes(y=NLCI, color=Type), linetype="dashed")+
  #geom_line( aes(y=NUCI, color=Type), linetype="dashed")+
  geom_smooth(aes(ymin=lci,ymax=uci,color=Method,fill =Method),
              stat = "identity")+ ylim(c(160,220))+
  scale_fill_manual(breaks=c("Crude","Bias-corrected"),values=c(cbPalette[c(2,6,2,6)])) +
  scale_color_manual(breaks=c("Crude","Bias-corrected"),values=c(cbPalette[c(2,6,2,6)])) +
  geom_line( aes(y=mean,color=Method), size=1.5)+
  # 
  # scale_y_discrete(breaks=as.character(lx),limits=rev(as.character(lx)))+
  #geom_abline(slope=0,intercept = 0,aes(color="red"), size=1)+ ylim(c(-1,1))+
  #scale_color_manual(breaks=c("Crude Output Analysis","With Bias Correction"),values=c(cbPalette[c(2,8,2,8)])) + 
  theme( panel.background = element_rect(fill = "transparent",colour = NA),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         # Change axis line
         axis.line = element_line(colour = "black"))+
  
  #scale_colour_manual(values=c(cbPalette[c(2,8)]))+
  #labs(x="Sample Size",y= "Normalized Performance") 
  theme(legend.text = element_text(size=30),
        legend.title = element_text(size=30,face="bold"),
        axis.text = element_text(size=20,face="bold"))

 p= p +
   geom_line(aes(y=c(rep(trcost$Exp[c(1,1,2,2,3,3,4,4)],2))), size=1.5)
p= p+facet_wrap(~n, scales="free",strip.position = "top" )+theme(strip.background = element_rect(),
                                                                     strip.placement = "outside", strip.text = element_text(size = 30,face="bold"))
p
dev.off()

##### corrupt
library(forcats)
# setwd("C:/Users/Administrator/OneDrive - North Carolina State University/Google Drive/Kimia/results/OPTE results/plots")
test=read.csv("./SimCI_bias_corrupt.csv")
trcost=data.frame(matrix(c("Sc1","Sc2","Sc3","Sc4",183,176,180,201), nrow = 4, ncol = 2))
colnames(trcost)=c("Sc","Exp")

trcost$Exp=as.numeric(trcost$Exp)
#trcost$n= rep("10",4)
#trcost=rbind(trcost,trcost)
#trcost$n[5:8]=rep("50",4)
test$n=as.character(test$n)
test$Method=as.factor(test$Method)
lx=factor(c("Crude","Bias-corrected"))
# lx=relevel(lx, ref = "IU-inflated")
lx=relevel(lx, ref = "Crude")
lx
colnames(test)
png("./SimCI_corrupt_1123.png", units="cm",width=60,height=20,res=350 )
p<-ggplot(test 
          , aes(x=Sc, y=mean,group=Method))+ 
  geom_point(aes(color=Method)) + 
  #geom_line(aes(y=NLCI, color=Type), linetype="dashed")+
  #geom_line( aes(y=NUCI, color=Type), linetype="dashed")+
  geom_smooth(aes(ymin=lci,ymax=uci,color=Method,fill =Method),
              stat = "identity")+ ylim(c(160,220))+
  scale_fill_manual(breaks=c("Crude","Bias-corrected"),values=c(cbPalette[c(2,6,2,6)])) +
  scale_color_manual(breaks=c("Crude","Bias-corrected"),values=c(cbPalette[c(2,6,2,6)])) +
  geom_line( aes(y=mean,color=Method), size=1.5)+
  # 
  # scale_y_discrete(breaks=as.character(lx),limits=rev(as.character(lx)))+
  #geom_abline(slope=0,intercept = 0,aes(color="red"), size=1)+ ylim(c(-1,1))+
  #scale_color_manual(breaks=c("Crude Output Analysis","With Bias Correction"),values=c(cbPalette[c(2,8,2,8)])) + 
  theme( panel.background = element_rect(fill = "transparent",colour = NA),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         # Change axis line
         axis.line = element_line(colour = "black"))+
  
  #scale_colour_manual(values=c(cbPalette[c(2,8)]))+
  #labs(x="Sample Size",y= "Normalized Performance") 
  theme(legend.text = element_text(size=30),
        legend.title = element_text(size=30,face="bold"),
        axis.text = element_text(size=20,face="bold"))

p= p +
  geom_line(aes(y=c(rep(trcost$Exp[c(1,1,2,2,3,3,4,4)],2))), size=1.5)
p= p+facet_wrap(~n, scales="free",strip.position = "top" )+theme(strip.background = element_rect(),
                                                                 strip.placement = "outside", strip.text = element_text(size = 30,face="bold"))
p
dev.off()



