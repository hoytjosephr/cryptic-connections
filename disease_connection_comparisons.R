rm(list=ls())

##READ THIS!###
#meta-data relevant to this file:
#species = which spp was sampled
###MYLU = M. lucifugus
###PESU = P. subflavus
###MYSE = M. septentrionalis
#newsite = coded site names
#pdate2 = pdate2 is the partial date (units are months) begining with the first date of sampling in November and ending with last date of sampling in March. Compare to Ldate for actual date visted.

############################ Pd invasion data ####################
 rm(list=ls())
library(igraph)
library(reshape2)
library(ggplot2)
library(scales)
library(lme4)
library(tidyverse)
library(reshape2)

#load in disease data (0/1)
dd2 = read.csv("disease_dat.csv")

#run mixed model
lr1=glmer(gd~species*pdate2 +(1|newsite),data=dd2,family="binomial");summary(lr1)#

#create dataframe for smooth lines
dd2.new=expand.grid(seq(from=min(dd$pdate2),to=max(dd$pdate2),by=.1),factor(c("MYLU","MYSE","PESU")),factor(c("MI_BC","MI_MC","WI_SJ","WI_SB","MI_GA","WI_SP","WI_ST","MI_TM")))
colnames(dd2.new)=c("pdate2","species","site")
#predict model
dd2.new$phat=predict(lr1,dd2.new,re.form=NA,type="response")

#use full names instead of abbreviations
dd2$new_spec[dd2$species=="MYSE"]="Myotis septentrionalis"
dd2$new_spec[dd2$species=="MYLU"]="Myotis lucifugus"
dd2$new_spec[dd2$species=="PESU"]="Perimyotis subflavus"

dd2.new$new_spec[dd2.new$species=="MYSE"]="Myotis septentrionalis"
dd2.new$new_spec[dd2.new$species=="MYLU"]="Myotis lucifugus"
dd2.new$new_spec[dd2.new$species=="PESU"]="Perimyotis subflavus"

#color palette
lifeisbeautiful=c("#801637","#047878","#FFB733","#F57336","#C22121","darkblue", "black", "tan")

##summarize data - plot summarized points instead of 0/1s
r1=aggregate(gd~newsite+pdate+species,FUN=mean,data=dd);r1=r1[order(r1$newsite,r1$pdate,r1$species),];r1 
r2=aggregate(c~newsite+pdate+species,FUN=sum,data=dd);r2=r2[order(r2$newsite,r2$pdate,r2$species),];r2 
r1$N=r2$c
r1

#add error bars
r1$se=(r1$gd*(1-r1$gd)/r1$N)^.5
Pmax.all.zeros<-function(N,Pr=0.05){ 1 - exp( log(Pr) / N )}
r1$se[!is.na(r1$se==0|r1$se==1)]=Pmax.all.zeros(r1$N[!is.na(r1$se==0|r1$se==1)],.05)/2
r1=subset(r1, N > 3) #remove very small sample sizes
r1$pdate2=r1$pdate-min(r1$pdate) #rescale date to match 0/1 model fit data

#use full names instead of abbreviations
r1$new_spec[r1$species=="MYSE"]="Myotis septentrionalis"
r1$new_spec[r1$species=="MYLU"]="Myotis lucifugus"
r1$new_spec[r1$species=="PESU"]="Perimyotis subflavus"


#plot predicted lines
#remove any NA species
dd2=dd2[!(is.na(dd2$new_spec)),]
r1.noepfu=r1[!(is.na(r1$new_spec)),]

p<-ggplot(dd2,aes(x=pdate2,y=gd))+
  facet_wrap(~new_spec,scales = "fixed")+
  geom_point(data=r1.noepfu,aes(x=pdate2,y=gd,color=newsite),size=4)+
  geom_errorbar(data=r1.noepfu,aes(ymin=gd-se, ymax=gd+se, color=newsite), width=.1)+
  labs(color='Site') +
  geom_line(data=dd2.new,aes(x=pdate2,y=phat),size=1)+
  ylab(expression(~italic("P. destructans")~"prevalence"))+
  xlab("Month")+
  scale_x_continuous(labels = c("Nov","Dec", "Jan", "Feb", "Mar"))+
  scale_color_manual(values=lifeisbeautiful)+
  coord_cartesian(ylim=c(-0.1,1.1))+
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),strip.text.x = element_text (face="italic",size = 15,hjust = 0.5, vjust = 0.5),axis.title=element_text(size=25),axis.text=element_text(size=15),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.text = element_text(size=14), axis.text.x=element_text(angle=0))
p

#################Pd and dust comparison#################
#read in aggregated disease and dust file
prev_all2 = read.csv("prev_all2.csv")

#below is the calculation for change in dust and fungus - already in the file, but explicitly here...
prev_all2$d_dust=(prev_all2$dust_late-prev_all2$dust_early)/(1-prev_all2$dust_early)
prev_all2$d_gd=(prev_all2$gd_late-prev_all2$gd_early)/(1-prev_all2$gd_early)

#analysis
k3=lmer(d_gd~d_dust + (1|newsite), data=prev_all2);summary(k3)
2*(1-pt(5.496,8))
#get a p-value from the t value
library(MuMIn)
r.squaredGLMM(k3)
hist(resid(k3))
shapiro.test(resid(k3))
#passes normality test

#predict model
preds=predict(k3,re.form=NA,type="response", se=T)
preds$ci = 1.96*preds$se.fit
prev_all2 = cbind(prev_all2,preds)

c <- ggplot(prev_all2, aes(x=d_dust, y=d_gd)) + 
  geom_line(aes(y = fit),col="black", size=1)+
  geom_ribbon(aes(ymin=fit-ci, ymax=fit+ci),alpha=0.1)+
  geom_jitter(aes(col= species),size=5, shape=20, width=.05,stroke=3)+
  scale_color_manual(name="", values=c("lightblue","red","orange"))+
  ylab("Change in pathogen prevalence")+
  xlab("Change in UVF-dust prevalence")+
  geom_abline(intercept = 0,slope=1,linetype="dashed")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.title=element_text(size=18),axis.text=element_text(size=15),axis.text.x=element_text(face="italic"), axis.line=element_line(),legend.position="top",legend.text = element_text(size=12,face="italic"))
c


#############Pd and shared group comparison###############
dm = read.csv("sim_comb_df6.csv")
prev_all3 = read.csv("cluster_pd_comparison.csv")

#nonlinear model
require(nlme)
e4=d_gd~1 / (1 + exp(-c1*(dcluster_size-c2)))

f3b <- nlme(e4, 
           data = prev_all3,
           fixed = c(c1 ~ 1,c2~1),
           random = c1~ 1|site,
           start = c(c1=1,c2=.1),
           control = nlmeControl(maxIter = 100))
summary(f3b)

1-var(resid(f3b))/var(prev_all2$d_gd) #~R2 for nlme model (similar to R2 for a regular nls b/c random effects are very small)
fitted(f3b)
ranef(f3b)#random effects - very small as indicated by SD of random effect in summary

#create new dataframe for smooth lines
xseq=seq(from=min(prev_all3$dcluster_size), to=max(prev_all3$dcluster_size), by=.001)
newdat=data.frame(dcluster_size=xseq, d_gd=NA)

##nonlinear
newdat$yhat=predict(f3b,newdata = newdat, level = 0)


#change spp names from abbreviations
prev_all3$nspecies[prev_all3$species=="MYLU"]="Myotis lucifugus"
prev_all3$nspecies[prev_all3$species=="MYSE"]="Myotis septentrionalis"
prev_all3$nspecies[prev_all3$species=="PESU"]="Perimyotis subflavus"

g2 <- ggplot(prev_all3, aes(x=dcluster_size, y=d_gd)) + 
  ylab("Change in pathogen prevalence")+
  xlab("Change in shared group connections")+#M. lucifugus - Social
  coord_cartesian(ylim=c(0,1),xlim=c(0,1))+
  geom_line(data=newdat,aes(x=dcluster_size,y=yhat),size = 1)+
  geom_point(aes(colour = nspecies), size=5.5)+
  scale_color_manual(name="", values=c("lightblue","red","orange"))+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.title=element_text(size=18),axis.text=element_text(size=15),axis.text.x=element_text(face="italic"), axis.line=element_line(),legend.position="top",legend.text = element_text(size=12,face="italic"))
g2

ggsave(file="/Users/klangwig/Dropbox/Contact rate MS/Figures/23AUG2018_ChangeInPdWithChangeinSharedGroup.pdf",width=5.5,height=4,units="in",dpi=300,useDingbats = F)

