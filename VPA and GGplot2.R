## VPA and GGplot2
#prepare sample-species and sample-factors table.txt
library(vegan)
library(ggplot2)
setwd('')
sp <-read.table(file=file.choose(),sep="\t",header=T,row.names=1)
sp
se <-read.table(file=file.choose(),sep="\t",header=T,row.names=1)
se
##Analysis strategy selection: >4 CCA Unimodal model, 3~4 CCA and RDA, <3 RDA Linear model
decorana(sp)
sp0 <- rda(sp ~ 1, se)  
sp0
plot(sp0)
C.whole <- rda(sp ~ ., se)  
C.whole
plot(C.whole)  ##-----this is over of cca and rda, if vpa {##step vpa}
##Extract coordinate elements
new<-C.whole$CCA
new
samples<-data.frame(sample=row.names(new$u),RDA1=new$u[,1],RDA2=new$u[,2])
samples
species<-data.frame(spece=row.names(new$v),RDA1=new$v[,1],RDA2=new$v[,2])
species
envi<-data.frame(en=row.names(new$biplot),RDA1=new$biplot[,1],RDA2=new$biplot[,2])
envi
line_x =c(0,envi[1,2],0,envi[2,2],0,envi[3,2],0,envi[4,2],0,envi[5,2],0,envi[6,2])
line_x
line_y =c(0,envi[1,3],0,envi[2,3],0,envi[3,3],0,envi[4,3],0,envi[5,3],0,envi[6,3])
line_y
line_g =c("pH","pH","T","T","S2","S2","NH4","NH4","NO2","NO2","Fe2","Fe2") ##This is your double factors
line_g
line_data =data.frame(x=line_x,y=line_y,group=line_g)
line_data
p=ggplot(data=samples,aes(RDA1,RDA2)) 
##onestep graph
p1=p +geom_point(aes(color=sample),size=2) +
  
  geom_text(data=envi,aes(label=en),color="blue") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  geom_line(data=line_data,aes(x=x,y=y,group=group),color="green")+
  theme_bw() + theme(panel.grid=element_blank())
p1

##-----this is vpa based on cca and rda-----##
x.sig = anova(C.whole)
x.p = x.sig$Pr[1]
plot(C.whole,dis=c('wa','cn'),main = paste( "(p=", x.p,")"))
inf_factor = vif.cca(C.whole)
max_env = which(inf_factor == max(inf_factor))
env.st3 = env.st2
max_env = which(inf_factor == max(inf_factor))
env.st3 = env.st2
while ( inf_factor[max_env] > 20){		
  env.st3 = env.st3[,-max_env]		
  C.reduced = cca(fg.dat2, env.st3)		
  inf_factor = vif.cca(C.reduced)		
  max_env = which(inf_factor == max(inf_factor))		
}
inf_factor
colnames(env.st3)
ind.p = array(0,dim=c(1,ncol(env.st2)))
ind.F = array(0,dim=c(1,ncol(env.st2)))
for(j in 1:ncol(env.st2)){		
  ind.cca = cca(fg.dat2, env.st2[,j]) #ind.cca = cca(fg.dat, env.st[,j], env.st[,-j])  #		
  ind.sig = anova(ind.cca,step=1000)		
  ind.p[1,j] = ind.sig$Pr[1]		
  ind.F[1,j] = ind.sig$F[1]		
}
colnames(ind.p) = colnames(env.st2)
t(rbind(ind.F,ind.p))
###two group,
#pH,Ce,H2Rec,RT
env.grp1 = c("pH","Ce")      
env.grp2 = c("H2Rec","RT")  

#pCCA
env.selected = cbind( env.st2[ ,env.grp1], env.st2[,env.grp2]) ;env.selected
C.grp1.par  = cca(otu, env.selected[,env.grp1], env.selected[,c(env.grp2)])
C.grp2.par  = cca(otu, env.selected[,env.grp2], env.selected[,c(env.grp1)])
C.grp1_2.par = C.whole

str(C.whole)
total.chi = C.whole$tot.chi            ;total.chi         
total.constrained = C.whole$CCA$tot.chi;total.constrained 
str(C.grp1.par)
grp1.par.chi  = C.grp1.par$CCA$tot.chi ;grp1.par.chi
grp2.par.chi  = C.grp2.par$CCA$tot.chi ;grp2.par.chi
grp1_2.chi = C.grp1_2.par$CCA$tot.chi  ;grp1_2.chi

overlap.grp1_2.chi = grp1_2.chi - (grp1.par.chi + grp2.par.chi) ;overlap.grp1_2.chi

grp1.percent = grp1.par.chi / total.chi  ;grp1.percent
grp2.percent = grp2.par.chi / total.chi  ;grp2.percent
grp1_2.percent = overlap.grp1_2.chi/total.chi; if(grp1_2.percent < 0 ){grp1_2.percent = 0} ;grp1_2.percent

unexplained.percent = (total.chi - total.constrained) / total.chi

grp1.percent
grp2.percent
grp1_2.percent
unexplained.percent


###three group
# for pCCA
env.grp1 = c("pH","Ce")      
env.grp2 = "H2Rec" 
env.grp3 = "RT"    

env.selected = cbind(env.st2[,env.grp1], env.st2[,env.grp2], env.st2[,env.grp3])
head(env.selected)
colnames(env.selected)[3] = "H2Rec"
colnames(env.selected)[4] = "RT" 
head(env.selected)


C.grp1.par  = cca(otu, env.selected[,env.grp1], env.selected[,c(env.grp2,env.grp3)])
C.grp2.par  = cca(otu, env.selected[,env.grp2], env.selected[,c(env.grp1,env.grp3)])
C.grp3.par  = cca(otu, env.selected[,env.grp3], env.selected[,c(env.grp1,env.grp2)])
C.grp1_2.par = cca(otu, env.selected[,c(env.grp1,env.grp2)], env.selected[,env.grp3])
C.grp1_3.par = cca(otu, env.selected[,c(env.grp1,env.grp3)], env.selected[,env.grp2])
C.grp2_3.par = cca(otu, env.selected[,c(env.grp2,env.grp3)], env.selected[,env.grp1])
C.grp1_2_3.par = C.whole

str(C.whole)

total.chi = C.whole$tot.chi ; total.chi
total.constrained = C.whole$CCA$tot.chi ; total.constrained
grp1.par.chi  = C.grp1.par$CCA$tot.chi  ; grp1.par.chi
grp2.par.chi  = C.grp2.par$CCA$tot.chi  ; grp2.par.chi
grp3.par.chi  = C.grp3.par$CCA$tot.chi  ; grp3.par.chi
grp1_2.chi = C.grp1_2.par$CCA$tot.chi   ; grp1_2.chi 
grp1_3.chi = C.grp1_3.par$CCA$tot.chi   ; grp1_3.chi
grp2_3.chi = C.grp2_3.par$CCA$tot.chi   ; grp2_3.chi

overlap.grp1_2.chi = grp1_2.chi - (grp1.par.chi + grp2.par.chi)
overlap.grp1_3.chi = grp1_3.chi - (grp1.par.chi + grp3.par.chi)
overlap.grp2_3.chi = grp2_3.chi - (grp2.par.chi + grp3.par.chi)
overlap.all = total.constrained-(grp1.par.chi+grp2.par.chi+grp3.par.chi+overlap.grp1_2.chi+overlap.grp1_3.chi+overlap.grp2_3.chi)

grp1.percent = grp1.par.chi / total.chi
grp2.percent = grp2.par.chi / total.chi
grp3.percent = grp3.par.chi / total.chi
grp1_2.percent = overlap.grp1_2.chi/total.chi; if(grp1_2.percent < 0 ){grp1_2.percent = 0}
grp1_3.percent = overlap.grp1_3.chi/total.chi; if(grp1_3.percent < 0 ){grp1_3.percent = 0}
grp2_3.percent = overlap.grp2_3.chi/total.chi; if(grp2_3.percent < 0 ){grp2_3.percent = 0}
overlap.percent = overlap.all/total.chi;if(overlap.percent < 0){overlap.percent = 0}

unexplained.percent = (total.chi - total.constrained) / total.chi

grp1.percent
grp2.percent
grp3.percent
grp1_2.percent
grp1_3.percent
grp2_3.percent
overlap.percent
unexplained.percent