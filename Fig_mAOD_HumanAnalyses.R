#   *************************************************************************
#   ***           Scripts to perform the human analyses                   ***
#   ***   from the paper Transition from background selection to          ***
#   ***  associative overdominance promotes diversity in regions of       ***
#   ***                low recombination                                  ***
#   ***   by Gilbert, Pouyet, Excoffier and Peischl                       ***
#   ***   Theses scripts were prepared by Fanny Pouyet on 2019.           ***
#   ***              Last update 05.10.2019                               ***
#   *************************************************************************

library(scales)
source("populations.R")
library(gplots)

#####Prepare data
########

where="PDF/" # where to put the figures
pop.names=c("YRI", "LWK","GBR","IBS","PJL","BEB","KHV", "JPT","CLM","PEL")
order=c(8,10,1,3,6,7,4,9,2,5)
who = "window2Mb-step0p5Mb" 
sub.who = "2Mb Windows - 500kb step"
popcol=c("orangered","orange", "royalblue","deepskyblue", "plum","hotpink", "navyblue","slateblue","limegreen","forestgreen")

# Upload data for DAFi and for PI
nbsnp = read.table("NbSNP_LR_MR0p1-0p5_perwindow.txt",h=T)

pi.low = read.table(paste("/storage/pouyet/1000G_data/strictMask/nucdiv.lowRec0p05.windows2Mb-step0p5Mb.pi",sep=""),h=T)
pi.low[,"BIN_MEAN"]<- (pi.low[,2]+pi.low[,3])/2 # get mean position of the bin
pi.medium = read.table(paste("/storage/pouyet/1000G_data/strictMask/nucdiv.RR0p1-0p5.windows2Mb-step0p5Mb.pi",sep=""),h=T)
pi.medium[,"BIN_MEAN"]<- (pi.medium[,2]+pi.medium[,3])/2 #get mean position of the bin

#Step 1 . Make sure I have enough SNPs per windows for both low and medium RR i.e. 50 SNPs (nota bene: I always have the same number of SNP for all population, see Pouyet et al 2018 for an explanation)
not.enough.SNPs = which( rowSums(is.na(pi.low[,(4+1*0:9)])) > 0 |rowSums(is.na(pi.medium[,(4+1*0:9)])) > 0 | (nbsnp[,1])<50 | (nbsnp[,2])<50  ) 
pi.medium[(not.enough.SNPs),(4:(ncol(pi.medium)-1))]<-NA
pi.low[(not.enough.SNPs),(4:(ncol(pi.low)-1))]<-NA


#Step 2. Find the outliers

# a. initialisation
outliers = list()
for(pop in pop.names){
  outliers[[pop]] = list();
  outliers[[pop]][["Pi"]]=list()
}
mm.pi=ssd.pi = as.data.frame(matrix(0,ncol=length(pop.names),nrow=22))
#b. compute outliers for each population, an outlier has to have a pi higher than the mean +2 sd
for(chr in c(1:22)){
  sub.pi.low<-pi.low[(pi.low[,1]==chr),]
  sub.pi.medium<-pi.medium[(pi.medium[,1]==chr),]
  p=0
  for(pop in pop.names){
    p=p+1
    mm.pi[chr,p] = mean( (sub.pi.low[,(3+1*p)]- sub.pi.medium[,(3+1*p)]) / sub.pi.medium[,(3+1*p)]  ,na.rm=T)
    ssd.pi[chr,p] =  sd((sub.pi.low[,(3+1*p)]- sub.pi.medium[,(3+1*p)]) / sub.pi.medium[,(3+1*p)] , na.rm=T)
    outliers[[pop]][["Pi"]][[chr]] <-  which( ((sub.pi.low[,(3+1*p)]- sub.pi.medium[,(3+1*p)]) / sub.pi.medium[,(3+1*p)]) > (mm.pi[chr,p]+2*ssd.pi[chr,p]) )
  }
}

#c. They have to be outlier for all 10 pop
nb.pop.outliers=10
outliers.all=list()
for(chr in c(1:22)){
  outliers.all[[chr]]=-100 #initialize, I don't care if it's negative
  tmp.chr=c()
  for(pop in pop.names){
    pi.outliers = outliers[[pop]][["Pi"]][[chr]]
    #concatenate windows that are outliers for  Pi per population. I will then count the number of windows that I see nb.pop.outliers times
    tmp.chr<-c(tmp.chr,names(which(table(pi.outliers)==1)))
  }
  outliers.all[[chr]]<-tmp.chr
}
#Get all AOD candidates regions ( i.e. candidate in at least nb.pop.outliers populations)
outliers.pionly =list()
for(chr in c(1:22)){
  outliers.pionly[[chr]] = as.numeric(names(which(table(outliers.all[[chr]])>=nb.pop.outliers)))
}


if(!exists("genot.1000G")) genot.1000G=read.table("/storage/pouyet/1000G_data/strictMask/1000GPAN.genot_annot.autosome.polymorphic.noCpG.strictMask.bed")
#Column name for recombination of YRI Hapmap
colrec=104
if(!exists("genot.rec"))  genot.rec<-genot.1000G[(genot.1000G[,colrec]>0 ),]  

# If I want to count how many SNPs I have per window. It's very slow.
 #  nbsnp <- pi.low[,1:3] # get chr and positions of the bins
 #  nbsnp[,4:5]<-0 # nb snps for LR and MR regions. Inititalize at 0.
 #  for(i in 1:nrow(nbsnp)){
 #    sub.tab = genot.rec[which(genot.rec[,101] == nbsnp[i,1] & genot.rec[,102] >=nbsnp[i,2] & genot.rec[,102]<=nbsnp[i,3]),]
 #    nbsnp[i,4]<-sum(sub.tab[,colrec]<=0.05)
 #    nbsnp[i,5]<-sum(sub.tab[,colrec]>=0.1 & sub.tab[,colrec]<0.5)
 # # 
 #  }
 #  colnames(nbsnp)=c("CHROM", "POSSTART","POSEND","NbSNP_LR","NbSNP_MR")
 #  write.table(nbsnp[,4:5], "NbSNP_LR_MR0p1-0p5_perwindow.txt",quote=F, col.names=T,row.names=F,sep="\t")
 # 


#Step 3. Make the plots

################################################
#############   SupItem 2A      ###############
#############    Pi scans        ###############
################################################

make.index1pdf=function(outliers.fig, main.outliers=""){
  
  pdf(paste("SupItem2A.pdf",sep=""),width=48, height=90)
  par(mar=c(5.1,4.1,4.1,4.1),mfrow=c(6,1))
  for(chr in c(1:22)){
  par(mar=c(5.1,4.1,4.1,4.1))
    print(paste("CHR",chr,sep=" "))
    #get subtable per chromosome of pi LR and pi MR, make the normalised difference
    sub.pi.low<-pi.low[(pi.low[,1]==chr),]
    sub.pi.medium<-pi.medium[(pi.medium[,1]==chr),]
    #plot pop 1 and initialize the plot
    plot(  ((sub.pi.low[,4]- sub.pi.medium[,4]) / sub.pi.medium[,4])~sub.pi.low[,"BIN_MEAN"], col=popcol[1],lwd=1,lty=1, t="l",
         main=paste("CHROMOSOME ",chr,"\n Pi",sep=""),
         xlab =paste("position ",sub.who,sep=""), ylab = "Normalised difference in Pi",
         ylim=range( (sub.pi.low[,(3+1*1:10)]- sub.pi.medium[,(3+1*1:10)]) / sub.pi.medium[,(3+1*1:10)], na.rm=T),
         cex.lab=5, cex.main=6,cex.axis=6)
    #add pop 1 its mean (mm.pi) and mean+2sd (ssd.pi) per chromosome: chromosome are in row and populations in columns
    abline(h=mm.pi[chr,1], col=popcol[1])
    abline(h=mm.pi[chr,1]+2*ssd.pi[chr,1], lty=2,col=popcol[1])
    #if there is an outlier window on the chromosome, add a grey polygon to show it
    if(length(outliers.fig[[chr]])>0){
      for(window in outliers.fig[[chr]]){
        polygon(y=c(-10,-10,20,20),
                x=c(sub.pi.low[window,"BIN_START"],sub.pi.low[window,"BIN_END"],sub.pi.low[window,"BIN_END"],sub.pi.low[window,"BIN_START"]),col=alpha("grey",0.6),border=FALSE)
      }
    }
    #thenadd the other populations: normalised difference, mean and mean+2sd
    for(pop in 2:10){
      lines(((sub.pi.low[,(3+1*pop)]- sub.pi.medium[,(3+1*pop)]) / sub.pi.medium[,(3+1*pop)])~sub.pi.low[,"BIN_MEAN"], col=popcol[pop],lwd=1,
            lty=1, t="l")
      abline(h=mm.pi[chr,pop],col=popcol[pop])
      abline(h=mm.pi[chr,pop]+2*ssd.pi[chr,pop],lty=2,col=popcol[pop])
    }
    legend("topleft", legend=c(pop.names), col=c(popcol,grey(0.4)), lty=c(rep(1,10),2),lwd=1.2,bty="n",cex=1.2)
}
  dev.off()
  
  
}
make.index1pdf(outliers.pionly)


################################################
#############   SupItem 2B       ###############
#############    Table S1        ###############
#    SFS   of all AOD window candidate regions #
################################################

#Useful for making SFS
color.set = c( "blue","orange")
legend.set = c("LR", "Any R")
pop.name = c("YRI","LWK","GBR","IBS","BEB","PJL","KHV","JPT","PEL","CLM")
numpop = length(pop.name)
numSFSentries = 21
pop.ind=list()
pop.ind[["YRI"]] = c(71:77,91:93) #AFR
pop.ind[["JPT"]] = c(78:84,88:90) #SAS
pop.ind[["GBR"]] <-1:10 #EUR
pop.ind[["CLM"]] <-11:20 #AMR
pop.ind[["IBS"]] <- 21:30 #EUR
pop.ind[["KHV"]] <-  c(31:33,44:48,51,52) #SAS
pop.ind[["PEL"]] <- 34:43 #AMR
pop.ind[["PJL"]] <- c(49,50,53:60) #EAS
pop.ind[["BEB"]] <- 61:70 #EAS
pop.ind[["LWK"]] <- c(85:87,94:100) #AFR



###Functions to compute sfs and make the plots, ds refers to downsampling
making.sfs = function(genot.table){#low is [,,1] vs all. [,,2]
  # SFS are counting the number of SNPs per population that are at frequency n/20 (10 indiv per pop, diploid) with n rqnging from 0 to 20.
  sfs = array(0,dim = c(numpop, numSFSentries, 2)) #2 stands for LR & all sites
  for (i in 1:numpop) {
    # I add 0:20 to make the table because otherwise it might miss one of the SFS values
    sfs[i,,1] = table(c(0:20,rowSums(genot.table[which(genot.table[,colrec]<=0.05),pop.ind[[pop.name[i]]]])))
    sfs[i,,2] = table(c(0:20,rowSums(genot.table[,pop.ind[[pop.name[i]]]])))
    sfs[i,,1] = sfs[i,,1] -1 #remove the 0:20 I added for the table() function
    sfs[i,,2] = sfs[i,,2] -1 
  }
  return(sfs)
}

making.sfs.ds = function(genot.table,sampling.size,nb.replicate=100){ #as making.sfs but where I downsample with replacement 
  low <- genot.table[which(genot.table[,colrec]<=0.05),] 
  sfs = array(0,dim = c(numpop, numSFSentries, 2,nb.replicate)) #2 stands for LR & all sites
  for(rep in 1:nb.replicate){
    # downsampling of sampling.size
    tmp.table <- genot.table[sample(1:nrow(genot.table), sampling.size,replace=T),]
    tmp.low <- low[sample(1:nrow(low), sampling.size,replace=T),]
    # same as making.sfs
    for (i in 1:numpop){
      sfs[i,,1,rep] = table(c(0:20,rowSums(  tmp.low[, pop.ind[[pop.name[i]]]])))
      sfs[i,,2,rep] = table(c(0:20,rowSums(tmp.table[, pop.ind[[pop.name[i]]]])))
      sfs[i,,1,rep] = sfs[i,,1,rep] -1 
      sfs[i,,2,rep] = sfs[i,,2,rep] -1 
    }
  }
  return(sfs)
} 

make.sfs.plots=function(aod,aod.around,main.set,my.title ){
  # make sfs plots for the AOD candidate region and for the surrounding region
  normal =  (1/1:19/sum(1/1:19))
  yylim=12.5
  sfs.aod = making.sfs(aod)
  sfs.aod.around = making.sfs(aod.around)
  nb.replicate=100
  rec=1
  for(pop in 1:length(pop.name)){
    sampling = nrow(aod[which(aod[,colrec]<=0.05 & rowSums(aod[,pop.ind[[pop.name[pop]]]])>0 & rowSums(aod[,pop.ind[[pop.name[pop]]]])<20),])
    if(rec==1){
     #Downsampling to make the sfs 
      sfs.aod.ds = making.sfs.ds(aod[which(aod[,colrec]<=0.05 & rowSums(aod[,pop.ind[[pop.name[pop]]]])>0 & rowSums(aod[,pop.ind[[pop.name[pop]]]])<20),],sampling,nb.replicate)
      sfs.aod.around.ds = making.sfs.ds(aod.around[which(aod.around[,colrec]<=0.05 & rowSums(aod.around[,pop.ind[[pop.name[pop]]]])>0 & rowSums(aod.around[,pop.ind[[pop.name[pop]]]])<20),],sampling,nb.replicate)
    }else{
      sfs.aod.ds = making.sfs.ds(aod[which( rowSums(aod[,pop.ind[[pop.name[pop]]]])>0 & rowSums(aod[,pop.ind[[pop.name[pop]]]])<20),],sampling,nb.replicate)
      sfs.aod.around.ds = making.sfs.ds(aod.around[which( rowSums(aod.around[,pop.ind[[pop.name[pop]]]])>0 & rowSums(aod.around[,pop.ind[[pop.name[pop]]]])<20),],sampling,nb.replicate)
      
    }
    #figures I will show the mean of the downsampling and the cinfidence intervals that are +2sd from it
    plot(NULL, xlim=c(1,19)/20, ylim = c(0,yylim), 
        main=paste( my.title, "\n", pop.name[pop]," (",sampling,")", sep=""),
           cex.main=6, cex.lab=6,cex.axis=6 , ylab = "Freq norm to neut expec", xlab = "Der allele freq")
    for(chr in 1:length(main.set)){
      if(chr==1){ cur.sfs = sfs.aod.ds}else{cur.sfs= sfs.aod.around.ds}
      m=std=c()
      ff= as.data.frame(matrix(0, ncol=nb.replicate, nrow = 19))
      for(rep in 1:nb.replicate) { ff[,rep] = cur.sfs[pop,2:20,rec,rep]/sum(cur.sfs[pop,2:20,rec,rep])/normal}
      m = rowMeans(ff); std = apply(ff,1,sd)
      polygon(x=c(1:19/20,19:1/20), y=c(m+2*std, rev(m-2*std)), col = alpha(color.set[chr],0.2),border=F)
      lines(y=m,x = 1:19/20, col=color.set[chr],lwd = 4)
    }
    if(pop==1) legend("top", main.set, bty="n", col = color.set,lwd = 4, cex = 6)
  }
}

pdf("SupItem2B.pdf", width=cm(30),height=cm(13))
layout.matrix= matrix(c(1,3,5,7,9,2,4,6,8,10),ncol=5,nrow=2,byrow=T) 11,13,15,17,19,12,14,16,18,20,21,23,25,27,29,22,24,26,28,30,31,33,35,37,39,32,34,36,38,40),ncol=5,nrow=8, byrow = T)
layout(layout.matrix)
par(mar=c(14,14,14,1), mgp=c(10,3,0),las=1)
#also I take as infos the chr, position start and end and average pi for all pop, mean DAFi LR and mean DAFi MR and number of SNP (nrow) LR and MR to then create tableS1
nrow.LR=nrow.MR=tab.chr=tab.start=tab.end =tab.dafi.LR = tab.dafi.MR= tab.pidiff = c()
jj=0
mm=0
for(chr in 1:22){
  #Are there outliers in that chromosome ? If so I will go on
  if(length(outliers.pionly[[chr]])>0){
    windows.chr = pi.low[(pi.low[,"CHROM"] == chr),1:3]
    positions.chr = windows.chr[sort(outliers.pionly[[chr]]),]
    genot.chr = genot.rec[genot.rec[,101]==chr,]
    pi.diff = rowMeans((pi.low[pi.low[,"CHROM"]==chr,][sort(outliers.pionly[[chr]]),4:13]-
           pi.medium[pi.medium[,"CHROM"]==chr,][sort(outliers.pionly[[chr]]),4:13])/ 
           pi.medium[pi.medium[,"CHROM"]==chr,][sort(outliers.pionly[[chr]]),4:13] )
    #get the posistion and concatenate if some windows are overlapping!!
    if(nrow(positions.chr)>1){
      k=1
      bound.position.chr=positions.chr[1,]
      mm=mm+1
      tab.pidiff[mm] = mean(pi.diff[1])
      for(pos in 2:nrow(positions.chr)){
        if(positions.chr[pos,2]<positions.chr[(pos-1),3]){
          bound.position.chr[nrow(bound.position.chr),3]<-positions.chr[pos,3]
          tab.pidiff[mm] = mean(rep(tab.pidiff[mm],k), pi.diff[pos] )
          k=k+1
        }else{
          k=1
          mm=mm+1
          bound.position.chr[nrow(bound.position.chr)+1,]<-positions.chr[pos,]
          tab.pidiff[mm] = pi.diff[pos]
        }
      }
    }else{
      mm=mm+1
      bound.position.chr=positions.chr
      tab.pidiff[mm] = pi.diff
    }
    print(bound.position.chr)
    #bound.position.chr gives the boundary of outlier window per chromosome

    for(pos in 1:nrow(bound.position.chr)){
      jj=jj+1
      start.pos = bound.position.chr[pos,2]
      end.pos = bound.position.chr[pos,3]
      #get subtable of the candidate region and it's surrounding area. For the surrounding I need to make sure not to overlap with an other candidate region..
      aod = genot.chr[(genot.chr[,102]>= start.pos & genot.chr[,102]<=end.pos),]
      aod.surrounding = genot.chr[((genot.chr[,102]>= max(0,start.pos-5000000,bound.position.chr[(pos-1),3] ) &
        genot.chr[,102]<start.pos) | 
        (genot.chr[,102]>end.pos & genot.chr[,102]<=min(end.pos+5000000,bound.position.chr[(pos+1),2],na.rm=T)) ),]
      #info for tab S1
      tab.chr[jj] = chr
      tab.start[jj]= start.pos
      tab.end[jj] = end.pos
      nrow.LR[jj] = sum(aod[,colrec]<=0.05)
      nrow.MR[jj] = sum(aod[,colrec]<0.5 & aod[,colrec]>=0.1)
      LR= aod[which(aod[,colrec]<=0.05),]
      tab.dafi.LR[jj] = mean(colMeans(LR[,1:100]))/2
      MR= aod[which(aod[,colrec]<0.5 & aod[,colrec]>=0.1),]
      tab.dafi.MR[jj] = mean(colMeans(MR[,1:100]))/2

      #plot for Item 2B
      make.sfs.plots(aod,aod.surrounding,
                     main.set = c("AOD Candidate", "Surrounding"),
                     my.title=paste("CHR",chr,":",start.pos,"-",end.pos,sep=""))
      }
    }
}

dev.off()


tabS1<-cbind(tab.chr,tab.start,tab.end,nrow.LR,nrow.MR,tab.pidiff, tab.dafi.LR,tab.dafi.MR)
write.table(tabS1, "tableS1.tab",col.names=T,row.names=F, sep="\t",quote=F)




################################################
#############   SupItem 2C       ###############
# HEATMAPS of all AOD window candidate regions #
################################################
#First, create temporary png files, 1 per heatmap and per population, then concatenate them

k=0
for(chr in 1:22){
  #Are there outliers in that chromosome ?
  if(length(outliers.pionly[[chr]])>0){
    windows.chr = pi.low[(pi.low[,"CHROM"] == chr),1:3]
    positions.chr = windows.chr[sort(outliers.pionly[[chr]]),]
    genot.chr = genot.rec[genot.rec[,101]==chr,]
    #get the posistion and concatenate if some windows are overlapping
    if(nrow(positions.chr)>1){
      bound.position.chr=positions.chr[1,]
      for(pos in 2:nrow(positions.chr)){
        if(positions.chr[pos,2]<positions.chr[(pos-1),3]){
          bound.position.chr[nrow(bound.position.chr),3]<-positions.chr[pos,3]
        }else{
          bound.position.chr[nrow(bound.position.chr)+1,]<-positions.chr[pos,]
        }
      }
    }else{
      bound.position.chr=positions.chr
    }
    #bound.position.chr gives the boundary of outlier window per chromosome
    for(pos in 1:nrow(bound.position.chr)){
      #for each of these region I do the heatmap of SNP that are LR and polymorphic among the population (between 1 and 19 derived alleles in total)
      k=k+1
      start.pos = bound.position.chr[pos,2]
      end.pos = bound.position.chr[pos,3]
      aod = genot.chr[(genot.chr[,102]>= start.pos & genot.chr[,102]<=end.pos),]
      for(pop in 1:10){
        pop.heatmap = aod[which(aod[,colrec]<=0.05 & rowSums(aod[,pop.ind[[pop.names[[pop]]]]])> 0 & rowSums(aod[,pop.ind[[pop.names[[pop]]]]]) < 20),pop.ind[[pop.names[[pop]]]]]
        colnames(pop.heatmap) <- paste("Ind ",1:10,sep="")
        mylmat = as.matrix(rbind(4:3,2:1)+4*(pop-1))
        png(paste("tmpHeatmap",pop,"_",k,".png",sep=""),width=600,height=600,res=300)
        heatmap.2(t(as.matrix(pop.heatmap)),
                  ylab="IND",xlab="LOCI",col=c("yellow", "orange", "red"),
                  breaks=c(-0.5,0.5,1.5,2.5),trace="none", main=paste("CHR",chr,":", start.pos,"-",end.pos,"\n", pop.names[[pop]],sep=""),dendrogram="row",cexRow = 10,cex.main=10,  labRow=F, labCol = F)
        dev.off()
      }
    }
  }
}
}

library(raster)
library(png)
#Concatenation of the figures per candidate region of 10 populations. It's ugly but I couldn't find another way to have several heatmaps on the same page
#I represent 4 candidate regions per page.
pdf(paste("SupItem2C.pdf",sep=""),height=60,width=48)
layout.matrix= matrix(c(1,3,5,7,9,2,4,6,8,10, 
                        11,13,15,17,19,12,14,16,18,20,
                        21,23,25,27,29,22,24,26,28,30,
                        31,33,35,37,39,32,34,36,38,40),ncol=5,nrow=8, byrow = T)
layout(layout.matrix)
for(k in 1:21){ 
  for(pop in 1:10){
    img <- readPNG( paste("tmpHeatmap",pop,"_",k,".png",sep=""))
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n",ylab="",xlab="")
    rasterImage(img,0,0,1,1)
  }
}
dev.off()

