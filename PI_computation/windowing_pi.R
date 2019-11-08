#   *************************************************************************
#   ***           Scripts to perform the human analyses                   ***
#   ***   from the paper Transition from background selection to          ***
#   ***  associative overdominance promotes diversity in regions of       ***
#   ***                low recombination                                  ***
#   ***   by Gilbert, Pouyet, Excoffier and Peischl                       ***
#   ***   Theses scripts were prepared by Fanny Pouyet on 2019.           ***
#   ***              Last update 22.08.2019                               ***
#   *************************************************************************

#Merging outputs vcftools
rm(list=ls())
begining.colnames = c("CHROM","BIN_START","BIN_END")
pop.names =  c("YRI","LWK","GBR","IBS","BEB","PJL","KHV","JPT","CLM","PEL")
window.range = read.table("windows2Mb_step0p5Mb",h=T)

for(rec in c("lowRec0p05","RR0p1-0p5")){
  print(rec)
  bed.file.low.rec = read.table(paste("strict_mask_YRI.recomb.",rec,".noCpG.bed",sep=""))
  covered.low.rec=list()
  for(chrom in 1:22){
    sub.window.range = window.range[(window.range[,1]==chrom),]
    covered.low.rec[[chrom]]=rep(0,sum(window.range[,1]==chrom))
    # create output table
    for(i in 1:sum(window.range[,1]==chrom)){
         low.rec = bed.file.low.rec[which(bed.file.low.rec[,1] == sub.window.range[i,1] & bed.file.low.rec[,2]<=sub.window.range[i,3] & bed.file.low.rec[,3] >=sub.window.range[i,2]),]
         covered.low.rec[[chrom]][i] = sum(low.rec[,3] -low.rec[,2])
    }
  }

  print("start_compute")
  for(chrom in c(1,10:19,2,20:22, 3:9)){
    pi.sites = read.table(paste("chr",chrom,".",rec,".sites.pi",sep=""),h=T)
    #initialisation
    window.pi <- cbind(window.range[window.range[,1]==chrom ,], as.data.frame(matrix(0,ncol=10,nrow= sum(window.range[,1]==chrom))) )    
    # create output table
    for(i in 1:nrow(window.pi)){
      if(length(which(pi.sites[,2]<=window.pi[i,3] & pi.sites[,2]>=window.pi[i,2])) > 0){
       window.pi[i,4:ncol(window.pi)] <- colSums(pi.sites[which(pi.sites[,2]<=window.pi[i,3] & pi.sites[,2]>=window.pi[i,2]),4:13]) / covered.low.rec[[chrom]][i] 
      }
    }
    if(!exists("pi.output")){
      pi.output<-window.pi
    }else{
      pi.output<- rbind(pi.output,window.pi)
    }
  }
  print("wrintingoutput")
  colnames(pi.output)<- colnames(pi.sites)
  write.table(pi.output, paste("nucdiv.",rec,".windows2Mb-step0p5Mb.pi", sep=""),col.names = T,quote=F,row.names = F,sep="\t")
  rm(pi.output)
}
