setwd("~/Dropbox/GenomicAnalysisMusDom/20161130files/BEDfiles/")

####Next few sections are to make custom genome files so that intervals are only shuffled in regions of the genome with coverage

#Make exclude file to exclude ends of chromosomes with no coverage
genome.exclude=read.table("mm10_noRandom.bed",header=F,stringsAsFactors=F)
colnames(genome.exclude)=c("chr","length")
genome.exclude$window.start=0
genome.exclude$window.stop=0
#this was done for set of sliding windows - equivalent would be window.start=position of 1st SNP on each chromosome, window.stop - last SNP on each chr
#note: X is set to 20 in "chr" column
for(i in 1:20){
  genome.exclude$window.start[i]=min(window.stats[window.stats$chr==genome.exclude$chr[i],]$Window.start)
  genome.exclude$window.stop[i]=max(window.stats[window.stats$chr==genome.exclude$chr[i],]$Window.end)
}

#Make new genome file with chr lengths set to end of windows
mm10_max.windows.bed=genome.exclude[,c("chr","window.stop")]
write.table(mm10_max.windows.bed,file="mm10_max.windows.bed",row.names=F,col.names=F,sep="\t",quote=F,options(scipen=10))

#Make exclude file for beginnings 
mm10_exclude.beginnings.bed=genome.exclude[,1:3]
mm10_exclude.beginnings.bed[,2]=1
mm10_exclude.beginnings.bed[,3]=mm10_exclude.beginnings.bed[,3]-1
write.table(mm10_exclude.beginnings.bed,file="mm10_exclude.beginnings.bed",row.names=F,col.names=F,sep="\t",quote=F,options(scipen=10))
#for autosomes exclude whole X - make manual file 
mm10_exclude.beginnings.wholeX.bed=mm10_exclude.beginnings.bed
mm10_exclude.beginnings.wholeX.bed[20,3]=mm10_max.windows.bed$window.stop[20]
write.table(mm10_exclude.beginnings.wholeX.bed,file="mm10_exclude.beginnings.wholeX.bed",row.names=F,col.names=F,sep="\t",quote=F,options(scipen=10))

#functions for shuffling and counting overlap using bedtools
#if you want to compare peak SNPs overlap probably need a new function for counting that
bedTools.Shuffle.auto<-function(functionstring="bedtools shuffle",RegBed,GenomeBed,ExcludeBed)
{
  #create temp files
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  # create the command string and call the command using system()
  command=paste(functionstring,"-i",RegBed,"-g", GenomeBed,"-excl",ExcludeBed,"-noOverlapping",">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table(out,header=F)
  unlink(out)  
  return(res)
}

bedTools.Shuffle.X<-function(functionstring="bedtools shuffle",RegBed,GenomeBed,ExcludeBed)
{
  #create temp files
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  # create the command string and call the command using system()
  command=paste(functionstring,"-i",RegBed,"-g", GenomeBed,"-excl",ExcludeBed,"-chrom","-noOverlapping",">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table(out,header=F)
  unlink(out)  
  return(res)
}


bedTools.count.overlap<-function(functionstring="bedtools intersect",RegCounts,RegCompare)
{
  #create temp files
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",RegCounts,"-b", RegCompare,"-wao",">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table(out,header=F)
  unlink(out)  
  return(res)
}

#for overlap use column 6 instead of 7 because overlap bp for single bp=0
get.overlap.count=function(shuffled.perm.row,count.perm) sum(count.perm[count.perm[,8]==shuffled.perm.row[4],6]>0)
get.overlap.bp=function(shuffled.perm.row,count.perm) sum(count.perm[count.perm[,8]==shuffled.perm.row[4],7])

#read in bed file with all significant intervals and split into autosomal and X - will permute separately
bed.file=read.table(paste(outlier.p99.stats$Stat[i],".outlier.p99.regions.20170419.bed",sep=""),stringsAsFactors=F)
write.table(bed.file[bed.file[,1]!="chrX",],file="outliers.auto.temp.bed",row.names=F,col.names=F,sep="\t",quote=F,options(scipen=10))
write.table(bed.file[bed.file[,1]=="chrX",],file="outliers.X.temp.bed",row.names=F,col.names=F,sep="\t",quote=F,options(scipen=10))


#need bed file for the regions you are comparing shuffled to "RegCompare" - here is "GWAScollapsed.bed"
#Make vectors for output of 10,000 permutations for # overlapped regions and bp overlap
overlap.counts=rep(-1,10000)
overlap.bp=rep(-1,10000)
reg.count.check=rep(-1,10000) 
for (j in 7478:10000){
  auto.shuffled.perm=bedTools.Shuffle.auto(RegBed="outliers.auto.temp.bed",GenomeBed="mm10_max.windows.20170412.bed",ExcludeBed="mm10_exclude.beginnings.wholeX.20170412.bed")
  x.shuffled.perm=bedTools.Shuffle.X(RegBed="outliers.X.temp.bed",GenomeBed="mm10_max.windows.20170412.bed",ExcludeBed="mm10_exclude.beginnings.20170412.bed")
  write.table(auto.shuffled.perm,file="Outliers.shuffled.temp.bed",row.names=F,col.names=F,sep="\t",quote=F,options(scipen=10))
  write.table(x.shuffled.perm,file="Outliers.shuffled.temp.bed",row.names=F,col.names=F,sep="\t",quote=F,options(scipen=10))
  count.perm=bedTools.count.overlap(RegCounts="Outliers.shuffled.temp.bed",RegCompare="GWAScollapsed.bed")
  overlap.counts[j]=sum(count.perm[,7]>0)
  overlap.bp[j]=sum(count.perm[,7])
  reg.count.check[j]=nrow(count.perm)
  unlink("Outliers.shuffled.temp.bed")
}
unlink("outliers.auto.temp.bed")
unlink("outliers.X.temp.bed")
table(reg.count.check==outlier.p99.stats$NumRegions[i]) #should all be TRUE
#copy the values into output table for all traits/stats
outlier.p99.stats$perm10000.NumOverlapGWAS[i]=sum(overlap.counts>=outlier.p99.stats$NumOverlapGWAS[i])
outlier.p99.stats$perm10000.bpOverlapGWAS[i]=sum(overlap.bp>=outlier.p99.stats$bpOverlapGWAS[i])

