#! /usr/bin/env Rscript
# repeat masker summary
# analysis of *.out file
# input arguments: <rmsk_result.out> <reads.fas.
# calculates totoal proportion of matched repetetive families and averagee score and print the table 

opt = list()
opt$fasta = commandArgs(T)[1]
opt$repeat_masker = commandArgs(T)[2]
opt$output = commandArgs(T)[3]


RMfiles=system(paste("ls ",opt$repeat_masker,sep=""),intern=T)
fasta=system(paste("ls ",opt$fasta,sep=""),intern=T)



getsubdir=function(x){
	xx=gsub(".*/","",x)
  sdir=gsub(xx,"",x,fixed=T)
	sdir
}

fasta.summary=function(reads){
	suppressPackageStartupMessages(library(Biostrings,quiet=T))
	if (exists("readDNAStringSet")){
		seqs=readDNAStringSet(reads)	
	}else{
		seqs=read.DNAStringSet(reads)
	}
	
	N=length(seqs)
	Ls=nchar(seqs)
	L=sum(Ls)
	M=median(Ls)
	A=mean(Ls)
	output=c(N,L,M,A)	
	names(output)=c("N","length","median","mean")
	output
}

repCont=function(repname){   # uses external variable!
	pos=repname==classfamily
	RC=unname(sum(lengths[pos])/totalLength)
	RC
}
repScore=function(repname){   # uses external variable!
	pos=repname==classfamily
	S=unname(mean(score[pos]))
	S
}

N=length(RMfiles)

for (i in seq_along(RMfiles)){
	cat(paste((RMfiles[i]),"\n"))
	
	sdir=getsubdir(RMfiles[i])
	seqs.summary=fasta.summary(fasta[i])
	totalLength=seqs.summary["length"]
	
	total.counts=data.frame(c("All_Reads_Length","All_Reads_Number"),c(NA,NA),c(NA,NA),c(totalLength,seqs.summary['N']),c(totalLength,seqs.summary['N']),c(NA,NA))
	colnames(total.counts)=c("class/fammily","class","family","hits","content","Mean_Score")
	

	# get RM data if not empty
	rmsk_empty=length(grep("There were no repetitive sequences detected",scan(RMfiles[i],what=character(),n=1,sep="\n",quiet=T),fixed=T))>0
	rmsk_notperformed = length(grep("No Repbase used with RepetMasker",scan(RMfiles[i],what=character(),n=1,sep="\n",quiet=T),fixed=T))>0
	if (!rmsk_empty & !rmsk_notperformed) {
		rmsk=read.table(RMfiles[i],header=F,as.is=T,skip=2,comment.char="*",fill=T)
		score=rmsk[,1]
		Ids=rmsk[,5]
		starts=rmsk[,6]
		ends=rmsk[,7]
		lengths=ends-starts
		classfamily=rmsk[,11]
		out.table=table(classfamily)
		out.table=data.frame(names(out.table),c(out.table),stringsAsFactors=F)
		contents=sapply(out.table[,1],repCont)*100
		meanScore=sapply(out.table[,1],repScore)
		class=gsub("/.*","",out.table[,1])
		families=gsub(".*/","",out.table[,1])
		out.table=data.frame(out.table[,1],class,families,out.table[,-1],contents,meanScore,stringsAsFactors=F)
		colnames(out.table)=c("class/fammily","class","family","hits","content","Mean_Score")
		out.table=out.table[order(out.table$content,decreasing=T),]
		out.table=rbind(out.table,total.counts)
	}else{
		out.table=total.counts
	}
	out.table=out.table[order(out.table$content,decreasing=T),]
	
	write.table(out.table,file=paste(sdir,opt$output,sep=""),row.names=F,quote=F,sep="\t")

	
	# merge all tables:
	
	colnames(out.table)=paste(colnames(out.table),c("",rep(sdir,5)))

	if (i==1) {
		mergedTable=out.table[,c(1,4,5,6)]
	}else{
		
		mergedTable=merge(mergedTable,out.table[,c(1,4,5,6)],by=1,all=T)
	}
	
}
	mergedTable=mergedTable[order(rowSums(mergedTable[,seq(2,N*3,by=3),drop=FALSE],na.rm=T),decreasing=T),]
	mergedTable[is.na(mergedTable)]<-0
	
	write.table(t(mergedTable),file=paste(opt$output,"summary.csv",sep=""),col.names=F,row.names=T,quote=F,sep="\t")

#save.image("tmp.RData")








