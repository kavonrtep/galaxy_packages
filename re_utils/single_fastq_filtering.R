#!/usr/bin/env Rscript



library(optparse,quiet=TRUE)
matricesSum = function(m1,m2){
    mfin = matrix(0, nrow = nrow(m2), ncol = max(ncol(m1), ncol(m2)))
    rownames(mfin) = rownames(m2)
    if (sum(m1)>0){
        mfin[,1:ncol(m1)] = mfin[,1:ncol(m1)] + m1
    }
    if (sum(m2)>0){
        mfin[,1:ncol(m2)] = mfin[,1:ncol(m2)] + m2
    }

    return(mfin)
}

plotFrequencies = function(x){
    par(xaxs = 'i', yaxs = 'i')
    plot(NULL, xlim = c(1,ncol(x)), ylim = c(0,1),type ='n', xlab="position", ylab = "proportion")
    col = c(A="red", "T"="orange", C="blue", G="green")
    for (n in c("A","C", "T","G")){
        points(x[n,]/colSums(x), type = 'l', col = col[n])
    }
    grid(ny = 50, nx = ncol(x) - 1  )
    legend("topright", col = col, lwd = 2, legend = names(col))
}

megablast = function(fx,params="-F \"m D\" -D 4 -p 90 ",database){
  # assume blastn and makeblastdb in $PATH
	tmp_in = tempfile()
	tmp_out = tempfile()
	writeFasta(fx,file=tmp_in)
  ## check if db is present:
  dbfiles = paste0(database, ".nhr")
  if (!file.exists(dbfiles)){
    cat("running makeblastdb\n")
    system(paste0("makeblastdb -dbtype nucl -in ", database))
  }
  params = '-num_alignments 1 -num_threads 4  -outfmt "6 qseqid qcovs pident" '
	cmd=paste("blastn", " -db", database, "-query", tmp_in,
            "-out", tmp_out, " ", params)
  status=system(cmd,intern=TRUE)

	if (file.info(tmp_out)$size != 0){
    results=read.table(tmp_out,sep="\t",as.is=TRUE,comment.char="")
    colnames(results) = c("qseqid", "qcovs", "pident")
	}else{
		results=NULL
	}
  unlink(tmp_in)
	unlink(tmp_out)
	return(results)
}



option_list = list(
  make_option(c('-a', '--fastqA'),action='store',type='character',help='fastq file A',default=NA),
  make_option(c('-x', '--fastqX'),action='store',type='character',help='output fastq file X',default=NA),
  make_option(c('-c', '--cut_off'),action='store',type='numeric',help='Quality cut-off value [default %default]',default=10),
  make_option(c('-r', '--rnd_seed'),action='store',type='numeric',help='seed for random number generator [default %default]',default=123),
  make_option(c('-p', '--percent_above'),action='store',type='numeric',
              help='Percent of bases in sequence that must have quality equal to / higher than cut-off value [default %default]',default=95),
  make_option(c("-G", "--png_file_output"), action = "store", type = "character", default=NA),
  make_option(c('-e', '--trim_end'), action='store',type='numeric',help="trimming - end position [no trimming by default]", default=NA),
  make_option(c('-s', '--trim_start'), action='store',type='numeric',help="triming position - start  [no trimming by default]", default=NA),

  make_option(c('-n', '--sample_size'),action='store',type='numeric',help='requested sample size (number of pairs)[no sampling by default]',default=NULL),
  make_option(c('-N', '--Max_N'),action='store',type='integer',help='maximum number of Ns in sequences [default %default]',default=0),
  make_option(c('-f', '--format'),action='store',type='character',help='format of output - fastq or fasta [default %default] ',default="fasta"),
  make_option(c('-C', '--cutadapt_options'),action='store',type='character',help='file specifying cutadapt options',default=NULL),
  make_option(c('-F', '--filter_seq'),action='store',type='character',help='file specifying sequences for filtering (e.g. plastid DNA)',default=NULL),
  make_option(c('-j', '--chunk_size'),action='store',type='numeric',help='Number of sequences processed in single step. This option affect speed of processing and memory usage [default %default]',default=1000000)





)

description=paste(strwrap('Script for filterinq fastq file f
						'),collapse="\n")
# arguments whic hcould be modifed
cutadapt_arguments=paste(
		" --anywhere='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' ",
		" --anywhere='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'  ",
		" --anywhere='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' --anywhere='ATCTCGTATGCCGTCTTCTGCTTG' ",
		" --anywhere='CAAGCAGAAGACGGCATACGAGAT' --anywhere='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC' ",
		"--error-rate=0.05 --times=1 --overlap=15 --discard ",
		" ",sep='')


epilogue=paste(strwrap(c('Description:\n
This tool is designed to make memory efficient preprocessing of single fastq files. Input files can be in GNU zipped archive (.gz extension).
Reads are filtered based on the quality, presence of N bases and adapters.
Cutaddapt program must be installed to use adapter filtering.
The input files are process in chunks,
All reads which pass the quality filter fill be writen into output file.
If --sample size is specified, only sample of sequences will be returned.
By default cutadapt us run with this options:
\n ',
cutadapt_arguments,
"\n\nIf you want to use different adapter sequences specify them in separate file and use -C --cutadapt options.
Consult cutadapt manual for usage. ")
),collapse="\n")

epilogue=paste(epilogue,"


# example:
               TODO

",sep='\n')


parser=OptionParser(
		option_list=option_list,
		epilogue=epilogue,
		description=description,
)



opt = parse_args(parser, args=commandArgs(TRUE))
if  (!(opt$format %in% c("fasta","fastq"))){
	stop("wrong output format")
}

if (any(is.na(c(opt$fastqA, opt$fastqX )))){
	cat("\nInput ond output files must be specified!\n\n")
	print_help(parser)
	q()
}


if (!is.null(opt$cutadapt_options)){
	cutadapt_argumenns=scan(opt$cutadapt_options,what="character", comment.char="#")
	#remove "EOF"
	cutadapt_arguments=paste(" ",cutadapt_arguments," ",collapse=' ',sep=" ")
}


cutadapt_cmd=function(input,output,cutadapt_arguments){
	cmds=paste("cutadapt --format=",
			"fastq",
			cutadapt_arguments,
			" ",
			c(input),
			" -o ",
			c(output),
			sep="")

}

set.seed(opt$rnd_seed)

PROP=opt$percent_above/100
MINQ=opt$cut_off
TRIM_END=opt$trim_end
TRIM_START=opt$trim_start

CHUNK_SIZE=opt$chunk_size

MIN_LENGTH=TRIM_END - ifelse(is.na(TRIM_START),1,TRIM_START)+1
if (is.na(MIN_LENGTH)){
	MIN_LENGTH=1
}

suppressPackageStartupMessages(library(ShortRead))
Qfilter=srFilter(function(x){
			rowSums(as(quality(x),'matrix')>MINQ,na.rm=TRUE)/nchar(sread(x))>=PROP
		}
)

Min_Length_filter=srFilter(function(x){
			nchar(sread(x))>=MIN_LENGTH
		}
)


fin_filter=compose(Qfilter,nFilter(opt$Max_N),Min_Length_filter)

if (opt$format=="fasta"){
	writeFun=writeFasta
}else{
	writeFun=writeFastq
}


f1out=opt$fastqX

#find wich character specify pair: - it could be

# check size of sequences - must be same...
n1=countLines(opt$fastqA)/4

cat("number of sequences in ",opt$fastqA,":",n1,"\n")
cat("--------\n")


number_of_chunks=round(n1/CHUNK_SIZE)
CHUNK_SIZE = round(n1/number_of_chunks)

if (number_of_chunks==0){
	CHUNK_SIZE=n1
	number_of_chunks=1
}
if (!is.null(opt$sample_size)){
	sample_size_in_chunk=round(opt$sample_size/number_of_chunks)
  n_missing = opt$sample_size - sample_size_in_chunk * number_of_chunks
}else{
	sample_size_in_chunk=CHUNK_SIZE
  n_missing = 0
}

# adjust the chunk size to get exact count of sequences:
CHUNK_SIZE = ceiling(n1/number_of_chunks)
save.image("tmp.RData")


f1 <- FastqStreamer(opt$fastqA, CHUNK_SIZE)
total=0
nucleotideFrequenciesForward = matrix(0)
while(TRUE){
	print (Sys.time())
	fq1 <- yield(f1)
	print (Sys.time())
	if (length(fq1)==0){
		break
	}
	
	
	
	
	if (!is.na(TRIM_END)){
		fq1=narrow(fq1,end=ifelse((TRIM_END)>width(fq1),width(fq1),TRIM_END))
	}
	if (!is.na(TRIM_START)){
		fq1=narrow(fq1,start=ifelse((TRIM_START)<width(fq1),TRIM_START,width(fq1)))
	}
	


	fqF1=fq1[fin_filter(fq1)]

	cat("running cutadapt on chunk\n")
	# cut addapt here: # filter parts:
	tmp_in1=tempfile(fileext=".fastq")
	tmp_out1=tempfile(fileext=".fastq")
        if (is.null(formals(writeFastq)$compress)){
            writeFastq(fqF1,file=tmp_in1)
        }else{
            writeFastq(fqF1,file=tmp_in1, compress = FALSE)
        }


	cmd1=cutadapt_cmd(tmp_in1,tmp_out1, cutadapt_arguments)

	status=system(cmd1,intern=TRUE)
	print (Sys.time())
	# collect results
	ftmp1=FastqFile(tmp_out1)

	fqF1=readFastq(ftmp1)
	close(ftmp1)

  # remove sequences similar to filter database (e.g. plastid DNA)
  if (!is.null(opt$filter_seq)){
		blast_results =  megablast(fqF1, database=opt$filter_seq)
		if (!is.null(blast_results)){
			exclude=with(blast_results, unique(qseqid[pident >= 90 & qcovs >=90]))
      ## note - blast will truncate seq ids - anything beyond space is omitted
			seq_to_exclude = gsub(" .*","",id(fqF1)) %in% exclude
			excluded_contamination=signif(sum(seq_to_exclude)/length(seq_to_exclude)*100,3)
      cat(excluded_contamination,"% filtered out after blast\n" )
			fqF1 = fqF1[!seq_to_exclude]
		}
	}



	# clean up
	unlink(tmp_out1)
	unlink(tmp_in1)


  

	# filter complete pairs again:

	# create new id - last character must differentiate pair - for interlacig
	if (length(fqF1)>(sample_size_in_chunk + n_missing)){
		smp=sort(sample(seq_along(fqF1),sample_size_in_chunk + n_missing))
    n_missing = 0
		writeFun(fqF1[smp],file=f1out,mode='a')
    nfrq1 = alphabetByCycle(sread(fqF1[smp]))

	}else{
      writeFun(fqF1,file=f1out,mode='a')
      nfrq1 = alphabetByCycle(sread(fqF1))
	}
  nucleotideFrequenciesForward = matricesSum(nucleotideFrequenciesForward, nfrq1)

}

if( !is.na(opt$png_file_output)){
    png(opt$png_file_output, width = 1000, height = 500)
    plotFrequencies(nucleotideFrequenciesForward)
    dev.off()
}

close(f1)




