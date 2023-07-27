#!/usr/bin/env Rscript
# TODO - add triming from start and end, removing short reads
library(optparse, quiet = TRUE)
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

option_list = list(make_option(c("-a", "--fastqA"), action = "store", type = "character", 
    help = "fastq file A", default = NA), make_option(c("-b", "--fastqB"), action = "store", 
    type = "character", help = "fastq file B", default = NA), make_option(c("-x", 
    "--fastqX"), action = "store", type = "character", help = "output fastq file X", 
    default = NA), make_option(c("-y", "--fastqY"), action = "store", type = "character", 
    help = "output fastq file Y", default = NA), make_option(c("-c", "--cut_off"), 
    action = "store", type = "numeric", help = "Quality cut-off value [default %default]", 
    default = 10), make_option(c("-r", "--rnd_seed"), action = "store", type = "numeric", 
                               help = "seed for random number generator [default %default]", default = 123),
    make_option(c("-G", "--png_file_output"), action = "store", type = "character", default=NA),
    make_option(c("-p", "--percent_above"), action = "store", type = "numeric", help = "Percent of bases in sequence that must have quality equal to / higher than cut-off value [default %default]", 
        default = 95), make_option(c("-e", "--trim_end"), action = "store", type = "numeric", 
        help = "trimming - end position [no trimming by default]", default = NA), 
    make_option(c("-s", "--trim_start"), action = "store", type = "numeric", help = "triming position - start  [no trimming by default]", 
        default = NA), make_option(c("-n", "--sample_size"), action = "store", type = "numeric", 
        help = "requested sample size (number of pairs)[no sampling by default]", 
        default = NULL), make_option(c("-N", "--Max_N"), action = "store", type = "integer", 
        help = "maximum number of Ns in sequences [default %default]", default = 0), 
    make_option(c("-R", "--rename"), action = "store_true", type = "logical", help = "Rename sequences", 
        default = FALSE), make_option(c("-f", "--format"), action = "store", type = "character", 
        help = "format of output - fastq or fasta [default %default] ", default = "fasta"), 
    make_option(c("-C", "--cutadapt_options"), action = "store", type = "character", 
        help = "file specifying cutadapt options", default = NULL), make_option(c("-j", 
        "--chunk_size"), action = "store", type = "numeric", help = "Number of sequences processed in single step. This option affect speed of processing and memory usage [default %default]", 
        default = 1000000),
    make_option(c('-F', '--filter_seq'),action='store',type='character',help='file specifying sequences for filtering (e.g. plastid DNA)',default=NULL)
    )

description = paste(strwrap("Script for filterinq fastq files generated from paired end sequencing
\t\t\t\t\t\t"), 
    collapse = "\n")
# arguments whic hcould be modifed
cutadapt_arguments = paste(" --anywhere='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' ", 
    " --anywhere='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'  ", 
    " --anywhere='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' --anywhere='ATCTCGTATGCCGTCTTCTGCTTG' ", 
    " --anywhere='CAAGCAGAAGACGGCATACGAGAT' --anywhere='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC' ", 
    "--error-rate=0.05 --times=1 --overlap=15 --discard ", " ", sep = "")


epilogue = paste(strwrap(c("Description:\n
This tool is designed to make memory efficient preprocessing of two fastq files. Input files can be in GNU zipped archive (.gz extension).
Reads are filtered based on the quality, presence of N bases and adapters.
Two input fastq files are procesed in parallel. Cutaddapt program must be installed to use adapter filtering.
Only complete pair are kept. As the input files are process in chunks,
it is required that pair reads are complete and in the same order in both input files. All reads which pass the quality
filter fill be writen into ouytput files. If --sample size is specified, only sample of sequences will be returned.
By default cutadapt us run with this options:
\n ", 
    cutadapt_arguments, "\n\nIf you want to use different adapter sequences specify them in separate file and use -C --cutadapt options.
Consult cutadapt manual for usage. ")), 
    collapse = "\n")

epilogue = paste(epilogue, "
fastqA                    fastqB
   |                         |
   V                         V
----------------------------------
|      Trimming (optional)       |
----------------------------------
   |                         |
   |                         |
   V                         V
----------------------------------
|        Filter by quality       |
----------------------------------
   |                         |
   |                         |
   V                         V
---------------------------------
|        Discard single reads   |--->
|        Keep complete pairs    |
---------------------------------
   |                         |
   V                         V
----------------------------------
|        cutadapt filtering    |
----------------------------------
   |                         |
   V                         V
---------------------------------
|        Discard single reads   |--->
|        Keep complete pairs    |
---------------------------------
   |                         |
   |                         |
   V                         V
----------------------------------
|      sample (optional)         |
----------------------------------
   |                         |
   |                         |
   V                         V
 FastqX                    FastqY
(or fastaX)             (or fastaY)

# example:

", 
    sep = "\n")


parser = OptionParser(option_list = option_list, epilogue = epilogue, description = description, 
    )



opt = parse_args(parser, args = commandArgs(TRUE))
if (!(opt$format %in% c("fasta", "fastq"))) {
    stop("wrong output format")
}
if (any(is.na(c(opt$fastqA, opt$fastqB, opt$fastqX, opt$fastqY)))) {
    cat("\nInput ond output files must be specified!\n\n")
    print_help(parser)
    q()
}


if (!is.null(opt$cutadapt_options)) {
    if (file.exists(opt$cutadapt_options)) {
        cutadapt_arguments = scan(opt$cutadapt_options, what = "character", comment.char = "#")
        
    } else {
        cutadapt_arguments = opt$cutadapt_options
    }
    ## remove 'EOF'
    cutadapt_arguments = paste(" ", cutadapt_arguments, " ", collapse = " ", sep = " ")
}



cutadapt_cmd = function(input, output, cutadapt_arguments) {
    cmds = paste("cutadapt --format=", "fastq", cutadapt_arguments, " ", c(input), 
        " -o ", c(output), sep = "")
}

set.seed(opt$rnd_seed)
PROP = opt$percent_above/100
MINQ = opt$cut_off
TRIM_END = opt$trim_end
TRIM_START = opt$trim_start

CHUNK_SIZE = opt$chunk_size

MIN_LENGTH = 1 + TRIM_END - ifelse(is.na(TRIM_START), 1, TRIM_START)
if (is.na(MIN_LENGTH)) {
    MIN_LENGTH = 1
}

suppressPackageStartupMessages(library(ShortRead))
Qfilter = srFilter(function(x) {
    rowSums(as(quality(x), "matrix") > MINQ, na.rm = TRUE)/nchar(sread(x)) >= PROP
})

Min_Length_filter = srFilter(function(x) {
    nchar(sread(x)) >= MIN_LENGTH
})




fin_filter = compose(Qfilter, nFilter(opt$Max_N), Min_Length_filter)

if (opt$format == "fasta") {
    writeFun = writeFasta
} else {
    writeFun = writeFastq
}

f1out = opt$fastqX
f2out = opt$fastqY
# are output files empty?:


# find wich character specify pair: - it could be
lA = readLines(opt$fastqA, n = 1)
lB = readLines(opt$fastqB, n = 1)
pair_tag_position = which(!strsplit(lA, "")[[1]] == strsplit(lB, "")[[1]])

if ((length(pair_tag_position)) != 1 & !opt$rename) {
    warning("unable to determine pair tag")
}

# pair_regex=paste(substring(lA,pair_tag_position-1,pair_tag_position-1),'.+',sep='')
# assume 1/ same order 2/ all id must be same

# check size of sequences - must be same...
n1 = countLines(opt$fastqA)/4


number_of_chunks = round(n1/CHUNK_SIZE)
## adjust chunk size to make last chunk of the same size is all other
## this is to avoid small last chunk
CHUNK_SIZE = round(n1/number_of_chunks)

if (number_of_chunks == 0) {
    CHUNK_SIZE = n1
    number_of_chunks = 1
}
if (!is.null(opt$sample_size)) {
  sample_size_in_chunk = round(opt$sample_size/number_of_chunks)
  n_missing = opt$sample_size - sample_size_in_chunk * number_of_chunks
} else {
  sample_size_in_chunk = CHUNK_SIZE
  n_missing = 0
}

cat("number chunks ", number_of_chunks, "\n")
cat("chunks size ", CHUNK_SIZE, "\n")
# adjust the chunk size to get exact count of sequences:
CHUNK_SIZE = ceiling(n1/number_of_chunks)
F_id = ifelse(opt$rename, "/1", "1")
R_id = ifelse(opt$rename, "/2", "2")

f1 <- FastqStreamer(opt$fastqA, CHUNK_SIZE)
f2 <- FastqStreamer(opt$fastqB, CHUNK_SIZE)
total = 0
chunk = 0
nucleotideFrequenciesForward = nucleotideFrequenciesReverse = matrix(0)
while (TRUE) {
    chunk = chunk + 1
    cat("chunk number ", chunk, "\n")
    fq1 <- yield(f1)
    fq2 <- yield(f2)
    if (length(fq1) == 0) {
        break
    }
    cat("chunk number ", chunk, " imported\n")
    cat("chunk size", length(fq1), "\n")
    ## rename
    chunk_id = sprintf(paste0("%0", round(log10(number_of_chunks)) + 1, "d"), chunk)
    cat("chunk id ", chunk_id, "\n")
    fmt = paste0("%0", round(log10(length(fq1))) + 1, "d")
    
    ## either do not rename
    if (opt$rename) {
        ## or use index
        fq1@id = fq2@id = BStringSet(paste0(chunk_id, sprintf(fmt, 1:length(fq1))))
    } else {
        ## or use largest common substring
        ind = mapply(function(x, y) which(x != y)[1], strsplit(as.character(id(fq1)), 
            split = ""), strsplit(as.character(id(fq2)), split = ""))
        fq1@id = fq2@id = subseq(id(fq1), 1, ind - 1)
        
    }
    
    if (!is.na(TRIM_END)) {
        fq1 = narrow(fq1, end = ifelse((TRIM_END) > width(fq1), width(fq1), TRIM_END))
        fq2 = narrow(fq2, end = ifelse((TRIM_END) > width(fq2), width(fq2), TRIM_END))
    }
    if (!is.na(TRIM_START)) {
        fq1 = narrow(fq1, start = ifelse((TRIM_START) < width(fq1), TRIM_START, width(fq1)))
        fq2 = narrow(fq2, start = ifelse((TRIM_START) < width(fq2), TRIM_START, width(fq2)))
        
    }
    
    
    fqF1 = fq1[fin_filter(fq1)]
    fqF2 = fq2[fin_filter(fq2)]
    
    inc1 = id(fqF1) %in% id(fqF2)
    inc2 = id(fqF2) %in% id(fqF1)
    
    cat("running cutadapt on chunk\n")
    # cut addapt here: # filter parts:
    tmp_in1 = tempfile(fileext = ".fastq")
    tmp_out1 = tempfile(fileext = ".fastq")
    tmp_in2 = tempfile(fileext = ".fastq")
    tmp_out2 = tempfile(fileext = ".fastq")
    
    
    if (is.null(formals(writeFastq)$compress)) {
        writeFastq(fqF1[inc1], file = tmp_in1)
        writeFastq(fqF2[inc2], file = tmp_in2)
    } else {
        writeFastq(fqF1[inc1], file = tmp_in1, compress = FALSE)
        writeFastq(fqF2[inc2], file = tmp_in2, compress = FALSE)
    }
    
    
    cmd1 = cutadapt_cmd(tmp_in1, tmp_out1, cutadapt_arguments)
    cmd2 = cutadapt_cmd(tmp_in2, tmp_out2, cutadapt_arguments)
    # this should run cutadapt in parallel
    
    status = system(paste(cmd1, "& \n", cmd2, "&\nwait"), intern = TRUE)
    # collect results
    ftmp1 = FastqFile(tmp_out1)
    
    fqF1 = readFastq(ftmp1)
    close(ftmp1)
    
    ftmp2 = FastqFile(tmp_out2)
    fqF2 = readFastq(ftmp2)
    close(ftmp2)
    ## clean up
    unlink(tmp_out1)
    unlink(tmp_in1)
    unlink(tmp_out2)
    unlink(tmp_in2)
    


    ## remove sequences similar to filter database (e.g. plastid DNA)
    if (!is.null(opt$filter_seq)){
      blast_results1 =  megablast(fqF1, database=opt$filter_seq)
      blast_results2 =  megablast(fqF2, database=opt$filter_seq)
      if (!is.null(blast_results1) &  !is.null(blast_results2)){
        exclude1=with(blast_results1, unique(qseqid[pident >= 90 & qcovs >=90]))
        exclude2=with(blast_results2, unique(qseqid[pident >= 90 & qcovs >=90]))
        ## note - blast will truncate seq ids - anything beyond space is omitted
        seq_to_exclude1 = gsub(" .*","",id(fqF1)) %in% exclude1
        seq_to_exclude2 = gsub(" .*","",id(fqF2)) %in% exclude2
        excluded_contamination1=signif(sum(seq_to_exclude1)/length(seq_to_exclude1)*100,3)
        excluded_contamination2=signif(sum(seq_to_exclude2)/length(seq_to_exclude2)*100,3)
        cat(excluded_contamination1,"% filtered out after blast - forward reads\n" )
        cat(excluded_contamination2,"% filtered out after blast - reverse reads\n" )
        fqF1 = fqF1[!seq_to_exclude1]
        fqF2 = fqF2[!seq_to_exclude2]
      }
    }



    
    # filter complete pairs again: id1=gsub(pair_regex,'',id(fqF1))
    # id2=gsub(pair_regex,'',id(fqF2))
    inc1 = id(fqF1) %in% id(fqF2)
    inc2 = id(fqF2) %in% id(fqF1)
    total = sum(inc1) + total
    ## create new id - last character must differentiate pair - for interlacig
    
    fqF1@id = BStringSet(paste0(id(fqF1), F_id))
    fqF2@id = BStringSet(paste0(id(fqF2), R_id))
    
    
    if (sum(inc1) > (sample_size_in_chunk + n_missing)) {
        smp = sort(sample(sum(inc1), sample_size_in_chunk + n_missing))
        n_missing = 0  ## this was to correct rounding error
        writeFun(fqF1[inc1][smp], file = f1out, mode = "a")
        writeFun(fqF2[inc2][smp], file = f2out, mode = "a")
        nfrq1 = alphabetByCycle(sread(fqF1[inc1][smp]))
        nfrq2 = alphabetByCycle(sread(fqF2[inc2][smp]))
    } else {
        writeFun(fqF1[inc1], file = f1out, mode = "a")
        writeFun(fqF2[inc2], file = f2out, mode = "a")
        nfrq1 = alphabetByCycle(sread(fqF1[inc1]))
        nfrq2 = alphabetByCycle(sread(fqF2[inc2]))
    }
    nucleotideFrequenciesForward = matricesSum(nucleotideFrequenciesForward, nfrq1)
    nucleotideFrequenciesReverse = matricesSum(nucleotideFrequenciesReverse, nfrq2)
    
}
if( !is.na(opt$png_file_output)){
    png(opt$png_file_output, width = 1000, height = 1000)
    par(mfrow=c(2,1))
    plotFrequencies(nucleotideFrequenciesForward)
    mtext("forward reads")
    plotFrequencies(nucleotideFrequenciesReverse)
    mtext("reverse reads")
    dev.off()
}
close(f1)
close(f2)
 
