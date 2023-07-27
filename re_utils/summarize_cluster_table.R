#!/usr/bin/env Rscript
library(optparse)
option_list <- list( 
  make_option(c("-c", "--cluster_table"), default=NA, type = "character",
              help="file from RepeatExplorer2 clustering - CLUSTER_TABLE.csv"),

  make_option(c("-m", "--comparative_counts"),default = NA,type = "character",
              help="file from RepeatExplorer2 output - COMPARATIVE_ANALYSIS_COUNTS.csv"),
  make_option(c("-o", "--output"), type="character",
              help="output file name")
)


opt = parse_args(OptionParser(option_list = option_list))

## for testing
cluster_annotation = opt$cluster_table
header_line = grep(".*Cluster.*Supercluster.*Size", readLines(cluster_annotation))
annot = read.table(cluster_annotation, sep="\t",header=TRUE,as.is=TRUE, skip = header_line - 1)


input_read_counts = as.numeric(strsplit(
  grep("Number_of_analyzed_reads",
       readLines(con=cluster_annotation, n=header_line),
       value=TRUE)
 ,split="\t")[[1]][2]
)

## complete classification table:
unique_groups = sort(unique(annot$Final_annotatio))

groups_to_remove = grep("contamination|organelle", unique_groups, value=TRUE)
groups_to_keep =  unique_groups[!(unique_groups %in% groups_to_remove)]

if (length(groups_to_remove)>0){
  input_count_reads_corrected = input_read_counts - sum(annot$Size_adjusted[annot$Final_annotation %in% groups_to_remove])

}else{
  input_count_reads_corrected = input_read_counts 
}

proportion = numeric()
sum_of_reads = numeric()
for (g in groups_to_keep){
  sum_of_reads[g] = sum(annot$Size_adjusted[annot$Final_annotation %in% g])
  proportion[g] = sum_of_reads[g] / input_count_reads_corrected
}



summary_table = data.frame(Annotation = groups_to_keep,
                           Number_of_reads = sum_of_reads,
                           "Proportion[%]" = proportion * 100 , check.names = FALSE)

print(opt$output)
write.table(summary_table,file = opt$output,
            row.names = FALSE, col.names = TRUE, sep="\t")
