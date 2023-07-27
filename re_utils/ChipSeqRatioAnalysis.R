#!/usr/bin/env Rscript
library(R2HTML, quietly=T)
library(base64enc, quietly=T)


htmlheader=
"	<html xmlns:mml=\"http://www.w3.org/1998/Math/MathML\">
  <head>
  <title> ChIP-Seq Mapper Output </title>
 <style>
html,body{font-family:Verdana,sans-serif;font-size:15px;line-height:1.5}

table {
  border-collapse: collapse;
  border: 1px solid black;
  width: 1000pt
}
table, th, td {
  border: 1px solid black;
}
</style>
  
  </head>



"


                                        #arguments
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
HTMLfile <- args[2]
thr = 5
threshld <- thr/(thr+1)

inputN=as.numeric(args[3])
chipN=as.numeric(args[4])
                                        #dataframe preprocessing and table creation
df <- read.delim(input, comment.char="#")

df$"Ratio Chip/Input"=df$Chip_Hits/df$Input_Hits
df$"Normalized ratio Chip/Input"=(df$Chip_Hits/chipN)/(df$Input_Hits/inputN)

df$"Ratio Chip/(Chip+Input)"=df$Chip_Hits/(df$Chip_Hits + df$Input_Hits)
df$"Normalized ratio Chip/(Chip+Input)"=(df$Chip_Hits/chipN)/((df$Input_Hits/inputN)+(df$Chip_Hits/chipN))

outputTable = df[df$"Normalized ratio Chip/(Chip+Input)" > threshld,
                 ]
outputTable = outputTable[!is.na(outputTable$Cluster),
                          c('Cluster',	'Chip_Hits',	'Input_Hits',
                            'Normalized ratio Chip/Input','Normalized ratio Chip/(Chip+Input)')]
save.image("tmp.RData")                                        #Plot creation
pngfile <- tempfile()
png(pngfile, width = 1000, height = 1200, pointsize=20)
par(mfrow=c(2,1))
lims=range(df$"Normalized ratio Chip/Input"[df$"Normalized ratio Chip/Input">0], finite = TRUE)
suppressWarnings(plot(df$Cluster,df$"Normalized ratio Chip/Input", log="y", xlab="Cluster Nr.", ylab="Normalized ChiP/Input ratio", pch=20, ylim=lims))
abline(h=1,col='#00000080', lwd = 2)
abline(h=thr,col='#FF000080', lwd = 2)


suppressWarnings(plot(df$Cluster,df$"Normalized ratio Chip/(Chip+Input)", xlab="Cluster Nr.", ylab="Normalized Chip/(Chip+Input)", pch=20))
abline(h=0.5,col='#00000080', lwd = 2)
abline(h=threshld,col='#FF000080', lwd = 2)


dev.off()
graph <- paste('<img src="data:image/png;base64 ,',
               base64encode(pngfile),
               '" alt="image" />'
)

                                        #HMTL report creation + writing final output
directory=dirname(HTMLfile)
filename=basename(HTMLfile)
## create HTML header
cat(htmlheader, file = filename)


HTML(graph, file=filename)
if (nrow(outputTable)>0){
  HTML(outputTable, file=filename, classtable = "dataframe",
       row.names=FALSE, align='left', caption=paste("Clusters with Normalized ChIP/Input ratio >", thr), captionalign="top")
}
HTMLEndFile(filename) 
# file.rename(from=filename, to=HTMLfile)
system(sprintf("cp -r ./%s %s", filename, HTMLfile))
write.table(df, file=input, sep="\t", row.names = FALSE)
