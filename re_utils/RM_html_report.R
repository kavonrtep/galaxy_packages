#!/usr/bin/env Rscript
### this script is expected to run from clustering directory! ######

## assume RM-custom_output_tablesummary.csv file in active directory

suppressPackageStartupMessages(library(R2HTML))
######################################################################################
htmlheader="
		<html xmlns:mml=\"http://www.w3.org/1998/Math/MathML\">
		<head>
		<title> Clustering summary </title>
		<style>
		<!--
		table { background:#FFFFFF;
		border:1px solid gray;
		border-collapse:collapse;
		color:#fff;
		font:normal 10px verdana, arial, helvetica, sans-serif;
		}
		caption { border:1px solid #5C443A;
		color:#5C443A;
		font-weight:bold;
		font-size:20pt
		padding:6px 4px 8px 0px;
		text-align:center;
		
		}
		td, th { color:#363636;
		padding:.4em;
		}
		tr { border:1px dotted gray;
		}
		thead th, tfoot th { background:#5C443A;
		color:#FFFFFF;
		padding:3px 10px 3px 10px;
		text-align:left;
		text-transform:uppercase;
		}
		tbody td a { color:#3636FF;
		text-decoration:underline;
		}
		tbody td a:visited { color:gray;
		text-decoration:line-through;
		}
		tbody td a:hover { text-decoration:underline;
		}
		tbody th a { color:#3636FF;
		font-weight:normal;
		text-decoration:none;
		}
		tbody th a:hover { color:#363636;
		}
		tbody td+td+td+td a { background-image:url('bullet_blue.png');
		background-position:left center;
		background-repeat:no-repeat;
		color:#FFFFFF;
		padding-left:15px;
		}
		tbody td+td+td+td a:visited { background-image:url('bullet_white.png');
		background-position:left center;
		background-repeat:no-repeat;
		}
		tbody th, tbody td { text-align:left;
		vertical-align:top;
		}
		tfoot td { background:#5C443A;
		color:#FFFFFF;
		padding-top:3px;
		}
		.odd { background:#fff;
		}
		tbody tr:hover { background:#EEEEEE;
		border:1px solid #03476F;
		color:#000000;
		}
		-->
		</style>
		
		</head>
		
		"
######################################################################################
######################################################################################



#basic statistics:
# Number of reads used for clustering

RM=read.table("RM-custom_output_tablesummary.csv",sep="\t",header=TRUE,as.is=TRUE,check.names=FALSE)

#Any hits to RM database?
N=NA

# convert to legible format:
RM2=data.frame(
		'total length [bp]'=RM$All_Reads_Length[c(T,F,F)],
		'number of reads'=RM$All_Reads_Number[c(T,F,F)],
		check.names=FALSE,stringsAsFactors=FALSE
)

RMpart1=RM[c(T,F,F),-c(1:3),drop=FALSE] #counts
RMpart2=RM[c(F,T,F),-c(1:3),drop=FALSE] #percent

RMjoined=list()

for (i in colnames(RMpart1)){
	RMjoined[[i]]=paste(RMpart1[,i],"hits, ",signif(RMpart2[,i],3),"%",sep='')
}



if (ncol(RM)>3){  # not emppty	
	RM2=cbind(cluster=paste("CL",1:nrow(RM2),sep=''),
			RM2,
			"Genome proportion[%]"=signif(RM2$'number of reads'/N*100,3),
			"cumulative GP [%]"=signif(cumsum(RM2$'number of reads'/N*100),3),
			as.data.frame(RMjoined,stringsAsFactors=FALSE))
	
	##### RM2 formating for html output: #####
	##########################################
	bold=RMpart2>3
	for (i in 6:ncol(RM2)){
		rmcol=RM2[,i]
		RM2[,i]=ifelse(bold[,i-5],paste("<b>",rmcol,"</b>",sep=''),rmcol)
	}
	
	# join hits to one  column
	RMstring=character(nrow(RM2))
	for (i in 1:nrow(RM2)){
		x=ifelse(RMpart2[i,]>0,paste(colnames(RM2[,-(1:5),drop=FALSE])," (",RM2[i,-(1:5),drop=FALSE],")",sep=''),"")
		# reorder based on GR
		x=x[order(RMpart2[i,],decreasing=TRUE)]
		
		RMstring[i]=paste(x[x!=''],collapse="<br />")
		if (nchar(RMstring[i])>240){
			RMstring[i]=paste(substring(RMstring[i],1,220),"......",sep='')
		}
		
	}
}else{  # no RM hits
	RM2=cbind(cluster=paste("CL",1:nrow(RM2),sep=''),
			RM2,
			"Genome proportion[%]"=signif(RM2$'number of reads'/N*100,3),
			"cumulative GP [%]"=signif(cumsum(RM2$'number of reads'/N*100),3))
	RMstring=rep("",nrow(RM)/3)
}


# RM2 add link to subpage


RM2=data.frame(RM2[,1:3],'Repeat Masker'=RMstring,check.names=FALSE)


##################################################################################################
####################                              HTML output                                #####
##################################################################################################


htmlout=commandArgs(T)[1]  # full absolute path

cat(htmlheader,file=htmlout)

HTML.title("RepeatMasker search against custom database",file=htmlout,HR=1)

HTML(RM2,file=htmlout,align='left',caption="",captionalign='')
HTMLEndFile(htmlout)

