library("optparse")
library(plyr)
library(lubridate)
library("forcats")
library(ggplot2)
library(scales) # labels = percent()
#library(ggalt)

option_list = list(
			make_option(c("-i", "--input"), type="character", default="",
              help="path to file", metavar="character"),
			make_option(c("-l", "--lineages_list"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-o", "--outfolder"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-d", "--dates"), type="character", default="", 
              help="leaf_dates.csv", metavar="character"),
			make_option(c("-e", "--exact"), action = "store_true", default = FALSE),
			make_option(c("-f", "--from"), type="character", default="", 
              help="date in yyyy-mm-dd format", metavar="character"),
            make_option(c("-t", "--to"), type="character", default="", 
              help="date in yyyy-mm-dd format", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
options(stringsAsFactors = FALSE)


is_important = function(lins, exact){
	if(exact){
		res = lins %in% important_lin$lineage
	}else{
		res = sapply(lins, function(lin){
		print(lin)
		splitter = unlist(strsplit(lin,".", fixed = TRUE))
		print(rev(seq(1,length(splitter))))
		for (i in rev(seq(1,length(splitter)))){
			print(paste(splitter[0:i], collapse="."))
			if( paste(splitter[0:i], collapse=".") %in% important_lin$lineage){
				return(TRUE)
			}
		}
		return(FALSE)
	})
	}
	return(res)
}


print("Reading dates..")
dates = read.csv2(opt$dates, header = F, sep = "\t")
colnames(dates) = c("taxon", "date")
dates$date = as.Date(dates$date)
print("Reading pangolineages..")
pangolined = read.csv2(opt$input, header = T, sep = ",")
pangolined = pangolined[,c("taxon", "lineage")]
print("Joining dates and pangolineages..")
pangolined = join(pangolined, dates, by = "taxon", match = "first")
print(head(pangolined))
important_lin = read.csv2(opt$lineages_list, header = T, sep = ",", comment.char = '@') 
print(important_lin)
all_lin = unique(pangolined$lineage)
print("Counting samples..")
from = as.Date(opt$from, format = "%Y-%m-%d")
to = as.Date(opt$to, format = "%Y-%m-%d")
counts = sapply(all_lin, function(l){nrow(pangolined[pangolined$lineage == l & pangolined$date>=from & pangolined$date <=to,])})
print(head(counts))
print("Deciding what is important..")
desc = unlist(sapply(all_lin, function(l){
  d = important_lin[important_lin$lineage == l, "desc"][1]
  if (is.na(d)){d = ""}
  d
  })) 
counts_df = data.frame(lineage = names(counts), count = counts,  desc = desc ) #names(counts) %in% important_lin$lineage 
counts_df$important = is_important(names(counts), opt$exact)
ordered_counts_df = counts_df[order( counts_df[,"important"], counts_df[,"count"], decreasing = T ),]
write.table(ordered_counts_df, paste(c(opt$outfolder,"/","from_",as.character(from), "_to_", as.character(to),"_counts.csv"),collapse = ""), sep="\t", quote = F, row.names = F)







