library("optparse")
library(plyr)
library(lubridate)
library(ggplot2)
library(ggalt)
library(scales) # function percent()


option_list = list(
			make_option(c("--dates_stats"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("--variants_file"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("--output"), type="character", default="", 
              help="", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

dates_stats = read.csv2(file = opt$dates_stats,header = T, sep = "\t", stringsAsFactors=FALSE)
variants = read.csv2(file = opt$variants_file,header = T, sep = ",", stringsAsFactors=FALSE)

melted = ldply(variants$entry, function(e){
	row = dates_stats[dates_stats$entry == e, c("close_foreign_dates","close_foreign_strains")]
	dates = unlist(strsplit(row$close_foreign_dates, ","))
	strains = unlist(strsplit(row$close_foreign_strains, ","))
	data.frame(entry = rep(e,length(strains)), date = dates, strain = strains)
})

write.csv2(file = opt$output, melted, row.names = FALSE, quote = FALSE)
