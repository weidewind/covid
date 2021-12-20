library("optparse")
library(ggplot2)
library(scales) # function percent()

option_list = list(
			make_option(c("-i", "--infections"), type="character", default=NULL, 
              help="path to file", metavar="character"), # 
			make_option(c("-m", "--melted_entries_file"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("-o", "--output_folder"), type="character", default="", 
              help="", metavar="character")	  
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

df = read.csv2(opt$melted_entries_file, header = T, stringsAsFactors=FALSE, colClasses = c("numeric", "character", "character", "character", "Date", "logical", "Date"))
df = df[!is.na(df$date), ]

infections = read.table(opt$infections, header = T, stringsAsFactors=FALSE, sep = "\t")
