library("optparse")
library(ggplot2)
library(scales) # function percent()

option_list = list(
			make_option(c("-i", "--input"), type="character", default=NULL, 
              help="path to file", metavar="character"), # transmission_lineages.withduplicates.outstats
			make_option(c("-o", "--output_folder"), type="character", default="", 
              help="", metavar="character")	  
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

fs = read.table(file = opt$input, header = FALSE, sep = " ")
stats = data.frame(count = as.numeric(fs[,2]), size = as.numeric(substr(fs[,6], 0, nchar(fs[,6])-1)))
g = ggplot(stats, aes(x = size, y = count*10)) + geom_col() +xlab("clade size") + ylab("count")+
 scale_y_log10(labels = function(x) x/10)+ scale_x_log10()+theme_minimal()
ggsave(g, file=paste(c(opt$output_folder, "/", "clade_size_hist.png"), collapse=""), width=10, height=7)	