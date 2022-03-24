library("optparse")
library(plyr)
library(lubridate)


option_list = list(
			make_option(c("-o", "--output"), type="character", default="", 
              help="", metavar="character"),
			 make_option(c("-d", "--dates_stats"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-p", "--pangolined"), type="character", default="", 
              help="", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
options(stringsAsFactors = FALSE)


dstats = read.table(opt$dates_stats, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
pango = read.table(opt$pangolined, header = T, stringsAsFactors=FALSE, sep = ",")

output = as.data.frame(do.call(rbind, lapply(1:nrow(dstats), function(x){
	entry = dstats[x, "entry"]
	foreigners = unlist(strsplit(dstats[x, "close_foreign_strains"],","))
	pangolins = sapply(foreigners, function(f){
		pango[pango$taxon == f,"lineage"][1]
	})
	st = apply(count(pangolins),1, function(r){paste(c(r[1],":",r[2]),collapse = "")})
	c(entry = entry, foreign_pango_stat = paste(st, collapse = ";"))
})))


write.table(output, opt$output, sep  = "\t")