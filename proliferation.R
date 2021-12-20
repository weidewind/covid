library("optparse")
library(plyr)
library(lubridate)
library(ggplot2)


option_list = list(
			make_option(c("-i", "--important_lineages"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("-e", "--melted_entries_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-o", "--output"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-w", "--window"), type="numeric", default="", 
              help="", metavar="numeric")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

print("entry_strains_file:")
print(opt$entry_strains_file)

important_lin = read.csv2(opt$important_lineages, header = T, sep = ",", stringsAsFactors=FALSE) 

df = read.csv2(opt$melted_entries_file, header = T, stringsAsFactors=FALSE, colClasses = c("numeric", "character", "character", "character", "Date", "logical", "Date"))

strains_per_entry <- function(df, window, date, lineage){
	entries = df[df$date >= date - window/2 & df$date <= date + window/2 & df$lineage == lineage & !df$is_stem, "entry"]
	entries_per_date = length(unique(df[df$entry_date == date & df$lineage == lineage & !df$is_stem,"entry"])) # how many entries happened on that day? entry_count - how many entries are alive on that day?
	data.frame(date = date, entries_per_date = entries_per_date, strain_count = length(entries), entry_count = length(unique(entries)), strains_per_entry = length(entries)/length(unique(entries)))
}

window = opt$window
print("window/2:")
print(window/2)
print("mindate")
print(min(df$date, na.rm = TRUE))
print("will start from")
print(min(df$date, na.rm = TRUE)+window/2)

sapply(important_lin$lineage, function(lin){
	spe = ldply(seq(min(df$date, na.rm = TRUE)+window/2, max(df$date, na.rm = TRUE)-window/2, 1), function(d){
		strains_per_entry(df, window = window, date = d, lineage = lin)
	})
	p = ggplot(spe, aes(x=date,y=strains_per_entry))+geom_line() + geom_line(aes(y=entry_count), color = "blue") + geom_line(aes(y=strain_count), color = "gray") + geom_line(aes(y=entries_per_date), color = "red") +
	theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_date(date_breaks = "3 months") +
	labs(title = "", x = "", y = paste(c("Число сэмплов за ", window, " дней/число завозов"), collapse=""))
	ggsave(p, file=paste(c(opt$output, "_", lin, "_", window, "days.nonstem.png"), collapse = ""), width=15, height=5, limitsize = FALSE)
	write.csv2(spe, paste(c(opt$output, "_", lin, "_", window, "days.nonstem.csv"), collapse = ""), row.names = F, quote = F)
})

