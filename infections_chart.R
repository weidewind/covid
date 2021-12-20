library("optparse")
library(plyr)
library(lubridate)
library("forcats")
library(ggplot2)
#library(ggalt)

option_list = list(
			make_option(c("-i", "--input"), type="character", default="",
              help="path to file", metavar="character"),
			make_option(c("-o", "--outfolder"), type="character", default="", 
              help="", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)
print("Reading input..")
df = read.csv2(opt$input, header = T, sep = "\t")
print(head(df))
df$date = as.Date(df$date, format = "%d.%m.%Y")

dfsum = ldply(unique(df$date), function(d){
	data.frame(date = d, infections = sum(df[df$date == d, "infections"]))
})
print(head(dfsum))


by_week <- function(x,n=1){
  seq(min(x,na.rm=T),max(x,na.rm=T),by=paste0(n," weeks"))
}

hplot = ggplot(df,aes(x=date)) +
	scale_x_date(date_breaks = "1 month", limits = c(as.Date("2020-02-01", format = "%Y-%m-%d"), max(df$date))) +
	geom_freqpoly(breaks = by_week(df$date), aes(weight = infections))+ labs(title = "", x = "", y = "Infections per week",stat="identity") +
	theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(hplot, file=paste(c(opt$outfolder, "/", "infections.png"), collapse=""),  width=15, height=5)
