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
			make_option(c("-t", "--title"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-o", "--outfolder"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-d", "--dates"), type="character", default="", 
              help="leaf_dates.csv", metavar="character"),
			 make_option(c("-e", "--exact"), action = "store_true", default = FALSE)
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)
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
counts = sapply(all_lin, function(l){nrow(pangolined[pangolined$lineage == l,])})
print(head(counts))
print("Deciding what is important..")
desc = unlist(sapply(all_lin, function(l){
  d = important_lin[important_lin$lineage == l, "desc"][1]
  if (is.na(d)){d = ""}
  d
  }))
  
counts_df = data.frame(lineage = names(counts), count = counts,  desc = desc ) #names(counts) %in% important_lin$lineage 

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

counts_df$important = is_important(names(counts), opt$exact)

print(head(counts_df))
ordered_counts_df = counts_df[order( counts_df[,"important"], counts_df[,"count"], decreasing = T ),]
write.table(ordered_counts_df, paste(c(opt$outfolder,"/",opt$title,"_counts.csv"),collapse = ""), sep="\t", quote = F, row.names = F)

df = data.frame(pangolined, important = is_important(pangolined$lineage, opt$exact))
 pangolined_counts = sapply(pangolined$lineage, function(lin){
	 counts_df[counts_df$lineage == lin, "count"]
 })
df$count = pangolined_counts
print(head(df))

if(opt$exact){
	extag = "exact"
}else{
	extag = "generalized"
}
# print("Plotting..")
# p = ggplot(df, aes(lineage, fill = important)) + scale_fill_manual(values=c("darkblue", "darkred"))+ geom_bar(aes(x=forcats::fct_infreq(lineage)))  + theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_text(data = df[df$important == TRUE,], stat = "count", aes(label = after_stat(count)), vjust = -1, hjust = 0) + ggtitle(opt$title)
# output = paste(c(opt$outfolder,"/",opt$title, "_", extag, "_barplot.png"),collapse = "")
# ggsave(p, file=output, width=15, height=10)

# p = ggplot(df[df$count>5 | df$important,], aes(lineage, fill = important)) + scale_fill_manual(values=c("darkblue", "darkred"))+ geom_bar(aes(x=forcats::fct_infreq(lineage)))  + theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='none') + geom_text(data = df[df$important == TRUE,], stat = "count", aes(label = after_stat(count)), vjust = -1, hjust = 0) + ggtitle(opt$title)
# output = paste(c(opt$outfolder,"/",opt$title,"_", extag, "_important_barplot.png"),collapse = "")
# ggsave(p, file=output, width=9, height=6)


by_week <- function(x,n=1){
	seq(min(x,na.rm=T),max(x,na.rm=T),by=paste0(n," weeks"))
}

# legend_lineage = apply(df, 1, function(row){
	# if(row[["important"]]){
		# if ( substr(row[["lineage"]],0,2) =="AY" & !row[["lineage"]] %in% c("AY.4", "AY.12")){
		# return("other delta")} else {return(row[["lineage"]])}
	# }else{return("other")}
# })

get_legend_lineage = function(row, exact){
	if (row[["lineage"]] %in% important_lin$lineage){return(row[["lineage"]])}
	else{
		if(exact){
			return("other")
		}else{
			if(row[["important"]]){
				splitter = unlist(strsplit(row[["lineage"]],".", fixed = TRUE))
				for (i in rev(seq(1,length(splitter)))){
					if( paste(splitter[0:i], collapse=".") %in% important_lin$lineage){
						return(paste(splitter[0:i], collapse="."))
					}
				}
			}else{return("other")}
		}
	}
}

legend_lineage = apply(df, 1, function(row){get_legend_lineage(row, opt$exact)})
df$legend_lineage = legend_lineage



stacked_plot = function(dframe, output, percent = FALSE){ #percent = TRUE does not work as expected
	dframe = dframe[!is.na(dframe$date),]
	
	if(percent){
		hplot = ggplot(dframe,aes(date, fill = legend_lineage)) +
			geom_histogram(breaks = by_week(dframe$date),position = "fill") + 
			scale_y_continuous(labels = percent(c(0, 0.25, 0.5, 0.75, 1))) +
			labs(title = "", x = "", y = "Percentage", fill = "Pangolin lineage\n")
	}else{
		hplot = ggplot(dframe,aes(date, fill = legend_lineage)) +
			geom_histogram(breaks = by_week(dframe$date))+ labs(title = "", x = "", y = "Samples per week", fill = "Pangolin lineage\n")
	}
	hplot = hplot + theme_minimal() +scale_fill_manual(breaks = c(important_lin$lineage, "other"), values=c(important_lin$color,"lightgray"))+
	scale_x_date(date_breaks = "1 month", limits = c(as.Date("2020-02-01", format = "%Y-%m-%d"), max(df$date))) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
	ggsave(hplot, file=output,  width=15, height=5)
}



df_delta = df[df$lineage == "B.1.617.2" | substr(df$lineage,0,2) == "AY",]
output_delta = paste(c(opt$outfolder,"/",opt$title,"_", extag, "_delta_timeseries_stacked.png"),collapse = "") 
stacked_plot(df_delta, output_delta)
output_delta_p = paste(c(opt$outfolder,"/",opt$title,"_", extag, "_delta_timeseries_stacked_percent.png"),collapse = "") 
stacked_plot(df_delta, output_delta_p, percent = TRUE)
#df_major = df[df$count>10,]
df_major = df
output = paste(c(opt$outfolder,"/",opt$title,"_", extag, "_major_timeseries_stacked.png"),collapse = "") 
stacked_plot(df_major, output)
output_p = paste(c(opt$outfolder,"/",opt$title,"_", extag, "_major_timeseries_stacked_percent.png"),collapse = "") 
stacked_plot(df_major, output_p, percent = TRUE)

