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
              help="", metavar="numeric"),
			make_option(c("-x", "--exact"), action = "store_true", default = FALSE),
			make_option(c("-r", "--restrict"), action = "store_true", default = FALSE,help = "cut off the part of the timeline with no lineages of interest")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

strains_per_entry <- function(df, window, date, lineage, normalized = FALSE){
	entries = df[df$date >= date - window/2 & df$date <= date + window/2 & df$legend_lineage == lineage & !df$is_stem, "entry"]
	entries_per_date = length(unique(df[df$entry_date == date & df$legend_lineage == lineage & !df$is_stem,"entry"])) # how many entries happened on that day? entry_count - how many entries are alive on that day?
	out = data.frame(date = date, entries_per_date = entries_per_date, strain_count = length(entries), entry_count = length(unique(entries)), strains_per_entry = length(entries)/length(unique(entries)))
	if (normalized){
		total_strain_count = nrow(df[df$date >= date - window/2 & df$date <= date + window/2 & !df$is_stem,])
		out$strains_per_entry = out$strains_per_entry/total_strain_count
	}
	return(out)
}

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

get_legend_lineage_delta_aggr = function(row, exact){
	if (substr(row[["lineage"]],0,3) == "AY." | row[["lineage"]] == "B.1.617.2")
		return("delta")
	else
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


is_important = function(lins, exact){
	if(exact){
		res = lins %in% important_lin$lineage
	}else{
		res = sapply(lins, function(lin){
		splitter = unlist(strsplit(lin,".", fixed = TRUE))
		for (i in rev(seq(1,length(splitter)))){
			if( paste(splitter[0:i], collapse=".") %in% important_lin$lineage){
				return(TRUE)
			}
		}
		return(FALSE)
	})
	}
	return(res)
}

important_lin = read.csv2(opt$important_lineages, header = T, sep = ",", stringsAsFactors=FALSE) 
df = read.csv2(opt$melted_entries_file, header = T, stringsAsFactors=FALSE, colClasses = c("numeric", "character", "character", "character", "Date", "logical", "Date"))
df = df[!is.na(df$date), ]

window = opt$window
print("window/2:")
print(window/2)
print("mindate")
print(min(df$date, na.rm = TRUE))
print("will start from")
print(min(df$date, na.rm = TRUE)+window/2)

df$important = is_important(df$lineage, opt$exact)
legend_lineage = apply(df, 1, function(row){get_legend_lineage_delta_aggr(row, opt$exact)})
df$legend_lineage = legend_lineage

plot_multilin = function(df, window, normalized = FALSE){
	multilin_spe = ldply(unique(df$legend_lineage), function(lin){
		spe = ldply(seq(min(df$date, na.rm = TRUE)+window/2, max(df$date, na.rm = TRUE)-window/2, 1), function(d){
			strains_per_entry(df, window = window, date = d, lineage = lin, normalized = normalized)
		})
		data.frame(legend_lineage = lin, spe)
	})
	print("multilin_spe:")
	print(head(multilin_spe))
	y_lab = paste(c("log(Samples per ", window, " days/number of entries)"), collapse="")
	norm_tag = ""
	if (normalized){
		y_lab = paste(c("log(Normalized samples per ", window, " days/number of entries)"), collapse="")
		norm_tag = "normalized."
	}
	mindate = as.Date("2020-02-01", format = "%Y-%m-%d")
	if(opt$restrict){
		mindate = min(df[df$important == TRUE,"date"])
	}
	p = ggplot(multilin_spe, aes(x=date,y=strains_per_entry, color = legend_lineage))+geom_line() +
	theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
	scale_color_manual(breaks = c(important_lin$lineage, "delta", "other"), values=c(important_lin$color,"red", "lightgray")) +
	scale_x_date(date_breaks = "1 month", limits = c(mindate, max(df$date))) +
	scale_y_continuous(trans='log10')+
	labs(title = "", x = "", y = y_lab, color = "Pangolin lineage\n")
	ggsave(p, file=paste(c(opt$output,  "_", window, "days.nonstem.",norm_tag,"png"), collapse = ""), width=15, height=5, limitsize = FALSE)
	write.csv2(multilin_spe, paste(c(opt$output,  "_", window, "days.nonstem.",norm_tag,"csv"), collapse = ""), row.names = F, quote = F)


}

plot_multilin(df=df, window = window, normalized = FALSE)
plot_multilin(df=df, window = window, normalized = TRUE)



