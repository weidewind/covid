# libproj and/or proj.h/proj_api.h not found in standard search locations
# Install PROJ library and if necessary set PKG_CPPFLAGS/PKG_LIBS accordingly
# Solution:
# devtools::install_github("hrbrmstr/ggalt", ref = "noproj", lib="~/R/library/")

library("optparse")
library(plyr)
library(lubridate)
library(ggplot2)
library(ggalt)
library(scales) # function percent()

option_list = list(
			make_option(c("-i", "--important_lineages"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("-p", "--pangolined_file"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("-o", "--output_folder"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-r", "--corrected_dates_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-d", "--entry_dates_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-e", "--entry_strains_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-c", "--cluster_dates_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-t", "--tag"), type="character", default="", 
              help="gentag to start output file names with ", metavar="character"),
			make_option(c("--coverage"), type="numeric", default=NULL, help="only used for tagging"),
			make_option(c("-w", "--entry_weight_file"), type="character", default=NULL, 
              help="path to file produced by estimate_coverage.R", metavar="character"), 	
			make_option(c("--mutations_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("--melted_entries_file"), type="character", default="", 
              help="", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

important_lin = read.csv2(opt$important_lineages, header = T, sep = ",", stringsAsFactors=FALSE, comment.char = '@') 
pangolined = read.csv2(opt$pangolined_file, header = T, sep = ",")
cluster_strains = read.csv2(file = opt$entry_strains_file,sep = "\t", header = F, stringsAsFactors=FALSE)
print(head(cluster_strains))
cluster_strains = data.frame(entry=as.character(cluster_strains$V1), strains=as.character(cluster_strains$V2), strain_count = sapply(cluster_strains$V2, function(s){length(unlist(strsplit(s,";")))}))
cluster_dates = read.csv2(file = opt$cluster_dates_file,header = T, sep = "\t", stringsAsFactors=FALSE)
cluster_dates = data.frame(entry = cluster_dates$entry, min_cluster_date = as.Date(cluster_dates$min_date, format = "%Y-%m-%d"), max_cluster_date = as.Date(cluster_dates$max_date, format = "%Y-%m-%d"), all_dates = cluster_dates$all_dates)
cordates_df <- read.csv2(opt$corrected_dates_file, stringsAsFactors=FALSE)
entrydates = read.csv2(file = opt$entry_dates_file,sep = "\t", header = T, stringsAsFactors=FALSE)
if(opt$mutations_file){
	mutations_df = read.csv2(file = opt$mutations_file,header = T, sep = ",", stringsAsFactors=FALSE)
}

cordates_df <- data.frame(entry = cordates_df$entry, 
                          median = as.Date(cordates_df$corrected_median_foreign, format = "%Y-%m-%d"), mean = as.Date(cordates_df$corrected_mean_foreign, format = "%Y-%m-%d"),
                          min = as.Date(cordates_df$corrected_min_foreign, format = "%Y-%m-%d"),
                          mean_all =as.Date(cordates_df$mean_all, format = "%Y-%m-%d"), 
                          date_count = cordates_df$date_count, 
                          corrected_median_all =as.Date(cordates_df$corrected_median_all, format = "%Y-%m-%d"),
                          corrected_mean_all =as.Date(cordates_df$corrected_mean_all, format = "%Y-%m-%d"),
                          corrected_grouped_mean_all = as.Date(cordates_df$corrected_grouped_mean_all, format = "%Y-%m-%d"),
                          min_foreign_date = as.Date(cordates_df$min_foreign, format = "%Y-%m-%d"),
                          median_foreign_date =  as.Date(cordates_df$median_foreign, format = "%Y-%m-%d"))



dates <-merge(cordates_df, cluster_dates, by = "entry")
dates <-merge(dates, cluster_strains, by= "entry")
dates <-merge(dates, data.frame(entry=entrydates$entry, 
	tree_distance_to_closest_foreign_strain = as.numeric(entrydates$tree_distance_to_closest_foreign_strain), 
	tree_distance_to_closest_rus_strain = as.numeric(entrydates$tree_distance_to_closest_rus_strain)), by= "entry")

by_month <- function(x,n=1){
  seq(min(x,na.rm=T),max(x,na.rm=T),by=paste0(n," months"))
}

by_week <- function(x,n=1){
  seq(min(x,na.rm=T),max(x,na.rm=T),by=paste0(n," weeks"))
}

## Adding pangolines
print("Adding pangolin data..")
print(str(dates))
lin_stat = sapply(dates$strains, function(s){
  strains = unlist(strsplit(s, ";"))
  lineages = sapply(strains, function(st){
   pangolined[pangolined$taxon == st, "lineage"]
  })
  print("Str:")
  print(str(lineages))
  print("Lineages:")
  print(lineages)
 tbl = as.data.frame(table(lineages))
 kv = apply(tbl, 1, function(row){
   paste(c(row["lineages"], as.numeric(row["Freq"])), collapse = ":")
 })
 paste(kv, collapse = ";")
})

major_lins = sapply(dates$strains, function(s){
    lins = sapply(unlist(strsplit(s, ";")), function(st){
       pangolined[pangolined$taxon == st, "lineage"]
      })
    tbl = as.data.frame(table(lins))
    if (max(tbl$Freq)/sum(tbl$Freq) > 0.5){major_lin = as.character(tbl[which.max(tbl$Freq), "lins"])}else{major_lin = NA}
    major_lin
})
names(major_lins) = NULL

major_lins_pruned = sapply(major_lins, function(l){
  if (l %in% important_lin$lineage){l}else{
      allnom = unlist(strsplit(l, ".", fixed = TRUE))
      parentnom = paste(head(allnom, length(allnom)-1), collapse = "." )
      if (parentnom %in% important_lin$lineage) {parentnom}else{"other"}
    }
})

entries = data.frame(entry = dates$entry, lineage_stat = lin_stat, major_lineage = major_lins, major_lineage_pruned = major_lins_pruned, important = major_lins %in% important_lin$lineage)
row.names(entries) = NULL
dd=merge(dates,entries, by = "entry")
if(!is.null(opt$entry_weight_file)){
	ew = read.csv2(file = opt$entry_weight_file, header = T, stringsAsFactors=FALSE)
	colnames(ew) = c("entry", "entry_weight")
	dd <-merge(dd, ew, by = "entry")
}else{
	dd$entry_weight = rep(1, nrow(dd))
	
}
###



fortnight_plot<-function(dd, pngpath, percent = FALSE, weight = FALSE){
		if(weight == FALSE){
			dd$entry_weight = rep(1, nrow(dd))
		}
		if (percent){
			hplot = ggplot(dd,aes(corrected_grouped_mean_all, fill = major_lineage_pruned)) +
				geom_histogram(breaks = by_week(dd$corrected_grouped_mean_all), position = "fill",aes(weight = entry_weight)) +
				labs(title = "", x = "", y = "Entries percentage", fill = "Pangolin lineage\n")	+
				scale_y_continuous(labels = percent(c(0,0.25,0.5,0.75,1)))
			ggsave(hplot, file=pngpath, width=15, height=5)
			}else{
			hplot = ggplot(dd,aes(corrected_grouped_mean_all, fill = major_lineage_pruned)) +
				geom_histogram(breaks = by_week(dd$corrected_grouped_mean_all),aes(weight = entry_weight)) +
				labs(title = "", x = "", y = "Entries per week", fill = "Pangolin lineage\n")
			}
			hplot = hplot+scale_x_date(date_breaks = "1 month", limits = c(as.Date("2020-02-01", format = "%Y-%m-%d"), max(cluster_dates$max_cluster_date))) +
				theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
				scale_fill_manual(breaks = c(important_lin$lineage, "other"), values=c(important_lin$color,"lightgray"))
			ggsave(hplot, file=pngpath, width=15, height=5)			

}


rug_plot_selection<-function(dates, pngpath, value, mutations = NULL){
		dates$entry_sorted <- forcats::fct_reorder(dates$entry, as.Date(dates[[value]], format = "%Y-%m-%d"), .desc = TRUE)
		if(mutations){
			dates = dates[dates$entry %in% mutations$entry,]
		}
		print("Head dates:")
		print(head(dates))
		print("Head melted_cluster:")
		melted_cluster = melt_cluster_dates(dates)
		melted_cluster$mutated = sapply(melted_cluster)
		print(head(melted_cluster))
		p<-ggplot(dates) +  
			geom_dumbbell(data = dates,aes(x = get(value), xend = max_cluster_date, y=entry_sorted),size = 1.2,
						size_x = 3, 
						size_xend = 3,
						colour = "gray",
						colour_xend = "gray",
						colour_x = "blue") +
			geom_point(data = dates,aes(x=min_cluster_date, y=entry_sorted), color = "red", size = 3)+
			geom_point(data = melted_cluster, aes(x=date, y=entry_sorted), color = "gray", size = 3, alpha = 0.3)+
			geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(x=as.Date("2019-12-15", format = "%Y-%m-%d"), y=entry_sorted, label="strain count", vjust=-0.5),size=1, fontface = "bold") +
			geom_text(data = dates, aes(label= strain_count, y=entry_sorted,  x= as.Date("2019-12-15", format = "%Y-%m-%d"),  size=1, fontface = "bold", family="Calibri")) +
			theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			scale_y_discrete(expand=expansion(add=3)) 
			ggsave(p, file=pngpath, width=25, height=85, limitsize = FALSE)
}



rug_plot<-function(dates, pngpath, value){
		dates$entry_sorted <- forcats::fct_reorder(dates$entry, as.Date(dates[[value]], format = "%Y-%m-%d"), .desc = TRUE)
		print("Head dates:")
		print(head(dates))
		print("Head melted_cluster:")
		melted_cluster = melt_cluster_dates(dates)
		print(head(melted_cluster))
		#color = log(strain_count)
		p<-ggplot(dates) +  
			geom_dumbbell(data = dates,aes(x = get(value), xend = max_cluster_date, y=entry_sorted),size = 1.2,
						size_x = 3, 
						size_xend = 3,
						colour = "gray",
						colour_xend = "gray",
						colour_x = "blue") +
			geom_point(data = dates, aes(x=median, y=entry_sorted), color = "lightblue") +
			geom_point(data = dates,aes(x=median_foreign_date, y=entry_sorted), color = "green") +
			geom_point(data = dates,aes(x=min_cluster_date, y=entry_sorted), color = "red", size = 3)+
			geom_point(data = melted_cluster, aes(x=date, y=entry_sorted), color = "gray", size = 3, alpha = 0.3)+
			geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(x=as.Date("2019-12-15", format = "%Y-%m-%d"), y=entry_sorted, label="strain count", vjust=-0.5),size=1, fontface = "bold") +
			geom_text(data = dates, aes(label= strain_count, y=entry_sorted,  x= as.Date("2019-12-15", format = "%Y-%m-%d"),  size=1, fontface = "bold", family="Calibri")) +
		   # geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(x=as.Date("2019-11-15", format = "%Y-%m-%d"), y=entry_sorted, label="foreign dates count",size=2, vjust=-0.5)) +
		  #  geom_text(aes(label= date_count, y=entry_sorted,  x= as.Date("2019-11-15", format = "%Y-%m-%d"),  size=2, family="Calibri"), color="black") +
			#geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(label= "dist to closest foreign", y=entry_sorted,  x= as.Date("2021-04-01", format = "%Y-%m-%d"), size=2, vjust=-0.5)) +
			#geom_text(aes(label= round(as.numeric(tree_distance_to_closest_foreign_strain), 6), y=entry_sorted,  x= as.Date("2021-04-01", format = "%Y-%m-%d"),  size=3, family="Calibri"),  color="black") +
			theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			scale_y_discrete(expand=expansion(add=3)) 
			ggsave(p, file=pngpath, width=25, height=85, limitsize = FALSE)
}

melt_cluster_dates = function(df_sorted){
	splitted <- strsplit(as.character(df_sorted$all_dates), ",")
	data.frame(entry_sorted = rep.int(df_sorted$entry_sorted, sapply(splitted, length)), date = as.Date(unlist(splitted), format = "%Y-%m-%d"))
}

worst05 = quantile(dd$min_cluster_date-dd$corrected_grouped_mean_all, probs = 0.05, na.rm = TRUE)
print("left 05 quantile is ")
print(worst05)
dd_pruned05 = dd[dd$min_cluster_date-dd$corrected_grouped_mean_all > worst05, ]

stem_dd = dd[dd$tree_distance_to_closest_foreign_strain == 0 & dd$tree_distance_to_closest_rus_strain == 0, ]
print("stem_dd:")
print(head(stem_dd))
nonstem_dd = dd[!(dd$tree_distance_to_closest_foreign_strain == 0 & dd$tree_distance_to_closest_rus_strain == 0), ]
print("nonstem_dd:")
print(str(nonstem_dd))
print(head(nonstem_dd))
pngpath = paste(c(opt$output_folder, "/", opt$tag, "_entries_per_fortnight_pangolin_colored"), collapse = "")
print("Drawing pangolin-colored barplot with entries per fortnight..")
fortnight_plot(stem_dd, pngpath=paste(c(pngpath, "_stem.png"), collapse = ""))
fortnight_plot(nonstem_dd, pngpath=paste(c(pngpath, "_nonstem.png"), collapse = ""))
fortnight_plot(stem_dd, pngpath=paste(c(pngpath, "_stem_percent.png"), collapse = ""), percent = TRUE)
fortnight_plot(nonstem_dd, pngpath=paste(c(pngpath, "_nonstem_percent.png"), collapse = ""), percent = TRUE)

nonstem_dd_pruned05 = dd_pruned05[!(dd_pruned05$tree_distance_to_closest_foreign_strain == 0 & dd_pruned05$tree_distance_to_closest_rus_strain == 0), ]
fortnight_plot(nonstem_dd_pruned05, pngpath=paste(c(pngpath, "_nonstem_pruned05.png"), collapse = ""))
fortnight_plot(dd_pruned05, pngpath=paste(c(pngpath, "_all_pruned05.png"), collapse = ""))
fortnight_plot(nonstem_dd_pruned05, pngpath=paste(c(pngpath, "_nonstem_pruned05_percent.png"), collapse = ""), percent = TRUE)
fortnight_plot(dd_pruned05, pngpath=paste(c(pngpath, "_all_pruned05_percent.png"), collapse = ""), percent = TRUE)

fortnight_plot(dd, pngpath=paste(c(pngpath, "_coverage_", opt$coverage, "_all_with_weights.png"), collapse = ""), weight = TRUE)
fortnight_plot(dd_pruned05, pngpath=paste(c(pngpath, "_coverage_", opt$coverage, "_all_pruned05_with_weights.png"), collapse = ""), weight = TRUE)
fortnight_plot(nonstem_dd_pruned05, pngpath=paste(c(pngpath, "_coverage_", opt$coverage, "_nonstem_pruned05_with_weights.png"), collapse = ""), weight = TRUE)
fortnight_plot(dd, pngpath=paste(c(pngpath, "_all_testnoweights.png"), collapse = ""))



print("Drawing the rug..")
rug_plot(stem_dd, paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_stem.png"),collapse=""), "corrected_grouped_mean_all")
rug_plot(nonstem_dd, paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_nonstem.png"),collapse=""), "corrected_grouped_mean_all")
rug_plot(nonstem_dd_pruned05, paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_nonstem_pruned05.png"),collapse=""), "corrected_grouped_mean_all")
rug_plot(nonstem_dd_pruned05, paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_nonstem_pruned05.png"),collapse=""), "corrected_grouped_mean_all")


