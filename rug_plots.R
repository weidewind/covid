library("optparse")
library(plyr)
library(lubridate)
library(ggplot2)
library(ggalt)

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
              help="gentag to start output file names with ", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

important_lin = read.csv2(opt$important_lineages, header = T, sep = ",", stringsAsFactors=FALSE) 
pangolined = read.csv2(opt$pangolined_file, header = T, sep = ",")
cluster_strains = read.csv2(file = opt$entry_strains_file,sep = "\t", header = F, stringsAsFactors=FALSE)
print(head(cluster_strains))
cluster_strains = data.frame(entry=as.character(cluster_strains$V1), strains=as.character(cluster_strains$V2), strain_count = sapply(cluster_strains$V2, function(s){length(unlist(strsplit(s,";")))}))
cluster_dates = read.csv2(file = opt$cluster_dates_file,header = T, sep = "\t", stringsAsFactors=FALSE)
cluster_dates = data.frame(entry = cluster_dates$entry, min_cluster_date = as.Date(cluster_dates$min_date, format = "%Y-%m-%d"), max_cluster_date = as.Date(cluster_dates$max_date, format = "%Y-%m-%d"))
cordates_df <- read.csv2(opt$corrected_dates_file, stringsAsFactors=FALSE)
entrydates = read.csv2(file = opt$entry_dates_file,sep = "\t", header = T, stringsAsFactors=FALSE)


cordates_df <- data.frame(entry = cordates_df$entry, 
                          median = as.Date(cordates_df$median, format = "%Y-%m-%d"), mean = as.Date(cordates_df$mean, format = "%Y-%m-%d"),
                          min = as.Date(cordates_df$min, format = "%Y-%m-%d"),
                          mean_all =as.Date(cordates_df$mean_all, format = "%Y-%m-%d"), 
                          date_count = cordates_df$date_count, min_foreign_date = as.Date(cordates_df$min_foreign_date, format = "%Y-%m-%d"), 
                          corrected_median_all =as.Date(cordates_df$corrected_median_all, format = "%Y-%m-%d"),
                          corrected_mean_all =as.Date(cordates_df$corrected_mean_all, format = "%Y-%m-%d"),
                          corrected_grouped_mean_all = as.Date(cordates_df$corrected_grouped_mean_all, format = "%Y-%m-%d"),
                          min_foreign_date = as.Date(cordates_df$min_foreign_date, format = "%Y-%m-%d"),
                          median_foreign_date =  as.Date(cordates_df$median_foreign_date, format = "%Y-%m-%d"))



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
###



fortnight_plot<-function(dd, pngpath){
		hplot = ggplot(dd,aes(corrected_grouped_mean_all, fill = major_lineage_pruned)) +
		  geom_histogram(breaks = by_week(dd$corrected_grouped_mean_all)) +
		  scale_x_date(labels = scales::date_format("%Y-%b-%d"),
					   breaks = by_week(dd$corrected_grouped_mean_all,4), limits = c(min(dd$corrected_grouped_mean_all), max(dd$corrected_grouped_mean_all))) + labs(title = "", x = "", y = "Число завозов", fill = "Pangolin lineage\n")+
		  theme_minimal() +scale_fill_viridis_d(option = "plasma", direction = -1)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		  ggsave(hplot, file=pngpath, width=6, height=4)
}

rug_plot<-function(dates, pngpath, value){
		dates$entry_sorted <- forcats::fct_reorder(dates$entry, as.Date(dates[[value]], format = "%Y-%m-%d"), .desc = TRUE)
		p<-ggplot(dates, 
			    aes(x = get(value), xend = max_cluster_date, y=entry_sorted, color = log(strain_count))) +  
			geom_dumbbell(size = 1.2,
						size_x = 3, 
						size_xend = 3,
						colour = "gray",
						colour_xend = "gray",
						colour_x = "blue") +
			geom_point(aes(x=median, y=entry_sorted), color = "lightblue") +
			geom_point(aes(x=median_foreign_date, y=entry_sorted), color = "green") +
			geom_point(aes(x=min_cluster_date, y=entry_sorted), color = "red", size = 3)+
			geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(x=as.Date("2019-12-15", format = "%Y-%m-%d"), y=entry_sorted, label="strain count", vjust=-0.5),size=1, fontface = "bold") +
			geom_text(aes(label= strain_count, y=entry_sorted,  x= as.Date("2019-12-15", format = "%Y-%m-%d"),  size=1, fontface = "bold", family="Calibri")) +
		   # geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(x=as.Date("2019-11-15", format = "%Y-%m-%d"), y=entry_sorted, label="foreign dates count",size=2, vjust=-0.5)) +
		  #  geom_text(aes(label= date_count, y=entry_sorted,  x= as.Date("2019-11-15", format = "%Y-%m-%d"),  size=2, family="Calibri"), color="black") +
			#geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(label= "dist to closest foreign", y=entry_sorted,  x= as.Date("2021-04-01", format = "%Y-%m-%d"), size=2, vjust=-0.5)) +
			#geom_text(aes(label= round(as.numeric(tree_distance_to_closest_foreign_strain), 6), y=entry_sorted,  x= as.Date("2021-04-01", format = "%Y-%m-%d"),  size=3, family="Calibri"),  color="black") +
			theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			scale_y_discrete(expand=expansion(add=3)) 
			ggsave(p, file=pngpath, width=25, height=85, limitsize = FALSE)
}


stem_dd = dd[dd$tree_distance_to_closest_foreign_strain == 0 & dd$tree_distance_to_closest_rus_strain == 0, ]
print("stem_dd:")
print(head(stem_dd))
nonstem_dd = dd[!(dd$tree_distance_to_closest_foreign_strain == 0 & dd$tree_distance_to_closest_rus_strain == 0), ]
print("nonstem_dd:")
print(str(nonstem_dd))
print(head(nonstem_dd))
stem_pngpath=paste(c(opt$output_folder, "/", opt$tag, "_entries_per_fortnight_pangolin_colored_stem.png"), collapse = "")
nonstem_pngpath=paste(c(opt$output_folder, "/", opt$tag, "_entries_per_fortnight_pangolin_colored_nonstem.png"), collapse = "")

print("Drawing pangolin-colored barplot with entries per fortnight..")
fortnight_plot(stem_dd, pngpath=stem_pngpath)
fortnight_plot(nonstem_dd, pngpath=nonstem_pngpath)

print("Drawing the rug..")
rug_plot(stem_dd, paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_stem.png"),collapse=""), "corrected_grouped_mean_all")
rug_plot(nonstem_dd, paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_nonstem.png"),collapse=""), "corrected_grouped_mean_all")


