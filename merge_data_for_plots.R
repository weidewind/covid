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
			make_option(c("--coverage"), type="numeric", default=NULL, help="only used for tagging"),
			make_option(c("-w", "--entry_weight_file"), type="character", default=NULL, 
              help="path to file produced by estimate_coverage.R", metavar="character"),
			make_option(c("--dates_stats_file"), type="character", default=NULL, 
              help="", metavar="character"),
			make_option(c("--all_pangolined_file"), type="character", default=NULL, 
              help="", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

important_lin = read.csv2(opt$important_lineages, header = T, sep = ",", stringsAsFactors=FALSE, comment.char = '@') 
pangolined = read.csv2(opt$pangolined_file, header = T, sep = ",")
cluster_strains = read.csv2(file = opt$entry_strains_file,sep = "\t", header = F, stringsAsFactors=FALSE)
print(head(cluster_strains))
cluster_strains = data.frame(entry=as.character(cluster_strains$V1), strains=as.character(cluster_strains$V2), previous_export=as.character(cluster_strains$V3), strain_count = sapply(cluster_strains$V2, function(s){length(unlist(strsplit(s,";")))}))
cluster_dates = read.csv2(file = opt$cluster_dates_file,header = T, sep = "\t", stringsAsFactors=FALSE)
cluster_dates = data.frame(entry = cluster_dates$entry, min_cluster_date = as.Date(cluster_dates$min_date, format = "%Y-%m-%d"), max_cluster_date = as.Date(cluster_dates$max_date, format = "%Y-%m-%d"))
cordates_df <- read.csv2(opt$corrected_dates_file, stringsAsFactors=FALSE)
entrydates = read.csv2(file = opt$entry_dates_file,sep = "\t", header = T, stringsAsFactors=FALSE)

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


dstats = read.table(opt$dates_stats_file, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
pango = read.table(opt$all_pangolined_file, header = T, stringsAsFactors=FALSE, sep = ",")

foreign_pango_stat = as.data.frame(do.call(rbind, lapply(1:nrow(dstats), function(x){
	entry = dstats[x, "entry"]
	foreigners = unlist(strsplit(dstats[x, "close_foreign_strains"],","))
	pangolins = sapply(foreigners, function(f){
		pango[pango$taxon == f,"lineage"][1]
	})
	st = apply(count(pangolins),1, function(r){paste(c(r[1],":",r[2]),collapse = "")})
	c(entry = entry, foreign_pango_stat = paste(st, collapse = ";"))
})))

dates <-merge(cordates_df, cluster_dates, by = "entry")
dates <-merge(dates, foreign_pango_stat, by= "entry")
dates <-merge(dates, cluster_strains, by= "entry")
dates <-merge(dates, data.frame(entry=entrydates$entry, 
	tree_distance_to_closest_foreign_strain = as.numeric(entrydates$tree_distance_to_closest_foreign_strain), 
	tree_distance_to_closest_rus_strain = as.numeric(entrydates$tree_distance_to_closest_rus_strain)), by= "entry")


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

write.table(dd, file = opt$output, row.names = F, quote = F, sep="\t")


