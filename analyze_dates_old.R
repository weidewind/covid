#!/usr/bin/env Rscript

# libproj and/or proj.h/proj_api.h not found in standard search locations
# Install PROJ library and if necessary set PKG_CPPFLAGS/PKG_LIBS accordingly
# Solution:
# devtools::install_github("hrbrmstr/ggalt", ref = "noproj", lib="~/R/library/")

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

clock_rate = 0.0011 # sub per site per year
worst_case_count = 2978

entrydates = read.csv2(file = opt$entry_dates_file,sep = "\t", header = T, stringsAsFactors=FALSE)
important_lin = read.csv2(opt$important_lineages, header = T, sep = ",", stringsAsFactors=FALSE) 
pangolined = read.csv2(opt$pangolined_file, header = T, sep = ",")
cluster_strains = read.csv2(file = opt$entry_strains_file,sep = "\t", header = F, stringsAsFactors=FALSE)
print(head(cluster_strains))
cluster_strains = data.frame(entry=as.character(cluster_strains$V1), strains=as.character(cluster_strains$V2), strain_count = sapply(cluster_strains$V2, function(s){length(unlist(strsplit(s,";")))}))
cluster_dates = read.csv2(file = opt$cluster_dates_file,header = T, sep = "\t", stringsAsFactors=FALSE)
cluster_dates = data.frame(entry = cluster_dates$entry, min_cluster_date = as.Date(cluster_dates$min_date, format = "%Y-%m-%d"), max_cluster_date = as.Date(cluster_dates$max_date, format = "%Y-%m-%d"))



mean_all = apply(entrydates, 1, function(row){
 as.character(mean.Date(as.Date(c(substr(as.character(unlist(strsplit(row["close_foreign_dates"], ","))), 0,10),substr(as.character(unlist(strsplit(row["close_rus_dates"], ","))), 0,10)), format = "%Y-%m-%d")))
})

# Estimate entry date by close foreign strains dates and lengths of corresponding branches
# Columns close_foreign_dates, tree_distance_foreign_to_mrca, tree_distance_mrca_to_entry contain equally ordered ","-separated data,
# thus we can convert each row to dataframe, where one row corresponds to a single foreign strain. 
print("Calculating corrected dates..")
corrected_dates = apply(entrydates, 1, function(row){
	print(row["entry"])
  if (row["number_of_close_foreign_strains"] == worst_case_count){
      r = substr(as.character(unlist(strsplit(row["close_foreign_dates"], ","))), 0,10)
      names(r) <-NULL
      r
  }else{
      df = data.frame(ftm = unlist(strsplit(row["tree_distance_foreign_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_foreign_dates"], ","))))
      df = df[!is.na(df$fd),]
      # Now from each foreign strain (ie, for each row in the dataframe df) we compute the date of entry
      r <- apply(df, 1, function(r){
        res = ymd(r["fd"]) - dyears(as.numeric(r["ftm"])/clock_rate) + dyears(as.numeric(r["mte"])/clock_rate)
      # convert to character, because otherwise it is auto-converted to numeric  
        substr(as.character(res), 0,10)
      })
      names(r) <-NULL
      r
  }
})

# Estimate entry date by close foreign strains AND close rus strains dates and lengths of corresponding branches
print("Calculating all corrected dates..")
all_corrected_dates = apply(entrydates, 1, function(row){
  print(row["entry"])
  fr = NULL  
  if (row["number_of_close_foreign_strains"] == worst_case_count){
      r = substr(as.character(unlist(strsplit(row["close_foreign_dates"], ","))), 0,10)
      names(r) <-NULL
      fr = r
  }else {
      df = data.frame(ftm = unlist(strsplit(row["tree_distance_foreign_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_foreign_dates"], ","))))
      df = df[!is.na(df$fd) ,]
      # Now from each foreign strain (ie, for each row in the dataframe df) we compute the date of entry
      r <- apply(df, 1, function(r){
        res = ymd(r["fd"]) - dyears(as.numeric(r["ftm"])/clock_rate) + dyears(as.numeric(r["mte"])/clock_rate)
      # convert to character, because otherwise it is auto-converted to numeric  
        substr(as.character(res), 0,10)
      })
      names(r) <-NULL
      fr = r
  }
      
      rdf = data.frame(rdf = ymd(unlist(strsplit(row["close_rus_dates"], ","))), rstr = unlist(strsplit(row["close_rus_strains"], ",")))
      rdf = rdf[!is.na(rdf$rdf),]
      rr <- apply(rdf, 1, function(r){
        res =  ymd(r["rdf"]) - dyears(as.numeric(row["tree_distance_to_closest_rus_strain"])/clock_rate)
        substr(as.character(res), 0,10)
      })
      names(rr) <-NULL            
      c(fr,rr)
})

# Estimate entry date SEPARATELY by close foreign strains dates, corrected by branch lengths, and close rus strains, corrected by branch length.
# This results in 2 median values; return take their mean
print("Calculating all grouped corrected dates..")
all_grouped_corrected_mean = apply(entrydates, 1, function(row){
  print(row["entry"])
  frm = NULL  
  if (row["number_of_close_foreign_strains"] == worst_case_count){
      r = substr(as.character(unlist(strsplit(row["close_foreign_dates"], ","))), 0,10)
      names(r) <-NULL
      fr = r
  }else {
      df = data.frame(ftm = unlist(strsplit(row["tree_distance_foreign_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_foreign_dates"], ","))))
      df = df[!is.na(df$fd) ,]
      # Now from each foreign strain (ie, for each row in the dataframe df) we compute the date of entry
      r <- apply(df, 1, function(r){
        res = ymd(r["fd"]) - dyears(as.numeric(r["ftm"])/clock_rate) + dyears(as.numeric(r["mte"])/clock_rate)
        
      # convert to character, because otherwise it is auto-converted to numeric  
        substr(as.character(res), 0,10)
      })
      names(r) <-NULL
      fr = r
  }
      
      rdf = data.frame(rdf = ymd(unlist(strsplit(row["close_rus_dates"], ","))), rstr = unlist(strsplit(row["close_rus_strains"], ",")))
      rdf = rdf[!is.na(rdf$rdf),]
      rr <- apply(rdf, 1, function(r){
        res =  ymd(r["rdf"]) - dyears(as.numeric(row["tree_distance_to_closest_rus_strain"])/clock_rate)
        substr(as.character(res), 0,10)
      })
      names(rr) <-NULL
      
      frm = median(as.Date(fr,format = "%Y-%m-%d"), na.rm = TRUE) 
      rrm = median(as.Date(rr,format = "%Y-%m-%d"), na.rm = TRUE) 
      as.character(mean(c(frm, rrm),na.rm = TRUE))
})


all_corrected_medians<-sapply(all_corrected_dates, function(list){
  if (length(list) > 0){
  as.character(median(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})

all_corrected_means<-sapply(all_corrected_dates, function(list){
  if (length(list) > 0){
  as.character(mean.Date(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})

means<-sapply(corrected_dates, function(list){
  if (length(list) > 0){
  as.character(mean.Date(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})

medians<-sapply(corrected_dates, function(list){
  if (length(list) > 0){
  as.character(median(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})

mins<-sapply(corrected_dates, function(list){
  if (length(list) > 0){
  as.character(min(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})
  
counts<-sapply(corrected_dates, length)

corrected_dates_str  = sapply(corrected_dates, function(list){
  paste(list, collapse = ",")
})
cordates_df <- data.frame(entry = entrydates$entry, corrected_dates = corrected_dates_str, 
                          median = as.Date(medians, format = "%Y-%m-%d"), mean = as.Date(means, format = "%Y-%m-%d"),
                          min = as.Date(mins, format = "%Y-%m-%d"), mean_all =as.Date(mean_all, format = "%Y-%m-%d"), 
                          date_count = counts, min_foreign_date = as.Date(entrydates$min_foreign_date, format = "%Y-%m-%d"), 
                          corrected_median_all =as.Date(all_corrected_medians, format = "%Y-%m-%d"),
                          corrected_mean_all =as.Date(all_corrected_means, format = "%Y-%m-%d"),
                          corrected_grouped_mean_all = as.Date(all_grouped_corrected_mean, format = "%Y-%m-%d"),
                          date_count = counts, min_foreign_date = as.Date(entrydates$min_foreign_date, format = "%Y-%m-%d"),
                          median_foreign_date =  as.Date(entrydates$median_foreign_date, format = "%Y-%m-%d"))

# write it down for later use
print("Writing corrected entry dates to file..")
write.csv2(cordates_df, paste(c(opt$output_folder, "/", opt$tag, "_corrected_entry_dates.csv"), collapse = ""), quote = F, row.names = F)

dates <-merge(cordates_df, cluster_dates, by = "entry")
dates <-merge(dates, cluster_strains, by= "entry")
dates <-merge(dates, entrydates[,c("entry", "tree_distance_to_closest_foreign_strain")])

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

###

#dd = dates
#dd$corrected_grouped_mean_all = lubridate::ymd(dd$corrected_grouped_mean_all)
dd=merge(dates,entries, by = "entry")

print("Drawing pangolin-colored barplot with entries per fortnight..")
hplot = ggplot(dd,aes(corrected_grouped_mean_all, fill = major_lineage_pruned)) +
  geom_histogram(breaks = by_week(dd$corrected_grouped_mean_all)) +
  scale_x_date(labels = scales::date_format("%Y-%b-%d"),
               breaks = by_week(dd$corrected_grouped_mean_all,4), limits = c(min(dd$corrected_grouped_mean_all), max(dd$corrected_grouped_mean_all))) + labs(title = "", x = "", y = "Число завозов", fill = "Pangolin lineage\n")+
  theme_minimal() +scale_fill_viridis_d(option = "plasma", direction = -1)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(hplot, file=paste(c(opt$output_folder, "/", opt$tag, "_entries_per_fortnight_pangolin_colored.png"), collapse = ""), width=6, height=4)


print("Drawing the rug..")
dates$entry_sorted <- forcats::fct_reorder(dates$entry, as.Date(dates$corrected_grouped_mean_all, format = "%Y-%m-%d"), .desc = TRUE)
p<-ggplot(dates, 
       aes(x = corrected_grouped_mean_all, xend = max_cluster_date, y=entry_sorted, color = log(strain_count))) +  
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
	ggsave(p, file=paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date.png"),collapse=""), width=25, height=85, limitsize = FALSE)
