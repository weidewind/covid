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
			make_option(c("-o", "--output_folder"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-d", "--entry_dates_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-t", "--tag"), type="character", default="", 
              help="gentag to start output file names with ", metavar="character"),
			make_option(c("-b", "--branchlength"), type="character", default="", 
              help="subs or subs_per_site? ", metavar="character")
			make_option(c("-s", "--seqlength"), type="character", default="", 
              help="if branchlength is in subs, we need seqlength", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

clock_rate = 0.0011 # sub per site per year
if (opt$branchlength == "subs"){
	seqlength = opt$seqlength
}else if (opt$branchlength == "subs_per_site"){
	seqlength = 1
}else{
	stop(paste(c("--branchlength must be either 'subs' or 'subs_per_site'! Got ", opt$branchlength), collapse=""))
}

entrydates = read.csv2(file = opt$entry_dates_file,sep = "\t", header = T, stringsAsFactors=FALSE)


mean_all = apply(entrydates, 1, function(row){
 as.character(mean.Date(as.Date(c(substr(as.character(unlist(strsplit(row["close_foreign_dates"], ","))), 0,10),substr(as.character(unlist(strsplit(row["close_rus_dates"], ","))), 0,10)), format = "%Y-%m-%d")))
})

# Estimate entry date by close foreign strains dates and lengths of corresponding branches
# Columns close_foreign_dates, tree_distance_foreign_to_mrca, tree_distance_foreign_mrca_to_entry contain equally ordered ","-separated data,
# thus we can convert each row to dataframe, where one row corresponds to a single foreign strain. 

print("Calculating corrected dates..")

corrected_dates = apply(entrydates, 1, function(row){
	print(row["entry"])
      df = data.frame(ftm = unlist(strsplit(row["tree_distance_foreign_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_foreign_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_foreign_dates"], ","))))
      df = df[!is.na(df$fd),]
      # Now from each foreign strain (ie, for each row in the dataframe df) we compute the date of entry
	  print(head(df))
      r <- apply(df, 1, function(r){
        res = ymd(r["fd"]) + dyears((as.numeric(r["mte"])- as.numeric(r["ftm"]))/(seqlength*clock_rate))
      # convert to character, because otherwise it is auto-converted to numeric  
        substr(as.character(res), 0,10)
      })
      names(r) <-NULL
      r
})

if(is.matrix(corrected_dates)){
	corrected_dates = split(corrected_dates, rep(1:ncol(corrected_dates), each = nrow(corrected_dates)))
}
print("Corrected dates:")
print(str(corrected_dates))
#print(head(corrected_dates))

# Estimate entry date by close foreign strains AND close rus strains dates and lengths of corresponding branches
print("Calculating all corrected dates..")
all_corrected_dates = apply(entrydates, 1, function(row){
  print(row["entry"])

  df = data.frame(ftm = unlist(strsplit(row["tree_distance_foreign_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_foreign_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_foreign_dates"], ","))))
  df = df[!is.na(df$fd) ,]
  # Now from each foreign strain (ie, for each row in the dataframe df) we compute the date of entry
  fr <- apply(df, 1, function(r){
	res = ymd(r["fd"]) + dyears((as.numeric(r["mte"])- as.numeric(r["ftm"]))/(seqlength*clock_rate))
  # convert to character, because otherwise it is auto-converted to numeric  
	substr(as.character(res), 0,10)
  })
  names(fr) <-NULL

  rdf = data.frame(ftm = unlist(strsplit(row["tree_distance_rus_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_rus_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_rus_dates"], ","))))
  rdf = rdf[!is.na(rdf$fd) ,]
  # Now from each foreign strain (ie, for each row in the dataframe df) we compute the date of entry
  rr <- apply(rdf, 1, function(r){
	res = ymd(r["fd"]) + dyears((as.numeric(r["mte"])- as.numeric(r["ftm"]))/(seqlength*clock_rate))
  # convert to character, because otherwise it is auto-converted to numeric  
	substr(as.character(res), 0,10)
  })
  names(rr) <-NULL
  
  c(fr,rr)
})

if(is.matrix(all_corrected_dates)){
	all_corrected_dates = split(all_corrected_dates, rep(1:ncol(all_corrected_dates), each = nrow(all_corrected_dates)))
}
print("All Corrected dates:")
print(str(all_corrected_dates))
#print(head(all_corrected_dates))

# Estimate entry date SEPARATELY by close foreign strains dates, corrected by branch lengths, and close rus strains, corrected by branch length.
# This results in 2 median values; return their mean
print("Calculating all grouped corrected dates..")
all_grouped_corrected_mean = apply(entrydates, 1, function(row){
  print(row["entry"])
  frm = NULL  

      df = data.frame(ftm = unlist(strsplit(row["tree_distance_foreign_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_foreign_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_foreign_dates"], ","))))
      df = df[!is.na(df$fd) ,]
      # Now from each foreign strain (ie, for each row in the dataframe df) we compute the date of entry
      fr <- apply(df, 1, function(r){
        res = ymd(r["fd"]) + dyears((as.numeric(r["mte"])- as.numeric(r["ftm"]))/(seqlength*clock_rate))
        
      # convert to character, because otherwise it is auto-converted to numeric  
        substr(as.character(res), 0,10)
      })
      names(fr) <-NULL


      rdf = data.frame(rtm = unlist(strsplit(row["tree_distance_rus_to_mrca"], ",")), mte = unlist(strsplit(row["tree_distance_rus_mrca_to_entry"], ",")), fd = ymd(unlist(strsplit(row["close_rus_dates"], ","))))
      rdf = rdf[!is.na(rdf$rdf),]
      rr <- apply(rdf, 1, function(r){
        res = ymd(r["fd"]) + dyears((as.numeric(r["mte"])- as.numeric(r["ftm"]))/(seqlength*clock_rate))
        substr(as.character(res), 0,10)
      })
      names(rr) <-NULL
      
      frm = median(as.Date(fr,format = "%Y-%m-%d"), na.rm = TRUE) 
      rrm = median(as.Date(rr,format = "%Y-%m-%d"), na.rm = TRUE) 
      as.character(mean(c(frm, rrm),na.rm = TRUE))
})

print("All Grouped Corrected mean:")
print(str(all_grouped_corrected_mean))
print(head(all_grouped_corrected_mean))

all_corrected_medians<-sapply(all_corrected_dates,  function(list){
  if (length(list) > 0){
  as.character(median(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})

all_corrected_means<-sapply(all_corrected_dates,  function(list){
  if (length(list) > 0){
  as.character(mean.Date(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})

print("All corrected means:")
print(str(all_corrected_means))


means<-sapply(corrected_dates,  function(list){
  if (length(list) > 0){
  as.character(mean.Date(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})
print("corrected_dates means:")
print(str(means))


medians<-sapply(corrected_dates, function(list){
  if (length(list) > 0){
  as.character(median(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})
print("corrected_dates medians:")
print(str(medians))


mins<-sapply(corrected_dates, function(list){
  if (length(list) > 0){
  as.character(min(as.Date(list, format = "%Y-%m-%d"), na.rm = TRUE))
  }else{
    NA
  }
})
  
counts<-sapply(all_corrected_dates, length)


cordates_df <- data.frame(entry = entrydates$entry, 
                          median = as.Date(medians, format = "%Y-%m-%d"), mean = as.Date(means, format = "%Y-%m-%d"),
                          min = as.Date(mins, format = "%Y-%m-%d"),
                          mean_all =as.Date(mean_all, format = "%Y-%m-%d"), 
                          date_count = counts, min_foreign_date = as.Date(entrydates$min_foreign_date, format = "%Y-%m-%d"), 
                          corrected_median_all =as.Date(all_corrected_medians, format = "%Y-%m-%d"),
                          corrected_mean_all =as.Date(all_corrected_means, format = "%Y-%m-%d"),
                          corrected_grouped_mean_all = as.Date(all_grouped_corrected_mean, format = "%Y-%m-%d"),
                          date_count = counts, min_foreign_date = as.Date(entrydates$min_foreign_date, format = "%Y-%m-%d"),
                          median_foreign_date =  as.Date(entrydates$median_foreign_date, format = "%Y-%m-%d"))

# write it down for later use
print("Writing corrected entry dates to file..")
write.csv2(cordates_df, paste(c(opt$entry_dates_file, ".corrected.", opt$tag, ".csv"), collapse = ""), quote = F, row.names = F)

