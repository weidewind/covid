# libproj and/or proj.h/proj_api.h not found in standard search locations
# Install PROJ library and if necessary set PKG_CPPFLAGS/PKG_LIBS accordingly
# Solution:
# devtools::install_github("hrbrmstr/ggalt", ref = "noproj", lib="~/R/library/")

library("optparse")
library(plyr)
library(lubridate)
library(ggplot2)
library(ggalt)
library("RColorBrewer")
library(scales) # function percent()

option_list = list(
			make_option(c("-o", "--output_folder"), type="character", default="", 
              help="", metavar="character"),
			 make_option(c("-d", "--data_file"), type="character", default="", 
              help="produced by merge_data_for_plots.R", metavar="character"),
			make_option(c("-t", "--tag"), type="character", default="", 
              help="gentag to start output file names with ", metavar="character"),
			make_option(c("--coverage"), type="numeric", default=NULL, help="only used for tagging"),
			make_option(c("--mutations_file"), type="character", default=NULL, 
              help="", metavar="character"),
			make_option(c("--foreign_mutations_file"), type="character", default=NULL, 
              help="", metavar="character"),
			make_option(c("--melted_entries_file"), type="character", default="", 
              help=",", metavar="character"),
			make_option(c("--districts"), type="character", default="", 
              help=",", metavar="character"),
			make_option(c("--lineage"), type="character", default="", 
              help=",", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

dd = read.table(opt$data_file, stringsAsFactors = FALSE, header = TRUE, sep = "\t", colClasses =c("character","Date","Date","Date","Date","numeric","Date","Date","Date","Date","Date","Date","Date","character","numeric","numeric","numeric","character","character","character","logical","numeric"))
print("DD str:")
print(str(dd))
# sapply(c("median","mean","min","mean_all","date_count","corrected_median_all","corrected_mean_all","corrected_grouped_mean_all","min_foreign_date","median_foreign_date","min_cluster_date","max_cluster_date"), function(colname){
	# dd[[colname]] = as.Date(dd[[colname]], format = "%Y-%m-%d")
# })



melted_entries = read.csv2(opt$melted_entries_file, header = T, stringsAsFactors=FALSE, colClasses = c("numeric", "character", "character", "character", "Date", "logical", "Date"))
if(!is.null(opt$mutations_file) & !is.null(opt$foreign_mutations_file)){
	mutations_df = read.csv2(file = opt$mutations_file,header = T, sep = ",", stringsAsFactors=FALSE)
	mutations_df$is_double_mut = as.factor(mutations_df$is_double_mut)
	f_mutations_df = read.csv2(file = opt$foreign_mutations_file,header = T, sep = ",", stringsAsFactors=FALSE)
	f_mutations_df$is_double_mut = as.factor(f_mutations_df$is_double_mut)
	f_mutations_df$date = as.Date(f_mutations_df$date, format = "%Y-%m-%d")
}


by_month <- function(x,n=1){
  seq(min(x,na.rm=T),max(x,na.rm=T),by=paste0(n," months"))
}

by_week <- function(x,n=1){
  seq(min(x,na.rm=T),max(x,na.rm=T),by=paste0(n," weeks"))
}


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


rug_plot_selected_lineage<-function(dates, pngpath, value, lineage = NULL, districts = NULL){
		dates$entry_sorted <- forcats::fct_reorder(dates$entry, as.Date(dates[[value]], format = "%Y-%m-%d"), .desc = TRUE)
		print("Head dates:")
		print(head(dates))
		print(str(dates))
		melted_cluster = melt_cluster_dates(dates$entry_sorted)
		print("Head melted_cluster:")
		print(head(melted_cluster))
		print("melted_cluster where lineage is NA:")
		print(melted_cluster[is.na(melted_cluster$lineage),])
		if (!is.null(lineage)){
			melted_cluster_lin = melted_cluster[is_lineage(melted_cluster$lineage, lineage),]
			lin_entries = unique(melted_cluster_lin[,"entry"])
			dates = dates[dates$entry %in% lin_entries,]
			print("head dates for lineage:")
			print(head(dates))
			melted_cluster_nonlin = melted_cluster[melted_cluster$entry %in% lin_entries & !is_lineage(melted_cluster$lineage, lineage),] 
			print("Head melted_cluster_lin:")
			print(head(melted_cluster_lin))
		}

		p<-ggplot(dates) +  
			geom_dumbbell(data = dates,aes(x = get(value), xend = max_cluster_date, y=entry_sorted),size = 1,
						size_x = 2, 
						size_xend = 2,
						colour = "gray",
						colour_xend = "gray",
						colour_x = "blue") +
			geom_point(data = dates,aes(x=median_foreign_date, y=entry_sorted), color = "green") +
			geom_hline(yintercept = dates[dates$strain_count >100, "entry_sorted"], color = "red") + 
			geom_count(data = melted_cluster_nonlin, aes(x=date, y=entry_sorted), color = "black", alpha = 0.5) +
			#geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(x=as.Date("2019-12-15", format = "%Y-%m-%d"), y=entry_sorted, label="strain count", vjust=-0.5),size=1, fontface = "bold") +
			geom_text(data = dates, aes(label= strain_count, y=entry_sorted,  x= as.Date("2019-12-15", format = "%Y-%m-%d"), fontface = "bold", family="Calibri")) +
			theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			scale_y_discrete(expand=expansion(add=3)) 
			if (is.null(districts)){
				p<-p+geom_count(data = melted_cluster_lin, aes(x=date, y=entry_sorted, color="gray"), alpha = 0.5)
			}else{
				melted_cluster_lin = merge(melted_cluster_lin, districts, by = "strain", all.x = T)
				print("Head melted_cluster_lin with districts:")
				print(head(melted_cluster_lin))
				p<-p+geom_count(data = melted_cluster_lin, aes(x=date, y=entry_sorted, color=district), alpha = 0.5)
			}
			param = 100/nrow(dates)
			if(param < 1){param = 1}
			ggsave(p, file=pngpath, width=25, height=85/param, limitsize = FALSE)
}


rug_plot_selected_with_foreign<-function(dates, pngpath, value, mutations, foreign_mutations = NULL, districts = NULL){
		dates = dates[dates$entry %in% mutations$entry,]
		dates$entry_sorted <- forcats::fct_reorder(dates$entry, as.Date(dates[[value]], format = "%Y-%m-%d"), .desc = TRUE)
		print("Head dates:")
		print(head(dates))
		print(str(dates))
		print("Head mutations:")
		print(head(mutations))
		melted_cluster = melt_cluster_dates(dates$entry_sorted)
		print("Head melted_cluster:")
		print(head(melted_cluster))
		melted_cluster_nondelta = melted_cluster[!(melted_cluster$strain %in% mutations$taxon),]
		melted_cluster_delta = melted_cluster[melted_cluster$strain %in% mutations$taxon,]
		melted_cluster_delta$mutated = unlist(sapply(melted_cluster_delta$strain, function(s){
			mutations[mutations$taxon == s, "is_double_mut"][1]
		}))
		print("Head melted_cluster_delta:")
		print(head(melted_cluster_delta))
		
		if(!is.null(foreign_mutations)){
			foreign_mutations = resort(entry_sorted = dates$entry_sorted, mydf = foreign_mutations)
			foreign_mutations$mutated = foreign_mutations$is_double_mut
			print("foreign resorted:")
			print(head(foreign_mutations))
			print(str(foreign_mutations))
			non_delta_foreign = foreign_mutations[foreign_mutations$lineage != "B.1.617.2" & substr(foreign_mutations$lineage, 0,2) != "AY",]
			delta_foreign = foreign_mutations[foreign_mutations$lineage == "B.1.617.2" | substr(foreign_mutations$lineage, 0,2) == "AY",]
		}
		p<-ggplot(dates) +  
			geom_dumbbell(data = dates,aes(x = get(value), xend = max_cluster_date, y=entry_sorted),size = 1,
						size_x = 2, 
						size_xend = 2,
						colour = "gray",
						colour_xend = "gray",
						colour_x = "blue") +
			geom_hline(yintercept = dates[dates$strain_count >100, "entry_sorted"], color = "red") + 
			geom_count(data = melted_cluster_nondelta, aes(x=date, y=entry_sorted), color = "black", alpha = 0.5) +
			geom_text(data = dates[dates$entry_sorted == tail(levels(dates$entry_sorted),1),], aes(x=as.Date("2019-12-15", format = "%Y-%m-%d"), y=entry_sorted, label="strain count", vjust=-0.5), fontface = "bold") +
			geom_text(data = dates, aes(label= strain_count, y=entry_sorted,  x= as.Date("2019-12-15", format = "%Y-%m-%d"), fontface = "bold", family="Calibri")) +
			theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			scale_y_discrete(expand=expansion(add=3)) 
			if (!is.null(foreign_mutations)){
				p<-p + geom_count(data = non_delta_foreign, aes(x=date, y=entry_sorted), color = "black", shape = 2, alpha = 0.5)
			}
			if (is.null(districts)){
				p<-p+geom_count(data = melted_cluster_delta, aes(x=date, y=entry_sorted, color=mutated), alpha = 0.5)
			}else{
				melted_cluster_delta = merge(melted_cluster_delta, districts, by = "strain", all.x = T)
				p<-p+geom_count(data = melted_cluster_delta[melted_cluster_delta$mutated == 1,], aes(x=date, y=entry_sorted, color=district), alpha = 0.5)
			}
			ggsave(p, file=pngpath, width=25, height=85, limitsize = FALSE)
}



rug_plot<-function(dates, pngpath,value){
		dates$entry_sorted <- forcats::fct_reorder(dates$entry, as.Date(dates[[value]], format = "%Y-%m-%d"), .desc = TRUE)
		print("Head dates:")
		print(head(dates))
		print("Head melted_cluster:")
		melted_cluster = melt_cluster_dates(dates$entry_sorted)
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
			theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			scale_y_discrete(expand=expansion(add=3)) 
			ggsave(p, file=pngpath, width=25, height=85, limitsize = FALSE)
}

melt_cluster_dates = function(entry_sorted){
	entry_sorted = entry_sorted[!is.na(entry_sorted)]
	ldply(entry_sorted, function(e){
		df = melted_entries[melted_entries$entry == e,]
		df$entry_sorted = e
		df
	})
}

resort = function(entry_sorted, mydf){
	ldply(entry_sorted, function(e){
		df = mydf[mydf$entry == e,]
		if(nrow(df) == 0){
			print("no such entry at mydf: ")
			print(e)
			}else{
			df$entry_sorted = e
			df
		}
	})
}

is_lineage<-function(char_vector, lineage, tag = NULL){
  sapply(char_vector, function(string){
    if (lineage == "delta"){
      if( string == "B.1.617.2" | substr(string, 0,2) == "AY"){ return(TRUE) }else{ return(FALSE) }
    }
    if (lineage == "P.1" | lineage == "B.1.1.28.1"){
      if (string == "P.1" | string == "B.1.1.28.1") {return(TRUE)} else{return(FALSE)}
    }
	if (lineage == "BA" | lineage == "B.1.1.529" | lineage == "omicron"){
		if (string == "B.1.1.529" | substr(string,0,2) == "BA"){ return(TRUE) }else{ return(FALSE) }
	}
    if (substr(string,0,nchar(lineage)) == lineage) {return(TRUE)} else{return(FALSE)}
  })
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
nonstem_dd_pruned05 = dd_pruned05[!(dd_pruned05$tree_distance_to_closest_foreign_strain == 0 & dd_pruned05$tree_distance_to_closest_rus_strain == 0), ]

if(!is.null(opt$districts)){
	districts = read.table(opt$districts, header = F, sep = "\t")
	colnames(districts) = c("strain", "district")
}
if(!is.null(opt$mutations_file) & !is.null(opt$foreign_mutations_file)){
	print("Drawing the rug..")
	if(is.null(opt$districts)){
		rug_plot_selected_with_foreign(dd, pngpath = paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_all_delta_doublemut.png"),collapse=""), value= "corrected_grouped_mean_all",  mutations = mutations_df, foreign_mutations = f_mutations_df)
	}else{
		rug_plot_selected_with_foreign(dd, pngpath = paste(c(opt$output_folder, "/", opt$tag, "_corrected_grouped_mean_all_with_uncorr_median_foreign_date_all_delta_doublemut_districts.png"),collapse=""), value= "corrected_grouped_mean_all",  mutations = mutations_df, foreign_mutations = f_mutations_df,districts = districts )
	}
} else if (!is.null(opt$lineage)){
	print("selected lineage is ")
	print(opt$lineage)
	path = paste(c(opt$output_folder, "/", opt$tag, "_rugplot_", opt$lineage, "_pruned05.png"),collapse="")
	rug_plot_selected_lineage(dd_pruned05, path, "corrected_grouped_mean_all", lineage = opt$lineage, districts = districts)
	path = paste(c(opt$output_folder, "/", opt$tag, "_rugplot_", opt$lineage, "_all.png"),collapse="")
	rug_plot_selected_lineage(dd, path, "corrected_grouped_mean_all", lineage = opt$lineage, districts = districts)

} else {
	pngpath = paste(c(opt$output_folder, "/", opt$tag, "_entries_per_fortnight_pangolin_colored"), collapse = "")
	print("Drawing pangolin-colored barplot with entries per fortnight..")
	fortnight_plot(stem_dd, pngpath=paste(c(pngpath, "_stem.png"), collapse = ""))
	fortnight_plot(nonstem_dd, pngpath=paste(c(pngpath, "_nonstem.png"), collapse = ""))
	fortnight_plot(stem_dd, pngpath=paste(c(pngpath, "_stem_percent.png"), collapse = ""), percent = TRUE)
	fortnight_plot(nonstem_dd, pngpath=paste(c(pngpath, "_nonstem_percent.png"), collapse = ""), percent = TRUE)
	
	fortnight_plot(dd, pngpath=paste(c(pngpath, "_all.png"), collapse = ""))
	fortnight_plot(dd, pngpath=paste(c(pngpath, "_all_percent.png"), collapse = ""), percent = TRUE)

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
}

