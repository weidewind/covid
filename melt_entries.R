library("optparse")
library(plyr)
library(lubridate)


option_list = list(
			make_option(c("-p", "--pangolined_file"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("-d", "--dates_stats_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-y", "--entry_dates_file"), type="character", default="", 
              help="", metavar="character"),			
			make_option(c("-e", "--entry_strains_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-l", "--leaf_dates_file"), type="character", default="", 
              help="", metavar="character"),
			 make_option(c("-o", "--output"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("-t", "--tag"), type="character", default="", 
              help="gentag to start output file names with ", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

print("options:")
print(opt)

pangolined = read.csv2(opt$pangolined_file, header = T, sep = ",")
cluster_strains = read.csv2(file = opt$entry_strains_file, sep = "\t", header = F, stringsAsFactors=FALSE)
cluster_strains = data.frame(entry=as.character(cluster_strains$V1), strains=as.character(cluster_strains$V2))
dates_stats = read.csv2(file = opt$dates_stats_file,sep = "\t", header = T, stringsAsFactors=FALSE)
entry_dates = read.csv2(file = opt$entry_dates_file,sep = ";", header = T, stringsAsFactors=FALSE)
entry_to_date = data.frame(entry = entry_dates$entry, entry_date = as.Date(entry_dates$corrected_grouped_mean_all, format = "%Y-%m-%d"))
stem_strains = dates_stats[as.numeric(dates_stats$tree_distance_to_closest_rus_strain) == 0 & as.numeric(dates_stats$tree_distance_to_closest_foreign_strain) == 0,"close_rus_strains"]
print("stem_strains:")
print(str(stem_strains))
#print(stem_strains)

stem_strains.list <-unlist(lapply(stem_strains, function(ss){
	strains = unlist(strsplit(ss, ","))
}))
print("stem_strains.list:")
print(str(stem_strains.list))
#print(stem_strains.list)
stem_strains.unique = unique(stem_strains.list)


print("cluster_strains:")
head(cluster_strains)

strain_to_date = read.csv2(file = opt$leaf_dates_file,sep = "\t", header = F, stringsAsFactors=FALSE)
strain_to_date = data.frame(strain = strain_to_date$V1, date = as.Date(strain_to_date$V2, format = "%Y-%m-%d"))
is_stem = sapply(strain_to_date$strain, function(s){
	if(s %in% stem_strains.unique){TRUE}else{FALSE}
})
strain_to_date$is_stem = is_stem

print("strain_to_date:")
head(strain_to_date)

cluster_strains.list <- split(cluster_strains, seq(nrow(cluster_strains)))
strain_to_entry = ldply(cluster_strains.list, function(row){
	print("str row:")
	print(str(row))
	print("row:")
	print(row)
	print("row['strains']")
	print(row["strains"])
	print(str(row["strains"]))
	strains = unlist(strsplit(row[["strains"]], ";"))
	print("strains:")
	print(str(strains))
	lineages  = unlist(sapply(strains, function(st){
		l = pangolined[pangolined$taxon == st, "lineage"]
		if (length(l) == 0){
			print("No lineage for taxon ")
			print(st)
		}
		l
	}))
	print("lineages:")
	print(str(lineages))
	subdf = data.frame(entry = rep(row[["entry"]], length(strains)), strain = strains, lineage = lineages)
	print("subdf head:")
	print(head(subdf))
	subdf
})

print("strain_to_entry")
print(head(strain_to_entry))

df = join(strain_to_entry, strain_to_date, by = "strain")
df = join(df, entry_to_date, by = "entry")
write.csv2(df, opt$output, quote = FALSE, row.names = FALSE)