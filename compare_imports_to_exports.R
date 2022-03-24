library("optparse")
library(plyr)
library(lubridate)
library(ggplot2)
library(ggalt)
library("RColorBrewer")
library(scales) # function percent()

option_list = list(
			make_option(c("-o", "--output_file"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("--imports_data"), type="character", default="", 
              help="", metavar="character"),
			 make_option(c("--imports_dates_stats"), type="character", default="", 
              help="", metavar="character"),
			 make_option(c("--exports_dates_stats"), type="character", default="", 
              help="", metavar="character"),
			 make_option(c("--entries_file"), type="numeric", default="",
			 help="transmission_lineages.withduplicates.out.broken"),
			 make_option(c("--reverse"), type="numeric", default=0,
			 help="1 or 0"),
			 make_option(c("--lineages"), type="character", default=NULL, help="", metavar="character")
			
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)
print("options: ")
print(opt)

imports = read.csv2(file = opt$imports_dates_stats,sep = "\t", header = T, stringsAsFactors=FALSE)
exports = read.csv2(file = opt$exports_dates_stats,sep = "\t", header = T, stringsAsFactors=FALSE)
if(opt$reverse == 1){
	t = data.frame(entry = exports$entry, 
	closest_ingroup_nodes = exports$closest_outgroup_nodes,
	tree_distance_to_closest_rus_strain = exports$tree_distance_to_closest_foreign_strain,
	closest_outgroup_nodes = exports$closest_ingroup_nodes,
	tree_distance_to_closest_foreign_strain = exports$tree_distance_to_closest_rus_strain)
	exports = t
}
print("Imports:")
print(str(imports))
print(head(imports))
print("Exports:")
print(str(exports))
print(head(exports))

is_secondary = read.csv2(file = opt$entries_file,sep = "\t", header = F, stringsAsFactors=FALSE)
names(is_secondary) = c("entry", "cluster_strains", "previous_export")

imports = merge(imports, is_secondary[,c("entry", "previous_export")], by = "entry", all.x = TRUE)
exports$previous_export = rep(NA,nrow(exports))
imports$type = sapply(imports$previous_export, function(ex){
	if (ex == "" || is.na(ex)){return("primary import")}else{return("secondary import")}
})
exports$type = rep("export", nrow(exports))
dat = rbind(imports, exports)
print("Dat:")
print(str(dat))
print(head(dat))
dat$tree_distance_to_closest_rus_strain = as.double(dat$tree_distance_to_closest_rus_strain)
dat$tree_distance_to_closest_foreign_strain = as.double(dat$tree_distance_to_closest_foreign_strain)
dat$distance_difference = dat$tree_distance_to_closest_foreign_strain - dat$tree_distance_to_closest_rus_strain
print(str(dat))
print(head(dat))
print(tail(dat))
print("number of exports:")
print(nrow(dat[dat$type == "export",]))
print("number of primary imports:")
print(nrow(dat[dat$type == "primary import",]))
print("number of secondary imports:")
print(nrow(dat[dat$type == "secondary import",]))

is_lineage<-function(char_vector, lineage, tag = NULL){
  sapply(char_vector, function(string){
	if (is.na(string)) {return(FALSE)}else{
		if (lineage == string){return(TRUE)}else{
			if (lineage == "delta"){
			  if( string == "B.1.617.2" | substr(string, 0,2) == "AY"){ return(TRUE) }else{ return(FALSE) }
			}
			if (lineage == "P.1" | lineage == "B.1.1.28.1"){
			  if (string == "P.1" | string == "B.1.1.28.1") {return(TRUE)} else{return(FALSE)}
			}
			if (lineage == "BA" | lineage == "B.1.1.529" | lineage == "omicron"){
				if (string == "B.1.1.529" | substr(string,0,2) == "BA"){ return(TRUE) }else{ return(FALSE) }
			}
			if (substr(string,0,nchar(lineage)+1) == paste(c(lineage, "."),collapse = "")) {return(TRUE)} else{return(FALSE)}
		}
	}
  })
}

plot_dist = function(dat, pngpath){
	p <- ggplot(dat, aes(x=tree_distance_to_closest_foreign_strain, fill = type)) +
		geom_histogram(alpha = 0.5, position="identity",color = "black") +
		scale_fill_manual(values=c("deepskyblue1","deepskyblue4","coral"), 
		breaks = c("primary import","secondary import","export"), 
		labels =c("primary import","secondary import","export"), name="") +
		xlab("distance to closest foreign strain")
		theme_minimal()

	ggsave(p, file= pngpath, width=7, height=7, limitsize = FALSE)
}

plot_delta_dist = function(dat, pngpath){	
	p <- ggplot(dat, aes(x=distance_difference, fill = type)) +
		geom_histogram(alpha = 0.5, position="identity",color = "black") +
		scale_fill_manual(values=c("deepskyblue1","deepskyblue4","coral"), 
		breaks = c("primary import","secondary import","export"), 
		labels =c("primary import","secondary import","export"), name="") +
		xlab("distance to foreign - distance to russian")
		theme_minimal()

	ggsave(p, file= pngpath, width=7, height=7, limitsize = FALSE)
}


pngpath = paste(c(opt$output_file, ".png"), collapse = "")
plot_dist(dat,pngpath)
pngpath = paste(c(opt$output_file, "_diff.png"), collapse = "")
plot_delta_dist(dat,pngpath)
if (!is.null(opt$lineages)){
	impdat = read.csv2(file = opt$imports_data,sep = "\t", header = T, stringsAsFactors=FALSE)
	dat = merge(dat, impdat[, c("entry", "major_lineage")], by = "entry")
	lins = unlist(strsplit(opt$lineages, ","))
	sapply(lins, function(lin){
		pngpath = paste(c(opt$output_file, "_", lin, ".png"), collapse = "")
		plot_dist(dat[is_lineage(dat$major_lineage,lin), ], pngpath)
		pngpath = paste(c(opt$output_file, "_", lin, "_diff.png"), collapse = "")
		plot_delta_dist(dat[is_lineage(dat$major_lineage,lin), ], pngpath)
	})
}