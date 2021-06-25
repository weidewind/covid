#!/usr/bin/env Rscript

library("optparse")
library("forcats")
library("ggplot2")
 
option_list = list(
  make_option(c("-i", "--important_lineages"), type="character", default=NULL, 
              help="path to file", metavar="character"),
  make_option(c("-p", "--pangolined_file"), type="character", default=NULL, 
              help="path to file", metavar="character"),
  make_option(c("-o", "--output_folder"), type="character", default="", 
              help="", metavar="character"),
make_option(c("-e", "--entry_file"), type="character", default="", 
              help="", metavar="character"),
  make_option(c("-t", "--tag"), type="character", default="", 
              help="gentag to start output file names with ", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(error = quote({dump.frames(to.file=TRUE); q()}))

pangolined = read.csv2(opt$pangolined_file, header = T, sep = ",", stringsAsFactors=FALSE)
important_lin = read.csv2(opt$important_lineages, header = T, sep = ",", stringsAsFactors=FALSE) 


## For each unique pangolin lineage, count strains attributed to it. 
## Write counts to file, ordered by lineage importance and then by count
all_lin = unique(pangolined$lineage)
counts = sapply(all_lin, function(lin){nrow(pangolined[pangolined$lineage == lin,])})
print(str(important_lin))
print(head(important_lin))
print(str(all_lin))
print(head(all_lin))
desc = unlist(sapply(all_lin, function(lin){
  d = important_lin[important_lin$lineage == lin, "desc"][1]
  if (is.na(d)){d = ""}
  d
  }))
counts_df = data.frame(lineage = names(counts), count = counts, important = names(counts) %in% important_lin$lineage, desc = desc ) #names(counts) %in% important_lin$lineage
ordered_counts_df = counts_df[
  order( counts_df[,"important"], counts_df[,"count"], decreasing = T ),
]
write.csv2(ordered_counts_df, paste(c(opt$output_folder,"/", opt$tag, "_pangolin_counts.csv"),collapse = ""), quote = F, row.names = F)


## Present the same data as a bar plot
print("Plotting counts..")
df = data.frame(pangolined, important = pangolined$lineage %in% important_lin$lineage)
p = ggplot(df, aes(lineage, fill = important)) + scale_fill_manual(values=c("darkblue", "darkred"))+ geom_bar(aes(x=forcats::fct_infreq(lineage)))  + theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_text(data = df[df$important == TRUE,], stat = "count", aes(label = after_stat(count)), vjust = -1, hjust = 0)+ ggtitle(opt$tag)
output = paste(c(opt$output_folder,"/", opt$tag, "_barplot.png"),collapse = "")
ggsave(p, file=output, width=15, height=10)


## For each entry, analyze its pangolin content and write it to file. 
## If >50% of entry content is attributed to lineage x, x is called the major_lineage for this entry
entries = read.csv2(opt$entry_file, sep = "\t", header = F, stringsAsFactors=FALSE)

colnames(entries) = c("entry", "strains")

print("Collecting lineage counts..")
lin_stat = sapply(entries$strains, function(s){
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

print("Searching for major lineages..")
major_lins = sapply(entries$strains, function(s){
    lins = sapply(unlist(strsplit(s, ";")), function(st){
       pangolined[pangolined$taxon == st, "lineage"]
      })
    tbl = as.data.frame(table(lins))
    if (max(tbl$Freq)/sum(tbl$Freq) > 0.5){major_lin = as.character(tbl[which.max(tbl$Freq), "lins"])}else{major_lin = NA}
    major_lin
})
names(major_lins) = NULL

entries = data.frame(entry = entries$entry, lineage_stat = lin_stat, major_lineage = major_lins, important = major_lins %in% important_lin$lineage)
row.names(entries) = NULL
write.table(entries, paste(c(opt$output_folder,"/", opt$tag, "_pangolin_counts_for_entry.csv"),collapse = ""),sep=",", row.names = F, quote = F)


## Present the same data as anordered barplot
print("Preparing to plot..")
mls = unique(entries$major_lineage)
ml_counts = sapply(mls, function(ml){
  if(is.na(ml)){nrow(entries[is.na(entries$major_lineage),])}else{
  nrow(entries[!is.na(entries$major_lineage) & entries$major_lineage == ml,])
  }
})
names(ml_counts) = NULL
imp =  mls %in% important_lin$lineage

df = data.frame(lineage = mls, count = ml_counts, important = imp, proportion = ml_counts/sum(ml_counts))
df = df[order(df$count),]

p = ggplot(entries, aes(major_lineage, fill = important)) + scale_fill_manual(values=c("darkblue", "darkred"))+ geom_bar(aes(x=forcats::fct_infreq(major_lineage)))  + theme_minimal() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_text(data = entries[entries$important == TRUE,], stat = "count", aes(label = after_stat(count)), vjust = -1, hjust = 0)+ ggtitle(opt$tag) + xlab("major entry lineage")
output = paste(c(opt$output_folder,"/", opt$tag, "_entry_barplot.png"),collapse = "")
ggsave(p, file=output, width=15, height=10)

p = ggplot(df, aes(x=lineage, y = proportion, fill = important)) + scale_fill_manual(values=c("darkblue", "darkred"))+ geom_bar(aes(x=reorder(lineage, -count)),stat = "identity")  + theme_minimal()  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(opt$tag) + xlab("major entry lineage")+ ylab("proportion") + geom_text(data = df[df$important == TRUE,], stat = "identity", aes(label = round(proportion,2)), vjust = -0.2, hjust = 0)
output = paste(c(opt$output_folder,"/", opt$tag, "_entry_barplot_percent.png"),collapse = "")
ggsave(p, file=output, width=15, height=10)

