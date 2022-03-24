library(ggplot2)
library(plyr)
library("optparse")


option_list = list(
			make_option(c("--folder"), type="character", default="", 
              help="", metavar="character"),
			make_option(c("--output"), type="character", default="", 
              help="", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)
folder=opt$folder
file=paste(c(folder,"country_stat.include_inner.1.country_stats"), collapse = "")
exfile=paste(c(folder,"country_stat.include_inner.exports.1.country_stats"), collapse = "")
infile=paste(c(folder,"country_stat.include_inner.inner.1.country_stats"), collapse = "")
pfile=paste(c(folder,"country_stat.include_inner.parents.1.1.country_stats"), collapse = "")
e_to_ex_file=paste(c(folder,"transmission_lineages.withduplicates.out.broken.entries.exports"), collapse = "")
distfile=paste(c(folder,"translin.inner.withdist"), collapse = "")
pdistfile=paste(c(folder,"translin.parents.1.withdist"), collapse = "")

e_to_ex=read.csv(e_to_ex_file, sep = "\t", header = FALSE)
dists = read.csv(distfile, sep = ",", header = TRUE)
pdists = read.csv(pdistfile, sep = ",", header = TRUE)
pdists$node = as.character(pdists$node)
colnames(e_to_ex) = c("entry", "export")
cstat=read.csv(file, sep = "\t", header = TRUE)
cstat$type = rep("import", nrow(cstat))
ecstat=read.csv(exfile, sep = "\t", header = TRUE)
ecstat$type = rep("export", nrow(ecstat))
icstat=read.csv(infile, sep = "\t", header = TRUE)
icstat$type = rep("internal", nrow(icstat))
pcstat=read.csv(pfile, sep = "\t", header = TRUE)
pcstat$type = rep("parent", nrow(pcstat))
allcstat = rbind(cstat,ecstat,icstat,pcstat)
allcstat$prop = allcstat$rus_count/(allcstat$rus_count+allcstat$foreign_count)
allcstat$dist = sapply(allcstat$entry, function(e){
  if(allcstat[allcstat$entry == e, "type"] == "import"){
     0
  }else{
  d = dists[dists$node == e, "distance"]
  if(length(d) == 0){
    d = mean(pdists[pdists$node == e, "distance"]) # mean!
    if(length(d) == 0){
      0 #import
    }else{d}
  }else{
    d
  }
  }
})

# best = c("156357", "157318", "157647")
# prop_for_best = ldply(best, function(e){
  # data.frame(entry = e, prop = allcstat[allcstat$entry == e, "prop"])
# })

# prop_for_best_exports = ldply(e_to_ex[e_to_ex$entry %in% best, "export"], function(e){
  # data.frame(entry = e, prop = allcstat[allcstat$entry == e, "prop"])
# })

# ggplot(allcstat) +geom_histogram(aes(x=prop,fill = type)) +
  # scale_fill_manual(values=c("deepskyblue3","deepskyblue4","coral","gray"), breaks = c("parent","import","export","internal"), labels =c("parent","import", "export","internal"), name="")   +
  # geom_vline(data = prop_for_best, aes(xintercept = prop),color = "deepskyblue4") +
  # geom_vline(data = prop_for_best_exports, aes(xintercept = prop), color = "coral")+
  # xlab("rus/(rus+foreign)")+ theme_minimal()

p <- ggplot(allcstat) + geom_count(aes(x = dist, y = prop, color = type)) +
  scale_color_manual(values=c("deepskyblue3","deepskyblue4","coral","gray"), breaks = c("parent","import","export","internal"), labels =c("parent","import", "export","internal"), name="") +
  ylab("rus/(rus+foreign)")+ xlab("distance from entry")+ theme_minimal()
ggsave(p, file=paste(c(opt$output, "_rus_prop", ".png"), collapse = ""), width=10, height=10, limitsize = FALSE)

p <- ggplot(allcstat[allcstat$type != "internal",]) + geom_count(aes(x = dist, y = prop, color = type)) +
  scale_color_manual(values=c("deepskyblue3","deepskyblue4","coral"), breaks = c("parent","import","export"), labels =c("parent","import", "export"), name="") +
  ylab("rus/(rus+foreign)")+ xlab("distance from entry")+ theme_minimal()
  ggsave(p, file=paste(c(opt$output, "_rus_prop_nointernal", ".png"), collapse = ""), width=10, height=10, limitsize = FALSE)

