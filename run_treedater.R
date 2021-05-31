library(ape)
library(treedater)
library(optparse)
library(lubridate)
library(plyr)
 
option_list = list(
	make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="tree file", metavar="character"),
	make_option(c("-d", "--dates"), type="character", default=NULL, 
              help="file with dates", metavar="character"),
	make_option(c("-u", "--upper"), type="character", default=NULL, 
	            help="file with max dates", metavar="character"),
	make_option(c("-l", "--lower"), type="character", default=NULL, 
	            help="file with min dates", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
	make_option(c("-s", "--seqlength"), type="integer", default="0", 
	            help="length of the alignment used to produce the tree", metavar="character") # 29903 for far20
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tre <- ape::read.tree( opt$tree )
dates <- read.csv2(opt$dates, header=T, sep = ",")
mindates <- read.csv2(opt$lower, header=T, sep = ",")
maxdates <- read.csv2(opt$upper, header=T, sep = ",")

#tre <- ape::read.tree("C:/Users/ִלטענטי/Documents/CMD/covid/gennady_tree/for_iqtree_far20.fasta.treefile")
#dates <- read.csv2("C:/Users/ִלטענטי/Documents/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.mean", header=T, sep = ",")
#mindates <- read.csv2("C:/Users/ִלטענטי/Documents/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.min", header=T, sep = ",")
#maxdates <- read.csv2("C:/Users/ִלטענטי/Documents/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.max", header=T, sep = ",")

# one leaf in our tree can correspond to several strains. here we assign to every leaf a pre-computed mean date
tipdates <- as.Date( dates$date, format="%Y-%m-%d")
sts <- lubridate::decimal_date( tipdates )
names(sts) <- dates$name
sts<-sts[!is.na(sts)]

# treedater can also use time intervals, so here we also assign corresponding pre-computed time intervals [min date, max date] to leafs
# if a mean date is also present, treedater will use it as the initial value for further optimization
datebounds <- data.frame(lower = lubridate::decimal_date(as.Date(mindates$date, format="%Y-%m-%d")), upper = lubridate::decimal_date(as.Date(maxdates$date, format="%Y-%m-%d")))
row.names(datebounds) <- dates$name

# if no date is given for a leaf, we assign to it the time period [global min date, global max date]
globmin <- min(datebounds$lower[!is.na(datebounds$lower)])
globmax <- max(datebounds$upper[!is.na(datebounds$upper)])
ext_lower <- sapply(datebounds$lower, function(d){
  if(is.na(d)){globmin}else{d}
})
ext_upper <- sapply(datebounds$upper, function(d){
  if(is.na(d)){globmax}else{d}
})
ext_datebounds <-data.frame(lower = ext_lower, upper = ext_upper)
row.names(ext_datebounds) <- row.names(datebounds)


# `tre` is an `ape::phylo` phylogeny, 
# `sts` is a named vector of sample times for each tip in `tre`
# `s` is the length of the genetic sequences used to estimate `tre`
#' @param estimateSampleTimes If some sample times are not known with certainty,
#'         bounds can be provided with this option. This should take the
#'         form of a data frame with columns 'lower' and 'upper'
#'         providing the sample time bounds for each uncertain tip. Row
#'         names of the data frame should correspond to elements in
#'         tip.label of the input tree. Tips with sample time bounds in
#'         this data frame do not need to appear in the *sts* argument,
#'         however if they are included in *sts*, that value will be
#'         used as a starting condition for optimisation.
dtr <- dater( tre, sts, s=opt$seqlength, estimateSampleTimes=ext_datebounds, clock = 'strict', ncpu = 20)
summary(dtr)
write.tree( dtr, file=opt$output)
outliers <- outlierTips( dtr, alpha = 0.20)
summary(outliers)

