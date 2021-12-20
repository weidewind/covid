library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
library("optparse")

option_list = list(
melted_entries_file
its = 1000
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
              help="gentag to start output file names with ", metavar="character"),
			make_option(c("-w", "--entry_weight_file"), type="character", default=NULL, 
              help="path to file produced by estimate_coverage.R", metavar="character")			  
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

options(stringsAsFactors = FALSE)

## read melted_entries.withdates.csv, cut timeline to bins ("samples"); 
## produce a list of vectors v = (X1, .., Xn), where Xi is a "species" sample count

by_week <- function(x,n=1){
	last = max(x,na.rm=T)
	m = seq(min(x,na.rm=T),last,by=paste0(n," weeks"))
	if (m[length(m)] != last)
		m = c(m, last+1)
	else
	  m[length(m)] = m[length(m)]+1
	return(m)
}


## Formulae 4a
sample_coverage <- function(df){
  N = sum(df$n)
  f1 = nrow(df[df$n == 1,])
  f2 = nrow(df[df$n == 2,])
  return(1-(f1/N)*(N-1)*f1/((N-1)*f1+2*f2))
  
}

## Formulae 4b
rarefied_coverage <- function(df, m){
  N = sum(df$n)
  t = sapply(df$n, function(x){
    (x*choose((N-x),m))/(N*choose((N-1),m))
  })
  return(1-sum(unlist(t)))
}

test_rarefied_coverage <-function(df,m){
  N = sum(df$n)
  t = sapply(df$n, function(x){
    p = x/N
    p*(1-p)^m
  })
  return(1-sum(t))
}

## Formulae 9a: how many samples is needed to reach the specified coverage?
m_asterisk <- function(df, coverage){
  N = sum(df$n)
  f1 = nrow(df[df$n == 1,])
  f2 = nrow(df[df$n == 2,])
  m_a = log(N*(1-coverage)/f1)/log((N-1)*f1/((N-1)*f1+2*f2)) - 1
  return(m_a)
}


## Pick m (the size of rarified sample) for the specified coverage
m_rarefied <- function(df, coverage, init.coverage = NULL){
  if (is.null(init.coverage))
    init.coverage = sample_coverage(df)
  m = sum(df$n)
  new_coverage = init.coverage
  while (new_coverage > coverage){
    m = m-1
    new_coverage = rarefied_coverage(df, m = m)
  }
 data.frame(m=m, rarefied_coverage = new_coverage)
}


# from a subset of melted dataframe

melted = read.csv2(melted_entries_file, header = T, stringsAsFactors=FALSE, colClasses = c("numeric", "character", "character", "character", "Date", "logical", "Date"))
melted = melted[!is.na(melted$date), ]
melted$bin = cut(melted$date, breaks = by_week(melted$date))
counts = melted %>% count(bin, entry)
print(head(counts))
bins = unique(melted$bin)
coverages = sapply(bins, function(bin){
  subset = counts[counts$bin == bin,]
  print(subset)
  sample_coverage(subset)
})
ns = sapply(bins, function(bin){
  sum(counts[counts$bin == bin, "n"])
})

processed = data.frame(bin = bins, N = ns, coverage = coverages)
print(head(processed))
ggplot(processed, aes(x = bin)) + geom_histogram(aes(y = coverage),stat="identity")

ctn = ldply(seq(0,1,0.05), function(cov){
  data.frame(coverage = cov, num_bins = nrow(processed[processed$coverage >= cov,]))
})
plot(ctn$coverage, ctn$num_bins) 

## the majority of bins reach 50% coverage? then cov_threshold = 0.5
cov_threshold = 0.5
m_a = sapply(bins, function(bin){
  subset = counts[counts$bin == bin,]
  m_asterisk(subset, cov_threshold)
})
processed$m_a = m_a
# how many samples we need to add to reach this coverage?
ggplot(processed, aes(x = bin)) + geom_histogram(aes(y = m_a),stat="identity")

# how many samples chould we take to downsample to this coverage?
m_r = ldply(bins, function(bin){
  subset = counts[counts$bin == bin,]
  m_rarefied(subset, coverage = cov_threshold, init.coverage = processed[processed$bin == bin, "coverage"])
})
processed$m = m_r$m
processed$rarefied_coverage = m_r$rarefied_coverage
head(processed)

# test: compute coverage for m (formulae (2))

t_rarefied_coverage = sapply(bins, function(bin){
  subset = counts[counts$bin == bin,]
  test_rarefied_coverage(subset, processed[processed$bin == bin, "m"])
})
processed$test_rarefied_coverage = t_rarefied_coverage






