library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
library(optparse)

option_list = list(
			make_option(c("-m", "--melted_entries_file"), type="character", default=NULL, 
              help="path to file", metavar="character"),
			make_option(c("-i", "--iterations"), type="integer", default=1000, 
              help="number of iterations", metavar="character"),
			make_option(c("-c", "--coverage_threshold"), type="numeric", default="", 
              help="coverage to rarefy to", metavar="character"),
			make_option(c("-o", "--output"), type="character", default="", 
              help="path to output file", metavar="character")	  
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

sample_coverage <- function(df){
  N = sum(df$n)
  f1 = nrow(df[df$n == 1,])
  f2 = nrow(df[df$n == 2,])
  res = 1-(f1/N)*(N-1)*f1/((N-1)*f1+2*f2)
  # if (is.nan(res)){
    # print("Debug sample_coverage: df ")
	# print(df)
  # }
  return(res)
  
}

rarefied_coverage <- function(df, m){
  N = sum(df$n)
  t = sapply(df$n, function(x){
    # equivalent for (x*choose((N-x),m))/(N*choose((N-1),m)):
	res =  x/N
    t = m
    while (t >= 1){
      res = res*(N-x-t+1)/(N-t)
      t = t-1
    }
	# if (is.nan(res)){
		# ve = c(x, N, m, choose((N-x),m),choose((N-1),m))
		# names(ve) = c("x", "N", "m", "choose((N-x),m)","choose((N-1),m)")
		# print(ve)
	# }
	res
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

m_asterisk <- function(df, coverage){
  N = sum(df$n)
  f1 = nrow(df[df$n == 1,])
  f2 = nrow(df[df$n == 2,])
  m_a = log(N*(1-coverage)/f1)/log((N-1)*f1/((N-1)*f1+2*f2)) - 1
  return(m_a)
}


m_rarefied <- function(df, coverage, init.coverage = NULL){
  # print("Debug: init.coverage")
  # print(init.coverage)
  # print("Debug: sample coverage")
  # print(sample_coverage(df))
  if (is.null(init.coverage))
    init.coverage = sample_coverage(df)
  m = sum(df$n)
  new_coverage = init.coverage
  while (new_coverage > coverage){
    m = m-1
    new_coverage = rarefied_coverage(df, m = m)
	# print("Debug: m, new_coverage")
	# print(m)
	# print(new_coverage)
	# if(is.nan(new_coverage)){
		# print("Debug: df")
		# print(df)
	# }
  }
 data.frame(m=m, rarefied_coverage = new_coverage)
}


melted = read.csv2(opt$melted_entries_file, header = T, stringsAsFactors=FALSE, colClasses = c("numeric", "character", "character", "character", "Date", "logical", "Date"))
melted = melted[!is.na(melted$date), ]
melted$bin = cut(melted$date, breaks = by_week(melted$date, 2))
counts = melted %>% count(bin, entry)
# print("Debug: head(counts))")
# print(head(counts))
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
gp = ggplot(processed, aes(x = bin)) + geom_histogram(aes(y = coverage),stat="identity") + geom_hline(yintercept=opt$coverage_threshold, color="red")
ggsave(gp, file=paste(c(opt$output, ".png"), collapse=""), width=10, height=5)	

m_r = ldply(bins, function(bin){
  subset = counts[counts$bin == bin,]
  m_rarefied(subset, coverage = opt$coverage_threshold, init.coverage = processed[processed$bin == bin, "coverage"])
})
processed$m = m_r$m
processed$rarefied_coverage = m_r$rarefied_coverage
head(processed)

tot_its = c()
for (i in 1:opt$iterations){
  iteration = unique(unlist(sapply(bins, function(bin){
    subset = melted[melted$bin == bin, ]
    m = processed[processed$bin == bin, "m"]
    subset[sample(nrow(subset),m), "entry"]
  })))
  tot_its = c(tot_its, iteration)
}
entry_freq = as.data.frame(table(tot_its))
entry_freq$Freq = entry_freq$Freq/opt$iterations
write.csv2(entry_freq,file = opt$output, row.names = FALSE, quote = FALSE)
