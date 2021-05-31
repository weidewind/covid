#!/usr/bin/env Rscript

# DELTRAN maximum parsimony reconstruction

mylib = file.path("~", "R","library", fsep = .Platform$file.sep)

#req = c( "ape", "optparse")
##install.packages("pak", lib=mylib)
##pkg_install(pkg, lib = mylib, upgrade = TRUE, ask = F)
# miss <- setdiff(req,  installed.packages(lib.loc = mylib)[,"Package"]) #installed.packages(lib.loc = mylib)[,"Package"]
# if (length(miss)) {
# print ("Installing packages: ")
# print (miss)
# install.packages(req, lib=mylib) #, lib.loc=mylib
# }

library(ape, lib.loc = mylib)
library(optparse, lib.loc = mylib)
library(plyr, lib.loc = mylib)
library(castor, lib.loc = mylib)
library(phylotools, lib.loc = mylib)


option_list = list(
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="path to tree file", metavar="character"),
  make_option(c("-s", "--states"), type="character", default=NULL, 
              help="path to states file", metavar="character"),	
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="path to output file with ancestral states", metavar="character")			  
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

with (opt,{

  print("Parsing tree..\n")
  myTree <- ape::read.tree(file=tree)
  #myTree <-multi2di(myTree)
  
  print("Reading states..\n")
  myStates <- read.csv(states, header=FALSE, sep="\t",colClasses=c('character','integer','character'))
  nstates <- 2
  levels <- seq(nstates)
  tip_states <- sapply(myTree$tip.label, function(label){
    as.integer(myStates[myStates[,1] == label,2])
  })
  
  print("Starting parsimony reconstruction..")
  castorRes<-asr_max_parsimony(myTree, tip_states, Nstates=nstates,
                               transition_costs="all_equal",
                               edge_exponent=0, weight_by_scenarios=TRUE,
                               check_input=TRUE)
  print("Parsimony reconstruction finished")
  
  castorProbs<-data.frame(castorRes$ancestral_likelihoods)
  myTree$node.label <- seq(length(myTree$tip.label)+1,length(myTree$tip.label)+nrow(castorProbs))
  names(castorProbs)<-seq(nstates)
  row.names(castorProbs) <-myTree$node.label
  
  tipdat<-ldply(seq(1, length(tip_states)), function(i){
	    sapply(seq(nstates), function(state){
			if (tip_states[i] == state){1.0}else{0.0}
	    })
	})
  names(tipdat)<-seq(nstates)
  row.names(tipdat)<-myTree$tip.label
  allProbs<-rbind(tipdat, castorProbs)
  
  write.table(allProbs,file=paste(c(output, ".probs"),collapse=""), quote = F, row.names = T, col.names = F)
  write.tree(myTree, file=paste(c(output, ".newick"),collapse=""))

  })