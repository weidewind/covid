#!/usr/bin/env Rscript

# DELTRAN maximum parsimony reconstruction

mylib = file.path("~", "R","library", fsep = .Platform$file.sep)

#req = c("igraph", "phangorn", "ape", "optparse")
##install.packages("pak", lib=mylib)
##pkg_install(pkg, lib = mylib, upgrade = TRUE, ask = F)
# miss <- setdiff(req,  installed.packages(lib.loc = mylib)[,"Package"]) #installed.packages(lib.loc = mylib)[,"Package"]
# if (length(miss)) {
	# print ("Installing packages: ")
	# print (miss)
	# install.packages(req, lib=mylib) #, lib.loc=mylib
# }
library(phangorn, lib.loc = mylib)
library(ape, lib.loc = mylib)
library(optparse, lib.loc = mylib)


 library(phangorn)
 library(ape)
 library(optparse)
 library(plyr)
 library(castor)
 library(phylotools)


option_list = list(
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="path to tree file", metavar="character"),
  make_option(c("-s", "--states"), type="character", default=NULL, 
              help="path to states file", metavar="character"),	
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="path to output file with ancestral states", metavar="character")			  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# with (opt,{
	# print ("Reading data..")
	# s <- read.phyDat(states, format = "fasta") #phydat
	# t <- read.tree(tree)
	# print ("MP ancestral reconstruction..")
	# anc.mpr <- ancestral.pars(tree, s, "MPR") #tree - pml
	# write.phyDat(anc.mpr, file, format = "fasta")
# })


myTree <- ape::read.tree(text='((((A:1,B:1,C:1):1,((D:1,E:1,(L:1,M:1):1):1,(F:1,G:1):1):1,(H:1,I:1,J:1):1):1,K:1):1,N:1);')
fafile <- "C:/Users/ִלטענטי/Documents/workspace/covid/data/mptest.fa"

bwTree <-  ape::read.tree(text='(((A:1,B:1,C:1,D:1):1,(E:1,(O:1,P:1,(Q:1,R:1):1):1,G:1,H:1):1):1,(I:1,(J:1,(S:1,T:1,(U:1,V:1, (W:1,X:1):1):1):1,(Y:1,Z:1):1,M:1):1));')
bwfafile <- "C:/Users/ִלטענטי/Documents/workspace/covid/data/mptest_bw.fa"   


#----------- playing with FANGORN


myTree <- bwTree
fafile <- bwfafile
lev = c(1,2)
myTree <-multi2di(myTree)

myDat <- read.phyDat(fafile, format = "fasta", type="USER", levels=lev)

anc.mpr <- ancestral.pars(myTree, myDat, "ACCTRAN", return="prob") #return="phydat" gives nice output with as.data.frame(anc.mpr), but does not show ambiguities
stat.pars <- parsimony(multi2di(myTree), myDat) # some edge statistics?
plotAnc(myTree, anc.mpr, 1, col=c("red", "yellow", "blue"))
write.phyDat(anc.mpr, "C:/Users/ִלטענטי/Documents/workspace/covid/data/mptest_out.nexus",  format = "nexus", colsep = "", nbcol = -1)

stateProbs<-ldply(c(1:length(anc.mpr)), function(nnode){
  row =c(names(anc.mpr)[nnode],unlist(lapply(c(1:length(lev)), function(nstate){anc.mpr[[nnode]][nstate]})))
})
stateProbs

myTree$node.label <- names(anc.mpr)[length(myTree$tip.label)+1:length(names(anc.mpr))]
write.csv(stateProbs,file="C:/Users/ִלטענטי/Documents/workspace/covid/data/mptest.out.probs", quote = F, row.names = F)
write.tree(myTree, file="C:/Users/ִלטענטי/Documents/workspace/covid/data/mptest.out.newick")

#-----------CASTOR

myFa <- read.fasta(fafile)
nstates <- 2
levels <- seq(nstates)
#myTree <-multi2di(myTree)

tip_states <- sapply(myTree$tip.label, function(label){
  as.integer(myFa[myFa$seq.name == label,2])
})
castorRes<-asr_max_parsimony(myTree, tip_states, Nstates=nstates,
                  transition_costs="all_equal",
                  edge_exponent=0, weight_by_scenarios=TRUE,
                  check_input=TRUE)
castorProbs<-data.frame(castorRes$ancestral_likelihoods)
myTree$node.label <- seq(length(myTree$tip.label)+1,length(myTree$tip.label)+nrow(castorProbs))
names(castorProbs)<-seq(nstates)
row.names(castorProbs) <-myTree$node.label
write.tree(myTree, file="C:/Users/ִלטענטי/Documents/workspace/covid/data/castor.out.newick")

# merge state data for leaves and internal nodes
tipdat<-ldply(seq(1, length(tip_states)), function(i){
  sapply(seq(nstates), function(state){
    if (tip_states[i] == state){1.0}else{0.0}
  })
})
names(tipdat)<-seq(nstates)
row.names(tipdat)<-myTree$tip.label
allProbs<-rbind(tipdat, castorProbs)
write.table(allProbs,file="C:/Users/ִלטענטי/Documents/workspace/covid/data/castor.out.probs", quote = F, row.names = T, col.names = F)

#

##

# ----- Convert castor output to phy.dat-like object for plot printing
tiplist<-lapply(seq(1, length(tip_states)), function(i){
  matrix(sapply(seq(1,nstates), function(state){
    if (tip_states[i] == state){1.0}else{0.0}
  }), nrow=1)
})
names(tiplist)<-names(tip_states)
nodelist<-split(castorProbs, seq(nrow(castorProbs)))
nodelist_of_matrices<-lapply(nodelist,function(r){matrix(unlist(r), nrow=1)})
names(nodelist_of_matrices)<-row.names(castorProbs)
castorProbsPhy<-c(tiplist, nodelist_of_matrices)
att<-attributes(anc.mpr)
att$names<-attributes(castorProbsPhy)$names
attributes(castorProbsPhy)<-att
plotAnc(myTree, castorProbsPhy, 1, col=c("red", "yellow", "blue"))
##



#----------------

#reorder(myTree, order = "postorder")
parent = myTree$edge[,1]
child = myTree$edge[,2]
res <- numeric(max(myTree$edge))
res[1:Ntip(myTree)] <- 1L
for(i in postorder(myTree)){ #rev(1:nrow(myTree$edge))
 # if (allProbs[child[i],1] == 0){print (child[i])}
  #print(paste(c(parent[i], child[i], allProbs[parent[i],1],allProbs[child[i],1])))
  res[parent[i]] = res[parent[i]] + res[child[i]]  
}


