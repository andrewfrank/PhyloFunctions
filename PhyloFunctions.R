# ~~~~ Andrew Frank (c) 2016

library(ape)
library(spider)
library(phangorn)
library(phytools)
library(ggtree)

library(rwty)

library(ggplot2)
library(parallel)
library(plyr)
library(reshape2)
library(gridExtra)


#### Functions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### nexus2DNAbin

nexus2DNAbin <- function(x) {
	tmp <- do.call(rbind,x)
	tmp.DNAbin <- as.DNAbin(tmp)
	return(tmp.DNAbin)
}

## Description

# A function that takes an object given by read.nexus.data and returns a dataframe similar
# to the output of read.DNA; useful for getting around ape's weird way of reading in
# nexus files.

## Arguments

# x - an object of class list containing DNA sequence data, like those given by
# 	  read.nexus.data

## Source

# Andrew Frank (c) 2016

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### cbind.fill

cbind.fill <- function(x) {
    nm <- lapply(x, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x)
    	rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

## Description

# Concatenate data frames by column

## Arguments

# x - an object of class list containing matrices or data frames

## Source

# http://stackoverflow.com/questions/7962267/cbind-a-df-with-an-empty-df-cbind-fill

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### concat.alignmentList

concat.alignmentList <- function(x) {
	tmp <- lapply(x,as.character)
	tmp.align <- cbind.fill(tmp)
	tmp.DNAbin <- as.DNAbin(tmp.align)
	return(tmp.DNAbin)
}

## Description

# A function that takes a list of alignments (of class DNAbin) and concatenates them
# together into a single DNAbin object. If a taxa is missing data for that character, the
# function inserts NA.

# This function needs to be rewritten so that NAs are not inserted *literally* into the
# alignment, which causes alignment issues (because inserting NA takes up two sites).

## Arguments

# x - an object of class list containing DNA sequence alignments in DNAbin format

## Source

# Andrew Frank (c) 2016

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

random.allele.alignment <- function(x) {
	taxa.num <- nrow(x)
	taxa.IDs <- seq(from = 1, to = taxa.num, by = 2)
	selected.alleles <- sapply(
		taxa.IDs,
		function(y){
			random.num <- sample(1:2,1)
			if (random.num == 1) {
				return(y)
			} else if (random.num == 2) {
				return(y+1)
			}
		}
	)
	random.alignment <- x[selected.alleles,]
	return(random.alignment)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### phylo.sites

phylo.sites <- function(x,phased) {

	variable.sites <- seg.sites(x)
	variable.alignment <- as.character(x[,variable.sites])
	colnames(variable.alignment) <- variable.sites
	taxa.num <- nrow(variable.alignment)

	phylo.sites.check <- apply(
		variable.alignment,
		2,
		function(y) {
			if (
				taxa.num %in% table(y) ||
				(taxa.num - 1) %in% table(y)
			) {
				return(FALSE)
			} else if (
				phased == TRUE &&
				(taxa.num - 2) %in% table(y)
			) {
				repeated.base <- names(which(table(y) == 2))
				poss.ind <- names(which(y == repeated.base))
				equality.check <- all.equal(poss.ind[1],poss.ind[2])
				return(equality.check != "1 string mismatch")
			} else return(TRUE)
		}
	)

	phylo.sites <- which(phylo.sites.check == TRUE)
	return(as.numeric(names(phylo.sites)))
}

## Description

#

## Arguments

# x - an object of class DNAbin

## Source

# modified from seg.sites from ape package

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### write.pfConfig

write.pfConfig <- function(x,file,branchlengths,models,selection,alignment.length,search) {
	sink(file = file)
	cat("## ALIGNMENT FILE ##")
	cat("\n")
	cat(paste0("alignment = ",x,";"))
	cat("\n")
	cat("\n")
	cat("## BRANCHLENGTHS: linked | unlinked ##")
	cat("\n")
	cat(paste0("branchlengths = ",branchlengths,";"))
	cat("\n")
	cat("\n")
	cat("## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | beast | <list> ##")
	cat("\n")
	cat("##              for PartitionFinderProtein: all_protein | <list> ##")
	cat("\n")
	cat(paste0("models = ",models,";"))
	cat("\n")
	cat("\n")
	cat("# MODEL SELECTION: AIC | AICc | BIC #")
	cat("\n")
	cat(paste0("model_selection = ",selection,";"))
	cat("\n")
	cat("\n")
	cat("## DATA BLOCKS: see manual for how to define ##")
	cat("\n")
	cat("[data_blocks]")
	cat("\n")
	cat(paste0("pos1 = 1-",alignment.length,";"))
	cat("\n")
	cat("\n")
	cat("## SCHEMES, search: all | user | greedy ##")
	cat("\n")
	cat("[schemes]")
	cat("\n")
	cat(paste0("search = ",search,";"))
	cat("\n")
	cat("\n")
	cat("#user schemes go here if search=user. See manual for how to define.#")
	cat("\n")
	sink()
}

## Description

# A function that takes a character object, x (presumed to match name(s) of intended input
# alignments to PartitionFinder), and writes a configuration file for PartitionFinder
# based on user selected branchlength parameters, model set, selection method, and search
# method.

# Currently incapable of automatically writing codon data_blocks

## Arguments

# x - an object of class character

# branchlengths - argument taking values "linked" or "unlinked"

# models - argument taking values "all", "raxml", "mrbayes", "beast", or a vector of these
#			possible values

# selection - argument taking values "AIC", "AICc", or "BIC"

# search - argument taking values "all", "user", or "greedy"

## Source

# Andrew Frank (c) 2016

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### read.pfResults

read.pfResults <- function(file) {

	results <- read.table(
		file = file.path(file),
		sep="|",
		header = TRUE,
		fill = TRUE,
		stringsAsFactors = FALSE,
		nrows = 1,
		skip = 18,
		strip.white = TRUE,
		blank.lines.skip = FALSE
	)

	return(results$Best.Model)

}

## Description

## Arguments

## Source

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### write.MrBayes.block

write.MrBayes.block <- function(x,file,model) {

	sink(file = file)
	cat("begin mrbayes;\n")
	cat(paste0("	exe ",x,".nex;\n"))
	cat("	set autoclose=yes;\n")
	cat(paste0("	log start file=",x,".txt replace;\n"))

	if (model == "GTR") {
		cat("	lset nst=6;\n")
	} else if (model == "GTR+I") {
		cat("	lset nst=6 rates=propinv;\n")
	} else if (model == "GTR+G") {
		cat("	lset nst=6 rates=gamma;\n")
	} else if (model == "GTR+I+G") {
		cat("	lset nst=6 rates=invgamma;\n")

	} else if (model == "SYM") {
		cat("	lset nst=6;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "SYM+I") {
		cat("	lset nst=6 rates=propinv;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "SYM+G") {
		cat("	lset nst=6 rates=gamma;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "SYM+I+G") {
		cat("	lset nst=6 rates=invgamma;\n")
		cat("	prset statefreqpr=fixed(equal);\n")

	} else if (model == "HKY") {
		cat("	lset nst=2;\n")
	} else if (model == "HKY+I") {
		cat("	lset nst=2 rates=propinv;\n")
	} else if (model == "HKY+G") {
		cat("	lset nst=2 rates=gamma;\n")
	} else if (model == "HKY+I+G") {
		cat("	lset nst=2 rates=invgamma;\n")

	} else if (model == "K80") {
		cat("	lset nst=2;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "K80+I") {
		cat("	lset nst=2 rates=propinv;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "K80+G") {
		cat("	lset nst=2 rates=gamma;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "K80+I+G") {
		cat("	lset nst=2 rates=invgamma;\n")
		cat("	prset statefreqpr=fixed(equal);\n")

	} else if (model == "F81") {
		cat("	lset nst=1;\n")
	} else if (model == "F81+I") {
		cat("	lset nst=1 rates=propinv;\n")
	} else if (model == "F81+G") {
		cat("	lset nst=1 rates=gamma;\n")
	} else if (model == "F81+I+G") {
		cat("	lset nst=1 rates=invgamma;\n")

	} else if (model == "JC") {
		cat("	lset nst=1;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "JC+I") {
		cat("	lset nst=1 rates=propinv;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "JC+G") {
		cat("	lset nst=1 rates=gamma;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else if (model == "JC+I+G") {
		cat("	lset nst=1 rates=invgamma;\n")
		cat("	prset statefreqpr=fixed(equal);\n")
	} else return("ERROR")

	cat("	mcmc ngen=1000000 samplefreq=100 printfreq=1000 nchains=1 nruns=1 savebrlens=yes;\n")
	cat("	sumt relburnin = yes burninfrac = 0.1 conformat=simple;\n")
	cat("	quit;\n")
	cat("end;")

	sink()

}

## Description

# Apparently this doesn't play well with lapply or sapply, so I use it in a for-loop

## Arguments

## Source

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

adv.root.test <- function(tree,outgroup) {

	tree <- multi2di(tree, random = FALSE)
	tree <- ladderize(tree)
	outgroup.MRCA <- getMRCA(tree,outgroup)

	if (length(outgroup.MRCA) == 0) {
#		rooted.tree <- root(tree, resolve.root = TRUE, outgroup = outgroup[which(outgroup %in% tree$tip.label)])
		return(NULL)
	} else if (is.monophyletic(tree,outgroup) == TRUE) {
		rooted.tree <- reroot(
			tree,
			node.number = outgroup.MRCA,
			position = tree$edge.length[
				which(tree$edge[,2] == outgroup.MRCA)
			]/2
		)
	} else if (is.monophyletic(tree,outgroup) == FALSE) {
		rooted.tree <- root(tree, resolve.root = TRUE, outgroup = outgroup[1])
	}

	splitsB <- lapply(
		prop.part(tree),
		function(x) {
			if (any(x==1)) {tree$tip.label[x]
			} else {tree$tip.label[-x]}
		}
	)
	splitsA <- lapply(
		prop.part(rooted.tree),
		function(x) {
			if (any(x==1)) {rooted.tree$tip.label[x]
			} else {rooted.tree$tip.label[-x]}
		}
	)
	rooted.tree$node.label <- tree$node.label[match(splitsA,splitsB)]
	tree <- rooted.tree
	tree <- di2multi(tree)
	return(tree)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

adv.root <- function(tree,outgroup,keep.unrooted) {

	outgroup.MRCA <- getMRCA(tree,outgroup)

	if (length(outgroup.MRCA) == 0) {
		new.outgroup <- outgroup[which(outgroup %in% tree$tip.label)]
		if (length(new.outgroup) == 0 && keep.unrooted == TRUE) {
			return(tree)
		} else if (length(new.outgroup) == 0 && keep.unrooted == FALSE) {
			return(NULL)
		} else rooted.tree <- root(tree, resolve.root = TRUE, outgroup = new.outgroup)
	} else if (is.monophyletic(tree,outgroup) == TRUE) {
		rooted.tree <- reroot(
			tree,
			node.number = outgroup.MRCA,
			position = tree$edge.length[
				which(tree$edge[,2] == outgroup.MRCA)
			]/2
		)
	} else if (is.monophyletic(tree,outgroup) == FALSE) {
		rooted.tree <- root(tree, resolve.root = TRUE, outgroup = outgroup[1])
	}

	return(rooted.tree)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


lappend <- function (lst, ...) {
	lst <- c(lst, list(...))
	return(lst)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


extract.simple.subtree <- function(x,node,taxa) {

	subtree <- extract.clade(x,node)
	if (length(which(taxa %in% subtree$tip.label)) == 0) {
		return(NULL)
	} else
	subtree$edge.length <- NULL
	subtree$node.label <- NULL
	write.tree(subtree,file="")

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


check.species.monophyly <- function(x,taxa) {

	tip.num <- which(x$tip.label %in% taxa)
	mono.topo <- is.monophyletic(x,taxa)
	mrca <- getMRCA(x,taxa)

	mono.status <- "Monophyletic"
	support <- x$node.label[which(unique(x$edge[,1]) == mrca)]
	nonmono.subtrees.report <- NULL

	if (mono.topo == FALSE) {

		mono.status <- "Non-monophyletic"
		descendants <- unlist(Descendants(x,mrca))

		repeat {

			siblings <- list()
			for (i in seq_along(descendants)) {
				temp <- Siblings(x,descendants[i],include.self=TRUE)
				siblings <- lappend(siblings,temp)
			}
			siblings <- unique(siblings)
			names(siblings) <- lapply(siblings,mrca.phylo,x=x)

			mono.test <- lapply(
				siblings,
				function(y) {
					members <- y[which(y %in% tip.num)]
					if (! FALSE %in% (y <= length(x$tip.label))) {
						if (length(members) == 0) {
							return("Alien")
						} else if (length(members) == length(y)) {
							return("Member")
						} else return("Non-monophyetic")
					} else
						clades <- y[which(y > length(x$tip.label))]
						clade.descendants <- unlist(sapply(clades,Descendants,x=x))
						clade.members <- clade.descendants[which(clade.descendants %in% tip.num)]
						total.members <- c(members,clade.members)
						if (length(total.members) == 0) {
							return("Alien")
						} else if (length(total.members) == (length(y) - length(clades) + length(clade.descendants))) {
							return("Member")
						} else return("Non-monophyetic")
				}
			)

			nonmono.clades <- as.numeric(names(mono.test)[which(mono.test == "Non-monophyetic")])
			mono.clades <- as.numeric(names(mono.test)[c(which(mono.test == "Alien"), which(mono.test == "Member"))])

			if (TRUE %in% (nonmono.clades %in% mrca)) {
				nonmono.clades <- nonmono.clades[-which(nonmono.clades %in% mrca)]
			}

			if (length(nonmono.clades) != 0) {

				nonmono.clades.descendants <- Descendants(x,nonmono.clades,type="all")
				names(nonmono.clades.descendants) <- nonmono.clades
				if (class(nonmono.clades.descendants) != "list") {nonmono.clades.descendants <- list(nonmono.clades.descendants)}

				nonmono.subtrees <- lapply(
					nonmono.clades.descendants,
					function (y) {
						embedded.clades <- y[which(y %in% nonmono.clades)]
						if (length(embedded.clades) == 0) {
							subtree <- extract.simple.subtree(x,mrca.phylo(x,y),taxa)
							return(subtree)
						} else
							embedded.clades.descendants <- unlist(Descendants(x,embedded.clades))
							remaining.tips <- y[-c(which(y %in% embedded.clades),which(y %in% embedded.clades.descendants))]
							members <- remaining.tips[which(remaining.tips %in% tip.num)]
							if (length(members) == 0 || length(members) == length(remaining.tips)) {
								return(NULL)
							} else
								subtree <- extract.simple.subtree(x,mrca.phylo(x,y),taxa)
								return(subtree)
					}
				)
				nonmono.subtrees.report <- c(nonmono.subtrees.report,unlist(nonmono.subtrees))
				names(nonmono.subtrees.report) <- NULL

				nonmono.subtrees.mod <- gsub(nonmono.subtrees.report,pattern=";",replacement="")
				for (i in seq_along(nonmono.subtrees.mod)) {
					subtree.presence <- agrep(nonmono.subtrees.mod[-i],pattern=nonmono.subtrees.mod[i],0.0001)
					if (length(subtree.presence) != 0) {
						nonmono.subtrees.report <- nonmono.subtrees.report[-i]
					}
				}

				mono.clades.inside <- sapply(
					nonmono.clades.descendants,
					function (y) {
						mono.clades[which(mono.clades %in% y)]
					}
				)
				mono.clades.inside <- which(mono.clades %in% unique(unlist(mono.clades.inside)))
				if (length(mono.clades.inside) != 0) {mono.clades <- mono.clades[-mono.clades.inside]}

			}

			if (length(mono.clades) == 0) {break}
			descendants <- mono.clades

		}
	}

	if (is.null(nonmono.subtrees.report) && mono.topo == TRUE) {
		nonmono.subtrees.report <- "No non-monophyletic subtrees"
	} else if (is.null(nonmono.subtrees.report) && mono.topo == FALSE) {
		nonmono.subtrees.report <- "Ingroup polytomy"
	}

	nonmono.subtrees.report <- paste(nonmono.subtrees.report,collapse="")

	report <- list(mono.status,support,nonmono.subtrees.report)
	names(report) <- c("Monophyly","Support","Subtrees")
	return(report)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


color.tips <- function(x,color.list) {

	tipnames <- x$tip.label

	tip.colors <- c()
	for (i in seq_along(color.list)) {
		tip.colors[grep(pattern = names(color.list)[i],tipnames)] <- color.list[i]
	}

	return(unlist(tip.colors))

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read.galax.results <- function(file) {

	raw.results <- readLines(file)
	start.line <- which(raw.results == "Totals") + 11

	results <- read.table(
		file = file,
		header = TRUE,
		fill = TRUE,
		stringsAsFactors = FALSE,
		nrows = 1,
		skip = start.line,
		strip.white = TRUE,
		blank.lines.skip = FALSE
	)

	return(results$Ipct)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Corrects names from Alan and Emily Lemmons dataset to a more reasonable form... also
# corrects the mislabeling of JQR576 as a gilberti, and correctly lists it as a
# skiltonianus

shorten.name <- function(x) {
	short.names <- gsub(
		x$tip.label,
		pattern = "I\\d+_(\\w+)_Squamata_Scincidae_Plestiodon_(\\w+)_seq\\d",
		replacement = "\\1_\\2"
	)
	x$tip.label <- short.names

	corrected.names <- sub(
		x$tip.label,
		pattern = "JQR576_gilberti",
		replacement = "JQR576_skiltonianus"
	)
	x$tip.label <- corrected.names

	return(x)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Assumes alleles occur in numerical order, i.e. skink_seq1, skink_seq2, lizard_seq1,
# lizard seq2, etc.

pick.random.alleles <- function(tree,seed) {

	set.seed(seed)

	tips.num <- length(tree$tip.label)

	taxa.IDs <- seq(from = 1, to = tips.num, by = 2)
	selected.alleles <- sapply(
		taxa.IDs,
		function(y){
			random.num <- sample(1:2,1)
			if (random.num == 1) {
				return(y)
			} else if (random.num == 2) {
				return(y+1)
			}
		}
	)
	random.tree <- drop.tip(tree, tip = selected.alleles)
	return(random.tree)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create.pruned.tree <- function(tree,focus.taxa) {

	taxa <- tree$tip.label
	taxa.test <- focus.taxa %in% taxa

	if (FALSE %in% taxa.test) {
		return(NULL)
	} else
		dropped.taxa <- taxa[which(! taxa %in% focus.taxa)]
		pruned.tree <- drop.tip(
			tree,
			dropped.taxa)
		return(pruned.tree)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

remove.polytomy.quartet <- function(tree) {

	if (tree$Nnode < 3) {return(NULL)}
	else return(tree)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

remove.erroneous.roots <- function(tree,outgroup.taxon) {

	taxa <- tree$tip.label
	ingroup <- taxa[- which(taxa == outgroup.taxon)]

	if (outgroup.taxon %in% taxa == FALSE) {
		return(NULL)
	} else if (getMRCA(tree,taxa) == getMRCA(tree,ingroup)) {
		return(NULL)
	}	else return(tree)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read.nexus.mrbayes <- function(file,branchsupport) {
	tree <- read.nexus(file)
	if (branchsupport == TRUE) {
		majrule.tree <- tree[[1]]
	} else if (branchsupport == FALSE) {
		majrule.tree <- tree[[2]]
	}
	return(majrule.tree)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Would like to implement some kind of burn-in option to this function

# By default, MrBayes outputs a nexus file containing all the trees from the
# posterior, which is very difficult to parse by APE's read.nexus. This function
# uses bash commands on the original text file to get around reading every
# posterior tree into memory. It outputs a multiPhylo object with the number of
# trees sampled from the posterior indicated by posterior.sampling.num.

read.nexus.mrbayesposterior <- function(file,posterior.sampling.num) {

	# Get line numbers and put together bash commands to prepare for reading
	# the nexus file

	command.getdivide <- paste0(					# Put together bash command to find the
		"grep -n \"_seq\\d;\" \"",					# line where the final taxa occurs in
		file,																# taxa block
		"\" | cut -f1 -d:"
		)
	divide.line <- system(								# Run the bash command to find the line
		command.getdivide,									# where the final taxa occurs in the
		intern = TRUE												# taxa block
		)

	sampling.start <- as.numeric(					# Specify the line to start reading in
		divide.line) +											# trees, which is the dividing line
		posterior.sampling.num							# plus the posterior sampling number

	command.getbeginning <- paste0(				# Put together bash command to read
		"gsed -n 1,",												# the nexus file from the first line
		divide.line,												# to the last line of the taxa block,
		"p \"",															# using gsed
		file,
		"\""
		)
	command.gettrees <- paste0(						# Put togther bash command to read the
		"gsed -n ",													# nexus file from the sampling starting
		sampling.start,											# line, and then each line according to
		"~",																# the posterior sampling number, using
		posterior.sampling.num,							# gsed
		"p \"",
		file,
		"\""
		)

	# Execute bash commands to read in taxa block and tree blocks, as well making
	# a final end command for the nexus file

	beginningblock <- system(
		command.getbeginning,
		intern = TRUE
		)
	treeblock <- system(
		command.gettrees,
		intern = TRUE
		)
	endblock <- "end;"

	# Assemble blocks together and process trees into a DNAbin object

	assembled.nexus <- c(									# Concatenate each block read in by the
		beginningblock,											# bash commands (and the end command)
		treeblock,
		endblock)

	temporary.nexus <- tempfile()					# Establish temporary file to write to

	writeLines(														# Write each line from assembled.nexus
		assembled.nexus, 										# to the temporary file
		con = temporary.nexus
		)
	trees <- read.nexus(									# Read in the assembled nexus file as
		file = temporary.nexus)							# trees, which makes them into DNAbin
																				# objects

	# In order to identify trees by the gene they came from, paste the
	# file basename to each tree's name

	newnames <- paste(										# Create vector of new gene tree names
		basename(file),											# by pasting the basename of the tree
		names(trees),												# file to the name of each tree
		sep = "_"
		)
	names(trees) <- newnames							# Assign the new name vector to the gene
																				# trees

	return(trees)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

concat.multiPhylo <- function(trees) {

	temporary.nexus <- tempfile()

	mclapply(
		trees,
		mc.cores = 8,
		write.tree,
		file = temporary.nexus,
		tree.names = TRUE,
		append = TRUE
		)

	concat.trees <- read.tree(temporary.nexus)

	return(concat.trees)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.genetreemap <- function(trees) {

	genetreenames.col <- paste(
		names(trees),
		collapse = ","
		)
	genetreemap <- paste0(
		"{",
		genetreenames.col,
		"}"
		)
	return(genetreemap)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.phylonet.nexus <- function(
	files,										# vector of paths to tree files, if tree.type is a
														# posterior option, then vector of paths to
														# posterior tree files
	tree.names,								# vector of tree names matching the tree list
	tree.type,								# "MrBayes" | "MrBayesPosterior", "BEAST"
	outgroup,
	posterior.sampling.num,		# enumeration of how often to sample posterior trees
														# e.g. 1000 to sample 10 trees from the posterior
														# posterior trees from 10000 total trees
	allele.sampling,					# "all", "full" only include trees w/ all taxa
	allele.mapping, 					# "unphased" | "mapped" | "random"
	allele.map,
	analysis.type,						# "MP","ML","ML_CV","ML_Bootstrap","MPL"
	ret.num,									# num of reticulations allowed
	output.path								# output directory for nexus files generated
) {

	if (tree.type == "MrBayes") {

		trees <- lapply(
			files,
			read.nexus.mrbayes
		)
		names(trees) <- basename(files)

	} else if (tree.type == "MrBayesPosterior") {

		posterior.trees <- lapply(
			files,
			read.nexus.mrbayesposterior,
			posterior.sampling.num
			)
		names(posterior.trees) <- basename(files)

		genetree.map <- lapply(
			posterior.trees,
			write.genetreemap
			)

		genetree.map.concat <- paste0(
			"(",
			paste(genetreemap,collapse = ","),
			")"
			)

		trees <- concat.multiPhylo(posterior.trees)

	} else if (tree.type == "manual") {

		trees <- read.tree(files)
		names(trees) <- seq(1,length(trees))

	} else return("ERROR: BEAST trees not implemented yet")

	if (allele.sampling == "full") {
		tips.num <- sapply(
			trees,
			function (x) {
				length(x$tip.label)
			}
		)
		trees <- trees[which(tips.num == max(tips.num))]
	}

	if (allele.mapping == "random") {

		trees <- lapply(
			trees,
			shorten.name
		)

		tip.nums <- sapply(
			trees,
			function(x) length(x$tip.label))

		trees.unphased <- which(
		  tip.nums <= 0.5*max(tip.nums))

		trees.phased <- trees[-trees.unphased]

		trees.pruned <- lapply(
			trees.phased,
			pick.random.alleles,
			seed = 12345
		)

		trees <- c(trees.pruned,trees[trees.unphased])

		root.presence <- sapply(
			trees,
			function(x) outgroup %in% x$tip.label)

		trees.noRoot <- which(root.presence == FALSE)
		trees.root <- trees[-trees.noRoot]

		trees.rooted <- lapply(
			trees.root,
			root,
			outgroup = outgroup
		)

		trees <- c(trees.rooted,trees[trees.noRoot])

	} else if (allele.mapping == "mapped") {
		allele.map <- allele.map

		trees <- lapply(
			trees,
			adv.root,
			outgroup = outgroup,
			keep.unrooted = TRUE
		)

	}

	analysis <- paste0("InferNetwork_",analysis.type)

	file.name <- paste(
#		tree.type,
#		allele.sampling,
#		allele.mapping,
#		analysis.type,
		"garli_MLTrees_HighQuality_Newick",
		paste0(ret.num,"ret"),
		sep="_"
	)
	full.file.name <- file.path(
		output.path,
		paste0(file.name,".nex")
	)
	output.file.name <- paste0(file.name,".txt")

	sink(full.file.name)
	cat("#NEXUS\n\n")
	cat("BEGIN TREES;\n\n")
	for (i in seq_along(trees)) {
		cat(paste0("tree ",names(trees[i])," = "))
		cat(write.tree(trees[[i]],file=""),"\n")
	}
	cat("\nEND;\n\n")
	cat("BEGIN PHYLONET;\n\n")
	if (tree.type == "MrBayesPosterior") {
		cat(analysis, genetreemap.concat, "")
	} else cat(analysis, "(all) ")
	cat(ret.num,"-x 10","-di -pl 8")
#	if (tree.type == "BEAST") {
#		cat(" -bl")
#		if (analysis.type == "ML_Bootstrap") {
#			cat(" -ms /home/afrank/develop/msdir/ms")
#		}
#	}
	if (allele.mapping == "mapped") {cat(allele.map," ")}
	cat(output.file.name)
	cat(";")
	cat("\n\nEND;")
	sink()


}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.DendroscopeNetworks <- function(networks,commands) {

	temporary.newick <- tempfile()

	writeLines(
		networks,
		con = temporary.newick)

	commands <- paste0(
		"<open file=",temporary.newick,";> ",
		commands
		)

	command.dendroscope <- paste0(
		"/Users/Andrew/Develop/Dendroscope/Dendroscope.app/Contents/MacOS/JavaApplicationStub ",
		"-g ",
		commands)

	system(command.dendroscope)


}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DNAbin2xml <- function(x,file) {

	sink(file = file)
	for(i in seq_along(rownames(x))) {
		taxon <- rownames(x)[i]
		sequence <- paste(x[i,],collapse="")
		cat(paste0("<sequence id=\"seq_",taxon,"\" "))
		cat(paste0("taxon=\"",taxon,"\" "))
		cat("totalcount=\"4\" ")
		cat(paste0("value=\"",sequence,"\"/>"))
		cat("\n")
	}
	sink()

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# x is trees as phylo object

calculate.T1T2 <- function(x) {

	node.distances <- branching.times(x)
	T2 <- node.distances[3]/node.distances[1]
	T1 <- node.distances[2]/node.distances[1]
	values <- c(T1,T2)
	names(values) <- c("T1","T2")

	return(values)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

T1T2.table <- function(trees) {

		T1T2.values <- lapply(
			trees,
			calculate.T1T2)
		T1T2.df <- as.data.frame(t(
			do.call(cbind.data.frame,T1T2.values)))

		unique.quartets <- unique.multiPhylo(trees)
		raw.topology.ids <- attr(
			unique.quartets,"old.index")
		raw.topology.counts <- table(raw.topology.ids)
		names(unique.quartets) <- names(raw.topology.counts)

		topology.counts <- sort(
			raw.topology.counts,
			decreasing = TRUE)
		topology.order <- names(topology.counts)

		names(topology.counts) <- mapvalues(
			names(topology.counts),
			from = topology.order,
			to = c("A","B","C"))
		topology.ids <- mapvalues(
			raw.topology.ids,
			from = topology.order,
			to = c("A","B","C"))

		T1T2.df$Topology.IDs <- topology.ids

		names(unique.quartets) <- mapvalues(
			names(unique.quartets),
			from = topology.order,
			to = c("A","B","C"))
		attr(unique.quartets, "old.index") <- NULL
		unique.quartets <- unique.quartets[
			c("A","B","C")]

		returned.list <- list(
			T1T2.df,
			topology.counts,
			unique.quartets)
		names(returned.list) <- c(
			"T1T2.table",
			"Topology.counts",
			"Unique.quartets")

		return(returned.list)

	}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# x is the combined t1 or t2 dataframe, y indicates whether you chose t1 or t2

calc.T1T2.signifiance <- function(x,y) {

	if (y == "T1") {
		aov.out <- aov(T1 ~ Topologies, data=x)
	} else if (y == "T2") {
		aov.out <- aov(T2 ~ Topologies, data=x)
	}
	p.values <- TukeyHSD(aov.out)$Topologies[10:12]
	sig.p.values <- which(p.values <= 0.05)

	if (length(sig.p.values) == 0) {
		sigdiffs <- c("x","x","x")
	} else if (length(sig.p.values) == 3) {
		sigdiffs <- c("x","y","z")
	} else if (length(sig.p.values) == 1) {
		if (sig.p.values == 1) {
			sigdiffs <- c("x","y","x,y")
		} else if (sig.p.values == 2) {
			sigdiffs <- c("x","x,y","y")
		} else if (sig.p.values == 3) {
			sigdiffs <- c("x,y","x","y")
		}
	} else if (length(sig.p.values) == 2) {
		if (paste(sig.p.values,collapse="") == "12") {
			sigdiffs <- c("x","y","y")
		} else if (paste(sig.p.values,collapse="") == "13") {
			sigdiffs <- c("x","y","x")
		} else if (paste(sig.p.values,collapse="") == "23") {
			sigdiffs <- c("x","y","y")
		}
	}
	return(sigdiffs)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot.T1T2 <- function(T1T2.object,file,taxa.annotations.df, outgroup.taxon) {

	display.trees <- lapply(
		T1T2.object$Unique.quartets,
		drop.tip,
		outgroup.taxon)
	names(display.trees) <- paste(
		names(T1T2.object$Unique.quartets),
		as.character(T1T2.object$Topology.counts),
		sep = ", ")

	anno.df <- taxa.annotations.df[
		which(taxa.annotations.df$tip.label %in%
			display.trees$tip.label),]

	trees <- lapply(
		names(display.trees),
		function(x) {
			tree <- ggtree(
				display.trees[[x]],
				branch.length = "none") %<+%
				anno.df +
				geom_tippoint(
					aes(color = species),
					size = 5,
					alpha = 0.75) +
				ggtitle(names(display.trees[x]))
			return(tree)
		}
		)
	tmp.tree <- trees[[1]] +
		theme(
			legend.position = "bottom",
			legend.box = "horizontal") +
		guides(col = guide_legend(nrow = 1))
	trees.legend <- g_legend(tmp.tree)

	T1T2.melted <- melt(T1T2.object$T1T2.table)
	boxplot <- ggplot(
			T1T2.melted,
			aes(factor(Topology.IDs), value)) +
		geom_boxplot() +
		facet_grid(. ~ variable) +
		ggtitle(paste(display.trees$tip.label, collapse = ", "))

	pdf(file, paper = "USr", width = 11)
	grid.arrange(
		trees[[1]], trees[[2]], trees[[3]], trees.legend, boxplot,
		ncol = 2,
		layout_matrix = cbind(c(1,2,3,4), c(5,5,5,4)),
		widths = c(0.25,1),
		heights = c(1,1,1,0.3))
	dev.off()

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
