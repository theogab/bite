bindTwoTrees <- function(acceptor, donnor, nb_spec){
		# Binds two different trees and keeps them ultrametrice
		# User should provide number acceptor and donnor trees and number of species
		#
		# Args:
		# 	acceptor:	acceptor tree (to which to bind)
		#	donnor:  	donnor tree (which to bind)
		#	nb_spec:	total number of species
		#
		# Returns:
		#	An object of class phylo (a new binded tree).

		donnor_stem_age <- branching.times(donnor)[1] + donnor$root.edge # total age of the tree
		br.ti <- branching.times(acceptor)
		labels <- c(acceptor$tip.label, acceptor$node.label)
		
		# if the total age of the donor tree is larger than the time of mrca of the acceptor tree, then we bind both trees at their root
		if (br.ti[1] < donnor_stem_age) samp <- acceptor$node.label[1] 
		else{#otherwise we randomly choose an edge to plug the donnor tree, that corresponds to the age of the donor tree
		br.ti <- c(br.ti, setNames(rep(0, Ntip(acceptor)), acceptor$tip.label))
		edges <- cbind(labels[acceptor$edge[,1]], labels[acceptor$edge[,2]])
		samp_vect <- edges[br.ti[edges[,2]]<donnor_stem_age & br.ti[edges[,1]]>donnor_stem_age, 2]
		samp <- sample(samp_vect, 1)
		}
		
		k <- (which(labels==samp))
		#we bind it at the right position (branching time of the edge selected minus age of the donor tree) so that we have a ultrametric tree
		new_tree <- bind.tree(acceptor, donnor, where = k, position = donnor_stem_age - br.ti[samp])
		new_tree$node.label[is.na(new_tree$node.label)] <- paste("n", nb_spec, sep="") #we rename the node label that has been erased
		names(new_tree$edge.length) <- NULL
		if (is.ultrametric(new_tree)) return(new_tree)
		else print("Error") #this shouldn't happen
		
}
                 

sortBindAllTrees <- function(list_of_trees, param){
		# Binds all different trees together
		# User should provide a list of tree and parameters (?)
		#
		# Args:
		# 	list_of_trees:	list of trees to bind
		#	param:  		parameters of ... 
		#
		# Returns:
		#	An object of class phylo (a new binded tree).

		if (length(list_of_trees)==1){
		final_tree <- read.tree(text=write.tree(list_of_trees[[1]][[1]]))
		break_node <- grep("s1", final_tree$node.label)+Ntip(final_tree)
		param <- cbind(param, break_node=break_node)
		return(list(tree=final_tree, param=param)) #if there is only one tree
		}
		root_times <- c()
		nb_spec <- 0
		br_times <- c()
		for (i in 1:length(list_of_trees)) {#we define the root edge as the difference between the total age of the tree and the time of mrca
		list_of_trees[[i]][[1]]$root.edge <- list_of_trees[[i]][[2]]-branching.times(list_of_trees[[i]][[1]])[1]
		root_times <- c(root_times, setNames(list_of_trees[[i]][[2]], i))
		br_times <- c(br_times, setNames(branching.times(list_of_trees[[i]][[1]])[1], i))
		nb_spec <- nb_spec+Ntip(list_of_trees[[i]][[1]])
		}
		ordered <- sort(root_times, decreasing = TRUE) #we order the trees by decreasing total age
		print(sort(br_times, decreasing = TRUE))
		final_tree <- list_of_trees[[as.numeric(names(ordered)[1])]][[1]]
		for (i in 1:(length(ordered)-1)){# we begin by the most ancient and sequentially bind the young trees to it
		final_tree <- bindTwoTrees(final_tree, list_of_trees[[as.numeric(names(ordered)[i+1])]][[1]], nb_spec+i+1)
		}
		final_tree <- read.tree(text=write.tree(final_tree))
		print(ordered)
		if (!grepl("s", final_tree$node.label[1])){
		tmp <- final_tree$node.label[1]
		changing_lab <- paste("s", names(ordered)[1], sep="")
		final_tree$node.label[1] <- changing_lab
		final_tree$node.label[duplicated(final_tree$node.label)] <- tmp
		param[changing_lab, "root_times"] <- NA
		param[changing_lab, "mrca_times"] <- branching.times(final_tree)[1]
		print("replaced!")
		}
		break_node <- c()
		for (i in 1:(length(list_of_trees))) break_node <- c(break_node, grep(paste("s", i, sep=""), final_tree$node.label)+Ntip(final_tree))
		param <- cbind(param, break_node=break_node)
		return(list(tree=final_tree, param=param))
}
