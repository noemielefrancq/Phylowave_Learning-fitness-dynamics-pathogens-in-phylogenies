## Small tool functions
#' @title Compute number of descendants of all nodes (internal and external)
#' @import ape
#' @param tree Tree (object of class "phylo" or "multiPhylo")
#' @return Vector containing the number of descendants of each node
#' @export
number.descendants.all.nodes <-  function(tree) {
  res = numeric(max(tree$edge))
  res[1:Ntip(tree)] <- 1L
  for (i in postorder(tree)) {
    tmp = tree$edge[i,1]
    res[tmp] = res[tmp] + res[tree$edge[i, 2]]
  }
  return(res)
}

#' @title Output the descendants of all nodes (internal and external)
#' @import ape
#' @param tree Tree (object of class "phylo" or "multiPhylo")
#' @return Vector containing the descendant of each node
#' @export
descendants.all.nodes <-  function(tree) {
  res = as.list(1:(2*Ntip(tree)-1))
  for (i in postorder(tree)) {
    tmp = tree$edge[i,1]
    res[[tmp]] = c(res[[tmp]], res[[tree$edge[i, 2]]])
  }
  return(res)
}

#' @title Compute the path from the tree root to all nodes
#' @import ape
#' @param tree Tree (object of class "phylo" or "multiPhylo").
#' @return List with each element being a vector containing the path fromß the root to the node.
#' @export
path.root.to.all.nodes <-  function(tree) {
  res = as.list(1:(2*Ntip(tree)-1))
  for (i in rev(postorder(tree))) {
    tmp_p = tree$edge[i,1]
    tmp_c = tree$edge[i,2]
    res[[tmp_c]] = c(res[[tmp_p]], res[[tmp_c]])
  }
  return(res)
}

#' @title Compute cumulative index of all nodes
#' @import ape
#' @param tree Tree (object of class "phylo" or "multiPhylo")
#' @param index Index of each node
#' @param time Sampling (or inferred) time of each node
#' @return List containing one vector per node, with the  descendant of each node
#' @export
cumulative.index.all.nodes <-  function(tree, index, time) {
  index[which(is.na(index) | is.nan(index) | is.infinite(index))] = 0
  res = index
  for (i in rev(postorder(tree))) {
    tmp_p = tree$edge[i,1]
    tmp_c = tree$edge[i,2]
    res[tmp_c] = res[tmp_p] + res[tmp_c] * (time[tmp_c] - time[tmp_p])
  }
  return(res)
}

test.nested.groups = function(tree, root_groups, metadata){
  potential_parents=rep(root_groups, each=length(root_groups)-1)
  potential_descendants=c()
  for (i in 1: length(root_groups)) potential_descendants=c(potential_descendants, root_groups[-i])
  nested=c()
  for (i in 1: length(potential_parents)){
    nested=c(nested, potential_descendants[i] %in% getDescendants(tree, potential_parents[i]))
  }
  parents=potential_parents[nested]
  descendants=potential_descendants[nested]
  timing_parents=metadata$time[parents]
  nested_dataframe=as.data.frame(cbind(parents, descendants, timing_parents))
  if (length(timing_parents)>1) nested_dataframe =nested_dataframe[order(-nested_dataframe$timing_parents),]#start by the most recent parent
  nested_dataframe = nested_dataframe[!duplicated(nested_dataframe), ]
  return(nested_dataframe)
}

update.root.groups = function(tree, root_groups, nroot, parent, descendant){
  other_roots=root_groups[-which(root_groups==parent)]
  #get a new ancestral group root : sister node of the descendent one
  path_desc=nodepath(tree, from=parent, to=descendant)
  parent_new_root_desc=tail(path_desc,2)[1]#we get the parent of the descendant root
  desc_parent_new_root_desc=tree$edge[,2][tree$edge[,1]==parent_new_root_desc]
  new_root_desc=desc_parent_new_root_desc[desc_parent_new_root_desc!=descendant]
  #get a new ancestral group root : sister node of the ancestral one
  path_ansc=nodepath(tree, from=nroot, to=parent)
  parent_new_root_ansc=tail(path_ansc,2)[1]#we get the parent of the ancestral root
  desc_parent_new_root_ansc=tree$edge[,2][tree$edge[,1]==parent_new_root_ansc]
  new_root_ansc=desc_parent_new_root_ansc[desc_parent_new_root_ansc!=parent]
  root_groups=c(other_roots, new_root_desc, new_root_ansc)
  root_groups=root_groups[!duplicated(root_groups)]
  return(root_groups)
}

compute.root.structural.groups=function(tree, metadata, genetic_distance_mat, nroot, nodes_height){
  mds=cmdscale(genetic_distance_mat)
  mds=mds[-nroot,] # remove root of the tree
  threshold_date=ceiling(min(nodes_height)) #we only keep older nodes, where computation of the index might be unreliable
  nb_nodes=sum(nodes_height<(threshold_date))
  while (nb_nodes<=2){
    threshold_date=threshold_date+1
    nb_nodes=sum(nodes_height>(threshold_date-1) & nodes_height<threshold_date)#while we don't get at least 2 nodes in the same year
  }
  oldest_nodes=names(nodes_height[nodes_height<threshold_date])
  reduced_mat_edge=tree$edge[tree$edge[,1] %in% oldest_nodes,]
  nodes_and_tips_kept=unique(as.vector(reduced_mat_edge))
  mds=mds[rownames(mds) %in% nodes_and_tips_kept,]
  diag_best_nb_groups=fviz_nbclust(mds, kmeans, method = "silhouette", k.max=min(10,dim(mds)[1]-1))#we find the best number of groups
  best_nb_groups=which.max(diag_best_nb_groups$data$y)
  groups=kmeans(mds, centers = best_nb_groups)
  root_groups=c()
  for (i in 1:length(groups$size)){
    nodes_and_tips_group=names(groups$cluster[groups$cluster==i])
    time_nodes_group=nodes_height[names(nodes_height) %in% nodes_and_tips_group]
    root=names(which.min(time_nodes_group))
    root_groups=c(root_groups, root)
  }
  root_groups=as.numeric(root_groups)
  #quality check
  nested_dataframe=test.nested.groups(tree, root_groups, metadata)
  while(dim(nested_dataframe)[1]>0){ # while there are nested groups
    parent=nested_dataframe$parents[1]
    descendant=nested_dataframe$descendants[1]
    if (dim(nested_dataframe)[1] ==2 & nested_dataframe$parents[1] == nested_dataframe$parents[2]){
      if (length(nodepath(tree, from = nested_dataframe$parents[1], to =nested_dataframe$descendants[1])) == 2 & length(nodepath(tree, from = nested_dataframe$parents[2], to =nested_dataframe$descendants[2])) == 2){
        root_groups=root_groups[-which(root_groups==nested_dataframe$parents[1])]
        nested_dataframe=test.nested.groups(tree, root_groups, metadata)
      }
    }
    if (dim(nested_dataframe)[1]>0){
      root_groups=update.root.groups(tree,root_groups, nroot, parent, descendant)
      nested_dataframe=test.nested.groups(tree, root_groups, metadata)
    }
  }
  return(root_groups)
}

#' Merge groups
#' @param timed_tree Hamming distance matrix
#' @param metadata presence/absence matrix
#' @param initial_splits potential splits found previously
#' @param group_count_threshold Minimum number of sequences per group
#' @param group_freq_threshold Minimum frequency of group, at least once in the time series
#' @return Groups
#' @export
merge.groups = function(timed_tree, metadata, structural_splits = NULL, initial_splits, group_count_threshold, group_freq_threshold){
  groups = numeric(nrow(metadata))
  metadata_tmp = metadata
  
  ## For each selected node, find the offsping and label it as being part of the new group 
  splits = unique(c(initial_splits, structural_splits))
  for(n in 1:length(splits)){
    ## Mrca of the group
    mrca = splits[n]
    
    ## Labels of all tips and nodes in the tree
    new_group = rep(0, nrow(metadata_tmp))
    new_group[c(mrca, getDescendants(timed_tree, mrca))] = splits[n] ## Members of the group
    
    ## Add column with new info
    metadata_tmp = cbind(metadata_tmp, new_group)
    colnames(metadata_tmp) = c(colnames(metadata_tmp)[1:(ncol(metadata_tmp)-1)], splits[n])
  }
  
  ## From the above columns, create new labels: 
  # Group_nodes: which nodes were included in that group
  metadata_tmp$group_nodes = apply(metadata_tmp[,which(is.na(match(colnames(metadata_tmp), splits)) == F)], MARGIN = 1, function(x)paste0(x,collapse = '_'))
  # Group: just a numbering of groups
  metadata_tmp$group = as.numeric(factor(apply(metadata_tmp[,which(is.na(match(colnames(metadata_tmp), splits)) == F)] == 0, MARGIN = 1, function(x)paste0(x,collapse = '_'))))
  root = as.character(max(metadata_tmp$group))
  if(is.null(structural_splits) == F){
    structural_roots=as.character(seq(as.numeric(root)-length(structural_splits),as.numeric(root)-1))
  }else{
    structural_roots = NULL
  }
  ## Count groups, overall and by time
  # Overall
  count_overall = table(metadata_tmp$group)
  # Remove the root from this list
  count_overall = count_overall[-which(names(count_overall) %in% c(root, structural_roots))]
  
  # By time
  count_time = freq_time = matrix(NA, 
                                  nrow = length(unique(metadata_tmp$group)), 
                                  ncol = length(unique(round(metadata_tmp$time))))
  groups = unique(metadata_tmp$group)
  years = sort(unique(round(metadata_tmp$time)), decreasing = F)
  
  # For each year, fill the matrices
  for(i in 1:length(years)){
    tmp = table(factor(metadata_tmp$group[which(round(metadata_tmp$time) == years[i])], levels = as.numeric(groups)))
    count_time[,i] = tmp
  }
  for(i in 1:length(groups)){freq_time[i,] = count_time[i,]/colSums(count_time)}
  rownames(count_time) = rownames(freq_time) = groups
  
  # Remove the root from this matrix, to avoid removing it
  freq_time = freq_time[-which(rownames(freq_time) %in% c(root, structural_roots)),]
  
  ## Find the first groups to remove
  group_to_remove = NULL
  if(length(years) == 1){
    tmp = names(which(length(which((freq_time>group_freq_threshold) == 0)) == length(years)))
  }else{
    tmp = names(which(apply(freq_time>group_freq_threshold, MARGIN = 1, function(x)length(which(x == 0))) == length(years)))
  }
  if(min(count_overall) < group_count_threshold){
    # By count
    group_to_remove = names(which.min(count_overall))
  }else if (length(tmp) > 0){
    # By freq
    group_to_remove = min(tmp)
  }
  
  k=0
  while(length(group_to_remove) > 0 & !is.null(group_to_remove)){
    # Edge group label
    edges = timed_tree$edge
    edge_group = timed_tree$edge 
    edge_group[,1] = as.numeric(metadata_tmp$group[edges[,1]])
    edge_group[,2] = as.numeric(metadata_tmp$group[edges[,2]])
    edge_group_nodes = timed_tree$edge
    edge_group_nodes[,1] = metadata_tmp$group_nodes[edges[,1]]
    edge_group_nodes[,2] = metadata_tmp$group_nodes[edges[,2]]
    edge_group_times = timed_tree$edge
    edge_group_times[,1] = metadata_tmp$time[edges[,1]]
    edge_group_times[,2] = metadata_tmp$time[edges[,2]]
    
    ## Find closest group in time and replace the group to remove by it
    tmp = matrix(edge_group[which(edge_group[,2] == group_to_remove | edge_group[,1] == group_to_remove),], ncol = 2)
    tmp_times = matrix(edge_group_times[which(edge_group[,2] == group_to_remove | edge_group[,1] == group_to_remove),], ncol = 2)
    tmp_times = matrix(tmp_times[which(tmp[,1] - tmp[,2] != 0),], ncol = 2)
    tmp = matrix(tmp[which(tmp[,1] - tmp[,2] != 0),], ncol = 2)
    tmp = tmp[which.min(tmp_times[,2] - tmp_times[,1]),]
    tmp = tmp[tmp != group_to_remove]
    metadata_tmp$group[which(metadata_tmp$group ==  group_to_remove)] = tmp
    metadata_tmp$group_nodes[which(metadata_tmp$group ==  group_to_remove)] = tmp
    
    ## Update count and freq matrices
    freq_time[which(rownames(freq_time) == tmp),] = freq_time[which(rownames(freq_time) == tmp),] + freq_time[which(rownames(freq_time) == group_to_remove),]
    freq_time = freq_time[-which(rownames(freq_time) == group_to_remove),]
    count_overall[which(names(count_overall) == tmp)] = count_overall[which(names(count_overall) == tmp)] + count_overall[which(names(count_overall) == group_to_remove)]
    count_overall = count_overall[-which(names(count_overall) == group_to_remove)]
    
    ## Find new group to remove, if any
    group_to_remove = NULL
    if(length(years) == 1){
      tmp = names(which(length(which((freq_time>group_freq_threshold) == 0)) == length(years)))
    }else{
      print(dim(freq_time))
      tmp = names(which(apply(freq_time>group_freq_threshold, MARGIN = 1, function(x)length(which(x == 0))) == length(years)))
    }
    if(min(count_overall) < group_count_threshold){
      # By count
      group_to_remove = names(which.min(count_overall))
    }else if (length(tmp) > 0){
      # By freq
      group_to_remove = min(tmp)
    }
    # print(group_to_remove)
    k=k+1
  }
  # ##delete smaller groups 
  # for (i in 1:max(metadata_tmp$group)){#remove groups that are too small
  #   num_samples_per_group=sum(metadata_tmp$group==i, na.rm = T)
  #   if (num_samples_per_group<group_count_threshold){
  #     metadata_tmp$group[metadata_tmp$group==i]=NA
  #   }
  # }
  
  ## Rename groups so that they follow each other 
  metadata_tmp$group = as.numeric(factor(metadata_tmp$group, levels = sort(as.numeric(unique(metadata_tmp$group)))))
  
  ## Create lineage final tree
  tree_tmp = timed_tree
  groups = metadata_tmp$group
  edges = tree_tmp$edge
  edge_group = tree_tmp$edge
  edge_group[,1] = as.numeric(metadata_tmp$group[edges[,1]])
  edge_group[,2] = as.numeric(metadata_tmp$group[edges[,2]])
  
  edge_group_nodes= tree_tmp$edge
  edge_group_nodes[,1] = metadata_tmp$group_nodes[edges[,1]]
  edge_group_nodes[,2] = metadata_tmp$group_nodes[edges[,2]]
  
  to_remove = which(edge_group[,1] - edge_group[,2] == 0)[1]
  k=1
  while(!is.na(to_remove)){
    ## Find the edge to remove
    tmp = tree_tmp$edge[to_remove,] 
    
    # Collapse that edge (i.e. remove node tmp[2])
    tree_tmp$edge[which(tree_tmp$edge[,1] == tmp[2]),1] = tmp[1]
    tree_tmp$edge = tree_tmp$edge[-to_remove,]
    edge_group = edge_group[-to_remove,]
    edge_group_nodes = edge_group_nodes[-to_remove,]
    
    # Next edge to remove: edges where both nodes are from the same group
    to_remove = which(edge_group[,1] - edge_group[,2] == 0)[1]
    k=k+1
  }
  
  ## Add a root to the tmp tree
  tree_tmp$edge = rbind(c(tree_tmp$edge[1,1], 1), tree_tmp$edge)
  edge_group = rbind(c(edge_group[1,1], NA), edge_group)
  edge_group_nodes = rbind(c(edge_group_nodes[1,1], NA), edge_group_nodes)
  
  ## Find "tips", i.e. nodes than do not have offspring, i.e. nodes present only once in tree_tmp$edge
  tips = as.numeric(names(table(tree_tmp$edge))[which(table(tree_tmp$edge) == 1)])
  tree_tmp$tip.label = edge_group[match(tips, tree_tmp$edge[,2]),2]
  tree_tmp$edge[match(tips, tree_tmp$edge)] =  1:length(tips)
  nodes = unique(tree_tmp$edge[which(tree_tmp$edge > length(tips))])
  tree_tmp$Nnode = length(nodes)
  
  ## Renumber nodes
  for(i in 1:tree_tmp$Nnode){
    tree_tmp$edge[which(!is.na(match(tree_tmp$edge, nodes[i])))] = as.integer(length(tips)+i)
  }
  
  ## Remove edge lengths from tree abject
  tree_tmp = tree_tmp[-which(names(tree_tmp) == "edge.length")]
  if(length(which(names(tree_tmp) == "node.label")) > 0) { # Remove node labels, if any
    tree_tmp = tree_tmp[-which(names(tree_tmp) == "node.label")]
  }
  
  ## Make sure tree_tmp is a tree object
  attributes(tree_tmp) = list('class' = 'phylo', 'names' = names(tree_tmp))
  tree_tmp = as.phylo(tree_tmp)
  # checkValidPhylo(tree_tmp) ## Check everything is fine
  
  ## Compute state of each node
  states = as.character(edge_group[match(1:(length(tips)+tree_tmp$Nnode), tree_tmp$edge)])
  names(states) = 1:(length(tips)+tree_tmp$Nnode)
  
  return(list('groups' = metadata_tmp$group,
              'lineage_tree' = tree_tmp,
              'tip_and_nodes_groups' = states))
}

#' Sub function to run the generalised additive model (gam)
run_gam_model = function(metadata_sub, clade_data, weights_sub, k_smooth, coefs, sp, log_y = T){
  ## Set copy of the clade_data to set 'by' in model
  clade_data_group = clade_data
  
  if(length(dim(clade_data)) == 1){stop('Wrong clade_data format')}
  
  if(log_y == TRUE){
    ## Write formula
    f = paste("log(metadata_sub$index) ~ s(clade_data$time, bs = 'cr', k = k_smooth)") ## log(y)
    if(ncol(clade_data) > 1){
      for(c in 2:ncol(clade_data)){
        clade_data[,c] = clade_data[,c]*metadata_sub$time
        clade_data_group[,c] = as.numeric(clade_data[,c]>0)
        f = paste0(f, "+ s(clade_data[,", c, "], bs = 'cr', k = k_smooth, by = clade_data_group[,", c, "])", collapse = '')
      }
    }
    f = as.formula(f)
    mod = try(mgcv::bam(f, weights = weights_sub, method = 'GCV.Cp', coef = coefs, sp = sp, nthreads=1, family = gaussian(link = identity)), silent=TRUE)
  }else if (log_y == FALSE){
    ## Write formula
    f = paste("metadata_sub$index ~ s(clade_data$time, bs = 'cr', k = k_smooth)") ## No log(y)
    if(ncol(clade_data) > 1){
      for(c in 2:ncol(clade_data)){
        clade_data[,c] = clade_data[,c]*metadata_sub$time
        clade_data_group[,c] = as.numeric(clade_data[,c]>0)
        f = paste0(f, "+ s(clade_data[,", c, "], bs = 'cr', k = k_smooth, by = clade_data_group[,", c, "])", collapse = '')
      }
    }
    f = as.formula(f)
    mod = try(mgcv::bam(f, weights = weights_sub, method = 'fREML', coef = coefs, sp = sp, nthreads=1, family = gaussian(link = log)), silent=TRUE)
  }
  
  return(mod)
}

#' Sub function to test a node
test_node_gam = function(i, metadata_sub, metadata_sub_rownames, 
                         nodes_to_test_tmp, descendants_all_nodes, 
                         min_group_size, clade_membership_recorded, 
                         clade_data_null, weights_sub, k_smooth, mod_null, log_y = T){
  sub = descendants_all_nodes[[as.numeric(nodes_to_test_tmp[i])]]
  sub = sub[which(is.na(match(as.character(sub), metadata_sub_rownames)) == F)]
  l_sub = length(sub)
  
  if(l_sub >= min_group_size & nrow(metadata_sub)-l_sub >= min_group_size){
    m = match(metadata_sub_rownames, as.character(sub))        
    clade_membership_new = rep(NA, nrow(metadata_sub))
    clade_membership_new[which(is.na(m) == T)] = 0
    clade_membership_new[which(is.na(m) == F)] = 1
    
    if(is.null(clade_membership_recorded)){ ## No previously recorded node
      ## Add new node to the clade_data_null, to test the fit
      clade_data_tmp = cbind(clade_data_null, clade_membership_new)
      colnames(clade_data_tmp) = c(colnames(clade_data_null), nodes_to_test_tmp[i])
      
      ## Compute groups
      groups = clade_membership_new
      
      ## Fit model
      if(k_smooth != -1){
        mod = run_gam_model(metadata_sub = metadata_sub, clade_data = clade_data_tmp, weights_sub = weights_sub, k_smooth = k_smooth, 
                            coefs = c(mod_null$coefficients, rep(0, k_smooth)), sp = NULL, log_y = log_y)
      }else{
        mod = run_gam_model(metadata_sub = metadata_sub, clade_data = clade_data_tmp, weights_sub = weights_sub, k_smooth = k_smooth, 
                            coefs = NULL, sp = NULL, log_y = log_y)
      }
    } else if(!is.null(clade_membership_recorded)){
      ## Add new node to the clade_data_null, to test the fit
      clade_data_tmp = cbind(clade_data_null, clade_membership_new)
      colnames(clade_data_tmp) = c(colnames(clade_data_null), nodes_to_test_tmp[i])
      
      ## Compute groups
      groups = apply(clade_data_tmp[,2:ncol(clade_data_tmp)], MARGIN = 1, function(x)paste0(x,collapse = '_'))
      
      if(all(table(groups) >= min_group_size)){
        ## Fit model
        if(k_smooth != -1){
          mod = run_gam_model(metadata_sub = metadata_sub, clade_data = clade_data_tmp, weights_sub = weights_sub, k_smooth = k_smooth, 
                              coefs = c(mod_null$coefficients, rep(0, k_smooth)), sp = NULL, log_y = log_y)
        }else{
          mod = run_gam_model(metadata_sub = metadata_sub, clade_data = clade_data_tmp, weights_sub = weights_sub, k_smooth = k_smooth, 
                              coefs = NULL, sp = NULL, log_y = log_y)
        }
      }else{ ## Some group(s) are too small
        mod =  NULL
      }
    }
    if("try-error" %in% class(mod) | is.null(mod)) {
      return(list('mod' = NA,
                  'groups' = NA))
    }else{
      return(list('mod' = mod,
                  'groups' = groups))
    }
  }else{
    return(list('mod' = NA,
                'groups' = NA))
  }
}


#' Find groups defined by index dynamics
#'
#' This code is fitting a generalised additive model (gam) to the index dynamics, trying to find the smallest number of lineages to explain best the index dynamics.
#' Model targets are AIC, explained deviance or maximum number of groups
#' The script is run until a target is met, and the user can choose to record all the different states (keep_track), or not
#'
#' Potential nodes can be filtered out by support and threshold_node_support. This can be a minimum number of mutations on the branch leading to the node, or a bootstrap support.
#'
#' @import mgcv
#' @importFrom phytools getDescendants getParent
#' @param timed_tree Timed tree
#' @param metadata Metadata dataframe
#' @param node_support Numeric value of support for each node (e.g. mutations on the branch leading to the node, or bootstrap support)
#' @param threshold_node_support Threshold on the node support for the nodes to be considered in the detection algorithm
#' @param time_window_initial Time in years of when to start the screen 
#' @param time_window_increment Window
#' @param min_descendants_per_tested_node To start the analysis, start from nodes that have this minimum number of sequences
#' @param min_group_size Minimum group size, when creating a new potential split
#' @param p_value_smooth p-value threshold to consider a new spline as significant.
#' @param stepwise_deviance_explained_threshold  Change in explained deviance, recommend 0.005 (0.5 percent): By how much the new split should increase the percentage of deviance explained, at the very least.
#' @param stepwise_AIC_threshold AIC threshold, commonly 7: By how much the new split should increase the goodness of fit, at the very least.
#' @param weight_by_time NULL or numeric (in years): Size of the window of time on which to compute the weights
#' @param weighting_transformation NULL, inv_freq, inv_sqrt, or inv_log: Type of weighting to use.
#' @param k_smooth Number of knots for the splines
#' @param plot_screening TRUE or FALSE, to plot the discovered dynamics as soon as they're found
#' @param parallelize_code TRUE or FALSE, to parallelise code
#' @param number_cores Integer, if parallelised, how many cores to use
#' @param max_groups_found Integer, Maximum number of groups to find
#' @param keep_track TRUE or FALSE, whether user wants to output all the steps
#' @param log_y TRUE or FALSE, whether user wants to fit the index values on a log scale (TRUE), or with a log link (FALSE). The algorithm with the log link is much slower.
#' @return Potential split nodes at target, and if keep_track==T, all the different nodes at each iteration
#' @export
find.groups.by.index.dynamics = function(timed_tree, 
                                         metadata, 
                                         node_support = NULL,
                                         threshold_node_support = 0.5,
                                         time_window_initial = 1, 
                                         time_window_increment = 0.5, 
                                         min_descendants_per_tested_node = 100, 
                                         min_group_size = 25,
                                         p_value_smooth = 0.01,
                                         stepwise_deviance_explained_threshold = 0.005,
                                         stepwise_AIC_threshold = 7,
                                         weight_by_time = NULL,
                                         weighting_transformation = c(NULL, 'inv_freq', 'inv_sqrt', 'inv_log'),
                                         k_smooth = 5,
                                         plot_screening = F,
                                         parallelize_code = F,
                                         number_cores = 2,
                                         max_groups_found = 100,
                                         keep_track = F,
                                         log_y = T){
  ## Store paths from root to all tips and nodes
  paths_root_to_all_nodes = path.root.to.all.nodes(timed_tree)
  
  ## Store descendants for all nodes
  descendants_all_nodes = descendants.all.nodes(timed_tree)
  
  ## Store number of descendants for all nodes
  nb_descendants_all_nodes = number.descendants.all.nodes(timed_tree)
  
  ## Remove potential NAs
  idx = which(is.infinite(metadata$index) | is.nan(metadata$index))
  if(length(idx)>0) metadata = metadata[-idx,]
  
  ## Remove potential 1s
  idx = which(metadata$index >= 1)
  if(length(idx)>0) metadata$index[idx] = 0.99
  
  ## Row names
  metadata_rownames = metadata$ID
  
  ## List all nodes
  nodes_to_test = metadata_rownames[which(metadata$is.node == 'yes')]
  ## Remove the root 
  nodes_to_test = nodes_to_test[which(nodes_to_test != timed_tree$Nnode + 2)]
  
  ## Remove nodes with node support below threshold
  if(!is.null(node_support)){
    idx = which(as.numeric(node_support) < threshold_node_support)
    nodes_to_test = nodes_to_test[which(is.na(match(nodes_to_test, idx + timed_tree$Nnode+1)))] 
  }
  
  ## Remove nodes that have less than minimum number of descendants
  nodes_to_test = nodes_to_test[which(nb_descendants_all_nodes[as.numeric(nodes_to_test)] >= min_descendants_per_tested_node)]
  
  ## Set weights, if user selected weighting by time
  weights = NULL
  if(!is.numeric(weight_by_time) & !is.null(weight_by_time)) stop('Weight_by_time is not NULL or numeric')
  if(!is.null(weight_by_time)){
    weights = rep(NA, length(metadata$time))
    window_weights = seq(round(min(metadata$time), digits = 1)-1e4, round(max(metadata$time), digits = 1)+1e4, weight_by_time)
    freq = hist(metadata$time, breaks = window_weights, plot = F, include.lowest = T)
    freq$counts[which(freq$counts == 0)] = NA
    if(length(weighting_transformation) != 1) stop('Weighting_transformation parameter is not of length 1, please supply one weighting transformation')
    freqs_weighted = switch(weighting_transformation, 
                            'inv_freq' = {sum(freq$counts, na.rm = T)/freq$counts}, 
                            'inv_sqrt' = {sum(sqrt(freq$counts), na.rm = T)/sqrt(freq$counts)}, 
                            'inv_log'  = {sum(log(freq$counts), na.rm = T)/log(freq$counts)}
    )
    for(i in 1:(length(window_weights)-1)){
      weights[which(metadata$time> window_weights[i] &
                      metadata$time<= window_weights[i+1])] = freqs_weighted[i]
    }
  }
  
  ## To store iterations
  best_dev_explained = NULL
  best_AIC = NULL
  best_BIC = NULL
  best_summary = NULL
  best_mod = NULL
  best_groups = NULL
  best_nodes_names = NULL
  
  ## Initialize window
  w = 0 
  window = 0
  
  ## Store split nodes
  nodes_recorded = NULL
  while(window < max(metadata$time)){
    
    ## Set argument for while loop
    go_to_next_window = F
    n_tested_within_window = 0 ## To prevent endless loop
    n_groups = 0
    
    ## Set new window
    window = time_window_initial + w*time_window_increment
    print(paste0('Considering tree up to ', window, '-------------------------->'))
    
    ## Set deviance explained in the new window
    dev_explained_window = 0
    
    while(go_to_next_window == F & n_tested_within_window < 1000 &  n_groups < max_groups_found){
      ## Get the metadata for this time window
      metadata_sub = metadata[which(metadata$time <= window),]
      
      if(nrow(metadata_sub) >= min_group_size*2){
        ## Get weights for isolates in this window 
        weights_sub = weights[which(metadata$time <= window)]
        ## Isolated names
        metadata_sub_rownames = metadata_sub$ID
        ## Get the list of nodes to test for this time window
        nodes_to_test_tmp = nodes_to_test[which(!is.na(match(as.numeric(nodes_to_test), metadata_sub$ID)))]
        
        ## Compute null fit - fit with only the recorded nodes
        mod_null_significant = F
        while(mod_null_significant == F){
          if(length(nodes_recorded) == 0){
            clade_data_null = matrix(NA, nrow = length(metadata_sub_rownames), ncol = 1)
            colnames(clade_data_null) = c('time')
            clade_data_null = as.data.frame(clade_data_null)
            clade_data_null$time = metadata_sub$time
            
            clade_membership_recorded = NULL
            
            ## Fit model
            mod_null = run_gam_model(metadata_sub = metadata_sub, clade_data = clade_data_null, weights_sub = weights_sub, k_smooth = k_smooth, coefs = NULL, sp = NULL, log_y = log_y)
            
            ## Retrieve summary and names of variable
            s_null = summary(mod_null)
            nodes_names_null = 'time'
            
            ## Compute fit error and AIC
            dev_explained_null = s_null$dev.expl
            AIC_null = AIC(mod_null)
            BIC_null = BIC(mod_null)
            
            dev_explained_null_init = s_null$dev.expl
          } else {
            ## First get the clade membership of the nodes found in previous iterations
            clade_membership_recorded = NULL
            for(r in 1:length(nodes_recorded)){
              sub_tmp = descendants_all_nodes[[as.numeric(nodes_recorded[r])]]
              sub_tmp = sub_tmp[which(is.na(match(as.character(sub_tmp), metadata_sub_rownames)) == F)]
              m = match(metadata_sub_rownames, as.character(sub_tmp))
              clade_membership_recorded[[r]] = rep(NA, nrow(metadata_sub))
              clade_membership_recorded[[r]][which(is.na(m) == T)] = 0
              clade_membership_recorded[[r]][which(is.na(m) == F)] = 1
              names(clade_membership_recorded)[r] = nodes_recorded[r]
            }
            nodes_to_test_tmp = nodes_to_test_tmp[which(is.na(match(nodes_to_test_tmp, nodes_recorded)))]
            
            ## Put all data in one dataframe, making sure the columns are ordered by time
            clade_data_null = matrix(NA, nrow = length(metadata_sub_rownames), ncol = length(clade_membership_recorded) + 1)
            colnames(clade_data_null) = c('time', names(clade_membership_recorded))
            clade_data_null = as.data.frame(clade_data_null)
            clade_data_null[,(1:length(clade_membership_recorded))+1] = do.call(cbind, clade_membership_recorded)
            clade_data_null$time = metadata_sub$time
            
            ## Compute groups
            groups = apply(clade_data_null[,2:ncol(clade_data_null)], MARGIN = 1, function(x)paste0(x,collapse = '_'))
            
            # Fit model
            # mod_null = run_gam_model(metadata_sub = metadata_sub, clade_data = clade_data_null, weights_sub = weights_sub, k_smooth = k_smooth, coefs = NULL, sp = NULL)
            
            ## Best model so far is the last one
            mod_null = best_mod[[length(best_mod)]]
            
            ## Retrieve summary and names of variable (in correct order)
            s_null = summary(mod_null)
            nodes_names_null = colnames(clade_data_null)
            
            ## Compute fit error and AIC
            dev_explained_null = s_null$dev.expl
            AIC_null = AIC(mod_null)
            BIC_null = BIC(mod_null)
          }
          
          if(!all(s_null$s.pv <= p_value_smooth)){
            mod_null_significant = F
            if(length(nodes_recorded) == 1){
              nodes_recorded = NULL
              print('Problem: no node recorded')
            }else{
              index = which(s_null$s.pv > p_value_smooth)-1
              index = index[which(index !=0)]
              if(length(index) > 0){
                nodes_recorded = nodes_recorded[-index]
                print('Problem: no node recorded')
              }else{
                mod_null_significant = T
              }
            }
          }else{
            mod_null_significant = T
          }
        }
        
        if(n_tested_within_window == 0){ ## If this is the first iteration within the window, set a null deviance explained
          dev_explained_window = dev_explained_null
        }
        
        print(paste0('Testing ', length(nodes_to_test_tmp), ' nodes'))
        if(parallelize_code == T){
          cl = parallel::makeCluster(number_cores)
          parallel::clusterExport(cl=cl, list("metadata_sub", "metadata_sub_rownames", "nodes_to_test_tmp", 
                                              "descendants_all_nodes", "min_group_size",
                                              "clade_membership_recorded", "clade_data_null", 
                                              "weights_sub", "k_smooth", "mod_null", "test_node_gam", "run_gam_model"),
                                  envir=environment())
          mods_tmp = parallel::parLapplyLB(cl = cl, X = 1:length(nodes_to_test_tmp),
                                           fun = function(X)test_node_gam(X, metadata_sub, metadata_sub_rownames,
                                                                          nodes_to_test_tmp, descendants_all_nodes,
                                                                          min_group_size, clade_membership_recorded,
                                                                          clade_data_null, weights_sub, k_smooth, mod_null, log_y = log_y))
          parallel::stopCluster(cl)
        }else{
          mods_tmp = lapply(X = 1:length(nodes_to_test_tmp), 
                            FUN = function(X)test_node_gam(X, metadata_sub, metadata_sub_rownames, 
                                                           nodes_to_test_tmp, descendants_all_nodes, 
                                                           min_group_size, clade_membership_recorded, 
                                                           clade_data_null, weights_sub, k_smooth, mod_null, log_y = log_y))
        }
        
        ## Set elements to store results
        dev_explained_tmp = rep(NA, length(nodes_to_test_tmp))
        AICs_tmp = rep(NA, length(nodes_to_test_tmp))
        BICs_tmp = rep(NA, length(nodes_to_test_tmp))
        summaries_tmp = as.list(rep(NA, length(nodes_to_test_tmp)))
        groups_tmp = as.list(rep(NA, length(nodes_to_test_tmp)))
        
        for(i in 1:length(mods_tmp)){
          if(!is.na(mods_tmp[[i]]$mod[1])){
            # Fill in results
            s = summary(mods_tmp[[i]]$mod)
            if(all(!is.na(s$s.pv) & !is.na(s$p.pv))){
              if(all(s$s.pv < p_value_smooth, na.rm = T) & all(s$p.pv < p_value_smooth)){
                dev_explained_tmp[i] = s$dev.expl
                AICs_tmp[i] = AIC(mods_tmp[[i]]$mod)
                BICs_tmp[i] = BIC(mods_tmp[[i]]$mod)
                summaries_tmp[[i]] = s
                groups_tmp[[i]] = mods_tmp[[i]]$groups
              }else{
                dev_explained_tmp[i] = NA
                AICs_tmp[i] = NA
                BICs_tmp[i] = NA
                summaries_tmp[[i]] = NA
                groups_tmp[[i]] = NA
              }
            }else{
              dev_explained_tmp[i] = NA
              AICs_tmp[i] = NA
              BICs_tmp[i] = NA
              summaries_tmp[[i]] = NA
              groups_tmp[[i]] = NA
            }
          }else{
            dev_explained_tmp[i] = NA
            AICs_tmp[i] = NA
            BICs_tmp[i] = NA
            summaries_tmp[[i]] = NA
            groups_tmp[[i]] = NA
          }
        }
        ## Find the model that maximizes the explained variance
        if(!all(is.na(dev_explained_tmp))){
          diff_dev_explained = max(dev_explained_tmp, na.rm = T) -  dev_explained_window
          max_dev_explained = max(dev_explained_tmp, na.rm = T)
          BIC_diff = BICs_tmp - BIC_null
          AIC_diff = AICs_tmp - AIC_null
          min_BIC_diff = min(BIC_diff, na.rm = T)
          min_AIC_diff = min(AIC_diff, na.rm = T)
        }else{
          max_dev_explained = -1000
          diff_dev_explained = -1000
          min_BIC_diff = -1000
          min_AIC_diff = -1000
        }
        
        if(min_AIC_diff < -stepwise_AIC_threshold & diff_dev_explained > stepwise_deviance_explained_threshold){
          dev_explained_window = max(dev_explained_tmp, na.rm = T)
          index = which.max(dev_explained_tmp)
          nodes_recorded = c(nodes_recorded, nodes_to_test_tmp[index])
          n_groups = n_groups + 1
          print('Better fit found')
          print(paste0('Deviance explained by ',  round(diff_dev_explained*100, digits = 3), '% more'))
          print(nodes_recorded)
          
          if(plot_screening == T){
            ## Find group colours
            colors_names = MetBrewer::met.brewer(name="Cross", n=length(levels(as.factor(groups_tmp[[index]]))), type="continuous")
            colors = as.factor(groups_tmp[[index]])
            levels(colors) = colors_names
            colors = as.character(colors)
            
            ## Plot
            plot(metadata_sub$time, metadata_sub$index, pch = 16, bty = 'n', col = colors)
            if(log_y == T){
              points(metadata_sub$time, exp(predict(mods_tmp[[index]]$mod)), pch = 16, cex = 0.5, col = colors)
            }else if(log_y == F){
              points(metadata_sub$time, predict(mods_tmp[[index]]$mod), pch = 16, cex = 0.5, col = colors)
            }
          }
          
          ## Save best result
          idx = length(best_dev_explained)+1
          best_dev_explained = c(best_dev_explained, dev_explained_window)
          best_AIC = c(best_AIC, AICs_tmp[[index]])
          best_BIC = c(best_BIC, BICs_tmp[[index]])
          best_summary[[idx]] = summaries_tmp[[index]]
          best_mod[[idx]] = mods_tmp[[index]]$mod
          best_groups[[idx]] = groups_tmp[[index]]
          best_nodes_names[[idx]] = nodes_recorded
        }else{
          go_to_next_window = T
          print('No better fit found')
        }
        n_tested_within_window = n_tested_within_window + 1
      }else{
        go_to_next_window = T
      }
    }
    w = w +1
  }
  if(keep_track == F){
    return('potential_splits' =  as.numeric(nodes_recorded))
  }else{
    return(list('potential_splits' =  as.numeric(nodes_recorded),
                'best_dev_explained' =  as.numeric(best_dev_explained),
                'first_dev' = dev_explained_null_init,
                'best_AIC' = best_AIC,
                'best_BIC'= best_BIC,
                'best_summary' = best_summary,
                'best_mod' = best_mod,
                'best_groups' = best_groups,
                'best_nodes_names' = best_nodes_names))
  }
}

#' Sub function to ouput a fit, for a selected number of nodes
compute_fit_for_node_set = function(timed_tree, 
                                    metadata, 
                                    node_set,
                                    weight_by_time = NULL,
                                    weighting_transformation = c(NULL, 'inv_freq', 'inv_sqrt', 'inv_log'),
                                    k_smooth = 5,
                                    plot = T, 
                                    log_y = T){
  ## Make sure split nodes are character strings
  node_set = as.character(node_set)
  
  ## Store paths from root to all tips and nodes
  paths_root_to_all_nodes = path.root.to.all.nodes(timed_tree)
  
  ## Store descendants for all nodes
  descendants_all_nodes = descendants.all.nodes(timed_tree)
  
  ## Store number of descendants for all nodes
  nb_descendants_all_nodes = number.descendants.all.nodes(timed_tree)
  
  ## Remove potential NAs
  idx = which(is.infinite(metadata$index) | is.nan(metadata$index))
  if(length(idx)>0) metadata = metadata[-idx,]
  
  ## Remove potential 1s
  idx = which(metadata$index >= 1)
  if(length(idx)>0) metadata$index[idx] = 0.99
  
  ## Row names
  metadata_rownames = metadata$ID
  
  ## Set weights, if user selected weighting by time
  weights = NULL
  if(!is.numeric(weight_by_time) & !is.null(weight_by_time)) stop('Weight_by_time is not NULL or numeric')
  if(!is.null(weight_by_time)){
    weights = rep(NA, length(metadata$time))
    window_weights = seq(round(min(metadata$time), digits = 1)-1e4, round(max(metadata$time), digits = 1)+1e4, weight_by_time)
    freq = hist(metadata$time, breaks = window_weights, plot = F, include.lowest = T)
    freq$counts[which(freq$counts == 0)] = NA
    if(length(weighting_transformation) != 1) stop('Weighting_transformation parameter is not of length 1, please supply one weighting transformation')
    freqs_weighted = switch(weighting_transformation, 
                            'inv_freq' = {sum(freq$counts, na.rm = T)/freq$counts}, 
                            'inv_sqrt' = {sum(sqrt(freq$counts), na.rm = T)/sqrt(freq$counts)}, 
                            'inv_log'  = {sum(log(freq$counts), na.rm = T)/log(freq$counts)}
    )
    for(i in 1:(length(window_weights)-1)){
      weights[which(metadata$time> window_weights[i] &
                      metadata$time<= window_weights[i+1])] = freqs_weighted[i]
    }
  }
  
  ## Compute final, best fit
  if(length(node_set) == 0){
    clade_data_best = matrix(NA, nrow = length(metadata_rownames), ncol = 1)
    colnames(clade_data_best) = c('time')
    clade_data_best = as.data.frame(clade_data_best)
    clade_data_best$time = metadata$time
    
    ## Fit model
    mod_best = run_gam_model(metadata_sub = metadata, clade_data = clade_data_best, weights_sub = weights, k_smooth = k_smooth, sp = NULL, log_y = log_y)
    
    ## Retrieve summary and names of variable
    s_best = summary(mod_best)
    nodes_names_best = 'time'
    
    ## Compute fit error and AIC
    dev_explained_best = s_best$dev.expl
    AIC_best = AIC(mod_best)
    BIC_best = BIC(mod_best)
  } else {
    ## First get the clade membership of the nodes found in previous iterations
    clade_membership_recorded = NULL
    for(r in 1:length(node_set)){
      sub_tmp = descendants_all_nodes[[as.numeric(node_set[r])]]
      sub_tmp = sub_tmp[which(is.na(match(as.character(sub_tmp), metadata_rownames)) == F)]
      m = match(metadata_rownames, as.character(sub_tmp))
      clade_membership_recorded[[r]] = rep(NA, nrow(metadata))
      clade_membership_recorded[[r]][which(is.na(m) == T)] = 0
      clade_membership_recorded[[r]][which(is.na(m) == F)] = 1
      names(clade_membership_recorded)[r] = node_set[r]
    }
    
    ## Put all data in one dataframe, making sure the columns are ordered by time
    clade_data_best = matrix(NA, nrow = length(metadata_rownames), ncol = length(clade_membership_recorded) + 1)
    colnames(clade_data_best) = c('time', names(clade_membership_recorded))
    clade_data_best = as.data.frame(clade_data_best)
    clade_data_best[,(1:length(clade_membership_recorded))+1] = do.call(cbind, clade_membership_recorded)
    clade_data_best$time = metadata$time
    
    ## Compute groups
    groups_best = apply(clade_data_best[,2:ncol(clade_data_best)], MARGIN = 1, function(x)paste0(x,collapse = '_'))
    
    # Fit model
    mod_best = run_gam_model(metadata_sub = metadata, clade_data = clade_data_best, weights_sub = weights, k_smooth = k_smooth, sp = NULL, log_y = log_y)
    
    ## Retrieve summary and names of variable (in correct order)
    s_best = summary(mod_best)
    nodes_names_best = colnames(clade_data_best)
    
    ## Compute fit error and AIC
    dev_explained_best = s_best$dev.expl
    AIC_best = AIC(mod_best)
    BIC_best = BIC(mod_best)
  }
  colors_names = MetBrewer::met.brewer(name="Cross", n=length(levels(as.factor(groups_best))), type="continuous")
  colors_names = sample(colors_names, replace = F)
  colors = as.factor(groups_best)
  levels(colors) = colors_names
  colors = as.character(colors)
  
  if(plot == T){
    plot(metadata$time, metadata$index, pch = 1, bty = 'n', col = adjustcolor(colors, alpha.f = 0.5) , cex = 0.5, 
         ylim = c(0, 1),
         yaxt = 'n',
         xlab = 'Time',
         ylab = 'Index')
    axis(2, las = 2)
    if(log_y == T){
      points(metadata$time, exp(predict(mod_best)), pch = 16, cex = 0.5, col = colors)
    }else if(log_y == F){
      points(metadata$time, predict(mod_best), pch = 16, cex = 0.5, col = colors)
    }
    
  }
  
  return(list('mod' = mod_best,
              'summary' = s_best,
              'weights' = weights,
              'deviance_explained' = dev_explained_best,
              'AIC' = AIC_best,
              'groups' = groups_best,
              'colors' = colors))
}



