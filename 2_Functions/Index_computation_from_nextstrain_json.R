#!/usr/bin/env Rscript

# Install packages if needed

library(argparse)
library(stringr)
library(ape)
library(treeio)
library(lubridate)
library(MetBrewer)
library(extrafont)

########################################################################################################################################
## Argument parser
########################################################################################################################################
parser <- ArgumentParser(description = "Compute index dynamics from command line")

parser$add_argument(
  "--tree_file",
  type = "character",
  default = "",
  help = "Path to the tree file. Has to be a json tree from a Nextstrain build."
)

parser$add_argument(
  "--mu",
  type = "double",
  default = 2.5e-7,
  help = "Mutation rate (mutations/site/year)"
)

parser$add_argument(
  "--genome_length",
  type = "integer",
  default = 4086189,
  help = "Genome length in base pairs"
)

parser$add_argument(
  "--timescale",
  type = "double",
  default = 2,
  help = "Timescale parameter for index"
)

parser$add_argument(
  "--wind",
  type = "double",
  default = 1,
  help = "Window length (years)"
)

parser$add_argument(
  "--plot_res",
  type = "logical",
  default = F,
  help = "Plot results to PDF (TRUE/FALSE)"
)

args <- parser$parse_args()
########################################################################################################################################

########################################################################################################################################
## Use the arguments
########################################################################################################################################
cat("Tree file:", args$tree_file, "\n")
cat("Mutation rate:", args$mu, "\n")
cat("Genome length:", args$genome_length, "\n")
cat("Timescale:", args$timescale, "\n")
cat("Window:", args$wind, "\n")
cat("Plot results:", args$plot_res, "\n")

tree_file = args$tree_file
mu = args$mu
genome_length = args$genome_length
timescale = args$timescale
wind = args$wind
plot_res = args$plot_res
########################################################################################################################################

########################################################################################################################################
## Useful functions
########################################################################################################################################
## Load index functions
compute.index = function(time_distance_mat, timed_tree, time_window, metadata, mutation_rate, timescale, genome_length){
  ## TO DO: ADD CHECKS
  
  ## Matrices computation
  matrices = matrices.computation(time_distance_mat, timed_tree, time_window, metadata, mutation_rate, genome_length)
  
  ## Index computation
  # Compute index bandwidth
  bandwidth = index.bandwidth(timescale, genome_length, mutation_rate)
  # Check zero diagonal
  stopifnot(all(diag(matrices$hamming_mat) == 0))
  # Check window matrix has the same dimensions as H
  stopifnot(all(dim(matrices$hamming_mat) == dim(matrices$window_mat)))
  # Let B <- b^h
  B <- bandwidth^matrices$hamming_mat
  # Take columns sums, of the corrected matrix,  subtract 1 (= b^0) to ignore diagonal
  Bsums <- colSums(B, na.rm = T)-1
  
  return(Bsums/(colSums(matrices$window_mat)-1))
}
index.bandwidth = function(timescale, genome_length, mutation_rate) {
  if (timescale == 0) return(NA)
  h = genome_length*(1 - exp(-2 * mutation_rate * timescale))
  objfun = function(b) (1/2 * (1 - b^genome_length) - (1 - b^h))^2
  opt = optimise(objfun, c(1e-09, 1 - 1e-09))
  return(opt$minimum)
}
matrices.computation = function(time_distance_mat, timed_tree, time_window, metadata, mutation_rate, genome_length){
  ## Hamming distance matrix computation by window
  ## Retrieve sequence names from tree
  names_seqs = timed_tree$tip.label
  
  ## Construct a new distance matrix, that compare tips and nodes to the population circulating within a window of time (including branches, i.e. unsampled individuals in the population)
  timed_distance_mat_by_window = time_distance_mat
  
  ## Compute branches birth and death
  branches_times = timed_tree$edge
  branches_times[,1] = branches_times[,2] = NA
  times_tips_nodes = metadata$time
  branches_times[,1] = times_tips_nodes[match(timed_tree$edge[,1], metadata$ID)]
  branches_times[,2] = times_tips_nodes[match(timed_tree$edge[,2], metadata$ID)]
  
  ## Loop through all tips and nodes to correct the matrix for each individual
  for(sample in colnames(time_distance_mat)){
    ## Metadata
    metadata_tmp = metadata[which(metadata$name_seq == sample),]
    
    ## Distances between these chosen tips and nodes
    i = which(colnames(time_distance_mat) == sample)
    time_distance_mat_tmp = time_distance_mat[,i] 
    
    ## Filter branches alive at sampling time
    t_min = metadata_tmp$time - time_window
    t_max = metadata_tmp$time + time_window
    tmp = which((branches_times[,1] <= t_min & branches_times[,2] >= t_max) | ## branches alive, with no sampled individuals 
                  (branches_times[,1] <= t_min & (branches_times[,2] >= t_min & branches_times[,2] <= t_max)) | ## branches born before interval and died within interval
                  ((branches_times[,1] >= t_min & branches_times[,1] <= t_max) & branches_times[,2] >= t_max) | ## branches born in interval and died after interval
                  ((branches_times[,1] >= t_min & branches_times[,1] <= t_max) & (branches_times[,2] >= t_min & branches_times[,2] <= t_max))) ## branches born and died within interval
    
    ## For each sampling above, get distances
    edge_lengths = branches_times[tmp,]
    edges = timed_tree$edge[tmp,]
    
    if(is.null(dim(edges))){
      ## List all nodes in the tables
      all_nodes_and_tips = sort(c(edges[2])) ## Choice: consider one individual per branch, and the lastest sample (offspring)
      time_distance_mat_tmp[is.na(match(metadata$ID, unique(all_nodes_and_tips)))] = NA ## Only keep the relevant nodes
      
      ## Self distance should be 0
      time_distance_mat_tmp[i] = 0
    }
    if(!is.null(dim(edges))){
      ## List all nodes in the tables
      all_nodes_and_tips = sort(c(edges[,2])) ## Choice: consider one individual per branch, and the lastest sample (offspring)
      time_distance_mat_tmp[is.na(match(metadata$ID, unique(all_nodes_and_tips)))] = NA ## Only keep the relevant nodes
      
      ## Remove connecting nodes: those are messing up the distance matrix
      tbl = table(c(edges[,1], edges[,2])) 
      ## Any node that is present 3 times should be removed: it's a connecting node - which is screwing up the distance matrix later
      ## Any node that is present 2 times is a mrca, but is already not taken into account is the computation, so it's fine
      connecting_nodes = as.numeric(names(tbl)[which(tbl == 3)])
      idx = !is.na(match(metadata$ID, connecting_nodes))
      time_distance_mat_tmp[idx] = NA ## Remove connecting nodes
      
      if(length(which(connecting_nodes == i)) > 0){ ## The node i is a connecting node
        idx = match(metadata$ID, edges[which(edges[,1] == i),2])
        idx = which(!is.na(idx))
        time_distance_mat_tmp[idx] = NA ## Remove tips descending from connecting nodes
      }
      
      ## Correct this vector for the nodes that have a sampling time >metadata_tmp$time
      too_late = which(edge_lengths[,2] > metadata_tmp$time)
      m = match(edges[too_late,2], metadata$ID)
      time_distance_mat_tmp[m] = time_distance_mat_tmp[m]  - abs(edge_lengths[too_late,2] - metadata_tmp$time)
      
      ## Correct this vector for the nodes that have a sampling time <metadata_tmp$time
      too_early = which(edge_lengths[,2] < metadata_tmp$time)
      m = match(edges[too_early,2], metadata$ID)
      time_distance_mat_tmp[m] = time_distance_mat_tmp[m]  + abs(edge_lengths[too_early,2] - metadata_tmp$time)
      
      ## Self distance should be 0
      time_distance_mat_tmp[i] = 0
    }
    ## Put this vector in the big matrix 
    timed_distance_mat_by_window[,i] = time_distance_mat_tmp
  }
  
  ## Isolates within the same range
  same_year_mat = !is.na(timed_distance_mat_by_window)
  
  return(list('hamming_mat' = timed_distance_mat_by_window*mutation_rate*genome_length/2, 
              'window_mat' = same_year_mat))
}
dist.nodes.with.names = function(timed_tree){
  # Check tree is binary
  stopifnot(is.binary.phylo(timed_tree))
  names_seqs = timed_tree$tip.label
  genetic_distance_mat = ape::dist.nodes(timed_tree) ## Pairwise time of differences
  colnames(genetic_distance_mat) = c(names_seqs, length(names_seqs)+(1:(length(names_seqs)-1)))
  rownames(genetic_distance_mat) = c(names_seqs, length(names_seqs)+(1:(length(names_seqs)-1)))
  return(genetic_distance_mat)
}
axisPhylo_NL = function (side = 1, root.time = NULL, backward = TRUE, at_axis = NULL, lab_axis = NULL, ...){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  type <- lastPP$type
  if (type == "unrooted")
    stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
  if (type == "radial")
    stop("axisPhylo() not meaningful for this type of plot")
  if (is.null(root.time))
    root.time <- lastPP$root.time
  if (type %in% c("phylogram", "cladogram")) {
    xscale <- if (lastPP$direction %in% c("rightwards", "leftwards"))
      range(lastPP$xx)
    else range(lastPP$yy)
    tmp <- lastPP$direction %in% c("leftwards", "downwards")
    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp))
      tscale <- tscale[2:1]
    if (!is.null(root.time)) {
      tscale <- tscale + root.time
      if (backward)
        tscale <- tscale - xscale[2]
    }
    beta <- diff(xscale)/diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]
    if(is.null(at_axis) == T){
      x <- beta * lab + alpha
      lab <- pretty(tscale)
    }
    # if(is.null(at_axis) != F){
    x <- at_axis
    lab <- lab_axis
    # }
    axis(side = side, at = x, labels = lab, ...)
  }
  else {
    n <- lastPP$Ntip
    xx <- lastPP$xx[1:n]
    yy <- lastPP$yy[1:n]
    r0 <- max(sqrt(xx^2 + yy^2))
    alpha <- sort(setNames(rect2polar(xx, yy)$angle, 1:n))
    angles <- c(diff(alpha), 2 * pi - alpha[n] + alpha[1L])
    j <- which.max(angles)
    i <- if (j == 1L)
      n
    else j - 1L
    firstandlast <- as.integer(names(angles[c(i, j)]))
    theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
    x0 <- r0 * cos(theta0)
    y0 <- r0 * sin(theta0)
    inc <- diff(pretty(c(0, r0))[1:2])
    srt <- 360 * theta0/(2 * pi)
    coef <- -1
    if (abs(srt) > 90) {
      srt <- srt + 180
      coef <- 1
    }
    len <- 0.025 * r0
    r <- r0
    while (r > 1e-08) {
      x <- r * cos(theta0)
      y <- r * sin(theta0)
      if (len/r < 1) {
        ra <- sqrt(len^2 + r^2)
        thetaa <- theta0 + coef * asin(len/r)
        xa <- ra * cos(thetaa)
        ya <- ra * sin(thetaa)
        segments(xa, ya, x, y)
        text(xa, ya, r0 - r, srt = srt, adj = c(0.5,
                                                1.1), ...)
      }
      r <- r - inc
    }
    segments(x, y, x0, y0)
  }
}
########################################################################################################################################

########################################################################################################################################
## Read tree
########################################################################################################################################
tree_json = read.nextstrain.json(tree_file)

tree = tree_json@phylo
metadata = as.data.frame(tree_json@data)

## Rescale the tree to make sure it is on the same scale as metadata
tree_edge = tree$edge
tree_edge_times = tree_edge
for(i in 1:length(tree$tip.label)){
  idx = which(tree_edge == i, arr.ind = T)
  idx2 = which(metadata$node == i, arr.ind = T)
  tree_edge_times[idx] = metadata$num_date[idx2]
}
for(i in 1:(length(tree$tip.label)-1)){
  idx = which(tree_edge == i + length(tree$tip.label), arr.ind = T)
  idx2 = which(metadata$node == i + length(tree$tip.label), arr.ind = T)
  tree_edge_times[idx] = metadata$num_date[idx2]
}
tree$edge.length = tree_edge_times[,2]-tree_edge_times[,1]
########################################################################################################################################

########################################################################################################################################
## Prepare data
########################################################################################################################################
tree$tip.label = metadata$accession[match(1:length(tree$tip.label), metadata$node)]
names_seqs = tree$tip.label

## Create dataset with names of each sequence, time sampling, lineage
dataset_tips = data.frame('ID' = 1:length(names_seqs),
                          'name_seq' = tree$tip.label,
                          'time' = metadata$num_date[match(1:length(tree$tip.label), metadata$node)],
                          'clade_membership' = metadata$clade_membership[match(1:length(tree$tip.label), metadata$node)])
dataset_tips$time = as.numeric(dataset_tips$time)
########################################################################################################################################

########################################################################################################################################
## Preparation data nodes
########################################################################################################################################
## Compute distance between each pair of sequences AND NODES in the tree
genetic_distance_mat = dist.nodes.with.names(tree)

## Get the time each node
nroot = length(tree$tip.label)+1 ## Checked and it's the root 
distance_to_root = genetic_distance_mat[nroot,]
root_height = dataset_tips$time[which(dataset_tips$name_seq == names(distance_to_root[20]))] - distance_to_root[20]  ## Take one tip, doesn't matter which tip is used
nodes_height = root_height + distance_to_root[length(names_seqs)+(1:(length(names_seqs)-1))]

# Meta-data with nodes 
dataset_with_nodes = data.frame('ID' = c(1:length(names_seqs), length(names_seqs)+(1:(length(names_seqs)-1))),
                               'name_seq' = c(names_seqs, length(names_seqs)+(1:(length(names_seqs)-1))),
                               'time' = c(dataset_tips$time, nodes_height),
                               'is.node' = c(rep('no', length(names_seqs)), rep('yes', (length(names_seqs)-1))),
                               'clade_membership' = metadata$clade_membership[match(c(1:length(tree$tip.label),length(names_seqs)+(1:(length(names_seqs)-1))), metadata$node)]) 
########################################################################################################################################

########################################################################################################################################
## Compute index of every tip and node
########################################################################################################################################
dataset_with_nodes$index = compute.index(time_distance_mat = genetic_distance_mat, 
                                        timed_tree = tree, 
                                        time_window = wind,
                                        metadata = dataset_with_nodes, 
                                        mutation_rate = mu,
                                        timescale = timescale,
                                        genome_length = genome_length)
########################################################################################################################################

########################################################################################################################################
## Generate colors needed
########################################################################################################################################
lev = as.character(unique(dataset_with_nodes$clade_membership))
n = length(lev)
colors_genotypes = c(met.brewer(name="Cross", n=n, type="continuous"))

dataset_with_nodes$clade_membership_color = dataset_with_nodes$clade_membership
dataset_with_nodes$clade_membership_color = as.factor(dataset_with_nodes$clade_membership_color)
labels = levels(dataset_with_nodes$clade_membership_color)
levels(dataset_with_nodes$clade_membership_color) = colors_genotypes[match(labels, lev)]
dataset_with_nodes$clade_membership_color = as.character(dataset_with_nodes$clade_membership_color)
dataset_with_nodes$clade_membership_color[which(is.na(dataset_with_nodes$clade_membership_color))] = 'grey60'
########################################################################################################################################

########################################################################################################################################
## Plot result
########################################################################################################################################
if(plot_res == T){
  pdf(paste0('Index_dynamics_nextstran_json_', Sys.Date(), '.pdf'), width = 6/2.54, height = 5/2.54)
  par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
      cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)
  
  min_year = min(dataset_with_nodes$time)
  max_year = max(dataset_with_nodes$time)
  max_index = max(dataset_with_nodes$index)
  
  ## Plot tree
  tree_to_plot = tree
  plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.2, edge.color = 'grey', 
       x.lim = c(min_year, max_year)-root_height)
  nodelabels(pch = 16, col = dataset_with_nodes$clade_membership_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.2)
  tiplabels(pch = 16, col = dataset_with_nodes$clade_membership_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.2)
  axisPhylo_NL(side = 1, root.time = root_height, backward = F,
               at_axis = seq(round(min_year), round(max_year), 10)-root_height,
               lab_axis = seq(round(min_year), round(max_year), 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
  legend('topleft', title = 'Colours', col = colors_genotypes, legend = lev, pch = 16, cex = 0.2, box.col = 'white', bg = 'white')
  
  ## Plot index
  plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
       col = adjustcolor(dataset_with_nodes$clade_membership_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
       bty = 'n', xlim = c(min_year, max_year), cex = 0.2, 
       pch = 16, bty = 'n', ylim = c(0, max_index), 
       ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
  points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
         dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
         col = adjustcolor(dataset_with_nodes$clade_membership_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
         cex = 0.2, pch = 16)
  axis(1, at = seq(round(min_year), round(max_year), 10), labels = seq(round(min_year), round(max_year), 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
  axis(2, las = 2, tck=-0.01, lwd = 0.5)
  title(main="Index dynamics", line=-0.25, outer = F)
  title(ylab="Index", line=0.5, outer = F)
  title(xlab="Time (years)", line=0, outer = F)
  dev.off()
}
########################################################################################################################################

########################################################################################################################################
## Save csv with results
########################################################################################################################################
write.csv(dataset_with_nodes, paste0('Index_dynamics_nextstran_json_', Sys.Date(), '.csv'))
########################################################################################################################################

########################################################################################################################################
## Print done
########################################################################################################################################
cat('Done')
########################################################################################################################################