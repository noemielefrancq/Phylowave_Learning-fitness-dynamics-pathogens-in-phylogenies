#' @title Compute the diversity index from pairwise evolutionary time differences between all nodes and tips
#' @param time_distance_mat Hamming distance matrix
#' @param timed_tree Timed tree
#' @param time_window Time window for the index computation
#' @param metadata Metadata for all tips and nodes (see details)
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param timescale Index timescale (in years)
#' @param genome_length Genome length (in bp)
#' @return The index of each sequence and node in the tree
#' @references Wirth, T, Wong, V, Vandenesch, F, Rasigade, J-P. Applied phyloepidemiology: Detecting drivers of pathogen transmission from genomic signatures using density measures. Evol Appl. 2020; 13: 1513– 1525. https://doi.org/10.1111/eva.12991
#' @export
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

#' @title Compute the index bandwidth from a set timescale, genome length and mutation rate 
#' @param timescale Index timescale (in years)
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param genome_length Genome length (in bp)
#' @return Index bandwidth
#' @references Wirth, T, Wong, V, Vandenesch, F, Rasigade, J-P. Applied phyloepidemiology: Detecting drivers of pathogen transmission from genomic signatures using density measures. Evol Appl. 2020; 13: 1513– 1525. https://doi.org/10.1111/eva.12991
#' @export
index.bandwidth = function(timescale, genome_length, mutation_rate) {
  if (timescale == 0) return(NA)
  h = genome_length*(1 - exp(-2 * mutation_rate * timescale))
  objfun = function(b) (1/2 * (1 - b^genome_length) - (1 - b^h))^2
  opt = optimise(objfun, c(1e-09, 1 - 1e-09))
  return(opt$minimum)
}

#' @title Computes a matrix of pairwise distances for all pairs of isolates sampled at within the same window of time
#' @param time_distance_mat Hamming distance matrix
#' @param timed_tree Timed tree
#' @param time_window Time window for the index computation
#' @param metadata Metadata for all tips and nodes (see details)
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param genome_length Genome length (in bp)
#' @export
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

#' @title Compute distance between all tips and nodes in the timed tree
#' @import ape
#' @param timed_tree Timed tree
#' @return Distance matrix
#' @export
dist.nodes.with.names = function(timed_tree){
  # Check tree is binary
  stopifnot(is.binary.phylo(timed_tree))
  names_seqs = timed_tree$tip.label
  genetic_distance_mat = ape::dist.nodes(timed_tree) ## Pairwise time of differences
  colnames(genetic_distance_mat) = c(names_seqs, length(names_seqs)+(1:(length(names_seqs)-1)))
  rownames(genetic_distance_mat) = c(names_seqs, length(names_seqs)+(1:(length(names_seqs)-1)))
  return(genetic_distance_mat)
}


#' @title Compute best timescale
#' @importFrom stats lm
#' @param genome_length Genome length (in bp)
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param time_window Time window for gathering sequences in groups
#' @param Ne Effective population size
#' @param r_squared_threshold Threshold to find the best timescale
#' @export
compute.timescale=function(genome_length, mutation_rate, time_window, Ne=1E5, r_squared_threshold=0.99){
  t_seq =seq(0.1,2,0.1) ## test for different values of timescale  
  b_seq=c()
  for (t in t_seq){#compute the bandwidth for each  timescale
    b <- index.bandwidth(timescale = t, genome_length, mutation_rate) 
    b_seq=c(b_seq,b)
  }
  times = seq(0,20,0.1)
  wind=c()
  for (k in 1:length(t_seq)){#for each timescale
    b=b_seq[k]
    index_constant_pop = theoretical_index(b, genome_length, mu, Ne, times)
    len_wind=5
    linear_model=lm(log(index_constant_pop[1:len_wind])~times[1:len_wind])
    R_2=summary(linear_model)$r.squared
    while (R_2>r_squared_threshold & len_wind<length(times)){
      len_wind=len_wind+1#increase the window as long as it is not linear enough
      linear_model=lm(log(index_constant_pop[1:len_wind])~times[1:len_wind])
      R_2=summary(linear_model)$r.squared
    }
    if (len_wind==length(times)){#we reached a threshold
      break
    }
    wind=c(wind,times[len_wind])
  }
  slope=lm(wind~t_seq[1:length(wind)])$coefficients[2]
  timescale=time_window/slope
  return(timescale)
}


#' @title Compute the theoretical index for a population of constant size Ne during times
#' @param bandwidth Index timescale (in years)
#' @param genome_length Genome length (in bp)
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param Ne Effective population size
#' @param time Time at which to compute the theoretical index (can be a double or vector)
#' @export
theoretical_index = function(bandwidth, genome_length, mutation_rate, Ne, times) {
  return((1-(bandwidth*exp(-1/(genome_length*mutation_rate*Ne)))^(genome_length*mutation_rate*times +1)) /  (1-(bandwidth*exp(-1/(genome_length*mutation_rate*Ne)))) / (  (1-(exp(-1/(genome_length*mutation_rate*Ne)))^(genome_length*mutation_rate*times +1)) /  (1-(exp(-1/(genome_length*mutation_rate*Ne))))))
}


#' @title Equivalent of axisPhylo() in ape, with the option of relabeling the axis
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


#' Compute the diversity index from a distance matrix (no tree required)
#'
#' Computes the diversity index directly from a pairwise distance matrix between
#' samples, without requiring a phylogenetic tree. This is useful when users have
#' precomputed genetic distances (e.g., from sequence alignments) without a tree.
#'
#' @param distance_mat Pairwise distance matrix (e.g., Hamming distances) between samples
#' @param time_window Time window for the index computation (positive number, in years)
#' @param metadata Data frame with columns \code{ID} and \code{time}
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param timescale Index timescale (in years)
#' @param genome_length Genome length (in bp)
#' @return Named vector of index values for each sample
#' @export
compute.index.from.distance = function(distance_mat, time_window, metadata, mutation_rate, timescale, genome_length){
  ## Input validation
  if (!is.matrix(distance_mat)) {
    stop("'distance_mat' must be a matrix.")
  }
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data.frame.")
  }
  required_cols <- c("ID", "time")
  missing_cols <- setdiff(required_cols, names(metadata))
  if (length(missing_cols) > 0) {
    stop("'metadata' is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(time_window) || length(time_window) != 1 || time_window <= 0) {
    stop("'time_window' must be a single positive number.")
  }
  if (!is.numeric(mutation_rate) || length(mutation_rate) != 1 || mutation_rate <= 0) {
    stop("'mutation_rate' must be a single positive number.")
  }
  if (!is.numeric(timescale) || length(timescale) != 1 || timescale <= 0) {
    stop("'timescale' must be a single positive number.")
  }
  if (!is.numeric(genome_length) || length(genome_length) != 1 || genome_length <= 0) {
    stop("'genome_length' must be a single positive integer.")
  }
  
  ## Match metadata to matrix
  sample_ids <- colnames(distance_mat)
  times <- metadata$time[match(sample_ids, metadata$ID)]
  if (any(is.na(times))) {
    stop("Some sample IDs in 'distance_mat' are not found in 'metadata$ID'.")
  }
  n <- length(sample_ids)
  
  ## Time-window filtering and temporal correction
  hamming_mat <- distance_mat
  window_mat <- matrix(TRUE, nrow = n, ncol = n)
  
  for (j in seq_len(n)) {
    t_j <- times[j]
    for (i in seq_len(n)) {
      if (i == j) next
      t_i <- times[i]
      ## Filter: only keep pairs within the time window
      if (abs(t_i - t_j) > time_window) {
        hamming_mat[i, j] <- NA
        window_mat[i, j] <- FALSE
      } else {
        ## Temporal correction: adjust distance by time difference
        hamming_mat[i, j] <- hamming_mat[i, j] + abs(t_i - t_j)
      }
    }
  }
  
  ## Scale distances
  hamming_mat <- hamming_mat * mutation_rate * genome_length / 2
  
  ## Index computation
  bandwidth <- index.bandwidth(timescale, genome_length, mutation_rate)
  B <- bandwidth^hamming_mat
  Bsums <- colSums(B, na.rm = TRUE) - 1
  index <- Bsums / (colSums(window_mat) - 1)
  names(index) <- sample_ids
  
  return(index)
}


