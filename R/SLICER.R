#' @importFrom igraph V shortest_paths distances graph_from_adjacency_matrix clusters count_components get.shortest.paths
#' @importFrom lle lle
#' @importFrom alphahull ahull areaahull
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics legend lines plot points segments
#' @importFrom stats dist rnorm var
NULL

adaptive_knn_graph = function(traj_dist, k)
{
  adj_mat = matrix(0, nrow = nrow(traj_dist), ncol = ncol(traj_dist))
  knn = t(apply(traj_dist, 1, order))
  for (i in 1:nrow(traj_dist))
  {
    adj_mat[i, knn[i, 2:(k[i] + 1)]] = traj_dist[i, knn[i, 2:(k[i] + 1)]]
  }
  return (adj_mat)
}

#' Helper function for k selection
#' @param k Nearest neighbors
#' @param traj_exp Cell expression matrix
#' @return Width of LLE embedding
#' @export

width_k = function(k, traj_exp)
{
  n = nrow(traj_exp)
  traj_lle = lle(traj_exp, m = 2, k = k)$Y
  traj_graph = conn_knn_graph(traj_lle, k)
  dists = distances(traj_graph)
  len = max(dists)
  alpha = len / 10
  alpha_hull = ahull(traj_lle, alpha = alpha)
  return(areaahull(alpha_hull) / len)
}

#' Select the number of nearest neighbors for LLE to use
#'
#' \code{select_k} uses the alpha-hull to determine which value
#' of k yields an embedding that most resembles a trajectory.
#'
#' @param exp_mat Matrix of expression levels
#' @param kmin Smallest value of k to try
#' @param kmax Largest value of k to try
#' @param by Increment
#' @return The optimal value of k
#' @export
#' @examples
#' \dontrun{
#' genes = select_genes(traj)
#' k = select_k(traj[,genes])
#' }
select_k = function(exp_mat,
                    kmin = 5,
                    kmax = 50,
                    by = 5)
{
  min_k = kmin
  min_width = Inf
  w = Inf
  for (k in seq(kmin, kmax, by))
  {
    tryCatch({
      w = width_k(k, exp_mat)
    },
    error = function(cond)
    {
      w = Inf
    })
    if (is.na(w))
    {
      w = Inf
    }
    if (w < min_width)
    {
      min_k = k
      min_width = w
    }
  }
  return (min_k)
}

dev_ij = function(i, j, traj_exp, adj_mat)
{
  return(sum((traj_exp[i, j] - traj_exp[which(adj_mat[i, ] > 0), j]) ^ 2))
}

selection_val = function(j, traj_exp, adj_mat)
{
  n = nrow(traj_exp)
  dev = sapply(1:n, dev_ij, j, traj_exp, adj_mat)
  k = sum(adj_mat[1, ] > 0)
  return (var(traj_exp[, j]) / (sum(dev) / (n * k - 1)))
}

min_conn_k = function(traj_exp)
{
  traj_dist = as.matrix(dist(traj_exp))
  conn_comp = 2
  k = 0
  while (conn_comp > 1)
  {
    k = k + 1
    adj_mat = adaptive_knn_graph(traj_dist, rep(k, nrow(traj_exp)))
    traj_graph = graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted =
                                               TRUE)
    conn_comp = count_components(traj_graph)
    
  }
  return (k)
}

#' Identify clusters corresponding to putative cell types
#'
#' \code{detect_cell_types} divides the k-nearest neighbor
#' graph (built from the LLE embedding) into connected
#' components. These connected components represent clusters of cells corresponding to putative cell types.
#'
#' @param embedding Low-dimensional LLE embedding of cells
#' @param k Number of nearest neighbors to use when detecting clusters
#' @return Vector containing a numerical cluster assignment for each cell
#' @export

detect_cell_types = function(embedding, k)
{
  traj_dist = as.matrix(dist(embedding))
  adj_mat = adaptive_knn_graph(traj_dist, rep(k, nrow(embedding)))
  traj_graph = graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted =
                                             TRUE)
  cluster_assignments = clusters(traj_graph)$membership
  plot(embedding[, 1], embedding[, 2], pch = 16, col = cluster_assignments)
  return(cluster_assignments)
}

#' Construct a k-nearest neighbor graph that is fully connected
#'
#' This function constructs a k-nearest neighbor graph using an
#' LLE embedding, then adds the minimum number of edges needed to
#' make the graph fully connected.
#'
#' @param embedding Low-dimensional LLE embedding of cells
#' @param k Number of nearest neighbors
#' @return An igraph object corresponding to the k-NN graph
#' @export
#' @examples
#' genes=1:200
#' cells=sample(1:500,30)
#' k=10
#' traj_lle = lle::lle(traj[cells,genes],m=2,k)$Y
#' traj_graph = conn_knn_graph(traj_lle,5)

conn_knn_graph = function(embedding, k)
{
  if (sum(duplicated(embedding)) > 0)
  {
    dup = which(duplicated(embedding))
    for (i in dup)
    {
      embedding[i, ] = embedding[i, ] + rnorm(ncol(embedding), 0, 0.00001)
    }
  }
  traj_dist = as.matrix(dist(embedding))
  adj_mat = adaptive_knn_graph(traj_dist, rep(k, nrow(embedding)))
  traj_graph = graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted =
                                             TRUE)
  conn_comp = count_components(traj_graph)
  n = nrow(embedding)
  while (conn_comp > 1)
  {
    cluster_assignments = clusters(traj_graph)$membership
    min_dist = Inf
    min_i = 0
    min_j = 0
    for (i in 2:n)
    {
      for (j in 1:(i - 1))
      {
        if (cluster_assignments[i] != cluster_assignments[j])
        {
          if (traj_dist[i, j] < min_dist)
          {
            min_dist = traj_dist[i, j]
            min_i = i
            min_j = j
          }
        }
      }
    }
    adj_mat[min_i, min_j] = min_dist
    adj_mat[min_j, min_i] = min_dist
    traj_graph = graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted =
                                               TRUE)
    conn_comp = count_components(traj_graph)
  }
  return(traj_graph)
}

#' Select genes to use in building a cell trajectory
#'
#' This function uses "neighborhood variance" to identify genes
#' that vary smoothly, rather than fluctuating randomly, across
#' the set of cells. Genes selected in this way can then be used
#' to construct a trajectory.
#'
#' @param embedding Low-dimensional LLE embedding of cells
#' @return Vector containing indices of selected genes
#' @export
#' @examples
#' \dontrun{
#' genes = select_genes(traj)
#' }
select_genes = function(embedding)
{
  k = min_conn_k(embedding)
  n = nrow(embedding)
  m = ncol(embedding)
  traj_dist = as.matrix(dist(embedding))
  adj_mat = adaptive_knn_graph(traj_dist, rep(k, n))
  sel_vals = sapply(1:m, selection_val, embedding, adj_mat)
  genes = which(sel_vals > 1)
  return(genes)
}

#' Determine the position of each cell within the trajectory
#'
#' This function calculates the geodesic distance from the start
#' cell to each other cell. This value corresponds to the
#' distance a cell has migrated through the process described by
#' the cell trajectory.
#'
#' @param traj_graph Nearest neighbor graph built from LLE embedding
#' @param start Index of starting cell
#' @return Vector of distances
#' @export
#' @examples
#' genes=1:200
#' cells=sample(1:500,30)
#' k=10
#' traj_lle = lle::lle(traj[cells,genes],m=2,k)$Y
#' traj_graph = conn_knn_graph(traj_lle,5)
#' start = 1
#' dists = process_distance(traj_graph,start)
process_distance = function(traj_graph, start)
{
  geodesic_dists = distances(traj_graph, v = start)
  return(geodesic_dists)
}

#' Identify candidate start cells for the trajectory
#'
#' Plots the embedding generated by LLE and highlights
#' potential starting cells for the trajectory. The candidates
#' are chosen based on the longest shortest path through the
#' nearest neighbor graph.
#'
#' @param traj_graph Nearest neighbor graph built from LLE embedding
#' @param embedding Low-dimensional LLE embedding of cells
#' @return Indices of potential starting cells
#' @export
#' @examples
#' genes=1:200
#' cells=sample(1:500,30)
#' k=10
#' traj_lle = lle::lle(traj[cells,genes],m=2,k)$Y
#' traj_graph = conn_knn_graph(traj_lle,5)
#' find_extreme_cells(traj_graph,traj_lle)
find_extreme_cells = function(traj_graph, embedding)
{
  dists = distances(traj_graph)
  starts = unique(max.col(dists))
  plot(
    embedding[, 1],
    embedding[, 2],
    pch = 16,
    xlab = "Manifold Dim 1",
    ylab = "Manifold Dim 2"
  )
  points(embedding[starts[1], 1], embedding[starts[1], 2], pch = 15, col =
           "Red")
  points(embedding[starts[2], 1], embedding[starts[2], 2], pch = 15, col =
           "Yellow")
  legend(
    "left",
    pch = c(15, 15),
    legend = starts,
    col = c("Red", "Yellow")
  )
  return(starts)
}

#' Sort cells according to their progress through a process
#'
#' Uses the values computed by \code{process_distance} to order cells.
#'
#' @param traj_graph Nearest neighbor graph built from LLE embedding
#' @param start Index of starting cell
#' @return Sorted vector of cell indices
#' @export
#' @examples
#' genes=1:200
#' cells=sample(1:500,30)
#' data(traj)
#' k=10
#' traj_lle = lle::lle(traj[cells,genes],m=2,k)$Y
#' traj_graph = conn_knn_graph(traj_lle,5)
#' start=1
#' cells_ordered = cell_order(traj_graph,start)
cell_order = function(traj_graph, start)
{
  geodesic_dists = distances(traj_graph, v = start)
  cell_orders = order(geodesic_dists)
  return(cell_orders)
}

#' Plot trajectory colored by expression level of a gene
#'
#' This function plots the embedding produced by LLE, coloring
#' cells by their expression levels of a gene of interest.
#'
#' @param exp_mat Matrix of expression levels
#' @param embedding Low-dimensional LLE embedding of cells
#' @param samples Indices of cells to include in the plot
#' @param gene_ind Index of gene to use
#' @param cell_symbols Symbols to use for plotting each cell
#' @param title Plot title
#' @return None
#' @export
#' @examples
#' \dontrun{
#' graph_gene(traj,traj_lle,1:nrow(traj),1)
#' }
graph_gene = function(exp_mat,
                      embedding,
                      samples,
                      gene_ind,
                      cell_symbols = 16,
                      title = "Gene Expression")
{
  gene_exp = log(exp_mat[gene_ind, samples] + 1)
  col_scl <-
    (gene_exp - min(gene_exp, na.rm = T)) / (max(gene_exp, na.rm = T) - min(gene_exp, na.rm =
                                                                              T))
  plotclr <-
    colorRampPalette(c("black", "red", "yellow"), space = "rgb")(50)
  color_scl = round(col_scl * length(plotclr))
  color_scl[color_scl == 0] = 1
  plot(
    embedding[, 1],
    embedding[, 2],
    pch = cell_symbols,
    col = plotclr[color_scl],
    xlab = "Manifold Dim 1",
    ylab = "Manifold Dim 2",
    main = title,
    bg = plotclr[color_scl]
  )
  #points(seq(-1,-0.5,length.out=50),rep(0,50),col=plotclr,pch=16,cex=2)
  #legend("bottomleft",pch=c(23,16,15,17),legend=c("Embryonic Day 14.5","Embryonic Day 16.5","Embryonic Day 18.5","Postnatal Day 107"),pt.bg="Black")
}

#' Plot trajectory colored by process distance
#'
#' This function plots the embedding produced by LLE, coloring
#' cells by their progress through a process.
#'
#' @param traj_graph Nearest neighbor graph built from LLE embedding
#' @param embedding Low-dimensional LLE embedding of cells
#' @param start Index of start cell
#' @param cell_symbols Symbols to use for plotting each cell
#' @return None
#' @export
#' @examples
#' genes=1:200
#' cells=sample(1:500,30)
#' k=10
#' traj_lle = lle::lle(traj[cells,genes],m=2,k)$Y
#' traj_graph = conn_knn_graph(traj_lle,5)
#' start=1
#' graph_process_distance(traj_graph,traj_lle,start)
graph_process_distance = function(traj_graph,
                                  embedding,
                                  start,
                                  cell_symbols = 16)
{
  geodesic_dists = process_distance(traj_graph, start)
  col_scl <-
    (geodesic_dists - min(geodesic_dists, na.rm = T)) / (max(geodesic_dists, na.rm =
                                                               T) - min(geodesic_dists, na.rm = T))
  plotclr <-
    colorRampPalette(c("black", "red", "yellow"), space = "rgb")(50)
  color_scl = round(col_scl * length(plotclr))
  color_scl[color_scl == 0] = 1
  plot(
    embedding[, 1],
    embedding[, 2],
    pch = cell_symbols,
    col = plotclr[color_scl],
    xlab = "Manifold Dim 1",
    ylab = "Manifold Dim 2"
  )
  for (i in 1:length(V(traj_graph)))
  {
    if (i == start) {
      i = i + 1
    }
    path = get.shortest.paths(traj_graph, start, i)[[1]]
    path_inds = as.numeric(path[[1]])
    start_inds = path_inds[1:length(path_inds) - 1]
    end_inds = path_inds[2:length(path_inds)]
    segments(
      embedding[start_inds, 1],
      embedding[start_inds, 2],
      embedding[end_inds, 1],
      embedding[end_inds, 2],
      col = rep("black", length(start_inds)),
      lwd = 1
    )
  }
  points(
    seq(-1, -0.5, length.out = 50),
    rep(-1.5, 50),
    col = plotclr,
    pch = 16,
    cex = 2
  )
  #legend("bottomleft",pch=c(23,16,15,17),legend=c("Embryonic Day 14.5","Embryonic Day 16.5","Embryonic Day 18.5","Postnatal Day 107"),pt.bg="Black")
}

#' Compute the geodesic entropy profile of a trajectory
#'
#' The geodesic entropy of a trajectory can be used to detect
#' branches. This function computes geodesic entropy and
#' produces a plot that can be used to visually confirm
#' the branches detected by \code{assign_branches}.
#'
#' @param traj_graph Nearest neighbor graph built from LLE embedding
#' @param start Index of start cell
#' @return Vector of geodesic entropy values. Item k is the
#' geodesic entropy k steps away from the start cell.
#' @export
#' @examples
#' \dontrun{
#' traj_lle = lle(traj[,genes],m=2,k)$Y
#' traj_graph = conn_knn_graph(traj_lle,5)
#' start=1
#' compute_geodesic_entropy(traj_graph,start)
#' }
compute_geodesic_entropy = function(traj_graph, start)
{
  smart_table = function(x, n)
  {
    tab = rep(0, n)
    tab[x] = 1
    return (tab)
  }
  n = length(V(traj_graph))
  paths = shortest_paths(traj_graph, from = start)
  path_table = sapply(paths$vpath, smart_table, n)
  path_counts = rowSums(path_table)
  
  vertex_freqs = matrix(0, n, n)
  for (i in 1:length(paths$vpath)) {
    for (j in 1:length(paths$vpath))
    {
      if (length(paths$vpath[[j]]) >= i) {
        vertex_freqs[i, paths$vpath[[j]][i]] = vertex_freqs[i, paths$vpath[[j]][i]] + 1
      }
    }
  }
  vertex_probs = vertex_freqs / rowSums(vertex_freqs)
  entropy_i = function(vertex_probs, i)
  {
    return (-sum(vertex_probs[i, which(vertex_probs[i, ] > 0)] * log2(vertex_probs[i, which(vertex_probs[i, ] > 0)])))
  }
  max_path_length = max(sapply(paths$vpath, function(x)
    length(x)))
  entropies = rep(Inf, max_path_length)
  for (i in 1:max_path_length)
  {
    entropies[i] = entropy_i(vertex_probs, i)
  }
  plot(
    1:max_path_length,
    entropies,
    type = "l",
    xlab = "Steps from Start Cell",
    ylab = "Geodesic Entropy"
  )
  lines(1:max_path_length, rep(1, max_path_length), lty = 2)
  return (entropies)
}

calc_vertex_freqs <- function(path_counts, path_len) {
  n <- length(path_len)
  vertex_freqs <- matrix(0, n, n)
  for (i in seq_len(n)) {
    vertex_freqs[path_len[i], i] <- path_counts[i]
  }
  vertex_freqs
}

calc_path_counts <- function(vpath, which_ind) {
  counts <- integer(length(which_ind))
  for (i in seq_along(vpath)) {
    index <- which_ind[vpath[[i]]]
    counts[index] <- counts[index] + 1L
  }
  counts
}

as_seq <- function(vec) {
  f <- factor(vec)
  levels(f) <- seq_along(levels(f))
  as.integer(f)
}

#' Detect branches in the trajectory and assign cells to branches
#'
#' This function uses geodesic entropy to automatically determine
#' the number and location of branches in the trajectory.
#' Each cell is then assigned to the corresponding branch.
#'
#' @param traj_graph Nearest neighbor graph built from LLE embedding
#' @param start Index of start cell
#' @param min_branch_len Minimum number of cells required to call a branch
#' @param cells List of indices indicating which cells to assign
#' to branches (used for recursive calls; not intended to be set
#' by users).
#' @return Vector of integers assigning each cell to a branch
#' @export
#' @examples
#' \dontrun{
#' traj_lle = lle::lle(traj[,genes],m=2,k)$Y
#' traj_graph = conn_knn_graph(traj_lle,5)
#' start=1
#' branches = assign_branches(traj_graph,start)
#' plot(traj_lle,pch=16,col=branches)
#' }
assign_branches <- function (traj_graph, start, min_branch_len = 10, cells = V(traj_graph))
{

  n = length(cells)
  which_ind <- integer(length(V(traj_graph)))
  which_ind[cells] <- 1:n
  paths = shortest_paths(traj_graph, from = start, to = cells)
  path_counts <- calc_path_counts(paths$vpath, which_ind)
  path_len <- vapply(paths$vpath, length, FUN.VALUE = integer(1))
  vertex_freqs <- calc_vertex_freqs(path_counts, path_len)
  vertex_probs = vertex_freqs/rowSums(vertex_freqs)
  max_path_length = max(path_len)
  entropies <- vapply(1:max_path_length,
                      function(i, vertex_probs) {
                        probs <- vertex_probs[i,]
                        val <- probs[probs > 0]
                        -sum( val * log2(val))
                      }, FUN.VALUE = numeric(1),
                      vertex_probs = vertex_probs)
  
  branch_point = which(entropies > 1)[1] - 1
  num_branches = round(2^(entropies[branch_point + 1]))
  if (is.na(num_branches) | is.nan(num_branches)) {
    return(rep(1, n))
  }
  crit_verts = order(vertex_probs[branch_point + 1, ], decreasing = TRUE)[1:num_branches]
  num_branches <- num_branches - sum(vertex_freqs[branch_point + 1, crit_verts] < min_branch_len)

  if (num_branches < 2) {
    return(rep(1, n))
  }
  assign_cell = function(x, branch_point, crit_verts) {
    if (length(x) <= branch_point) {
      return(1)
    }
    for (i in 1:length(crit_verts)) {
      if (which_ind[x[branch_point + 1]] == crit_verts[i]) {
        return(i + 1)
      }
    }
    return(1)
  }
  path_long <- path_len >= branch_point + 1
  branch_assignments <- rep(1L, n)
  branch_assignments[path_long] <- cell_assignments[path_long] +
    vapply(paths$vpath[path_long], function(path, table, i) {
      match(path[i], table = table, nomatch = 0L)
    }, FUN.VALUE = integer(1), table = crit_verts, i = branch_point + 1)
  if (max(branch_assignments)-1 > num_branches) {
    # compute once
    cell_dist <- distances(traj_graph, v = cells, to = cells)
    
    while ((max(branch_assignments)- 1) > num_branches) {
      # min_dist <- Inf
      # min_i <- 0
      # min_j <- 0
      branch_splits <- split(1:n, branch_assignments)
      combos <- combn(unique(branch_assignments), 2, simplify = F) |>
        Filter(f = function(x) {! 1 %in% x}, x = _)
      mean_dists <- vapply(
        combos,
        FUN = function(id, dist, split_id){
          mean(dist[split_id[[id[1]]],
                    split_id[[id[2]]]],)
        },
        FUN.VALUE = numeric(1),
        dist = cell_dist,
        split_id = branch_splits)
      selected <- combos[[which.min(mean_dists)]]
      branch_assignments[branch_assignments==max(selected)] <- min(selected)
      branch_assignments <- as_seq(branch_assignments)
    }
    
  }
  while ((max(branch_assignments) - 1) > num_branches) {
    min_dist = Inf
    min_i = 0
    min_j = 0
    browser()
    for (i in 3:max(branch_assignments)) {
      for (j in 2:(i - 1)) {
        mean_dist = mean(distances(traj_graph,
                                   v = cells[branch_assignments == i],
                                   to = cells[branch_assignments == j]))
        if (!is.nan(mean_dist) && mean_dist < min_dist) {
          min_dist = mean_dist
          min_i = i
          min_j = j
        }
      }
    }
    branch_assignments[branch_assignments == min_i] = min_j
  }
  branch_assignments <- as_seq(branch_assignments)
  geodesic_dists = process_distance(traj_graph, start)
  for (i in 1:num_branches) {
    branch_i = which(branch_assignments == (i + 1))
    recurse_br = assign_branches(traj_graph,
                                 start = cells[branch_i[which.min(geodesic_dists[cells[branch_i]])]],
                                 min_branch_len, cells[branch_i])
    rec_num_branches = max(recurse_br)
    if (rec_num_branches > 1) {
      add_this = max(branch_assignments)
      #
      not_one <- recurse_br > 1
      recurse_br[not_one] <- recurse_br[not_one] + add_this
      # for (j in 2:rec_num_branches) {
      #   recurse_br[recurse_br == j] = add_this + j -
      #     1
      # }
      recurse_br[recurse_br == 1] = max(branch_assignments)
      branch_assignments[branch_i] = recurse_br
    }
  }
  as_seq(branch_assignments)
}
