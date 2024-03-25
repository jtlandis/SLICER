## SLICER 0.2.0

## Background

SLICER is an algorithm for constructing trajectories that describe gene expression changes during a sequential biological process. SLICER can capture highly nonlinear gene expression changes, automatically select genes related to the process, and detect multiple branch and loop features in the trajectory. 

## Installation

You can also install SLICER from GitHub using devtools
```{r,eval=FALSE}
library("devtools")
install_github("jtlandis/SLICER")
```

## Sample Data and Code
A sample dataset containing 500 simulated "cells" each expressing 300 "genes" is included with the SLICER R package. The example below shows how to run SLICER on this sample data. Note that documentation for each function is available from within R.

```{r,eval=FALSE}
library(SLICER)
genes = select_genes(traj)
k = select_k(traj[,genes], kmin=5)
traj_lle = lle(traj[,genes], m=2, k)$Y
traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = 1
cells_ordered = cell_order(traj_graph, start)
branches = assign_branches(traj_graph,start)
```

## A Few Notes on Using SLICER for Trajectory Construction
1. The select_k function returns a value of k, the number of nearest neighbors to use in dimensionality reduction by locally linear embedding. Although this method of selecting k generally gives good results, we have found that, in some cases, it is necessary to manually tune the value of this parameter. Thus, best practice is to visually examine LLE plots for a range of k values. Also, note that there is a separate nearest neighbor parameter that determines the number of edges in the k-nearest neighbor graph that SLICER builds in the low-dimensional LLE space. The select_k function does not select the value of this parameter. We generally set this parameter at 5 (as in the code snippet above), but occasionally, we find that tweaking it slightly improves the results. SLICER is fundamentally an exploratory, unsupervised analysis tool, so setting of parameter values should always be guided by careful consideration of the biological sensibility of results.

2. SLICER does not restrict the dimensionality of the low-dimensional projection to 2. We have found that in some cases, particularly when the dataset involves many cell fates, using a more high-dimensional projection can improve results.

3. SLICER was initially developed with single cell RNA-seq in mind, but we have also found that it can give good results when applied to other types of data, including bulk RNA-seq, single cell qPCR data, and single cell epigenomic data.

4. A key assumption of SLICER is that a possibly branching, intrinsically one-dimensional process is the dominant source of variation among samples. It is important to think carefully about whether this assumption is reasonable for any given dataset. 

[![](http://cranlogs.r-pkg.org/badges/grand-total/SLICER)](http://cran.rstudio.com/web/packages/SLICER/index.html)
