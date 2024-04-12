library(random)
#library(caret)
#library(extraTrees) #uses rJava
#library(tsne)
#library(factoextra)
#library(philentropy)
#related packages: FNN in cran.r-project.org/web/packages/FNN/FNN.pdf

#
#Updates and Comments:
#
#08.15.23 "rjaccard" method added for RELATIVE count-based Jaccard similarity in get_dist()
#          this method is better if features can have drastically disparate counts (some very high, some very low)
#          as it ensures equally-scaled (between 0..1) contribution of each hitting feature into the overall metric

#06.12.23 SE_sample_by_dist() added, which si version of SE_sample() with distance matrix input (faster for smaller matrices)
#
#06.30.22 removed M_OPTI parameter flag from get_dist() as it gives virtually no speed benefits (<1%)
#         tested alternative code for count-based jaccard
#         count_sim() added for count-based similarity metric, but it is actually ~50% slower than the corresponding Jaccard-based distance in get_dist()
#         NB: various metrics were speed-tested in get_dist() with vectorial input (non-vevctorial input is much slower due to unlist() conversion )
#         from fast to slow:  canberra > pearson ~ minkowsky > spearman ~ cosine ~ cosinev >>> tversky > cjaccard >> kendall
#
#06.29.22 verbose mode added to SE_sample() and some speed tweaks
#
#06.28.22 "cjaccard" method added for count-based Jaccard similarity in get_dist()
#         added for SE_sample()  MODE == 3 which samples fixed #nearest neighbors per each seed (may avoid bias toward seeds in dense areas)
#
#10.25.21 embed_SPE_D() added that uses distance matrix as input, but same algorithm as in embed_SPE()
#02.28.21 fix for SE_sample() for MODE == 2
#01.20.21 pull4spot() fix when pulled point matches the spot exactly and causes NAs due to floating precision
#01.16.21 hex2droll() - generator of triangular 2d grid of circular shape
#         TODO: OptiS2Lattice() to code!
#
#12.17.20 push2spot(), pull4spot() - added to locally transform spherical 3D layout
#11.11.20 SE_sample() updated for multi-point seeds and extra mode added
#09.15.20 SE_sample() added for sampling data points by distance
#05.18.20 auxiliary f() added: scale_0center_XY()
#05.14.20 fix for mapXYtoS2() default XY data scaling estimates, gap parameter introduced (theta)

#02.06.20 update for embed_SPE() and expansion for Euclidean only version (embed_SPE_EU() )
#         added OptiS2Lattice() and auxiliary functions for optimizing S2 lattice based on cluster centroids
#
#01.31.20 GetRotMat() - rotation matrix around arbitrary vector u and angle theta
#         mapXYtoS2spot() 
#
#01.30.20 InitS2Lattice() - spherical grid for default clusters form 2 to 32 sizes
#04.16.19 fix in SK_recursive_clust() for 0-memberships when first dirty clustering gives fewer than K clusters
#04.22.19 updated mapXYtoS2() to create customized 3D surface based on 2D input
#


count_sim<-function(V1, V2)
#jaccard similarity for count-based values, V1 and V2 should be vectors
#NB!: this is actually slower than corresponding distance metric inside get_dist() that uses apply()
{
  nn <- min( length(V1), length(V2) )
  ms <- 0 #minimum of corresponding counts
  bs <- 0 #maximum of corresponding counts
  for (ix in 1:nn)
  {
    if (V1[ix] < V2[ix])
    {
      bs <- bs + V2[ix]
      ms <- ms + V1[ix]
    }
    else
    {
      bs <- bs + V1[ix]
      ms <- ms + V2[ix]
    }
  }
  if (bs == 0) return (0) #no overlap in features, so defaults to 0 similarity
  
  return  ( ms / bs  )
}

get_edist2 <- function(V1, V2)
#special case of squared Euclidean distance (faster to calculate)
{
  Z <- V1 - V2
  return ( sum ( Z * Z ) )
}

get_dist <- function(V1, V2, D_TYPE = "minkowski", A_COFF = 2)
  #calculates diverse distance metrics between vectors V1 and V2
  #
  #D_TYPE - the metric type, : "minkowski", "canberra", "tversky", "rjaccard", "cjaccard", "pearson", "cosine", "cosinev", "spearman", "kendall"
  #A_COFF - positive, metric parameter (not used for canberra, as a fraction for Tversky, as a power coefficient for the rest)
  #       - e.g., Minkowski with 1: Manhattan, with 2: Euclidean; Tversky with 1: Tanimoto, with 2: Dice
  #
{
  REV_A <- 1/A_COFF
  
  if (!is.vector(V1)) V1 <- unlist(V1)
  if (!is.vector(V2)) V2 <- unlist(V2)
  if (length(V1) != length(V2)) return (NULL)
  
  type <- tolower(D_TYPE)
  
  if (type == "euclidean2") { # Added for consistency with get_ad.exe default
    return(get_edist2(V1,V2))
  }
  
  if (type == "minkowski")
  {#Manhattan and Euclidean metrics are special cases of this
    Z <- abs(V1 - V2)
    S <- sum ( Z ^ A_COFF )
    return ( S ^ REV_A )
  }
  
  if (type == "canberra")
  {#NB: good for sparse multi-dimensional vectors, because mutually null indices are skipped
    Z <- abs(V1 - V2)
    S <- abs(V1) + abs(V2)
    i <- S > 0
    return ( sum (Z[i] / S[i]) )
  }

  if (type == "rjaccard") #relative count-based jaccard index (extension of continuousTanimoto)
  {
    V12 <- rbind(V1,V2)
    va <- apply(V12, 2, min)
    vb <- apply(V12, 2, max)
    hits <- (vb > 0)
    MS <- length( vb[hits] )
    if (MS <= 0) return (1.0)
    
    S <- sum (va[hits] / vb[hits] ) / MS
    return (1.0 - S)
  } #rjaccard
  
  if (type == "cjaccard") #count-based jaccard index (extension of binary Tanimoto/Jaccard, as below)
  {
    V12 <- rbind(V1,V2)
    MS <- sum ( apply(V12, 2, max) )
    if (MS <= 0) return (1.0)
    
    S <- sum ( apply(V12, 2, min) ) / MS
    return (1.0 - S)
  } #cjaccard
  
  
  if (type == "tversky")
  { #A_COFF = 1 is for Tanimoto/Jaccard, 2 for Dice coefficients. 
    #Assumes binary vetors, but will work for normalized continuous values as well.
    #calculates specific version of Tversky where a = b = (1 / A_COFF)
    #https://en.wikipedia.org/wiki/Jaccard_index#Other_definitions_of_Tanimoto_distance
    V12 <- V1 - V2
    V21 <- -V12
    V12[V12 < 0] <- 0
    V21[V21 < 0] <- 0
    ###V12 and V21 now store complements, i.e. specific features for each vector
    
    VZ   <- V1* V2
    VZ[VZ < 0] <- 0
    SVZ <- sum (VZ)
    
    S <- SVZ / (REV_A*(sum(V12*V12) + sum(V21*V21)) + SVZ)
    
    #alternative (avoids complement calculation), assumes non-negative vectors
    #S <- SVZ / ( REV_A*(sum(V1*V1) + sum(V2*V2)) + SVZ*(1 - 2*REV_A) ) 
    #may need to enforce non-negativiy on V1 and V2 first, as it is enforced on VZ above
    #NB: for arbitrary valued vectors, cosine correlation metric can be used instead
    
    return ( 1 - S )
  }
  
  #remaining are correlation-based distance metrics
  crf <- cor(V1, V2, method = "pearson")
  if (type == "cosine") crf <- 2*sum (V1 * V2) / ( sum(V1*V1) + sum(V2*V2) )
  if (type == "cosinev") crf <- sum (V1 * V2) / sqrt( sum(V1*V1)*sum(V2*V2) ) #vector cosine similarity
  
  if (type == "kendall")  crf <- cor(V1, V2, method = "kendall") #very slow
  if (type == "spearman")  crf <- cor(V1, V2, method = "spearman")
  
  #NB, distance equation below becomes 0 for perfect anti-correlation (= -1), which is assumed desirable behavior
  #to override that, A_COFF can be set to odd value (e.g., 1)
  return ( (1 - crf^A_COFF)*A_COFF/2 ) 
  
} #get_dist()




#non-linear lower dimentional space projection
#Agrafiotis DK et al. stochastic proximity embedding. Molecular Informatics, 2010, 29(11), 758-770
#10.1002/minf.201000134
#
#NB: does not preserve distant order, so large clusters will be scrambled, but closest neighbors are mostly kept
#NB: pretty slow now, but can replace some inner for{} loops with apply() to cure that
#

embed_SPE <-function(FM, LO_D = 2, NCYC = 100, LRN = 0.5, RSEED = NULL, DIST = "Minkowski", DPAR = 2.0)
  #LO_D - 2 or 3 coordinates for lo-dimension projection space; LRN - learning rate
  #FM - feature matrix in high-dimensional space; (DIST, DPAR) define distance metric to be used in FM
  #returns matrix with projected coordinates
{
  N <- dim(FM)[1]
  D <- dim(FM)[2]
  
  if (D <= LO_D) return (NULL)
  
  EPS <- 0.01
  DL <- LRN / NCYC
  if (!is.null(RSEED))
  {
    set.seed(RSEED)
  }
  
  PM <- matrix(0, N, LO_D)
  PM[,1] <- runif(N)
  PM[,2] <- runif(N)
  if (LO_D > 2)  PM[,3] <- runif(N)
  
  
  for (cycles in 1:NCYC)
  {
    for(s in 1:N)
    {
      ij <- sample.int(N, 2)
      Vi <- PM[ij[1],]
      Vj <- PM[ij[2],]
      
      dij <- get_dist(Vi, Vj) #Euclidean for lo-D space
      rij <- get_dist(FM[ij[1],], FM[ij[2],], DIST, DPAR) #specified metric for hi-D space
      #---
      
      if (abs(dij - rij) < EPS) next
      
      K <- LRN*(rij - dij)*(Vj - Vi)/(dij + EPS)
      PM[ij[1],] <- Vi -  K
      PM[ij[2],] <- Vj +  K
    }
    
    LRN <- LRN - DL
  }
  
  return (PM)
} #embed_SPE()


embed_SPE_D <-function(DM, LO_D = 2, NCYC = 100, LRN = 0.5, RSEED = NULL)
  #a version of SPE embedding based on distance matrix as input
  #LO_D - 2 or 3 coordinates for lo-dimension projection space; LRN - learning rate
  #DM - distance matrix (normally, from a in high-dimensional feature space)
  #returns matrix with projected coordinates
{
  N <- dim(DM)[1]
  
  EPS <- 0.01
  DL <- LRN / NCYC
  if (!is.null(RSEED))
  {
    set.seed(RSEED)
  }
  
  PM <- matrix(0, N, LO_D)
  PM[,1] <- runif(N)
  PM[,2] <- runif(N)
  if (LO_D > 2)  PM[,3] <- runif(N)
  
  for (cycles in 1:NCYC)
  {
    for(s in 1:N)
    {
      ij <- sample.int(N, 2)
      Vi <- PM[ij[1],]
      Vj <- PM[ij[2],]
      
      dij <- get_edist2(Vi, Vj)
      rij <- DM[ij[1], ij[2]]
      
      
      if (abs(dij - rij) < EPS) next
      
      K <- LRN*(rij - dij)*(Vj - Vi)/(dij + EPS)
      PM[ij[1],] <- Vi -  K
      PM[ij[2],] <- Vj +  K
    }
    
    LRN <- LRN - DL
    
  }
  
  return (PM)
} #embed_SPE_D()

embed_SPE_EU <-function(FM, LO_D = 2, NCYC = 100, LRN = 0.5, RSEED = NULL, INMEM = TRUE)
  #LO_D - 2 or 3 coordinates for lo-dimension projection space; LRN - learning rate
  #FM - feature matrix in high-dimensional space, Euclidean metric is used
  #returns matrix with projected coordinates
{
  N <- dim(FM)[1]
  D <- dim(FM)[2]
  
  if (D <= LO_D) return (NULL)
  
  EPS <- 0.01
  DL <- LRN / NCYC
  if (!is.null(RSEED))
  {
    set.seed(RSEED)
  }
  
  PM <- matrix(0, N, LO_D)
  PM[,1] <- runif(N)
  PM[,2] <- runif(N)
  if (LO_D > 2)  PM[,3] <- runif(N)
  
  DM <- NULL
  if (INMEM) 
  {
    DM <- as.matrix( dist(FM) ) #takes up memory,  but saves on recalculations
    DM <- DM * DM
  }
  
  for (cycles in 1:NCYC)
  {
    for(s in 1:N)
    {
      ij <- sample.int(N, 2)
      Vi <- PM[ij[1],]
      Vj <- PM[ij[2],]
      
      dij <- get_edist2(Vi, Vj)
      
      #---
      #NB: get_dist() is slow but memory-sparing
      if (INMEM) rij <- DM[ij[1], ij[2]] else   rij <- get_edist2(FM[ij[1],], FM[ij[2],])
      #---
      
      if (abs(dij - rij) < EPS) next
      
      K <- LRN*(rij - dij)*(Vj - Vi)/(dij + EPS)
      PM[ij[1],] <- Vi -  K
      PM[ij[2],] <- Vj +  K
    }
    
    LRN <- LRN - DL
    
    #---
    #plot(PM[,1], PM[,2]) #dbg
    #Sys.sleep(1) #dbg
    #---
  }
  
  return (PM)
} #embed_SPE_EU()

SK_clust_findmin <- function(V, FMref, DIST = "Minkowski", DPAR = 2.0)
  #auxiliary for SK_clust()
{
  #dis <- apply(FMref, 1, get_edist2, V)
  dis <- apply(FMref, 1, get_dist, V, DIST, DPAR)
  return ( which.min(dis) )
}

#----
#K-means stochastic clustering with supplied or random-seeded starting POINTS
SK_clust <- function(FM, POINTS = 10, DIST = "Minkowski", DPAR = 2.0, RSEED = NULL)
  #POINTS can be a vector of indices in FM or a number to sample randomly
  #returns a vector of cluster assignments (from 1 to K)
{
  N <- dim(FM)[1]
  rs <- POINTS
  K  <- length(rs)
  
  if (K == 1)
  {
    if (!is.null(RSEED)) set.seed(RSEED)
    K <- POINTS
    rs <- sample.int(N, K)
  }
  
  K_List <- apply(FM, 1, SK_clust_findmin, FM[rs,], DIST, DPAR)
  
  return (K_List)
}

SK_clust_wseeds <- function(FM, seeds, DIST = "Minkowski", DPAR = 2.0, RSEED = NULL)
  #version with explicit seeds supplied
  #returns a vector of cluster assignments (from 1 to K)
{
  K_list <- NULL
  if (dim(FM)[2] == dim(seeds)[2])
  {
    K_List <- apply(FM, 1, SK_clust_findmin, seeds, DIST, DPAR)
  }
  
  return (K_List)
}

SK_fix_seeds <- function(FM, cluster_list)
{#calculates central points for each cluster in cluster_list[]
  membs <- table(cluster_list)
  K <- length(membs)
  seeds <- matrix(0, K, dim(FM)[2])
  for (i in 1:K)
  {
    seeds[i,] <- apply(FM[cluster_list == i,], 2, median)
  }
  return (seeds)
}

SK_recursive_clust <- function(FM, K = 24, Kmin = 4, nbatch = 64, MaxDepth = NULL, DIST = "Minkowski", DPAR = 2.0, RSEED = NULL)
  #performs recursive paritioning of FM into at most K (and approx. at least Kmin) clusters in each slice
  #FM is a feature matrix with rows as samples to be clustered
  #nbatch is default sampling for each cluster's centroid estimate
  #returns a matrix of clusters (columns) arranged into slices with individual row ids propagated as cluster ids
  #MaxDepth is upper limit for number of levels (slices), if not specified, imputed from the data
{
  N <- dim(FM)[1]
  M <- dim(FM)[2]
  Kup <- 2*K             #the minimum set size to undergo clustering  
  
  if (!is.null(RSEED)) set.seed(RSEED)
  
  Smax <- MaxDepth
  if (is.null(MaxDepth))    
    Smax <- round(1 + log(N) / log(Kmin)) #maximum number of clustering levels to undertake
  
  Scur <- 2 #current level of clustering (level 1 is the supercluster)
  
  Slices <- matrix(0, N, Smax + 1) #stores cluster ids
  Slices[,1] <- 1              #default super-cluster id
  
  Glob <- c(1:N)               #used for global indexing of samples in the dataset
  
  while (Scur < Smax)
  {
    prev <- unique(Slices[,Scur-1])
    for( pid in prev)
    {
      p <- Glob[Slices[,Scur-1] == pid] #p has the current work cluster points mapped onto global ids
      
      if (length(p) < Kup)   
      {#treat as leafs
        Slices[p, Scur] <- p
        next
      }
      
      TopK <- SK_clust(FM[p,], K, DIST=DIST, DPAR=DPAR, RSEED=RSEED) #initial quick and dirty clustering
      tt <- sort(table(TopK))
      
      ttl <- length(tt) #fix 041619, before it was assumed to equal K
      too_small <- min( median(tt), tt[ max(0,ttl-Kmin)+1] ) 
      
      if (length(tt[tt > too_small]) < Kmin)
      {
        tt <- tt[(max(0,ttl-Kmin)+1):ttl]
      }else
      {
        tt <- tt[tt > too_small]
      }
      
      ok_ids <- as.numeric(rownames(as.matrix(tt)))
      K1 <- length(ok_ids)
      
      cmeds <- matrix(0.0, K1, M) #stores K1 median points for the clusters
      for (i in 1:K1)
      {
        id <- ok_ids[i] #cluster ids that point to flags in TopK
        pK1 <- p[TopK == id]
        tsample <- min(length(pK1), nbatch)
        minibatch <- sample(x = pK1, size = tsample)
        if (tsample == 1)
        {
          cmeds[i,] <- unlist( FM[minibatch,] )
        } else 
        {
          #if (tsample < 1) print (paste("Err: K1=", K1, "id=", id, "Scur=", Scur)) #dbg, does not happen
          cmeds[i, ] <- apply( FM[minibatch,], 2, median )
        }
      }
      
      #refinement clustering, just one iteration
      TopK1 <- SK_clust_wseeds(FM[p,], cmeds, DIST=DIST, DPAR=DPAR, RSEED=RSEED)
      
      
      fine_ids <- unique(TopK1)
      K1f <- length(fine_ids) #can get smaller than K1 of the first iteration (e.g., when some small clusters disappear)
      
      medoids <- rep(0, K1f)     #stores K1 cluster centroid ids
      for (i in 1:K1f)
      { 
        id <- fine_ids[i]
        pK1f <- p[TopK1 == id]
        tsample <- min(length(pK1f), nbatch)
        minibatch <- sample(x = pK1f, size = tsample)
        
        if (tsample == 1)
        {
          #cmeds[i,] <- unlist( FM[minibatch,] )
          medoids[i] <- minibatch
        }else 
        {
          cmeds[i,] <- apply(FM[minibatch,], 2, median)
          medoids[i] <- minibatch[ SK_clust_findmin(cmeds[i,], FM[minibatch,], DIST, DPAR) ]  
        }
        
        Slices[pK1f, Scur] <- medoids[i]     #these serve as unique ids for clusters of clusters
      } #for i
      
    } #per pid
    
    tt <- table(Slices[,Scur])
    Scur <- Scur + 1
    
    if (length(tt) == N) break
    
    #if (max(tt) < Kup)
    if (median(tt) < Kmin)  break
  } #while
  
  
  #-------------------------------------------------------------------------
  Scur <- Scur - 1
  ##assign all the nodes as leafs in Scur+1 layer
  prev <- unique(Slices[,Scur])
  if (length(prev) < N)
  {
    for (pid in prev)
    {
      p <- Glob[Slices[,Scur] == pid]
      Slices[p, Scur+1] <- p
    }#for pid
  } else Scur <- Scur-1
  
  
  #-----------dbg 041319 - 041619
  #write.table(Slices, "DSSTOX_clusters_dbg_041619.txt", sep='\t', quote = FALSE)
  #Slices <- read.table("DSSTOX_clusters_dbg_041319.txt", sep='\t')
  #-----------dbg end
  
  #-------------------------------------------------------------------------
  #Slices now contain: [,1] for supercluster and up to Scur+1 for singletons
  
  #The final span of slices will be in reverse order of layers
  fslices <- matrix(0, N, Scur)
  rownames(fslices) <- rownames(FM)
  colnames(fslices) <- paste("clayer#", c(1:Scur), sep='')
  fslices[,1] <- Slices[,Scur+1] 
  
  
  #-------------------------------------------------------------------------
  #now, back-propagate cluster ids for consistency between slices
  for (s in Scur:2)
  {
    us <- unique(Slices[, s]) #iterate over cluster ids in the current slice "s"
    for (p in us)
    {#from parent to child
      p_us <- Slices[,s] == p #get members of cluster p
      sp <- fslices[p_us, Scur - s + 1] #get prev.layer cluster ids for these members
      
      #c <- unique(sp)[1]  #use first member's cluster id
      c <- as.numeric(rownames(as.matrix(sort(table(sp), decreasing = TRUE))))[1] #use most dominant member
      
      fslices[p_us, Scur - s + 2] <- c  #propagate chosen child's cluster id as the next layer's cluster id
    } #p
  } #s
  
  return (fslices)
} #end of SK_recursive_clust()

XYZ_layout <- function(FM, Slices, NDIMS = 3, C_ISO = 0.1, UNI_R = 1000.0, MAX_ITER = 100)
  #creates progressively coordinates for the hierarchy of clusters from the slices matrix
  #NDIMS = 2 or 3, defines if 2D or 3D coordinates are created
  #UNI_R = maximum ditance between points (essentially defines the span of coordinates)
  #C_ISO controls for fractional overlap between clusters (0 - no control, the higher, the more isolated each cluster is)
  #MAX_ITER controls maximum number of iterations for optimizing low-dimensional coordinates to reflect distances in FM space
  #
  #NB: uses embed_SPE() method with Euclidean distances, so may not be applicable if other distance metrics are needed
{
  N <- dim(Slices)[1]
  NS <- dim(Slices)[2]
  if (dim(FM)[1] != N) return (NULL)
  
  Rscale <- 0.1
  if (abs(C_ISO) < 1.0) Rscale <- C_ISO
  
  gXYZr <- matrix(0.0, N, NDIMS + 1) #stores coordinates and cluster radius
  Glob <- c(1:N) #used for global indexing of samples in the dataset
  
  tops <- unique(Slices[,NS])
  
  n_iters <- min(3*length(tops), 100)
  XYZ <- embed_SPE_EU( FM[tops,], NDIMS, n_iters) #generate coordinates
  
  f <- (2.0 - Rscale) * UNI_R / get_dist(apply(XYZ, 2, max), apply(XYZ, 2, min))
  centr <- apply(XYZ, 2, mean)
  for (d in 1:NDIMS)
  {
    gXYZr[tops,d] <- f*(XYZ[,d] - centr[d])
  }
  
  DD <- as.matrix( dist(gXYZr[tops,1:NDIMS]) )
  DD[DD == 0] <- max(DD)
  gXYZr[tops,NDIMS+1] <- apply(DD, 1, min) / (2.0 + Rscale)
  
  #recursively update other layers
  for (sl in NS:2)
  {
    curr <- unique(Slices[,sl])
    for( pid in curr)
    {
      l1 <- Glob[Slices[,sl] == pid] 
      p <- unique(Slices[l1,sl-1]) #p has the current work cluster points mapped onto global ids
      if (length(p) == 1) next
      
      n_iters <- min(3*length(p), 100)
      XYZ <- embed_SPE_EU( FM[p,], NDIMS, n_iters) #generate coordinates
      
      #shifting & scaling
      estim_range <- get_dist(apply(XYZ, 2, max), apply(XYZ, 2, min))
      f <- (2.0 - Rscale)*gXYZr[pid,NDIMS+1] / estim_range
      t <- gXYZr[pid,1:NDIMS] - f*apply(XYZ, 2, mean)
      
      for (d in 1:NDIMS)
      {
        gXYZr[p,d] <- f*XYZ[,d] + t[d]  
      }
      
      DD <- as.matrix( dist(gXYZr[p,1:NDIMS]) )
      DD[DD == 0] <- max(DD)
      gXYZr[p,NDIMS+1] <- apply(DD, 1, min) / (2.0 + Rscale)
      
      #      points3d(x=gXYZr[p,1], y=gXYZr[p,2], z=gXYZr[p,3], type = "s", radius = gXYZr[p,NDIMS+1], xlab = 'x', ylab = 'y', zlab = 'z',
      #             col = hcolr[pid], lwd = 2, Add = FALSE) 
      
      #symbols(x=gXYZr[p,1], y=gXYZr[p,2], circles = gXYZr[p,NDIMS+1], 
      #        bg = bcolr[sl-1], fg = hcolr[p], add = TRUE, inches = FALSE)
      
    } #for pid
    print ( paste("Current level:", sl) ) #---dbg only
  } #for sl
  
  colnames(gXYZr) <- c( c("X", "Y", "Z")[1:NDIMS], "R")
  return (gXYZr[,1:NDIMS])
} #end of XYZ_layout()

KDT_SPLIT <- function(M, cdim, points)
  #auxiliary function for KDT_PCA:
  #M - feature matrix, cdim - current dimension
  #points - indices of rows in FM that are to be split by the median
{
  n <- length(points)
  if (n == 0) return (NULL) #empty branch
  if (n == 1) return ( list(dim=0, pivot = points, left = NULL, right = NULL) ) #leaf node
  
  ix <- sort.list( M[points,cdim] )  
  piv <- ceiling(length(points) / 2)
  
  for (z in piv:n)
  {#adjustments
    if (M[ix[piv],cdim] < M[ix[piv+1],cdim]) break
    
    piv <- piv + 1  #makes sure right branch is strictly greater than pivot, 
    if (piv == n) 
    { #prevents empty right-branch, at the expense of seeding it with a single leaf equal to the pivot,which is ok
      #this helps with inseparable points that are complete duplicates
      piv <- piv - 1
      break
    }
  } #z iterations
  
  return ( list(dim = cdim, pivot = ix[piv], left = KDT_SPLIT(M, cdim+1, ix[1:piv]), right = KDT_SPLIT(M, cdim+1, ix[-1:-piv]) ) )
}


KDT_PCA <-function(FM, minv = 0.95, norm = c("BoxCox", "center", "scale") )
  #attempts to construct a k-d tree in the principal component subspace
  #minv - fraction of total variance that needs to be explained by retained PCs
  #norm - data normalization mode, if NULL - skipped, 
  #NB: caret package Preprocessing() is used for PCA and normalization steps
  #
  #returns Caret object, that can be used in KDT_PCA_NN() for finding approximate nearest neighbors
{
  meth <- c(norm, "pca")
  
  #P <- prcomp(FM) #too bothersome calculating variance proportions from st.dev data in prcomp object
  CaretP = preProcess(FM,  method=meth, thresh = minv)
  
  N <- dim(FM)[1]
  if (log2(N) < CaretP$numComp)
  {
    cat("Need ", 2^CaretP$numComp, " data points to cover ", CaretP$numComp, " PCA dimensions! (Reduce minv threshold?)")
    return (NULL)  
  }
  
  #proceed with tree construction
  TFM <- predict(CaretP, FM)
  return ( list(kdt = KDT_SPLIT(TFM, 1, 1:N), pca = CaretP, mtr = FM) )
}

KDT_PCA_NN <- function(TREE, V, N = 1)
  #finds n nearest neighbor(s), using OBJ, a kdt object generated by KDT_PCA(), and a search point vector V
{
  if (length(V) != dim(TREE$mtr)[2])
  {
    cat("Dimensionality mismatch, expected: ", dim(TREE$mtr)[2], "supplied: ", length(V))
    return (NULL)
  }
  M_PCA <- predict(TREE$pca, TREE$mtr)
  V_PCA <- predict(TREE$pca, V)
  
  if (V_PCA[TREE$kdt$dim] > M_PCA[TREE$kdt$pivot, TREE$kdt$dim])
  {
    TREE$kdt$right
    TREE$kdt$left
  }
}

GET_NNS <- function(V, D, K1, dist_type = "Minkowski", dist_par = 2.0)
{#gets nearest neighbors for vector V in feature matrix D
  xd <- apply(D, 1, get_dist, V, dist_type, dist_par)
  return ( sort(xd, index.return = TRUE)$ix[1:K1])
}

GET_NNS_FROM <- function(MD, RD, K, dist_type = "Minkowski", dist_par = 2.0)
{#for each row in feature matrix MD finds its K nearest neighbor rows FROM feature matrix RD
  ns <- dim(MD)[1]
  VALS <- matrix(0, ns, 2*K) #results
  nix <- 1:K
  for (x in 1:ns)
  {
    nns <- as.vector( GET_NNS(MD[x,], RD, K, dist_type, dist_par) )
    VALS[x, nix] <- nns
    VALS[x, nix+K] <- as.vector( apply(RD[nns,], 1, get_dist, MD[x,], dist_type, dist_par)  )
  }
  
  colnames(VALS) <- c(paste("NN", nix, sep=''), paste("D_", nix, sep=''))
  return ( VALS )
}

GET_NNS_SELF <- function(MD, K = 5, maxchunk = 0, dist_type = "Minkowski", dist_par = 2.0)
{#for each row in MD finds its K nearest neighbor rows
  
  ns <- dim(MD)[1]
  
  out <- rep(0, ns)
  
  if ( (ns > maxchunk) && (maxchunk > 0) )  
  {
    nsk <- max(3, ceiling( ns / maxchunk ) ) #minimum 3 clusters
    out <- SK_clust(MD, nsk, dist_type, dist_par)
  }
  
  VALS <- matrix(0, ns, 2*K) #results
  nix <- 1:K
  for (x in 1:ns)
  {
    rset <- (out[x] == out)
    rset[x] <- FALSE
    
    nns <- as.vector( GET_NNS(MD[x,], MD[rset,], K, dist_type, dist_par) )
    
    remap <- c(1:ns)[rset]
    nns <- remap[nns] #fix for id shift)
    
    VALS[x, nix] <- nns
    VALS[x, nix+K] <- as.vector( apply(MD[nns,], 1, get_dist, MD[x,], dist_type, dist_par)  )
  }
  
  colnames(VALS) <- c(paste("NN", nix, sep=''), paste("D_", nix, sep=''))
  return ( VALS )
}




GET_NN_STATS_aux <- function (V, FM, K1, dist_type = "Minkowski", dist_par = 2.0)
{
  xd <- apply(FM, 1, get_dist, V, dist_type, dist_par)
  return ( sort(xd)[2:K1] )
}

GET_NN_STATS <- function(MD, K = 5, maxchunk = 1024, dist_type = "Minkowski", dist_par = 2.0)
  #M - normalized matrix
  #K - number of nearest neighbors
  #MAXNS - maximum sample size to estimate nearest neighbor variance, not used if <= 0
  #dist_type and dist_par are supplied for get_dist() call for distance metric calculations
  #output is a list with the array of K-nearest neighbor distances along with derivative stats
{
  ns <- dim(MD)[1]
  
  out <- rep(0, ns)
  clst <- 0
  
  if ( (ns > maxchunk) && (maxchunk > 0) )  
  {
    nsk <- max(3, ceiling( ns / maxchunk ) ) #minimum 3 clusters
    out <- SK_clust(MD, nsk, dist_type, dist_par)
    tout <- table(out)
    #tout <- cbind(1:nsk, table(out))
    #tout <- sort.int(table(out), index.return = TRUE)
    
    #1.preset
    cs <- max(0.6*maxchunk, min(tout))
    ce <- min(maxchunk,max(tout))
    
    #2.refine
    cs <- max(tout[tout <= cs])
    ce <- min(tout[tout >= ce])
    
    #3.get qualifying clusters
    clst <- (1:nsk)[(tout >= cs) & (tout <= ce)]
    
    #rset <- sample.int(ns, MAXNS ) #biased
  }
  
  VALS <- NULL
  #process clusters
  for (c in clst)
  {
    rset <- which(out == c)
    V <- as.vector( apply(MD[rset,], 1, GET_NN_STATS_aux, MD[rset,], K+1, dist_type, dist_par)  )
    VALS <- c(VALS, V)
  }
  
  return (   list(s = sd(VALS), m = mean(VALS), v = VALS, k = K, d = dist_type, p = dist_par)  )
}

GET_AD <- function(TESTM, BASEM, K = 5, NN_STATS = NULL)
  #estimates distance-based applicability domain (coverage) for TESTM by BASEM, where
  # TESTM - test samples matrix
  # BASEM - dataset matrix
  # K = number of nearest neighbors to use, normally
  # NN_STATS - nearest neighbor distribution of BASEM to use for TESTM coverage calculation, will be recalculated if NULL
  # --
  # returns a matrix of Z scores for TESTM (with K columns and rows matching TESTM)
{
  if (is.null(NN_STATS))    NN_STATS <- GET_NN_STATS(BASEM, K)
  
  nt <- dim(TESTM)[1]
  if ( is.null(nt) )
  {
    nt <- 1
    TESTM <- matrix(TESTM, nrow = nt, ncol = length(TM))
  }
  
  Zs <- matrix(0, nrow=nt, ncol = K)
  for (c in 1:nt)
  {
    x <- apply(BASEM, 1, get_dist, TESTM[c,], NN_STATS$d, NN_STATS$p)
    
    if (K == 1) x <- min(x) else x <- sort(x)[1:K] 
    Zs[c,] <- (x - NN_STATS$m) / NN_STATS$s
  }
  
  return (Zs)
}



hclust2tree <- function(h, FM)
  #creates binary tree with exlpicit coordinates based on h - hclust() object
{
  nodes
  n <- dim(h$merge)[1] + 1
  #n <- dim(A2$x)[1]
  max( abs( c(h$merge[,1], h$mege[,2]) ) )
}

scale_0center_XY <- function(XY, scale = 1.0)
{#zero-centers and max.scales relative to origin, the input XY coordinates
  nXY <- XY
  nXY[,1] <- XY[,1] - mean(XY[,1])
  nXY[,2] <- XY[,2] - mean(XY[,2])
  
  d <- sqrt(nXY[,1]*nXY[,1] + nXY[,2]*nXY[,2])
  xy_r <- scale / max ( d )
  
  #first, rescale to max length 1
  nXY <- nXY * xy_r
  return (nXY)
}

mapXYtoS2 <-function(XY, r = 1.0, theta = 0.1)
  #returns 3D coordinates that are projections of the 2D plane XY onto a 3D sphere surface
  # XY data are recentered to 0 and rescaled so as max.distance from center is 2r/theta, while av.dist is r
  #projection options (conceptually): 
  #         stereographic (quadratic compression of distances), 
  #         azimouthal (extra parameter for height)
  #         tangent based (slightly less compressive) - current
  #   P1 ( x,y ) projects into: 2r^2/(x^2+y^2+r^2)*( x, y, (x^2+y^2-r^2)/(2r) )
  #
  #   the stretch/compression ratio is 2f/(1+f^2) where f is a coefficient to represent the original 2D distance as f*r
  #     e.g., if there is a very small cluster at a distance f*r from the center of 2D plane, 
  #     then its size will not change if f = 1, stretched if f < 1, squshed if f > 1
#
#   r is sphere's radius, the input data is rescaled to comply with r:
#     the uncovered part of a half-meridian after the projection is assumed to be theta*r, and since
#     the uncovered segment is ~2r^2/L, where L is the radius of the 2D circle containing all XY points, hence
#     L => 2r/theta, while distance compression at L is ~2/(pi*theta) fold
#
#   if theta is 0, no scaling/centering is done on XY and it is projected onto a sphere as is
{
  if (theta > 0)
  {#center and stretch input XY, otherwize skip
    
    XY <- scale_0center_XY(XY)
    
    d_2d <- sqrt(XY[,1]*XY[,1] + XY[,2]*XY[,2]) #recalculate new radial distances
    
    #now, progressively stretch to the new max distance L = 2*r / theta
    pL <- log(2/theta)
    pR0 <- 0.5*r*theta
    stre <- r*exp( pL *(2*d_2d - 1.0) ) / d_2d
    stre[d_2d < 0.5] <- 2*r #lower hemisphere points should not be drastically squashed, so it will be proportional to their 2d distance from center
    #rescale to L
    XY[,1] <- XY[,1] * stre
    XY[,2] <- XY[,2] * stre
  } #end of transforming
  
  
  r2 <- r*r
  
  XYZ <- matrix(0, dim(XY)[1], 3)
  for (i in 1:dim(XY)[1])
  {
    l2 <- XY[i,1]*XY[i,1] + XY[i,2]*XY[i,2]
    cf <- 2*r2/(l2+r2)
    XYZ[i,1] <- cf*XY[i,1]
    XYZ[i,2] <- cf*XY[i,2]
    XYZ[i,3] <- 0.5*cf*(l2-r2)/r
  } #for i
  
  return (XYZ)
}

InitS2Lattice <- function(N = 2)
  #creates a lattice of N points on a 0-centered sphere of radius 1
  #N has to be in [2;32] range
{
  if (N < 2) return (NULL)
  if (N > 32) return (NULL)
  
  p <- matrix(0, N, 3) #X by x,y,z
  colnames(p) <- c("x", "y", "z")
  
  if (N == 2)
  {
    p[1,3] <- -1
    p[2,3] <- 1
  }
  
  if (N == 3)
  {
    a <- sqrt(3)/2
    p[1,] <- c(1, 0, 0)
    p[2,] <- c(-0.5, a, 0)
    p[3,] <- c(-0.5, -a, 0)
  }
  
  if (N == 4)
  {
    a <- sqrt(2)
    b <- a / sqrt(3)
    c <- 1 / 3
    p[1,] <- c(2*a*c, 0, -c)
    p[2,] <- c(-a*c, b, -c)
    p[3,] <- c(-a*c, -b, -c)
    p[4,] <- c(0, 0, 1)
  }
  
  if (N == 5)
  {
    p[1:3,] <- InitSphereLattice(3)
    p[4,] <- c(0, 0, -1)
    p[5,] <- c(0, 0, 1)
  }
  
  if (N == 6)
  {
    p[1,] <- c(1, 0, 0)
    p[2,] <- c(-1, 0, 0)
    p[3,] <- c(0, 1, 0)
    p[4,] <- c(0, -1, 0)
    p[5,] <- c(0, 0, 1)
    p[6,] <- c(0, 0, -1)
  }
  
  if (N == 7) return (InitSphereLattice(8)[-8,])
  
  if (N == 8)
  {#cube
    a <- 1/sqrt(3) #0.57735
    p[1,] <- c(a, a, a)
    p[2,] <- c(-a, a, a)
    p[3,] <- c(a, -a, a)
    p[4,] <- c(-a, -a, a)
    
    p[5,] <- c(a, a, -a)
    p[6,] <- c(-a, a, -a)
    p[7,] <- c(a, -a, -a)
    p[8,] <- c(-a, -a, -a)
  }
  
  if (N == 9) return (InitSphereLattice(12)[c(-1,-9,-12),])
  
  if (N == 10) return (InitSphereLattice(12)[c(-9,-12),])
  
  if (N == 11) return (InitSphereLattice(12)[-12,])
  
  if (N == 12)
  {#icosahedron
    phi <- (1 + sqrt(5))/2
    sc <- sqrt(phi*phi + 1)
    
    p[1,] <- c(0, 1, phi)
    p[2,] <- c(0, 1, -phi)
    p[3,] <- c(0, -1, phi)
    p[4,] <- c(0, -1, -phi)
    
    p[5,] <- c( 1,  phi, 0)
    p[6,] <- c( 1, -phi, 0)
    p[7,] <- c(-1,  phi, 0)
    p[8,] <- c(-1, -phi, 0)
    
    p[9,] <-  c( phi, 0, 1)
    p[10,] <- c(-phi, 0, 1)
    p[11,] <- c( phi, 0, -1)
    p[12,] <- c(-phi, 0, -1)
    p <- p / sc
  }
  
  if (N == 13) return ( InitSphereLattice(14)[-14,] )
  
  if (N == 14)
  {#cube with FCC
    p[1:8,] <- InitSphereLattice(8)
    p[ 9,] <- c( 1, 0, 0)
    p[10,] <- c(-1, 0, 0)
    p[11,] <- c( 0, 1, 0)
    p[12,] <- c( 0,-1, 0)
    p[13,] <- c( 0, 0, 1)
    p[14,] <- c( 0, 0,-1)
  }
  
  #15 - 19, may need to revisit these grids:
  if (N == 15) return ( InitSphereLattice(20)[c(-9,-13,-16,-17,-20),] )
  if (N == 16) return ( InitSphereLattice(20)[c(-13,-16,-17,-20),] )
  if (N == 17) return ( InitSphereLattice(20)[c(-16,-17,-20),] )
  if (N == 18) return ( InitSphereLattice(20)[c(-17,-20),] )
  if (N == 19) return ( InitSphereLattice(20)[-20,] )
  
  if (N == 20)
  {#dodecahedron
    phi <- (1 + sqrt(5))/2
    ihp <- 1/phi
    
    p[1:8,] <- sqrt(3)*InitSphereLattice(8)
    
    p[9 ,] <- c(0, phi, ihp)
    p[10,] <- c(0, phi, -ihp)
    p[11,] <- c(0, -phi, ihp)
    p[12,] <- c(0, -phi, -ihp)
    
    p[13,] <- c( phi, ihp, 0)
    p[14,] <- c( phi, -ihp, 0)
    p[15,] <- c(-phi, ihp, 0)
    p[16,] <- c(-phi, -ihp, 0)
    
    p[17,] <- c( ihp, 0, phi)
    p[18,] <- c(-ihp, 0, phi)
    p[19,] <- c( ihp, 0, -phi)
    p[20,] <- c(-ihp, 0, -phi)
    
    p <- p / sqrt(3)
  }
  
  if (N == 21)  return ( InitSphereLattice(24)[c(-16,-21,-24),] )
  
  if (N == 22)  return ( InitSphereLattice(24)[c(-21,-24),] )
  
  if (N == 23)  return ( InitSphereLattice(24)[-24,] )
  
  if (N == 24)
  {#truncated octahedron
    a <- 1/sqrt(5)
    a2 <- 2*a
    p[1,] <- c(0, a, a2)
    p[2,] <- c(0, -a, a2)
    p[3,] <- c(0, a, -a2)
    p[4,] <- c(0, -a, -a2)
    
    p[5,] <- c(a, 0, a2)
    p[6,] <- c(-a, 0, a2)
    p[7,] <- c(a, 0, -a2)
    p[8,] <- c(-a, 0, -a2)
    
    p[ 9,] <- c(a, a2, 0)
    p[10,] <- c(-a, a2, 0)
    p[11,] <- c(a, -a2, 0)
    p[12,] <- c(-a, -a2, 0)
    
    p[13,] <- c(0, a2, a)
    p[14,] <- c(0, -a2, a)
    p[15,] <- c(0, a2, -a)
    p[16,] <- c(0, -a2, -a)
    
    p[17,] <- c(a2, 0, a)
    p[18,] <- c(-a2, 0, a)
    p[19,] <- c(a2, 0, -a)
    p[20,] <- c(-a2, 0, -a)
    
    p[21,] <- c(a2, a, 0)
    p[22,] <- c(-a2, a, 0)
    p[23,] <- c(a2, -a, 0)
    p[24,] <- c(-a2, -a, 0)
  }
  
  if (N == 25)
  {
    #todo
  }
  
  if (N == 26)  return (InitSphereLattice(32)[c(-22,-23, -26,-27, -30,-31),])
  
  if (N == 27)  return (InitSphereLattice(32)[c(-23, -26,-27, -30,-31),])
  
  if (N == 28)  return (InitSphereLattice(32)[c(-26,-27, -30,-31),])
  
  if (N == 29)  return (InitSphereLattice(32)[c(-27, -30,-31),])
  
  if (N == 30)  return (InitSphereLattice(32)[c(-30,-31),])
  
  if (N == 31)  return (InitSphereLattice(32)[-31,])
  
  if (N == 32)
  {#dodecahedron with FCC icosahedron
    p[1:20,] <- InitSphereLattice(20)
    p[21:32,] <- InitSphereLattice(12)
  }
  
  return ( p )
} #InitSphereLattice


GetRotMat <- function(ux,uy,uz, theta)
{
  #rotation matrix around a unit vector u (ux,uy,uz) by theta angle
  #NB counter-clockwise rotation in right-handed coordinate system  
  
  co <- cos(theta)
  si <- sin(theta)
  ve <- 1 - co  #versine
  
  R <- matrix(0, 3, 3)
  uve <- c(ux,uy,uz) * ve
  
  #R[1,] <- c( ux*ux*ve + co, ux*uy*ve - uz*si, ux*uz*ve + uy*si)
  R[1,] <- ux*uve + c( co, - uz*si, uy*si)
  R[2,] <- uy*uve + c( uz*si, co, -ux*si)
  R[3,] <- uz*uve + c( -uy*si, ux*si, co)
  
  return (R)
}

vec_crossprod <- function(a1,a2,a3, b1,b2,b3)
{#vector product of 3d vectors a and b
  return ( c(a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1) )
}

vec_renorm <- function( v )
{#renormalizes vector v to unit length or returns as is vector if length is non-positive
  sv <- sum(v*v)
  if (sv > 0)  return ( v / sqrt(sv) )
  
  return ( v )
}

mapXYtoS2spot <- function(XY, x, y, z)
{
  #simply projects 2D points onto 3D sphere (of unit radius) at (x,y,z) point
  n <- dim(XY)[1]
  v1 <- c(0,0,-1)
  v2 <- vec_renorm ( c(x,y,z) )
  ax  <- vec_renorm ( vec_crossprod(0,0,-1, x, y, z) )
  RM <- GetRotMat(ax[1], ax[2], ax[3], acos( -v2[3]) )
  
  tXY <- matrix(-1, n, 3)
  tXY[,1:2] <- XY
  XYZ <- t(RM %*% t(tXY))
  for (p in 1:n)
  {
    XYZ[p,] <- vec_renorm(XYZ[p,])
  }
  
  return (XYZ)
}



dist3Dsq <- function(x1,y1,z1, x2,y2,z2)
{
  dx <- (x1 - x2)
  dy <- (y1 - y2)
  dz <- (z1 - z2)
  return (dx*dx + dy*dy + dz*dz)
}

dist3D <- function(x1,y1,z1, x2,y2,z2)
{
  return ( sqrt( dist3Dsq(x1,y1,z1, x2,y2,z2) ) )
}


push2spot <- function (U0S2_XYZ, SPOT, d, lim = pi)
{
  #given 3D points on a unit-radius zero-centered sphere, and an (x,y,z) vector SPOT on it with arc-based size d,
  #pushes all points toward the spot proportionally to their latitude on meridian formed by (x,y,z) and its opposite )
  #
  #NB: default arc limit for pulling is pi (i.e., opposite point to spot is unchanged)
  #NB: if lim is negative (e.g., -1), then the input data point with largest arc will be used as a limit
  np <- dim(U0S2_XYZ)[1]
  upd <- U0S2_XYZ 
  gs <- acos( U0S2_XYZ %*% as.matrix(SPOT, 3, 1) ) #angles between points and spot
  
  nadir <- lim
  if (lim < 0) nadir <- max(gs, na.rm = TRUE)
  if (nadir == 0) return (upd) #no change
  
  sc <- d/(lim-d) #angular scaling
  for (p in 1:np)
  {
    px <- U0S2_XYZ[p,1]
    py <- U0S2_XYZ[p,2]
    pz <- U0S2_XYZ[p,3]
    
    if (nadir < gs[p]) next
    
    vt <- vec_crossprod(SPOT[1],SPOT[2],SPOT[3], px,py,pz)
    if (sum(vt * vt) == 0) next # do nothing for special cases
    vt <- vec_renorm(vt)
    
    g <- (gs[p] - lim)*sc
    if (gs[p] < d) #should not happen normally (point within spot's d area)
    {
      g <- gs[p] - d
    }
    
    RM <- GetRotMat(vt[1],vt[2],vt[3], g) #counter-clockwise rotation, hence angle is inverted
    
    pusht <- t(RM %*% as.matrix(U0S2_XYZ[p,], 3, 1))
    upd[p,] <- pusht[1,]
  }
  
  return (upd)
}


pull4spot <- function (U0S2_XYZ, SPOT, maxd)
{
  #given 3D points on a unit-radius zero-centered sphere, and an (x,y,z) SPOT vector,
  #pulls points away from the spot proportionally to their latitude on meridian formed by (x,y,z) and its opposite pole
  #NB: the pull is relative to arc distance from spot to supplied points, so farthest ones will be pulled to maxd
  
  np <- dim(U0S2_XYZ)[1]
  upd <- U0S2_XYZ 
  tt <- U0S2_XYZ %*% as.matrix(SPOT, 3, 1)
  tt[tt > 1] <- 1
  tt[tt < -1] <- -1
  
  #there should be no missing values due to acos!
  gs <- acos( tt )
  nadir <- max(gs)
  if (nadir == 0) return (upd) #no change
  
  sc <- (maxd - nadir)/nadir #scaling
  
  for (p in 1:np)
  {
    px <- U0S2_XYZ[p,1]
    py <- U0S2_XYZ[p,2]
    pz <- U0S2_XYZ[p,3]
    g <- gs[p] #gamma, angle between current point and spot center
    
    vt <- vec_crossprod(SPOT[1],SPOT[2],SPOT[3], px,py,pz)
    if (sum(vt * vt) == 0) next # do nothing for special cases
    
    vt <- vec_renorm(vt)
    RM <- GetRotMat(vt[1],vt[2],vt[3], g*sc) #counter-clockwise rotation, hence angle is inverted
    
    pusht <- t(RM %*% as.matrix(U0S2_XYZ[p,], 3, 1))
    upd[p,] <- pusht[1,]
  }
  
  return (upd)
}


GetS2LatticeStrain <- function()
{
  #Sum over distances, (Dab - D2ab)/d2ab * sqrt(Na*Nb)
}


OptiS2Lattice <- function(S3D, NSIZ, ORID)
{
  #S3D contains X,Y,Z points on a unit sphere
  #NSIZ contains weights for S3D points, such as cluster sizes
  #ORID - distance matrix in original space
  #DO <- as.matrix(ORID)
}


hex2droll <- function(n = 7)
  #returns x,y coordinates for n points 
  #that unroll from (0,0) counterclockwise on equilateral triangular grid
{
  if (n < 2)
  {
    return (NULL)
  }
  
  #the full d-hex has 1 + 3d(d+1) points
  d <- ceiling(sqrt((n -1)/3)) #maximum hex size we need for this n
  maxn <- 1+3*d*(d+1) + (d+2) #last term is extra edge in addition to filled d-hex, due to roll-pattern
  
  hx_tpl <- c(-1, 1, -2, 0, -1, -1, 1, -1, 2, 0, 1, 1) #template pattern for 5 + 1 edges (with  a single horizontal shift before last edge)
  #--- dx and dy are cell steps for the template hx_tpl()
  dx <- 0.5
  dy <- 0.5*sqrt(3)
  #---
  
  xy <- matrix(0, maxn, 2)
  xy[2,1] <- 2*dx 
  #xy[1,] and xy[2,] store first two points (origin and initial shift)
  #after that, template can be applied iteratively
  
  curr <- 2
  for (hex in 1:d)
  {
    for (i in 1:6) #template application
    {
      ddx <- hx_tpl[2*i-1]*dx
      ddy <- hx_tpl[2*i]*dy
      
      for (ihex in 1:hex)
      {
        xy[curr+1,1] <- xy[curr,1] + ddx
        xy[curr+1,2] <- xy[curr,2] + ddy
        curr <- curr + 1
      }
      
      if (i == 5)
      {
        #shift
        xy[curr+1,1] <- xy[curr,1] + 2*dx
        xy[curr+1,2] <- xy[curr,2]
        curr <- curr + 1
        
      }
    }#for i
    
  } #for hex, keep working while needed
  
  return (xy[1:n,]) #return first n points
}




SE_sample <- function(FM, N = 0.1, MODE = 0, ixSTART = 1, DIST = "Minkowski", DPAR = 2.0, verbose = FALSE)
#exaustive search to sample N  points as rows from FM feature matrix by a family of sphere exclusion algorithms
#Golbraikh & Tropsha 2002 - https://pubmed.ncbi.nlm.nih.gov/12489684/
#NB: get_dist() is used for distance metric calculations
#  
#FM - feature matrix to sample rows from
#N - an exact number of points to sample or a fraction of total #rows (default is 10%)
# ixSTART - index of the first point row (will be always included), 
#           can also be a vector of row ids
#  
# MODE  = 0 for max of min (sampling by diversity), 
#           e.g.,  i-th candidate has largest value for the minimum of its distances to previous i-1 selected points
#       = 1 for min of max (sampling by class similarity)
#           e.g.,  i-th candidate has to have smallest value for the maximum of its distances to previous i-1 selected points
#       = 2 for min of min (sampling by individiual samples similarity)
#           e.g., i-th candidate has to have smallest value for the  minimum of its distances to previous i-1 selected points
#
#       = 3 special mode, to sample N nearest neighbors for each seed by min() dist
#         ixSTART has to be a vector and  if N is a fraction, N is set to 1
#NB: in case of ties, lowest row id is picked
#
#verbose flag enables progress info print out (warnings will be shown anyway)
#-----------------------
#example usage:
#prep a random 2D feature space with 2000 points:
#FF <- matrix(runif(2000), 1000, 2)
#plot(FF[,1], FF[,2])
#
# hah <- SE_sample(FF, N=100) #sample by default (diversity)
# points(FF[hah,1], FF[hah,2], col = "red")
# hah2 <- SE_sample(FF, N=30, MODE = 1, ixSTART = 200) #sample by similarity
# points(FF[hah2,1], FF[hah2,2], col = "blue")
# hah3 <- SE_sample(FF, N=1, MODE = 3, ixSTART = hah2) #sample by individual seeds'similarity
# points(FF[hah3,1], FF[hah3,2], col = "green")
# points(FF[hah4,1], FF[hah4,2], col = "red")
{
  if (verbose)
  {
    print("Verbose mode ON.") #dbg / logging mode
  }
  
  nr <- dim(FM)[1]
  SMP <- 0
  if ( (N > 1) && (N < nr) ) SMP <- N
  if ((N < 1) && (N > 0)) SMP <- round(nr * N, 0)
  
  
  if ((SMP == 0) && (MODE < 3)) return (NULL)
  
  ix <- ixSTART
  li <- length(ix)
  setN <- ix
  
  #----------------------
  if (MODE == 3)
  {#reinterpret SMP as #nearest neighbors to be sampled (def. 1) per each seed point
    
    if ( (N < 1) || (SMP == 0) ) SMP <- 1
    candidates <- (1:nr)[-setN]
    for (z in ix)
    {
      if ((z > nr) || (z < 1))
      {
        print(("Invalid row index in the seed(s).")) #dbg
        return (NULL)
      }
      
      dis <- apply(FM[candidates,], 1, get_dist, FM[z,], DIST, DPAR)  
      
      if (SMP == 1)
      {
        setN <- c(setN, candidates[ which.min(dis) ])
      }
      else
      {
        NNs <- sort(dis, index.return = TRUE)
        setN <- c(setN, candidates[ NNs$ix[1:SMP] ])  
      }
      if (verbose)    print( paste0(z, " seed processed") ) #-- dbg
    }#for z
    
    setN <- unique(setN)
    SMP <- length(setN)
    
    if (verbose)    print( paste0("sampled (with seeds) ", SMP, " points") ) #dbg
    return (setN)
  }
  #end of MODE == 3
  
  
  
  if (li+1 > SMP) return (NULL)
  
  if (li == 1)
  {
    print("single seed...") #dbg
    
    #reset the seed's row index to the default, if invalid
    if (ixSTART > nr) ix <- 1 
    if (ixSTART < 1) ix <- 1
    
    dis <- apply(FM, 1, get_dist, FM[ix,], DIST, DPAR)
    dis[ix] <- min(dis) - 1.0
    cd <- which.max(dis)
    if ( (MODE == 1) || (MODE == 2) )
    { 
      dis[ix] <- 2*max(dis)
      cd <- which.min(dis)
    }
    
    setN <- c(ix, cd)  
  }
  
  li <- length(setN)+1
  
  for (i in li:SMP)
  {
    candidates <- (1:nr)[-setN]
    xx <- NULL
    for (ss in candidates)
    {
      dis <- apply(FM[setN,], 1, get_dist, FM[ss,], DIST, DPAR)
      xdis <- min(dis)
      if (MODE == 1) xdis <- max(dis)
      xx <- c (xx, xdis)
    }
    cdx <- which.max(xx)
    if ( (MODE == 1)||(MODE == 2) ) cdx <- which.min(xx)
    
    setN <- c(setN, candidates[cdx])
    if (verbose)    print( paste0(i, " added: ", candidates[cdx]) ) #-- dbg
  } #for i
  
  return ( setN )
} #end of SE_sample()


SE_sample_by_dist <- function(DM, N = 0.1, MODE = 0, ixSTART = 1, verbose = FALSE)
  #version of the above function that uses distance matrix instead of feature matrix
{
  if (verbose)
  {
    print("Verbose mode ON.") #dbg / logging mode
  }
  
  nr <- dim(DM)[1]
  SMP <- 0
  if ( (N > 1) && (N < nr) ) SMP <- N
  if ((N < 1) && (N > 0)) SMP <- round(nr * N, 0)
  
  
  if ((SMP == 0) && (MODE < 3)) return (NULL)
  
  ix <- ixSTART
  li <- length(ix)
  setN <- ix
  
  #----------------------
  if (MODE == 3)
  {#reinterpret SMP as #nearest neighbors to be sampled (def. 1) per each seed point
    
    if ( (N < 1) || (SMP == 0) ) SMP <- 1
    candidates <- (1:nr)[-setN]
    for (z in ix)
    {
      if ((z > nr) || (z < 1))
      {
        print(("Invalid row index in the seed(s).")) #dbg
        return (NULL)
      }
      
      dis <- DM[candidates,z]
      
      if (SMP == 1)
      {
        setN <- c(setN, candidates[ which.min(dis) ])
      }
      else
      {
        NNs <- sort(dis, index.return = TRUE)
        setN <- c(setN, candidates[ NNs$ix[1:SMP] ])  
      }
      if (verbose)    print( paste0(z, " seed processed") ) #-- dbg
    }#for z
    
    setN <- unique(setN)
    SMP <- length(setN)
    
    if (verbose)    print( paste0("sampled (with seeds) ", SMP, " points") ) #dbg
    return (setN)
  }
  #end of MODE == 3
  
  
  
  if (li+1 > SMP) return (NULL)
  
  if (li == 1)
  {
    print("single seed...") #dbg
    
    #reset the seed's row index to the default, if invalid
    if (ixSTART > nr) ix <- 1 
    if (ixSTART < 1) ix <- 1
    
    dis <- DM[ix,]
    dis[ix] <- min(dis) - 1.0
    cd <- which.max(dis)
    if ( (MODE == 1) || (MODE == 2) )
    { 
      dis[ix] <- 2*max(dis)
      cd <- which.min(dis)
    }
    
    setN <- c(ix, cd)  
  }
  
  li <- length(setN)+1
  
  for (i in li:SMP)
  {
    candidates <- (1:nr)[-setN]
    xx <- NULL
    for (ss in candidates)
    {
      dis <- DM[setN,ss]
      xdis <- min(dis)
      if (MODE == 1) xdis <- max(dis)
      xx <- c (xx, xdis)
    }
    cdx <- which.max(xx)
    if ( (MODE == 1)||(MODE == 2) ) cdx <- which.min(xx)
    
    setN <- c(setN, candidates[cdx])
    if (verbose)    print( paste0(i, " added: ", candidates[cdx]) ) #-- dbg
  } #for i
  
  return ( setN )
} #end of SE_sample_by_dist(), distance matrix version