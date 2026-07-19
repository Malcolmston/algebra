// Package clustering provides a broad, self-contained toolkit of unsupervised
// clustering algorithms and their supporting machinery, implemented using only
// the Go standard library.
//
// Partitioning methods include Lloyd's k-means with k-means++ seeding
// ([KMeans], [KMeansPlusPlusInit]), a mini-batch variant ([MiniBatchKMeans]),
// bisecting k-means ([BisectingKMeans]), Partitioning Around Medoids
// ([KMedoids]) and fuzzy c-means ([FuzzyCMeans]). Density-based methods include
// [DBSCAN] and [OPTICS] with reachability extraction. [MeanShift] performs
// mode-seeking with flat or Gaussian kernels. Hierarchical agglomerative
// clustering ([LinkageCluster]) supports single, complete, average, Ward,
// centroid, median and weighted linkage via the Lance-Williams recurrence, and
// produces a [Dendrogram] that can be cut by cluster count ([Dendrogram.CutTree])
// or by height ([Dendrogram.CutTreeByHeight]). Model-based clustering is
// provided by Gaussian mixture EM ([FitGMM]) with full, diagonal and spherical
// covariances, soft responsibilities and AIC/BIC scoring. [SpectralClustering]
// embeds points with the eigenvectors of the normalised graph Laplacian, using
// a built-in symmetric Jacobi eigensolver ([JacobiEigenSymmetric]).
//
// A rich set of distance metrics is available — Euclidean, squared Euclidean,
// Manhattan, Chebyshev, Minkowski, cosine, correlation, Canberra, Bray-Curtis,
// Hamming, Jaccard, weighted Euclidean and Mahalanobis — all satisfying the
// [Metric] function type so they can be passed interchangeably to the
// algorithms.
//
// Internal ([SilhouetteScore], [DaviesBouldinIndex], [CalinskiHarabaszIndex],
// [DunnIndex], [Inertia], [ElbowPoint], [GapStatistic]) and external
// ([RandIndex], [AdjustedRandIndex], [NormalizedMutualInformation],
// [AdjustedMutualInformation], [FowlkesMallowsIndex], [VMeasure], [Purity],
// [ClusteringAccuracy]) validation measures let clusterings be scored and
// compared. Preprocessing helpers ([Standardize], [MinMaxScale],
// [CovarianceMatrix]) and synthetic dataset generators ([MakeBlobs],
// [MakeMoons], [MakeCircles]) round out the package.
//
// Randomised routines accept a *math/rand.Rand; passing nil uses a deterministic
// source so results are reproducible. Every other routine is deterministic and
// depends only on the standard library.
package clustering
