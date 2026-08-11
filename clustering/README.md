# clustering — provides a broad, self-contained toolkit of unsupervised clustering algorithms

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/clustering.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/clustering)

Package clustering provides a broad, self-contained toolkit of unsupervised
clustering algorithms and their supporting machinery, implemented using only
the Go standard library.

Partitioning methods include Lloyd's k-means with k-means++ seeding (`KMeans`,
`KMeansPlusPlusInit`), a mini-batch variant (`MiniBatchKMeans`), bisecting
k-means (`BisectingKMeans`), Partitioning Around Medoids (`KMedoids`) and
fuzzy c-means (`FuzzyCMeans`). Density-based methods include `DBSCAN` and
`OPTICS` with reachability extraction. `MeanShift` performs mode-seeking with
flat or Gaussian kernels. Hierarchical agglomerative clustering
(`LinkageCluster`) supports single, complete, average, Ward, centroid, median
and weighted linkage via the Lance-Williams recurrence, and produces a
`Dendrogram` that can be cut by cluster count (`Dendrogram.CutTree`) or by
height (`Dendrogram.CutTreeByHeight`). Model-based clustering is provided by
Gaussian mixture EM (`FitGMM`) with full, diagonal and spherical covariances,
soft responsibilities and AIC/BIC scoring. `SpectralClustering` embeds points
with the eigenvectors of the normalised graph Laplacian, using a built-in
symmetric Jacobi eigensolver (`JacobiEigenSymmetric`).

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/clustering"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AdjustedMutualInformation(labelsTrue, labelsPred []int) float64` | AdjustedMutualInformation returns the mutual information adjusted for chance using the expected mutual information under a hypergeometric model of… |
| `func AdjustedRandIndex(labelsTrue, labelsPred []int) float64` | AdjustedRandIndex returns the Rand index corrected for chance, in the range [-1, 1] (typically [0, 1]); 1 indicates identical clusterings and 0… |
| `func AgglomerativeClustering(data [][]float64, k int, linkage Linkage, metric Metric) ([]int, error)` | AgglomerativeClustering is a convenience wrapper that clusters data with the given linkage and metric and directly returns a flat clustering into k… |
| `func AssignClusters(data, centroids [][]float64, metric Metric) []int` | AssignClusters returns, for each sample, the index of the nearest centroid under the given metric. |
| `func AssignClustersInto(data, centroids [][]float64, metric Metric, labels []int) bool` | AssignClustersInto assigns each sample to its nearest centroid, writing the result into labels (which must have length len(data)). |
| `func BestLabelMatching(labelsTrue, labelsPred []int) map[int]int` | BestLabelMatching returns a permutation that relabels labelsPred to best match labelsTrue by greedily pairing predicted clusters to true classes in… |
| `func BetweenClusterSumOfSquares(data [][]float64, labels []int) float64` | BetweenClusterSumOfSquares returns the between-cluster sum of squares: the weighted squared distance of each cluster centroid to the overall mean. |
| `func BrayCurtis(a, b []float64) float64` | BrayCurtis returns the Bray-Curtis dissimilarity between a and b. |
| `func CalinskiHarabaszIndex(data [][]float64, labels []int) float64` | CalinskiHarabaszIndex returns the Calinski-Harabasz index (variance ratio criterion): the ratio of between-cluster dispersion to within-cluster… |
| `func Canberra(a, b []float64) float64` | Canberra returns the Canberra distance between a and b, a weighted version of the Manhattan distance that is sensitive to small changes near zero. |
| `func Centroid(points [][]float64) []float64` | Centroid returns the coordinate-wise mean of the given points. |
| `func Chebyshev(a, b []float64) float64` | Chebyshev returns the Chebyshev (L-infinity, maximum coordinate) distance between a and b. |
| `func Cholesky(a [][]float64) ([][]float64, error)` | Cholesky returns the lower-triangular Cholesky factor L of a symmetric positive-definite matrix a, such that a = L*L^T. It returns an error if a is… |
| `func CloneMatrix(m [][]float64) [][]float64` | CloneMatrix returns a deep copy of the matrix m. |
| `func CloneVector(v []float64) []float64` | CloneVector returns a deep copy of the vector v. |
| `func ClusterIndices(labels []int) map[int][]int` | ClusterIndices groups sample indices by their cluster label, returning a map from label to the sorted list of sample indices assigned to it. |
| `func ClusteringAccuracy(labelsTrue, labelsPred []int) float64` | ClusteringAccuracy returns the fraction of samples correctly labelled after optimally matching predicted clusters to true classes with… |
| `func ColumnMeans(data [][]float64) []float64` | ColumnMeans returns the mean of each column (feature) of data. |
| `func ColumnMinMax(data [][]float64) (min, max []float64)` | ColumnMinMax returns the per-column minimum and maximum of data. |
| `func ColumnStdDevs(data [][]float64) []float64` | ColumnStdDevs returns the population standard deviation of each column of data. |
| `func Completeness(labelsTrue, labelsPred []int) float64` | Completeness returns the completeness score: 1 when all members of a true class are assigned to the same cluster. |
| `func ComputeInertia(data, centroids [][]float64, labels []int, metric Metric) float64` | ComputeInertia returns the sum over all samples of the distance from each sample to its assigned centroid, using the given metric. |
| `func CondensedDistanceMatrix(data [][]float64, metric Metric) []float64` | CondensedDistanceMatrix returns the pairwise distances between rows of data in condensed upper-triangular form: entries for pairs (i, j) with i < j… |
| `func CondensedIndex(n, i, j int) int` | CondensedIndex returns the position in a condensed distance matrix (as produced by CondensedDistanceMatrix) of the pair (i, j) for n points. |
| `func ContingencyMatrix(labelsTrue, labelsPred []int) (table [][]int, trueLabels, predLabels []int)` | ContingencyMatrix returns the contingency table between two labellings of the same samples, together with the sorted distinct labels of each. |
| `func CorrelationDistance(a, b []float64) float64` | CorrelationDistance returns 1 - PearsonCorrelation(a, b), in the range [0, 2]. |
| `func Cosine(a, b []float64) float64` | Cosine returns the cosine distance 1 - CosineSimilarity(a, b), in the range [0, 2]. |
| `func CosineSimilarity(a, b []float64) float64` | CosineSimilarity returns the cosine of the angle between a and b, in the range [-1, 1]. |
| `func CountClusters(labels []int) int` | CountClusters returns the number of distinct non-negative cluster labels, ignoring the noise label -1 used by density-based methods. |
| `func CovarianceMatrix(data [][]float64) [][]float64` | CovarianceMatrix returns the sample covariance matrix of data (dividing by n-1). |
| `func CovarianceMatrixMLE(data [][]float64) [][]float64` | CovarianceMatrixMLE returns the maximum-likelihood covariance matrix of data (dividing by n). |
| `func DaviesBouldinIndex(data [][]float64, labels []int) float64` | DaviesBouldinIndex returns the Davies-Bouldin index of a clustering, the mean over clusters of the worst-case ratio of within-cluster scatter to… |
| `func DegreeMatrix(affinity [][]float64) [][]float64` | DegreeMatrix returns the diagonal degree matrix of an affinity matrix, where each diagonal entry is the row sum of affinities. |
| `func Determinant(a [][]float64) (float64, error)` | Determinant returns the determinant of the square matrix a. |
| `func DistanceMatrix(data [][]float64, metric Metric) [][]float64` | DistanceMatrix returns the full symmetric n x n matrix of pairwise distances between the rows of data using the given metric. |
| `func Dot(a, b []float64) float64` | Dot returns the dot (inner) product of a and b. |
| `func DunnIndex(data [][]float64, labels []int, metric Metric) float64` | DunnIndex returns the Dunn index of a clustering: the ratio of the minimum inter-cluster distance to the maximum intra-cluster diameter. |
| `func ElbowInertias(data [][]float64, kValues []int, rng *rand.Rand) ([]float64, error)` | ElbowInertias runs k-means for each k in kValues and returns the resulting inertia for each, the sequence used to locate an "elbow" in the inertia… |
| `func ElbowPoint(kValues []int, inertias []float64) int` | ElbowPoint returns the index into kValues of the elbow of the inertia curve, chosen as the point of maximum distance to the straight line joining the… |
| `func EstimateBandwidth(data [][]float64, quantile float64, metric Metric) float64` | EstimateBandwidth returns a bandwidth estimate for mean-shift based on the average distance to the k-th nearest neighbour, using the given quantile… |
| `func Euclidean(a, b []float64) float64` | Euclidean returns the Euclidean (L2) distance between a and b. |
| `func FlatKernel(distance, bandwidth float64) float64` | FlatKernel is a uniform (flat) kernel: it returns 1 within the bandwidth and 0 outside, giving each in-window point equal weight. |
| `func ForgyInit(data [][]float64, k int, rng *rand.Rand) [][]float64` | ForgyInit is an alias for RandomInit: it picks k random data points as the initial centroids (the Forgy method). |
| `func FowlkesMallowsIndex(labelsTrue, labelsPred []int) float64` | FowlkesMallowsIndex returns the Fowlkes-Mallows index between two clusterings: the geometric mean of the pairwise precision and recall. |
| `func GapStatistic(data [][]float64, k, nRefs int, rng *rand.Rand) (gap, stderr float64, err error)` | GapStatistic computes the gap statistic for a given k: the difference between the expected log within-cluster dispersion under a uniform reference… |
| `func GaussianKernel(distance, bandwidth float64) float64` | GaussianKernel returns the (unnormalised) Gaussian kernel weight exp(-distance^2 / (2*bandwidth^2)). |
| `func Hamming(a, b []float64) float64` | Hamming returns the number of coordinates in which a and b differ. |
| `func Homogeneity(labelsTrue, labelsPred []int) float64` | Homogeneity returns the homogeneity score: 1 when each cluster contains only members of a single true class. |
| `func Identity(n int) [][]float64` | Identity returns the n x n identity matrix. |
| `func Inertia(data [][]float64, labels []int) float64` | Inertia returns the total within-cluster sum of squared distances of samples to their cluster centroid (the k-means objective). |
| `func Inverse(a [][]float64) ([][]float64, error)` | Inverse returns the inverse of the square matrix a, or ErrSingularMatrix if a is singular. |
| `func JaccardBinary(a, b []float64) float64` | JaccardBinary returns the Jaccard distance between two binary vectors, treating any nonzero coordinate as set. |
| `func JacobiEigenSymmetric(a [][]float64) (eigenvalues []float64, eigenvectors [][]float64, err error)` | JacobiEigenSymmetric computes all eigenvalues and eigenvectors of a real symmetric matrix using the cyclic Jacobi rotation method. |
| `func KDistances(data [][]float64, k int, metric Metric) []float64` | KDistances returns, for each sample, the distance to its k-th nearest neighbour (excluding the point itself). |
| `func KMeansPlusPlusInit(data [][]float64, k int, rng *rand.Rand, metric Metric) [][]float64` | KMeansPlusPlusInit selects k initial centroids from data using the k-means++ seeding strategy, which spreads the initial centres according to squared… |
| `func KNearestAffinityMatrix(data [][]float64, k int, metric Metric) [][]float64` | KNearestAffinityMatrix returns a symmetric k-nearest-neighbour affinity matrix: entry `i``j` is 1 if either i is among the k nearest neighbours of j… |
| `func KNearestNeighbors(point []float64, data [][]float64, k int, metric Metric) []int` | KNearestNeighbors returns the indices of the k rows of data closest to point under the given metric, ordered from nearest to farthest. |
| `func LabelCounts(labels []int) map[int]int` | LabelCounts returns a map from cluster label to the number of samples with that label. |
| `func Mahalanobis(a, b []float64, invCov [][]float64) float64` | Mahalanobis returns the Mahalanobis distance between a and b given the inverse covariance matrix invCov, defined as sqrt((a-b)^T invCov (a-b)). |
| `func MakeBlobs(n, k int, centers [][]float64, std float64, rng *rand.Rand) ([][]float64, []int)` | MakeBlobs generates n samples drawn from k isotropic Gaussian blobs whose centres are placed at the given locations (one row per centre). |
| `func MakeCircles(n int, factor, noise float64, rng *rand.Rand) ([][]float64, []int)` | MakeCircles generates two concentric circles with n samples total, an inner circle of the given factor (0 < factor < 1) of the outer radius, and… |
| `func MakeMoons(n int, noise float64, rng *rand.Rand) ([][]float64, []int)` | MakeMoons generates the classic two interleaving half-circles ("moons") dataset with n samples total and the given Gaussian noise level. |
| `func Manhattan(a, b []float64) float64` | Manhattan returns the Manhattan (L1, taxicab) distance between a and b. |
| `func MatMul(a, b [][]float64) ([][]float64, error)` | MatMul returns the matrix product a*b. |
| `func MatVecMul(m [][]float64, v []float64) ([]float64, error)` | MatVecMul returns the product of matrix m and column vector v. |
| `func Medoid(data [][]float64, indices []int, metric Metric) int` | Medoid returns the index of the point in data that minimises the total distance to all other points in the given index subset, using the supplied… |
| `func MinMaxScale(data [][]float64) [][]float64` | MinMaxScale returns a copy of data with each column linearly rescaled to the range [0, 1]. |
| `func Minkowski(a, b []float64, p float64) float64` | Minkowski returns the Minkowski distance of order p between a and b. |
| `func MultivariateNormalPDF(x, mean []float64, cov [][]float64) (float64, error)` | MultivariateNormalPDF returns the density of a multivariate normal distribution with the given mean and covariance evaluated at x. |
| `func MutualInformation(labelsTrue, labelsPred []int) float64` | MutualInformation returns the mutual information (in nats) between two clusterings viewed as discrete random variables. |
| `func NearestNeighbor(point []float64, data [][]float64, metric Metric) (int, float64)` | NearestNeighbor returns the index of the row of data closest to point under the given metric, together with that distance. |
| `func Norm(v []float64) float64` | Norm returns the Euclidean (L2) norm of the vector v. |
| `func NormP(v []float64, p float64) float64` | NormP returns the Lp norm of the vector v for p > 0. |
| `func NormalizedHamming(a, b []float64) float64` | NormalizedHamming returns the fraction of coordinates in which a and b differ, in the range [0, 1]. |
| `func NormalizedLaplacian(affinity [][]float64) [][]float64` | NormalizedLaplacian returns the symmetric normalised graph Laplacian L_sym = I - D^{-1/2} W D^{-1/2} for the affinity matrix W. |
| `func NormalizedMutualInformation(labelsTrue, labelsPred []int) float64` | NormalizedMutualInformation returns the mutual information normalised by the arithmetic mean of the two cluster entropies, in the range [0, 1]. |
| `func PartitionCoefficient(membership [][]float64) float64` | PartitionCoefficient returns the fuzzy partition coefficient of a membership matrix, in the range [1/c, 1]; values near 1 indicate crisp (well… |
| `func PartitionEntropy(membership [][]float64) float64` | PartitionEntropy returns the fuzzy partition entropy of a membership matrix; lower values indicate crisper clustering. |
| `func PearsonCorrelation(a, b []float64) float64` | PearsonCorrelation returns the Pearson product-moment correlation coefficient between a and b, in the range [-1, 1]. |
| `func Purity(labelsTrue, labelsPred []int) float64` | Purity returns the clustering purity: for each predicted cluster the count of its most common true label is summed and divided by the total number of… |
| `func Quantile(values []float64, q float64) float64` | Quantile returns the q-quantile (0 <= q <= 1) of the values using linear interpolation. |
| `func RBFAffinityMatrix(data [][]float64, gamma float64) [][]float64` | RBFAffinityMatrix returns the Gaussian (RBF) affinity matrix of data, where entry `i``j` = exp(-gamma * \|\|x_i - x_j\|\|^2). |
| `func RandIndex(labelsTrue, labelsPred []int) float64` | RandIndex returns the Rand index between two clusterings: the fraction of sample pairs that are either together in both clusterings or separate in… |
| `func RandomInit(data [][]float64, k int, rng *rand.Rand) [][]float64` | RandomInit selects k distinct data points uniformly at random as initial centroids. |
| `func RegionQuery(data [][]float64, point []float64, eps float64, metric Metric) []int` | RegionQuery returns the indices of samples within distance eps of the query point (including any coincident points), using the given metric. |
| `func RegionQueryAll(data [][]float64, eps float64, metric Metric) [][]int` | RegionQueryAll returns, for each sample, the indices of all samples within distance eps (including the point itself), using the given metric. |
| `func SilhouetteSamples(data [][]float64, labels []int, metric Metric) []float64` | SilhouetteSamples returns the silhouette coefficient of every sample. |
| `func SilhouetteScore(data [][]float64, labels []int, metric Metric) float64` | SilhouetteScore returns the mean silhouette coefficient over all samples, a value in [-1, 1] where higher is better. |
| `func SolveLinearSystem(a [][]float64, b []float64) ([]float64, error)` | SolveLinearSystem solves a*x = b for x, where a is a square nonsingular matrix. |
| `func SortEigen(eigenvalues []float64, eigenvectors [][]float64) ([]float64, [][]float64)` | SortEigen returns the eigenvalues and eigenvectors reordered by ascending eigenvalue. |
| `func SortedKDistances(data [][]float64, k int, metric Metric) []float64` | SortedKDistances returns the k-distances of all points sorted in ascending order, ready to plot as a k-distance graph for eps selection. |
| `func SpectralClustering(data [][]float64, k int, gamma float64, rng *rand.Rand) ([]int, error)` | SpectralClustering partitions data into k clusters using the normalised spectral clustering algorithm of Ng, Jordan and Weiss: it builds an RBF… |
| `func SpectralClusteringFromAffinity(affinity [][]float64, k int, rng *rand.Rand) ([]int, error)` | SpectralClusteringFromAffinity runs normalised spectral clustering starting from a precomputed affinity matrix. |
| `func SquaredEuclidean(a, b []float64) float64` | SquaredEuclidean returns the squared Euclidean distance between a and b. |
| `func Standardize(data [][]float64) [][]float64` | Standardize returns a copy of data with each column transformed to zero mean and unit variance (z-score normalisation). |
| `func TotalSumOfSquares(data [][]float64) float64` | TotalSumOfSquares returns the total sum of squared distances of all samples to the overall mean. |
| `func Trace(a [][]float64) float64` | Trace returns the sum of the diagonal entries of the square matrix a. |
| `func Transpose(m [][]float64) [][]float64` | Transpose returns the transpose of the matrix m. |
| `func UniqueLabels(labels []int) []int` | UniqueLabels returns the sorted distinct labels present in labels. |
| `func UnnormalizedLaplacian(affinity [][]float64) [][]float64` | UnnormalizedLaplacian returns the graph Laplacian L = D - W for the affinity matrix W. |
| `func UpdateCentroids(data [][]float64, labels []int, k int, prev [][]float64) [][]float64` | UpdateCentroids recomputes each centroid as the mean of the samples assigned to it. |
| `func VMeasure(labelsTrue, labelsPred []int) float64` | VMeasure returns the V-measure: the harmonic mean of homogeneity and completeness. |
| `func VectorAdd(a, b []float64) []float64` | VectorAdd returns the element-wise sum a+b. |
| `func VectorAddInPlace(a, b []float64) []float64` | VectorAddInPlace adds b into a element-wise, modifying and returning a. |
| `func VectorScale(v []float64, s float64) []float64` | VectorScale returns v scaled by the scalar s. |
| `func VectorSub(a, b []float64) []float64` | VectorSub returns the element-wise difference a-b. |
| `func WeightedCentroid(points [][]float64, weights []float64) []float64` | WeightedCentroid returns the weighted coordinate-wise mean of points using the given nonnegative weights. |
| `func WeightedEuclidean(a, b, w []float64) float64` | WeightedEuclidean returns the Euclidean distance between a and b with each coordinate weighted by the corresponding entry of w. |
| `func WithinClusterSumOfSquares(data [][]float64, labels []int) map[int]float64` | WithinClusterSumOfSquares returns the sum of squared distances to the centroid for each cluster, keyed by label. |
| `func Zeros(rows, cols int) [][]float64` | Zeros returns a rows x cols matrix filled with zeros. |

### Types

| Type | What it is |
| --- | --- |
| `CovarianceType` | CovarianceType selects the structure imposed on the component covariance matrices of a Gaussian mixture model. |
| `DBSCANResult` | DBSCANResult holds the outcome of a DBSCAN run. |
| `Dendrogram` | Dendrogram is the result of agglomerative clustering: the ordered list of merges and the number of original observations. |
| `FuzzyCMeansResult` | FuzzyCMeansResult holds the outcome of fuzzy c-means clustering. |
| `GMM` | GMM is a fitted Gaussian mixture model. |
| `GMMOptions` | GMMOptions configures Gaussian-mixture EM fitting. |
| `GaussianComponent` | GaussianComponent describes a single multivariate normal component of a mixture: its mixing weight, mean and covariance matrix. |
| `KMeansOptions` | KMeansOptions configures a k-means run. |
| `KMeansResult` | KMeansResult holds the outcome of a k-means clustering run. |
| `KMedoidsResult` | KMedoidsResult holds the outcome of a k-medoids (PAM) clustering run. |
| `Linkage` | Linkage identifies an inter-cluster linkage criterion for agglomerative clustering. |
| `MeanShiftResult` | MeanShiftResult holds the outcome of a mean-shift clustering run. |
| `MergeStep` | MergeStep records a single agglomeration in a dendrogram. |
| `Metric` | Metric is a distance function between two equal-length real vectors. |
| `OPTICSResult` | OPTICSResult holds the ordering and reachability produced by the OPTICS algorithm. |

<details>
<summary><code>DBSCANResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DBSCAN(data [][]float64, eps float64, minPts int, metric Metric) (*DBSCANResult, error)` | DBSCAN runs density-based spatial clustering (DBSCAN) on data. |

</details>

<details>
<summary><code>Dendrogram</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LinkageCluster(data [][]float64, linkage Linkage, metric Metric) (*Dendrogram, error)` | Linkage computes an agglomerative hierarchical clustering of data using the Lance-Williams update formula for the given linkage, and returns the… |
| `func LinkageFromDistanceMatrix(dm [][]float64, linkage Linkage) (*Dendrogram, error)` | LinkageFromDistanceMatrix computes agglomerative clustering from a precomputed symmetric distance matrix. |
| `func (dg *Dendrogram) CopheneticCorrelation(originalCondensed []float64) float64` | CopheneticCorrelation returns the Pearson correlation between the original pairwise distances (condensed) and the cophenetic distances from the… |
| `func (dg *Dendrogram) CopheneticDistances() [][]float64` | CopheneticDistances returns the cophenetic distance matrix: the height at which each pair of original observations is first joined in the dendrogram. |
| `func (dg *Dendrogram) CutTree(k int) ([]int, error)` | CutTree cuts the dendrogram to produce exactly k flat clusters and returns a label per original observation in 0..k-1. |
| `func (dg *Dendrogram) CutTreeByHeight(height float64) []int` | CutTreeByHeight cuts the dendrogram at the given linkage height and returns a flat clustering: merges with distance strictly less than height are… |
| `func (dg *Dendrogram) Leaves(id int) []int` | Leaves returns the original observation indices belonging to the cluster with the given dendrogram id (0..N-1 for leaves, N..2N-2 for internal nodes). |
| `func (dg *Dendrogram) MergeDistances() []float64` | MergeDistances returns the merge distances in the order the merges occurred. |

</details>

<details>
<summary><code>FuzzyCMeansResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FuzzyCMeans(data [][]float64, c int, m float64, maxIter int, tol float64, rng *rand.Rand) (*FuzzyCMeansResult, error)` | FuzzyCMeans runs the fuzzy c-means (soft k-means) algorithm on data with c clusters and fuzziness exponent m > 1 (a common value is 2). |

</details>

<details>
<summary><code>GMM</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FitGMM(data [][]float64, k int) (*GMM, error)` | FitGMM fits a Gaussian mixture model with k components to data using the EM algorithm and default options. |
| `func FitGMMWithOptions(data [][]float64, k int, opts GMMOptions) (*GMM, error)` | FitGMMWithOptions fits a Gaussian mixture model with k components using the supplied options. |
| `func (g *GMM) AIC(data [][]float64) float64` | AIC returns the Akaike information criterion of the fitted mixture on data. |
| `func (g *GMM) BIC(data [][]float64) float64` | BIC returns the Bayesian information criterion of the fitted mixture on data. |
| `func (g *GMM) NumParameters() int` | NumParameters returns the number of free parameters of the fitted mixture, used by the AIC and BIC criteria. |
| `func (g *GMM) Predict(newData [][]float64) []int` | Predict returns the most probable component for each row of newData. |
| `func (g *GMM) PredictProba(newData [][]float64) [][]float64` | PredictProba returns the posterior probability (responsibility) of each component for every row of newData. |
| `func (g *GMM) Score(newData [][]float64) float64` | Score returns the mean per-sample log-likelihood of newData. |
| `func (g *GMM) ScoreSamples(newData [][]float64) []float64` | ScoreSamples returns the per-sample log-likelihood of newData under the fitted mixture. |

</details>

<details>
<summary><code>KMeansResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BisectingKMeans(data [][]float64, k int, rng *rand.Rand) (*KMeansResult, error)` | BisectingKMeans clusters data into k clusters top-down: it starts with all points in one cluster and repeatedly splits the cluster with the largest… |
| `func KMeans(data [][]float64, k int) (*KMeansResult, error)` | KMeans runs k-means clustering on data with the default options, returning the best result over the default number of restarts. |
| `func KMeansWithOptions(data [][]float64, k int, opts KMeansOptions) (*KMeansResult, error)` | KMeansWithOptions runs k-means clustering on data using the supplied options. |
| `func MiniBatchKMeans(data [][]float64, k, batchSize, iterations int, rng *rand.Rand) (*KMeansResult, error)` | MiniBatchKMeans runs a mini-batch variant of k-means, updating centroids from random subsets of size batchSize for the given number of iterations. |
| `func (r *KMeansResult) Predict(newData [][]float64) []int` | KMeansPredict assigns each row of newData to the nearest centroid of a fitted result, returning the cluster labels. |
| `func (r *KMeansResult) Transform(newData [][]float64) [][]float64` | Transform returns, for each row of newData, the Euclidean distance to every centroid of the fitted result. |

</details>

<details>
<summary><code>KMedoidsResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func KMedoids(data [][]float64, k int, metric Metric, rng *rand.Rand) (*KMedoidsResult, error)` | KMedoids runs the Partitioning Around Medoids (PAM) algorithm on data using the given metric (nil means Euclidean). |
| `func KMedoidsFromDistances(dm [][]float64, k int) (*KMedoidsResult, error)` | KMedoidsFromDistances runs PAM directly on a precomputed symmetric distance matrix, avoiding recomputation of pairwise distances. |
| `func (r *KMedoidsResult) Predict(newData [][]float64, metric Metric) []int` | Predict assigns each row of newData to the nearest medoid of a fitted result using the given metric. |

</details>

<details>
<summary><code>Linkage</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LinkageByName(name string) (Linkage, bool)` | LinkageByName returns the Linkage for a name such as "single", "complete", "average"/"upgma", "ward", "centroid", "median" or "weighted". |
| `func (l Linkage) String() string` | String returns the canonical name of the linkage. |

</details>

<details>
<summary><code>MeanShiftResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MeanShift(data [][]float64, bandwidth float64, metric Metric) (*MeanShiftResult, error)` | MeanShift runs the mean-shift clustering algorithm on data using a flat kernel of the given bandwidth and the given metric (nil means Euclidean). |
| `func MeanShiftWithKernel(data [][]float64, bandwidth float64, metric Metric, kernel func(dist, bw float64) float64, maxIter int, tol float64) (*MeanShiftResult, error)` | MeanShiftWithKernel runs mean-shift with a caller-supplied kernel, iteration cap and convergence tolerance. |

</details>

<details>
<summary><code>Metric</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MahalanobisMetric(invCov [][]float64) Metric` | MahalanobisMetric returns a Metric computing the Mahalanobis distance using the given inverse covariance matrix. |
| `func MetricByName(name string) (Metric, error)` | MetricByName returns a Metric for a named distance function. |
| `func MinkowskiMetric(p float64) Metric` | MinkowskiMetric returns a Metric computing the Minkowski distance of order p. |
| `func NewMahalanobisMetric(cov [][]float64) (Metric, error)` | NewMahalanobisMetric builds a Mahalanobis Metric from a covariance matrix by inverting it. |
| `func WeightedEuclideanMetric(w []float64) Metric` | WeightedEuclideanMetric returns a Metric computing the weighted Euclidean distance with the given per-coordinate weights. |

</details>

<details>
<summary><code>OPTICSResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func OPTICS(data [][]float64, maxEps float64, minPts int, metric Metric) (*OPTICSResult, error)` | OPTICS computes the OPTICS ordering of data, generalising DBSCAN across a range of eps values. |
| `func (o *OPTICSResult) ExtractDBSCAN(eps float64) []int` | ExtractDBSCAN extracts a flat DBSCAN-style clustering from an OPTICS result at the given eps threshold, returning a label per original sample (noise… |

</details>

### Constants

`NoiseLabel`

### Variables

`ErrDimensionMismatch`, `ErrEmptyData`, `ErrInvalidK`, `ErrSingularMatrix`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/clustering).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
