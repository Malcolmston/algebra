package clustering

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool {
	return math.Abs(a-b) <= eps
}

// twoBlobs returns a small well-separated two-cluster dataset.
func twoBlobs() [][]float64 {
	return [][]float64{
		{0, 0}, {0.1, 0.1}, {-0.1, 0.05}, {0.05, -0.1},
		{10, 10}, {10.1, 9.9}, {9.9, 10.1}, {10.05, 10.05},
	}
}

func TestDistanceMetrics(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{4, 6, 3}
	tests := []struct {
		name string
		got  float64
		want float64
		eps  float64
	}{
		{"euclidean", Euclidean(a, b), 5, tol},
		{"sqeuclidean", SquaredEuclidean(a, b), 25, tol},
		{"manhattan", Manhattan(a, b), 7, tol},
		{"chebyshev", Chebyshev(a, b), 4, tol},
		{"minkowski_p1", Minkowski(a, b, 1), 7, tol},
		{"minkowski_p2", Minkowski(a, b, 2), 5, tol},
		{"dot", Dot(a, b), 4 + 12 + 9, tol},
		{"norm", Norm([]float64{3, 4}), 5, tol},
		{"normp_inf", NormP([]float64{-3, 4, 1}, math.Inf(1)), 4, tol},
		{"hamming", Hamming(a, b), 2, tol},
		{"normhamming", NormalizedHamming(a, b), 2.0 / 3.0, tol},
		{"cosine_same", Cosine([]float64{1, 0}, []float64{2, 0}), 0, tol},
		{"cosine_orth", Cosine([]float64{1, 0}, []float64{0, 1}), 1, tol},
		{"cossim_orth", CosineSimilarity([]float64{1, 0}, []float64{0, 1}), 0, tol},
		{"pearson_perfect", PearsonCorrelation([]float64{1, 2, 3}, []float64{2, 4, 6}), 1, tol},
		{"pearson_neg", PearsonCorrelation([]float64{1, 2, 3}, []float64{3, 2, 1}), -1, tol},
		{"canberra", Canberra([]float64{1, 0}, []float64{0, 1}), 2, tol},
		{"braycurtis", BrayCurtis([]float64{1, 2}, []float64{3, 4}), 4.0 / 10.0, tol},
		{"jaccard", JaccardBinary([]float64{1, 1, 0, 0}, []float64{1, 0, 1, 0}), 1 - 1.0/3.0, tol},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if !approx(tc.got, tc.want, tc.eps) {
				t.Errorf("%s = %v, want %v", tc.name, tc.got, tc.want)
			}
		})
	}
}

func TestWeightedAndMahalanobis(t *testing.T) {
	a := []float64{0, 0}
	b := []float64{3, 4}
	if got := WeightedEuclidean(a, b, []float64{1, 1}); !approx(got, 5, tol) {
		t.Errorf("weighted euclidean = %v, want 5", got)
	}
	// Identity inverse covariance -> Euclidean distance.
	inv := Identity(2)
	if got := Mahalanobis(a, b, inv); !approx(got, 5, tol) {
		t.Errorf("mahalanobis(identity) = %v, want 5", got)
	}
	// Scaled covariance.
	cov := [][]float64{{4, 0}, {0, 4}}
	m, err := NewMahalanobisMetric(cov)
	if err != nil {
		t.Fatal(err)
	}
	if got := m(a, b); !approx(got, 2.5, tol) {
		t.Errorf("mahalanobis scaled = %v, want 2.5", got)
	}
}

func TestMetricByName(t *testing.T) {
	for _, name := range []string{"euclidean", "Manhattan", "city_block", "cosine", "hamming"} {
		if _, err := MetricByName(name); err != nil {
			t.Errorf("MetricByName(%q) error: %v", name, err)
		}
	}
	if _, err := MetricByName("nope"); err == nil {
		t.Error("expected error for unknown metric")
	}
}

func TestLinearAlgebra(t *testing.T) {
	a := [][]float64{{4, 3}, {6, 3}}
	det, err := Determinant(a)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(det, -6, tol) {
		t.Errorf("det = %v, want -6", det)
	}
	inv, err := Inverse(a)
	if err != nil {
		t.Fatal(err)
	}
	prod, _ := MatMul(a, inv)
	id := Identity(2)
	for i := range prod {
		for j := range prod[i] {
			if !approx(prod[i][j], id[i][j], 1e-9) {
				t.Errorf("A*A^-1 not identity at %d,%d: %v", i, j, prod[i][j])
			}
		}
	}
	// Solve linear system: [[2,1],[1,3]] x = [3,5] -> x = [0.8, 1.4].
	x, err := SolveLinearSystem([][]float64{{2, 1}, {1, 3}}, []float64{3, 5})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(x[0], 0.8, 1e-9) || !approx(x[1], 1.4, 1e-9) {
		t.Errorf("solve = %v, want [0.8 1.4]", x)
	}
	// Cholesky reconstruction.
	spd := [][]float64{{4, 2}, {2, 3}}
	l, err := Cholesky(spd)
	if err != nil {
		t.Fatal(err)
	}
	recon, _ := MatMul(l, Transpose(l))
	for i := range spd {
		for j := range spd[i] {
			if !approx(recon[i][j], spd[i][j], 1e-9) {
				t.Errorf("cholesky recon mismatch at %d,%d", i, j)
			}
		}
	}
	if _, err := Inverse([][]float64{{1, 2}, {2, 4}}); err == nil {
		t.Error("expected singular matrix error")
	}
}

func TestKMeans(t *testing.T) {
	data := twoBlobs()
	res, err := KMeans(data, 2)
	if err != nil {
		t.Fatal(err)
	}
	if res.K != 2 {
		t.Fatalf("K = %d", res.K)
	}
	// First four and last four points should share a cluster respectively.
	if res.Labels[0] != res.Labels[1] || res.Labels[0] != res.Labels[2] || res.Labels[0] != res.Labels[3] {
		t.Errorf("first blob not grouped: %v", res.Labels)
	}
	if res.Labels[4] != res.Labels[5] || res.Labels[4] != res.Labels[6] || res.Labels[4] != res.Labels[7] {
		t.Errorf("second blob not grouped: %v", res.Labels)
	}
	if res.Labels[0] == res.Labels[4] {
		t.Errorf("blobs should be in different clusters")
	}
	if res.Inertia > 0.1 {
		t.Errorf("inertia too high for tight blobs: %v", res.Inertia)
	}
	// Predict on a new point near the origin blob.
	pred := res.Predict([][]float64{{0, 0}})
	if pred[0] != res.Labels[0] {
		t.Errorf("predict mismatch")
	}
}

func TestKMeansErrors(t *testing.T) {
	if _, err := KMeans(nil, 2); err == nil {
		t.Error("expected empty data error")
	}
	if _, err := KMeans(twoBlobs(), 0); err == nil {
		t.Error("expected invalid k error")
	}
	if _, err := KMeans(twoBlobs(), 100); err == nil {
		t.Error("expected invalid k error for k>n")
	}
}

func TestKMeansKnownCentroids(t *testing.T) {
	// Two clusters on a line; centroids should be the means of each group.
	data := [][]float64{{0}, {2}, {10}, {12}}
	res, err := KMeansWithOptions(data, 2, KMeansOptions{NInit: 20, Rand: rand.New(rand.NewSource(7))})
	if err != nil {
		t.Fatal(err)
	}
	cents := []float64{res.Centroids[0][0], res.Centroids[1][0]}
	sort.Float64s(cents)
	if !approx(cents[0], 1, 1e-9) || !approx(cents[1], 11, 1e-9) {
		t.Errorf("centroids = %v, want [1 11]", cents)
	}
	if !approx(res.Inertia, 4, 1e-9) {
		t.Errorf("inertia = %v, want 4", res.Inertia)
	}
}

func TestMiniBatchKMeans(t *testing.T) {
	data := twoBlobs()
	res, err := MiniBatchKMeans(data, 2, 4, 200, rand.New(rand.NewSource(3)))
	if err != nil {
		t.Fatal(err)
	}
	if CountClusters(res.Labels) != 2 {
		t.Errorf("expected 2 clusters, got %d", CountClusters(res.Labels))
	}
}

func TestKMedoids(t *testing.T) {
	data := twoBlobs()
	res, err := KMedoids(data, 2, Euclidean, rand.New(rand.NewSource(1)))
	if err != nil {
		t.Fatal(err)
	}
	if len(res.MedoidIndices) != 2 {
		t.Fatalf("expected 2 medoids")
	}
	// Medoids should be actual data points, one per blob.
	m0 := res.MedoidIndices[0]
	m1 := res.MedoidIndices[1]
	if (m0 < 4) == (m1 < 4) {
		t.Errorf("medoids not split across blobs: %v", res.MedoidIndices)
	}
	if res.Labels[0] == res.Labels[4] {
		t.Errorf("blobs not separated by k-medoids")
	}
}

func TestHierarchicalLinkages(t *testing.T) {
	data := twoBlobs()
	for _, lk := range []Linkage{SingleLinkage, CompleteLinkage, AverageLinkage, WardLinkage, CentroidLinkage, MedianLinkage, WeightedLinkage} {
		dg, err := LinkageCluster(data, lk, Euclidean)
		if err != nil {
			t.Fatal(err)
		}
		if len(dg.Merges) != len(data)-1 {
			t.Errorf("%v: expected %d merges, got %d", lk, len(data)-1, len(dg.Merges))
		}
		labels, err := dg.CutTree(2)
		if err != nil {
			t.Fatal(err)
		}
		if CountClusters(labels) != 2 {
			t.Errorf("%v: expected 2 clusters, got %d", lk, CountClusters(labels))
		}
		if labels[0] == labels[4] {
			t.Errorf("%v: blobs not separated", lk)
		}
	}
}

func TestHierarchicalKnownMerge(t *testing.T) {
	// Points at 0,1,4 on a line. Single linkage merges 0&1 (dist 1) first.
	data := [][]float64{{0}, {1}, {4}}
	dg, err := LinkageCluster(data, SingleLinkage, Euclidean)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(dg.Merges[0].Distance, 1, tol) {
		t.Errorf("first merge distance = %v, want 1", dg.Merges[0].Distance)
	}
	if !approx(dg.Merges[1].Distance, 3, tol) {
		t.Errorf("second merge distance = %v, want 3", dg.Merges[1].Distance)
	}
	// Cut into 2 clusters: {0,1} and {4}.
	labels, _ := dg.CutTree(2)
	if labels[0] != labels[1] || labels[0] == labels[2] {
		t.Errorf("cut labels = %v", labels)
	}
}

func TestCutTreeByHeight(t *testing.T) {
	data := [][]float64{{0}, {1}, {10}, {11}}
	dg, _ := LinkageCluster(data, SingleLinkage, Euclidean)
	labels := dg.CutTreeByHeight(5)
	if CountClusters(labels) != 2 {
		t.Errorf("expected 2 clusters cutting at height 5, got %d: %v", CountClusters(labels), labels)
	}
}

func TestCopheneticCorrelation(t *testing.T) {
	data := twoBlobs()
	dg, _ := LinkageCluster(data, AverageLinkage, Euclidean)
	cond := CondensedDistanceMatrix(data, Euclidean)
	c := dg.CopheneticCorrelation(cond)
	if c < 0.5 || c > 1.0001 {
		t.Errorf("cophenetic correlation out of expected range: %v", c)
	}
}

func TestDBSCAN(t *testing.T) {
	data := [][]float64{
		{0, 0}, {0.2, 0}, {0, 0.2}, {0.2, 0.2},
		{10, 10}, {10.2, 10}, {10, 10.2},
		{50, 50}, // outlier
	}
	res, err := DBSCAN(data, 0.5, 3, Euclidean)
	if err != nil {
		t.Fatal(err)
	}
	if res.NClusters != 2 {
		t.Errorf("expected 2 clusters, got %d: %v", res.NClusters, res.Labels)
	}
	if res.Labels[7] != NoiseLabel {
		t.Errorf("outlier should be noise, got %d", res.Labels[7])
	}
	if res.Labels[0] == res.Labels[4] {
		t.Errorf("distinct clusters merged")
	}
}

func TestOPTICS(t *testing.T) {
	data := [][]float64{
		{0, 0}, {0.2, 0}, {0, 0.2},
		{10, 10}, {10.2, 10}, {10, 10.2},
	}
	res, err := OPTICS(data, math.Inf(1), 2, Euclidean)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Order) != len(data) {
		t.Errorf("order length = %d", len(res.Order))
	}
	labels := res.ExtractDBSCAN(0.5)
	if CountClusters(labels) != 2 {
		t.Errorf("OPTICS extraction expected 2 clusters, got %d: %v", CountClusters(labels), labels)
	}
}

func TestMeanShift(t *testing.T) {
	data := twoBlobs()
	res, err := MeanShift(data, 3.0, Euclidean)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.ClusterCenters) != 2 {
		t.Errorf("expected 2 centres, got %d", len(res.ClusterCenters))
	}
	if res.Labels[0] == res.Labels[4] {
		t.Errorf("blobs not separated by mean-shift")
	}
}

func TestGMM(t *testing.T) {
	data := twoBlobs()
	g, err := FitGMM(data, 2)
	if err != nil {
		t.Fatal(err)
	}
	if len(g.Components) != 2 {
		t.Fatalf("expected 2 components")
	}
	labels := g.Predict(data)
	if labels[0] == labels[4] {
		t.Errorf("GMM did not separate blobs: %v", labels)
	}
	// Weights sum to 1.
	var sum float64
	for _, c := range g.Components {
		sum += c.Weight
	}
	if !approx(sum, 1, 1e-6) {
		t.Errorf("weights sum = %v, want 1", sum)
	}
	// Responsibilities per row sum to 1.
	proba := g.PredictProba(data)
	for i, p := range proba {
		var s float64
		for _, v := range p {
			s += v
		}
		if !approx(s, 1, 1e-6) {
			t.Errorf("row %d responsibilities sum = %v", i, s)
		}
	}
	if g.BIC(data) == 0 {
		t.Errorf("BIC should be nonzero")
	}
}

func TestMultivariateNormalPDF(t *testing.T) {
	// Standard bivariate normal at the mean: 1/(2*pi).
	pdf, err := MultivariateNormalPDF([]float64{0, 0}, []float64{0, 0}, Identity(2))
	if err != nil {
		t.Fatal(err)
	}
	if !approx(pdf, 1/(2*math.Pi), 1e-9) {
		t.Errorf("pdf at mean = %v, want %v", pdf, 1/(2*math.Pi))
	}
}

func TestSilhouette(t *testing.T) {
	data := twoBlobs()
	labels := []int{0, 0, 0, 0, 1, 1, 1, 1}
	s := SilhouetteScore(data, labels, Euclidean)
	if s < 0.9 {
		t.Errorf("silhouette for clean split too low: %v", s)
	}
	// A single sample per cluster edge case handled without panic.
	_ = SilhouetteSamples(data, labels, Euclidean)
}

func TestEvaluationIndices(t *testing.T) {
	data := twoBlobs()
	good := []int{0, 0, 0, 0, 1, 1, 1, 1}
	bad := []int{0, 1, 0, 1, 0, 1, 0, 1}
	if DaviesBouldinIndex(data, good) >= DaviesBouldinIndex(data, bad) {
		t.Errorf("Davies-Bouldin should be lower for good clustering")
	}
	if CalinskiHarabaszIndex(data, good) <= CalinskiHarabaszIndex(data, bad) {
		t.Errorf("Calinski-Harabasz should be higher for good clustering")
	}
	if DunnIndex(data, good, Euclidean) <= DunnIndex(data, bad, Euclidean) {
		t.Errorf("Dunn index should be higher for good clustering")
	}
	// Inertia decomposition: total = within + between.
	tot := TotalSumOfSquares(data)
	within := Inertia(data, good)
	between := BetweenClusterSumOfSquares(data, good)
	if !approx(tot, within+between, 1e-6) {
		t.Errorf("SS decomposition failed: total=%v within=%v between=%v", tot, within, between)
	}
}

func TestElbow(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	data, _ := MakeBlobs(90, 3, nil, 0.5, rng)
	ks := []int{1, 2, 3, 4, 5, 6}
	inertias, err := ElbowInertias(data, ks, rand.New(rand.NewSource(2)))
	if err != nil {
		t.Fatal(err)
	}
	// Inertia must be non-increasing in k.
	for i := 1; i < len(inertias); i++ {
		if inertias[i] > inertias[i-1]+1e-6 {
			t.Errorf("inertia increased at k=%d: %v", ks[i], inertias)
		}
	}
	elbow := ElbowPoint(ks, inertias)
	if ks[elbow] < 2 || ks[elbow] > 4 {
		t.Errorf("elbow at k=%d, expected near 3", ks[elbow])
	}
}

func TestGapStatistic(t *testing.T) {
	rng := rand.New(rand.NewSource(5))
	data, _ := MakeBlobs(60, 3, nil, 0.4, rng)
	gap, se, err := GapStatistic(data, 3, 5, rand.New(rand.NewSource(6)))
	if err != nil {
		t.Fatal(err)
	}
	if se < 0 {
		t.Errorf("stderr negative: %v", se)
	}
	if math.IsNaN(gap) {
		t.Errorf("gap is NaN")
	}
}

func TestComparisonMetrics(t *testing.T) {
	a := []int{0, 0, 1, 1}
	b := []int{1, 1, 0, 0} // same partition, relabelled
	if !approx(RandIndex(a, b), 1, tol) {
		t.Errorf("RandIndex identical partition = %v, want 1", RandIndex(a, b))
	}
	if !approx(AdjustedRandIndex(a, b), 1, tol) {
		t.Errorf("ARI identical = %v, want 1", AdjustedRandIndex(a, b))
	}
	if !approx(NormalizedMutualInformation(a, b), 1, 1e-9) {
		t.Errorf("NMI identical = %v, want 1", NormalizedMutualInformation(a, b))
	}
	if !approx(FowlkesMallowsIndex(a, b), 1, tol) {
		t.Errorf("FMI identical = %v, want 1", FowlkesMallowsIndex(a, b))
	}
	if !approx(VMeasure(a, b), 1, 1e-9) {
		t.Errorf("VMeasure identical = %v, want 1", VMeasure(a, b))
	}
	if !approx(Purity(a, b), 1, tol) {
		t.Errorf("Purity identical = %v, want 1", Purity(a, b))
	}
	if !approx(ClusteringAccuracy(a, b), 1, tol) {
		t.Errorf("accuracy identical = %v, want 1", ClusteringAccuracy(a, b))
	}
	// A random-vs-truth partition: ARI near 0.
	c := []int{0, 1, 0, 1}
	if ari := AdjustedRandIndex(a, c); ari > 0.5 {
		t.Errorf("ARI for poor match too high: %v", ari)
	}
}

func TestContingencyMatrix(t *testing.T) {
	table, tl, pl := ContingencyMatrix([]int{0, 0, 1, 1}, []int{0, 0, 1, 1})
	if len(tl) != 2 || len(pl) != 2 {
		t.Fatalf("labels wrong")
	}
	if table[0][0] != 2 || table[1][1] != 2 || table[0][1] != 0 {
		t.Errorf("contingency table wrong: %v", table)
	}
}

func TestSpectralClustering(t *testing.T) {
	data := twoBlobs()
	labels, err := SpectralClustering(data, 2, 0.1, rand.New(rand.NewSource(1)))
	if err != nil {
		t.Fatal(err)
	}
	if labels[0] == labels[4] {
		t.Errorf("spectral did not separate blobs: %v", labels)
	}
}

func TestJacobiEigen(t *testing.T) {
	// Symmetric matrix with known eigenvalues 2 and 4 (from [[3,1],[1,3]]).
	a := [][]float64{{3, 1}, {1, 3}}
	vals, vecs, err := JacobiEigenSymmetric(a)
	if err != nil {
		t.Fatal(err)
	}
	sorted := append([]float64(nil), vals...)
	sort.Float64s(sorted)
	if !approx(sorted[0], 2, 1e-9) || !approx(sorted[1], 4, 1e-9) {
		t.Errorf("eigenvalues = %v, want [2 4]", sorted)
	}
	// Verify A v = lambda v for each column.
	for c := 0; c < 2; c++ {
		v := []float64{vecs[0][c], vecs[1][c]}
		av, _ := MatVecMul(a, v)
		for i := range v {
			if !approx(av[i], vals[c]*v[i], 1e-8) {
				t.Errorf("eigenvector %d failed: Av=%v lambda*v=%v", c, av, VectorScale(v, vals[c]))
			}
		}
	}
}

func TestBisectingKMeans(t *testing.T) {
	rng := rand.New(rand.NewSource(9))
	data, _ := MakeBlobs(60, 3, nil, 0.3, rng)
	res, err := BisectingKMeans(data, 3, rand.New(rand.NewSource(10)))
	if err != nil {
		t.Fatal(err)
	}
	if res.K != 3 {
		t.Errorf("expected 3 clusters, got %d", res.K)
	}
	if CountClusters(res.Labels) != 3 {
		t.Errorf("expected 3 distinct labels")
	}
}

func TestFuzzyCMeans(t *testing.T) {
	data := twoBlobs()
	res, err := FuzzyCMeans(data, 2, 2, 300, 1e-6, rand.New(rand.NewSource(1)))
	if err != nil {
		t.Fatal(err)
	}
	// Membership rows sum to 1.
	for i, row := range res.Membership {
		var s float64
		for _, v := range row {
			s += v
		}
		if !approx(s, 1, 1e-6) {
			t.Errorf("row %d membership sum = %v", i, s)
		}
	}
	if res.Labels[0] == res.Labels[4] {
		t.Errorf("fuzzy c-means did not separate blobs")
	}
	pc := PartitionCoefficient(res.Membership)
	if pc < 0.5 || pc > 1.0001 {
		t.Errorf("partition coefficient out of range: %v", pc)
	}
	if PartitionEntropy(res.Membership) < 0 {
		t.Errorf("partition entropy negative")
	}
}

func TestPreprocessing(t *testing.T) {
	data := [][]float64{{1, 10}, {2, 20}, {3, 30}}
	means := ColumnMeans(data)
	if !approx(means[0], 2, tol) || !approx(means[1], 20, tol) {
		t.Errorf("means = %v", means)
	}
	std := Standardize(data)
	sm := ColumnMeans(std)
	if !approx(sm[0], 0, 1e-9) || !approx(sm[1], 0, 1e-9) {
		t.Errorf("standardized means not zero: %v", sm)
	}
	scaled := MinMaxScale(data)
	if !approx(scaled[0][0], 0, tol) || !approx(scaled[2][0], 1, tol) {
		t.Errorf("minmax scale wrong: %v", scaled)
	}
	cov := CovarianceMatrix(data)
	if !approx(cov[0][0], 1, 1e-9) {
		t.Errorf("cov[0][0] = %v, want 1", cov[0][0])
	}
	if !approx(Quantile([]float64{1, 2, 3, 4}, 0.5), 2.5, 1e-9) {
		t.Errorf("median quantile wrong")
	}
}

func TestCentroidAndMedoid(t *testing.T) {
	pts := [][]float64{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	c := Centroid(pts)
	if !approx(c[0], 1, tol) || !approx(c[1], 1, tol) {
		t.Errorf("centroid = %v, want [1 1]", c)
	}
	m := Medoid(pts, nil, Euclidean)
	if m < 0 || m >= len(pts) {
		t.Errorf("medoid index invalid: %d", m)
	}
}

func TestKDistances(t *testing.T) {
	data := [][]float64{{0}, {1}, {2}, {10}}
	kd := KDistances(data, 1, Euclidean)
	// Nearest neighbour distances: 1,1,1,8.
	want := []float64{1, 1, 1, 8}
	for i := range want {
		if !approx(kd[i], want[i], tol) {
			t.Errorf("kd[%d] = %v, want %v", i, kd[i], want[i])
		}
	}
	sorted := SortedKDistances(data, 1, Euclidean)
	if sorted[len(sorted)-1] != 8 {
		t.Errorf("sorted k-distances max = %v, want 8", sorted[len(sorted)-1])
	}
}

func TestNearestNeighbors(t *testing.T) {
	data := [][]float64{{0, 0}, {1, 0}, {5, 5}}
	idx, d := NearestNeighbor([]float64{0.1, 0}, data, Euclidean)
	if idx != 0 || !approx(d, 0.1, 1e-9) {
		t.Errorf("nearest = %d dist %v", idx, d)
	}
	knn := KNearestNeighbors([]float64{0, 0}, data, 2, Euclidean)
	if knn[0] != 0 || knn[1] != 1 {
		t.Errorf("knn = %v, want [0 1]", knn)
	}
}

func TestGenerators(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	data, labels := MakeBlobs(30, 3, nil, 0.5, rng)
	if len(data) != 30 || len(labels) != 30 {
		t.Errorf("MakeBlobs wrong sizes")
	}
	md, ml := MakeMoons(40, 0.05, rng)
	if len(md) != 40 || len(ml) != 40 {
		t.Errorf("MakeMoons wrong sizes")
	}
	cd, cl := MakeCircles(40, 0.5, 0.02, rng)
	if len(cd) != 40 || len(cl) != 40 {
		t.Errorf("MakeCircles wrong sizes")
	}
}

// ExampleKMeans demonstrates clustering a small dataset into two groups.
func ExampleKMeans() {
	data := [][]float64{
		{1, 1}, {1.5, 2}, {1, 1.5},
		{8, 8}, {8.5, 8}, {8, 8.5},
	}
	res, err := KMeansWithOptions(data, 2, KMeansOptions{NInit: 10, Rand: rand.New(rand.NewSource(1))})
	if err != nil {
		panic(err)
	}
	// Relabel so the cluster containing point 0 is always cluster A.
	first := res.Labels[0]
	for i, l := range res.Labels {
		if l == first {
			fmt.Printf("point %d -> A\n", i)
		} else {
			fmt.Printf("point %d -> B\n", i)
		}
	}
	// Output:
	// point 0 -> A
	// point 1 -> A
	// point 2 -> A
	// point 3 -> B
	// point 4 -> B
	// point 5 -> B
}

// ExampleDBSCAN demonstrates density-based clustering with noise detection.
func ExampleDBSCAN() {
	data := [][]float64{
		{0, 0}, {0.1, 0.1}, {0, 0.1},
		{5, 5}, {5.1, 5}, {5, 5.1},
		{100, 100}, // noise
	}
	res, _ := DBSCAN(data, 0.5, 3, Euclidean)
	fmt.Println("clusters:", res.NClusters)
	fmt.Println("last point noise:", res.Labels[6] == NoiseLabel)
	// Output:
	// clusters: 2
	// last point noise: true
}

// ExampleSilhouetteScore shows scoring a clean two-cluster split.
func ExampleSilhouetteScore() {
	data := [][]float64{{0, 0}, {0, 1}, {10, 0}, {10, 1}}
	labels := []int{0, 0, 1, 1}
	score := SilhouetteScore(data, labels, Euclidean)
	fmt.Printf("%.3f\n", score)
	// Output:
	// 0.900
}
