package clustering

import (
	"math"
	"math/rand"
	"sort"
)

// JacobiEigenSymmetric computes all eigenvalues and eigenvectors of a real
// symmetric matrix using the cyclic Jacobi rotation method. It returns the
// eigenvalues and a matrix whose columns are the corresponding orthonormal
// eigenvectors. The input matrix is not modified.
func JacobiEigenSymmetric(a [][]float64) (eigenvalues []float64, eigenvectors [][]float64, err error) {
	n := len(a)
	if n == 0 || len(a[0]) != n {
		return nil, nil, ErrDimensionMismatch
	}
	m := CloneMatrix(a)
	v := Identity(n)
	const maxSweeps = 100
	for sweep := 0; sweep < maxSweeps; sweep++ {
		off := offDiagonalNorm(m)
		if off < 1e-14 {
			break
		}
		for p := 0; p < n-1; p++ {
			for q := p + 1; q < n; q++ {
				if math.Abs(m[p][q]) < 1e-300 {
					continue
				}
				theta := (m[q][q] - m[p][p]) / (2 * m[p][q])
				t := sign(theta) / (math.Abs(theta) + math.Sqrt(theta*theta+1))
				if theta == 0 {
					t = 1
				}
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				jacobiRotate(m, v, p, q, c, s, n)
			}
		}
	}
	eigenvalues = make([]float64, n)
	for i := 0; i < n; i++ {
		eigenvalues[i] = m[i][i]
	}
	return eigenvalues, v, nil
}

func offDiagonalNorm(m [][]float64) float64 {
	var s float64
	for i := range m {
		for j := i + 1; j < len(m); j++ {
			s += 2 * m[i][j] * m[i][j]
		}
	}
	return math.Sqrt(s)
}

func sign(x float64) float64 {
	if x < 0 {
		return -1
	}
	return 1
}

func jacobiRotate(m, v [][]float64, p, q int, c, s float64, n int) {
	for i := 0; i < n; i++ {
		mip := m[i][p]
		miq := m[i][q]
		m[i][p] = c*mip - s*miq
		m[i][q] = s*mip + c*miq
	}
	for i := 0; i < n; i++ {
		mpi := m[p][i]
		mqi := m[q][i]
		m[p][i] = c*mpi - s*mqi
		m[q][i] = s*mpi + c*mqi
	}
	for i := 0; i < n; i++ {
		vip := v[i][p]
		viq := v[i][q]
		v[i][p] = c*vip - s*viq
		v[i][q] = s*vip + c*viq
	}
}

// SortEigen returns the eigenvalues and eigenvectors reordered by ascending
// eigenvalue. eigenvectors is a matrix whose columns are eigenvectors; the
// returned matrix preserves the column-eigenvector convention.
func SortEigen(eigenvalues []float64, eigenvectors [][]float64) ([]float64, [][]float64) {
	n := len(eigenvalues)
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(a, b int) bool { return eigenvalues[idx[a]] < eigenvalues[idx[b]] })
	vals := make([]float64, n)
	vecs := make([][]float64, len(eigenvectors))
	for i := range vecs {
		vecs[i] = make([]float64, n)
	}
	for newCol, oldCol := range idx {
		vals[newCol] = eigenvalues[oldCol]
		for r := range eigenvectors {
			vecs[r][newCol] = eigenvectors[r][oldCol]
		}
	}
	return vals, vecs
}

// RBFAffinityMatrix returns the Gaussian (RBF) affinity matrix of data, where
// entry [i][j] = exp(-gamma * ||x_i - x_j||^2). Diagonal entries are 1.
func RBFAffinityMatrix(data [][]float64, gamma float64) [][]float64 {
	n := len(data)
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		a[i][i] = 1
		for j := i + 1; j < n; j++ {
			w := math.Exp(-gamma * SquaredEuclidean(data[i], data[j]))
			a[i][j] = w
			a[j][i] = w
		}
	}
	return a
}

// KNearestAffinityMatrix returns a symmetric k-nearest-neighbour affinity
// matrix: entry [i][j] is 1 if either i is among the k nearest neighbours of j
// or vice versa, and 0 otherwise. metric may be nil (Euclidean).
func KNearestAffinityMatrix(data [][]float64, k int, metric Metric) [][]float64 {
	if metric == nil {
		metric = Euclidean
	}
	n := len(data)
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		nbrs := KNearestNeighbors(data[i], data, k+1, metric)
		for _, j := range nbrs {
			if j == i {
				continue
			}
			a[i][j] = 1
			a[j][i] = 1
		}
	}
	return a
}

// DegreeMatrix returns the diagonal degree matrix of an affinity matrix, where
// each diagonal entry is the row sum of affinities.
func DegreeMatrix(affinity [][]float64) [][]float64 {
	n := len(affinity)
	d := Zeros(n, n)
	for i := 0; i < n; i++ {
		var s float64
		for j := 0; j < n; j++ {
			s += affinity[i][j]
		}
		d[i][i] = s
	}
	return d
}

// UnnormalizedLaplacian returns the graph Laplacian L = D - W for the affinity
// matrix W.
func UnnormalizedLaplacian(affinity [][]float64) [][]float64 {
	n := len(affinity)
	l := make([][]float64, n)
	for i := 0; i < n; i++ {
		l[i] = make([]float64, n)
		var deg float64
		for j := 0; j < n; j++ {
			deg += affinity[i][j]
		}
		for j := 0; j < n; j++ {
			if i == j {
				l[i][j] = deg - affinity[i][j]
			} else {
				l[i][j] = -affinity[i][j]
			}
		}
	}
	return l
}

// NormalizedLaplacian returns the symmetric normalised graph Laplacian
// L_sym = I - D^{-1/2} W D^{-1/2} for the affinity matrix W.
func NormalizedLaplacian(affinity [][]float64) [][]float64 {
	n := len(affinity)
	deg := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			deg[i] += affinity[i][j]
		}
	}
	dinv := make([]float64, n)
	for i := 0; i < n; i++ {
		if deg[i] > 0 {
			dinv[i] = 1 / math.Sqrt(deg[i])
		}
	}
	l := make([][]float64, n)
	for i := 0; i < n; i++ {
		l[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			var w float64
			if i == j {
				w = 1
			}
			l[i][j] = w - dinv[i]*affinity[i][j]*dinv[j]
		}
	}
	return l
}

// SpectralClustering partitions data into k clusters using the normalised
// spectral clustering algorithm of Ng, Jordan and Weiss: it builds an RBF
// affinity graph, forms the normalised Laplacian, embeds points using its k
// smallest eigenvectors, row-normalises the embedding and runs k-means on it.
func SpectralClustering(data [][]float64, k int, gamma float64, rng *rand.Rand) ([]int, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > len(data) {
		return nil, ErrInvalidK
	}
	affinity := RBFAffinityMatrix(data, gamma)
	return SpectralClusteringFromAffinity(affinity, k, rng)
}

// SpectralClusteringFromAffinity runs normalised spectral clustering starting
// from a precomputed affinity matrix.
func SpectralClusteringFromAffinity(affinity [][]float64, k int, rng *rand.Rand) ([]int, error) {
	n := len(affinity)
	if n == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > n {
		return nil, ErrInvalidK
	}
	lap := NormalizedLaplacian(affinity)
	vals, vecs, err := JacobiEigenSymmetric(lap)
	if err != nil {
		return nil, err
	}
	vals, vecs = SortEigen(vals, vecs)
	// Build embedding from the k smallest eigenvectors (columns 0..k-1).
	embed := make([][]float64, n)
	for i := 0; i < n; i++ {
		embed[i] = make([]float64, k)
		var norm float64
		for c := 0; c < k; c++ {
			embed[i][c] = vecs[i][c]
			norm += vecs[i][c] * vecs[i][c]
		}
		norm = math.Sqrt(norm)
		if norm > 0 {
			for c := 0; c < k; c++ {
				embed[i][c] /= norm
			}
		}
	}
	res, err := KMeansWithOptions(embed, k, KMeansOptions{NInit: 10, Rand: rng})
	if err != nil {
		return nil, err
	}
	return res.Labels, nil
}
