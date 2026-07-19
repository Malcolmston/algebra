package clustering

import (
	"math"
	"math/rand"
)

// MakeBlobs generates n samples drawn from k isotropic Gaussian blobs whose
// centres are placed at the given locations (one row per centre). std sets the
// standard deviation of each blob. It returns the samples and their true blob
// labels. If centers is nil, k centres are placed on a circle.
func MakeBlobs(n, k int, centers [][]float64, std float64, rng *rand.Rand) ([][]float64, []int) {
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	if std <= 0 {
		std = 1
	}
	if centers == nil {
		centers = make([][]float64, k)
		for i := 0; i < k; i++ {
			angle := 2 * math.Pi * float64(i) / float64(k)
			centers[i] = []float64{10 * math.Cos(angle), 10 * math.Sin(angle)}
		}
	}
	dim := len(centers[0])
	data := make([][]float64, n)
	labels := make([]int, n)
	for i := 0; i < n; i++ {
		c := i % len(centers)
		labels[i] = c
		point := make([]float64, dim)
		for d := 0; d < dim; d++ {
			point[d] = centers[c][d] + rng.NormFloat64()*std
		}
		data[i] = point
	}
	return data, labels
}

// MakeMoons generates the classic two interleaving half-circles ("moons")
// dataset with n samples total and the given Gaussian noise level. It returns
// the samples and their true class labels (0 or 1).
func MakeMoons(n int, noise float64, rng *rand.Rand) ([][]float64, []int) {
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	nOut := n / 2
	nIn := n - nOut
	data := make([][]float64, 0, n)
	labels := make([]int, 0, n)
	for i := 0; i < nOut; i++ {
		t := math.Pi * float64(i) / float64(nOut)
		x := math.Cos(t) + rng.NormFloat64()*noise
		y := math.Sin(t) + rng.NormFloat64()*noise
		data = append(data, []float64{x, y})
		labels = append(labels, 0)
	}
	for i := 0; i < nIn; i++ {
		t := math.Pi * float64(i) / float64(nIn)
		x := 1 - math.Cos(t) + rng.NormFloat64()*noise
		y := 0.5 - math.Sin(t) + rng.NormFloat64()*noise
		data = append(data, []float64{x, y})
		labels = append(labels, 1)
	}
	return data, labels
}

// MakeCircles generates two concentric circles with n samples total, an inner
// circle of the given factor (0 < factor < 1) of the outer radius, and Gaussian
// noise. It returns the samples and their true class labels (0 = outer,
// 1 = inner).
func MakeCircles(n int, factor, noise float64, rng *rand.Rand) ([][]float64, []int) {
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	if factor <= 0 || factor >= 1 {
		factor = 0.5
	}
	nOut := n / 2
	nIn := n - nOut
	data := make([][]float64, 0, n)
	labels := make([]int, 0, n)
	for i := 0; i < nOut; i++ {
		t := 2 * math.Pi * float64(i) / float64(nOut)
		data = append(data, []float64{
			math.Cos(t) + rng.NormFloat64()*noise,
			math.Sin(t) + rng.NormFloat64()*noise,
		})
		labels = append(labels, 0)
	}
	for i := 0; i < nIn; i++ {
		t := 2 * math.Pi * float64(i) / float64(nIn)
		data = append(data, []float64{
			factor*math.Cos(t) + rng.NormFloat64()*noise,
			factor*math.Sin(t) + rng.NormFloat64()*noise,
		})
		labels = append(labels, 1)
	}
	return data, labels
}
