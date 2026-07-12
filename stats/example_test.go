package stats_test

import (
	"fmt"

	"github.com/malcolmston/algebra/stats"
)

func ExampleMean() {
	xs := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	fmt.Printf("mean=%.2f stddev=%.4f\n", stats.Mean(xs), stats.StdDev(xs))
	// Output: mean=5.00 stddev=2.1381
}

func ExampleMedian() {
	fmt.Printf("odd=%.1f even=%.1f\n",
		stats.Median([]float64{1, 2, 3}),
		stats.Median([]float64{1, 2, 3, 4}))
	// Output: odd=2.0 even=2.5
}

func ExampleQuantile() {
	xs := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	fmt.Printf("Q1=%.2f Q3=%.2f IQR=%.2f\n",
		stats.Quantile(xs, 0.25), stats.Quantile(xs, 0.75), stats.IQR(xs))
	// Output: Q1=4.00 Q3=5.50 IQR=1.50
}

func ExampleCorrelation() {
	xs := []float64{1, 2, 3, 4, 5}
	ys := []float64{2, 4, 6, 8, 10}
	fmt.Printf("r=%.4f\n", stats.Correlation(xs, ys))
	// Output: r=1.0000
}

func ExampleChoose() {
	fmt.Printf("C(52,5)=%.0f P(5,2)=%.0f 5!=%.0f\n",
		stats.Choose(52, 5), stats.Perm(5, 2), stats.Factorial(5))
	// Output: C(52,5)=2598960 P(5,2)=20 5!=120
}

func ExampleNormal() {
	n := stats.Normal{Mu: 0, Sigma: 1}
	fmt.Printf("CDF(0)=%.4f CDF(1.96)=%.4f Q(0.975)=%.4f\n",
		n.CDF(0), n.CDF(1.96), n.Quantile(0.975))
	// Output: CDF(0)=0.5000 CDF(1.96)=0.9750 Q(0.975)=1.9600
}

func ExampleBinomial() {
	b := stats.Binomial{N: 5, P: 0.5}
	fmt.Printf("PMF(2)=%.4f mean=%.1f var=%.2f\n",
		b.PMF(2), b.Mean(), b.Variance())
	// Output: PMF(2)=0.3125 mean=2.5 var=1.25
}

func ExamplePoisson() {
	p := stats.Poisson{Lambda: 3}
	fmt.Printf("PMF(2)=%.4f CDF(2)=%.4f\n", p.PMF(2), p.CDF(2))
	// Output: PMF(2)=0.2240 CDF(2)=0.4232
}

func ExampleLinearRegression() {
	xs := []float64{1, 2, 3, 4, 5}
	ys := []float64{3, 5, 7, 9, 11}
	slope, intercept, r2 := stats.LinearRegression(xs, ys)
	fmt.Printf("y = %.1fx + %.1f (R2=%.2f)\n", slope, intercept, r2)
	// Output: y = 2.0x + 1.0 (R2=1.00)
}
