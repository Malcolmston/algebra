package stats

import (
	"math"
	"testing"
)

// statsMRApproxEqual reports whether a and b are within tol, treating two NaNs
// as equal so NaN sentinels can be checked in table tests.
func statsMRApproxEqual(a, b, tol float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return math.IsNaN(a) && math.IsNaN(b)
	}
	return math.Abs(a-b) <= tol
}

const statsMRTol = 1e-9

func TestMultipleLinearRegression(t *testing.T) {
	tests := []struct {
		name      string
		X         [][]float64
		y         []float64
		intercept bool
		want      []float64 // nil means expect nil coeffs + NaN r2
		wantR2    float64
	}{
		{
			name:      "exact fit with intercept y=1+2x1+3x2",
			X:         [][]float64{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 1}},
			y:         []float64{1, 3, 4, 6, 8},
			intercept: true,
			want:      []float64{1, 2, 3},
			wantR2:    1,
		},
		{
			name:      "exact fit no intercept y=2x",
			X:         [][]float64{{1}, {2}, {3}},
			y:         []float64{2, 4, 6},
			intercept: false,
			want:      []float64{2},
			wantR2:    1,
		},
		{
			name:      "intercept-only constant response",
			X:         [][]float64{{}, {}, {}},
			y:         []float64{5, 5, 5},
			intercept: true,
			want:      []float64{5},
			wantR2:    1,
		},
		{
			name:      "dimension mismatch rows vs y",
			X:         [][]float64{{1}, {2}, {3}},
			y:         []float64{1, 2},
			intercept: true,
			want:      nil,
			wantR2:    math.NaN(),
		},
		{
			name:      "ragged rows",
			X:         [][]float64{{1, 2}, {3}},
			y:         []float64{1, 2},
			intercept: false,
			want:      nil,
			wantR2:    math.NaN(),
		},
		{
			name:      "too few rows",
			X:         [][]float64{{1, 2}},
			y:         []float64{1},
			intercept: false,
			want:      nil,
			wantR2:    math.NaN(),
		},
		{
			name:      "singular collinear predictors",
			X:         [][]float64{{1, 1}, {2, 2}, {3, 3}},
			y:         []float64{1, 2, 3},
			intercept: false,
			want:      nil,
			wantR2:    math.NaN(),
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			coeffs, r2 := MultipleLinearRegression(tc.X, tc.y, tc.intercept)
			if tc.want == nil {
				if coeffs != nil {
					t.Fatalf("expected nil coeffs, got %v", coeffs)
				}
				if !math.IsNaN(r2) {
					t.Fatalf("expected NaN r2, got %v", r2)
				}
				return
			}
			if len(coeffs) != len(tc.want) {
				t.Fatalf("coeffs length = %d, want %d (%v)", len(coeffs), len(tc.want), coeffs)
			}
			for i := range tc.want {
				if !statsMRApproxEqual(coeffs[i], tc.want[i], statsMRTol) {
					t.Errorf("coeffs[%d] = %v, want %v", i, coeffs[i], tc.want[i])
				}
			}
			if !statsMRApproxEqual(r2, tc.wantR2, statsMRTol) {
				t.Errorf("r2 = %v, want %v", r2, tc.wantR2)
			}
		})
	}
}

func TestMultipleLinearRegressionR2Partial(t *testing.T) {
	// Data with noise so R2 is strictly between 0 and 1.
	X := [][]float64{{1}, {2}, {3}, {4}, {5}}
	y := []float64{1, 2, 1.3, 3.75, 2.25}
	coeffs, r2 := MultipleLinearRegression(X, y, true)
	if coeffs == nil {
		t.Fatal("unexpected nil coeffs")
	}
	if r2 <= 0 || r2 >= 1 {
		t.Errorf("r2 = %v, want in (0,1)", r2)
	}
	// Cross-check the slope/intercept against simple LinearRegression.
	slope, intercept, r2b := LinearRegression([]float64{1, 2, 3, 4, 5}, y)
	if !statsMRApproxEqual(coeffs[0], intercept, 1e-9) {
		t.Errorf("intercept = %v, want %v", coeffs[0], intercept)
	}
	if !statsMRApproxEqual(coeffs[1], slope, 1e-9) {
		t.Errorf("slope = %v, want %v", coeffs[1], slope)
	}
	if !statsMRApproxEqual(r2, r2b, 1e-9) {
		t.Errorf("r2 = %v, want %v", r2, r2b)
	}
}

func TestRidgeRegression(t *testing.T) {
	// Single predictor, no intercept: XtX=14, Xty=28.
	X := [][]float64{{1}, {2}, {3}}
	y := []float64{2, 4, 6}

	// lambda=0 reduces to OLS.
	b0 := RidgeRegression(X, y, 0, false)
	if len(b0) != 1 || !statsMRApproxEqual(b0[0], 2, statsMRTol) {
		t.Fatalf("lambda=0 coeffs = %v, want [2]", b0)
	}

	// lambda=1: (14+1)b = 28 => b = 28/15.
	b1 := RidgeRegression(X, y, 1, false)
	if len(b1) != 1 || !statsMRApproxEqual(b1[0], 28.0/15.0, statsMRTol) {
		t.Fatalf("lambda=1 coeffs = %v, want [%v]", b1, 28.0/15.0)
	}

	// Larger lambda shrinks the coefficient toward zero.
	b2 := RidgeRegression(X, y, 10, false)
	if !(math.Abs(b2[0]) < math.Abs(b1[0])) {
		t.Errorf("expected more shrinkage: |%v| should be < |%v|", b2[0], b1[0])
	}

	// Negative lambda is rejected.
	if got := RidgeRegression(X, y, -1, false); got != nil {
		t.Errorf("negative lambda = %v, want nil", got)
	}
	// NaN lambda is rejected.
	if got := RidgeRegression(X, y, math.NaN(), false); got != nil {
		t.Errorf("NaN lambda = %v, want nil", got)
	}
	// Dimension mismatch is rejected.
	if got := RidgeRegression([][]float64{{1}, {2}}, []float64{1}, 1, false); got != nil {
		t.Errorf("mismatch = %v, want nil", got)
	}
}

func TestRidgeRegressionEqualsOLS(t *testing.T) {
	X := [][]float64{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 1}, {1, 2}}
	y := []float64{1.2, 3.1, 3.9, 6.2, 7.8, 8.1}
	ols, _ := MultipleLinearRegression(X, y, true)
	ridge := RidgeRegression(X, y, 0, true)
	if len(ols) != len(ridge) {
		t.Fatalf("length mismatch: %d vs %d", len(ols), len(ridge))
	}
	for i := range ols {
		if !statsMRApproxEqual(ols[i], ridge[i], 1e-9) {
			t.Errorf("coeff[%d]: ridge(0)=%v, ols=%v", i, ridge[i], ols[i])
		}
	}
}

func TestRidgeInterceptNotPenalized(t *testing.T) {
	// If the intercept were penalized, a large lambda would pull it toward 0.
	// With an unpenalized intercept and heavy penalty on the slope, the fit
	// collapses toward a horizontal line at the mean of y.
	X := [][]float64{{1}, {2}, {3}, {4}}
	y := []float64{10, 12, 14, 16}
	b := RidgeRegression(X, y, 1e6, true)
	if b == nil {
		t.Fatal("unexpected nil")
	}
	meanY := Mean(y)
	if !statsMRApproxEqual(b[0], meanY, 1e-2) {
		t.Errorf("intercept = %v, want approx mean %v", b[0], meanY)
	}
	if math.Abs(b[1]) > 1e-2 {
		t.Errorf("slope = %v, want near 0 under heavy penalty", b[1])
	}
}

func TestPredict(t *testing.T) {
	tests := []struct {
		name      string
		coeffs    []float64
		x         []float64
		intercept bool
		want      float64
	}{
		{"with intercept", []float64{1, 2, 3}, []float64{4, 5}, true, 1 + 2*4 + 3*5},
		{"no intercept", []float64{2, 3}, []float64{4, 5}, false, 2*4 + 3*5},
		{"intercept only", []float64{7}, []float64{}, true, 7},
		{"wrong length with intercept", []float64{1, 2, 3}, []float64{4}, true, math.NaN()},
		{"wrong length no intercept", []float64{2, 3}, []float64{4, 5, 6}, false, math.NaN()},
		{"empty coeffs with intercept", []float64{}, []float64{}, true, math.NaN()},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := Predict(tc.coeffs, tc.x, tc.intercept)
			if !statsMRApproxEqual(got, tc.want, statsMRTol) {
				t.Errorf("Predict = %v, want %v", got, tc.want)
			}
		})
	}
}

func TestPredictRoundTrip(t *testing.T) {
	X := [][]float64{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 1}}
	y := []float64{1, 3, 4, 6, 8}
	coeffs, _ := MultipleLinearRegression(X, y, true)
	for i, row := range X {
		got := Predict(coeffs, row, true)
		if !statsMRApproxEqual(got, y[i], 1e-9) {
			t.Errorf("row %d: Predict = %v, want %v", i, got, y[i])
		}
	}
}

func TestCovarianceMatrix(t *testing.T) {
	x := []float64{1, 2, 3, 4}
	y := []float64{2, 4, 6, 8}
	cov := CovarianceMatrix([][]float64{x, y})
	want := [][]float64{
		{5.0 / 3.0, 10.0 / 3.0},
		{10.0 / 3.0, 20.0 / 3.0},
	}
	for i := range want {
		for j := range want[i] {
			if !statsMRApproxEqual(cov[i][j], want[i][j], 1e-9) {
				t.Errorf("cov[%d][%d] = %v, want %v", i, j, cov[i][j], want[i][j])
			}
		}
	}
	// Diagonal must equal stats.Variance and matrix must be symmetric.
	if !statsMRApproxEqual(cov[0][0], Variance(x), 1e-12) {
		t.Errorf("diag[0] = %v, want Variance = %v", cov[0][0], Variance(x))
	}
	if !statsMRApproxEqual(cov[1][1], Variance(y), 1e-12) {
		t.Errorf("diag[1] = %v, want Variance = %v", cov[1][1], Variance(y))
	}
	if cov[0][1] != cov[1][0] {
		t.Errorf("not symmetric: %v vs %v", cov[0][1], cov[1][0])
	}
	// Off-diagonal must equal stats.Covariance.
	if !statsMRApproxEqual(cov[0][1], Covariance(x, y), 1e-12) {
		t.Errorf("off-diag = %v, want Covariance = %v", cov[0][1], Covariance(x, y))
	}
}

func TestCovarianceMatrixInvalid(t *testing.T) {
	// No columns -> nil.
	if got := CovarianceMatrix(nil); got != nil {
		t.Errorf("nil input = %v, want nil", got)
	}
	// Length mismatch -> all NaN, shape preserved.
	mm := CovarianceMatrix([][]float64{{1, 2, 3}, {1, 2}})
	if len(mm) != 2 {
		t.Fatalf("shape = %d, want 2", len(mm))
	}
	for i := range mm {
		for j := range mm[i] {
			if !math.IsNaN(mm[i][j]) {
				t.Errorf("mm[%d][%d] = %v, want NaN", i, j, mm[i][j])
			}
		}
	}
	// n < 2 -> all NaN.
	small := CovarianceMatrix([][]float64{{1}, {2}})
	for i := range small {
		for j := range small[i] {
			if !math.IsNaN(small[i][j]) {
				t.Errorf("small[%d][%d] = %v, want NaN", i, j, small[i][j])
			}
		}
	}
}

func TestCorrelationMatrix(t *testing.T) {
	x := []float64{1, 2, 3, 4}
	y := []float64{2, 4, 6, 8} // perfectly correlated with x
	z := []float64{4, 3, 2, 1} // perfectly anti-correlated
	corr := CorrelationMatrix([][]float64{x, y, z})
	// Unit diagonal.
	for i := 0; i < 3; i++ {
		if !statsMRApproxEqual(corr[i][i], 1, 1e-12) {
			t.Errorf("diag[%d] = %v, want 1", i, corr[i][i])
		}
	}
	if !statsMRApproxEqual(corr[0][1], 1, 1e-9) {
		t.Errorf("corr(x,y) = %v, want 1", corr[0][1])
	}
	if !statsMRApproxEqual(corr[0][2], -1, 1e-9) {
		t.Errorf("corr(x,z) = %v, want -1", corr[0][2])
	}
	// Off-diagonal must equal stats.Correlation and be symmetric.
	if !statsMRApproxEqual(corr[0][1], Correlation(x, y), 1e-12) {
		t.Errorf("off-diag = %v, want Correlation = %v", corr[0][1], Correlation(x, y))
	}
	if corr[0][2] != corr[2][0] {
		t.Errorf("not symmetric: %v vs %v", corr[0][2], corr[2][0])
	}
}

func TestCorrelationMatrixInvalid(t *testing.T) {
	if got := CorrelationMatrix(nil); got != nil {
		t.Errorf("nil input = %v, want nil", got)
	}
	mm := CorrelationMatrix([][]float64{{1, 2, 3}, {1, 2}})
	if len(mm) != 2 {
		t.Fatalf("shape = %d, want 2", len(mm))
	}
	for i := range mm {
		for j := range mm[i] {
			if !math.IsNaN(mm[i][j]) {
				t.Errorf("mm[%d][%d] = %v, want NaN", i, j, mm[i][j])
			}
		}
	}
}

func TestSolveDense(t *testing.T) {
	// 2x + y = 5 ; x + 3y = 10  => x = 1, y = 3.
	a := [][]float64{{2, 1}, {1, 3}}
	b := []float64{5, 10}
	x, ok := statsSolveDense(a, b)
	if !ok {
		t.Fatal("unexpected singular")
	}
	if !statsMRApproxEqual(x[0], 1, 1e-9) || !statsMRApproxEqual(x[1], 3, 1e-9) {
		t.Errorf("solution = %v, want [1 3]", x)
	}
	// Singular system.
	sa := [][]float64{{1, 2}, {2, 4}}
	sb := []float64{3, 6}
	if _, ok := statsSolveDense(sa, sb); ok {
		t.Error("expected singular, got ok")
	}
}

// statsMRBenchData builds a deterministic n-by-p regression problem.
func statsMRBenchData(n, p int) ([][]float64, []float64) {
	X := make([][]float64, n)
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		row := make([]float64, p)
		sum := 1.0
		for j := 0; j < p; j++ {
			v := float64((i*7+j*13)%101) / 101.0
			row[j] = v
			sum += float64(j+1) * v
		}
		X[i] = row
		y[i] = sum
	}
	return X, y
}

func BenchmarkMultipleLinearRegression(b *testing.B) {
	X, y := statsMRBenchData(200, 8)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		coeffs, _ := MultipleLinearRegression(X, y, true)
		if coeffs == nil {
			b.Fatal("unexpected nil")
		}
	}
}

func BenchmarkRidgeRegression(b *testing.B) {
	X, y := statsMRBenchData(200, 8)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		coeffs := RidgeRegression(X, y, 0.5, true)
		if coeffs == nil {
			b.Fatal("unexpected nil")
		}
	}
}
