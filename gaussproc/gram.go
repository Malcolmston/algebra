package gaussproc

// GramMatrix returns the symmetric n-by-n Gram matrix K with
// K[i][j] = kernel.Eval(x[i], x[j]) for the data set x (a slice of vectors).
// Only the upper triangle is evaluated and mirrored, so the result is exactly
// symmetric.
func GramMatrix(kernel Kernel, x [][]float64) Matrix {
	n := len(x)
	k := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			v := kernel.Eval(x[i], x[j])
			k[i][j] = v
			k[j][i] = v
		}
	}
	return k
}

// CrossGramMatrix returns the len(x1)-by-len(x2) matrix K with
// K[i][j] = kernel.Eval(x1[i], x2[j]).
func CrossGramMatrix(kernel Kernel, x1, x2 [][]float64) Matrix {
	n, m := len(x1), len(x2)
	k := NewMatrix(n, m)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			k[i][j] = kernel.Eval(x1[i], x2[j])
		}
	}
	return k
}

// KernelDiagonal returns the vector of self-covariances kernel.Eval(x[i], x[i]).
func KernelDiagonal(kernel Kernel, x [][]float64) []float64 {
	d := make([]float64, len(x))
	for i := range x {
		d[i] = kernel.Eval(x[i], x[i])
	}
	return d
}

// NoisyGramMatrix returns GramMatrix(kernel, x) with the scalar noiseVar added
// to every diagonal entry, that is K + noiseVar·I. This is the covariance of
// noisy observations under a Gaussian process.
func NoisyGramMatrix(kernel Kernel, x [][]float64, noiseVar float64) Matrix {
	k := GramMatrix(kernel, x)
	for i := range k {
		k[i][i] += noiseVar
	}
	return k
}

// AddJitter returns a copy of square matrix k with the small positive value
// jitter added to each diagonal entry, improving numerical conditioning for the
// Cholesky factorisation.
func AddJitter(k Matrix, jitter float64) Matrix {
	return AddToDiagonal(k, jitter)
}

// IsKernelPSD reports whether the Gram matrix of kernel over the data set x is
// numerically positive definite, tested by attempting a Cholesky factorisation
// after adding the given jitter.
func IsKernelPSD(kernel Kernel, x [][]float64, jitter float64) bool {
	k := AddJitter(GramMatrix(kernel, x), jitter)
	_, err := Cholesky(k)
	return err == nil
}
