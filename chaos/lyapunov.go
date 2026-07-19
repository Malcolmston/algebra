package chaos

import (
	"math"
	"sort"
)

// LyapunovLog1D estimates the Lyapunov exponent of the one-dimensional map f
// from its analytic derivative df. After discarding transient iterates it
// averages log|df(x)| over n steps. For the tent map of slope mu the result is
// log(2*mu) and for the logistic map at r=4 it approaches log(2).
func LyapunovLog1D(f, df Map1D, x0 float64, transient, n int) float64 {
	x := x0
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	if n <= 0 {
		return 0
	}
	var sum float64
	for i := 0; i < n; i++ {
		d := math.Abs(df(x))
		if d < 1e-300 {
			d = 1e-300
		}
		sum += math.Log(d)
		x = f(x)
	}
	return sum / float64(n)
}

// Lyapunov1D estimates the Lyapunov exponent of the map f using a numerical
// derivative, so no analytic derivative is required.
func Lyapunov1D(f Map1D, x0 float64, transient, n int) float64 {
	df := func(x float64) float64 { return Multiplier1D(f, x) }
	return LyapunovLog1D(f, df, x0, transient, n)
}

// LyapunovSeparation1D estimates the Lyapunov exponent of the map f by the
// trajectory-separation method: two orbits start a distance d0 apart and the
// perturbation is rescaled after every step. It needs no derivative.
func LyapunovSeparation1D(f Map1D, x0, d0 float64, transient, n int) float64 {
	x := x0
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	y := x + d0
	var sum float64
	for i := 0; i < n; i++ {
		x = f(x)
		y = f(y)
		d := math.Abs(y - x)
		if d < 1e-300 {
			d = 1e-300
		}
		sum += math.Log(d / d0)
		// Rescale the perturbed orbit back to distance d0.
		y = x + d0*(y-x)/d
	}
	return sum / float64(n)
}

// BenettinLargestMap estimates the largest Lyapunov exponent of the map F using
// the Benettin single-perturbation method: a tangent vector is evolved by the
// numerical Jacobian and renormalised each step, and the average log growth is
// returned. The result is per iteration.
func BenettinLargestMap(F MapN, x0 Vec, transient, n int) float64 {
	x := x0.Clone()
	for i := 0; i < transient; i++ {
		x = F(x)
	}
	v := make(Vec, len(x))
	v[0] = 1
	v = v.Normalize()
	var sum float64
	for i := 0; i < n; i++ {
		J := JacobianMap(F, x, 1e-7)
		v = J.MulVec(v)
		g := v.Norm()
		if g < 1e-300 {
			g = 1e-300
		}
		sum += math.Log(g)
		v = v.Scale(1 / g)
		x = F(x)
	}
	return sum / float64(n)
}

// BenettinLargestFlow estimates the largest Lyapunov exponent of the flow with
// field f by evolving a tangent vector along the trajectory with the
// variational equation, renormalising every step. The result is per unit time.
func BenettinLargestFlow(f Field, x0 Vec, h float64, transient, n int) float64 {
	x := x0.Clone()
	for i := 0; i < transient; i++ {
		x = StepRK4(f, x, h)
	}
	v := make(Vec, len(x))
	v[0] = 1
	v = v.Normalize()
	var sum float64
	for i := 0; i < n; i++ {
		// Evolve state and tangent vector together over one step.
		x, v = stepFlowTangent(f, x, v, h)
		g := v.Norm()
		if g < 1e-300 {
			g = 1e-300
		}
		sum += math.Log(g)
		v = v.Scale(1 / g)
	}
	return sum / (float64(n) * h)
}

// stepFlowTangent advances the state x and a single tangent vector v by one
// RK4 step using the numerical Jacobian for the variational part.
func stepFlowTangent(f Field, x, v Vec, h float64) (Vec, Vec) {
	jv := func(state, w Vec) Vec {
		return JacobianField(f, state, 1e-7).MulVec(w)
	}
	k1x := f(x)
	k1v := jv(x, v)
	x2 := x.AddScaled(h/2, k1x)
	v2 := v.AddScaled(h/2, k1v)
	k2x := f(x2)
	k2v := jv(x2, v2)
	x3 := x.AddScaled(h/2, k2x)
	v3 := v.AddScaled(h/2, k2v)
	k3x := f(x3)
	k3v := jv(x3, v3)
	x4 := x.AddScaled(h, k3x)
	v4 := v.AddScaled(h, k3v)
	k4x := f(x4)
	k4v := jv(x4, v4)
	nx := x.Clone()
	nv := v.Clone()
	for i := range nx {
		nx[i] += h / 6 * (k1x[i] + 2*k2x[i] + 2*k3x[i] + k4x[i])
		nv[i] += h / 6 * (k1v[i] + 2*k2v[i] + 2*k3v[i] + k4v[i])
	}
	return nx, nv
}

// LyapunovSpectrumMap estimates the full Lyapunov spectrum of the map F by the
// QR (Benettin) method: an orthonormal frame is evolved by the Jacobian and
// reorthonormalised each step, and the log of the diagonal of R is averaged.
// The exponents are returned in descending order, per iteration.
func LyapunovSpectrumMap(F MapN, x0 Vec, transient, n int) Vec {
	d := len(x0)
	x := x0.Clone()
	for i := 0; i < transient; i++ {
		x = F(x)
	}
	Q := Eye(d)
	sums := make(Vec, d)
	for i := 0; i < n; i++ {
		J := JacobianMap(F, x, 1e-7)
		Z := J.Mul(Q)
		Qn, R := QR(Z)
		for k := 0; k < d; k++ {
			r := math.Abs(R[k][k])
			if r < 1e-300 {
				r = 1e-300
			}
			sums[k] += math.Log(r)
		}
		Q = Qn
		x = F(x)
	}
	for k := range sums {
		sums[k] /= float64(n)
	}
	sort.Sort(sort.Reverse(sort.Float64Slice(sums)))
	return sums
}

// LyapunovSpectrumFlow estimates the full Lyapunov spectrum of the flow with
// field f by the QR method: an orthonormal frame is evolved with the
// variational equation and reorthonormalised every reortho steps, and the
// growth of the R diagonal is averaged. Exponents are returned in descending
// order, per unit time.
func LyapunovSpectrumFlow(f Field, x0 Vec, h float64, transient, n, reortho int) Vec {
	d := len(x0)
	if reortho < 1 {
		reortho = 1
	}
	x := x0.Clone()
	for i := 0; i < transient; i++ {
		x = StepRK4(f, x, h)
	}
	Q := Eye(d)
	sums := make(Vec, d)
	steps := 0
	for i := 0; i < n; i++ {
		// Evolve state and each basis column over one step.
		cols := make([]Vec, d)
		for k := 0; k < d; k++ {
			cols[k] = Q.Col(k)
		}
		nx := x
		for k := 0; k < d; k++ {
			xx, vv := stepFlowTangent(f, x, cols[k], h)
			cols[k] = vv
			nx = xx
		}
		x = nx
		// Assemble matrix with columns cols.
		Z := NewMat(d, d)
		for k := 0; k < d; k++ {
			for r := 0; r < d; r++ {
				Z[r][k] = cols[k][r]
			}
		}
		steps++
		if steps%reortho == 0 || i == n-1 {
			Qn, R := QR(Z)
			for k := 0; k < d; k++ {
				r := math.Abs(R[k][k])
				if r < 1e-300 {
					r = 1e-300
				}
				sums[k] += math.Log(r)
			}
			Q = Qn
		} else {
			Q = Z
		}
	}
	t := float64(n) * h
	for k := range sums {
		sums[k] /= t
	}
	sort.Sort(sort.Reverse(sort.Float64Slice(sums)))
	return sums
}

// SumExponents returns the sum of the Lyapunov exponents, which equals the
// average divergence of the flow (the phase-volume contraction rate).
func SumExponents(spectrum Vec) float64 {
	return spectrum.Sum()
}

// KaplanYorkeDimension returns the Kaplan-Yorke (Lyapunov) dimension implied by
// a Lyapunov spectrum. The exponents may be given in any order; they are
// sorted internally in descending order. The dimension is
//
//	D_KY = j + (sum_{i<=j} lambda_i) / |lambda_{j+1}|
//
// where j is the largest index for which the partial sum of exponents is still
// non-negative.
func KaplanYorkeDimension(spectrum Vec) float64 {
	if len(spectrum) == 0 {
		return 0
	}
	s := spectrum.Clone()
	sort.Sort(sort.Reverse(sort.Float64Slice(s)))
	var partial float64
	j := 0
	for j < len(s) {
		if partial+s[j] < 0 {
			break
		}
		partial += s[j]
		j++
	}
	if j == 0 {
		return 0
	}
	if j >= len(s) {
		return float64(len(s))
	}
	return float64(j) + partial/math.Abs(s[j])
}

// MetricEntropy returns the Kolmogorov-Sinai entropy estimated by Pesin's
// identity as the sum of the positive Lyapunov exponents.
func MetricEntropy(spectrum Vec) float64 {
	var s float64
	for _, l := range spectrum {
		if l > 0 {
			s += l
		}
	}
	return s
}
