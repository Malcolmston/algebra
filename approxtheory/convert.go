package approxtheory

import "math"

// ChebyshevToMonomial converts Chebyshev coefficients on [-1, 1] into the
// coefficients of the equivalent polynomial in the monomial basis (ascending
// order). It uses the recurrence T_{k+1} = 2x T_k - T_{k-1} on coefficient
// vectors.
func ChebyshevToMonomial(cheb []float64) []float64 {
	n := len(cheb)
	if n == 0 {
		return nil
	}
	// Tprev = T_0 = [1]; Tcur = T_1 = [0,1].
	Tprev := []float64{1}
	out := PolyScale(Tprev, cheb[0])
	if n == 1 {
		return out
	}
	Tcur := []float64{0, 1}
	out = PolyAdd(out, PolyScale(Tcur, cheb[1]))
	for k := 2; k < n; k++ {
		// T_k = 2x*T_{k-1} - T_{k-2}
		shifted := append([]float64{0}, Tcur...) // multiply by x
		Tnext := PolySub(PolyScale(shifted, 2), Tprev)
		out = PolyAdd(out, PolyScale(Tnext, cheb[k]))
		Tprev, Tcur = Tcur, Tnext
	}
	return PolyNormalize(out)
}

// MonomialToChebyshev converts monomial coefficients of a polynomial on
// [-1, 1] into Chebyshev coefficients. It expresses each power x^k in the
// Chebyshev basis and accumulates, using x*T_j = (T_{|j+1|}+T_{|j-1|})/2 with
// the convention T_{-1} = T_1.
func MonomialToChebyshev(mono []float64) []float64 {
	n := len(mono)
	if n == 0 {
		return nil
	}
	out := make([]float64, n)
	// Represent x^k in Chebyshev basis by repeated multiplication using the
	// correct identity x*T_j = (T_{|j+1|}+T_{|j-1|})/2 with T_{-1}=T_1.
	cur := make([]float64, n)
	cur[0] = 1 // x^0
	for k := 0; k < n; k++ {
		for j := 0; j < n; j++ {
			out[j] += mono[k] * cur[j]
		}
		if k+1 < n {
			next := make([]float64, n)
			for j := 0; j < n; j++ {
				c := cur[j]
				if c == 0 {
					continue
				}
				up := j + 1
				if up < n {
					next[up] += 0.5 * c
				}
				down := j - 1
				if down < 0 {
					down = 1 // T_{-1} = T_1
				}
				if down < n {
					next[down] += 0.5 * c
				}
			}
			cur = next
		}
	}
	return out
}

// ChebProduct returns the Chebyshev series of the product of two series sharing
// the same domain. The product is computed exactly in the Chebyshev basis via
// 2 T_i T_j = T_{i+j} + T_{|i-j|}.
func ChebProduct(a, b *ChebSeries) *ChebSeries {
	na, nb := len(a.Coeffs), len(b.Coeffs)
	if na == 0 || nb == 0 {
		return &ChebSeries{Coeffs: []float64{0}, A: a.A, B: a.B}
	}
	out := make([]float64, na+nb-1)
	for i := 0; i < na; i++ {
		ai := a.Coeffs[i]
		if ai == 0 {
			continue
		}
		for j := 0; j < nb; j++ {
			bj := b.Coeffs[j]
			if bj == 0 {
				continue
			}
			p := 0.5 * ai * bj
			out[i+j] += p
			d := i - j
			if d < 0 {
				d = -d
			}
			out[d] += p
		}
	}
	return &ChebSeries{Coeffs: out, A: a.A, B: a.B}
}

// ChebFromMonomial builds a ChebSeries on [a, b] from a monomial polynomial
// defined on that same interval (the polynomial is first rescaled to [-1, 1]).
func ChebFromMonomial(mono []float64, a, b float64) *ChebSeries {
	// Substitute x = (A+B)/2 + (B-A)/2 t into the monomial polynomial to get a
	// polynomial in t, then convert to Chebyshev.
	shift := 0.5 * (a + b)
	half := 0.5 * (b - a)
	// build poly in t: sum mono[k] * (shift + half*t)^k
	polyT := []float64{}
	for k := len(mono) - 1; k >= 0; k-- {
		// Horner in the linear polynomial (shift + half*t)
		polyT = PolyMul(polyT, []float64{shift, half})
		polyT = PolyAdd(polyT, []float64{mono[k]})
	}
	cheb := MonomialToChebyshev(polyT)
	if math.IsNaN(cheb[0]) {
		cheb[0] = 0
	}
	return &ChebSeries{Coeffs: cheb, A: a, B: b}
}
