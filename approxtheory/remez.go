package approxtheory

import "math"

// RemezResult holds the outcome of a Remez minimax approximation: the monomial
// coefficients of the best polynomial (ascending order), the equioscillating
// error magnitude, the number of exchange iterations performed and whether the
// iteration converged to the requested tolerance.
type RemezResult struct {
	Coeffs     []float64 // monomial coefficients, ascending
	Error      float64   // minimax (leveled) error estimate
	Iterations int
	Converged  bool
}

// Eval evaluates the minimax polynomial at x.
func (r *RemezResult) Eval(x float64) float64 { return Polyval(r.Coeffs, x) }

// RemezPoly computes the degree-n minimax (best uniform) polynomial
// approximation of f on [a, b] using the Remez exchange algorithm. maxIter
// bounds the number of exchanges and tol is the relative leveling tolerance on
// the extrema. A degree-n polynomial has n+1 coefficients and the reference
// carries n+2 points.
func RemezPoly(f func(float64) float64, n int, a, b float64, maxIter int, tol float64) (*RemezResult, error) {
	if a >= b {
		return nil, ErrDimensionMismatch
	}
	if n < 0 {
		n = 0
	}
	if maxIter <= 0 {
		maxIter = 100
	}
	if tol <= 0 {
		tol = 1e-12
	}
	m := n + 2 // number of reference points

	// Initial reference: Chebyshev extrema on [a,b].
	ref := ChebPoints(n+1, a, b)
	if len(ref) != m {
		// Fall back to equispaced if the count is off.
		ref = Linspace(a, b, m)
	}

	coeffs := make([]float64, n+1)
	var lastErr float64
	converged := false
	iter := 0
	for ; iter < maxIter; iter++ {
		c, e, err := remezSolve(f, ref, n)
		if err != nil {
			return nil, err
		}
		coeffs = c
		lastErr = math.Abs(e)

		// Locate the extrema of the error on a fine grid.
		newRef, maxErr, minErr := remezExtrema(f, coeffs, a, b, m)
		if len(newRef) == m {
			ref = newRef
		}
		if maxErr > 0 && (maxErr-minErr) <= tol*maxErr {
			lastErr = maxErr
			converged = true
			iter++
			break
		}
		lastErr = maxErr
	}
	return &RemezResult{Coeffs: coeffs, Error: lastErr, Iterations: iter, Converged: converged}, nil
}

// remezSolve solves the linear system that levels the error with alternating
// signs at the reference points, returning the polynomial coefficients and the
// signed leveled error E.
func remezSolve(f func(float64) float64, ref []float64, n int) ([]float64, float64, error) {
	m := len(ref) // n+2
	A := make([][]float64, m)
	rhs := make([]float64, m)
	for i := 0; i < m; i++ {
		row := make([]float64, m)
		xp := 1.0
		for j := 0; j <= n; j++ {
			row[j] = xp
			xp *= ref[i]
		}
		sign := 1.0
		if i%2 == 1 {
			sign = -1.0
		}
		row[n+1] = sign
		A[i] = row
		rhs[i] = f(ref[i])
	}
	sol, err := solveLinear(A, rhs)
	if err != nil {
		return nil, 0, err
	}
	coeffs := make([]float64, n+1)
	copy(coeffs, sol[:n+1])
	return coeffs, sol[n+1], nil
}

// remezExtrema finds a fresh alternating reference of m = n+2 points from the
// local extrema of the error e(x) = f(x) - p(x) on a dense grid, and reports
// the largest and smallest extreme magnitudes.
func remezExtrema(f func(float64) float64, coeffs []float64, a, b float64, m int) ([]float64, float64, float64) {
	M := 20 * m
	if M < 400 {
		M = 400
	}
	grid := Linspace(a, b, M+1)
	e := make([]float64, len(grid))
	for i, x := range grid {
		e[i] = f(x) - Polyval(coeffs, x)
	}
	// Candidate extrema: endpoints plus interior local extrema of |e|.
	type cand struct {
		x    float64
		val  float64 // signed error
		mag  float64
		sign int
	}
	var cs []cand
	add := func(x, v float64) {
		s := 0
		if v > 0 {
			s = 1
		} else if v < 0 {
			s = -1
		}
		cs = append(cs, cand{x: x, val: v, mag: math.Abs(v), sign: s})
	}
	h := grid[1] - grid[0]
	add(grid[0], e[0])
	for i := 1; i < len(grid)-1; i++ {
		if (e[i]-e[i-1])*(e[i+1]-e[i]) <= 0 { // local extremum of e
			// Parabolic refinement of the extremum location and value.
			denom := e[i-1] - 2*e[i] + e[i+1]
			xr, vr := grid[i], e[i]
			if denom != 0 {
				delta := 0.5 * (e[i-1] - e[i+1]) / denom
				if delta > -1 && delta < 1 {
					xr = grid[i] + delta*h
					vr = e[i] - 0.25*(e[i-1]-e[i+1])*delta
				}
			}
			add(xr, vr)
		}
	}
	add(grid[len(grid)-1], e[len(grid)-1])

	// Collapse runs of equal sign, keeping the largest magnitude in each run.
	var alt []cand
	for _, c := range cs {
		if c.sign == 0 {
			continue
		}
		if len(alt) == 0 || alt[len(alt)-1].sign != c.sign {
			alt = append(alt, c)
		} else if c.mag > alt[len(alt)-1].mag {
			alt[len(alt)-1] = c
		}
	}
	// Trim to m points, dropping the smaller endpoint each time.
	for len(alt) > m {
		if alt[0].mag <= alt[len(alt)-1].mag {
			alt = alt[1:]
		} else {
			alt = alt[:len(alt)-1]
		}
	}
	out := make([]float64, 0, len(alt))
	maxMag, minMag := 0.0, math.Inf(1)
	for _, c := range alt {
		out = append(out, c.x)
		if c.mag > maxMag {
			maxMag = c.mag
		}
		if c.mag < minMag {
			minMag = c.mag
		}
	}
	if len(alt) == 0 {
		minMag = 0
	}
	return out, maxMag, minMag
}

// MinimaxError returns the achieved uniform error of a monomial polynomial p
// against f on [a, b], sampled on a fine grid.
func MinimaxError(f func(float64) float64, coeffs []float64, a, b float64) float64 {
	return MaxError(f, func(x float64) float64 { return Polyval(coeffs, x) }, a, b)
}
