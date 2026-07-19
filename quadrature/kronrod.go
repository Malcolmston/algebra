package quadrature

import "math"

// Result bundles a computed integral value with an estimate of its absolute
// error and the number of function evaluations spent producing it.
type Result struct {
	Value   float64 // the estimated integral
	AbsErr  float64 // an estimate of the absolute error
	Evals   int     // number of function evaluations
	Success bool    // whether the requested tolerance was met
}

// Gauss-Kronrod 7-15 abscissae and weights (QUADPACK constants). xgk holds the
// abscissae of the 15-point Kronrod rule on the positive half of [-1, 1],
// ordered from the outside inward with the centre last. Even indices are the
// points added to the embedded 7-point Gauss rule; odd indices are the Gauss
// nodes.
var (
	gk15xgk = [8]float64{
		0.991455371120813,
		0.949107912342759,
		0.864864423359769,
		0.741531185599394,
		0.586087235467691,
		0.405845151377397,
		0.207784955007898,
		0.000000000000000,
	}
	gk15wgk = [8]float64{
		0.022935322010529,
		0.063092092629979,
		0.104790010322250,
		0.140653259715525,
		0.169004726639267,
		0.190350578064785,
		0.204432940075298,
		0.209482141084728,
	}
	gk15wg = [4]float64{
		0.129484966168870,
		0.279705391489277,
		0.381830050505119,
		0.417959183673469,
	}
)

// GaussKronrod15 returns the abscissae and weights of the 15-point
// Gauss-Kronrod rule on [-1, 1] together with the weights of the embedded
// 7-point Gauss-Legendre rule (the gauss slice has length 7 and its weights
// correspond to the seven odd-indexed Kronrod abscissae). The abscissae are
// sorted ascending.
func GaussKronrod15() (nodes, kronrod, gauss []float64) {
	nodes = make([]float64, 15)
	kronrod = make([]float64, 15)
	gauss = make([]float64, 7)
	// Positions 0..14, centre at 7. gk15xgk[0] is the outermost abscissa.
	for i := 0; i < 7; i++ {
		nodes[i] = -gk15xgk[i]
		kronrod[i] = gk15wgk[i]
		nodes[14-i] = gk15xgk[i]
		kronrod[14-i] = gk15wgk[i]
	}
	nodes[7] = 0
	kronrod[7] = gk15wgk[7]
	// Gauss nodes are the odd-indexed Kronrod abscissae (1,3,5) plus centre.
	gi := 0
	for j := 0; j < 3; j++ {
		gauss[gi] = gk15wg[j]
		gi++
	}
	gauss[gi] = gk15wg[3]
	gi++
	for j := 2; j >= 0; j-- {
		gauss[gi] = gk15wg[j]
		gi++
	}
	return
}

// GaussKronrod15Eval applies the 15-point Gauss-Kronrod rule to f over [a, b],
// returning the Kronrod estimate together with a QUADPACK-style estimate of
// the absolute error obtained from the difference between the Kronrod and the
// embedded Gauss estimates.
func GaussKronrod15Eval(f Func, a, b float64) (value, absErr float64) {
	centr := 0.5 * (a + b)
	hlgth := 0.5 * (b - a)
	dhlgth := math.Abs(hlgth)

	var fv1, fv2 [7]float64
	fc := f(centr)
	resg := fc * gk15wg[3]
	resk := fc * gk15wgk[7]
	resabs := math.Abs(resk)

	for j := 0; j < 3; j++ {
		jtw := 2*j + 1
		absc := hlgth * gk15xgk[jtw]
		f1 := f(centr - absc)
		f2 := f(centr + absc)
		fv1[jtw] = f1
		fv2[jtw] = f2
		fsum := f1 + f2
		resg += gk15wg[j] * fsum
		resk += gk15wgk[jtw] * fsum
		resabs += gk15wgk[jtw] * (math.Abs(f1) + math.Abs(f2))
	}
	for j := 0; j < 4; j++ {
		jtwm1 := 2 * j
		absc := hlgth * gk15xgk[jtwm1]
		f1 := f(centr - absc)
		f2 := f(centr + absc)
		fv1[jtwm1] = f1
		fv2[jtwm1] = f2
		fsum := f1 + f2
		resk += gk15wgk[jtwm1] * fsum
		resabs += gk15wgk[jtwm1] * (math.Abs(f1) + math.Abs(f2))
	}

	reskh := resk * 0.5
	resasc := gk15wgk[7] * math.Abs(fc-reskh)
	for j := 0; j < 7; j++ {
		resasc += gk15wgk[j] * (math.Abs(fv1[j]-reskh) + math.Abs(fv2[j]-reskh))
	}

	value = resk * hlgth
	resabs *= dhlgth
	resasc *= dhlgth
	absErr = math.Abs((resk - resg) * hlgth)
	if resasc != 0 && absErr != 0 {
		absErr = resasc * math.Min(1, math.Pow(200*absErr/resasc, 1.5))
	}
	const uflow = 2.2250738585072014e-308
	const epmach = 2.220446049250313e-16
	if resabs > uflow/(50*epmach) {
		absErr = math.Max(epmach*50*resabs, absErr)
	}
	return value, absErr
}

// IntegrateGaussKronrod is an alias for GaussKronrod15Eval returning only the
// integral estimate over [a, b].
func IntegrateGaussKronrod(f Func, a, b float64) float64 {
	v, _ := GaussKronrod15Eval(f, a, b)
	return v
}

// AdaptiveGaussKronrod approximates the integral of f over [a, b] to the
// requested absolute tolerance by recursively subdividing the interval,
// applying the 15-point Gauss-Kronrod rule on each panel and bisecting panels
// whose local error estimate is too large. It returns the estimate, an error
// estimate and the number of subdivisions.
func AdaptiveGaussKronrod(f Func, a, b, tol float64) Result {
	if a == b {
		return Result{Value: 0, AbsErr: 0, Evals: 0, Success: true}
	}
	if tol <= 0 {
		tol = 1e-10
	}
	const maxDepth = 50
	var evals int
	var totalErr float64
	var rec func(a, b, tol float64, depth int) float64
	rec = func(a, b, tol float64, depth int) float64 {
		v, e := GaussKronrod15Eval(f, a, b)
		evals += 15
		if e <= tol || depth >= maxDepth {
			totalErr += e
			return v
		}
		m := 0.5 * (a + b)
		return rec(a, m, tol/2, depth+1) + rec(m, b, tol/2, depth+1)
	}
	v := rec(a, b, tol, 0)
	return Result{Value: v, AbsErr: totalErr, Evals: evals, Success: totalErr <= tol}
}
