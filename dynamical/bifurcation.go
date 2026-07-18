package dynamical

// BifurcationPoint holds the sampled long-term states of a one-dimensional map
// at a single parameter value. Param is the parameter; Values are the attractor
// samples collected after discarding transients. Plotting Values against Param
// for many parameters yields a bifurcation diagram.
type BifurcationPoint struct {
	Param  float64
	Values []float64
}

// Bifurcation samples the asymptotic behavior of a one-parameter family of
// one-dimensional maps. For each of steps parameter values evenly spaced from
// pmin to pmax (inclusive), it builds the map family(param), iterates from x0
// discarding transient points, then records the next samples iterates as the
// attractor samples for that parameter. It returns one [BifurcationPoint] per
// parameter value.
//
// steps must be at least 1; a single step samples pmin only.
func Bifurcation(family func(param float64) Map1D, pmin, pmax float64, steps int, x0 float64, transient, samples int) []BifurcationPoint {
	if steps < 1 {
		steps = 1
	}
	out := make([]BifurcationPoint, steps)
	for i := 0; i < steps; i++ {
		var p float64
		if steps == 1 {
			p = pmin
		} else {
			p = pmin + (pmax-pmin)*float64(i)/float64(steps-1)
		}
		f := family(p)
		x := x0
		for k := 0; k < transient; k++ {
			x = f(x)
		}
		vals := make([]float64, 0, samples)
		for k := 0; k < samples; k++ {
			x = f(x)
			vals = append(vals, x)
		}
		out[i] = BifurcationPoint{Param: p, Values: vals}
	}
	return out
}

// LogisticBifurcation is a convenience wrapper around [Bifurcation] for the
// logistic family x -> r*x*(1-x), sweeping r from rmin to rmax.
func LogisticBifurcation(rmin, rmax float64, steps int, x0 float64, transient, samples int) []BifurcationPoint {
	return Bifurcation(func(r float64) Map1D { return LogisticMap(r) }, rmin, rmax, steps, x0, transient, samples)
}
