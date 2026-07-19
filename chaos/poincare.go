package chaos

import "math"

// Section describes an oriented hyperplane {x : normal.(x-point)=0} used as a
// Poincare surface of section. Direction selects which crossings are kept:
// +1 for crossings with the field pointing along normal (increasing), -1 for
// the opposite orientation, and 0 for both.
type Section struct {
	// Normal is the (not necessarily unit) normal vector of the hyperplane.
	Normal Vec
	// Point is any point lying on the hyperplane.
	Point Vec
	// Direction selects the crossing orientation: +1, -1 or 0 (both).
	Direction int
}

// PlaneSection returns a Section for the coordinate hyperplane x[coord]=value,
// keeping crossings in the given direction.
func PlaneSection(dim, coord int, value float64, direction int) Section {
	n := make(Vec, dim)
	n[coord] = 1
	p := make(Vec, dim)
	p[coord] = value
	return Section{Normal: n, Point: p, Direction: direction}
}

// signedDistance returns normal.(x - point).
func (s Section) signedDistance(x Vec) float64 {
	return s.Normal.Dot(x.Sub(s.Point))
}

// PoincareSection integrates the flow with field f from x0 for n steps of size
// h and returns the sequence of intersection points with the section, found by
// detecting sign changes of the signed distance and linearly interpolating to
// the crossing. The transient steps are integrated first and discarded.
func PoincareSection(f Field, x0 Vec, h float64, transient, n int, s Section) ([]Vec, error) {
	x := x0.Clone()
	for i := 0; i < transient; i++ {
		x = StepRK4(f, x, h)
	}
	prev := x.Clone()
	dprev := s.signedDistance(prev)
	var crossings []Vec
	for i := 0; i < n; i++ {
		cur := StepRK4(f, prev, h)
		dcur := s.signedDistance(cur)
		if crossed(dprev, dcur, s.Direction) {
			// Linear interpolation for the crossing fraction.
			t := dprev / (dprev - dcur)
			pt := prev.AddScaled(t, cur.Sub(prev))
			crossings = append(crossings, pt)
		}
		prev, dprev = cur, dcur
	}
	if len(crossings) == 0 {
		return nil, ErrNoCrossing
	}
	return crossings, nil
}

// crossed reports whether the signed distance changed sign from d0 to d1 in a
// way consistent with the requested direction.
func crossed(d0, d1 float64, direction int) bool {
	if d0 == 0 {
		d0 = math.Copysign(1e-300, d1) // avoid double-count at exact zero
	}
	if (d0 < 0) == (d1 < 0) {
		return false
	}
	switch direction {
	case 1:
		return d0 < 0 && d1 >= 0
	case -1:
		return d0 > 0 && d1 <= 0
	default:
		return true
	}
}

// ReturnMap builds a first-return (Poincare) map from the crossings of a
// section: it returns the pairs (crossings[i], crossings[i+1]) as two parallel
// slices, projected onto coordinate coord. It is the discrete map whose
// fixed points are periodic orbits of the flow.
func ReturnMap(crossings []Vec, coord int) (in, out []float64) {
	if len(crossings) < 2 {
		return nil, nil
	}
	in = make([]float64, len(crossings)-1)
	out = make([]float64, len(crossings)-1)
	for i := 0; i+1 < len(crossings); i++ {
		in[i] = crossings[i][coord]
		out[i] = crossings[i+1][coord]
	}
	return in, out
}

// StroboscopicMap samples the trajectory of a periodically forced flow at
// times that are multiples of the forcing period, returning one point per
// period after the transient. It is the classical stroboscopic Poincare map
// used for forced oscillators; stepsPerPeriod steps of size h make up one
// period.
func StroboscopicMap(f Field, x0 Vec, h float64, stepsPerPeriod, transientPeriods, periods int) []Vec {
	x := x0.Clone()
	for i := 0; i < transientPeriods*stepsPerPeriod; i++ {
		x = StepRK4(f, x, h)
	}
	out := make([]Vec, 0, periods)
	for p := 0; p < periods; p++ {
		out = append(out, x.Clone())
		for i := 0; i < stepsPerPeriod; i++ {
			x = StepRK4(f, x, h)
		}
	}
	return out
}

// LocalMaxima returns the successive local maxima of the coordinate coord along
// the trajectory of f. These are the peaks used to build the classic Lorenz
// return map z_{n+1} vs z_n. The transient steps are discarded first.
func LocalMaxima(f Field, x0 Vec, h float64, transient, n, coord int) []float64 {
	x := x0.Clone()
	for i := 0; i < transient; i++ {
		x = StepRK4(f, x, h)
	}
	prev2 := x[coord]
	x = StepRK4(f, x, h)
	prev1 := x[coord]
	var maxima []float64
	for i := 0; i < n; i++ {
		x = StepRK4(f, x, h)
		cur := x[coord]
		if prev1 > prev2 && prev1 >= cur {
			maxima = append(maxima, prev1)
		}
		prev2, prev1 = prev1, cur
	}
	return maxima
}

// SuccessiveMaximaMap turns a sequence of successive maxima into the pairs
// (m[i], m[i+1]) of the return map, as parallel slices.
func SuccessiveMaximaMap(maxima []float64) (in, out []float64) {
	if len(maxima) < 2 {
		return nil, nil
	}
	in = make([]float64, len(maxima)-1)
	out = make([]float64, len(maxima)-1)
	for i := 0; i+1 < len(maxima); i++ {
		in[i] = maxima[i]
		out[i] = maxima[i+1]
	}
	return in, out
}
