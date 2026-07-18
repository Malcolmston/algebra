package diffgeo

// GeodesicState is a point on a curve in a surface's parameter domain together
// with its parametric velocity: the parameters (U, V) and their rates of change
// (DU, DV) with respect to the curve parameter.
type GeodesicState struct {
	U, V   float64
	DU, DV float64
}

// GeodesicAcceleration returns the parameter-space acceleration (d²u, d²v)
// dictated by the geodesic equations at state st on surface s:
//
//	d²u = −(Γ⁰₀₀ u'² + 2Γ⁰₀₁ u'v' + Γ⁰₁₁ v'²)
//	d²v = −(Γ¹₀₀ u'² + 2Γ¹₀₁ u'v' + Γ¹₁₁ v'²)
//
// where the Γ are the [Christoffel] symbols at (U, V). A curve satisfying these
// equations is a geodesic: it has zero tangential acceleration within the
// surface.
func GeodesicAcceleration(s Surface, st GeodesicState) (accU, accV float64) {
	c := ChristoffelSymbols(s, st.U, st.V)
	du, dv := st.DU, st.DV
	accU = -(c.Gamma[0][0][0]*du*du + 2*c.Gamma[0][0][1]*du*dv + c.Gamma[0][1][1]*dv*dv)
	accV = -(c.Gamma[1][0][0]*du*du + 2*c.Gamma[1][0][1]*du*dv + c.Gamma[1][1][1]*dv*dv)
	return accU, accV
}

// diffgeoGeoDeriv returns the time-derivative of a geodesic state: position
// rates equal the stored velocities and velocity rates come from
// [GeodesicAcceleration].
func diffgeoGeoDeriv(s Surface, st GeodesicState) GeodesicState {
	accU, accV := GeodesicAcceleration(s, st)
	return GeodesicState{U: st.DU, V: st.DV, DU: accU, DV: accV}
}

// diffgeoGeoAdd forms st + h·d component-wise, used inside the RK4 stages.
func diffgeoGeoAdd(st, d GeodesicState, h float64) GeodesicState {
	return GeodesicState{
		U:  st.U + h*d.U,
		V:  st.V + h*d.V,
		DU: st.DU + h*d.DU,
		DV: st.DV + h*d.DV,
	}
}

// GeodesicStep advances one geodesic [GeodesicState] by a step h of the curve
// parameter using the classical fourth-order Runge-Kutta method.
func GeodesicStep(s Surface, st GeodesicState, h float64) GeodesicState {
	k1 := diffgeoGeoDeriv(s, st)
	k2 := diffgeoGeoDeriv(s, diffgeoGeoAdd(st, k1, h/2))
	k3 := diffgeoGeoDeriv(s, diffgeoGeoAdd(st, k2, h/2))
	k4 := diffgeoGeoDeriv(s, diffgeoGeoAdd(st, k3, h))
	return GeodesicState{
		U:  st.U + h/6*(k1.U+2*k2.U+2*k3.U+k4.U),
		V:  st.V + h/6*(k1.V+2*k2.V+2*k3.V+k4.V),
		DU: st.DU + h/6*(k1.DU+2*k2.DU+2*k3.DU+k4.DU),
		DV: st.DV + h/6*(k1.DV+2*k2.DV+2*k3.DV+k4.DV),
	}
}

// GeodesicPath integrates a geodesic on surface s from the initial
// [GeodesicState] start over a curve-parameter length of tEnd, using steps
// equal-sized RK4 steps. It returns the sequence of states, of length steps+1,
// beginning with start. A geodesic is the surface analogue of a straight line:
// on a sphere the paths are great circles, on a plane they are straight lines.
func GeodesicPath(s Surface, start GeodesicState, tEnd float64, steps int) []GeodesicState {
	if steps < 1 {
		steps = 1
	}
	h := tEnd / float64(steps)
	path := make([]GeodesicState, 0, steps+1)
	path = append(path, start)
	st := start
	for i := 0; i < steps; i++ {
		st = GeodesicStep(s, st, h)
		path = append(path, st)
	}
	return path
}

// GeodesicLength returns the arc length in space of the geodesic sampled by
// [GeodesicPath], measured as the total length of the polyline through the
// image points s(U, V) of the path states. With enough steps it approaches the
// true intrinsic length between the endpoints.
func GeodesicLength(s Surface, path []GeodesicState) float64 {
	var total float64
	for i := 1; i < len(path); i++ {
		p0 := s(path[i-1].U, path[i-1].V)
		p1 := s(path[i].U, path[i].V)
		total += p0.Distance(p1)
	}
	return total
}

// diffgeoTransportState bundles a position/velocity geodesic state with the two
// contravariant components (W0, W1) of the vector being transported.
type diffgeoTransportState struct {
	g      GeodesicState
	W0, W1 float64
}

// diffgeoTransportDeriv returns the time-derivative of a transport state. The
// position and velocity evolve by the geodesic equations while the vector obeys
// the parallel-transport law dWᵏ/dt = −Γᵏᵢⱼ (dxⁱ/dt) Wʲ.
func diffgeoTransportDeriv(s Surface, ts diffgeoTransportState) diffgeoTransportState {
	c := ChristoffelSymbols(s, ts.g.U, ts.g.V)
	gd := diffgeoGeoDeriv(s, ts.g)
	du, dv := ts.g.DU, ts.g.DV
	// dWᵏ = −(Γᵏ₀ⱼ du + Γᵏ₁ⱼ dv) Wʲ
	dW0 := -((c.Gamma[0][0][0]*du+c.Gamma[0][1][0]*dv)*ts.W0 +
		(c.Gamma[0][0][1]*du+c.Gamma[0][1][1]*dv)*ts.W1)
	dW1 := -((c.Gamma[1][0][0]*du+c.Gamma[1][1][0]*dv)*ts.W0 +
		(c.Gamma[1][0][1]*du+c.Gamma[1][1][1]*dv)*ts.W1)
	return diffgeoTransportState{g: gd, W0: dW0, W1: dW1}
}

func diffgeoTransportAdd(ts, d diffgeoTransportState, h float64) diffgeoTransportState {
	return diffgeoTransportState{
		g:  diffgeoGeoAdd(ts.g, d.g, h),
		W0: ts.W0 + h*d.W0,
		W1: ts.W1 + h*d.W1,
	}
}

// ParallelTransport transports the tangent vector with parameter-space
// components (w0, w1) along the geodesic starting at start on surface s, over a
// curve-parameter length tEnd using steps RK4 steps. It returns the geodesic
// [GeodesicPath] and, aligned index-for-index with it, the transported vector
// components at each state.
//
// Parallel transport preserves the metric inner product; in particular the
// length √(E w0² + 2F w0 w1 + G w1²) of the transported vector is invariant
// along the path, a useful correctness check.
func ParallelTransport(s Surface, start GeodesicState, w0, w1, tEnd float64, steps int) ([]GeodesicState, [][2]float64) {
	if steps < 1 {
		steps = 1
	}
	h := tEnd / float64(steps)
	path := make([]GeodesicState, 0, steps+1)
	vecs := make([][2]float64, 0, steps+1)
	ts := diffgeoTransportState{g: start, W0: w0, W1: w1}
	path = append(path, ts.g)
	vecs = append(vecs, [2]float64{ts.W0, ts.W1})
	for i := 0; i < steps; i++ {
		k1 := diffgeoTransportDeriv(s, ts)
		k2 := diffgeoTransportDeriv(s, diffgeoTransportAdd(ts, k1, h/2))
		k3 := diffgeoTransportDeriv(s, diffgeoTransportAdd(ts, k2, h/2))
		k4 := diffgeoTransportDeriv(s, diffgeoTransportAdd(ts, k3, h))
		ts = diffgeoTransportState{
			g: GeodesicState{
				U:  ts.g.U + h/6*(k1.g.U+2*k2.g.U+2*k3.g.U+k4.g.U),
				V:  ts.g.V + h/6*(k1.g.V+2*k2.g.V+2*k3.g.V+k4.g.V),
				DU: ts.g.DU + h/6*(k1.g.DU+2*k2.g.DU+2*k3.g.DU+k4.g.DU),
				DV: ts.g.DV + h/6*(k1.g.DV+2*k2.g.DV+2*k3.g.DV+k4.g.DV),
			},
			W0: ts.W0 + h/6*(k1.W0+2*k2.W0+2*k3.W0+k4.W0),
			W1: ts.W1 + h/6*(k1.W1+2*k2.W1+2*k3.W1+k4.W1),
		}
		path = append(path, ts.g)
		vecs = append(vecs, [2]float64{ts.W0, ts.W1})
	}
	return path, vecs
}
