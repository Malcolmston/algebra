package fem

import "math"

// EvalP1_1D evaluates the P1 finite element function with the given nodal
// values at coordinate x. Points outside the mesh are clamped to the nearest
// element.
func EvalP1_1D(m *Mesh1D, u Vector, x float64) float64 {
	e := locateElement1D(m, x)
	i, j := m.ElementNodes(e)
	h := m.ElementLength(e)
	xi := (x - m.Nodes[i]) / h
	n := ShapeP1Interval(xi)
	return u[i]*n[0] + u[j]*n[1]
}

// EvalP2_1D evaluates the P2 finite element function with the given dof values
// (as returned by SolvePoisson1DP2) at coordinate x.
func EvalP2_1D(m *Mesh1D, u Vector, x float64) float64 {
	e := locateElement1D(m, x)
	conn := m.P2Connectivity()[e]
	i, _ := m.ElementNodes(e)
	h := m.ElementLength(e)
	xi := (x - m.Nodes[i]) / h
	n := ShapeP2Interval(xi)
	return u[conn[0]]*n[0] + u[conn[1]]*n[1] + u[conn[2]]*n[2]
}

func locateElement1D(m *Mesh1D, x float64) int {
	if x <= m.Nodes[0] {
		return 0
	}
	last := m.NumElements() - 1
	if x >= m.Nodes[m.NumNodes()-1] {
		return last
	}
	// binary search for the containing element
	lo, hi := 0, m.NumNodes()-1
	for hi-lo > 1 {
		mid := (lo + hi) / 2
		if m.Nodes[mid] <= x {
			lo = mid
		} else {
			hi = mid
		}
	}
	return lo
}

// L2Error1D returns the L2 norm of the difference between the P1 solution u and
// the exact function exact, integrated element by element.
func L2Error1D(m *Mesh1D, u Vector, exact func(float64) float64, gauss int) float64 {
	if gauss < 2 {
		gauss = 4
	}
	q := GaussLegendre1D(gauss)
	var sum float64
	for e := 0; e < m.NumElements(); e++ {
		i, j := m.ElementNodes(e)
		h := m.ElementLength(e)
		nodes, weights := q.MapToInterval(m.Nodes[i], m.Nodes[j])
		for k := range nodes {
			xi := (nodes[k] - m.Nodes[i]) / h
			sh := ShapeP1Interval(xi)
			uh := u[i]*sh[0] + u[j]*sh[1]
			d := uh - exact(nodes[k])
			sum += weights[k] * d * d
		}
	}
	return math.Sqrt(sum)
}

// H1SeminormError1D returns the H1 seminorm of the error (the L2 norm of the
// derivative error) for the P1 solution u.
func H1SeminormError1D(m *Mesh1D, u Vector, exactDeriv func(float64) float64, gauss int) float64 {
	if gauss < 2 {
		gauss = 4
	}
	q := GaussLegendre1D(gauss)
	var sum float64
	for e := 0; e < m.NumElements(); e++ {
		i, j := m.ElementNodes(e)
		h := m.ElementLength(e)
		duh := (u[j] - u[i]) / h
		nodes, weights := q.MapToInterval(m.Nodes[i], m.Nodes[j])
		for k := range nodes {
			d := duh - exactDeriv(nodes[k])
			sum += weights[k] * d * d
		}
	}
	return math.Sqrt(sum)
}

// H1Error1D returns the full H1 norm of the error, sqrt(L2^2 + seminorm^2).
func H1Error1D(m *Mesh1D, u Vector, exact, exactDeriv func(float64) float64, gauss int) float64 {
	l2 := L2Error1D(m, u, exact, gauss)
	semi := H1SeminormError1D(m, u, exactDeriv, gauss)
	return math.Sqrt(l2*l2 + semi*semi)
}

// L2Norm1D returns the L2 norm of a function over the meshed interval.
func L2Norm1D(m *Mesh1D, f func(float64) float64, gauss int) float64 {
	if gauss < 2 {
		gauss = 4
	}
	q := GaussLegendre1D(gauss)
	// composite integration element by element for accuracy
	var sum float64
	for e := 0; e < m.NumElements(); e++ {
		i, j := m.ElementNodes(e)
		en, ew := q.MapToInterval(m.Nodes[i], m.Nodes[j])
		for k := range en {
			v := f(en[k])
			sum += ew[k] * v * v
		}
	}
	return math.Sqrt(sum)
}

// L2Error2D returns the L2 norm of the difference between the P1 solution u and
// the exact function exact on the triangular mesh.
func L2Error2D(m *Mesh2D, u Vector, exact func(x, y float64) float64, degree int) float64 {
	if degree < 2 {
		degree = 4
	}
	rule := TriangleQuadrature(degree)
	var sum float64
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		area := TriangleArea(v1, v2, v3)
		tri := m.Triangles[t]
		for q, bc := range rule.Bary {
			x := bc[0]*v1[0] + bc[1]*v2[0] + bc[2]*v3[0]
			y := bc[0]*v1[1] + bc[1]*v2[1] + bc[2]*v3[1]
			uh := bc[0]*u[tri[0]] + bc[1]*u[tri[1]] + bc[2]*u[tri[2]]
			d := uh - exact(x, y)
			sum += rule.Weights[q] * area * d * d
		}
	}
	return math.Sqrt(sum)
}

// H1SeminormError2D returns the H1 seminorm of the error (the L2 norm of the
// gradient error) for the P1 solution u, given the exact gradient.
func H1SeminormError2D(m *Mesh2D, u Vector, exactGrad func(x, y float64) (gx, gy float64), degree int) float64 {
	if degree < 2 {
		degree = 4
	}
	rule := TriangleQuadrature(degree)
	var sum float64
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		grads, area := TriangleGradients(v1, v2, v3)
		a := area
		if a < 0 {
			a = -a
		}
		tri := m.Triangles[t]
		var ghx, ghy float64
		for i := 0; i < 3; i++ {
			ghx += u[tri[i]] * grads[i][0]
			ghy += u[tri[i]] * grads[i][1]
		}
		for q, bc := range rule.Bary {
			x := bc[0]*v1[0] + bc[1]*v2[0] + bc[2]*v3[0]
			y := bc[0]*v1[1] + bc[1]*v2[1] + bc[2]*v3[1]
			ex, ey := exactGrad(x, y)
			dx := ghx - ex
			dy := ghy - ey
			sum += rule.Weights[q] * a * (dx*dx + dy*dy)
		}
	}
	return math.Sqrt(sum)
}

// H1Error2D returns the full H1 norm of the error, sqrt(L2^2 + seminorm^2).
func H1Error2D(m *Mesh2D, u Vector, exact func(x, y float64) float64, exactGrad func(x, y float64) (gx, gy float64), degree int) float64 {
	l2 := L2Error2D(m, u, exact, degree)
	semi := H1SeminormError2D(m, u, exactGrad, degree)
	return math.Sqrt(l2*l2 + semi*semi)
}

// L2Norm2D returns the L2 norm of a function over the triangular mesh.
func L2Norm2D(m *Mesh2D, f func(x, y float64) float64, degree int) float64 {
	if degree < 2 {
		degree = 4
	}
	rule := TriangleQuadrature(degree)
	var sum float64
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		area := TriangleArea(v1, v2, v3)
		for q, bc := range rule.Bary {
			x := bc[0]*v1[0] + bc[1]*v2[0] + bc[2]*v3[0]
			y := bc[0]*v1[1] + bc[1]*v2[1] + bc[2]*v3[1]
			v := f(x, y)
			sum += rule.Weights[q] * area * v * v
		}
	}
	return math.Sqrt(sum)
}

// InterpolateNodal1D returns the nodal interpolant of f on the mesh (its value
// at each node).
func InterpolateNodal1D(m *Mesh1D, f func(float64) float64) Vector {
	u := make(Vector, m.NumNodes())
	for i, x := range m.Nodes {
		u[i] = f(x)
	}
	return u
}

// InterpolateNodal2D returns the nodal interpolant of f on the mesh.
func InterpolateNodal2D(m *Mesh2D, f func(x, y float64) float64) Vector {
	u := make(Vector, m.NumNodes())
	for i, p := range m.Nodes {
		u[i] = f(p[0], p[1])
	}
	return u
}
