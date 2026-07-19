package fem

// ElementStiffnessP1Interval returns the 2×2 P1 stiffness matrix for an
// interval element of length h, the integral of the product of basis-function
// derivatives.
func ElementStiffnessP1Interval(h float64) [2][2]float64 {
	c := 1 / h
	return [2][2]float64{
		{c, -c},
		{-c, c},
	}
}

// ElementMassP1Interval returns the 2×2 P1 mass matrix for an interval element
// of length h, the integral of products of basis functions.
func ElementMassP1Interval(h float64) [2][2]float64 {
	c := h / 6
	return [2][2]float64{
		{2 * c, c},
		{c, 2 * c},
	}
}

// ElementLoadP1Interval returns the 2-vector load contribution of the source f
// over the interval element [x0,x1], integrated with a Gauss rule of the given
// number of points (use at least 2 for smooth f).
func ElementLoadP1Interval(f func(float64) float64, x0, x1 float64, gauss int) [2]float64 {
	if gauss < 1 {
		gauss = 2
	}
	q := GaussLegendre1D(gauss)
	nodes, weights := q.MapToInterval(x0, x1)
	h := x1 - x0
	var out [2]float64
	for k := range nodes {
		xi := (nodes[k] - x0) / h
		n := ShapeP1Interval(xi)
		fv := f(nodes[k])
		out[0] += weights[k] * fv * n[0]
		out[1] += weights[k] * fv * n[1]
	}
	return out
}

// ElementStiffnessP2Interval returns the 3×3 P2 stiffness matrix for an
// interval element of length h in the node ordering (left, right, midpoint).
func ElementStiffnessP2Interval(h float64) [3][3]float64 {
	c := 1.0 / (3 * h)
	return [3][3]float64{
		{7 * c, c, -8 * c},
		{c, 7 * c, -8 * c},
		{-8 * c, -8 * c, 16 * c},
	}
}

// ElementMassP2Interval returns the 3×3 P2 mass matrix for an interval element
// of length h in the node ordering (left, right, midpoint).
func ElementMassP2Interval(h float64) [3][3]float64 {
	c := h / 30
	return [3][3]float64{
		{4 * c, -c, 2 * c},
		{-c, 4 * c, 2 * c},
		{2 * c, 2 * c, 16 * c},
	}
}

// ElementLoadP2Interval returns the 3-vector load contribution of f over the
// interval element [x0,x1] in the node ordering (left, right, midpoint).
func ElementLoadP2Interval(f func(float64) float64, x0, x1 float64, gauss int) [3]float64 {
	if gauss < 2 {
		gauss = 3
	}
	q := GaussLegendre1D(gauss)
	nodes, weights := q.MapToInterval(x0, x1)
	h := x1 - x0
	var out [3]float64
	for k := range nodes {
		xi := (nodes[k] - x0) / h
		n := ShapeP2Interval(xi)
		fv := f(nodes[k])
		for i := 0; i < 3; i++ {
			out[i] += weights[k] * fv * n[i]
		}
	}
	return out
}

// ElementStiffnessP1Triangle returns the 3×3 P1 stiffness matrix for the
// triangle with the given vertices.
func ElementStiffnessP1Triangle(v1, v2, v3 [2]float64) [3][3]float64 {
	grads, area := TriangleGradients(v1, v2, v3)
	a := area
	if a < 0 {
		a = -a
	}
	var k [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			k[i][j] = a * (grads[i][0]*grads[j][0] + grads[i][1]*grads[j][1])
		}
	}
	return k
}

// ElementMassP1Triangle returns the 3×3 P1 (consistent) mass matrix for the
// triangle with the given vertices.
func ElementMassP1Triangle(v1, v2, v3 [2]float64) [3][3]float64 {
	area := TriangleArea(v1, v2, v3)
	c := area / 12
	return [3][3]float64{
		{2 * c, c, c},
		{c, 2 * c, c},
		{c, c, 2 * c},
	}
}

// ElementLumpedMassP1Triangle returns the diagonal (row-sum lumped) P1 mass
// matrix for the triangle, each diagonal entry being area/3.
func ElementLumpedMassP1Triangle(v1, v2, v3 [2]float64) [3]float64 {
	c := TriangleArea(v1, v2, v3) / 3
	return [3]float64{c, c, c}
}

// ElementLoadP1Triangle returns the 3-vector load contribution of f over the
// triangle, integrated with a quadrature rule of the given degree.
func ElementLoadP1Triangle(f func(x, y float64) float64, v1, v2, v3 [2]float64, degree int) [3]float64 {
	rule := TriangleQuadrature(degree)
	area := TriangleArea(v1, v2, v3)
	var out [3]float64
	for q, bc := range rule.Bary {
		x := bc[0]*v1[0] + bc[1]*v2[0] + bc[2]*v3[0]
		y := bc[0]*v1[1] + bc[1]*v2[1] + bc[2]*v3[1]
		fv := f(x, y)
		w := rule.Weights[q] * area
		out[0] += w * fv * bc[0]
		out[1] += w * fv * bc[1]
		out[2] += w * fv * bc[2]
	}
	return out
}

// ElementStiffnessP2Triangle returns the 6×6 P2 stiffness matrix for the
// triangle with the given vertices, computed by quadrature.
func ElementStiffnessP2Triangle(v1, v2, v3 [2]float64) [6][6]float64 {
	rule := TriangleQuadrature(3)
	area := TriangleArea(v1, v2, v3)
	_, _, invT := TriangleJacobian(v1, v2, v3)
	var k [6][6]float64
	for q, bc := range rule.Bary {
		xi := bc[1]
		eta := bc[2]
		ref := GradShapeP2TriangleRef(xi, eta)
		var g [6][2]float64
		for i := 0; i < 6; i++ {
			g[i][0] = invT[0][0]*ref[i][0] + invT[0][1]*ref[i][1]
			g[i][1] = invT[1][0]*ref[i][0] + invT[1][1]*ref[i][1]
		}
		w := rule.Weights[q] * area
		for i := 0; i < 6; i++ {
			for j := 0; j < 6; j++ {
				k[i][j] += w * (g[i][0]*g[j][0] + g[i][1]*g[j][1])
			}
		}
	}
	return k
}

// ElementMassP2Triangle returns the 6×6 P2 mass matrix for the triangle with
// the given vertices, computed by quadrature.
func ElementMassP2Triangle(v1, v2, v3 [2]float64) [6][6]float64 {
	rule := TriangleQuadrature(5)
	area := TriangleArea(v1, v2, v3)
	var m [6][6]float64
	for q, bc := range rule.Bary {
		n := ShapeP2Triangle(bc)
		w := rule.Weights[q] * area
		for i := 0; i < 6; i++ {
			for j := 0; j < 6; j++ {
				m[i][j] += w * n[i] * n[j]
			}
		}
	}
	return m
}

// ElementLoadP2Triangle returns the 6-vector load contribution of f over the
// triangle, integrated with a degree-5 quadrature rule.
func ElementLoadP2Triangle(f func(x, y float64) float64, v1, v2, v3 [2]float64) [6]float64 {
	rule := TriangleQuadrature(5)
	area := TriangleArea(v1, v2, v3)
	var out [6]float64
	for q, bc := range rule.Bary {
		x := bc[0]*v1[0] + bc[1]*v2[0] + bc[2]*v3[0]
		y := bc[0]*v1[1] + bc[1]*v2[1] + bc[2]*v3[1]
		n := ShapeP2Triangle(bc)
		fv := f(x, y)
		w := rule.Weights[q] * area
		for i := 0; i < 6; i++ {
			out[i] += w * fv * n[i]
		}
	}
	return out
}
