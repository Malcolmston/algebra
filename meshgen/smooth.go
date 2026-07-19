package meshgen

// LaplacianSmooth returns a copy of the mesh after applying iterations passes of
// Laplacian smoothing. Each free vertex is moved a fraction relaxation toward
// the average of its edge-connected neighbours. When pinBoundary is true the
// boundary vertices are held fixed. The connectivity is preserved.
func LaplacianSmooth(m *Mesh, iterations int, relaxation float64, pinBoundary bool) *Mesh {
	out := m.Clone()
	if len(out.Vertices) == 0 || iterations <= 0 {
		return out
	}
	neigh := out.VertexNeighbors()
	fixed := boundaryMask(out, pinBoundary)
	for it := 0; it < iterations; it++ {
		next := make([]Vec2, len(out.Vertices))
		copy(next, out.Vertices)
		for v := range out.Vertices {
			if fixed[v] || len(neigh[v]) == 0 {
				continue
			}
			var sum Vec2
			for _, u := range neigh[v] {
				sum = sum.Add(out.Vertices[u])
			}
			avg := sum.Div(float64(len(neigh[v])))
			next[v] = out.Vertices[v].Lerp(avg, relaxation)
		}
		out.Vertices = next
	}
	return out
}

// LloydSmooth returns a copy of the mesh after applying iterations passes of
// centroidal (Lloyd) smoothing. Each free interior vertex is moved to the
// area-weighted centroid of the triangles incident to it, which drives the mesh
// toward a centroidal configuration. Boundary vertices are held fixed when
// pinBoundary is true.
func LloydSmooth(m *Mesh, iterations int, pinBoundary bool) *Mesh {
	out := m.Clone()
	if len(out.Vertices) == 0 || iterations <= 0 {
		return out
	}
	fixed := boundaryMask(out, pinBoundary)
	for it := 0; it < iterations; it++ {
		vtris := out.VertexTriangles()
		next := make([]Vec2, len(out.Vertices))
		copy(next, out.Vertices)
		for v := range out.Vertices {
			if fixed[v] || len(vtris[v]) == 0 {
				continue
			}
			var acc Vec2
			var wsum float64
			for _, ti := range vtris[v] {
				w := out.TriangleArea(ti)
				acc = acc.Add(out.TriangleCentroid(ti).Scale(w))
				wsum += w
			}
			if wsum > 0 {
				next[v] = acc.Div(wsum)
			}
		}
		out.Vertices = next
	}
	return out
}

// TaubinSmooth returns a copy of the mesh after applying iterations Taubin
// lambda/mu passes, a low-pass smoother that reduces the shrinkage of repeated
// Laplacian smoothing. Typical values are lambda = 0.5, mu = -0.53.
func TaubinSmooth(m *Mesh, iterations int, lambda, mu float64, pinBoundary bool) *Mesh {
	out := m.Clone()
	if len(out.Vertices) == 0 || iterations <= 0 {
		return out
	}
	neigh := out.VertexNeighbors()
	fixed := boundaryMask(out, pinBoundary)
	step := func(factor float64) {
		next := make([]Vec2, len(out.Vertices))
		copy(next, out.Vertices)
		for v := range out.Vertices {
			if fixed[v] || len(neigh[v]) == 0 {
				continue
			}
			var sum Vec2
			for _, u := range neigh[v] {
				sum = sum.Add(out.Vertices[u])
			}
			avg := sum.Div(float64(len(neigh[v])))
			delta := avg.Sub(out.Vertices[v])
			next[v] = out.Vertices[v].Add(delta.Scale(factor))
		}
		out.Vertices = next
	}
	for it := 0; it < iterations; it++ {
		step(lambda)
		step(mu)
	}
	return out
}

func boundaryMask(m *Mesh, pinBoundary bool) []bool {
	fixed := make([]bool, len(m.Vertices))
	if pinBoundary {
		for _, v := range m.BoundaryVertices() {
			fixed[v] = true
		}
	}
	return fixed
}
