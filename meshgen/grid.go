package meshgen

// GridMesh returns a triangulated regular grid of nx by ny nodes (so nx and ny
// must be at least two) with the given origin and spacing. Every rectangular
// cell is split into two counterclockwise triangles. It panics for degenerate
// dimensions.
func GridMesh(nx, ny int, x0, y0, dx, dy float64) *Mesh {
	if nx < 2 || ny < 2 {
		panic("meshgen: grid mesh needs at least 2x2 nodes")
	}
	verts := make([]Vec2, 0, nx*ny)
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			verts = append(verts, Vec2{x0 + float64(i)*dx, y0 + float64(j)*dy})
		}
	}
	idx := func(i, j int) int { return j*nx + i }
	var tris []Tri
	for j := 0; j < ny-1; j++ {
		for i := 0; i < nx-1; i++ {
			bl := idx(i, j)
			br := idx(i+1, j)
			tr := idx(i+1, j+1)
			tl := idx(i, j+1)
			tris = append(tris, Tri{bl, br, tr}, Tri{bl, tr, tl})
		}
	}
	return &Mesh{Vertices: verts, Triangles: tris}
}

// GridPoints returns the nodes of a regular nx by ny lattice in row-major
// order.
func GridPoints(nx, ny int, x0, y0, dx, dy float64) []Vec2 {
	pts := make([]Vec2, 0, nx*ny)
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			pts = append(pts, Vec2{x0 + float64(i)*dx, y0 + float64(j)*dy})
		}
	}
	return pts
}

// TriangulatePolygonFan returns a triangle-fan triangulation of a convex
// polygon given by its ordered vertices, as a Mesh. The polygon must have at
// least three vertices; the result is only a valid triangulation for convex
// polygons.
func TriangulatePolygonFan(poly []Vec2) (*Mesh, error) {
	if len(poly) < 3 {
		return nil, ErrNotEnoughPoints
	}
	verts := make([]Vec2, len(poly))
	copy(verts, poly)
	var tris []Tri
	for i := 1; i < len(poly)-1; i++ {
		t := Tri{0, i, i + 1}
		if Orient2D(verts[0], verts[i], verts[i+1]) < 0 {
			t = Tri{0, i + 1, i}
		}
		tris = append(tris, t)
	}
	return &Mesh{Vertices: verts, Triangles: tris}, nil
}

// TriangulatePolygonEarClip returns an ear-clipping triangulation of a simple
// (possibly non-convex) polygon given in counterclockwise or clockwise order.
// It returns ErrDegenerate if the polygon cannot be triangulated (for example
// if it self-intersects).
func TriangulatePolygonEarClip(poly []Vec2) (*Mesh, error) {
	n := len(poly)
	if n < 3 {
		return nil, ErrNotEnoughPoints
	}
	verts := make([]Vec2, n)
	copy(verts, poly)
	// Work with a CCW copy of the index ring.
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	if PolygonSignedArea(verts) < 0 {
		for i, j := 0, len(idx)-1; i < j; i, j = i+1, j-1 {
			idx[i], idx[j] = idx[j], idx[i]
		}
	}
	var tris []Tri
	guard := 0
	for len(idx) > 3 {
		guard++
		if guard > 10*n*n {
			return nil, ErrDegenerate
		}
		earFound := false
		m := len(idx)
		for k := 0; k < m; k++ {
			ip := idx[(k+m-1)%m]
			ic := idx[k]
			in := idx[(k+1)%m]
			a, b, c := verts[ip], verts[ic], verts[in]
			if Orient2D(a, b, c) <= 0 {
				continue // reflex or collinear
			}
			if pointInsideEar(verts, idx, ip, ic, in, a, b, c) {
				continue
			}
			tris = append(tris, Tri{ip, ic, in})
			idx = append(idx[:k], idx[k+1:]...)
			earFound = true
			break
		}
		if !earFound {
			return nil, ErrDegenerate
		}
	}
	tris = append(tris, Tri{idx[0], idx[1], idx[2]})
	return &Mesh{Vertices: verts, Triangles: tris}, nil
}

func pointInsideEar(verts []Vec2, idx []int, ip, ic, in int, a, b, c Vec2) bool {
	for _, vi := range idx {
		if vi == ip || vi == ic || vi == in {
			continue
		}
		if PointInTriangle(verts[vi], a, b, c, 1e-12) {
			return true
		}
	}
	return false
}
