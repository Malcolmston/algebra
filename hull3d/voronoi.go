package hull3d

import (
	"errors"
	"sort"
)

// VoronoiCell is the Voronoi region of a single site: the set of points closer
// to Site than to any other site, represented by its bounding half-spaces and,
// when bounded, its polytope. Cells on the outer boundary of the diagram are
// unbounded and have a nil Polytope.
type VoronoiCell struct {
	Site       int
	HalfSpaces []HalfSpace
	Polytope   *Polytope
	Bounded    bool
}

// VoronoiDiagram is the 3-D Voronoi diagram of a set of sites, the geometric
// dual of the Delaunay tetrahedralization.
type VoronoiDiagram struct {
	Sites []Vec3
	Cells []VoronoiCell
}

// VoronoiCells builds the Voronoi diagram of the given sites as the dual of the
// Delaunay tetrahedralization: each interior Delaunay vertex maps to a Voronoi
// cell whose vertices are the circumcentres of the incident tetrahedra. Cells on
// the convex hull are reported as unbounded. It returns an error if the sites are
// coplanar or too few.
//
// Each cell is additionally given by the bisector half-spaces separating its
// site from the sites it shares a Delaunay edge with; these define the cell
// exactly even when it is unbounded.
func VoronoiCells(sites []Vec3) (*VoronoiDiagram, error) {
	del, err := Delaunay(sites)
	if err != nil {
		return nil, err
	}
	pts := del.Points
	n := len(pts)

	// Circumcentre of each tetrahedron (a Voronoi vertex).
	circ := make([]Vec3, len(del.Tets))
	for i, t := range del.Tets {
		a, b, c, d := del.TetVertices(t)
		cc, _, err := Circumcenter(a, b, c, d)
		if err != nil {
			cc = a.Add(b).Add(c).Add(d).Scale(0.25)
		}
		circ[i] = cc
	}

	// For each site, gather incident tetrahedra and Delaunay neighbours.
	incident := make([][]int, n)
	neighbours := make([]map[int]bool, n)
	for i := range neighbours {
		neighbours[i] = map[int]bool{}
	}
	for ti, t := range del.Tets {
		vs := [4]int{t.A, t.B, t.C, t.D}
		for _, v := range vs {
			incident[v] = append(incident[v], ti)
			for _, w := range vs {
				if w != v {
					neighbours[v][w] = true
				}
			}
		}
	}

	// Sites on the convex hull have unbounded cells.
	onHull := make([]bool, n)
	for _, f := range del.BoundaryFaces() {
		onHull[f[0]] = true
		onHull[f[1]] = true
		onHull[f[2]] = true
	}

	diag := &VoronoiDiagram{Sites: pts, Cells: make([]VoronoiCell, n)}
	for s := 0; s < n; s++ {
		var hs []HalfSpace
		nbrs := make([]int, 0, len(neighbours[s]))
		for w := range neighbours[s] {
			nbrs = append(nbrs, w)
		}
		sort.Ints(nbrs)
		for _, w := range nbrs {
			mid := pts[s].Midpoint(pts[w])
			normal := pts[w].Sub(pts[s]) // outward: cell is where s is closer
			hs = append(hs, HalfSpaceFromPointNormal(mid, normal))
		}
		cell := VoronoiCell{Site: s, HalfSpaces: hs}
		if !onHull[s] {
			// Bounded interior cell: vertices are incident circumcentres.
			verts := make([]Vec3, 0, len(incident[s]))
			for _, ti := range incident[s] {
				verts = append(verts, circ[ti])
			}
			verts = DedupPoints(verts, 1e-9)
			if len(verts) >= 4 {
				if poly, err := ConvexHull(verts); err == nil {
					cell.Polytope = poly
					cell.Bounded = true
				}
			}
		}
		diag.Cells[s] = cell
	}
	return diag, nil
}

// Cell returns the Voronoi cell for site index i.
func (v *VoronoiDiagram) Cell(i int) (VoronoiCell, error) {
	if i < 0 || i >= len(v.Cells) {
		return VoronoiCell{}, errors.New("hull3d: site index out of range")
	}
	return v.Cells[i], nil
}

// NearestSite returns the index of the site nearest to q by brute-force search,
// which is the Voronoi cell containing q. It returns an error if there are no
// sites.
func (v *VoronoiDiagram) NearestSite(q Vec3) (int, error) {
	if len(v.Sites) == 0 {
		return -1, errors.New("hull3d: no sites")
	}
	best := 0
	bd := v.Sites[0].DistanceSq(q)
	for i := 1; i < len(v.Sites); i++ {
		if d := v.Sites[i].DistanceSq(q); d < bd {
			bd, best = d, i
		}
	}
	return best, nil
}

// Contains reports whether q lies in the Voronoi cell of its site, i.e. it
// satisfies all of the cell's bisector half-spaces within tolerance eps.
func (c VoronoiCell) Contains(q Vec3, eps float64) bool {
	for _, hs := range c.HalfSpaces {
		if !hs.Contains(q, eps) {
			return false
		}
	}
	return true
}

// NumBoundedCells returns the number of bounded (interior) Voronoi cells.
func (v *VoronoiDiagram) NumBoundedCells() int {
	c := 0
	for _, cell := range v.Cells {
		if cell.Bounded {
			c++
		}
	}
	return c
}
