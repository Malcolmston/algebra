package voronoi

import (
	"math"
	"sort"
)

// VoronoiDiagram is the dual of a Delaunay triangulation: it partitions the
// plane into cells, one per input site, where each cell is the set of points
// closer to that site than to any other. Its finite vertices are the
// circumcentres of the Delaunay triangles.
type VoronoiDiagram struct {
	tri      *Triangulation
	vertices []Point // circumcentre of triangle i
	cells    []VoronoiCell
	edges    []VoronoiEdge
}

// VoronoiCell is the region of the plane assigned to a single site. Vertices
// lists the finite cell vertices in counterclockwise order. Unbounded cells
// (those whose site is on the convex hull) have Bounded == false; their listed
// vertices form the finite portion of the cell boundary.
type VoronoiCell struct {
	Site     int
	SiteXY   Point
	Vertices []Point
	Bounded  bool
}

// VoronoiEdge is a finite edge of the Voronoi diagram, joining the
// circumcentres of two adjacent Delaunay triangles. Tri1 and Tri2 are the two
// triangle indices, and Dual is the Delaunay edge (between the two shared
// sites) that the Voronoi edge crosses.
type VoronoiEdge struct {
	A, B       Point
	Tri1, Tri2 int
	Dual       Edge
}

// Voronoi builds the Voronoi diagram of the given sites by first constructing
// their Delaunay triangulation. It returns the same errors as Triangulate.
func Voronoi(sites []Point) (*VoronoiDiagram, error) {
	tri, err := Triangulate(sites)
	if err != nil {
		return nil, err
	}
	return VoronoiFromTriangulation(tri), nil
}

// VoronoiFromTriangulation constructs the Voronoi diagram dual to an existing
// Delaunay triangulation, avoiding a second triangulation pass.
func VoronoiFromTriangulation(tri *Triangulation) *VoronoiDiagram {
	vd := &VoronoiDiagram{tri: tri}
	vd.vertices = tri.Circumcenters()
	vd.buildEdges()
	vd.buildCells()
	return vd
}

// buildEdges computes the finite Voronoi edges from shared Delaunay edges.
func (vd *VoronoiDiagram) buildEdges() {
	tri := vd.tri
	edgeToTri := make(map[Edge][]int)
	for ti, t := range tri.triangles {
		for _, e := range t.Edges() {
			c := e.Canonical()
			edgeToTri[c] = append(edgeToTri[c], ti)
		}
	}
	keys := make([]Edge, 0, len(edgeToTri))
	for e := range edgeToTri {
		keys = append(keys, e)
	}
	sortEdges(keys)
	for _, e := range keys {
		tris := edgeToTri[e]
		if len(tris) == 2 {
			vd.edges = append(vd.edges, VoronoiEdge{
				A:    vd.vertices[tris[0]],
				B:    vd.vertices[tris[1]],
				Tri1: tris[0],
				Tri2: tris[1],
				Dual: e,
			})
		}
	}
}

// buildCells assembles a polygon for each site from the circumcentres of its
// incident triangles, ordered counterclockwise about the site.
func (vd *VoronoiDiagram) buildCells() {
	tri := vd.tri
	incident := tri.VertexTriangles()

	// Boundary sites (on the convex hull) yield unbounded cells.
	boundarySite := make([]bool, len(tri.points))
	for _, e := range tri.BoundaryEdges() {
		boundarySite[e.A] = true
		boundarySite[e.B] = true
	}

	vd.cells = make([]VoronoiCell, len(tri.points))
	for s := range tri.points {
		site := tri.points[s]
		tset := incident[s]
		verts := make([]Point, 0, len(tset))
		for _, ti := range tset {
			verts = append(verts, vd.vertices[ti])
		}
		// Order the vertices counterclockwise about the site.
		sort.Slice(verts, func(i, j int) bool {
			return site.AngleTo(verts[i]) < site.AngleTo(verts[j])
		})
		vd.cells[s] = VoronoiCell{
			Site:     s,
			SiteXY:   site,
			Vertices: verts,
			Bounded:  !boundarySite[s] && len(verts) >= 3,
		}
	}
}

// Triangulation returns the underlying Delaunay triangulation.
func (vd *VoronoiDiagram) Triangulation() *Triangulation { return vd.tri }

// Sites returns the diagram's sites (the deduplicated input points).
func (vd *VoronoiDiagram) Sites() []Point { return vd.tri.points }

// NumSites returns the number of sites.
func (vd *VoronoiDiagram) NumSites() int { return len(vd.tri.points) }

// Vertices returns the finite Voronoi vertices (the Delaunay circumcentres),
// indexed by triangle.
func (vd *VoronoiDiagram) Vertices() []Point { return vd.vertices }

// NumVertices returns the number of finite Voronoi vertices.
func (vd *VoronoiDiagram) NumVertices() int { return len(vd.vertices) }

// Edges returns the finite Voronoi edges.
func (vd *VoronoiDiagram) Edges() []VoronoiEdge { return vd.edges }

// NumEdges returns the number of finite Voronoi edges.
func (vd *VoronoiDiagram) NumEdges() int { return len(vd.edges) }

// Cells returns all Voronoi cells, indexed by site.
func (vd *VoronoiDiagram) Cells() []VoronoiCell { return vd.cells }

// Cell returns the Voronoi cell of site index s.
func (vd *VoronoiDiagram) Cell(s int) VoronoiCell { return vd.cells[s] }

// BoundedCells returns only the cells that are finite (bounded) polygons.
func (vd *VoronoiDiagram) BoundedCells() []VoronoiCell {
	out := make([]VoronoiCell, 0, len(vd.cells))
	for _, c := range vd.cells {
		if c.Bounded {
			out = append(out, c)
		}
	}
	return out
}

// DelaunayEdges returns the edges of the dual Delaunay triangulation.
func (vd *VoronoiDiagram) DelaunayEdges() []Edge { return vd.tri.Edges() }

// NearestSite returns the index of the site closest to p and the distance to
// it. It returns -1 for an empty diagram.
func (vd *VoronoiDiagram) NearestSite(p Point) (int, float64) {
	best := -1
	bestD := math.Inf(1)
	for i, s := range vd.tri.points {
		if d := s.DistanceSq(p); d < bestD {
			bestD = d
			best = i
		}
	}
	if best < 0 {
		return -1, 0
	}
	return best, math.Sqrt(bestD)
}

// Area returns the area of the cell. Unbounded cells return positive infinity.
func (c VoronoiCell) Area() float64 {
	if !c.Bounded {
		return math.Inf(1)
	}
	return PolygonArea(c.Vertices)
}

// Centroid returns the centroid of a bounded cell. For an unbounded or
// degenerate cell it returns the vertex average (or the site if there are no
// vertices).
func (c VoronoiCell) Centroid() Point {
	if len(c.Vertices) == 0 {
		return c.SiteXY
	}
	if !c.Bounded {
		return Centroid(c.Vertices)
	}
	return PolygonCentroid(c.Vertices)
}

// Perimeter returns the perimeter of the finite portion of the cell boundary.
func (c VoronoiCell) Perimeter() float64 {
	if c.Bounded {
		return PolygonPerimeter(c.Vertices)
	}
	// Open polyline for unbounded cells.
	var s float64
	for i := 0; i+1 < len(c.Vertices); i++ {
		s += c.Vertices[i].Distance(c.Vertices[i+1])
	}
	return s
}

// NumVertices returns the number of finite vertices of the cell.
func (c VoronoiCell) NumVertices() int { return len(c.Vertices) }

// Length returns the Euclidean length of the Voronoi edge.
func (e VoronoiEdge) Length() float64 { return e.A.Distance(e.B) }

// Midpoint returns the midpoint of the Voronoi edge.
func (e VoronoiEdge) Midpoint() Point { return e.A.Midpoint(e.B) }
