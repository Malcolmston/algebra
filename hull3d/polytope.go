package hull3d

import (
	"errors"
	"fmt"
	"math"
	"sort"
)

// Face is a triangular face of a [Polytope], given as three indices into the
// polytope's vertex slice. For a closed convex polytope the vertices are ordered
// counter-clockwise as seen from outside, so that (B-A)×(C-A) is the outward
// normal.
type Face struct {
	A, B, C int
}

// Edge is an undirected edge between two vertex indices, stored with Lo <= Hi so
// that it can be used as a canonical map key.
type Edge struct {
	Lo, Hi int
}

// MakeEdge returns the canonical [Edge] joining vertex indices i and j.
func MakeEdge(i, j int) Edge {
	if i <= j {
		return Edge{i, j}
	}
	return Edge{j, i}
}

// Polytope is a boundary representation of a (typically convex) polyhedron: a
// list of vertices together with the oriented triangular faces of its surface.
type Polytope struct {
	Vertices []Vec3
	Faces    []Face
}

// NewPolytope returns a polytope with the given vertices and faces. The faces
// are assumed to reference valid vertex indices; use [Polytope.Validate] to
// check.
func NewPolytope(vertices []Vec3, faces []Face) *Polytope {
	return &Polytope{Vertices: vertices, Faces: faces}
}

// NumVertices returns the number of vertices.
func (p *Polytope) NumVertices() int { return len(p.Vertices) }

// NumFaces returns the number of triangular faces.
func (p *Polytope) NumFaces() int { return len(p.Faces) }

// Validate reports an error if any face references an out-of-range vertex index
// or is degenerate (repeats a vertex).
func (p *Polytope) Validate() error {
	n := len(p.Vertices)
	for i, f := range p.Faces {
		if f.A < 0 || f.A >= n || f.B < 0 || f.B >= n || f.C < 0 || f.C >= n {
			return fmt.Errorf("hull3d: face %d references vertex out of range", i)
		}
		if f.A == f.B || f.B == f.C || f.A == f.C {
			return fmt.Errorf("hull3d: face %d is degenerate", i)
		}
	}
	return nil
}

// FaceVertices returns the three vertex positions of face f.
func (p *Polytope) FaceVertices(f Face) (a, b, c Vec3) {
	return p.Vertices[f.A], p.Vertices[f.B], p.Vertices[f.C]
}

// FaceNormal returns the (unnormalised) outward normal (B-A)×(C-A) of face f.
func (p *Polytope) FaceNormal(f Face) Vec3 {
	a, b, c := p.FaceVertices(f)
	return b.Sub(a).Cross(c.Sub(a))
}

// FaceArea returns the area of triangular face f.
func (p *Polytope) FaceArea(f Face) float64 {
	return 0.5 * p.FaceNormal(f).Length()
}

// FaceCentroid returns the centroid (barycentre) of triangular face f.
func (p *Polytope) FaceCentroid(f Face) Vec3 {
	a, b, c := p.FaceVertices(f)
	return a.Add(b).Add(c).Scale(1.0 / 3.0)
}

// FacePlane returns the plane of face f, oriented so its normal is the outward
// normal. It returns an error if the face is degenerate.
func (p *Polytope) FacePlane(f Face) (Plane, error) {
	a, b, c := p.FaceVertices(f)
	return PlaneFromPoints(a, b, c)
}

// SurfaceArea returns the total surface area of the polytope, the sum of its
// face areas.
func (p *Polytope) SurfaceArea() float64 {
	var s float64
	for _, f := range p.Faces {
		s += p.FaceArea(f)
	}
	return s
}

// SignedVolume returns the signed volume enclosed by the surface via the
// divergence theorem. For a correctly outward-oriented closed surface it is
// positive and equal to the enclosed volume.
func (p *Polytope) SignedVolume() float64 {
	var v float64
	for _, f := range p.Faces {
		a, b, c := p.FaceVertices(f)
		v += Triple(a, b, c)
	}
	return v / 6.0
}

// Volume returns the (unsigned) volume enclosed by the surface.
func (p *Polytope) Volume() float64 { return math.Abs(p.SignedVolume()) }

// Centroid returns the centroid of the solid region enclosed by the surface,
// computed by integrating position over the tetrahedra formed with the origin.
// It falls back to the vertex centroid when the enclosed volume is negligible.
func (p *Polytope) Centroid() Vec3 {
	var vol float64
	var c Vec3
	for _, f := range p.Faces {
		a, b, cc := p.FaceVertices(f)
		w := Triple(a, b, cc)
		vol += w
		c = c.Add(a.Add(b).Add(cc).Scale(w))
	}
	if math.Abs(vol) < 1e-300 {
		return Centroid(p.Vertices)
	}
	// tetra centroid is quarter-sum; overall centroid = sum(w*centroid)/sum(w)
	return c.Scale(1.0 / (4.0 * vol))
}

// Edges returns the distinct undirected edges of the polytope's surface.
func (p *Polytope) Edges() []Edge {
	seen := make(map[Edge]struct{})
	var out []Edge
	add := func(i, j int) {
		e := MakeEdge(i, j)
		if _, ok := seen[e]; !ok {
			seen[e] = struct{}{}
			out = append(out, e)
		}
	}
	for _, f := range p.Faces {
		add(f.A, f.B)
		add(f.B, f.C)
		add(f.C, f.A)
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Lo != out[j].Lo {
			return out[i].Lo < out[j].Lo
		}
		return out[i].Hi < out[j].Hi
	})
	return out
}

// NumEdges returns the number of distinct undirected edges.
func (p *Polytope) NumEdges() int { return len(p.Edges()) }

// EulerCharacteristic returns V - E + F. For a simple closed polyhedron
// (topologically a sphere) it equals 2.
func (p *Polytope) EulerCharacteristic() int {
	return p.NumVertices() - p.NumEdges() + p.NumFaces()
}

// IsClosed reports whether every undirected edge is shared by exactly two faces,
// a necessary condition for the surface to bound a solid.
func (p *Polytope) IsClosed() bool {
	count := make(map[Edge]int)
	for _, f := range p.Faces {
		count[MakeEdge(f.A, f.B)]++
		count[MakeEdge(f.B, f.C)]++
		count[MakeEdge(f.C, f.A)]++
	}
	for _, c := range count {
		if c != 2 {
			return false
		}
	}
	return len(count) > 0
}

// BoundingBox returns the axis-aligned bounding box of the vertices. It returns
// an error if the polytope has no vertices.
func (p *Polytope) BoundingBox() (min, max Vec3, err error) {
	return BoundingBox(p.Vertices)
}

// VertexCentroid returns the arithmetic mean of the vertices.
func (p *Polytope) VertexCentroid() Vec3 { return Centroid(p.Vertices) }

// FacePlanes returns the outward-oriented supporting plane of every face,
// skipping degenerate faces. The i-th returned plane corresponds to the i-th
// non-degenerate face.
func (p *Polytope) FacePlanes() []Plane {
	var out []Plane
	for _, f := range p.Faces {
		if pl, err := p.FacePlane(f); err == nil {
			out = append(out, pl)
		}
	}
	return out
}

// Contains reports whether point q lies inside or on the polytope, treating it
// as the intersection of the half-spaces below each face plane. The polytope
// must be convex and outward-oriented. A point within tolerance eps of a face is
// considered contained.
func (p *Polytope) Contains(q Vec3, eps float64) bool {
	for _, f := range p.Faces {
		pl, err := p.FacePlane(f)
		if err != nil {
			continue
		}
		if pl.SignedDistance(q) > eps {
			return false
		}
	}
	return true
}

// Support returns the vertex of the polytope farthest in direction d (a support
// point of the convex hull of its vertices). It returns an error if there are no
// vertices.
func (p *Polytope) Support(d Vec3) (Vec3, error) {
	_, v, err := FarthestPoint(p.Vertices, d)
	return v, err
}

// Translate returns a copy of the polytope with every vertex shifted by t.
func (p *Polytope) Translate(t Vec3) *Polytope {
	vs := make([]Vec3, len(p.Vertices))
	for i, v := range p.Vertices {
		vs[i] = v.Add(t)
	}
	fs := append([]Face(nil), p.Faces...)
	return &Polytope{Vertices: vs, Faces: fs}
}

// Scaled returns a copy of the polytope scaled about the origin by factor s.
// Faces are reversed when s is negative so that outward orientation is
// preserved.
func (p *Polytope) Scaled(s float64) *Polytope {
	vs := make([]Vec3, len(p.Vertices))
	for i, v := range p.Vertices {
		vs[i] = v.Scale(s)
	}
	fs := make([]Face, len(p.Faces))
	for i, f := range p.Faces {
		if s < 0 {
			fs[i] = Face{f.A, f.C, f.B}
		} else {
			fs[i] = f
		}
	}
	return &Polytope{Vertices: vs, Faces: fs}
}

// Diameter returns the maximum distance between any two vertices (the geometric
// diameter of the vertex set). It returns 0 for fewer than two vertices.
func (p *Polytope) Diameter() float64 {
	var d2 float64
	for i := 0; i < len(p.Vertices); i++ {
		for j := i + 1; j < len(p.Vertices); j++ {
			if s := p.Vertices[i].DistanceSq(p.Vertices[j]); s > d2 {
				d2 = s
			}
		}
	}
	return math.Sqrt(d2)
}

// OrientOutward re-orders each face so its normal points away from the polytope
// centroid, making the surface consistently outward-oriented. It requires a
// convex polytope with an interior. The receiver is modified in place and
// returned.
func (p *Polytope) OrientOutward() *Polytope {
	c := Centroid(p.Vertices)
	for i, f := range p.Faces {
		a, b, cc := p.FaceVertices(f)
		n := b.Sub(a).Cross(cc.Sub(a))
		if n.Dot(c.Sub(a)) > 0 {
			p.Faces[i] = Face{f.A, f.C, f.B}
		}
	}
	return p
}

// Copy returns a deep copy of the polytope.
func (p *Polytope) Copy() *Polytope {
	return &Polytope{
		Vertices: append([]Vec3(nil), p.Vertices...),
		Faces:    append([]Face(nil), p.Faces...),
	}
}

// TetrahedronVolume returns the unsigned volume of the tetrahedron with the
// given four corners.
func TetrahedronVolume(a, b, c, d Vec3) float64 {
	return math.Abs(Triple(b.Sub(a), c.Sub(a), d.Sub(a))) / 6.0
}

// TriangleArea returns the area of the triangle with the given three corners.
func TriangleArea(a, b, c Vec3) float64 {
	return 0.5 * b.Sub(a).Cross(c.Sub(a)).Length()
}

// errNoInterior is returned when a construction cannot produce a full 3-D body.
var errNoInterior = errors.New("hull3d: points do not span three dimensions")
