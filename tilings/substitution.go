package tilings

import (
	"math"
)

// Affine2 is a general affine map of the plane, x -> M x + t, with a 2x2 linear
// part that need not be orthogonal (so it may include scaling and shear). It is
// used to place substitution tiles.
type Affine2 struct {
	M00, M01, M10, M11 float64
	Tx, Ty             float64
}

// IdentityAffine returns the identity affine map.
func IdentityAffine() Affine2 { return Affine2{M00: 1, M11: 1} }

// ScaleAffine returns the uniform scaling by factor s about the origin.
func ScaleAffine(s float64) Affine2 { return Affine2{M00: s, M11: s} }

// TranslateAffine returns the translation by (dx, dy).
func TranslateAffine(dx, dy float64) Affine2 {
	return Affine2{M00: 1, M11: 1, Tx: dx, Ty: dy}
}

// FromIsometry lifts an isometry to an affine map.
func FromIsometry(a Isometry) Affine2 {
	return Affine2{M00: a.M00, M01: a.M01, M10: a.M10, M11: a.M11, Tx: a.Tx, Ty: a.Ty}
}

// Apply returns the image of the point p under the affine map.
func (a Affine2) Apply(p Point) Point {
	return Point{a.M00*p.X + a.M01*p.Y + a.Tx, a.M10*p.X + a.M11*p.Y + a.Ty}
}

// Compose returns the affine map a∘inner (inner applied first).
func (a Affine2) Compose(inner Affine2) Affine2 {
	return Affine2{
		M00: a.M00*inner.M00 + a.M01*inner.M10,
		M01: a.M00*inner.M01 + a.M01*inner.M11,
		M10: a.M10*inner.M00 + a.M11*inner.M10,
		M11: a.M10*inner.M01 + a.M11*inner.M11,
		Tx:  a.M00*inner.Tx + a.M01*inner.Ty + a.Tx,
		Ty:  a.M10*inner.Tx + a.M11*inner.Ty + a.Ty,
	}
}

// Determinant returns the determinant of the linear part; its absolute value is
// the area scaling factor of the map.
func (a Affine2) Determinant() float64 { return a.M00*a.M11 - a.M01*a.M10 }

// ------------------------------------------------------------------
// Penrose tilings via Robinson triangles.
// ------------------------------------------------------------------

// RobinsonKind labels the two Robinson (half-rhombus) triangle types used by the
// Penrose deflation.
type RobinsonKind int

const (
	// Thin is a half of a thin Penrose rhombus.
	Thin RobinsonKind = iota
	// Thick is a half of a thick Penrose rhombus.
	Thick
)

// String returns the name of the Robinson triangle kind.
func (k RobinsonKind) String() string {
	if k == Thin {
		return "thin"
	}
	return "thick"
}

// RobinsonTriangle is an oriented Robinson triangle with vertices A, B and C. It
// is one half of a Penrose rhombus; two triangles of the same kind sharing their
// long edge form a rhombus.
type RobinsonTriangle struct {
	Kind    RobinsonKind
	A, B, C Point
}

// Vertices returns the three vertices A, B, C of the triangle.
func (t RobinsonTriangle) Vertices() []Point { return []Point{t.A, t.B, t.C} }

// Area returns the unsigned area of the triangle.
func (t RobinsonTriangle) Area() float64 {
	return math.Abs(t.B.Sub(t.A).Cross(t.C.Sub(t.A))) / 2
}

// Centroid returns the centroid of the triangle.
func (t RobinsonTriangle) Centroid() Point { return Centroid(t.A, t.B, t.C) }

// DeflateRobinson applies one step of the Penrose deflation (with inflation
// factor Phi) to each triangle, returning the finer set of triangles. This is
// the canonical Robinson-triangle subdivision underlying both the P2 and P3
// Penrose tilings.
func DeflateRobinson(tris []RobinsonTriangle) []RobinsonTriangle {
	out := make([]RobinsonTriangle, 0, len(tris)*3)
	for _, t := range tris {
		if t.Kind == Thin {
			p := t.A.Add(t.B.Sub(t.A).Scale(1 / Phi))
			out = append(out,
				RobinsonTriangle{Thin, t.C, p, t.B},
				RobinsonTriangle{Thick, p, t.C, t.A},
			)
		} else {
			q := t.B.Add(t.A.Sub(t.B).Scale(1 / Phi))
			r := t.B.Add(t.C.Sub(t.B).Scale(1 / Phi))
			out = append(out,
				RobinsonTriangle{Thick, r, t.C, t.A},
				RobinsonTriangle{Thick, q, r, t.B},
				RobinsonTriangle{Thin, r, q, t.A},
			)
		}
	}
	return out
}

// DeflateRobinsonN applies n rounds of Penrose deflation.
func DeflateRobinsonN(tris []RobinsonTriangle, n int) []RobinsonTriangle {
	for i := 0; i < n; i++ {
		tris = DeflateRobinson(tris)
	}
	return tris
}

// PenroseSun returns the "sun" seed patch: ten Robinson triangles arranged with
// five-fold symmetry around the origin, each with circumradius r. Deflating this
// patch generates a Penrose tiling of a disc.
func PenroseSun(r float64) []RobinsonTriangle {
	tris := make([]RobinsonTriangle, 0, 10)
	for i := 0; i < 10; i++ {
		b := Polar(r, float64(2*i-1)*math.Pi/10)
		c := Polar(r, float64(2*i+1)*math.Pi/10)
		if i%2 == 0 {
			b, c = c, b
		}
		tris = append(tris, RobinsonTriangle{Thin, Origin(), b, c})
	}
	return tris
}

// PenroseP3 returns the Robinson triangles of a Penrose P3 (rhombus) tiling
// obtained by deflating the sun seed of radius r a total of n times.
func PenroseP3(n int, r float64) []RobinsonTriangle {
	return DeflateRobinsonN(PenroseSun(r), n)
}

// PenroseP2 returns the Robinson triangles of a Penrose P2 (kite and dart)
// tiling obtained by deflating the sun seed of radius r a total of n times. The
// P2 and P3 tilings share the same Robinson-triangle deflation and differ only
// in how the triangles are grouped into tiles (see PenroseRhombi and
// PenroseKitesDarts).
func PenroseP2(n int, r float64) []RobinsonTriangle {
	return DeflateRobinsonN(PenroseSun(r), n)
}

// Rhombus is a Penrose rhombus formed from two Robinson triangles of the same
// kind sharing their long edge. Vertices are given in order around the boundary.
type Rhombus struct {
	Kind     RobinsonKind
	Vertices [4]Point
}

// PenroseRhombi assembles Robinson triangles into Penrose rhombi (the P3 tiles)
// by pairing triangles of the same kind that share their longest edge. Triangles
// without a partner within the patch are dropped.
func PenroseRhombi(tris []RobinsonTriangle) []Rhombus {
	type key struct {
		kind   RobinsonKind
		mx, my int64
	}
	byEdge := map[key][]RobinsonTriangle{}
	for _, t := range tris {
		mid := longestEdgeMidpoint(t)
		k := key{t.Kind, roundKey(mid.X), roundKey(mid.Y)}
		byEdge[k] = append(byEdge[k], t)
	}
	var out []Rhombus
	for k, group := range byEdge {
		if len(group) < 2 {
			continue
		}
		t1, t2 := group[0], group[1]
		a1, b1 := longestEdge(t1)
		apex1 := oppositeVertex(t1, a1, b1)
		apex2 := oppositeVertex(t2, a1, b1)
		out = append(out, Rhombus{Kind: k.kind, Vertices: [4]Point{a1, apex1, b1, apex2}})
	}
	return out
}

// PenroseKitesDarts is a synonym assembly returning the Robinson triangles
// grouped into rhombi; the P2 kite/dart tiles are locally derivable from these.
func PenroseKitesDarts(tris []RobinsonTriangle) []Rhombus {
	return PenroseRhombi(tris)
}

func roundKey(x float64) int64 { return int64(math.Round(x * 1e6)) }

func longestEdge(t RobinsonTriangle) (Point, Point) {
	ab := t.A.Distance2(t.B)
	bc := t.B.Distance2(t.C)
	ca := t.C.Distance2(t.A)
	if ab >= bc && ab >= ca {
		return t.A, t.B
	}
	if bc >= ab && bc >= ca {
		return t.B, t.C
	}
	return t.C, t.A
}

func longestEdgeMidpoint(t RobinsonTriangle) Point {
	p, q := longestEdge(t)
	return p.Midpoint(q)
}

func oppositeVertex(t RobinsonTriangle, p, q Point) Point {
	for _, v := range t.Vertices() {
		if !v.ApproxEqual(p, 1e-6) && !v.ApproxEqual(q, 1e-6) {
			return v
		}
	}
	return t.A
}

// ------------------------------------------------------------------
// Chair (L-tromino) rep-tile.
// ------------------------------------------------------------------

// chairCanonical is the six-vertex outline of the canonical L-tromino: the
// 2x2 square with the top-right unit cell removed.
var chairCanonical = []Point{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}

// chairLocalMaps place four copies of the canonical L-tromino into the twice-
// inflated tile; composed with a half-scaling they realise the deflation.
var chairLocalMaps = []Affine2{
	{M00: 1, M01: 0, M10: 0, M11: 1, Tx: 0, Ty: 0},  // identity
	{M00: 1, M01: 0, M10: 0, M11: 1, Tx: 1, Ty: 1},  // translate (1,1)
	{M00: 0, M01: -1, M10: 1, M11: 0, Tx: 4, Ty: 0}, // rotate +90, shift
	{M00: 0, M01: 1, M10: -1, M11: 0, Tx: 0, Ty: 4}, // rotate -90, shift
}

// Chair is one tile of the chair (L-tromino) substitution tiling, represented
// by the affine Frame that maps the canonical L-tromino into the plane.
type Chair struct {
	Frame Affine2
}

// NewChair returns the canonical chair tile (identity frame).
func NewChair() Chair { return Chair{Frame: IdentityAffine()} }

// Polygon returns the six boundary vertices of the chair tile in order.
func (c Chair) Polygon() []Point {
	out := make([]Point, len(chairCanonical))
	for i, p := range chairCanonical {
		out[i] = c.Frame.Apply(p)
	}
	return out
}

// Area returns the area of the chair tile (three times the square of its linear
// scale).
func (c Chair) Area() float64 {
	return math.Abs(c.Frame.Determinant()) * 3
}

// Substitute returns the four half-scale chair tiles that subdivide c.
func (c Chair) Substitute() []Chair {
	half := ScaleAffine(0.5)
	out := make([]Chair, 4)
	for i, lm := range chairLocalMaps {
		out[i] = Chair{Frame: c.Frame.Compose(half.Compose(lm))}
	}
	return out
}

// ChairSubstitute applies one deflation step to a set of chair tiles.
func ChairSubstitute(chairs []Chair) []Chair {
	out := make([]Chair, 0, len(chairs)*4)
	for _, c := range chairs {
		out = append(out, c.Substitute()...)
	}
	return out
}

// ChairTiling returns the chair tiles obtained by deflating the canonical chair
// n times.
func ChairTiling(n int) []Chair {
	chairs := []Chair{NewChair()}
	for i := 0; i < n; i++ {
		chairs = ChairSubstitute(chairs)
	}
	return chairs
}

// ------------------------------------------------------------------
// Pinwheel tiling.
// ------------------------------------------------------------------

// PinwheelTriangle is a right triangle similar to the (1, 2, sqrt5) prototile of
// the pinwheel tiling, with the right angle at A, the long leg ending at B and
// the short leg ending at C.
type PinwheelTriangle struct {
	A, B, C Point
}

// CanonicalPinwheel returns the prototile with the right angle at the origin,
// the long leg of length 2 along the x-axis and the short leg of length 1 along
// the y-axis.
func CanonicalPinwheel() PinwheelTriangle {
	return PinwheelTriangle{A: Point{0, 0}, B: Point{2, 0}, C: Point{0, 1}}
}

// Vertices returns the three vertices A, B, C.
func (t PinwheelTriangle) Vertices() []Point { return []Point{t.A, t.B, t.C} }

// Area returns the unsigned area of the triangle.
func (t PinwheelTriangle) Area() float64 {
	return math.Abs(t.B.Sub(t.A).Cross(t.C.Sub(t.A))) / 2
}

// Centroid returns the centroid of the triangle.
func (t PinwheelTriangle) Centroid() Point { return Centroid(t.A, t.B, t.C) }

// pinwheelNormalize builds a PinwheelTriangle from three points forming a
// (1,2,sqrt5)-similar right triangle, detecting the right-angle vertex and
// ordering the legs so that B ends the long leg and C the short leg.
func pinwheelNormalize(p, q, r Point) PinwheelTriangle {
	pts := [3]Point{p, q, r}
	// Find the right-angle vertex.
	for i := 0; i < 3; i++ {
		v := pts[i]
		u := pts[(i+1)%3]
		w := pts[(i+2)%3]
		du := u.Sub(v)
		dw := w.Sub(v)
		if math.Abs(du.Dot(dw)) <= 1e-9*du.Norm()*dw.Norm()+1e-12 {
			if du.Norm() >= dw.Norm() {
				return PinwheelTriangle{A: v, B: u, C: w}
			}
			return PinwheelTriangle{A: v, B: w, C: u}
		}
	}
	// Fallback: pick the vertex opposite the longest edge as the right angle.
	return PinwheelTriangle{A: p, B: q, C: r}
}

// Substitute returns the five sub-triangles, each similar to the prototile at
// scale 1/sqrt5, that exactly tile t. The prototile is split by the altitude
// from the right angle to the hypotenuse, and the larger of the two pieces is
// subdivided into four congruent triangles.
func (t PinwheelTriangle) Substitute() []PinwheelTriangle {
	a, b, c := t.A, t.B, t.C
	// Foot of the altitude from A onto line BC.
	cb := c.Sub(b)
	s := a.Sub(b).Dot(cb) / cb.Norm2()
	f := b.Add(cb.Scale(s))
	// Small triangle A,F,C (contains the short leg's continuation).
	small := pinwheelNormalize(a, f, c)
	// Large triangle A,F,B subdivided by its midpoints into four copies.
	mAF := a.Midpoint(f)
	mFB := f.Midpoint(b)
	mAB := a.Midpoint(b)
	q1 := pinwheelNormalize(a, mAF, mAB)
	q2 := pinwheelNormalize(mAF, f, mFB)
	q3 := pinwheelNormalize(mAB, mFB, b)
	q4 := pinwheelNormalize(mAF, mFB, mAB)
	return []PinwheelTriangle{small, q1, q2, q3, q4}
}

// PinwheelSubstitute applies one deflation step to a set of pinwheel triangles.
func PinwheelSubstitute(tris []PinwheelTriangle) []PinwheelTriangle {
	out := make([]PinwheelTriangle, 0, len(tris)*5)
	for _, t := range tris {
		out = append(out, t.Substitute()...)
	}
	return out
}

// PinwheelTiling returns the pinwheel triangles obtained by deflating the
// canonical prototile n times.
func PinwheelTiling(n int) []PinwheelTriangle {
	tris := []PinwheelTriangle{CanonicalPinwheel()}
	for i := 0; i < n; i++ {
		tris = PinwheelSubstitute(tris)
	}
	return tris
}
