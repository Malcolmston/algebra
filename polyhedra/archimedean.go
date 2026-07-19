package polyhedra

import "math"

// FaceType records that a solid has Count regular faces each with Sides sides.
type FaceType struct {
	Sides int
	Count int
}

// ArchimedeanSolid captures the combinatorial and metric data of one of the
// thirteen semiregular convex polyhedra (uniform, with regular polygon faces and
// a single vertex figure but more than one kind of face). Volume and
// circumradius are stored as coefficients for unit edge length; surface area is
// derived exactly from the regular face polygons.
type ArchimedeanSolid struct {
	Name       string     // common name
	VertexFig  string     // vertex configuration, e.g. "3.6.6"
	NumV       int        // number of vertices
	NumE       int        // number of edges
	NumF       int        // number of faces
	Faces      []FaceType // regular faces by polygon type
	VolumeUnit float64    // volume at edge length 1
	CircumUnit float64    // circumradius at edge length 1
}

// EulerCharacteristic returns V - E + F, which is 2 for every Archimedean solid.
func (s ArchimedeanSolid) EulerCharacteristic() int { return s.NumV - s.NumE + s.NumF }

// SurfaceArea returns the surface area of the solid at edge length a, summed
// exactly over its regular polygonal faces.
func (s ArchimedeanSolid) SurfaceArea(a float64) float64 {
	var area float64
	for _, ft := range s.Faces {
		area += float64(ft.Count) * RegularPolygonArea(ft.Sides, a)
	}
	return area
}

// Volume returns the volume of the solid at edge length a.
func (s ArchimedeanSolid) Volume(a float64) float64 { return s.VolumeUnit * a * a * a }

// Circumradius returns the circumradius (center-to-vertex distance) at edge
// length a.
func (s ArchimedeanSolid) Circumradius(a float64) float64 { return s.CircumUnit * a }

// Midradius returns the midradius at edge length a, derived from the
// circumradius and edge via ρ = sqrt(R² - a²/4). This holds for the uniform
// Archimedean solids, whose every edge subtends the same midsphere.
func (s ArchimedeanSolid) Midradius(a float64) float64 {
	r := s.Circumradius(a)
	m := r*r - a*a/4
	if m < 0 {
		return 0
	}
	return math.Sqrt(m)
}

// TotalFaceCount returns the number of faces implied by the Faces breakdown; it
// equals NumF for a consistent record and is useful as a self-check.
func (s ArchimedeanSolid) TotalFaceCount() int {
	n := 0
	for _, ft := range s.Faces {
		n += ft.Count
	}
	return n
}

// ---- The thirteen Archimedean solids ----

// NewTruncatedTetrahedron returns the truncated tetrahedron (vertex figure 3.6.6).
func NewTruncatedTetrahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"truncated tetrahedron", "3.6.6", 12, 18, 8,
		[]FaceType{{3, 4}, {6, 4}}, 23 * math.Sqrt2 / 12, math.Sqrt(22) / 4}
}

// NewCuboctahedron returns the cuboctahedron (vertex figure 3.4.3.4).
func NewCuboctahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"cuboctahedron", "3.4.3.4", 12, 24, 14,
		[]FaceType{{3, 8}, {4, 6}}, 5 * math.Sqrt2 / 3, 1}
}

// NewTruncatedCube returns the truncated cube (vertex figure 3.8.8).
func NewTruncatedCube() ArchimedeanSolid {
	return ArchimedeanSolid{"truncated cube", "3.8.8", 24, 36, 14,
		[]FaceType{{3, 8}, {8, 6}}, (21 + 14*math.Sqrt2) / 3, math.Sqrt(7+4*math.Sqrt2) / 2}
}

// NewTruncatedOctahedron returns the truncated octahedron (vertex figure 4.6.6).
func NewTruncatedOctahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"truncated octahedron", "4.6.6", 24, 36, 14,
		[]FaceType{{4, 6}, {6, 8}}, 8 * math.Sqrt2, math.Sqrt(10) / 2}
}

// NewRhombicuboctahedron returns the rhombicuboctahedron (vertex figure 3.4.4.4).
func NewRhombicuboctahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"rhombicuboctahedron", "3.4.4.4", 24, 48, 26,
		[]FaceType{{3, 8}, {4, 18}}, (12 + 10*math.Sqrt2) / 3, math.Sqrt(5+2*math.Sqrt2) / 2}
}

// NewTruncatedCuboctahedron returns the truncated cuboctahedron, also called the
// great rhombicuboctahedron (vertex figure 4.6.8).
func NewTruncatedCuboctahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"truncated cuboctahedron", "4.6.8", 48, 72, 26,
		[]FaceType{{4, 12}, {6, 8}, {8, 6}}, 22 + 14*math.Sqrt2, math.Sqrt(13+6*math.Sqrt2) / 2}
}

// NewSnubCube returns the snub cube (vertex figure 3.3.3.3.4). Its volume and
// circumradius involve the tribonacci constant and are given numerically.
func NewSnubCube() ArchimedeanSolid {
	return ArchimedeanSolid{"snub cube", "3.3.3.3.4", 24, 60, 38,
		[]FaceType{{3, 32}, {4, 6}}, 7.8894773999, 1.3437133737}
}

// NewIcosidodecahedron returns the icosidodecahedron (vertex figure 3.5.3.5).
func NewIcosidodecahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"icosidodecahedron", "3.5.3.5", 30, 60, 32,
		[]FaceType{{3, 20}, {5, 12}}, (45 + 17*math.Sqrt(5)) / 6, Phi}
}

// NewTruncatedDodecahedron returns the truncated dodecahedron (vertex figure 3.10.10).
func NewTruncatedDodecahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"truncated dodecahedron", "3.10.10", 60, 90, 32,
		[]FaceType{{3, 20}, {10, 12}}, (5.0 / 12) * (99 + 47*math.Sqrt(5)), math.Sqrt(74+30*math.Sqrt(5)) / 4}
}

// NewTruncatedIcosahedron returns the truncated icosahedron, the classic
// football/soccer-ball shape (vertex figure 5.6.6).
func NewTruncatedIcosahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"truncated icosahedron", "5.6.6", 60, 90, 32,
		[]FaceType{{5, 12}, {6, 20}}, (125 + 43*math.Sqrt(5)) / 4, math.Sqrt(58+18*math.Sqrt(5)) / 4}
}

// NewRhombicosidodecahedron returns the rhombicosidodecahedron (vertex figure 3.4.5.4).
func NewRhombicosidodecahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"rhombicosidodecahedron", "3.4.5.4", 60, 120, 62,
		[]FaceType{{3, 20}, {4, 30}, {5, 12}}, (60 + 29*math.Sqrt(5)) / 3, math.Sqrt(11+4*math.Sqrt(5)) / 2}
}

// NewTruncatedIcosidodecahedron returns the truncated icosidodecahedron, also
// called the great rhombicosidodecahedron (vertex figure 4.6.10).
func NewTruncatedIcosidodecahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"truncated icosidodecahedron", "4.6.10", 120, 180, 62,
		[]FaceType{{4, 30}, {6, 20}, {10, 12}}, 95 + 50*math.Sqrt(5), math.Sqrt(31+12*math.Sqrt(5)) / 2}
}

// NewSnubDodecahedron returns the snub dodecahedron (vertex figure 3.3.3.3.5).
// Its volume and circumradius involve algebraic constants and are given
// numerically.
func NewSnubDodecahedron() ArchimedeanSolid {
	return ArchimedeanSolid{"snub dodecahedron", "3.3.3.3.5", 60, 150, 92,
		[]FaceType{{3, 80}, {5, 12}}, 37.6166499018, 2.1558373751}
}

// ArchimedeanSolids returns all thirteen Archimedean solids in a canonical
// order (increasing face count, then vertices).
func ArchimedeanSolids() []ArchimedeanSolid {
	return []ArchimedeanSolid{
		NewTruncatedTetrahedron(),
		NewCuboctahedron(),
		NewTruncatedCube(),
		NewTruncatedOctahedron(),
		NewRhombicuboctahedron(),
		NewTruncatedCuboctahedron(),
		NewSnubCube(),
		NewIcosidodecahedron(),
		NewTruncatedDodecahedron(),
		NewTruncatedIcosahedron(),
		NewRhombicosidodecahedron(),
		NewTruncatedIcosidodecahedron(),
		NewSnubDodecahedron(),
	}
}
