package tilings

import (
	"math"
	"sort"
	"strings"
)

// InteriorAngle returns the interior angle of a regular n-gon in radians,
// (n-2)*pi/n. It panics for n < 3.
func InteriorAngle(n int) float64 {
	if n < 3 {
		panic("tilings: InteriorAngle requires n >= 3")
	}
	return float64(n-2) * math.Pi / float64(n)
}

// InteriorAngleDeg returns the interior angle of a regular n-gon in degrees.
func InteriorAngleDeg(n int) float64 { return Rad2Deg(InteriorAngle(n)) }

// ExteriorAngle returns the exterior angle of a regular n-gon in radians,
// 2*pi/n.
func ExteriorAngle(n int) float64 {
	if n < 3 {
		panic("tilings: ExteriorAngle requires n >= 3")
	}
	return 2 * math.Pi / float64(n)
}

// PolygonArea2 returns the area of a regular n-gon with circumradius r.
func PolygonAreaRegular(n int, r float64) float64 {
	if n < 3 {
		panic("tilings: PolygonAreaRegular requires n >= 3")
	}
	return 0.5 * float64(n) * r * r * math.Sin(2*math.Pi/float64(n))
}

// PolygonAreaFromSide returns the area of a regular n-gon with side length s.
func PolygonAreaFromSide(n int, s float64) float64 {
	if n < 3 {
		panic("tilings: PolygonAreaFromSide requires n >= 3")
	}
	return float64(n) * s * s / (4 * math.Tan(math.Pi/float64(n)))
}

// Circumradius returns the circumradius of a regular n-gon with side length s.
func Circumradius(n int, s float64) float64 {
	return s / (2 * math.Sin(math.Pi/float64(n)))
}

// Inradius returns the inradius (apothem) of a regular n-gon with side length s.
func Inradius(n int, s float64) float64 {
	return s / (2 * math.Tan(math.Pi/float64(n)))
}

// VertexConfiguration is the cyclic sequence of polygon sizes meeting at a
// vertex of a tiling; for example {4, 8, 8} for the truncated square tiling.
type VertexConfiguration []int

// AngleSum returns the sum of the interior angles of the polygons in the
// configuration, in radians.
func (v VertexConfiguration) AngleSum() float64 {
	var s float64
	for _, n := range v {
		s += InteriorAngle(n)
	}
	return s
}

// AngleSumDeg returns the interior-angle sum in degrees.
func (v VertexConfiguration) AngleSumDeg() float64 { return Rad2Deg(v.AngleSum()) }

// IsPlanar reports whether the configuration's interior angles sum to exactly
// 2*pi (360 degrees) to within eps, the condition for a flat vertex.
func (v VertexConfiguration) IsPlanar(eps float64) bool {
	return math.Abs(v.AngleSum()-2*math.Pi) <= eps
}

// Degree returns the number of polygons meeting at the vertex.
func (v VertexConfiguration) Degree() int { return len(v) }

// IsRegular reports whether the configuration consists of three or more copies
// of a single polygon type meeting at a flat vertex (a regular tiling vertex).
func (v VertexConfiguration) IsRegular(eps float64) bool {
	if len(v) < 3 || !v.IsPlanar(eps) {
		return false
	}
	for _, n := range v {
		if n != v[0] {
			return false
		}
	}
	return true
}

// Multiset returns the polygon sizes in sorted order, discarding the cyclic
// arrangement.
func (v VertexConfiguration) Multiset() []int {
	out := append([]int(nil), v...)
	sort.Ints(out)
	return out
}

// String returns the configuration in dotted notation, e.g. "3.4.6.4".
func (v VertexConfiguration) String() string {
	parts := make([]string, len(v))
	for i, n := range v {
		parts[i] = itoa(n)
	}
	return strings.Join(parts, ".")
}

// EqualMultiset reports whether two configurations use the same polygons the
// same number of times, ignoring order.
func (v VertexConfiguration) EqualMultiset(w VertexConfiguration) bool {
	a, b := v.Multiset(), w.Multiset()
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// UniformTiling names a regular or semiregular (Archimedean) tiling of the
// plane by its vertex configuration.
type UniformTiling struct {
	// Name is the common name of the tiling.
	Name string
	// Config is the vertex configuration.
	Config VertexConfiguration
	// Regular reports whether the tiling is one of the three regular tilings.
	Regular bool
}

// String returns the tiling's name and configuration.
func (t UniformTiling) String() string { return t.Name + " (" + t.Config.String() + ")" }

var regularTilings = []UniformTiling{
	{"triangular", VertexConfiguration{3, 3, 3, 3, 3, 3}, true},
	{"square", VertexConfiguration{4, 4, 4, 4}, true},
	{"hexagonal", VertexConfiguration{6, 6, 6}, true},
}

var semiregularTilings = []UniformTiling{
	{"trihexagonal", VertexConfiguration{3, 6, 3, 6}, false},
	{"snub square", VertexConfiguration{3, 3, 4, 3, 4}, false},
	{"snub trihexagonal", VertexConfiguration{3, 3, 3, 3, 6}, false},
	{"elongated triangular", VertexConfiguration{3, 3, 3, 4, 4}, false},
	{"rhombitrihexagonal", VertexConfiguration{3, 4, 6, 4}, false},
	{"truncated square", VertexConfiguration{4, 8, 8}, false},
	{"truncated hexagonal", VertexConfiguration{3, 12, 12}, false},
	{"truncated trihexagonal", VertexConfiguration{4, 6, 12}, false},
}

// RegularTilings returns the three regular tilings of the plane (triangular,
// square, hexagonal).
func RegularTilings() []UniformTiling {
	return append([]UniformTiling(nil), regularTilings...)
}

// SemiregularTilings returns the eight semiregular (Archimedean) tilings of the
// plane.
func SemiregularTilings() []UniformTiling {
	return append([]UniformTiling(nil), semiregularTilings...)
}

// ArchimedeanTilings returns all eleven uniform tilings of the plane (the three
// regular tilings followed by the eight semiregular ones).
func ArchimedeanTilings() []UniformTiling {
	out := make([]UniformTiling, 0, 11)
	out = append(out, regularTilings...)
	out = append(out, semiregularTilings...)
	return out
}

// TilingByName returns the uniform tiling with the given name and reports
// whether it was found.
func TilingByName(name string) (UniformTiling, bool) {
	for _, t := range ArchimedeanTilings() {
		if t.Name == name {
			return t, true
		}
	}
	return UniformTiling{}, false
}

// EnumerateVertexTypes returns every multiset of regular-polygon sizes (each at
// least 3) whose interior angles sum to exactly 2*pi, i.e. every combinatorial
// way that regular polygons can surround a vertex. The results are sorted
// non-decreasingly and returned in lexicographic order. There are 17 such
// vertex types.
func EnumerateVertexTypes() []VertexConfiguration {
	var out []VertexConfiguration
	// Angle contributed by an n-gon is (1 - 2/n) of pi; the sum over the vertex
	// must equal 2, i.e. sum(1 - 2/n_i) = 2. Enumerate non-decreasing sequences.
	const eps = 1e-9
	var rec func(start int, remaining float64, cur []int)
	rec = func(start int, remaining float64, cur []int) {
		if math.Abs(remaining) <= eps && len(cur) >= 3 {
			out = append(out, append(VertexConfiguration(nil), cur...))
			return
		}
		if remaining <= eps {
			return
		}
		for n := start; ; n++ {
			contrib := 1 - 2/float64(n)
			// Smallest remaining contribution uses the largest polygon; if even the
			// current polygon overshoots by too little to ever hit zero, stop.
			if contrib > remaining+eps {
				break
			}
			rec(n, remaining-contrib, append(cur, n))
			// An upper bound on n: a single remaining slot needs contrib==remaining;
			// beyond n where contrib exceeds remaining we already broke out.
			if n > 1000 {
				break
			}
		}
	}
	rec(3, 2, nil)
	sort.Slice(out, func(i, j int) bool {
		a, b := out[i], out[j]
		for k := 0; k < len(a) && k < len(b); k++ {
			if a[k] != b[k] {
				return a[k] < b[k]
			}
		}
		return len(a) < len(b)
	})
	return out
}

// IsArchimedeanConfig reports whether the given vertex configuration matches one
// of the eleven uniform tilings, comparing as a cyclic sequence up to rotation
// and reflection.
func IsArchimedeanConfig(v VertexConfiguration) bool {
	for _, t := range ArchimedeanTilings() {
		if cyclicMatch(v, t.Config) {
			return true
		}
	}
	return false
}

// cyclicMatch reports whether a equals b as a cyclic sequence, allowing
// rotation and reflection.
func cyclicMatch(a, b VertexConfiguration) bool {
	n := len(a)
	if n != len(b) {
		return false
	}
	if n == 0 {
		return true
	}
	try := func(seq []int) bool {
		for shift := 0; shift < n; shift++ {
			ok := true
			for i := 0; i < n; i++ {
				if seq[(i+shift)%n] != a[i] {
					ok = false
					break
				}
			}
			if ok {
				return true
			}
		}
		return false
	}
	if try(b) {
		return true
	}
	rev := make([]int, n)
	for i := 0; i < n; i++ {
		rev[i] = b[n-1-i]
	}
	return try(rev)
}

// DualTilingName returns the name of the dual (Laves/Catalan) tiling of a
// regular tiling: triangular and hexagonal are dual to each other and the square
// tiling is self-dual. It reports false for non-regular inputs.
func DualTilingName(name string) (string, bool) {
	switch name {
	case "triangular":
		return "hexagonal", true
	case "hexagonal":
		return "triangular", true
	case "square":
		return "square", true
	default:
		return "", false
	}
}
