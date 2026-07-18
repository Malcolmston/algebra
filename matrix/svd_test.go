package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

// --- test helpers ------------------------------------------------------------

func svdFloatsOf(t *testing.T, m *Matrix) [][]float64 {
	t.Helper()
	f, err := m.Floats()
	if err != nil {
		t.Fatalf("Floats: %v", err)
	}
	return f
}

// svdReconstruct forms U·diag(s)·Vᵀ from the thin factors (U is m×k, V is n×k).
func svdReconstruct(u [][]float64, s []float64, v [][]float64) [][]float64 {
	rows := len(u)
	cols := len(v)
	k := len(s)
	out := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		out[i] = make([]float64, cols)
		for j := 0; j < cols; j++ {
			var acc float64
			for c := 0; c < k; c++ {
				acc += u[i][c] * s[c] * v[j][c]
			}
			out[i][j] = acc
		}
	}
	return out
}

func svdApproxEqual(a, b [][]float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		for j := range a[i] {
			if math.Abs(a[i][j]-b[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// svdGram returns Mᵀ·M for an r×c matrix M given as a flat [][]float64.
func svdGram(m [][]float64) [][]float64 {
	rows := len(m)
	cols := 0
	if rows > 0 {
		cols = len(m[0])
	}
	g := make([][]float64, cols)
	for i := 0; i < cols; i++ {
		g[i] = make([]float64, cols)
		for j := 0; j < cols; j++ {
			var acc float64
			for r := 0; r < rows; r++ {
				acc += m[r][i] * m[r][j]
			}
			g[i][j] = acc
		}
	}
	return g
}

func svdIsIdentity(g [][]float64, tol float64) bool {
	for i := range g {
		for j := range g[i] {
			want := 0.0
			if i == j {
				want = 1.0
			}
			if math.Abs(g[i][j]-want) > tol {
				return false
			}
		}
	}
	return true
}

const svdTol = 1e-9

// --- tests -------------------------------------------------------------------

func TestSingularValues(t *testing.T) {
	cases := []struct {
		name string
		in   [][]float64
		want []float64
	}{
		{"diag2", [][]float64{{3, 0}, {0, 2}}, []float64{3, 2}},
		{"signed", [][]float64{{2, 0}, {0, -3}}, []float64{3, 2}},
		{"wide", [][]float64{{1, 0, 0}, {0, 1, 0}}, []float64{1, 1}},
		{"tall", [][]float64{{1, 0}, {0, 1}, {0, 0}}, []float64{1, 1}},
		{"rank1", [][]float64{{1, 2}, {2, 4}}, []float64{5, 0}},
		{"golden", [][]float64{{1, 1}, {0, 1}}, []float64{
			(1 + math.Sqrt(5)) / 2, (math.Sqrt(5) - 1) / 2,
		}},
		{"deficient3", [][]float64{{5, 0, 0}, {0, 4, 0}, {0, 0, 0}}, []float64{5, 4, 0}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := FromFloats(tc.in)
			s, err := m.SingularValues()
			if err != nil {
				t.Fatalf("SingularValues: %v", err)
			}
			if len(s) != len(tc.want) {
				t.Fatalf("len(s) = %d, want %d", len(s), len(tc.want))
			}
			for i := range s {
				if math.Abs(s[i]-tc.want[i]) > svdTol {
					t.Errorf("s[%d] = %v, want %v", i, s[i], tc.want[i])
				}
				if i > 0 && s[i] > s[i-1]+svdTol {
					t.Errorf("singular values not non-increasing: %v", s)
				}
			}
		})
	}
}

func TestSVDReconstructAndOrthonormal(t *testing.T) {
	cases := [][][]float64{
		{{3, 0}, {0, 2}},
		{{2, 0}, {0, -3}},
		{{1, 2}, {3, 4}},
		{{1, 0, 0}, {0, 1, 0}},   // wide (transpose path)
		{{1, 2}, {3, 4}, {5, 6}}, // tall
		{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}},
		{{1, 2}, {2, 4}}, // rank deficient
	}
	for idx, in := range cases {
		m := FromFloats(in)
		u, s, v, err := m.SVD()
		if err != nil {
			t.Fatalf("case %d: SVD: %v", idx, err)
		}
		uf := svdFloatsOf(t, u)
		vf := svdFloatsOf(t, v)
		got := svdReconstruct(uf, s, vf)
		if !svdApproxEqual(got, in, 1e-8) {
			t.Errorf("case %d: reconstruction mismatch\n got %v\nwant %v", idx, got, in)
		}
		// V always has orthonormal columns.
		if !svdIsIdentity(svdGram(vf), 1e-8) {
			t.Errorf("case %d: Vᵀ·V not identity", idx)
		}
		// Sign convention: largest-magnitude entry of each U column non-negative
		// (columns with nonzero singular value).
		k := len(s)
		for c := 0; c < k; c++ {
			if s[c] <= svdTol {
				continue
			}
			maxAbs, maxVal := 0.0, 0.0
			for r := range uf {
				if av := math.Abs(uf[r][c]); av > maxAbs {
					maxAbs = av
					maxVal = uf[r][c]
				}
			}
			if maxVal < 0 {
				t.Errorf("case %d: U column %d violates sign convention", idx, c)
			}
		}
	}
}

func TestRankNumeric(t *testing.T) {
	cases := []struct {
		name string
		in   [][]float64
		want int
	}{
		{"full2", [][]float64{{3, 0}, {0, 2}}, 2},
		{"rank1", [][]float64{{1, 2}, {2, 4}}, 1},
		{"zero", [][]float64{{0, 0}, {0, 0}}, 0},
		{"deficient3", [][]float64{{5, 0, 0}, {0, 4, 0}, {0, 0, 0}}, 2},
		{"wideFull", [][]float64{{1, 0, 0}, {0, 1, 0}}, 2},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := FromFloats(tc.in)
			r, err := m.RankNumeric(0)
			if err != nil {
				t.Fatalf("RankNumeric: %v", err)
			}
			if r != tc.want {
				t.Errorf("rank = %d, want %d", r, tc.want)
			}
		})
	}
}

func TestNullspace(t *testing.T) {
	// Rank-deficient square: one-dimensional kernel spanned by (2,-1).
	m := FromFloats([][]float64{{1, 2}, {2, 4}})
	ns, err := m.Nullspace(0)
	if err != nil {
		t.Fatalf("Nullspace: %v", err)
	}
	if ns.Cols() != 1 || ns.Rows() != 2 {
		t.Fatalf("nullspace shape = %dx%d, want 2x1", ns.Rows(), ns.Cols())
	}
	nf := svdFloatsOf(t, ns)
	// Column must be a unit vector.
	norm := math.Hypot(nf[0][0], nf[1][0])
	if math.Abs(norm-1) > svdTol {
		t.Errorf("nullspace column not unit: norm %v", norm)
	}
	// A·n must be ~0.
	af := svdFloatsOf(t, m)
	for i := 0; i < 2; i++ {
		v := af[i][0]*nf[0][0] + af[i][1]*nf[1][0]
		if math.Abs(v) > 1e-8 {
			t.Errorf("A·n component %d = %v, want 0", i, v)
		}
	}

	// Full-rank matrix: empty basis.
	full := FromFloats([][]float64{{3, 0}, {0, 2}})
	nsf, err := full.Nullspace(0)
	if err != nil {
		t.Fatalf("Nullspace(full): %v", err)
	}
	if nsf.Cols() != 0 {
		t.Errorf("full-rank nullspace cols = %d, want 0", nsf.Cols())
	}
}

func TestPinv(t *testing.T) {
	cases := [][][]float64{
		{{3, 0}, {0, 2}},
		{{1, 2}, {3, 4}},
		{{1, 0, 0}, {0, 1, 0}}, // wide
		{{1, 2}, {2, 4}},       // rank deficient
		{{1, 2}, {3, 4}, {5, 6}},
	}
	for idx, in := range cases {
		m := FromFloats(in)
		p, err := m.Pinv(0)
		if err != nil {
			t.Fatalf("case %d: Pinv: %v", idx, err)
		}
		af := svdFloatsOf(t, m)
		pf := svdFloatsOf(t, p)
		// Moore-Penrose property: A·A⁺·A == A.
		ap := svdMatMul(af, pf)
		apa := svdMatMul(ap, af)
		if !svdApproxEqual(apa, af, 1e-7) {
			t.Errorf("case %d: A·A⁺·A != A", idx)
		}
		// A⁺·A·A⁺ == A⁺.
		pa := svdMatMul(pf, af)
		pap := svdMatMul(pa, pf)
		if !svdApproxEqual(pap, pf, 1e-7) {
			t.Errorf("case %d: A⁺·A·A⁺ != A⁺", idx)
		}
	}
}

func svdMatMul(a, b [][]float64) [][]float64 {
	ra := len(a)
	ca := 0
	if ra > 0 {
		ca = len(a[0])
	}
	cb := 0
	if len(b) > 0 {
		cb = len(b[0])
	}
	out := make([][]float64, ra)
	for i := 0; i < ra; i++ {
		out[i] = make([]float64, cb)
		for j := 0; j < cb; j++ {
			var acc float64
			for k := 0; k < ca; k++ {
				acc += a[i][k] * b[k][j]
			}
			out[i][j] = acc
		}
	}
	return out
}

func TestCond2(t *testing.T) {
	m := FromFloats([][]float64{{3, 0}, {0, 2}})
	c, err := m.Cond2()
	if err != nil {
		t.Fatalf("Cond2: %v", err)
	}
	if math.Abs(c-1.5) > svdTol {
		t.Errorf("cond = %v, want 1.5", c)
	}

	// Rank-deficient: infinite condition number.
	sing := FromFloats([][]float64{{1, 2}, {2, 4}})
	c, err = sing.Cond2()
	if err != nil {
		t.Fatalf("Cond2(singular): %v", err)
	}
	if !math.IsInf(c, 1) {
		t.Errorf("cond = %v, want +Inf", c)
	}
}

func TestSVDDeterministic(t *testing.T) {
	in := [][]float64{{4, 1, 2}, {1, 3, 0}, {2, 0, 5}}
	m := FromFloats(in)
	_, s1, _, err := m.SVD()
	if err != nil {
		t.Fatalf("SVD: %v", err)
	}
	_, s2, _, err := m.SVD()
	if err != nil {
		t.Fatalf("SVD: %v", err)
	}
	for i := range s1 {
		if s1[i] != s2[i] {
			t.Errorf("non-deterministic singular value %d: %v vs %v", i, s1[i], s2[i])
		}
	}
}

func TestSVDSymbolic(t *testing.T) {
	m := FromExpr([][]algebra.Expr{
		{algebra.Sym("x"), algebra.Int(1)},
		{algebra.Int(1), algebra.Int(1)},
	})
	if _, _, _, err := m.SVD(); !errors.Is(err, ErrUnsupported) {
		t.Errorf("SVD symbolic err = %v, want ErrUnsupported", err)
	}
	if _, err := m.SingularValues(); !errors.Is(err, ErrUnsupported) {
		t.Errorf("SingularValues symbolic err = %v, want ErrUnsupported", err)
	}
	if _, err := m.Pinv(0); !errors.Is(err, ErrUnsupported) {
		t.Errorf("Pinv symbolic err = %v, want ErrUnsupported", err)
	}
	if _, err := m.RankNumeric(0); !errors.Is(err, ErrUnsupported) {
		t.Errorf("RankNumeric symbolic err = %v, want ErrUnsupported", err)
	}
}

// --- benchmarks --------------------------------------------------------------

func svdBenchMatrix(n int) *Matrix {
	vals := make([][]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			// Deterministic, well-mixed entries.
			vals[i][j] = math.Sin(float64(i*7+j*3+1)) + float64((i*j)%5)
		}
	}
	return FromFloats(vals)
}

func BenchmarkSVD(b *testing.B) {
	m := svdBenchMatrix(24)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, _, err := m.SVD(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkPinv(b *testing.B) {
	m := svdBenchMatrix(24)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.Pinv(0); err != nil {
			b.Fatal(err)
		}
	}
}
