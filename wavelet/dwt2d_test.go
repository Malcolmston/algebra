package wavelet

import (
	"math"
	"testing"
)

// makeImage builds a deterministic rows-by-cols test image.
func makeImage(rows, cols int) [][]float64 {
	m := make([][]float64, rows)
	for r := 0; r < rows; r++ {
		m[r] = make([]float64, cols)
		for c := 0; c < cols; c++ {
			m[r][c] = math.Sin(0.2*float64(r)) + math.Cos(0.3*float64(c)) + float64((r*c)%5)
		}
	}
	return m
}

func imageApproxEqual(a, b [][]float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for r := range a {
		if !sliceApproxEqual(a[r], b[r], tol) {
			return false
		}
	}
	return true
}

func TestDWT2DPerfectReconstruction(t *testing.T) {
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		for _, dims := range [][2]int{{4, 4}, {8, 6}, {16, 8}} {
			img := makeImage(dims[0], dims[1])
			res := DWT2D(img, w)
			rr, cc := res.Dims()
			if rr != dims[0]/2 || cc != dims[1]/2 {
				t.Errorf("%s %v: subband dims = (%d,%d)", w.Name(), dims, rr, cc)
			}
			back := IDWT2D(res, w)
			if !imageApproxEqual(back, img, 1e-9) {
				t.Errorf("%s %v: 2D reconstruction failed", w.Name(), dims)
			}
		}
	}
}

func TestDWT2DEnergyConservation(t *testing.T) {
	img := makeImage(8, 8)
	var in float64
	for _, row := range img {
		in += Energy(row)
	}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		res := DWT2D(img, w)
		if !approxEqual(res.Energy(), in, 1e-9) {
			t.Errorf("%s: 2D subband energy = %v want %v", w.Name(), res.Energy(), in)
		}
	}
}

func TestDWT2DConstantImage(t *testing.T) {
	// A constant image concentrates all energy in LL; the detail subbands vanish.
	rows, cols := 8, 8
	img := makeImage(rows, cols)
	for r := range img {
		for c := range img[r] {
			img[r][c] = 4
		}
	}
	res := DWT2D(img, DB2())
	for _, sb := range [][][]float64{res.LH, res.HL, res.HH} {
		for _, row := range sb {
			for _, v := range row {
				if math.Abs(v) > 1e-9 {
					t.Errorf("detail subband nonzero on constant image: %v", v)
				}
			}
		}
	}
}

func TestDWT2DPanicsOnRagged(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("ragged image should panic")
		}
	}()
	DWT2D([][]float64{{1, 2}, {3}}, Haar())
}
