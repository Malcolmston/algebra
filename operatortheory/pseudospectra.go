package operatortheory

import (
	"math"
	"math/cmplx"
)

// SmallestSingularValue returns the smallest singular value of the matrix.
func (m *Matrix) SmallestSingularValue() float64 {
	s := m.SingularValues()
	if len(s) == 0 {
		return 0
	}
	return s[len(s)-1]
}

// shiftedResolventTarget returns z*I - m for a square matrix.
func (m *Matrix) shiftedResolventTarget(z complex128) *Matrix {
	n := m.rows
	a := m.Scale(-1)
	for i := 0; i < n; i++ {
		a.data[i*n+i] += z
	}
	return a
}

// ResolventNorm returns the spectral norm of the resolvent at z,
// ||(z*I - m)^{-1}||_2 = 1 / sigma_min(z*I - m). It returns +Inf when z is (very
// close to) an eigenvalue. This quantity defines the pseudospectra: z belongs to
// the eps-pseudospectrum precisely when ResolventNorm(z) >= 1/eps.
func (m *Matrix) ResolventNorm(z complex128) float64 {
	if !m.IsSquare() {
		return math.NaN()
	}
	smin := m.shiftedResolventTarget(z).SmallestSingularValue()
	if smin <= 0 {
		return math.Inf(1)
	}
	return 1 / smin
}

// InPseudospectrum reports whether z lies in the eps-pseudospectrum of m, i.e.
// whether ResolventNorm(z) >= 1/eps. It requires eps > 0.
func (m *Matrix) InPseudospectrum(z complex128, eps float64) bool {
	if eps <= 0 {
		return false
	}
	return m.ResolventNorm(z) >= 1/eps
}

// PseudospectrumGrid evaluates the resolvent norm on a rectangular grid in the
// complex plane spanning [reMin,reMax] x [imMin,imMax] with nx columns and ny
// rows. The result is indexed as grid[row][col], where row 0 corresponds to
// imMin. It returns ErrInvalidArgument for non-positive grid sizes.
func (m *Matrix) PseudospectrumGrid(reMin, reMax, imMin, imMax float64, nx, ny int) ([][]float64, error) {
	if nx <= 0 || ny <= 0 {
		return nil, ErrInvalidArgument
	}
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	grid := make([][]float64, ny)
	for r := 0; r < ny; r++ {
		grid[r] = make([]float64, nx)
		var y float64
		if ny == 1 {
			y = imMin
		} else {
			y = imMin + (imMax-imMin)*float64(r)/float64(ny-1)
		}
		for c := 0; c < nx; c++ {
			var x float64
			if nx == 1 {
				x = reMin
			} else {
				x = reMin + (reMax-reMin)*float64(c)/float64(nx-1)
			}
			grid[r][c] = m.ResolventNorm(complex(x, y))
		}
	}
	return grid, nil
}

// PseudospectralAbscissa estimates the eps-pseudospectral abscissa, the maximum
// real part of any point in the eps-pseudospectrum. It scans a grid around the
// spectrum enlarged by a margin proportional to eps and refines the rightmost
// crossing. It requires eps > 0 and returns ErrInvalidArgument otherwise.
func (m *Matrix) PseudospectralAbscissa(eps float64, samples int) (float64, error) {
	if eps <= 0 {
		return 0, ErrInvalidArgument
	}
	if samples < 16 {
		samples = 16
	}
	vals := eigenvaluesQR(m)
	if len(vals) == 0 {
		return 0, ErrEmpty
	}
	reMin, reMax := real(vals[0]), real(vals[0])
	imMin, imMax := imag(vals[0]), imag(vals[0])
	for _, v := range vals {
		reMin = math.Min(reMin, real(v))
		reMax = math.Max(reMax, real(v))
		imMin = math.Min(imMin, imag(v))
		imMax = math.Max(imMax, imag(v))
	}
	margin := 3 * eps * (1 + m.OperatorNorm())
	reMin -= margin
	reMax += margin
	imMin -= margin
	imMax += margin
	best := math.Inf(-1)
	inv := 1 / eps
	for r := 0; r < samples; r++ {
		y := imMin + (imMax-imMin)*float64(r)/float64(samples-1)
		for c := 0; c < samples; c++ {
			x := reMin + (reMax-reMin)*float64(c)/float64(samples-1)
			if m.ResolventNorm(complex(x, y)) >= inv {
				if x > best {
					best = x
				}
			}
		}
	}
	if math.IsInf(best, -1) {
		// Fall back to the spectral abscissa.
		return m.SpectralAbscissa(), nil
	}
	return best, nil
}

// PseudospectralRadius estimates the eps-pseudospectral radius, the largest
// modulus of any point in the eps-pseudospectrum, by scanning a grid enclosing
// the spectrum enlarged by a margin proportional to eps.
func (m *Matrix) PseudospectralRadius(eps float64, samples int) (float64, error) {
	if eps <= 0 {
		return 0, ErrInvalidArgument
	}
	if samples < 16 {
		samples = 16
	}
	vals := eigenvaluesQR(m)
	if len(vals) == 0 {
		return 0, ErrEmpty
	}
	var rad float64
	for _, v := range vals {
		if a := cmplx.Abs(v); a > rad {
			rad = a
		}
	}
	margin := 3 * eps * (1 + m.OperatorNorm())
	lim := rad + margin
	inv := 1 / eps
	best := 0.0
	for r := 0; r < samples; r++ {
		y := -lim + 2*lim*float64(r)/float64(samples-1)
		for c := 0; c < samples; c++ {
			x := -lim + 2*lim*float64(c)/float64(samples-1)
			if m.ResolventNorm(complex(x, y)) >= inv {
				if a := math.Hypot(x, y); a > best {
					best = a
				}
			}
		}
	}
	return best, nil
}

// DistanceToSingularity returns the 2-norm distance from a square matrix to the
// nearest singular matrix, which equals the smallest singular value. It is the
// backward-error interpretation of the resolvent norm evaluated at 0.
func (m *Matrix) DistanceToSingularity() float64 {
	return m.SmallestSingularValue()
}
