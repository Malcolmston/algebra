package chaos

import (
	"math"
	"math/cmplx"
)

// StabilityKind classifies the linear stability of a fixed point or
// equilibrium.
type StabilityKind int

const (
	// Unknown indicates the classification could not be determined.
	Unknown StabilityKind = iota
	// StableNode: all eigenvalues real, negative (flow) / inside unit circle.
	StableNode
	// UnstableNode: all eigenvalues real, positive / outside unit circle.
	UnstableNode
	// Saddle: eigenvalues of mixed stability.
	Saddle
	// StableFocus: complex eigenvalues with contracting spiral.
	StableFocus
	// UnstableFocus: complex eigenvalues with expanding spiral.
	UnstableFocus
	// Center: purely imaginary eigenvalues / eigenvalues on the unit circle.
	Center
)

// String returns a human-readable name for the stability kind.
func (k StabilityKind) String() string {
	switch k {
	case StableNode:
		return "stable node"
	case UnstableNode:
		return "unstable node"
	case Saddle:
		return "saddle"
	case StableFocus:
		return "stable focus"
	case UnstableFocus:
		return "unstable focus"
	case Center:
		return "center"
	default:
		return "unknown"
	}
}

// FixedPoint1D refines a fixed point of the one-dimensional map f near the
// guess x0 by Newton's method applied to g(x)=f(x)-x, using a numerical
// derivative. It returns ErrNoConvergence if the tolerance is not met.
func FixedPoint1D(f Map1D, x0, tol float64, maxIter int) (float64, error) {
	x := x0
	const h = 1e-7
	for i := 0; i < maxIter; i++ {
		g := f(x) - x
		if math.Abs(g) < tol {
			return x, nil
		}
		dg := (f(x+h) - f(x-h)) / (2 * h) // f'
		dg -= 1                           // g' = f' - 1
		if dg == 0 {
			break
		}
		x -= g / dg
	}
	if math.Abs(f(x)-x) < tol {
		return x, nil
	}
	return x, ErrNoConvergence
}

// LogisticFixedPoints returns the two analytic fixed points of the logistic
// map with parameter r: 0 and 1 - 1/r.
func LogisticFixedPoints(r float64) []float64 {
	if r == 0 {
		return []float64{0}
	}
	return []float64{0, 1 - 1/r}
}

// Multiplier1D returns the derivative (multiplier) of the map f at the fixed
// point x, estimated by a central difference.
func Multiplier1D(f Map1D, x float64) float64 {
	const h = 1e-7
	return (f(x+h) - f(x-h)) / (2 * h)
}

// IsStable1D reports whether the fixed point x of the map f is linearly stable,
// i.e. whether the magnitude of its multiplier is less than one.
func IsStable1D(f Map1D, x float64) bool {
	return math.Abs(Multiplier1D(f, x)) < 1
}

// JacobianMap returns the numerical Jacobian of the n-dimensional map F at the
// point x, computed by central differences with step eps.
func JacobianMap(F MapN, x Vec, eps float64) Mat {
	n := len(x)
	fx := F(x)
	J := NewMat(len(fx), n)
	for j := 0; j < n; j++ {
		xp := x.Clone()
		xm := x.Clone()
		xp[j] += eps
		xm[j] -= eps
		fp := F(xp)
		fm := F(xm)
		for i := range fp {
			J[i][j] = (fp[i] - fm[i]) / (2 * eps)
		}
	}
	return J
}

// JacobianField returns the numerical Jacobian of the vector field f at x,
// computed by central differences with step eps.
func JacobianField(f Field, x Vec, eps float64) Mat {
	return JacobianMap(MapN(f), x, eps)
}

// FixedPointMap refines a fixed point of the n-dimensional map F near x0 by
// Newton's method on G(x)=F(x)-x with a numerical Jacobian.
func FixedPointMap(F MapN, x0 Vec, tol float64, maxIter int) (Vec, error) {
	x := x0.Clone()
	n := len(x)
	for it := 0; it < maxIter; it++ {
		g := F(x).Sub(x)
		if g.Norm() < tol {
			return x, nil
		}
		J := JacobianMap(F, x, 1e-7)
		for i := 0; i < n; i++ {
			J[i][i] -= 1 // Jacobian of G = J_F - I
		}
		delta, err := SolveLinear(J, g)
		if err != nil {
			return x, err
		}
		x = x.Sub(delta)
	}
	if F(x).Sub(x).Norm() < tol {
		return x, nil
	}
	return x, ErrNoConvergence
}

// Equilibrium refines an equilibrium of the vector field f (a root of f(x)=0)
// near x0 by Newton's method with a numerical Jacobian.
func Equilibrium(f Field, x0 Vec, tol float64, maxIter int) (Vec, error) {
	x := x0.Clone()
	for it := 0; it < maxIter; it++ {
		g := f(x)
		if g.Norm() < tol {
			return x, nil
		}
		J := JacobianField(f, x, 1e-7)
		delta, err := SolveLinear(J, g)
		if err != nil {
			return x, err
		}
		x = x.Sub(delta)
	}
	if f(x).Norm() < tol {
		return x, nil
	}
	return x, ErrNoConvergence
}

// LorenzEquilibria returns the analytic equilibria of the Lorenz system with
// parameters rho and beta: the origin, and (for rho>1) the symmetric pair
// C+/-.
func LorenzEquilibria(rho, beta float64) []Vec {
	eq := []Vec{{0, 0, 0}}
	if rho > 1 {
		s := math.Sqrt(beta * (rho - 1))
		eq = append(eq, Vec{s, s, rho - 1}, Vec{-s, -s, rho - 1})
	}
	return eq
}

// ClassifyFlow classifies the equilibrium of a continuous system from the
// eigenvalues of its Jacobian (2-by-2 or 3-by-3), by the sign of the real
// parts and the presence of an imaginary part.
func ClassifyFlow(J Mat) StabilityKind {
	eigs := smallEigs(J)
	if eigs == nil {
		return Unknown
	}
	var nPos, nNeg, nZero int
	hasImag := false
	for _, l := range eigs {
		re := real(l)
		if math.Abs(imag(l)) > 1e-9*(1+cmplx.Abs(l)) {
			hasImag = true
		}
		switch {
		case re > 1e-12:
			nPos++
		case re < -1e-12:
			nNeg++
		default:
			nZero++
		}
	}
	if nPos > 0 && nNeg > 0 {
		return Saddle
	}
	if nZero == len(eigs) || (nPos == 0 && nNeg == 0) {
		return Center
	}
	if hasImag {
		if nPos == 0 {
			return StableFocus
		}
		return UnstableFocus
	}
	if nPos == 0 {
		return StableNode
	}
	return UnstableNode
}

// ClassifyMap classifies the fixed point of a discrete system from the
// eigenvalues (multipliers) of its Jacobian, by their magnitude relative to
// the unit circle.
func ClassifyMap(J Mat) StabilityKind {
	eigs := smallEigs(J)
	if eigs == nil {
		return Unknown
	}
	var nIn, nOut, nOn int
	hasImag := false
	for _, l := range eigs {
		m := cmplx.Abs(l)
		if math.Abs(imag(l)) > 1e-9*(1+m) {
			hasImag = true
		}
		switch {
		case m > 1+1e-9:
			nOut++
		case m < 1-1e-9:
			nIn++
		default:
			nOn++
		}
	}
	if nIn > 0 && nOut > 0 {
		return Saddle
	}
	if nOn == len(eigs) {
		return Center
	}
	if hasImag {
		if nOut == 0 {
			return StableFocus
		}
		return UnstableFocus
	}
	if nOut == 0 {
		return StableNode
	}
	return UnstableNode
}

// smallEigs returns the eigenvalues of a 1-, 2- or 3-dimensional square
// matrix, or nil for other sizes.
func smallEigs(J Mat) []complex128 {
	switch J.Rows() {
	case 1:
		return []complex128{complex(J[0][0], 0)}
	case 2:
		e := Eigenvalues2(J)
		return []complex128{e[0], e[1]}
	case 3:
		e := Eigenvalues3(J)
		return []complex128{e[0], e[1], e[2]}
	default:
		return nil
	}
}

// IsHyperbolicFlow reports whether every eigenvalue of J has non-zero real
// part (a hyperbolic equilibrium of a flow).
func IsHyperbolicFlow(J Mat) bool {
	for _, l := range smallEigs(J) {
		if math.Abs(real(l)) <= 1e-12 {
			return false
		}
	}
	return true
}

// IsHyperbolicMap reports whether every eigenvalue of J has magnitude
// different from one (a hyperbolic fixed point of a map).
func IsHyperbolicMap(J Mat) bool {
	for _, l := range smallEigs(J) {
		if math.Abs(cmplx.Abs(l)-1) <= 1e-12 {
			return false
		}
	}
	return true
}
