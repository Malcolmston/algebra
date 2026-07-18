package hypercomplex

import "math"

// SplitComplex is a split-complex (hyperbolic) number z = Re + Im*j, where the
// unit j satisfies j*j = +1 (in contrast to the ordinary imaginary unit).
// Split-complex numbers model Lorentz boosts in 1+1 spacetime dimensions.
type SplitComplex struct {
	Re, Im float64
}

// Split constructs a split-complex number Re + Im*j.
func Split(re, im float64) SplitComplex {
	return SplitComplex{Re: re, Im: im}
}

// Add returns the sum z + w.
func (z SplitComplex) Add(w SplitComplex) SplitComplex {
	return SplitComplex{z.Re + w.Re, z.Im + w.Im}
}

// Sub returns the difference z - w.
func (z SplitComplex) Sub(w SplitComplex) SplitComplex {
	return SplitComplex{z.Re - w.Re, z.Im - w.Im}
}

// Scale returns z with both parts multiplied by the real factor s.
func (z SplitComplex) Scale(s float64) SplitComplex {
	return SplitComplex{z.Re * s, z.Im * s}
}

// Neg returns the additive inverse -z.
func (z SplitComplex) Neg() SplitComplex {
	return SplitComplex{-z.Re, -z.Im}
}

// Mul returns the product z*w. Using j*j = +1,
// (a+bj)(c+dj) = (ac+bd) + (ad+bc)j. Multiplication is commutative.
func (z SplitComplex) Mul(w SplitComplex) SplitComplex {
	return SplitComplex{
		Re: z.Re*w.Re + z.Im*w.Im,
		Im: z.Re*w.Im + z.Im*w.Re,
	}
}

// Conj returns the conjugate Re - Im*j.
func (z SplitComplex) Conj() SplitComplex {
	return SplitComplex{z.Re, -z.Im}
}

// ModulusSq returns the split-complex squared modulus Re^2 - Im^2, equal to the
// real part of z * conj(z). It is the quadratic form preserved by
// split-complex multiplication and may be negative or zero.
func (z SplitComplex) ModulusSq() float64 {
	return z.Re*z.Re - z.Im*z.Im
}

// Abs returns sqrt(|ModulusSq|), a non-negative magnitude that is multiplicative
// in absolute value.
func (z SplitComplex) Abs() float64 {
	return math.Sqrt(math.Abs(z.ModulusSq()))
}

// IsZeroDivisor reports whether z lies on the null cone Re = ±Im (with z not the
// zero element) to within the absolute tolerance tol, in which case z has no
// multiplicative inverse.
func (z SplitComplex) IsZeroDivisor(tol float64) bool {
	if math.Abs(z.Re) <= tol && math.Abs(z.Im) <= tol {
		return false
	}
	return math.Abs(z.ModulusSq()) <= tol
}

// Inverse returns the multiplicative inverse conj(z)/ModulusSq, together with a
// boolean that is false when z is a zero divisor (ModulusSq = 0), in which case
// the inverse does not exist and the zero element is returned.
func (z SplitComplex) Inverse() (SplitComplex, bool) {
	m := z.ModulusSq()
	if m == 0 {
		return SplitComplex{}, false
	}
	return z.Conj().Scale(1 / m), true
}

// Equal reports whether z and w agree in both parts to within the absolute
// tolerance tol.
func (z SplitComplex) Equal(w SplitComplex, tol float64) bool {
	return math.Abs(z.Re-w.Re) <= tol && math.Abs(z.Im-w.Im) <= tol
}

// Exp returns the split-complex exponential exp(z) = e^Re (cosh Im + j sinh Im).
func (z SplitComplex) Exp() SplitComplex {
	e := math.Exp(z.Re)
	return SplitComplex{
		Re: e * math.Cosh(z.Im),
		Im: e * math.Sinh(z.Im),
	}
}

// Argument returns the hyperbolic argument (rapidity) atanh(Im/Re) of z. It is
// well defined only in the right time-like sector Re > |Im|; elsewhere the
// result may be non-finite.
func (z SplitComplex) Argument() float64 {
	return math.Atanh(z.Im / z.Re)
}

// SplitFromModulusArgument returns the time-like split-complex number with the
// given non-negative modulus r and rapidity phi, namely
// r*(cosh phi + j sinh phi).
func SplitFromModulusArgument(r, phi float64) SplitComplex {
	return SplitComplex{
		Re: r * math.Cosh(phi),
		Im: r * math.Sinh(phi),
	}
}

// LightCone returns the coordinates (u, v) of z in the null (diagonal) basis,
// with u = Re + Im and v = Re - Im. In this basis multiplication acts
// component-wise: (u1,v1)*(u2,v2) = (u1*u2, v1*v2).
func (z SplitComplex) LightCone() (u, v float64) {
	return z.Re + z.Im, z.Re - z.Im
}

// SplitFromLightCone reconstructs a split-complex number from its null-basis
// coordinates (u, v), returning ((u+v)/2) + ((u-v)/2)*j.
func SplitFromLightCone(u, v float64) SplitComplex {
	return SplitComplex{Re: (u + v) / 2, Im: (u - v) / 2}
}
