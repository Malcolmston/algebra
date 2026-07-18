// Package special implements special mathematical functions in pure Go.
//
// The package collects the classical higher transcendental functions that do
// not appear in the standard math package, grouped by family:
//
//   - Bessel functions of the first and second kind (J0, J1, Jn, Y0, Y1, Yn),
//     the modified Bessel functions (I0, I1, In, K0, K1, Kn), the spherical
//     Bessel functions (SphericalJn, SphericalYn) and the Struve functions
//     (StruveH0, StruveH1).
//   - Airy functions Ai, Bi and their derivatives.
//   - Elliptic integrals and related quantities.
//   - The error function family, including the Fresnel integrals.
//   - The exponential, sine and cosine integrals.
//   - The Riemann zeta family (zeta, eta, beta).
//   - The Lambert W function.
//   - The gamma family (log-gamma, digamma, incomplete forms).
//
// Every routine is implemented with the Go standard library only, using series
// expansions for small arguments and asymptotic expansions for large ones so
// that results remain accurate across the full real line. The functions are
// deterministic and aim for correctness first: they are validated in the test
// suite against known closed-form and reference values.
package special

import "math"

// -----------------------------------------------------------------------------
// Unexported helpers (all prefixed with "special" to avoid collisions with
// sibling files in this package).
// -----------------------------------------------------------------------------

// specialEulerGamma is the Euler–Mascheroni constant γ.
const specialEulerGamma = 0.57721566490153286060651209008240243104215933593992

// -----------------------------------------------------------------------------
// Bessel functions of the first kind.
// -----------------------------------------------------------------------------

// BesselJ0 returns the Bessel function of the first kind of order zero, J0(x).
func BesselJ0(x float64) float64 {
	return math.J0(x)
}

// BesselJ1 returns the Bessel function of the first kind of order one, J1(x).
func BesselJ1(x float64) float64 {
	return math.J1(x)
}

// BesselJn returns the Bessel function of the first kind of integer order n,
// Jn(x). It is valid for any integer n (negative orders use the reflection
// Jn(-n, x) = (-1)^n Jn(n, x)).
func BesselJn(n int, x float64) float64 {
	return math.Jn(n, x)
}

// -----------------------------------------------------------------------------
// Bessel functions of the second kind.
// -----------------------------------------------------------------------------

// BesselY0 returns the Bessel function of the second kind of order zero, Y0(x).
// It is defined for x > 0.
func BesselY0(x float64) float64 {
	return math.Y0(x)
}

// BesselY1 returns the Bessel function of the second kind of order one, Y1(x).
// It is defined for x > 0.
func BesselY1(x float64) float64 {
	return math.Y1(x)
}

// BesselYn returns the Bessel function of the second kind of integer order n,
// Yn(x). It is defined for x > 0.
func BesselYn(n int, x float64) float64 {
	return math.Yn(n, x)
}

// -----------------------------------------------------------------------------
// Modified Bessel functions of the first kind.
// -----------------------------------------------------------------------------

// BesselI0 returns the modified Bessel function of the first kind of order
// zero, I0(x). The result is even in x and always ≥ 1.
func BesselI0(x float64) float64 {
	ax := math.Abs(x)
	if ax < 3.75 {
		t := x / 3.75
		t *= t
		return 1.0 + t*(3.5156229+t*(3.0899424+t*(1.2067492+
			t*(0.2659732+t*(0.0360768+t*0.0045813)))))
	}
	t := 3.75 / ax
	poly := 0.39894228 + t*(0.01328592+t*(0.00225319+t*(-0.00157565+
		t*(0.00916281+t*(-0.02057706+t*(0.02635537+t*(-0.01647633+
			t*0.00392377)))))))
	return (math.Exp(ax) / math.Sqrt(ax)) * poly
}

// BesselI1 returns the modified Bessel function of the first kind of order one,
// I1(x). The result is odd in x.
func BesselI1(x float64) float64 {
	ax := math.Abs(x)
	var ans float64
	if ax < 3.75 {
		t := x / 3.75
		t *= t
		ans = ax * (0.5 + t*(0.87890594+t*(0.51498869+t*(0.15084934+
			t*(0.02658733+t*(0.00301532+t*0.00032411))))))
	} else {
		t := 3.75 / ax
		poly := 0.02282967 + t*(-0.02895312+t*(0.01787654-t*0.00420059))
		poly = 0.39894228 + t*(-0.03988024+t*(-0.00362018+t*(0.00163801+
			t*(-0.01031555+t*poly))))
		ans = (math.Exp(ax) / math.Sqrt(ax)) * poly
	}
	if x < 0 {
		return -ans
	}
	return ans
}

// BesselIn returns the modified Bessel function of the first kind of integer
// order n, In(x). Negative orders satisfy In(-n, x) = In(n, x).
func BesselIn(n int, x float64) float64 {
	if n < 0 {
		n = -n
	}
	switch n {
	case 0:
		return BesselI0(x)
	case 1:
		return BesselI1(x)
	}
	if x == 0 {
		return 0
	}
	ax := math.Abs(x)
	const acc = 40.0
	const bigno = 1.0e10
	const bigni = 1.0e-10
	tox := 2.0 / ax
	bip, ans := 0.0, 0.0
	bi := 1.0
	// Downward recurrence (Miller's algorithm) for numerical stability.
	m := 2 * (n + int(math.Sqrt(acc*float64(n))))
	for j := m; j > 0; j-- {
		bim := bip + float64(j)*tox*bi
		bip = bi
		bi = bim
		if math.Abs(bi) > bigno {
			ans *= bigni
			bi *= bigni
			bip *= bigni
		}
		if j == n {
			ans = bip
		}
	}
	ans *= BesselI0(ax) / bi
	if x < 0 && (n&1) == 1 {
		return -ans
	}
	return ans
}

// -----------------------------------------------------------------------------
// Modified Bessel functions of the second kind.
// -----------------------------------------------------------------------------

// BesselK0 returns the modified Bessel function of the second kind of order
// zero, K0(x). It is defined for x > 0.
func BesselK0(x float64) float64 {
	if x <= 0 {
		return math.NaN()
	}
	if x <= 2.0 {
		t := x * x / 4.0
		return (-math.Log(x/2.0) * BesselI0(x)) + (-specialEulerGamma +
			t*(0.42278420+t*(0.23069756+t*(0.03488590+t*(0.00262698+
				t*(0.00010750+t*0.00000740))))))
	}
	t := 2.0 / x
	poly := 1.25331414 + t*(-0.07832358+t*(0.02189568+t*(-0.01062446+
		t*(0.00587872+t*(-0.00251540+t*0.00053208)))))
	return (math.Exp(-x) / math.Sqrt(x)) * poly
}

// BesselK1 returns the modified Bessel function of the second kind of order
// one, K1(x). It is defined for x > 0.
func BesselK1(x float64) float64 {
	if x <= 0 {
		return math.NaN()
	}
	if x <= 2.0 {
		t := x * x / 4.0
		return (math.Log(x/2.0) * BesselI1(x)) + (1.0/x)*(1.0+
			t*(0.15443144+t*(-0.67278579+t*(-0.18156897+t*(-0.01919402+
				t*(-0.00110404-t*0.00004686))))))
	}
	t := 2.0 / x
	poly := 1.25331414 + t*(0.23498619+t*(-0.03655620+t*(0.01504268+
		t*(-0.00780353+t*(0.00325614-t*0.00068245)))))
	return (math.Exp(-x) / math.Sqrt(x)) * poly
}

// BesselKn returns the modified Bessel function of the second kind of integer
// order n, Kn(x). It is defined for x > 0. Negative orders satisfy
// Kn(-n, x) = Kn(n, x).
func BesselKn(n int, x float64) float64 {
	if n < 0 {
		n = -n
	}
	if x <= 0 {
		return math.NaN()
	}
	switch n {
	case 0:
		return BesselK0(x)
	case 1:
		return BesselK1(x)
	}
	// Upward recurrence: K_{m+1} = K_{m-1} + (2m/x) K_m.
	tox := 2.0 / x
	bkm := BesselK0(x)
	bk := BesselK1(x)
	for j := 1; j < n; j++ {
		bkp := bkm + float64(j)*tox*bk
		bkm = bk
		bk = bkp
	}
	return bk
}

// -----------------------------------------------------------------------------
// Spherical Bessel functions.
// -----------------------------------------------------------------------------

// SphericalJ0 returns the spherical Bessel function of the first kind of order
// zero, j0(x) = sin(x)/x, with the removable singularity handled at x = 0.
func SphericalJ0(x float64) float64 {
	if x == 0 {
		return 1
	}
	return math.Sin(x) / x
}

// SphericalJ1 returns the spherical Bessel function of the first kind of order
// one, j1(x) = sin(x)/x^2 - cos(x)/x.
func SphericalJ1(x float64) float64 {
	if x == 0 {
		return 0
	}
	return math.Sin(x)/(x*x) - math.Cos(x)/x
}

// SphericalJn returns the spherical Bessel function of the first kind of order
// n, jn(x) = sqrt(pi/(2x)) J_{n+1/2}(x). It uses stable upward recurrence for
// x ≥ n and downward (Miller) recurrence otherwise.
func SphericalJn(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if x == 0 {
		if n == 0 {
			return 1
		}
		return 0
	}
	if n == 0 {
		return SphericalJ0(x)
	}
	if n == 1 {
		return SphericalJ1(x)
	}
	if x >= float64(n) {
		// Upward recurrence is stable when x >= n.
		jm := SphericalJ0(x)
		j := SphericalJ1(x)
		for k := 1; k < n; k++ {
			jp := float64(2*k+1)/x*j - jm
			jm = j
			j = jp
		}
		return j
	}
	// Downward recurrence for x < n, then normalise against j0.
	const start = 20
	m := n + start + int(math.Sqrt(float64(start*n)))
	jp := 0.0
	j := 1.0e-30
	var target float64
	for k := m; k >= 0; k-- {
		jm := float64(2*k+3)/x*j - jp
		jp = j
		j = jm
		if k == n {
			// After the update, j holds the (unnormalised) j_n(x).
			target = j
		}
	}
	return target * (SphericalJ0(x) / j)
}

// SphericalY0 returns the spherical Bessel function of the second kind of order
// zero, y0(x) = -cos(x)/x. It is defined for x > 0.
func SphericalY0(x float64) float64 {
	return -math.Cos(x) / x
}

// SphericalY1 returns the spherical Bessel function of the second kind of order
// one, y1(x) = -cos(x)/x^2 - sin(x)/x. It is defined for x > 0.
func SphericalY1(x float64) float64 {
	return -math.Cos(x)/(x*x) - math.Sin(x)/x
}

// SphericalYn returns the spherical Bessel function of the second kind of order
// n, yn(x). It is defined for x > 0 and uses stable upward recurrence.
func SphericalYn(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return SphericalY0(x)
	}
	if n == 1 {
		return SphericalY1(x)
	}
	ym := SphericalY0(x)
	y := SphericalY1(x)
	for k := 1; k < n; k++ {
		yp := float64(2*k+1)/x*y - ym
		ym = y
		y = yp
	}
	return y
}

// -----------------------------------------------------------------------------
// Airy functions.
// -----------------------------------------------------------------------------

// specialAiryValues computes Ai(x), Ai'(x), Bi(x) and Bi'(x) together using
// series near the origin and asymptotic expansions for large |x|.
func specialAiryValues(x float64) (ai, aip, bi, bip float64) {
	// Ai(0), Ai'(0), Bi(0), Bi'(0).
	const (
		ai0  = 0.355028053887817239260063186004183157  // 1/(3^(2/3) Γ(2/3))
		aip0 = -0.258819403792806798405183560189525631 // -1/(3^(1/3) Γ(1/3))
		bi0  = 0.614926627446000735150922369093613089  // 1/(3^(1/6) Γ(2/3))
		bip0 = 0.448288357353826357914823710398828391  // 3^(1/6)/Γ(1/3)
	)
	ax := math.Abs(x)
	if ax < 8.0 {
		// Power series: f = sum x^(3k)/((3k)!!-like), g = sum x^(3k+1)/...
		// Use f(x)=sum_{k} c_k x^{3k}, g(x)=sum_{k} d_k x^{3k+1}.
		var f, fp, g, gp float64
		// f series terms t_f satisfy t_{k+1} = t_k * x^3 / ((3k+2)(3k+3)).
		tf := 1.0
		tg := x
		f = tf
		g = tg
		fp = 0.0
		gp = 1.0
		x3 := x * x * x
		for k := 0; k < 40; k++ {
			// derivative accumulation before updating terms
			// f' contribution: d/dx tf = 3k * tf / x, handled via power tracking
			nf := tf * x3 / (float64(3*k+2) * float64(3*k+3))
			ng := tg * x3 / (float64(3*k+3) * float64(3*k+4))
			f += nf
			g += ng
			fp += float64(3*(k+1)) * nf / x
			gp += float64(3*(k+1)+1) * ng / x
			tf = nf
			tg = ng
			if math.Abs(nf)+math.Abs(ng) < 1e-18*(math.Abs(f)+math.Abs(g)) {
				break
			}
		}
		if x == 0 {
			fp = 0
			gp = 1
		}
		ai = ai0*f + aip0*g
		aip = ai0*fp + aip0*gp
		bi = bi0*f + bip0*g
		bip = bi0*fp + bip0*gp
		return
	}
	// Asymptotic expansions for large |x|.
	zeta := (2.0 / 3.0) * math.Pow(ax, 1.5)
	if x > 0 {
		sqx := math.Sqrt(ax)
		fac := 1.0 / (math.SqrtPi * math.Pow(x, 0.25))
		e := math.Exp(-zeta)
		ai = 0.5 * fac * e
		aip = -0.5 * fac * sqx * e
		bi = fac * math.Exp(zeta)
		bip = fac * sqx * math.Exp(zeta)
		return
	}
	// x < 0 oscillatory regime.
	fac := 1.0 / (math.SqrtPi * math.Pow(ax, 0.25))
	s, c := math.Sincos(zeta + math.Pi/4.0)
	sqx := math.Sqrt(ax)
	ai = fac * s
	bi = fac * c
	aip = fac * sqx * c
	bip = -fac * sqx * s
	return
}

// AiryAi returns the Airy function of the first kind, Ai(x).
func AiryAi(x float64) float64 {
	ai, _, _, _ := specialAiryValues(x)
	return ai
}

// AiryAiPrime returns the derivative of the Airy function of the first kind,
// Ai'(x).
func AiryAiPrime(x float64) float64 {
	_, aip, _, _ := specialAiryValues(x)
	return aip
}

// AiryBi returns the Airy function of the second kind, Bi(x).
func AiryBi(x float64) float64 {
	_, _, bi, _ := specialAiryValues(x)
	return bi
}

// AiryBiPrime returns the derivative of the Airy function of the second kind,
// Bi'(x).
func AiryBiPrime(x float64) float64 {
	_, _, _, bip := specialAiryValues(x)
	return bip
}

// -----------------------------------------------------------------------------
// Struve functions.
// -----------------------------------------------------------------------------

// StruveH0 returns the Struve function of order zero, H0(x). It is odd in x and
// uses a power series for small |x| and the relation to Y0 for large |x|.
func StruveH0(x float64) float64 {
	ax := math.Abs(x)
	var res float64
	if ax <= 20.0 {
		// H0(x) = (2/pi) sum_{k>=0} (-1)^k (x/2)^{2k+1} / [((2k+1)!!)^2 /? ]
		// Use term recurrence: t_0 = x, t_{k} = -t_{k-1} x^2 / ((2k+1)^2).
		term := ax
		sum := term
		x2 := ax * ax
		for k := 1; k < 100; k++ {
			term *= -x2 / (float64(2*k+1) * float64(2*k+1))
			sum += term
			if math.Abs(term) < 1e-17*math.Abs(sum) {
				break
			}
		}
		res = (2.0 / math.Pi) * sum
	} else {
		// Asymptotic: H0(x) ~ Y0(x) + (2/pi) sum ...
		res = BesselY0(ax) + specialStruveAsym(ax, 0)
	}
	if x < 0 {
		return -res
	}
	return res
}

// StruveH1 returns the Struve function of order one, H1(x). It is even in x and
// uses a power series for small |x| and the relation to Y1 for large |x|.
func StruveH1(x float64) float64 {
	ax := math.Abs(x)
	if ax <= 20.0 {
		// H1(x) = (2/pi) sum_{k>=0} (-1)^k (x/2)^{2k+2} / [ (2k+1)!! (2k+3)!! /?]
		// Term recurrence: t_0 = x^2/2? Use standard series
		// H1(x) = (2/pi) * sum_{k=0}^inf (-1)^k (x/2)^{2k+2} /
		//         (Gamma(k+3/2) Gamma(k+5/2)) * sqrt(pi)/2 ... simplify via ratio.
		half := ax / 2.0
		x2 := half * half
		term := 4.0 * x2 / 3.0 // k=0 term of the (2/pi)-scaled series
		sum := term
		for k := 1; k < 100; k++ {
			term *= -x2 / (float64(2*k+1) * float64(2*k+3) / 4.0)
			sum += term
			if math.Abs(term) < 1e-17*math.Abs(sum) {
				break
			}
		}
		return (2.0 / math.Pi) * sum
	}
	return BesselY1(ax) + specialStruveAsym(ax, 1)
}

// specialStruveAsym returns the (2/pi) * asymptotic tail added to Y_nu(x) in the
// large-argument expansion of the Struve function of order nu (0 or 1).
func specialStruveAsym(x float64, nu int) float64 {
	// H_nu(x) - Y_nu(x) ~ (1/pi) sum_{m=0}^{M} Gamma(m+1/2) (x/2)^{nu-2m-1} / Gamma(nu+1/2-m)
	// Evaluate a few terms.
	sum := 0.0
	for m := 0; m < 8; m++ {
		num := specialGammaHalfInt(2*m + 1) // Gamma(m+1/2) via helper: pass 2m+1
		den := specialGammaHalfInt(2*(nu-m) + 1)
		if math.IsInf(den, 0) || den == 0 {
			continue
		}
		term := num / den * math.Pow(x/2.0, float64(nu-2*m-1))
		sum += term
	}
	return sum / math.Pi
}

// specialGammaHalfInt returns Gamma(n/2) for a positive odd integer argument
// expressed as n (so Gamma((n)/2)); it supports the half-integer values that
// arise in the Struve asymptotic series and returns +Inf for non-positive
// half-integers.
func specialGammaHalfInt(n int) float64 {
	// n encodes the argument (n)/2 where n is odd.
	if n <= 0 {
		return math.Inf(1)
	}
	// Gamma((2m+1)/2) = (2m-1)!!/2^m * sqrt(pi); build iteratively.
	// Here argument = n/2 with n odd, so m = (n-1)/2.
	m := (n - 1) / 2
	res := math.SqrtPi
	for k := 1; k <= m; k++ {
		res *= (float64(2*k) - 1.0) / 2.0
	}
	return res
}

// -----------------------------------------------------------------------------
// Derivatives of the Bessel functions.
//
// These use the standard differentiation recurrences
//   J_n'(x) = (J_{n-1}(x) - J_{n+1}(x)) / 2
//   I_n'(x) = (I_{n-1}(x) + I_{n+1}(x)) / 2
//   K_n'(x) = -(K_{n-1}(x) + K_{n+1}(x)) / 2
// together with the special cases J_0' = -J_1, I_0' = I_1, K_0' = -K_1.
// -----------------------------------------------------------------------------

// BesselJ0Prime returns the derivative of J0, which equals -J1(x).
func BesselJ0Prime(x float64) float64 { return -math.J1(x) }

// BesselJ1Prime returns the derivative of J1(x).
func BesselJ1Prime(x float64) float64 { return 0.5 * (math.J0(x) - math.Jn(2, x)) }

// BesselJnPrime returns the derivative of the Bessel function Jn(x).
func BesselJnPrime(n int, x float64) float64 {
	if n == 0 {
		return -math.J1(x)
	}
	return 0.5 * (math.Jn(n-1, x) - math.Jn(n+1, x))
}

// BesselY0Prime returns the derivative of Y0, which equals -Y1(x). Defined for x > 0.
func BesselY0Prime(x float64) float64 { return -math.Y1(x) }

// BesselY1Prime returns the derivative of Y1(x). Defined for x > 0.
func BesselY1Prime(x float64) float64 { return 0.5 * (math.Y0(x) - math.Yn(2, x)) }

// BesselYnPrime returns the derivative of the Bessel function Yn(x). Defined for x > 0.
func BesselYnPrime(n int, x float64) float64 {
	if n == 0 {
		return -math.Y1(x)
	}
	return 0.5 * (math.Yn(n-1, x) - math.Yn(n+1, x))
}

// BesselI0Prime returns the derivative of I0, which equals I1(x).
func BesselI0Prime(x float64) float64 { return BesselI1(x) }

// BesselI1Prime returns the derivative of I1(x).
func BesselI1Prime(x float64) float64 { return 0.5 * (BesselI0(x) + BesselIn(2, x)) }

// BesselInPrime returns the derivative of the modified Bessel function In(x).
func BesselInPrime(n int, x float64) float64 {
	if n == 0 {
		return BesselI1(x)
	}
	return 0.5 * (BesselIn(n-1, x) + BesselIn(n+1, x))
}

// BesselK0Prime returns the derivative of K0, which equals -K1(x). Defined for x > 0.
func BesselK0Prime(x float64) float64 { return -BesselK1(x) }

// BesselK1Prime returns the derivative of K1(x). Defined for x > 0.
func BesselK1Prime(x float64) float64 { return -0.5 * (BesselK0(x) + BesselKn(2, x)) }

// BesselKnPrime returns the derivative of the modified Bessel function Kn(x). Defined for x > 0.
func BesselKnPrime(n int, x float64) float64 {
	if n == 0 {
		return -BesselK1(x)
	}
	return -0.5 * (BesselKn(n-1, x) + BesselKn(n+1, x))
}

// SphericalJnPrime returns the derivative of the spherical Bessel function
// jn(x), using jn'(x) = j_{n-1}(x) - (n+1)/x * jn(x).
func SphericalJnPrime(n int, x float64) float64 {
	if n == 0 {
		if x == 0 {
			return 0
		}
		return (x*math.Cos(x) - math.Sin(x)) / (x * x)
	}
	return SphericalJn(n-1, x) - float64(n+1)/x*SphericalJn(n, x)
}

// SphericalYnPrime returns the derivative of the spherical Bessel function
// yn(x), using yn'(x) = y_{n-1}(x) - (n+1)/x * yn(x). Defined for x > 0.
func SphericalYnPrime(n int, x float64) float64 {
	if n == 0 {
		return (x*math.Sin(x) + math.Cos(x)) / (x * x)
	}
	return SphericalYn(n-1, x) - float64(n+1)/x*SphericalYn(n, x)
}

// -----------------------------------------------------------------------------
// Exponentially scaled modified Bessel functions.
//
// These avoid overflow for large arguments by returning e^{-|x|} I_n(x) and
// e^{x} K_n(x) directly.
// -----------------------------------------------------------------------------

// BesselI0e returns the exponentially scaled I0, namely e^{-|x|} I0(x).
func BesselI0e(x float64) float64 {
	ax := math.Abs(x)
	if ax < 3.75 {
		t := x / 3.75
		t *= t
		return math.Exp(-ax) * (1.0 + t*(3.5156229+t*(3.0899424+t*(1.2067492+
			t*(0.2659732+t*(0.0360768+t*0.0045813))))))
	}
	t := 3.75 / ax
	poly := 0.39894228 + t*(0.01328592+t*(0.00225319+t*(-0.00157565+
		t*(0.00916281+t*(-0.02057706+t*(0.02635537+t*(-0.01647633+
			t*0.00392377)))))))
	return poly / math.Sqrt(ax)
}

// BesselI1e returns the exponentially scaled I1, namely e^{-|x|} I1(x).
func BesselI1e(x float64) float64 {
	ax := math.Abs(x)
	var ans float64
	if ax < 3.75 {
		t := x / 3.75
		t *= t
		ans = math.Exp(-ax) * ax * (0.5 + t*(0.87890594+t*(0.51498869+
			t*(0.15084934+t*(0.02658733+t*(0.00301532+t*0.00032411))))))
	} else {
		t := 3.75 / ax
		poly := 0.02282967 + t*(-0.02895312+t*(0.01787654-t*0.00420059))
		poly = 0.39894228 + t*(-0.03988024+t*(-0.00362018+t*(0.00163801+
			t*(-0.01031555+t*poly))))
		ans = poly / math.Sqrt(ax)
	}
	if x < 0 {
		return -ans
	}
	return ans
}

// BesselIne returns the exponentially scaled In, namely e^{-|x|} In(x).
func BesselIne(n int, x float64) float64 {
	return BesselIn(n, x) * math.Exp(-math.Abs(x))
}

// BesselK0e returns the exponentially scaled K0, namely e^{x} K0(x). Defined for x > 0.
func BesselK0e(x float64) float64 {
	if x <= 0 {
		return math.NaN()
	}
	if x <= 2.0 {
		return BesselK0(x) * math.Exp(x)
	}
	t := 2.0 / x
	poly := 1.25331414 + t*(-0.07832358+t*(0.02189568+t*(-0.01062446+
		t*(0.00587872+t*(-0.00251540+t*0.00053208)))))
	return poly / math.Sqrt(x)
}

// BesselK1e returns the exponentially scaled K1, namely e^{x} K1(x). Defined for x > 0.
func BesselK1e(x float64) float64 {
	if x <= 0 {
		return math.NaN()
	}
	if x <= 2.0 {
		return BesselK1(x) * math.Exp(x)
	}
	t := 2.0 / x
	poly := 1.25331414 + t*(0.23498619+t*(-0.03655620+t*(0.01504268+
		t*(-0.00780353+t*(0.00325614-t*0.00068245)))))
	return poly / math.Sqrt(x)
}

// BesselKne returns the exponentially scaled Kn, namely e^{x} Kn(x). Defined for x > 0.
func BesselKne(n int, x float64) float64 {
	return BesselKn(n, x) * math.Exp(x)
}

// -----------------------------------------------------------------------------
// Bessel functions of arbitrary real order (ascending series).
// -----------------------------------------------------------------------------

// BesselJnu returns the Bessel function of the first kind of arbitrary real
// order nu, J_nu(x), for x >= 0. It sums the convergent ascending power series
// and is accurate for small to moderate arguments.
func BesselJnu(nu, x float64) float64 {
	if x < 0 {
		return math.NaN()
	}
	if x == 0 {
		if nu == 0 {
			return 1
		}
		return 0
	}
	half := x / 2.0
	term := math.Pow(half, nu) / math.Gamma(nu+1)
	sum := term
	h2 := half * half
	for k := 1; k < 200; k++ {
		term *= -h2 / (float64(k) * (nu + float64(k)))
		sum += term
		if math.Abs(term) < 1e-17*math.Abs(sum) {
			break
		}
	}
	return sum
}

// BesselInu returns the modified Bessel function of the first kind of arbitrary
// real order nu, I_nu(x), for x >= 0. It sums the convergent ascending power
// series and is accurate for small to moderate arguments.
func BesselInu(nu, x float64) float64 {
	if x < 0 {
		return math.NaN()
	}
	if x == 0 {
		if nu == 0 {
			return 1
		}
		return 0
	}
	half := x / 2.0
	term := math.Pow(half, nu) / math.Gamma(nu+1)
	sum := term
	h2 := half * half
	for k := 1; k < 200; k++ {
		term *= h2 / (float64(k) * (nu + float64(k)))
		sum += term
		if math.Abs(term) < 1e-17*math.Abs(sum) {
			break
		}
	}
	return sum
}

// -----------------------------------------------------------------------------
// Kelvin functions.
// -----------------------------------------------------------------------------

// KelvinBer returns the Kelvin function ber(x) = Re[J0(x e^{3iπ/4})], summed
// from its convergent power series.
func KelvinBer(x float64) float64 {
	h := x / 2.0
	h4 := math.Pow(h, 4)
	term := 1.0
	sum := 1.0
	for k := 1; k < 100; k++ {
		d := float64(2*k-1) * float64(2*k) * float64(2*k-1) * float64(2*k)
		term *= -h4 / d
		sum += term
		if math.Abs(term) < 1e-18*math.Abs(sum) {
			break
		}
	}
	return sum
}

// KelvinBei returns the Kelvin function bei(x) = Im[J0(x e^{3iπ/4})], summed
// from its convergent power series.
func KelvinBei(x float64) float64 {
	h := x / 2.0
	h2 := h * h
	h4 := h2 * h2
	term := h2 // k=0 term
	sum := term
	for k := 1; k < 100; k++ {
		d := float64(2*k) * float64(2*k+1) * float64(2*k) * float64(2*k+1)
		term *= -h4 / d
		sum += term
		if math.Abs(term) < 1e-18*math.Abs(sum) {
			break
		}
	}
	return sum
}

// KelvinKer returns the Kelvin function ker(x) for x > 0, summed from its
// series with the leading logarithmic terms. It uses
//
//	ker(x) = -(ln(x/2)+γ) ber(x) + (π/4) bei(x)
//	         + Σ_{k≥0} (-1)^k φ(2k) / ((2k)!)^2 (x/2)^{4k}
//
// where φ(n) is the n-th harmonic number and φ(0) = 0.
func KelvinKer(x float64) float64 {
	if x <= 0 {
		return math.NaN()
	}
	h := x / 2.0
	h4 := math.Pow(h, 4)
	term := 1.0
	var sum, phi float64
	for k := 1; k < 100; k++ {
		d := float64(2*k-1) * float64(2*k) * float64(2*k-1) * float64(2*k)
		term *= -h4 / d
		phi += 1.0/float64(2*k-1) + 1.0/float64(2*k) // φ(2k)
		add := term * phi
		sum += add
		if math.Abs(add) < 1e-18*math.Abs(sum) {
			break
		}
	}
	return -(math.Log(h)+specialEulerGamma)*KelvinBer(x) + math.Pi/4.0*KelvinBei(x) + sum
}

// KelvinKei returns the Kelvin function kei(x) for x > 0, summed from its
// series with the leading logarithmic terms. It uses
//
//	kei(x) = -(ln(x/2)+γ) bei(x) - (π/4) ber(x)
//	         + Σ_{k≥0} (-1)^k φ(2k+1) / ((2k+1)!)^2 (x/2)^{4k+2}
//
// where φ(n) is the n-th harmonic number.
func KelvinKei(x float64) float64 {
	if x <= 0 {
		return math.NaN()
	}
	h := x / 2.0
	h2 := h * h
	h4 := h2 * h2
	term := h2  // k=0 term (x/2)^2 / (1!)^2
	phi := 1.0  // φ(1)
	sum := term // φ(1) = 1
	for k := 1; k < 100; k++ {
		d := float64(2*k) * float64(2*k+1) * float64(2*k) * float64(2*k+1)
		term *= -h4 / d
		phi += 1.0/float64(2*k) + 1.0/float64(2*k+1) // φ(2k+1)
		add := term * phi
		sum += add
		if math.Abs(add) < 1e-18*math.Abs(sum) {
			break
		}
	}
	return -(math.Log(h)+specialEulerGamma)*KelvinBei(x) - math.Pi/4.0*KelvinBer(x) + sum
}

// -----------------------------------------------------------------------------
// Modified Struve functions.
// -----------------------------------------------------------------------------

// StruveL0 returns the modified Struve function of order zero, L0(x). It is odd
// in x and uses the all-positive-term power series for moderate |x|.
func StruveL0(x float64) float64 {
	ax := math.Abs(x)
	term := ax
	sum := term
	x2 := ax * ax
	for k := 1; k < 200; k++ {
		term *= x2 / (float64(2*k+1) * float64(2*k+1))
		sum += term
		if math.Abs(term) < 1e-16*math.Abs(sum) {
			break
		}
	}
	res := (2.0 / math.Pi) * sum
	if x < 0 {
		return -res
	}
	return res
}

// StruveL1 returns the modified Struve function of order one, L1(x). It is even
// in x and uses the all-positive-term power series for moderate |x|.
func StruveL1(x float64) float64 {
	ax := math.Abs(x)
	half := ax / 2.0
	x2 := half * half
	term := 4.0 * x2 / 3.0
	sum := term
	for k := 1; k < 200; k++ {
		term *= x2 / (float64(2*k+1) * float64(2*k+3) / 4.0)
		sum += term
		if math.Abs(term) < 1e-16*math.Abs(sum) {
			break
		}
	}
	return (2.0 / math.Pi) * sum
}

// -----------------------------------------------------------------------------
// Modified spherical Bessel functions (Arfken normalisation).
// -----------------------------------------------------------------------------

// ModifiedSphericalI0 returns the modified spherical Bessel function of the
// first kind of order zero, i0(x) = sinh(x)/x.
func ModifiedSphericalI0(x float64) float64 {
	if x == 0 {
		return 1
	}
	return math.Sinh(x) / x
}

// ModifiedSphericalI1 returns the modified spherical Bessel function of the
// first kind of order one, i1(x) = (x cosh(x) - sinh(x)) / x^2.
func ModifiedSphericalI1(x float64) float64 {
	if x == 0 {
		return 0
	}
	return (x*math.Cosh(x) - math.Sinh(x)) / (x * x)
}

// ModifiedSphericalK0 returns the modified spherical Bessel function of the
// second kind of order zero, k0(x) = (pi/2) e^{-x}/x. Defined for x > 0.
func ModifiedSphericalK0(x float64) float64 {
	return math.Pi / 2.0 * math.Exp(-x) / x
}

// ModifiedSphericalK1 returns the modified spherical Bessel function of the
// second kind of order one, k1(x) = (pi/2) e^{-x} (1 + 1/x) / x. Defined for x > 0.
func ModifiedSphericalK1(x float64) float64 {
	return math.Pi / 2.0 * math.Exp(-x) * (1.0 + 1.0/x) / x
}

// -----------------------------------------------------------------------------
// Riccati-Bessel functions.
// -----------------------------------------------------------------------------

// RiccatiBesselPsi returns the Riccati-Bessel function psi_n(x) = x jn(x).
func RiccatiBesselPsi(n int, x float64) float64 {
	return x * SphericalJn(n, x)
}

// RiccatiBesselChi returns the Riccati-Bessel function chi_n(x) = -x yn(x).
func RiccatiBesselChi(n int, x float64) float64 {
	return -x * SphericalYn(n, x)
}

// -----------------------------------------------------------------------------
// Wronskians.
// -----------------------------------------------------------------------------

// BesselWronskianJY returns the Wronskian J_n(x) Y_n'(x) - J_n'(x) Y_n(x),
// which equals 2/(pi x) for every order n. Defined for x > 0.
func BesselWronskianJY(x float64) float64 {
	return 2.0 / (math.Pi * x)
}

// BesselWronskianIK returns the Wronskian I_n(x) K_n'(x) - I_n'(x) K_n(x),
// which equals -1/x for every order n. Defined for x > 0.
func BesselWronskianIK(x float64) float64 {
	return -1.0 / x
}

// -----------------------------------------------------------------------------
// Zeros.
// -----------------------------------------------------------------------------

// BesselJZero returns an approximation to the s-th positive zero (s >= 1) of the
// Bessel function J_n, using McMahon's asymptotic expansion. Accuracy improves
// rapidly with s.
func BesselJZero(n, s int) float64 {
	mu := 4.0 * float64(n) * float64(n)
	beta := (float64(s) + 0.5*float64(n) - 0.25) * math.Pi
	b8 := 8.0 * beta
	e := b8 * b8
	return beta - (mu-1)/b8 -
		4*(mu-1)*(7*mu-31)/(3*b8*e) -
		32*(mu-1)*(83*mu*mu-982*mu+3779)/(15*b8*e*e)
}

// BesselYZero returns an approximation to the s-th positive zero (s >= 1) of the
// Bessel function Y_n, using McMahon's asymptotic expansion. Accuracy improves
// rapidly with s.
func BesselYZero(n, s int) float64 {
	mu := 4.0 * float64(n) * float64(n)
	beta := (float64(s) + 0.5*float64(n) - 0.75) * math.Pi
	b8 := 8.0 * beta
	e := b8 * b8
	return beta - (mu-1)/b8 -
		4*(mu-1)*(7*mu-31)/(3*b8*e) -
		32*(mu-1)*(83*mu*mu-982*mu+3779)/(15*b8*e*e)
}

// AiryZeroAi returns an approximation to the s-th real zero (s >= 1) of the Airy
// function Ai, using the standard asymptotic expansion. The returned value is
// negative.
func AiryZeroAi(s int) float64 {
	t := 3.0 * math.Pi * (4.0*float64(s) - 1.0) / 8.0
	return -specialAiryZeroT(t)
}

// AiryZeroBi returns an approximation to the s-th real zero (s >= 1) of the Airy
// function Bi, using the standard asymptotic expansion. The returned value is
// negative.
func AiryZeroBi(s int) float64 {
	t := 3.0 * math.Pi * (4.0*float64(s) - 3.0) / 8.0
	return -specialAiryZeroT(t)
}

// specialAiryZeroT evaluates the asymptotic series T(t) used for the Airy zeros.
func specialAiryZeroT(t float64) float64 {
	t2 := t * t
	return math.Pow(t, 2.0/3.0) * (1.0 + 5.0/48.0/t2 - 5.0/36.0/(t2*t2) +
		77125.0/82944.0/(t2*t2*t2))
}
