# controltheory — provides classical linear control-systems primitives implemented with the Go standard

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/controltheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/controltheory)

Package controltheory provides classical linear control-systems primitives
implemented with the Go standard library only.

- Transfer functions G(s) = N(s)/D(s), where N and D are real polynomials
stored as `Poly` values. - State-space realizations (A, B, C, D) stored as
`StateSpace` values.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/controltheory"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func DampingRatio(pole complex128) float64` | DampingRatio returns the damping ratio of a complex pole, defined as -Re(p)/\|p\|. |
| `func IsHurwitzStable(p Poly) bool` | IsHurwitzStable reports whether the polynomial is Hurwitz stable, i.e. |
| `func LogSpace(startExp, endExp float64, n int) []float64` | LogSpace returns n logarithmically spaced values between 10^startExp and 10^endExp inclusive. |
| `func MagnitudeDB(value complex128) float64` | MagnitudeDB returns 20·log10(\|value\|) in decibels. |
| `func NaturalFrequency(pole complex128) float64` | NaturalFrequency returns the undamped natural frequency of a complex pole, which is its magnitude \|p\|. |
| `func NumRightHalfPlaneRoots(p Poly) int` | NumRightHalfPlaneRoots returns the number of roots of the polynomial with strictly positive real part, equal to the number of sign changes in the… |
| `func PhaseDeg(value complex128) float64` | PhaseDeg returns the phase angle of value in degrees. |

### Types

| Type | What it is |
| --- | --- |
| `BodePoint` | BodePoint is a single sample of a Bode plot at one angular frequency. |
| `NyquistPoint` | NyquistPoint is a single sample of a Nyquist plot at one angular frequency. |
| `PIDController` | PIDController is a discrete-time proportional-integral-derivative controller with configurable gains. |
| `Poly` | Poly represents a real polynomial in the Laplace variable s using the ascending-power convention: element i is the coefficient of s^i. |
| `RouthResult` | RouthResult holds the outcome of a Routh-Hurwitz stability analysis of a characteristic polynomial. |
| `StateSpace` | StateSpace is a SISO continuous-time state-space realization |
| `TransferFunction` | TransferFunction is a SISO continuous-time transfer function G(s) = Num(s) / Den(s) with real polynomial numerator and denominator stored in the… |

<details>
<summary><code>PIDController</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPIDController(kp, ki, kd float64) *PIDController` | NewPIDController returns a PIDController with the given proportional, integral, and derivative gains and zero internal state. |
| `func ZieglerNicholsPID(ku, pu float64) *PIDController` | ZieglerNicholsPID returns PID gains tuned by the classical Ziegler-Nichols ultimate-sensitivity rule from the ultimate gain ku (the proportional gain… |
| `func (c *PIDController) Integral() float64` | Integral returns the current accumulated integral of the error. |
| `func (c *PIDController) Reset()` | Reset clears the accumulated integral term and derivative history so the controller behaves as if freshly constructed. |
| `func (c *PIDController) TransferFunction() TransferFunction` | TransferFunction returns the ideal continuous-time transfer function of the controller, C(s) = Kp + Ki/s + Kd·s = (Kd·s^2 + Kp·s + Ki) / s. |
| `func (c *PIDController) Update(errValue, dt float64) float64` | Update advances the controller by one time step of length dt given the current error (setpoint minus measurement) and returns the control output |

</details>

<details>
<summary><code>Poly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPoly(coeffs ...float64) Poly` | NewPoly returns a Poly built from the given ascending-power coefficients. |
| `func PolyFromRoots(roots ...complex128) Poly` | PolyFromRoots returns the monic real polynomial whose roots are the given complex values. |
| `func (p Poly) Add(q Poly) Poly` | Add returns the sum p + q as a new polynomial. |
| `func (p Poly) Degree() int` | Degree returns the degree of the polynomial, i.e. |
| `func (p Poly) Derivative() Poly` | Derivative returns the derivative dp/ds as a new polynomial. |
| `func (p Poly) DivMod(q Poly) (quotient, remainder Poly)` | DivMod divides p by q and returns the quotient and remainder polynomials such that p = quotient*q + remainder with degree(remainder) < degree(q). |
| `func (p Poly) Eval(x float64) float64` | Eval evaluates the polynomial at the real point x using Horner's method. |
| `func (p Poly) EvalComplex(z complex128) complex128` | EvalComplex evaluates the polynomial at the complex point z using Horner's method, which is the form needed for frequency-response calculations. |
| `func (p Poly) LeadingCoeff() float64` | LeadingCoeff returns the coefficient of the highest nonzero power of s. |
| `func (p Poly) Mul(q Poly) Poly` | Mul returns the product p*q. |
| `func (p Poly) Roots() []complex128` | Roots returns all complex roots of the polynomial using the Durand-Kerner (Weierstrass) iteration with deterministic starting points. |
| `func (p Poly) Scale(k float64) Poly` | Scale returns the polynomial with every coefficient multiplied by k. |
| `func (p Poly) Sub(q Poly) Poly` | Sub returns the difference p - q as a new polynomial. |

</details>

<details>
<summary><code>RouthResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func RouthHurwitz(p Poly) RouthResult` | RouthHurwitz performs the Routh-Hurwitz stability test on the given characteristic polynomial (ascending-power `Poly`). |

</details>

<details>
<summary><code>StateSpace</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewStateSpace(a [][]float64, b, c []float64, d float64) StateSpace` | NewStateSpace constructs a StateSpace from the given matrices. |
| `func TransferFunctionToStateSpace(g TransferFunction) StateSpace` | TransferFunctionToStateSpace converts a proper transfer function into a state-space realization in controllable canonical form. |
| `func (s StateSpace) CharacteristicPolynomial() Poly` | CharacteristicPolynomial returns the characteristic polynomial det(sI - A) of the state matrix as an ascending-power monic `Poly`, computed with the… |
| `func (s StateSpace) ControllabilityMatrix() [][]float64` | ControllabilityMatrix returns the controllability matrix [B, AB, A^2 B, ..., A^(n-1) B], an n×n matrix for a SISO system. |
| `func (s StateSpace) ControllabilityRank() int` | ControllabilityRank returns the rank of the controllability matrix. |
| `func (s StateSpace) ImpulseResponse(times []float64) []float64` | ImpulseResponse returns the unit-impulse response C·e^{At}·B of the state-space system sampled at the given time points, obtained by simulating the… |
| `func (s StateSpace) IsControllable() bool` | IsControllable reports whether the system is completely state controllable, i.e. |
| `func (s StateSpace) IsObservable() bool` | IsObservable reports whether the system is completely observable, i.e. |
| `func (s StateSpace) ObservabilityMatrix() [][]float64` | ObservabilityMatrix returns the observability matrix [C; CA; CA^2; ...; CA^(n-1)], an n×n matrix for a SISO system. |
| `func (s StateSpace) ObservabilityRank() int` | ObservabilityRank returns the rank of the observability matrix. |
| `func (s StateSpace) Order() int` | Order returns the number of states n, the dimension of the A matrix. |
| `func (s StateSpace) Poles() []complex128` | Poles returns the poles of the realization, i.e. |
| `func (s StateSpace) Simulate(times []float64, x0 []float64, input func(t float64) float64) []float64` | Simulate integrates the state-space system forward from initial state x0 over the given monotonically increasing time points using the classical… |
| `func (s StateSpace) StepResponse(times []float64) []float64` | StepResponse returns the unit-step response of the state-space system from rest, sampled at the given time points. |
| `func (s StateSpace) TransferFunction() TransferFunction` | TransferFunction converts the state-space realization to an equivalent SISO transfer function G(s) = C(sI-A)^{-1}B + D using the Faddeev-LeVerrier… |

</details>

<details>
<summary><code>TransferFunction</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Feedback(g, h TransferFunction, sign int) TransferFunction` | Feedback returns the closed-loop transfer function of forward path g with feedback path h. |
| `func NewTransferFunction(num, den []float64) TransferFunction` | NewTransferFunction builds a TransferFunction from ascending-power numerator and denominator coefficient slices. |
| `func Parallel(g, h TransferFunction) TransferFunction` | Parallel returns the parallel connection G(s)+H(s) of two transfer functions, whose outputs are summed for a common input. |
| `func SecondOrderSystem(wn, zeta float64) TransferFunction` | SecondOrderSystem returns the canonical second-order transfer function wn^2 / (s^2 + 2·zeta·wn·s + wn^2) for natural frequency wn (rad/s) and… |
| `func Series(g, h TransferFunction) TransferFunction` | Series returns the cascade connection G(s)·H(s) of two transfer functions. |
| `func UnityFeedback(g TransferFunction) TransferFunction` | UnityFeedback returns the closed-loop transfer function of forward path g with unity negative feedback, G/(1+G). |
| `func (g TransferFunction) Bode(omegas []float64) []BodePoint` | Bode returns the Bode-plot samples (magnitude in dB and phase in degrees) of the transfer function at each supplied angular frequency. |
| `func (g TransferFunction) DCGain() float64` | DCGain returns the steady-state gain G(0) = Num(0)/Den(0). |
| `func (g TransferFunction) Evaluate(s complex128) complex128` | Evaluate returns the complex value of G(s) at the complex point s. |
| `func (g TransferFunction) FinalValue() float64` | FinalValue returns the steady-state value of the step response predicted by the final-value theorem, lim_{t->inf} y(t) = G(0), when the closed system… |
| `func (g TransferFunction) FrequencyResponse(omega float64) complex128` | FrequencyResponse returns G(jω), the complex response at angular frequency omega (in radians per second). |
| `func (g TransferFunction) GainCrossoverFrequency(omegas []float64) (float64, bool)` | GainCrossoverFrequency returns the lowest angular frequency in the supplied grid at which the open-loop magnitude \|G(jω)\| crosses unity (0 dB),… |
| `func (g TransferFunction) GainMargin(omegas []float64) (float64, bool)` | GainMargin returns the gain margin in decibels, defined as -20·log10\|G(jω)\| at the phase crossover frequency (where the phase is -180 degrees),… |
| `func (g TransferFunction) ImpulseResponse(times []float64) []float64` | ImpulseResponse returns the unit-impulse response of a strictly proper transfer function sampled at the given monotonically increasing time points. |
| `func (g TransferFunction) IsProper() bool` | IsProper reports whether the numerator degree does not exceed the denominator degree, the condition for a physically realizable system. |
| `func (g TransferFunction) IsStable() bool` | IsStable reports whether every pole has a strictly negative real part, the condition for asymptotic stability of a continuous-time system. |
| `func (g TransferFunction) IsStrictlyProper() bool` | IsStrictlyProper reports whether the numerator degree is strictly less than the denominator degree. |
| `func (g TransferFunction) Nyquist(omegas []float64) []NyquistPoint` | Nyquist returns the Nyquist-plot samples (real and imaginary parts of G(jω)) of the transfer function at each supplied angular frequency. |
| `func (g TransferFunction) Order() int` | Order returns the order of the system, i.e. |
| `func (g TransferFunction) PhaseCrossoverFrequency(omegas []float64) (float64, bool)` | PhaseCrossoverFrequency returns the lowest angular frequency in the supplied grid at which the open-loop phase crosses -180 degrees, refined by… |
| `func (g TransferFunction) PhaseMargin(omegas []float64) (float64, bool)` | PhaseMargin returns the phase margin in degrees, defined as 180 plus the open-loop phase at the gain crossover frequency, evaluated over the supplied… |
| `func (g TransferFunction) Poles() []complex128` | Poles returns the poles of the system, i.e. |
| `func (g TransferFunction) Series(h TransferFunction) TransferFunction` | Series returns the transfer function of g followed by h connected in series (cascade): the product G(s)·H(s). |
| `func (g TransferFunction) StepResponse(times []float64) []float64` | StepResponse returns the unit-step response of the transfer function sampled at the given monotonically increasing time points. |
| `func (g TransferFunction) Zeros() []complex128` | Zeros returns the zeros of the system, i.e. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/controltheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
