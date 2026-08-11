# physics — "science math" layer of the algebra module

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/physics.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/physics)

Package physics is the "science math" layer of the algebra module: a
standard-library-only collection of physical constants, unit conversions and
common physics formulas expressed in SI units.

- Physical constants. Exported float64 values in SI units, each documented
with its name, symbol, value and unit. Constants that are exact by definition
of the SI (2019 redefinition) are marked as such. The full set is also
available programmatically through `Constants` and `Lookup`, each entry being
a `Constant` record — this powers generated documentation.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/physics"
```

## Usage

This is the package's own `ExampleLookup`, so it compiles and its output is
asserted on every `go test ./physics/`.

```go
c, _ := physics.Lookup("c")
	fmt.Printf("%s = %g %s\n", c.Name, c.Value, c.Unit)
```

```
Speed of light = 2.99792458e+08 m/s
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AddInto(dst *Vec3, a, b Vec3)` | AddInto stores the component-wise sum a + b into *dst. |
| `func Beta(v float64) float64` | Beta returns the dimensionless velocity ratio β = v/c for a speed v (m/s), where c is the speed of light. |
| `func COPHeatPump(tCold, tHot float64) float64` | COPHeatPump returns the ideal (Carnot) coefficient of performance of a heat pump, COP = Th/(Th − Tc), delivering heat to a hot reservoir at tHot… |
| `func COPRefrigerator(tCold, tHot float64) float64` | COPRefrigerator returns the ideal (Carnot) coefficient of performance of a refrigerator, COP = Tc/(Th − Tc), moving heat from a cold reservoir at… |
| `func CapacitorCharge(c, v float64) float64` | CapacitorCharge returns the charge Q = C·V stored on a capacitor of capacitance c (F) at voltage v (V). |
| `func CapacitorEnergy(c, v float64) float64` | CapacitorEnergy returns the energy E = ½·C·V² stored in a capacitor of capacitance c (F) charged to voltage v (V). |
| `func CapacitorsParallel(c ...float64) float64` | CapacitorsParallel returns the equivalent capacitance of capacitors connected in parallel, the sum of the individual capacitances c (F). |
| `func CapacitorsSeries(c ...float64) float64` | CapacitorsSeries returns the equivalent capacitance of capacitors connected in series, 1/Σ(1/cᵢ), for capacitances c (F). |
| `func CarnotEfficiency(tCold, tHot float64) float64` | CarnotEfficiency returns the maximum thermal efficiency 1 − Tc/Th of a heat engine operating between a cold reservoir at tCold (K) and a hot… |
| `func CelsiusToFahrenheit(c float64) float64` | CelsiusToFahrenheit converts a temperature in degrees Celsius to degrees Fahrenheit. |
| `func CelsiusToKelvin(c float64) float64` | CelsiusToKelvin converts a temperature in degrees Celsius to kelvin. |
| `func CenterOfMass1D(masses, positions []float64) (float64, error)` | CenterOfMass1D returns the mass-weighted mean position (Σ mᵢ·xᵢ)/(Σ mᵢ) of a system of point masses arranged on a line, where masses`i` (kg)… |
| `func CentripetalAcceleration(v, r float64) float64` | CentripetalAcceleration returns the centripetal acceleration a = v²/r required to keep a body moving at speed v (m/s) on a circle of radius r (m). |
| `func Convert(value float64, fromUnit, toUnit string) (float64, error)` | Convert converts value from the unit fromUnit to the unit toUnit and returns the result. |
| `func CoulombForce(q1, q2, r float64) float64` | CoulombForce returns the magnitude and sign of the electrostatic force between two point charges q1 and q2 (C) separated by distance r (m): F =… |
| `func CurrentOhm(v, r float64) float64` | CurrentOhm returns the current I = V/R through a resistor of resistance r (Ω) with voltage v (V) across it, per Ohm's law. |
| `func DeBroglieWavelength(m, v float64) float64` | DeBroglieWavelength returns the de Broglie wavelength h/(m·v) of a particle of mass m (kg) moving at speed v (m/s), where h is the Planck constant. |
| `func DegToRad(deg float64) float64` | DegToRad converts an angle in degrees to radians. |
| `func Displacement(v0, a, t float64) float64` | Displacement returns the displacement s = v0·t + ½·a·t² of a body with initial velocity v0 (m/s) under constant acceleration a (m/s²) after… |
| `func Divergence(f func(Vec3) Vec3, p Vec3, h float64) float64` | Divergence returns a numerical approximation of the divergence ∇·f of the vector field f at the point p, using second-order central finite… |
| `func EVToJoules(ev float64) float64` | EVToJoules converts an energy in electronvolts to joules. |
| `func ElasticCollision1D(m1, v1, m2, v2 float64) (v1f, v2f float64)` | ElasticCollision1D returns the final velocities v1f and v2f (m/s) of two bodies after a one-dimensional perfectly elastic head-on collision, given… |
| `func ElectricField(q, r float64) float64` | ElectricField returns the magnitude of the electric field of a point charge, E = q/(4·π·ε0·r²), for a charge q (C) at distance r (m). |
| `func ElectricFieldPointCharge(q, r float64) float64` | ElectricFieldPointCharge returns the magnitude of the electric field E = k_e·q/r² produced by a point charge q (C) at distance r (m), where k_e is… |
| `func ElectricPotential(q, r float64) float64` | ElectricPotential returns the electric potential of a point charge, V = q/(4·π·ε0·r), for a charge q (C) at distance r (m), taking the potential… |
| `func ElectricPotentialPointCharge(q, r float64) float64` | ElectricPotentialPointCharge returns the electric potential V = k_e·q/r of a point charge q (C) at distance r (m), taking the potential to be zero… |
| `func EnergyMomentumRelation(m, p float64) float64` | EnergyMomentumRelation returns the total energy E of a body from the relativistic energy–momentum relation E = √((p·c)² + (m·c²)²), for rest… |
| `func EscapeVelocity(m, r float64) float64` | EscapeVelocity returns the escape speed from the surface of a body of mass m (kg) and radius r (m): v = √(2·G·m/r), where G is the Newtonian… |
| `func FahrenheitToCelsius(f float64) float64` | FahrenheitToCelsius converts a temperature in degrees Fahrenheit to degrees Celsius. |
| `func FahrenheitToKelvin(f float64) float64` | FahrenheitToKelvin converts a temperature in degrees Fahrenheit to kelvin. |
| `func FieldEnergyDensity(E float64) float64` | FieldEnergyDensity returns the energy density of an electric field in vacuum, u = ½·ε0·E², for a field magnitude E (V/m). |
| `func Force(m, a float64) float64` | Force returns the net force m·a on a body of mass m (kg) undergoing acceleration a (m/s²), per Newton's second law. |
| `func FreeFallDistance(t float64) float64` | FreeFallDistance returns the distance ½·g·t² fallen from rest under standard gravity after time t (s), ignoring air resistance. |
| `func GravitationalForce(m1, m2, r float64) float64` | GravitationalForce returns the magnitude of the gravitational attraction F = G·m1·m2/r² between two point masses m1 and m2 (kg) separated by… |
| `func GravitationalPotentialEnergy(m1, m2, r float64) float64` | GravitationalPotentialEnergy returns the gravitational potential energy U = −G·m1·m2/r of two point masses m1 and m2 (kg) separated by distance r… |
| `func HeatCapacity(mass, specificHeat, dT float64) float64` | HeatCapacity returns the heat Q = m·c·ΔT required to change the temperature of a mass (kg) with specific heat capacity specificHeat… |
| `func IdealGasMoles(P, V, T float64) float64` | IdealGasMoles returns the amount of substance of an ideal gas from n = P·V/(R·T), for pressure P (Pa) in volume V (m³) at temperature T (K), where… |
| `func IdealGasPressure(n, T, V float64) float64` | IdealGasPressure returns the pressure of an ideal gas from the ideal gas law P = n·R·T/V, for n moles at temperature T (K) in volume V (m³), where… |
| `func IdealGasTemperature(P, V, n float64) float64` | IdealGasTemperature returns the temperature of an ideal gas from T = P·V/(n·R), for pressure P (Pa) in volume V (m³) with n moles, where R is the… |
| `func IdealGasVolume(n, T, P float64) float64` | IdealGasVolume returns the volume of an ideal gas from V = n·R·T/P, for n moles at temperature T (K) and pressure P (Pa), where R is the molar gas… |
| `func ImpulseMomentum(force, dt float64) float64` | ImpulseMomentum returns the impulse J = force·dt delivered by a constant force (N) acting over a time interval dt (s). |
| `func InductorEnergy(l, i float64) float64` | InductorEnergy returns the energy stored in an inductor, E = ½·L·I², for an inductance l (H) carrying current i (A). |
| `func InelasticCollision1D(m1, v1, m2, v2 float64) float64` | InelasticCollision1D returns the common final velocity of two bodies that stick together after a one-dimensional perfectly inelastic collision, given… |
| `func JoulesToEV(j float64) float64` | JoulesToEV converts an energy in joules to electronvolts. |
| `func KelvinToCelsius(k float64) float64` | KelvinToCelsius converts a temperature in kelvin to degrees Celsius. |
| `func KelvinToFahrenheit(k float64) float64` | KelvinToFahrenheit converts a temperature in kelvin to degrees Fahrenheit. |
| `func KilogramsToPounds(kg float64) float64` | KilogramsToPounds converts a mass in kilograms to pounds (avoirdupois). |
| `func KineticEnergy(m, v float64) float64` | KineticEnergy returns the translational kinetic energy ½·m·v² of a body of mass m (kg) moving at speed v (m/s). |
| `func Laplacian(f func(Vec3) float64, p Vec3, h float64) float64` | Laplacian returns a numerical approximation of the Laplacian ∇²f of the scalar field f at the point p, using the second-order central… |
| `func LatentHeat(mass, L float64) float64` | LatentHeat returns the heat Q = m·L absorbed or released during a phase change of a mass (kg) with specific latent heat L (J/kg). |
| `func LengthContraction(properLength, v float64) float64` | LengthContraction returns the contracted length L = properLength/γ of an object measured in a frame in which it moves at speed v (m/s), where… |
| `func LorentzFactor(v float64) float64` | LorentzFactor returns the Lorentz factor γ = 1/√(1 − v²/c²) for a speed v (m/s), where c is the speed of light. |
| `func LorentzForceMag(q, E, v, B, theta float64) float64` | LorentzForceMag returns the combined electric and magnetic force on a charge as the scalar sum F = q·(E + v·B·sin(θ)), for a charge q (C) in an… |
| `func MagneticFieldStraightWire(current, r float64) float64` | MagneticFieldStraightWire returns the magnitude of the magnetic field around a long straight wire, B = μ0·I/(2·π·r), for a current current (A)… |
| `func MagneticForceOnCharge(q, v, B, theta float64) float64` | MagneticForceOnCharge returns the magnitude of the magnetic force on a moving charge, F = q·v·B·sin(θ), for a charge q (C) moving at speed v… |
| `func MassEnergy(m float64) float64` | MassEnergy returns the rest energy m·c² of a mass m (kg), where c is the speed of light. |
| `func MeanKineticEnergy(T float64) float64` | MeanKineticEnergy returns the mean translational kinetic energy of a gas particle, E = 3/2·k_B·T, at temperature T (K), where k_B is the Boltzmann… |
| `func MetersToMiles(m float64) float64` | MetersToMiles converts a length in metres to miles. |
| `func MilesToMeters(mi float64) float64` | MilesToMeters converts a length in miles to metres. |
| `func Momentum(m, v float64) float64` | Momentum returns the linear momentum m·v of a body of mass m (kg) moving at velocity v (m/s). |
| `func MostProbableSpeed(T, molarMass float64) float64` | MostProbableSpeed returns the most probable molecular speed of the Maxwell–Boltzmann distribution, v_p = √(2·R·T/M), for temperature T (K) and… |
| `func NumberDensity(P, T float64) float64` | NumberDensity returns the number of gas particles per unit volume from n = P/(k_B·T), for pressure P (Pa) at temperature T (K), where k_B is the… |
| `func OrbitalPeriod(m, r float64) float64` | OrbitalPeriod returns the period of a circular orbit of radius r (m) about a central mass m (kg), from Kepler's third law: T = 2π·√(r³/(G·m)),… |
| `func OrbitalVelocity(m, r float64) float64` | OrbitalVelocity returns the speed of a body in a circular orbit of radius r (m) about a central mass m (kg): v = √(G·m/r), where G is the… |
| `func ParallelCapacitance(c ...float64) float64` | ParallelCapacitance returns the equivalent capacitance of capacitors connected in parallel, the sum of the individual capacitances c (F), accumulated… |
| `func ParallelPlateCapacitance(area, separation, relPermittivity float64) float64` | ParallelPlateCapacitance returns the capacitance of a parallel-plate capacitor, C = ε0·εr·A/d, for plate area (m²), plate separation (m), and… |
| `func ParallelResistance(r ...float64) (float64, error)` | ParallelResistance returns the equivalent resistance of resistors connected in parallel, 1/Σ(1/rᵢ), for resistances r (Ω), accumulated in a… |
| `func PhotonEnergy(freq float64) float64` | PhotonEnergy returns the energy h·f of a photon of frequency freq (Hz), where h is the Planck constant. |
| `func PlanckSpectralRadiance(wavelength, T float64) float64` | PlanckSpectralRadiance returns the spectral radiance of a blackbody per unit wavelength from Planck's law, |
| `func PotentialEnergyGravity(m, h float64) float64` | PotentialEnergyGravity returns the gravitational potential energy m·g·h of a body of mass m (kg) at height h (m) near Earth's surface, using the… |
| `func PoundsToKilograms(lb float64) float64` | PoundsToKilograms converts a mass in pounds (avoirdupois) to kilograms. |
| `func Power(work, t float64) float64` | Power returns the average power work/time for work (J) done in time t (s). |
| `func PowerDissipated(i, r float64) float64` | PowerDissipated returns the electrical power dissipated by a resistor, P = I²·R, carrying current i (A) with resistance r (Ω). |
| `func ProjectileMaxHeight(v, theta float64) float64` | ProjectileMaxHeight returns the maximum height reached by a projectile launched from ground level with initial speed v (m/s) at launch angle theta… |
| `func ProjectilePosition(v, theta, t float64) (x, y float64)` | ProjectilePosition returns the horizontal and vertical position (x, y), in metres, of a projectile at time t (s), launched from the origin with… |
| `func ProjectileRange(v, theta float64) float64` | ProjectileRange returns the horizontal range of a projectile launched from and landing at the same height, with initial speed v (m/s) at launch angle… |
| `func ProjectileTimeOfFlight(v, theta float64) float64` | ProjectileTimeOfFlight returns the total time a projectile launched from and landing at the same height stays aloft, with initial speed v (m/s) at… |
| `func RCChargeVoltage(vSource, r, c, t float64) float64` | RCChargeVoltage returns the capacitor voltage while charging through a resistor, V = Vs·(1 − e^(−t/RC)), for a source voltage vSource (V)… |
| `func RCDischargeVoltage(v0, r, c, t float64) float64` | RCDischargeVoltage returns the capacitor voltage while discharging through a resistor, V = V0·e^(−t/RC), for an initial voltage v0 (V) across a… |
| `func RCTimeConstant(r, c float64) float64` | RCTimeConstant returns the time constant τ = R·C of an RC circuit, for a resistance r (Ω) and capacitance c (F). |
| `func RMSSpeed(T, molarMass float64) float64` | RMSSpeed returns the root-mean-square molecular speed of an ideal gas, v_rms = √(3·R·T/M), for temperature T (K) and molar mass molarMass… |
| `func RadToDeg(rad float64) float64` | RadToDeg converts an angle in radians to degrees. |
| `func RelativisticDopplerFactor(v float64) float64` | RelativisticDopplerFactor returns the relativistic Doppler frequency ratio f_observed/f_source for motion directly along the line of sight at speed v… |
| `func RelativisticEnergy(m, v float64) float64` | RelativisticEnergy returns the total relativistic energy E = γ·m·c² of a body of rest mass m (kg) moving at speed v (m/s), where γ is the… |
| `func RelativisticKineticEnergy(m, v float64) float64` | RelativisticKineticEnergy returns the relativistic kinetic energy K = (γ − 1)·m·c² of a body of rest mass m (kg) moving at speed v (m/s), where… |
| `func RelativisticMomentum(m, v float64) float64` | RelativisticMomentum returns the relativistic momentum p = γ·m·v of a body of rest mass m (kg) moving at speed v (m/s), where γ is the Lorentz… |
| `func RelativisticVelocityAddition(u, v float64) float64` | RelativisticVelocityAddition returns the combined velocity of two collinear velocities u and v (m/s) under special relativity: w = (u + v)/(1 +… |
| `func ResistanceOhm(v, i float64) float64` | ResistanceOhm returns the resistance R = V/I of a component with voltage v (V) across it carrying current i (A), per Ohm's law. |
| `func ResistorsParallel(r ...float64) float64` | ResistorsParallel returns the equivalent resistance of resistors connected in parallel, 1/Σ(1/rᵢ), for resistances r (Ω). |
| `func ResistorsSeries(r ...float64) float64` | ResistorsSeries returns the equivalent resistance of resistors connected in series, the sum of the individual resistances r (Ω). |
| `func ResonantFrequency(l, c float64) float64` | ResonantFrequency returns the resonant frequency of an LC circuit, f = 1/(2·π·√(L·C)), for an inductance l (H) and capacitance c (F). |
| `func ScaleInto(dst *Vec3, a Vec3, s float64)` | ScaleInto stores a scaled by the scalar s into *dst. |
| `func SeriesCapacitance(c ...float64) (float64, error)` | SeriesCapacitance returns the equivalent capacitance of capacitors connected in series, 1/Σ(1/cᵢ), for capacitances c (F), accumulated in a single… |
| `func SeriesResistance(r ...float64) float64` | SeriesResistance returns the equivalent resistance of resistors connected in series, the sum of the individual resistances r (Ω), accumulated in a… |
| `func StefanBoltzmannPower(emissivity, area, T float64) float64` | StefanBoltzmannPower returns the thermal power radiated by a grey body from the Stefan–Boltzmann law P = ε·σ·A·T⁴, for emissivity… |
| `func ThermalConduction(k, area, dT, thickness float64) float64` | ThermalConduction returns the steady-state rate of heat conduction through a slab, Q̇ = k·A·ΔT/L (Fourier's law), for thermal conductivity k… |
| `func TimeDilation(properTime, v float64) float64` | TimeDilation returns the dilated time interval Δt = γ·properTime measured in a frame moving at speed v (m/s) relative to the clock, where… |
| `func TimeToStop(v0, a float64) float64` | TimeToStop returns the time t = −v0/a for a body with initial velocity v0 (m/s) under constant acceleration a (m/s²) to reach zero velocity. |
| `func Velocity(v0, a, t float64) float64` | Velocity returns the velocity v = v0 + a·t of a body with initial velocity v0 (m/s) under constant acceleration a (m/s²) after time t (s). |
| `func VelocityFinal(v0, a, t float64) float64` | VelocityFinal returns the final velocity v = v0 + a·t of a body with initial velocity v0 (m/s) under constant acceleration a (m/s²) after time t… |
| `func VelocityFromDistance(v0, a, d float64) float64` | VelocityFromDistance returns the speed v = √(v0² + 2·a·d) reached by a body with initial velocity v0 (m/s) under constant acceleration a (m/s²)… |
| `func VisViva(mCentral, r, a float64) float64` | VisViva returns the orbital speed v = √(G·M·(2/r − 1/a)) of a body at distance r (m) from a central mass mCentral (kg) on an orbit of… |
| `func VoltageOhm(i, r float64) float64` | VoltageOhm returns the voltage V = I·R across a resistor carrying current i (A) with resistance r (Ω), per Ohm's law. |
| `func WaveSpeed(freq, wavelength float64) float64` | WaveSpeed returns the propagation speed of a wave, freq·wavelength, for a wave of frequency freq (Hz) and wavelength wavelength (m). |
| `func WienPeakWavelength(T float64) float64` | WienPeakWavelength returns the wavelength λ_max = b/T at which a blackbody's spectral radiance peaks (Wien's displacement law), for absolute… |
| `func Work(f, d float64) float64` | Work returns the work done by a constant force f (N) acting over a displacement d (m) parallel to the force: W = f·d. |

### Types

| Type | What it is |
| --- | --- |
| `Constant` | Constant is a self-describing physical constant record. |
| `Dim` | Dim is the exponent vector of the seven SI base dimensions. |
| `Measurement` | Measurement is a physical quantity paired with its standard (1σ) uncertainty. |
| `Quantity` | Quantity is a physical quantity: a numeric value together with its dimension. |
| `Vec3` | Vec3 is a three-dimensional Cartesian vector with float64 components. |

<details>
<summary><code>Constant</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Constants() []Constant` | Constants returns a copy of the table of physical constants known to the package, in a stable order. |
| `func ExtendedConstants() []Constant` | ExtendedConstants returns a copy of the table of additional physical constants defined in this file, in a stable order. |
| `func Lookup(symbol string) (Constant, bool)` | Lookup returns the `Constant` whose Symbol matches symbol exactly, and true if found. |
| `func LookupAny(symbol string) (Constant, bool)` | LookupAny returns the `Constant` whose Symbol matches symbol exactly, and true if found, searching both the base and the extended tables. |
| `func LookupByName(name string) (Constant, bool)` | LookupByName returns the `Constant` whose Name matches name exactly, and true if found, searching both the base and the extended tables. |

</details>

<details>
<summary><code>Dim</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (d Dim) String() string` | String renders the dimension as a middle-dot separated product of base-unit symbols with Unicode superscript exponents, for example "kg·m·s⁻²"… |

</details>

<details>
<summary><code>Measurement</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMeasurement(value, sigma float64) Measurement` | NewMeasurement returns a Measurement with the given value and 1σ standard uncertainty. |
| `func Propagate(f func([]float64) float64, means, sigmas []float64) (Measurement, error)` | Propagate propagates uncertainty through an arbitrary multi-variable function f evaluated at means, with per-input standard uncertainties sigmas. |
| `func (m Measurement) Add(o Measurement) Measurement` | Add returns the sum m + o. |
| `func (m Measurement) Apply(f func(float64) float64) Measurement` | Apply propagates uncertainty through an arbitrary single-variable function f. |
| `func (m Measurement) Div(o Measurement) Measurement` | Div returns the quotient m/o. |
| `func (m Measurement) Mul(o Measurement) Measurement` | Mul returns the product m·o. |
| `func (m Measurement) Pow(n float64) Measurement` | Pow returns m raised to the exact power n. |
| `func (m Measurement) Relative() float64` | Relative returns the fractional (relative) uncertainty Sigma/\|Value\|. |
| `func (m Measurement) Scale(k float64) Measurement` | Scale returns m multiplied by the exact constant k (a value carrying no uncertainty of its own). |
| `func (m Measurement) String() string` | String formats the measurement as "value ± sigma" using the shortest decimal representation that round-trips each component. |
| `func (m Measurement) Sub(o Measurement) Measurement` | Sub returns the difference m - o. |

</details>

<details>
<summary><code>Quantity</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewQuantity(value float64, dim Dim) Quantity` | NewQuantity returns a Quantity with the given value (in coherent SI base units) and dimension. |
| `func ParseQuantity(value float64, unitSymbol string) (Quantity, error)` | ParseQuantity returns the Quantity corresponding to value expressed in the unit named by unitSymbol. |
| `func (q Quantity) Add(r Quantity) (Quantity, error)` | Add returns the sum of q and r. |
| `func (q Quantity) Div(r Quantity) Quantity` | Div returns the quotient of q and r: values divide and dimension exponents subtract. |
| `func (q Quantity) IsDimensionless() bool` | IsDimensionless reports whether q has no physical dimension (all base exponents zero), such as a pure number or an angle in radians. |
| `func (q Quantity) Mul(r Quantity) Quantity` | Mul returns the product of q and r: values multiply and dimension exponents add. |
| `func (q Quantity) Pow(n int) Quantity` | Pow returns q raised to the integer power n: the value is raised to n and every dimension exponent is multiplied by n. |
| `func (q Quantity) Sub(r Quantity) (Quantity, error)` | Sub returns the difference q-r. |
| `func (q Quantity) To(unitSymbol string) (float64, error)` | To expresses q in the unit named by unitSymbol and returns the numeric value in that unit. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Curl(f func(Vec3) Vec3, p Vec3, h float64) Vec3` | Curl returns a numerical approximation of the curl ∇×f of the vector field f at the point p, using second-order central finite differences with… |
| `func Gradient(f func(Vec3) float64, p Vec3, h float64) Vec3` | Gradient returns a numerical approximation of the gradient ∇f of the scalar field f at the point p, using second-order central finite differences… |
| `func NewVec3(x, y, z float64) Vec3` | NewVec3 returns the vector with components (x, y, z). |
| `func SumSlice(v []Vec3) Vec3` | SumSlice returns the component-wise sum of every vector in v. |
| `func (a Vec3) Add(b Vec3) Vec3` | Add returns the component-wise sum a + b. |
| `func (a Vec3) Angle(b Vec3) float64` | Angle returns the unsigned angle in radians between a and b, in the range [0, π]. |
| `func (a Vec3) AngleBetween(b Vec3) float64` | AngleBetween returns the unsigned angle in radians between a and b, in the range [0, π]. |
| `func (a Vec3) Cross(b Vec3) Vec3` | Cross returns the vector (cross) product a×b, a vector orthogonal to both a and b whose magnitude equals the area of the parallelogram they span and… |
| `func (a Vec3) Distance(b Vec3) float64` | Distance returns the Euclidean distance \|a - b\| between the points a and b. |
| `func (a Vec3) Dot(b Vec3) float64` | Dot returns the scalar (inner) product a·b = aₓbₓ + a_yb_y + a_zb_z. |
| `func (a Vec3) Lerp(b Vec3, t float64) Vec3` | Lerp returns the linear interpolation (1-t)·a + t·b between a and b. |
| `func (a Vec3) Neg() Vec3` | Neg returns the vector pointing opposite to a, that is a scaled by -1. |
| `func (a Vec3) Norm() float64` | Norm returns the Euclidean length (magnitude) \|a\| = √(a·a) of the vector. |
| `func (a Vec3) Norm2() float64` | Norm2 returns the squared Euclidean length a·a. |
| `func (a Vec3) Normalize() Vec3` | Normalize returns the unit vector pointing in the direction of a. |
| `func (a Vec3) Scale(s float64) Vec3` | Scale returns the vector a with every component multiplied by the scalar s. |
| `func (a Vec3) Sub(b Vec3) Vec3` | Sub returns the component-wise difference a - b. |
| `func (a Vec3) Unit() Vec3` | Unit returns the unit vector (length 1) pointing in the direction of a. |

</details>

### Constants

`MolarMassConstant`, `ClassicalElectronRadius`, `ComptonWavelength`, `HartreeEnergy`, `MagneticFluxQuantum`, `ConductanceQuantum`, `JosephsonConstant`, `VonKlitzingConstant`, `WienDisplacement`, `FirstRadiation`, `SecondRadiation`, `ThomsonCrossSection`, `NuclearMagneton`, `BohrMagneton`, `MuonMass`, `TauMass`, `PlanckLength`, `PlanckTime`, `PlanckMass`, `PlanckTemperature`, `StandardAtmosphere`, `MolarVolumeSTP`, `SpeedOfLight`, `PlanckConstant`, `BoltzmannConstant`, `AvogadroConstant`, `ElementaryCharge`, `StandardGravity`, `ElectronVolt`, `GravitationalConstant`, `ElectronMass`, `ProtonMass`, `NeutronMass`, `VacuumPermittivity`, `VacuumPermeability`, `FineStructureConstant`, `GasConstant`, `StefanBoltzmann`, `FaradayConstant`, `RydbergConstant`, `BohrRadius`, `AtomicMassUnit`

### Variables

`DimLength`, `DimMass`, `DimTime`, `DimCurrent`, `DimTemperature`, `DimAmount`, `DimLuminosity`, `DimVelocity`, `DimAccel`, `DimForce`, `DimEnergy`, `DimPower`, `DimPressure`, `DimCharge`, `DimVoltage`, `CoulombConstant`, `ErrCenterOfMass`, `ReducedPlanck`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/physics).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
