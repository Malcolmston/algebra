package calendrical

import "math"

// Season identifies one of the four astronomical season boundaries.
type Season int

// The four astronomical seasons, named for the event that begins them in the
// northern hemisphere.
const (
	SpringEquinox Season = iota
	SummerSolstice
	AutumnEquinox
	WinterSolstice
)

// MoonPhase identifies one of the four principal lunar phases.
type MoonPhase int

// The four principal phases of the Moon.
const (
	NewMoon MoonPhase = iota
	FirstQuarter
	FullMoon
	LastQuarter
)

const degToRad = math.Pi / 180.0

func sinDeg(x float64) float64 { return math.Sin(x * degToRad) }
func cosDeg(x float64) float64 { return math.Cos(x * degToRad) }

// MeanSeasonJD returns the Julian Ephemeris Day of the mean (uncorrected)
// equinox or solstice for the given Gregorian year, using Meeus's polynomial
// (Astronomical Algorithms, chapter 27) valid for years -1000 to +3000.
func MeanSeasonJD(year int, s Season) float64 {
	if year >= 1000 {
		y := (float64(year) - 2000.0) / 1000.0
		switch s {
		case SpringEquinox:
			return 2451623.80984 + 365242.37404*y + 0.05169*y*y - 0.00411*y*y*y - 0.00057*y*y*y*y
		case SummerSolstice:
			return 2451716.56767 + 365241.62603*y + 0.00325*y*y + 0.00888*y*y*y - 0.00030*y*y*y*y
		case AutumnEquinox:
			return 2451810.21715 + 365242.01767*y - 0.11575*y*y + 0.00337*y*y*y + 0.00078*y*y*y*y
		default: // WinterSolstice
			return 2451900.05952 + 365242.74049*y - 0.06223*y*y - 0.00823*y*y*y + 0.00032*y*y*y*y
		}
	}
	y := float64(year) / 1000.0
	switch s {
	case SpringEquinox:
		return 1721139.29189 + 365242.13740*y + 0.06134*y*y + 0.00111*y*y*y - 0.00071*y*y*y*y
	case SummerSolstice:
		return 1721233.25401 + 365241.72562*y - 0.05323*y*y + 0.00907*y*y*y + 0.00025*y*y*y*y
	case AutumnEquinox:
		return 1721325.70455 + 365242.49558*y - 0.11677*y*y - 0.00297*y*y*y + 0.00074*y*y*y*y
	default: // WinterSolstice
		return 1721414.39987 + 365242.88257*y - 0.00769*y*y - 0.00933*y*y*y - 0.00006*y*y*y*y
	}
}

// seasonTerms holds the 24 periodic terms (A, B, C in degrees) used to correct
// the mean season to the apparent one.
var seasonTerms = [24][3]float64{
	{485, 324.96, 1934.136}, {203, 337.23, 32964.467}, {199, 342.08, 20.186},
	{182, 27.85, 445267.112}, {156, 73.14, 45036.886}, {136, 171.52, 22518.443},
	{77, 222.54, 65928.934}, {74, 296.72, 3034.906}, {70, 243.58, 9037.513},
	{58, 119.81, 33718.147}, {52, 297.17, 150.678}, {50, 21.02, 2281.226},
	{45, 247.54, 29929.562}, {44, 325.15, 31555.956}, {29, 60.93, 4443.417},
	{18, 155.12, 67555.328}, {17, 288.79, 4562.452}, {16, 198.04, 62894.029},
	{14, 199.76, 31436.921}, {12, 95.39, 14577.848}, {12, 287.11, 31931.756},
	{12, 320.81, 34777.259}, {9, 227.73, 1222.114}, {8, 15.45, 16859.074},
}

// SeasonJD returns the Julian Ephemeris Day of the apparent equinox or solstice
// for the given Gregorian year, applying Meeus's periodic corrections to the
// mean value. The result is accurate to roughly a minute for years near the
// present and to within about an hour over the -1000 to +3000 range.
func SeasonJD(year int, s Season) float64 {
	jd0 := MeanSeasonJD(year, s)
	t := (jd0 - 2451545.0) / 36525.0
	w := 35999.373*t - 2.47
	dl := 1.0 + 0.0334*cosDeg(w) + 0.0007*cosDeg(2*w)
	var sum float64
	for _, term := range seasonTerms {
		sum += term[0] * cosDeg(term[1]+term[2]*t)
	}
	return jd0 + (0.00001*sum)/dl
}

// SeasonFixed returns the RD fixed day (in Terrestrial Time, at Greenwich) of
// the apparent equinox or solstice for the given Gregorian year.
func SeasonFixed(year int, s Season) int {
	return FixedFromMoment(MomentFromJD(SeasonJD(year, s)))
}

// SpringEquinoxJD returns the Julian Ephemeris Day of the March equinox of the
// given Gregorian year.
func SpringEquinoxJD(year int) float64 { return SeasonJD(year, SpringEquinox) }

// SummerSolsticeJD returns the Julian Ephemeris Day of the June solstice of the
// given Gregorian year.
func SummerSolsticeJD(year int) float64 { return SeasonJD(year, SummerSolstice) }

// AutumnEquinoxJD returns the Julian Ephemeris Day of the September equinox of
// the given Gregorian year.
func AutumnEquinoxJD(year int) float64 { return SeasonJD(year, AutumnEquinox) }

// WinterSolsticeJD returns the Julian Ephemeris Day of the December solstice of
// the given Gregorian year.
func WinterSolsticeJD(year int) float64 { return SeasonJD(year, WinterSolstice) }

// MeanNewMoonJD returns the Julian Ephemeris Day of the mean New Moon for the
// given lunation index k, where k = 0 corresponds to the New Moon of 6 January
// 2000. Non-integer k of the form k+0.25, k+0.5, k+0.75 select the mean first
// quarter, full moon and last quarter respectively (Meeus, chapter 49).
func MeanNewMoonJD(k float64) float64 {
	t := k / 1236.85
	return 2451550.09766 + 29.530588861*k +
		0.00015437*t*t -
		0.000000150*t*t*t +
		0.00000000073*t*t*t*t
}

// LunationNumber returns the integer lunation index k whose mean New Moon is
// nearest to the given Julian Ephemeris Day. It is the inverse of MeanNewMoonJD
// rounded to the closest New Moon.
func LunationNumber(jd float64) int {
	return int(math.Round((jd - 2451550.09766) / 29.530588861))
}

// MoonPhaseJD returns the Julian Ephemeris Day of the requested principal phase
// for the lunation nearest index k. It applies the principal periodic
// corrections from Meeus chapter 49 and is accurate to a few minutes for dates
// near the present era.
func MoonPhaseJD(k int, phase MoonPhase) float64 {
	kk := float64(k)
	switch phase {
	case FirstQuarter:
		kk += 0.25
	case FullMoon:
		kk += 0.5
	case LastQuarter:
		kk += 0.75
	}
	t := kk / 1236.85
	jde := MeanNewMoonJD(kk)

	// Sun's mean anomaly, Moon's mean anomaly, Moon's argument of latitude.
	m := 2.5534 + 29.10535670*kk - 0.0000014*t*t - 0.00000011*t*t*t
	mp := 201.5643 + 385.81693528*kk + 0.0107582*t*t + 0.00001238*t*t*t - 0.000000058*t*t*t*t
	f := 160.7108 + 390.67050284*kk - 0.0016118*t*t - 0.00000227*t*t*t + 0.000000011*t*t*t*t
	omega := 124.7746 - 1.56375588*kk + 0.0020672*t*t + 0.00000215*t*t*t
	e := 1.0 - 0.002516*t - 0.0000074*t*t

	var corr float64
	switch phase {
	case NewMoon:
		corr = -0.40720*sinDeg(mp) +
			0.17241*e*sinDeg(m) +
			0.01608*sinDeg(2*mp) +
			0.01039*sinDeg(2*f) +
			0.00739*e*sinDeg(mp-m) -
			0.00514*e*sinDeg(mp+m) +
			0.00208*e*e*sinDeg(2*m) -
			0.00111*sinDeg(mp-2*f) -
			0.00057*sinDeg(mp+2*f) +
			0.00056*e*sinDeg(2*mp+m) -
			0.00042*sinDeg(3*mp) +
			0.00042*e*sinDeg(m+2*f) +
			0.00038*e*sinDeg(m-2*f) -
			0.00024*e*sinDeg(2*mp-m) -
			0.00017*sinDeg(omega) -
			0.00007*sinDeg(mp+2*m) +
			0.00004*sinDeg(2*mp-2*f) +
			0.00004*sinDeg(3*m) +
			0.00003*sinDeg(mp+m-2*f) +
			0.00003*sinDeg(2*mp+2*f) -
			0.00003*sinDeg(mp+m+2*f) +
			0.00003*sinDeg(mp-m+2*f) -
			0.00002*sinDeg(mp-m-2*f) -
			0.00002*sinDeg(3*mp+m) +
			0.00002*sinDeg(4*mp)
	case FullMoon:
		corr = -0.40614*sinDeg(mp) +
			0.17302*e*sinDeg(m) +
			0.01614*sinDeg(2*mp) +
			0.01043*sinDeg(2*f) +
			0.00734*e*sinDeg(mp-m) -
			0.00514*e*sinDeg(mp+m) +
			0.00209*e*e*sinDeg(2*m) -
			0.00111*sinDeg(mp-2*f) -
			0.00057*sinDeg(mp+2*f) +
			0.00056*e*sinDeg(2*mp+m) -
			0.00042*sinDeg(3*mp) +
			0.00042*e*sinDeg(m+2*f) +
			0.00038*e*sinDeg(m-2*f) -
			0.00024*e*sinDeg(2*mp-m) -
			0.00017*sinDeg(omega) -
			0.00007*sinDeg(mp+2*m) +
			0.00004*sinDeg(2*mp-2*f) +
			0.00004*sinDeg(3*m) +
			0.00003*sinDeg(mp+m-2*f) +
			0.00003*sinDeg(2*mp+2*f) -
			0.00003*sinDeg(mp+m+2*f) +
			0.00003*sinDeg(mp-m+2*f) -
			0.00002*sinDeg(mp-m-2*f) -
			0.00002*sinDeg(3*mp+m) +
			0.00002*sinDeg(4*mp)
	default: // FirstQuarter, LastQuarter
		corr = -0.62801*sinDeg(mp) +
			0.17172*e*sinDeg(m) -
			0.01183*e*sinDeg(mp+m) +
			0.00862*sinDeg(2*mp) +
			0.00804*sinDeg(2*f) +
			0.00454*e*sinDeg(mp-m) +
			0.00204*e*e*sinDeg(2*m) -
			0.00180*sinDeg(mp-2*f) -
			0.00070*sinDeg(mp+2*f) -
			0.00040*sinDeg(3*mp) -
			0.00034*e*sinDeg(2*mp-m) +
			0.00032*e*sinDeg(m+2*f) +
			0.00032*e*sinDeg(m-2*f) -
			0.00028*e*e*sinDeg(mp+2*m) +
			0.00027*e*sinDeg(2*mp+m) -
			0.00017*sinDeg(omega) -
			0.00005*sinDeg(mp-m-2*f) +
			0.00004*sinDeg(2*mp+2*f) -
			0.00004*sinDeg(mp+m+2*f) +
			0.00004*sinDeg(mp-2*m) +
			0.00003*sinDeg(mp+m-2*f) +
			0.00003*sinDeg(3*m) +
			0.00002*sinDeg(2*mp-2*f) +
			0.00002*sinDeg(mp-m+2*f) -
			0.00002*sinDeg(3*mp+m)
		w := 0.00306 - 0.00038*e*cosDeg(m) + 0.00026*cosDeg(mp) -
			0.00002*cosDeg(mp-m) + 0.00002*cosDeg(mp+m) + 0.00002*cosDeg(2*f)
		if phase == FirstQuarter {
			corr += w
		} else {
			corr -= w
		}
	}
	return jde + corr
}

// MoonPhaseFixed returns the RD fixed day (in Terrestrial Time) of the requested
// principal phase for the lunation nearest index k.
func MoonPhaseFixed(k int, phase MoonPhase) int {
	return FixedFromMoment(MomentFromJD(MoonPhaseJD(k, phase)))
}

// NextMoonPhaseJD returns the Julian Ephemeris Day of the first occurrence of
// the requested principal phase on or after the given Julian Ephemeris Day.
func NextMoonPhaseJD(jd float64, phase MoonPhase) float64 {
	k := LunationNumber(jd) - 2
	for {
		p := MoonPhaseJD(k, phase)
		if p >= jd {
			return p
		}
		k++
	}
}

// MeanLunarMonth is the mean length of the synodic month (New Moon to New Moon)
// in days.
const MeanLunarMonth = 29.530588861

// MeanTropicalYear is the mean length of the tropical year in days.
const MeanTropicalYear = 365.242189
