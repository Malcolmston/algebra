package calendrical

import "fmt"

// JulianDate is a proleptic Julian-calendar date. The Julian calendar has no
// year zero: the year before 1 CE is 1 BCE, represented here as -1. Use
// JulianYearToAstronomical to convert to a continuous numbering when needed.
type JulianDate struct {
	Year  int
	Month int
	Day   int
}

// JulianYearToAstronomical converts a Julian-calendar year that omits zero
// (…, -1 = 1 BCE, 1 = 1 CE, …) to a continuous astronomical year number.
func JulianYearToAstronomical(year int) int {
	if year < 0 {
		return year + 1
	}
	return year
}

// AstronomicalToJulianYear converts a continuous astronomical year number to a
// Julian-calendar year that omits zero.
func AstronomicalToJulianYear(year int) int {
	if year <= 0 {
		return year - 1
	}
	return year
}

// JulianLeapYear reports whether the given Julian-calendar year is a leap year.
// Every fourth year is a leap year; because there is no year zero, negative
// years are leap when year mod 4 == 3.
func JulianLeapYear(year int) bool {
	if year > 0 {
		return floorMod(year, 4) == 0
	}
	return floorMod(year, 4) == 3
}

// JulianDaysInMonth returns the number of days in the given month of the given
// Julian-calendar year.
func JulianDaysInMonth(year, month int) int {
	switch month {
	case February:
		if JulianLeapYear(year) {
			return 29
		}
		return 28
	case April, June, September, November:
		return 30
	case January, March, May, July, August, October, December:
		return 31
	default:
		return 0
	}
}

// JulianDaysInYear returns the number of days in the given Julian-calendar year
// (365 or 366).
func JulianDaysInYear(year int) int {
	if JulianLeapYear(year) {
		return 366
	}
	return 365
}

// FixedFromJulian converts a Julian-calendar year, month and day to an RD fixed
// day.
func FixedFromJulian(year, month, day int) int {
	y := year
	if year < 0 {
		y = year + 1
	}
	f := JulianEpoch - 1 +
		365*(y-1) +
		floorDiv(y-1, 4) +
		floorDiv(367*month-362, 12) +
		day
	switch {
	case month <= 2:
		return f
	case JulianLeapYear(year):
		return f - 1
	default:
		return f - 2
	}
}

// JulianFromFixed converts an RD fixed day to a Julian-calendar date.
func JulianFromFixed(fixed int) JulianDate {
	approx := floorDiv(4*(fixed-JulianEpoch)+1464, 1461)
	year := approx
	if approx <= 0 {
		year = approx - 1
	}
	priorDays := fixed - FixedFromJulian(year, January, 1)
	var correction int
	switch {
	case fixed < FixedFromJulian(year, March, 1):
		correction = 0
	case JulianLeapYear(year):
		correction = 1
	default:
		correction = 2
	}
	month := floorDiv(12*(priorDays+correction)+373, 367)
	day := fixed - FixedFromJulian(year, month, 1) + 1
	return JulianDate{Year: year, Month: month, Day: day}
}

// JulianValid reports whether year, month and day form a valid Julian-calendar
// date. Year zero is rejected.
func JulianValid(year, month, day int) bool {
	if year == 0 || month < 1 || month > 12 {
		return false
	}
	return day >= 1 && day <= JulianDaysInMonth(year, month)
}

// NewJulian returns a validated JulianDate, or an error if the fields do not
// describe a real date.
func NewJulian(year, month, day int) (JulianDate, error) {
	if !JulianValid(year, month, day) {
		return JulianDate{}, fmt.Errorf("%w: julian %d-%02d-%02d", ErrInvalidDate, year, month, day)
	}
	return JulianDate{Year: year, Month: month, Day: day}, nil
}

// Fixed returns the RD fixed day of the Julian-calendar date.
func (d JulianDate) Fixed() int { return FixedFromJulian(d.Year, d.Month, d.Day) }

// JDN returns the Julian Day Number of the Julian-calendar date.
func (d JulianDate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the Julian-calendar date.
func (d JulianDate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid Julian-calendar date.
func (d JulianDate) Valid() bool { return JulianValid(d.Year, d.Month, d.Day) }

// String formats the Julian date as YYYY-MM-DD.
func (d JulianDate) String() string {
	return fmt.Sprintf("%d-%02d-%02d", d.Year, d.Month, d.Day)
}

// JulianToGregorian converts a Julian-calendar date to the equivalent
// Gregorian date for the same instant.
func JulianToGregorian(year, month, day int) GregorianDate {
	return GregorianFromFixed(FixedFromJulian(year, month, day))
}

// GregorianToJulian converts a Gregorian date to the equivalent Julian-calendar
// date for the same instant.
func GregorianToJulian(year, month, day int) JulianDate {
	return JulianFromFixed(FixedFromGregorian(year, month, day))
}
