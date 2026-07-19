package calendrical

import "fmt"

// Gregorian month numbers.
const (
	January = iota + 1
	February
	March
	April
	May
	June
	July
	August
	September
	October
	November
	December
)

var gregorianMonthNames = [...]string{
	"", "January", "February", "March", "April", "May", "June",
	"July", "August", "September", "October", "November", "December",
}

// GregorianDate is a proleptic Gregorian date. Years use astronomical
// numbering, so the year before 1 is year 0 (1 BCE) and the year before that is
// -1 (2 BCE).
type GregorianDate struct {
	Year  int
	Month int
	Day   int
}

// GregorianLeapYear reports whether the given Gregorian year is a leap year.
func GregorianLeapYear(year int) bool {
	return floorMod(year, 4) == 0 && (floorMod(year, 100) != 0 || floorMod(year, 400) == 0)
}

// GregorianMonthName returns the English name of a Gregorian month (1 = January).
func GregorianMonthName(month int) string {
	if month < 1 || month > 12 {
		return ""
	}
	return gregorianMonthNames[month]
}

// GregorianDaysInMonth returns the number of days in the given Gregorian month
// of the given year.
func GregorianDaysInMonth(year, month int) int {
	switch month {
	case February:
		if GregorianLeapYear(year) {
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

// GregorianDaysInYear returns the number of days in the given Gregorian year
// (365 or 366).
func GregorianDaysInYear(year int) int {
	if GregorianLeapYear(year) {
		return 366
	}
	return 365
}

// FixedFromGregorian converts a Gregorian year, month and day to an RD fixed
// day.
func FixedFromGregorian(year, month, day int) int {
	f := GregorianEpoch - 1 +
		365*(year-1) +
		floorDiv(year-1, 4) -
		floorDiv(year-1, 100) +
		floorDiv(year-1, 400) +
		floorDiv(367*month-362, 12) +
		day
	switch {
	case month <= 2:
		return f
	case GregorianLeapYear(year):
		return f - 1
	default:
		return f - 2
	}
}

// GregorianYearFromFixed returns the Gregorian year containing an RD fixed day.
func GregorianYearFromFixed(fixed int) int {
	d0 := fixed - GregorianEpoch
	n400 := floorDiv(d0, 146097)
	d1 := floorMod(d0, 146097)
	n100 := floorDiv(d1, 36524)
	d2 := floorMod(d1, 36524)
	n4 := floorDiv(d2, 1461)
	d3 := floorMod(d2, 1461)
	n1 := floorDiv(d3, 365)
	year := 400*n400 + 100*n100 + 4*n4 + n1
	if n100 == 4 || n1 == 4 {
		return year
	}
	return year + 1
}

// GregorianFromFixed converts an RD fixed day to a Gregorian date.
func GregorianFromFixed(fixed int) GregorianDate {
	year := GregorianYearFromFixed(fixed)
	priorDays := fixed - FixedFromGregorian(year, January, 1)
	var correction int
	switch {
	case fixed < FixedFromGregorian(year, March, 1):
		correction = 0
	case GregorianLeapYear(year):
		correction = 1
	default:
		correction = 2
	}
	month := floorDiv(12*(priorDays+correction)+373, 367)
	day := fixed - FixedFromGregorian(year, month, 1) + 1
	return GregorianDate{Year: year, Month: month, Day: day}
}

// GregorianValid reports whether year, month and day form a valid Gregorian
// date.
func GregorianValid(year, month, day int) bool {
	if month < 1 || month > 12 {
		return false
	}
	return day >= 1 && day <= GregorianDaysInMonth(year, month)
}

// NewGregorian returns a validated GregorianDate, or an error if the fields do
// not describe a real date.
func NewGregorian(year, month, day int) (GregorianDate, error) {
	if !GregorianValid(year, month, day) {
		return GregorianDate{}, fmt.Errorf("%w: gregorian %04d-%02d-%02d", ErrInvalidDate, year, month, day)
	}
	return GregorianDate{Year: year, Month: month, Day: day}, nil
}

// Fixed returns the RD fixed day of the Gregorian date.
func (d GregorianDate) Fixed() int { return FixedFromGregorian(d.Year, d.Month, d.Day) }

// JDN returns the Julian Day Number of the Gregorian date.
func (d GregorianDate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the Gregorian date.
func (d GregorianDate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid Gregorian date.
func (d GregorianDate) Valid() bool { return GregorianValid(d.Year, d.Month, d.Day) }

// String formats the date in ISO extended form (YYYY-MM-DD), using a sign for
// non-positive years.
func (d GregorianDate) String() string {
	return fmt.Sprintf("%+05d-%02d-%02d", d.Year, d.Month, d.Day)
}

// GregorianDayNumber returns the ordinal day-of-year (1 = 1 January) of the
// given Gregorian date.
func GregorianDayNumber(year, month, day int) int {
	return FixedFromGregorian(year, month, day) - FixedFromGregorian(year-1, December, 31)
}

// GregorianDaysRemaining returns the number of days remaining in the year after
// the given Gregorian date.
func GregorianDaysRemaining(year, month, day int) int {
	return FixedFromGregorian(year, December, 31) - FixedFromGregorian(year, month, day)
}

// GregorianNewYear returns the RD fixed day of 1 January of the given year.
func GregorianNewYear(year int) int { return FixedFromGregorian(year, January, 1) }

// GregorianYearEnd returns the RD fixed day of 31 December of the given year.
func GregorianYearEnd(year int) int { return FixedFromGregorian(year, December, 31) }

// GregorianDateDifference returns the number of days from the first Gregorian
// date to the second (positive when the second date is later).
func GregorianDateDifference(y1, m1, d1, y2, m2, d2 int) int {
	return FixedFromGregorian(y2, m2, d2) - FixedFromGregorian(y1, m1, d1)
}
