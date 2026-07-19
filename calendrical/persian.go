package calendrical

import "fmt"

// Persian (Jalali) month numbers.
const (
	Farvardin = iota + 1
	Ordibehesht
	Khordad
	Tir
	Mordad
	Shahrivar
	Mehr
	Aban
	Azar
	Dey
	Bahman
	Esfand
)

var persianMonthNames = [...]string{
	"", "Farvardin", "Ordibehesht", "Khordad", "Tir", "Mordad", "Shahrivar",
	"Mehr", "Aban", "Azar", "Dey", "Bahman", "Esfand",
}

// PersianDate is a date in the arithmetical Persian (Jalali) calendar, which
// approximates the astronomical Persian calendar using a 2820-year intercalation
// cycle.
type PersianDate struct {
	Year  int
	Month int
	Day   int
}

// PersianLeapYear reports whether the given Persian year is a leap year in the
// arithmetical (2820-year) calendar.
func PersianLeapYear(year int) bool {
	var y int
	if year > 0 {
		y = year - 474
	} else {
		y = year - 473
	}
	y2 := floorMod(y, 2820) + 474
	return floorMod((y2+38)*31, 128) < 31
}

// PersianMonthName returns the transliterated name of a Persian month number.
func PersianMonthName(month int) string {
	if month < 1 || month > 12 {
		return ""
	}
	return persianMonthNames[month]
}

// PersianDaysInMonth returns the number of days in the given month of the given
// Persian year: 31 for the first six months, 30 for months 7–11, and 29 (or 30
// in a leap year) for Esfand.
func PersianDaysInMonth(year, month int) int {
	switch {
	case month < 1 || month > 12:
		return 0
	case month <= 6:
		return 31
	case month <= 11:
		return 30
	default:
		if PersianLeapYear(year) {
			return 30
		}
		return 29
	}
}

// PersianDaysInYear returns the number of days in the given Persian year (365 or
// 366).
func PersianDaysInYear(year int) int {
	if PersianLeapYear(year) {
		return 366
	}
	return 365
}

// FixedFromPersian converts a Persian (Jalali) year, month and day to an RD
// fixed day.
func FixedFromPersian(year, month, day int) int {
	var y int
	if year > 0 {
		y = year - 474
	} else {
		y = year - 473
	}
	y2 := floorMod(y, 2820) + 474
	var monthDays int
	if month <= 7 {
		monthDays = 31 * (month - 1)
	} else {
		monthDays = 30*(month-1) + 6
	}
	return PersianEpoch - 1 +
		1029983*floorDiv(y, 2820) +
		365*(y2-1) +
		floorDiv(31*y2-5, 128) +
		monthDays +
		day
}

// PersianYearFromFixed returns the Persian year containing an RD fixed day.
func PersianYearFromFixed(fixed int) int {
	d0 := fixed - FixedFromPersian(475, Farvardin, 1)
	n2820 := floorDiv(d0, 1029983)
	d1 := floorMod(d0, 1029983)
	var y2820 int
	if d1 == 1029982 {
		y2820 = 2820
	} else {
		y2820 = floorDiv(128*d1+46878, 46751)
	}
	year := 474 + 2820*n2820 + y2820
	if year > 0 {
		return year
	}
	return year - 1
}

// PersianFromFixed converts an RD fixed day to a Persian (Jalali) date.
func PersianFromFixed(fixed int) PersianDate {
	year := PersianYearFromFixed(fixed)
	dayOfYear := 1 + fixed - FixedFromPersian(year, Farvardin, 1)
	var month int
	if dayOfYear <= 186 {
		month = ceilDiv(dayOfYear, 31)
	} else {
		month = ceilDiv(dayOfYear-6, 30)
	}
	day := fixed - FixedFromPersian(year, month, 1) + 1
	return PersianDate{Year: year, Month: month, Day: day}
}

// ceilDiv returns the ceiling of a/b for positive b.
func ceilDiv(a, b int) int {
	return floorDiv(a+b-1, b)
}

// PersianValid reports whether year, month and day form a valid Persian date.
func PersianValid(year, month, day int) bool {
	if month < 1 || month > 12 {
		return false
	}
	return day >= 1 && day <= PersianDaysInMonth(year, month)
}

// NewPersian returns a validated PersianDate, or an error if the fields do not
// describe a real date.
func NewPersian(year, month, day int) (PersianDate, error) {
	if !PersianValid(year, month, day) {
		return PersianDate{}, fmt.Errorf("%w: persian %d-%02d-%02d", ErrInvalidDate, year, month, day)
	}
	return PersianDate{Year: year, Month: month, Day: day}, nil
}

// Fixed returns the RD fixed day of the Persian date.
func (d PersianDate) Fixed() int { return FixedFromPersian(d.Year, d.Month, d.Day) }

// JDN returns the Julian Day Number of the Persian date.
func (d PersianDate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the Persian date.
func (d PersianDate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid Persian date.
func (d PersianDate) Valid() bool { return PersianValid(d.Year, d.Month, d.Day) }

// String formats the Persian date as YYYY-MM-DD.
func (d PersianDate) String() string {
	return fmt.Sprintf("%d-%02d-%02d", d.Year, d.Month, d.Day)
}

// Nowruz returns the RD fixed day of Nowruz (the Persian new year, 1 Farvardin)
// occurring in the given Gregorian year.
func Nowruz(gYear int) int {
	persianYear := gYear - GregorianYearFromFixed(PersianEpoch) + 1
	return FixedFromPersian(persianYear, Farvardin, 1)
}
