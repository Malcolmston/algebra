package calendrical

import "fmt"

// Islamic month numbers (tabular calendar).
const (
	Muharram = iota + 1
	Safar
	RabiI
	RabiII
	JumadaI
	JumadaII
	Rajab
	Shaban
	Ramadan
	Shawwal
	DhuAlQadah
	DhuAlHijjah
)

var islamicMonthNames = [...]string{
	"", "Muharram", "Safar", "Rabi' I", "Rabi' II", "Jumada I", "Jumada II",
	"Rajab", "Sha'ban", "Ramadan", "Shawwal", "Dhu al-Qa'dah", "Dhu al-Hijjah",
}

// IslamicDate is a date in the tabular (arithmetical) Islamic calendar. This is
// the civil approximation to the observational lunar calendar; individual dates
// may differ from a local sighting by a day.
type IslamicDate struct {
	Year  int
	Month int
	Day   int
}

// IslamicLeapYear reports whether the given Islamic year is a leap year in the
// tabular calendar (11 leap years in each 30-year cycle).
func IslamicLeapYear(year int) bool {
	return floorMod(14+11*year, 30) < 11
}

// IslamicMonthName returns the English transliteration of an Islamic month
// number.
func IslamicMonthName(month int) string {
	if month < 1 || month > 12 {
		return ""
	}
	return islamicMonthNames[month]
}

// IslamicDaysInMonth returns the number of days in the given month of the given
// Islamic year: odd months have 30 days, even months 29, and the twelfth month
// has 30 in a leap year.
func IslamicDaysInMonth(year, month int) int {
	if month < 1 || month > 12 {
		return 0
	}
	if month%2 == 1 {
		return 30
	}
	if month == DhuAlHijjah && IslamicLeapYear(year) {
		return 30
	}
	return 29
}

// IslamicDaysInYear returns the number of days in the given Islamic year (354 or
// 355).
func IslamicDaysInYear(year int) int {
	if IslamicLeapYear(year) {
		return 355
	}
	return 354
}

// FixedFromIslamic converts an Islamic year, month and day to an RD fixed day.
func FixedFromIslamic(year, month, day int) int {
	return IslamicEpoch - 1 +
		(year-1)*354 +
		floorDiv(3+11*year, 30) +
		29*(month-1) +
		floorDiv(month, 2) +
		day
}

// IslamicFromFixed converts an RD fixed day to a tabular Islamic date.
func IslamicFromFixed(fixed int) IslamicDate {
	year := floorDiv(30*(fixed-IslamicEpoch)+10646, 10631)
	priorDays := fixed - FixedFromIslamic(year, Muharram, 1)
	month := floorDiv(11*priorDays+330, 325)
	day := fixed - FixedFromIslamic(year, month, 1) + 1
	return IslamicDate{Year: year, Month: month, Day: day}
}

// IslamicValid reports whether year, month and day form a valid tabular Islamic
// date.
func IslamicValid(year, month, day int) bool {
	if month < 1 || month > 12 {
		return false
	}
	return day >= 1 && day <= IslamicDaysInMonth(year, month)
}

// NewIslamic returns a validated IslamicDate, or an error if the fields do not
// describe a real date.
func NewIslamic(year, month, day int) (IslamicDate, error) {
	if !IslamicValid(year, month, day) {
		return IslamicDate{}, fmt.Errorf("%w: islamic %d-%02d-%02d", ErrInvalidDate, year, month, day)
	}
	return IslamicDate{Year: year, Month: month, Day: day}, nil
}

// Fixed returns the RD fixed day of the Islamic date.
func (d IslamicDate) Fixed() int { return FixedFromIslamic(d.Year, d.Month, d.Day) }

// JDN returns the Julian Day Number of the Islamic date.
func (d IslamicDate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the Islamic date.
func (d IslamicDate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid Islamic date.
func (d IslamicDate) Valid() bool { return IslamicValid(d.Year, d.Month, d.Day) }

// String formats the Islamic date as "Day MonthName Year AH".
func (d IslamicDate) String() string {
	return fmt.Sprintf("%d %s %d AH", d.Day, IslamicMonthName(d.Month), d.Year)
}

// IslamicNewYear returns the RD fixed day of 1 Muharram of the given Islamic
// year.
func IslamicNewYear(year int) int { return FixedFromIslamic(year, Muharram, 1) }

// RamadanStart returns the RD fixed day of 1 Ramadan of the given Islamic year
// in the tabular calendar.
func RamadanStart(year int) int { return FixedFromIslamic(year, Ramadan, 1) }
