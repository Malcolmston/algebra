package calendrical

import "fmt"

// Hebrew month numbers. The Hebrew year begins in Tishri (month 7); months 1–6
// (Nisan through Elul) fall in the second half of the civil year. In a leap
// year Adar is split into Adar I (month 12) and Adar II (month 13).
const (
	Nisan = iota + 1
	Iyyar
	Sivan
	Tammuz
	Av
	Elul
	Tishri
	Marheshvan
	Kislev
	Tevet
	Shevat
	Adar
	AdarII
)

var hebrewMonthNames = [...]string{
	"", "Nisan", "Iyyar", "Sivan", "Tammuz", "Av", "Elul", "Tishri",
	"Marheshvan", "Kislev", "Tevet", "Shevat", "Adar", "Adar II",
}

// HebrewDate is a date in the Hebrew (Jewish) calendar.
type HebrewDate struct {
	Year  int
	Month int
	Day   int
}

// HebrewLeapYear reports whether the given Hebrew year is embolismic (has a
// thirteenth month), under the 19-year Metonic cycle.
func HebrewLeapYear(year int) bool {
	return floorMod(7*year+1, 19) < 7
}

// HebrewMonthName returns the English transliteration of a Hebrew month number.
// In a leap year month 12 is "Adar I".
func HebrewMonthName(year, month int) string {
	if month < 1 || month > 13 {
		return ""
	}
	if month == Adar && HebrewLeapYear(year) {
		return "Adar I"
	}
	return hebrewMonthNames[month]
}

// LastMonthOfHebrewYear returns the number of the final month of the Hebrew
// year: 12 in a common year, 13 in a leap year.
func LastMonthOfHebrewYear(year int) int {
	if HebrewLeapYear(year) {
		return 13
	}
	return 12
}

// HebrewCalendarElapsedDays returns the number of days elapsed from the Hebrew
// epoch to the mean conjunction (molad) of Tishri of the given year, applying
// the traditional postponement (dehiyyah) of the new year.
func HebrewCalendarElapsedDays(year int) int {
	monthsElapsed := floorDiv(235*year-234, 19)
	partsElapsed := 12084 + 13753*monthsElapsed
	day := 29*monthsElapsed + floorDiv(partsElapsed, 25920)
	if floorMod(3*(day+1), 7) < 3 {
		return day + 1
	}
	return day
}

// HebrewNewYearDelay returns the year-length correction (0, 1 or 2 days) applied
// to Rosh ha-Shanah so that the year has a permissible length.
func HebrewNewYearDelay(year int) int {
	ny0 := HebrewCalendarElapsedDays(year - 1)
	ny1 := HebrewCalendarElapsedDays(year)
	ny2 := HebrewCalendarElapsedDays(year + 1)
	switch {
	case ny2-ny1 == 356:
		return 2
	case ny1-ny0 == 382:
		return 1
	default:
		return 0
	}
}

// HebrewNewYear returns the RD fixed day of Rosh ha-Shanah (1 Tishri) of the
// given Hebrew year.
func HebrewNewYear(year int) int {
	return HebrewEpoch + HebrewCalendarElapsedDays(year) + HebrewNewYearDelay(year)
}

// DaysInHebrewYear returns the length in days of the given Hebrew year (353,
// 354, 355, 383, 384 or 385).
func DaysInHebrewYear(year int) int {
	return HebrewNewYear(year+1) - HebrewNewYear(year)
}

// LongMarheshvan reports whether Marheshvan has 30 days in the given Hebrew
// year (a "complete" year).
func LongMarheshvan(year int) bool {
	n := DaysInHebrewYear(year)
	return n == 355 || n == 385
}

// ShortKislev reports whether Kislev has 29 days in the given Hebrew year (a
// "deficient" year).
func ShortKislev(year int) bool {
	n := DaysInHebrewYear(year)
	return n == 353 || n == 383
}

// LastDayOfHebrewMonth returns the number of days (29 or 30) in the given month
// of the given Hebrew year.
func LastDayOfHebrewMonth(year, month int) int {
	switch month {
	case Iyyar, Tammuz, Elul, Tevet, AdarII:
		return 29
	case Adar:
		if !HebrewLeapYear(year) {
			return 29
		}
		return 30
	case Marheshvan:
		if !LongMarheshvan(year) {
			return 29
		}
		return 30
	case Kislev:
		if ShortKislev(year) {
			return 29
		}
		return 30
	default:
		return 30
	}
}

// HebrewDaysInMonth is an alias for LastDayOfHebrewMonth.
func HebrewDaysInMonth(year, month int) int { return LastDayOfHebrewMonth(year, month) }

// FixedFromHebrew converts a Hebrew year, month and day to an RD fixed day.
func FixedFromHebrew(year, month, day int) int {
	total := HebrewNewYear(year) + day - 1
	if month < Tishri {
		for m := Tishri; m <= LastMonthOfHebrewYear(year); m++ {
			total += LastDayOfHebrewMonth(year, m)
		}
		for m := Nisan; m < month; m++ {
			total += LastDayOfHebrewMonth(year, m)
		}
	} else {
		for m := Tishri; m < month; m++ {
			total += LastDayOfHebrewMonth(year, m)
		}
	}
	return total
}

// HebrewFromFixed converts an RD fixed day to a Hebrew date.
func HebrewFromFixed(fixed int) HebrewDate {
	approx := floorDiv(98496*(fixed-HebrewEpoch), 35975351) + 1
	year := approx - 1
	for HebrewNewYear(year+1) <= fixed {
		year++
	}
	month := Tishri
	if fixed >= FixedFromHebrew(year, Nisan, 1) {
		month = Nisan
	}
	for fixed > FixedFromHebrew(year, month, LastDayOfHebrewMonth(year, month)) {
		month++
	}
	day := fixed - FixedFromHebrew(year, month, 1) + 1
	return HebrewDate{Year: year, Month: month, Day: day}
}

// HebrewValid reports whether year, month and day form a valid Hebrew date.
func HebrewValid(year, month, day int) bool {
	if month < 1 || month > LastMonthOfHebrewYear(year) {
		return false
	}
	if !HebrewLeapYear(year) && month == AdarII {
		return false
	}
	return day >= 1 && day <= LastDayOfHebrewMonth(year, month)
}

// NewHebrew returns a validated HebrewDate, or an error if the fields do not
// describe a real Hebrew date.
func NewHebrew(year, month, day int) (HebrewDate, error) {
	if !HebrewValid(year, month, day) {
		return HebrewDate{}, fmt.Errorf("%w: hebrew %d-%02d-%02d", ErrInvalidDate, year, month, day)
	}
	return HebrewDate{Year: year, Month: month, Day: day}, nil
}

// Fixed returns the RD fixed day of the Hebrew date.
func (d HebrewDate) Fixed() int { return FixedFromHebrew(d.Year, d.Month, d.Day) }

// JDN returns the Julian Day Number of the Hebrew date.
func (d HebrewDate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the Hebrew date.
func (d HebrewDate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid Hebrew date.
func (d HebrewDate) Valid() bool { return HebrewValid(d.Year, d.Month, d.Day) }

// String formats the Hebrew date as "Day MonthName Year".
func (d HebrewDate) String() string {
	return fmt.Sprintf("%d %s %d", d.Day, HebrewMonthName(d.Year, d.Month), d.Year)
}

// HebrewMonthsElapsed returns the number of lunar months from the Hebrew epoch
// to the start (molad) of the given month of the given Hebrew year.
func HebrewMonthsElapsed(year, month int) int {
	y := year
	if month < Tishri {
		y = year + 1
	}
	return (month - Tishri) + floorDiv(235*y-234, 19)
}

// HebrewMolad returns the RD moment (a real number of days since RD 0) of the
// mean lunar conjunction — the molad — that begins the given Hebrew month.
func HebrewMolad(year, month int) float64 {
	monthsElapsed := float64(HebrewMonthsElapsed(year, month))
	// Mean synodic month: 29 days, 12 hours, 793 parts (1 part = 1/1080 hour).
	const synodic = 29.0 + 12.0/24.0 + 793.0/25920.0
	return float64(HebrewEpoch) - 876.0/25920.0 + monthsElapsed*synodic
}

// HebrewYearFromGregorian returns the Hebrew year whose autumn new year
// (Rosh ha-Shanah) falls in the given Gregorian year.
func HebrewYearFromGregorian(gYear int) int {
	return gYear - GregorianYearFromFixed(HebrewEpoch) + 1
}

// RoshHaShanah returns the RD fixed day of Rosh ha-Shanah (1 Tishri) occurring
// in the given Gregorian year.
func RoshHaShanah(gYear int) int {
	return FixedFromHebrew(HebrewYearFromGregorian(gYear), Tishri, 1)
}

// YomKippur returns the RD fixed day of Yom Kippur (10 Tishri) occurring in the
// given Gregorian year.
func YomKippur(gYear int) int {
	return FixedFromHebrew(HebrewYearFromGregorian(gYear), Tishri, 10)
}

// Sukkot returns the RD fixed day of the first day of Sukkot (15 Tishri)
// occurring in the given Gregorian year.
func Sukkot(gYear int) int {
	return FixedFromHebrew(HebrewYearFromGregorian(gYear), Tishri, 15)
}

// HanukkahStart returns the RD fixed day of the first day of Hanukkah
// (25 Kislev) occurring in the given Gregorian year.
func HanukkahStart(gYear int) int {
	return FixedFromHebrew(HebrewYearFromGregorian(gYear), Kislev, 25)
}

// Passover returns the RD fixed day of the first day of Passover (15 Nisan)
// occurring in the given Gregorian year.
func Passover(gYear int) int {
	return FixedFromHebrew(HebrewYearFromGregorian(gYear)-1, Nisan, 15)
}

// Shavuot returns the RD fixed day of Shavuot (6 Sivan) occurring in the given
// Gregorian year.
func Shavuot(gYear int) int {
	return FixedFromHebrew(HebrewYearFromGregorian(gYear)-1, Sivan, 6)
}

// Purim returns the RD fixed day of Purim occurring in the given Gregorian
// year. Purim falls on 14 Adar in a common year and 14 Adar II in a leap year.
func Purim(gYear int) int {
	hy := HebrewYearFromGregorian(gYear) - 1
	return FixedFromHebrew(hy, LastMonthOfHebrewYear(hy), 14)
}

// TuBShevat returns the RD fixed day of Tu BiShvat (15 Shevat) occurring in the
// given Gregorian year.
func TuBShevat(gYear int) int {
	return FixedFromHebrew(HebrewYearFromGregorian(gYear)-1, Shevat, 15)
}
