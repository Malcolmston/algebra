package calendrical

import "fmt"

// Coptic month numbers. The Coptic year has twelve 30-day months followed by a
// short thirteenth month (the epagomenae) of five or six days.
const (
	Thoout = iota + 1
	Paope
	Athor
	Koiak
	Tobe
	Meshir
	Paremotep
	Parmoute
	Pashons
	Paone
	Epep
	Mesore
	Epagomenae
)

var copticMonthNames = [...]string{
	"", "Thoout", "Paope", "Athor", "Koiak", "Tobe", "Meshir", "Paremotep",
	"Parmoute", "Pashons", "Paone", "Epep", "Mesore", "Epagomenae",
}

var ethiopicMonthNames = [...]string{
	"", "Maskaram", "Teqemt", "Hedar", "Takhsas", "Ter", "Yakatit", "Magabit",
	"Miyazya", "Genbot", "Sane", "Hamle", "Nahase", "Paguemen",
}

// CopticDate is a date in the Coptic (Alexandrian) calendar.
type CopticDate struct {
	Year  int
	Month int
	Day   int
}

// EthiopicDate is a date in the Ethiopic calendar, which shares the Coptic
// structure but uses a different epoch.
type EthiopicDate struct {
	Year  int
	Month int
	Day   int
}

// CopticLeapYear reports whether the given Coptic year is a leap year (every
// fourth year, with the leap day added to the thirteenth month).
func CopticLeapYear(year int) bool {
	return floorMod(year, 4) == 3
}

// EthiopicLeapYear reports whether the given Ethiopic year is a leap year.
func EthiopicLeapYear(year int) bool {
	return floorMod(year, 4) == 3
}

// CopticMonthName returns the transliterated name of a Coptic month number.
func CopticMonthName(month int) string {
	if month < 1 || month > 13 {
		return ""
	}
	return copticMonthNames[month]
}

// EthiopicMonthName returns the transliterated name of an Ethiopic month number.
func EthiopicMonthName(month int) string {
	if month < 1 || month > 13 {
		return ""
	}
	return ethiopicMonthNames[month]
}

// CopticDaysInMonth returns the number of days in the given month of the given
// Coptic year (30 for months 1–12; 5 or 6 for the epagomenae).
func CopticDaysInMonth(year, month int) int {
	switch {
	case month < 1 || month > 13:
		return 0
	case month <= 12:
		return 30
	case CopticLeapYear(year):
		return 6
	default:
		return 5
	}
}

// EthiopicDaysInMonth returns the number of days in the given month of the given
// Ethiopic year.
func EthiopicDaysInMonth(year, month int) int {
	switch {
	case month < 1 || month > 13:
		return 0
	case month <= 12:
		return 30
	case EthiopicLeapYear(year):
		return 6
	default:
		return 5
	}
}

// CopticDaysInYear returns the number of days in the given Coptic year (365 or
// 366).
func CopticDaysInYear(year int) int {
	if CopticLeapYear(year) {
		return 366
	}
	return 365
}

// EthiopicDaysInYear returns the number of days in the given Ethiopic year.
func EthiopicDaysInYear(year int) int {
	if EthiopicLeapYear(year) {
		return 366
	}
	return 365
}

// FixedFromCoptic converts a Coptic year, month and day to an RD fixed day.
func FixedFromCoptic(year, month, day int) int {
	return CopticEpoch - 1 +
		365*(year-1) +
		floorDiv(year, 4) +
		30*(month-1) +
		day
}

// CopticFromFixed converts an RD fixed day to a Coptic date.
func CopticFromFixed(fixed int) CopticDate {
	year := floorDiv(4*(fixed-CopticEpoch)+1463, 1461)
	month := floorDiv(fixed-FixedFromCoptic(year, Thoout, 1), 30) + 1
	day := fixed - FixedFromCoptic(year, month, 1) + 1
	return CopticDate{Year: year, Month: month, Day: day}
}

// FixedFromEthiopic converts an Ethiopic year, month and day to an RD fixed day.
func FixedFromEthiopic(year, month, day int) int {
	return EthiopicEpoch + (FixedFromCoptic(year, month, day) - CopticEpoch)
}

// EthiopicFromFixed converts an RD fixed day to an Ethiopic date.
func EthiopicFromFixed(fixed int) EthiopicDate {
	c := CopticFromFixed(fixed + (CopticEpoch - EthiopicEpoch))
	return EthiopicDate{Year: c.Year, Month: c.Month, Day: c.Day}
}

// CopticValid reports whether year, month and day form a valid Coptic date.
func CopticValid(year, month, day int) bool {
	if month < 1 || month > 13 {
		return false
	}
	return day >= 1 && day <= CopticDaysInMonth(year, month)
}

// EthiopicValid reports whether year, month and day form a valid Ethiopic date.
func EthiopicValid(year, month, day int) bool {
	if month < 1 || month > 13 {
		return false
	}
	return day >= 1 && day <= EthiopicDaysInMonth(year, month)
}

// NewCoptic returns a validated CopticDate, or an error if the fields do not
// describe a real date.
func NewCoptic(year, month, day int) (CopticDate, error) {
	if !CopticValid(year, month, day) {
		return CopticDate{}, fmt.Errorf("%w: coptic %d-%02d-%02d", ErrInvalidDate, year, month, day)
	}
	return CopticDate{Year: year, Month: month, Day: day}, nil
}

// NewEthiopic returns a validated EthiopicDate, or an error if the fields do not
// describe a real date.
func NewEthiopic(year, month, day int) (EthiopicDate, error) {
	if !EthiopicValid(year, month, day) {
		return EthiopicDate{}, fmt.Errorf("%w: ethiopic %d-%02d-%02d", ErrInvalidDate, year, month, day)
	}
	return EthiopicDate{Year: year, Month: month, Day: day}, nil
}

// Fixed returns the RD fixed day of the Coptic date.
func (d CopticDate) Fixed() int { return FixedFromCoptic(d.Year, d.Month, d.Day) }

// JDN returns the Julian Day Number of the Coptic date.
func (d CopticDate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the Coptic date.
func (d CopticDate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid Coptic date.
func (d CopticDate) Valid() bool { return CopticValid(d.Year, d.Month, d.Day) }

// String formats the Coptic date as "Day MonthName Year".
func (d CopticDate) String() string {
	return fmt.Sprintf("%d %s %d", d.Day, CopticMonthName(d.Month), d.Year)
}

// Fixed returns the RD fixed day of the Ethiopic date.
func (d EthiopicDate) Fixed() int { return FixedFromEthiopic(d.Year, d.Month, d.Day) }

// JDN returns the Julian Day Number of the Ethiopic date.
func (d EthiopicDate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the Ethiopic date.
func (d EthiopicDate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid Ethiopic date.
func (d EthiopicDate) Valid() bool { return EthiopicValid(d.Year, d.Month, d.Day) }

// String formats the Ethiopic date as "Day MonthName Year".
func (d EthiopicDate) String() string {
	return fmt.Sprintf("%d %s %d", d.Day, EthiopicMonthName(d.Month), d.Year)
}
