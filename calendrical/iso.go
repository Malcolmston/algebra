package calendrical

import "fmt"

// ISODate is an ISO-8601 week date: a week-numbering year, an ISO week number
// (1–53) and a day of the week numbered 1 (Monday) through 7 (Sunday).
type ISODate struct {
	Year int
	Week int
	Day  int
}

// FixedFromISO converts an ISO week-date (year, week, weekday 1..7 with Monday
// = 1) to an RD fixed day.
func FixedFromISO(year, week, day int) int {
	return NthKDay(week, Sunday, FixedFromGregorian(year-1, December, 28)) + day
}

// ISOFromFixed converts an RD fixed day to an ISO week date.
func ISOFromFixed(fixed int) ISODate {
	approx := GregorianYearFromFixed(fixed - 3)
	year := approx
	if fixed >= FixedFromISO(approx+1, 1, 1) {
		year = approx + 1
	}
	week := 1 + floorDiv(fixed-FixedFromISO(year, 1, 1), 7)
	day := amod(fixed, 7)
	return ISODate{Year: year, Week: week, Day: day}
}

// ISOLongYear reports whether the given ISO week-numbering year is a "long"
// year of 53 weeks (rather than 52).
func ISOLongYear(year int) bool {
	jan1 := WeekdayFromFixed(GregorianNewYear(year))
	dec31 := WeekdayFromFixed(GregorianYearEnd(year))
	return jan1 == Thursday || dec31 == Thursday
}

// ISOWeeksInYear returns the number of ISO weeks (52 or 53) in the given ISO
// week-numbering year.
func ISOWeeksInYear(year int) int {
	if ISOLongYear(year) {
		return 53
	}
	return 52
}

// ISOValid reports whether year, week and day form a valid ISO week date.
func ISOValid(year, week, day int) bool {
	if day < 1 || day > 7 || week < 1 {
		return false
	}
	return week <= ISOWeeksInYear(year)
}

// NewISO returns a validated ISODate, or an error if the fields are not a real
// ISO week date.
func NewISO(year, week, day int) (ISODate, error) {
	if !ISOValid(year, week, day) {
		return ISODate{}, fmt.Errorf("%w: iso %04d-W%02d-%d", ErrInvalidDate, year, week, day)
	}
	return ISODate{Year: year, Week: week, Day: day}, nil
}

// Fixed returns the RD fixed day of the ISO week date.
func (d ISODate) Fixed() int { return FixedFromISO(d.Year, d.Week, d.Day) }

// JDN returns the Julian Day Number of the ISO week date.
func (d ISODate) JDN() int { return JDNFromFixed(d.Fixed()) }

// Weekday returns the day of the week of the ISO week date.
func (d ISODate) Weekday() Weekday { return WeekdayFromFixed(d.Fixed()) }

// Valid reports whether the receiver is a valid ISO week date.
func (d ISODate) Valid() bool { return ISOValid(d.Year, d.Week, d.Day) }

// String formats the ISO week date as YYYY-Www-D.
func (d ISODate) String() string {
	return fmt.Sprintf("%04d-W%02d-%d", d.Year, d.Week, d.Day)
}

// ISOWeekOfYear returns the ISO week number of an RD fixed day.
func ISOWeekOfYear(fixed int) int { return ISOFromFixed(fixed).Week }

// ISOYearOf returns the ISO week-numbering year of an RD fixed day, which may
// differ from the Gregorian year near a year boundary.
func ISOYearOf(fixed int) int { return ISOFromFixed(fixed).Year }
