package calendrical

// AddDays returns the RD fixed day n days after the given fixed day (n may be
// negative).
func AddDays(fixed, n int) int { return fixed + n }

// SubtractDays returns the RD fixed day n days before the given fixed day.
func SubtractDays(fixed, n int) int { return fixed - n }

// DayDifference returns the signed number of days from fixed1 to fixed2.
func DayDifference(fixed1, fixed2 int) int { return fixed2 - fixed1 }

// WeekDifference returns the signed number of whole weeks from fixed1 to fixed2
// (truncated toward zero).
func WeekDifference(fixed1, fixed2 int) int { return (fixed2 - fixed1) / 7 }

// PlusDays returns a new GregorianDate n days after the receiver.
func (d GregorianDate) PlusDays(n int) GregorianDate {
	return GregorianFromFixed(d.Fixed() + n)
}

// PlusMonths returns the GregorianDate n calendar months after the receiver.
// If the resulting month is shorter than the receiver's day, the day is clamped
// to the last valid day of that month (so 31 January plus one month is
// 28 February, or 29 February in a leap year).
func (d GregorianDate) PlusMonths(n int) GregorianDate {
	totalMonths := (d.Year*12 + (d.Month - 1)) + n
	year := floorDiv(totalMonths, 12)
	month := floorMod(totalMonths, 12) + 1
	day := d.Day
	if last := GregorianDaysInMonth(year, month); day > last {
		day = last
	}
	return GregorianDate{Year: year, Month: month, Day: day}
}

// PlusYears returns the GregorianDate n years after the receiver, clamping
// 29 February to 28 February in common years.
func (d GregorianDate) PlusYears(n int) GregorianDate {
	year := d.Year + n
	day := d.Day
	if last := GregorianDaysInMonth(year, d.Month); day > last {
		day = last
	}
	return GregorianDate{Year: year, Month: d.Month, Day: day}
}

// Before reports whether the receiver precedes the other Gregorian date.
func (d GregorianDate) Before(o GregorianDate) bool { return d.Fixed() < o.Fixed() }

// After reports whether the receiver follows the other Gregorian date.
func (d GregorianDate) After(o GregorianDate) bool { return d.Fixed() > o.Fixed() }

// Equal reports whether the receiver denotes the same day as the other date.
func (d GregorianDate) Equal(o GregorianDate) bool { return d.Fixed() == o.Fixed() }

// Sub returns the signed number of days from the other date to the receiver.
func (d GregorianDate) Sub(o GregorianDate) int { return d.Fixed() - o.Fixed() }

// GregorianYMD returns the year, month and day components of an RD fixed day as
// a convenience over GregorianFromFixed.
func GregorianYMD(fixed int) (year, month, day int) {
	g := GregorianFromFixed(fixed)
	return g.Year, g.Month, g.Day
}

// YearMonthDayDifference returns the difference between two RD fixed days broken
// into whole Gregorian years, months and days, with fixed1 <= fixed2 assumed
// for a non-negative result (the sign is preserved when reversed). It counts
// complete calendar periods, matching the usual notion of "n years, m months
// and d days between two dates".
func YearMonthDayDifference(fixed1, fixed2 int) (years, months, days int) {
	if fixed1 > fixed2 {
		y, m, d := YearMonthDayDifference(fixed2, fixed1)
		return -y, -m, -d
	}
	a := GregorianFromFixed(fixed1)
	b := GregorianFromFixed(fixed2)
	years = b.Year - a.Year
	months = b.Month - a.Month
	days = b.Day - a.Day
	if days < 0 {
		months--
		// borrow days from the month preceding b's month
		pm := b.Month - 1
		py := b.Year
		if pm < 1 {
			pm = 12
			py--
		}
		days += GregorianDaysInMonth(py, pm)
	}
	if months < 0 {
		years--
		months += 12
	}
	return years, months, days
}

// GregorianWeekday returns the day of the week of a Gregorian date.
func GregorianWeekday(year, month, day int) Weekday {
	return WeekdayFromFixed(FixedFromGregorian(year, month, day))
}

// NthWeekdayOfMonth returns the RD fixed day of the n-th given weekday in the
// specified Gregorian month. When n is positive it counts from the start of the
// month (n = 1 is the first such weekday); when n is negative it counts from the
// end (n = -1 is the last such weekday). It returns 0 and false if there is no
// such day (for example a fifth Monday that does not exist).
func NthWeekdayOfMonth(n int, k Weekday, year, month int) (int, bool) {
	if n == 0 {
		return 0, false
	}
	var result int
	if n > 0 {
		first := FixedFromGregorian(year, month, 1)
		result = KDayOnOrAfter(k, first) + 7*(n-1)
	} else {
		last := FixedFromGregorian(year, month, GregorianDaysInMonth(year, month))
		result = KDayOnOrBefore(k, last) + 7*(n+1)
	}
	if GregorianYearFromFixed(result) != year || GregorianFromFixed(result).Month != month {
		return 0, false
	}
	return result, true
}

// FirstWeekdayOfMonth returns the RD fixed day of the first given weekday in the
// specified Gregorian month.
func FirstWeekdayOfMonth(k Weekday, year, month int) int {
	return KDayOnOrAfter(k, FixedFromGregorian(year, month, 1))
}

// LastWeekdayOfMonth returns the RD fixed day of the last given weekday in the
// specified Gregorian month.
func LastWeekdayOfMonth(k Weekday, year, month int) int {
	return KDayOnOrBefore(k, FixedFromGregorian(year, month, GregorianDaysInMonth(year, month)))
}

// USMemorialDay returns the RD fixed day of United States Memorial Day (the last
// Monday in May) for the given Gregorian year.
func USMemorialDay(year int) int {
	return LastWeekdayOfMonth(Monday, year, May)
}

// USLaborDay returns the RD fixed day of United States Labor Day (the first
// Monday in September) for the given Gregorian year.
func USLaborDay(year int) int {
	return FirstWeekdayOfMonth(Monday, year, September)
}

// USThanksgiving returns the RD fixed day of United States Thanksgiving (the
// fourth Thursday in November) for the given Gregorian year.
func USThanksgiving(year int) int {
	d, _ := NthWeekdayOfMonth(4, Thursday, year, November)
	return d
}

// USElectionDay returns the RD fixed day of United States federal Election Day
// (the Tuesday after the first Monday in November) for the given Gregorian year.
func USElectionDay(year int) int {
	return KDayAfter(Tuesday, FirstWeekdayOfMonth(Monday, year, November))
}

// USIndependenceDay returns the RD fixed day of United States Independence Day
// (4 July) for the given Gregorian year.
func USIndependenceDay(year int) int {
	return FixedFromGregorian(year, July, 4)
}

// USThanksgivingCanada returns the RD fixed day of Canadian Thanksgiving (the
// second Monday in October) for the given Gregorian year.
func USThanksgivingCanada(year int) int {
	d, _ := NthWeekdayOfMonth(2, Monday, year, October)
	return d
}

// DayOfYearGregorian returns the ordinal day number within the year (1 = 1
// January) of an RD fixed day.
func DayOfYearGregorian(fixed int) int {
	year := GregorianYearFromFixed(fixed)
	return fixed - GregorianYearEnd(year-1)
}

// WeekOfYearGregorian returns the week number of an RD fixed day counting weeks
// from 1 January, where the first week is week 1 and weeks are Sunday-based.
func WeekOfYearGregorian(fixed int) int {
	year := GregorianYearFromFixed(fixed)
	jan1 := GregorianNewYear(year)
	return floorDiv(fixed-jan1+DayOfWeekFromFixed(jan1), 7) + 1
}

// IsWeekend reports whether an RD fixed day falls on a Saturday or Sunday.
func IsWeekend(fixed int) bool {
	w := WeekdayFromFixed(fixed)
	return w == Saturday || w == Sunday
}

// IsWeekday reports whether an RD fixed day falls on a Monday through Friday.
func IsWeekday(fixed int) bool { return !IsWeekend(fixed) }

// BusinessDaysBetween returns the number of Monday-to-Friday days in the
// half-open interval [start, end). The result is negative if end precedes
// start.
func BusinessDaysBetween(start, end int) int {
	if end < start {
		return -BusinessDaysBetween(end, start)
	}
	count := 0
	for d := start; d < end; d++ {
		if IsWeekday(d) {
			count++
		}
	}
	return count
}

// AddBusinessDays returns the RD fixed day reached by advancing n business days
// (Monday–Friday) from the given fixed day. n may be negative to move backward.
// If the starting day is a weekend, counting begins from the next business day
// in the direction of travel.
func AddBusinessDays(fixed, n int) int {
	step := 1
	if n < 0 {
		step = -1
		n = -n
	}
	d := fixed
	for n > 0 {
		d += step
		if IsWeekday(d) {
			n--
		}
	}
	return d
}
