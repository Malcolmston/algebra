package calendrical

// GregorianEaster returns the RD fixed day of Western (Gregorian) Easter Sunday
// in the given Gregorian year, computed via the Gregorian computus.
func GregorianEaster(year int) int {
	century := floorDiv(year, 100) + 1
	shiftedEpact := floorMod(14+11*floorMod(year, 19)-floorDiv(3*century, 4)+floorDiv(5+8*century, 25), 30)
	adjustedEpact := shiftedEpact
	if shiftedEpact == 0 || (shiftedEpact == 1 && 10 < floorMod(year, 19)) {
		adjustedEpact = shiftedEpact + 1
	}
	paschalMoon := FixedFromGregorian(year, April, 19) - adjustedEpact
	return KDayAfter(Sunday, paschalMoon)
}

// JulianEaster returns the RD fixed day of Orthodox (Julian) Easter Sunday in
// the given year, computed via the Julian computus. The result is an RD fixed
// day; convert with GregorianFromFixed to obtain the civil (Gregorian) date.
func JulianEaster(year int) int {
	shiftedEpact := floorMod(14+11*floorMod(year, 19), 30)
	jYear := year
	if year <= 0 {
		jYear = year - 1
	}
	paschalMoon := FixedFromJulian(jYear, April, 19) - shiftedEpact
	return KDayAfter(Sunday, paschalMoon)
}

// GregorianEasterDate returns Western Easter Sunday of the given year as a
// Gregorian date.
func GregorianEasterDate(year int) GregorianDate {
	return GregorianFromFixed(GregorianEaster(year))
}

// OrthodoxEasterDate returns Orthodox Easter Sunday of the given year as the
// equivalent Gregorian (civil) date.
func OrthodoxEasterDate(year int) GregorianDate {
	return GregorianFromFixed(JulianEaster(year))
}

// AshWednesday returns the RD fixed day of Ash Wednesday (46 days before Western
// Easter) in the given Gregorian year.
func AshWednesday(year int) int {
	return GregorianEaster(year) - 46
}

// PalmSunday returns the RD fixed day of Palm Sunday (the Sunday before Western
// Easter) in the given Gregorian year.
func PalmSunday(year int) int {
	return GregorianEaster(year) - 7
}

// GoodFriday returns the RD fixed day of Good Friday (the Friday before Western
// Easter) in the given Gregorian year.
func GoodFriday(year int) int {
	return GregorianEaster(year) - 2
}

// Pentecost returns the RD fixed day of Pentecost (49 days after Western Easter)
// in the given Gregorian year.
func Pentecost(year int) int {
	return GregorianEaster(year) + 49
}

// Ascension returns the RD fixed day of Ascension Thursday (39 days after
// Western Easter) in the given Gregorian year.
func Ascension(year int) int {
	return GregorianEaster(year) + 39
}

// AdventSunday returns the RD fixed day of the first Sunday of Advent (the
// fourth Sunday before Christmas) in the given Gregorian year.
func AdventSunday(year int) int {
	return KDayNearest(Sunday, FixedFromGregorian(year, November, 30))
}

// Christmas returns the RD fixed day of Christmas Day (25 December) in the given
// Gregorian year.
func Christmas(year int) int {
	return FixedFromGregorian(year, December, 25)
}

// EpiphanyGregorian returns the RD fixed day of Epiphany (6 January) in the
// given Gregorian year.
func EpiphanyGregorian(year int) int {
	return FixedFromGregorian(year, January, 6)
}
