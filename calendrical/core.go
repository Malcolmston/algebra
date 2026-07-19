package calendrical

import (
	"errors"
	"math"
)

// Epoch constants expressed as Rata Die (RD) fixed-day counts. RD 1 is the
// proleptic Gregorian date 1 January 1 (a Monday).
const (
	// GregorianEpoch is the RD of Gregorian 1 January 1.
	GregorianEpoch = 1
	// JulianEpoch is the RD of Julian 1 January 1 (proleptic).
	JulianEpoch = -1
	// ISOEpoch is the RD of the ISO week-date origin (same as Gregorian).
	ISOEpoch = 1
	// HebrewEpoch is the RD of Hebrew 1 Tishri, year 1 (Anno Mundi).
	HebrewEpoch = -1373427
	// IslamicEpoch is the RD of Islamic 1 Muharram, year 1.
	IslamicEpoch = 227015
	// PersianEpoch is the RD of Persian 1 Farvardin, year 1.
	PersianEpoch = 226896
	// CopticEpoch is the RD of Coptic 1 Thoout, year 1.
	CopticEpoch = 103605
	// EthiopicEpoch is the RD of Ethiopic 1 Maskaram, year 1.
	EthiopicEpoch = 2796
	// MJDEpoch is the RD of Modified Julian Date 0 (17 November 1858).
	MJDEpoch = 678576
	// UnixEpoch is the RD of the Unix epoch (1 January 1970).
	UnixEpoch = 719163
)

// JDEpoch is the RD moment corresponding to Julian Date 0.0. Because Julian
// days begin at noon, this value has a fractional part of one half.
const JDEpoch = -1721424.5

// Errors returned by validating constructors and parsers.
var (
	// ErrRange indicates a field is outside the valid range for its calendar.
	ErrRange = errors.New("calendrical: value out of range")
	// ErrInvalidDate indicates a syntactically or semantically invalid date.
	ErrInvalidDate = errors.New("calendrical: invalid date")
)

// floorDiv returns the floor of a/b, i.e. the quotient rounded toward negative
// infinity. b must be non-zero. Go's built-in division truncates toward zero,
// which differs for mixed-sign operands.
func floorDiv(a, b int) int {
	q := a / b
	if (a%b != 0) && ((a < 0) != (b < 0)) {
		q--
	}
	return q
}

// floorMod returns the non-negative remainder of a modulo b when b is positive
// (more generally, a value with the same sign as b). It satisfies
// a == b*floorDiv(a,b) + floorMod(a,b).
func floorMod(a, b int) int {
	m := a % b
	if m != 0 && ((m < 0) != (b < 0)) {
		m += b
	}
	return m
}

// amod returns the adjusted modulus of x by y: the unique value in the range
// [1, y] (for positive y) congruent to x modulo y. It is used by calendars
// whose components are numbered from 1 rather than 0.
func amod(x, y int) int {
	r := floorMod(x, y)
	if r == 0 {
		return y
	}
	return r
}

// FloorDiv returns the floor of a divided by b (rounded toward negative
// infinity). It panics if b is zero.
func FloorDiv(a, b int) int {
	if b == 0 {
		panic("calendrical: division by zero")
	}
	return floorDiv(a, b)
}

// FloorMod returns the mathematical modulus of a by b, carrying the sign of b.
// It panics if b is zero.
func FloorMod(a, b int) int {
	if b == 0 {
		panic("calendrical: division by zero")
	}
	return floorMod(a, b)
}

// AdjustedMod returns the value in [1, y] congruent to x modulo y, for y > 0.
// It panics if y is zero.
func AdjustedMod(x, y int) int {
	if y == 0 {
		panic("calendrical: division by zero")
	}
	return amod(x, y)
}

// --- Julian Day Number and related day counts ------------------------------

// JDNFromFixed converts an RD fixed-day count to a Julian Day Number, the
// integer count of days at noon since the start of the Julian Period.
func JDNFromFixed(fixed int) int { return fixed + 1721425 }

// FixedFromJDN converts a Julian Day Number to an RD fixed-day count.
func FixedFromJDN(jdn int) int { return jdn - 1721425 }

// MJDFromFixed converts an RD fixed-day count to a Modified Julian Date.
func MJDFromFixed(fixed int) int { return fixed - MJDEpoch }

// FixedFromMJD converts a Modified Julian Date to an RD fixed-day count.
func FixedFromMJD(mjd int) int { return mjd + MJDEpoch }

// UnixDayFromFixed converts an RD fixed-day count to the number of whole days
// since the Unix epoch (1 January 1970).
func UnixDayFromFixed(fixed int) int { return fixed - UnixEpoch }

// FixedFromUnixDay converts a day count since the Unix epoch to an RD fixed-day
// count.
func FixedFromUnixDay(day int) int { return day + UnixEpoch }

// UnixSecondsFromFixed converts an RD fixed-day count (interpreted at midnight
// UTC) to the corresponding number of Unix seconds.
func UnixSecondsFromFixed(fixed int) int64 {
	return int64(fixed-UnixEpoch) * 86400
}

// FixedFromUnixSeconds converts a count of Unix seconds to the RD fixed-day
// count of the day containing that instant.
func FixedFromUnixSeconds(sec int64) int {
	return int(floorDivInt64(sec, 86400)) + UnixEpoch
}

func floorDivInt64(a, b int64) int64 {
	q := a / b
	if (a%b != 0) && ((a < 0) != (b < 0)) {
		q--
	}
	return q
}

// MomentFromJD converts a (possibly fractional) Julian Date to an RD moment,
// a real number of days since RD 0 at midnight.
func MomentFromJD(jd float64) float64 { return jd + JDEpoch }

// JDFromMoment converts an RD moment to a Julian Date.
func JDFromMoment(t float64) float64 { return t - JDEpoch }

// MJDFromJD converts a Julian Date to a Modified Julian Date.
func MJDFromJD(jd float64) float64 { return jd - 2400000.5 }

// JDFromMJD converts a Modified Julian Date to a Julian Date.
func JDFromMJD(mjd float64) float64 { return mjd + 2400000.5 }

// FixedFromMoment returns the RD fixed day containing an RD moment (its floor).
func FixedFromMoment(t float64) int { return int(math.Floor(t)) }

// TimeFromMoment returns the fractional time-of-day component of an RD moment,
// in the half-open interval [0, 1).
func TimeFromMoment(t float64) float64 { return t - math.Floor(t) }
