// Package calendrical implements calendar and date arithmetic using only the
// Go standard library. It is a pure-arithmetic library: the core conversions
// do not depend on the time package and are deterministic — identical inputs
// always yield identical outputs.
//
// The package is organised around the Rata Die (RD) fixed-day count, an
// integer number of days elapsed since the proleptic Gregorian date
// 1 January 1 (RD 1, a Monday). Every calendar is converted to and from this
// common representation, following the algorithms described by Nachum
// Dershowitz and Edward M. Reingold in Calendrical Calculations. Julian Day
// Numbers (JDN), Modified Julian Dates (MJD) and Unix day counts are exposed
// as thin conversions on top of the fixed-day count.
//
// Supported calendars:
//
//   - Proleptic Gregorian calendar.
//   - Proleptic Julian calendar (no year zero).
//   - ISO-8601 week date (year, week, weekday).
//   - Hebrew (Jewish) calendar, including molad and holiday helpers.
//   - Islamic calendar in its tabular/arithmetical form.
//   - Persian (Jalali) calendar in its arithmetical (2820-year) form.
//   - Coptic and Ethiopic calendars.
//
// In addition the package provides day-of-week arithmetic (k-day search),
// leap-year predicates, date differences and offsets, Gregorian (Western) and
// Julian (Orthodox) Easter via the computus, and low-precision astronomical
// approximations for the equinoxes, solstices and the phases of the Moon based
// on Jean Meeus's Astronomical Algorithms.
//
// Unless a function name says otherwise, "fixed" arguments and results are RD
// day counts (int), month and day numbers are 1-based, and Gregorian years use
// astronomical numbering (there is a year 0). Astronomical routines return
// Julian Ephemeris Days (JDE, float64) in Terrestrial Time.
package calendrical
