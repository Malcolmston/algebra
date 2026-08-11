# calendrical — implements calendar and date arithmetic using only the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/calendrical.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/calendrical)

Package calendrical implements calendar and date arithmetic using only the Go
standard library. It is a pure-arithmetic library: the core conversions do not
depend on the time package and are deterministic — identical inputs always
yield identical outputs.

The package is organised around the Rata Die (RD) fixed-day count, an integer
number of days elapsed since the proleptic Gregorian date 1 January 1 (RD 1, a
Monday). Every calendar is converted to and from this common representation,
following the algorithms described by Nachum Dershowitz and Edward M. Reingold
in Calendrical Calculations. Julian Day Numbers (JDN), Modified Julian Dates
(MJD) and Unix day counts are exposed as thin conversions on top of the
fixed-day count.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/calendrical"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AddBusinessDays(fixed, n int) int` | AddBusinessDays returns the RD fixed day reached by advancing n business days (Monday–Friday) from the given fixed day. |
| `func AddDays(fixed, n int) int` | AddDays returns the RD fixed day n days after the given fixed day (n may be negative). |
| `func AdjustedMod(x, y int) int` | AdjustedMod returns the value in [1, y] congruent to x modulo y, for y > 0. |
| `func AdventSunday(year int) int` | AdventSunday returns the RD fixed day of the first Sunday of Advent (the fourth Sunday before Christmas) in the given Gregorian year. |
| `func Ascension(year int) int` | Ascension returns the RD fixed day of Ascension Thursday (39 days after Western Easter) in the given Gregorian year. |
| `func AshWednesday(year int) int` | AshWednesday returns the RD fixed day of Ash Wednesday (46 days before Western Easter) in the given Gregorian year. |
| `func AstronomicalToJulianYear(year int) int` | AstronomicalToJulianYear converts a continuous astronomical year number to a Julian-calendar year that omits zero. |
| `func AutumnEquinoxJD(year int) float64` | AutumnEquinoxJD returns the Julian Ephemeris Day of the September equinox of the given Gregorian year. |
| `func BusinessDaysBetween(start, end int) int` | BusinessDaysBetween returns the number of Monday-to-Friday days in the half-open interval [start, end). |
| `func Christmas(year int) int` | Christmas returns the RD fixed day of Christmas Day (25 December) in the given Gregorian year. |
| `func CopticChristmas(gYear int) []int` | CopticChristmas returns every RD fixed day in the given Gregorian year on which Coptic Christmas (29 Koiak) falls. |
| `func CopticDaysInMonth(year, month int) int` | CopticDaysInMonth returns the number of days in the given month of the given Coptic year (30 for months 1–12; 5 or 6 for the epagomenae). |
| `func CopticDaysInYear(year int) int` | CopticDaysInYear returns the number of days in the given Coptic year (365 or 366). |
| `func CopticInGregorian(cMonth, cDay, gYear int) []int` | CopticInGregorian returns every RD fixed day in the given Gregorian year on which the given Coptic month and day occur. |
| `func CopticLeapYear(year int) bool` | CopticLeapYear reports whether the given Coptic year is a leap year (every fourth year, with the leap day added to the thirteenth month). |
| `func CopticMonthName(month int) string` | CopticMonthName returns the transliterated name of a Coptic month number. |
| `func CopticValid(year, month, day int) bool` | CopticValid reports whether year, month and day form a valid Coptic date. |
| `func DayDifference(fixed1, fixed2 int) int` | DayDifference returns the signed number of days from fixed1 to fixed2. |
| `func DayOfWeekFromFixed(fixed int) int` | DayOfWeekFromFixed returns the day-of-week index (0 = Sunday) of an RD fixed day. |
| `func DayOfYearGregorian(fixed int) int` | DayOfYearGregorian returns the ordinal day number within the year (1 = 1 January) of an RD fixed day. |
| `func DaysInHebrewYear(year int) int` | DaysInHebrewYear returns the length in days of the given Hebrew year (353, 354, 355, 383, 384 or 385). |
| `func EpiphanyGregorian(year int) int` | EpiphanyGregorian returns the RD fixed day of Epiphany (6 January) in the given Gregorian year. |
| `func EthiopicDaysInMonth(year, month int) int` | EthiopicDaysInMonth returns the number of days in the given month of the given Ethiopic year. |
| `func EthiopicDaysInYear(year int) int` | EthiopicDaysInYear returns the number of days in the given Ethiopic year. |
| `func EthiopicInGregorian(eMonth, eDay, gYear int) []int` | EthiopicInGregorian returns every RD fixed day in the given Gregorian year on which the given Ethiopic month and day occur. |
| `func EthiopicLeapYear(year int) bool` | EthiopicLeapYear reports whether the given Ethiopic year is a leap year. |
| `func EthiopicMonthName(month int) string` | EthiopicMonthName returns the transliterated name of an Ethiopic month number. |
| `func EthiopicValid(year, month, day int) bool` | EthiopicValid reports whether year, month and day form a valid Ethiopic date. |
| `func FirstKDay(k Weekday, fixed int) int` | FirstKDay returns the RD of the first weekday k on or after the fixed anchor day. |
| `func FirstWeekdayOfMonth(k Weekday, year, month int) int` | FirstWeekdayOfMonth returns the RD fixed day of the first given weekday in the specified Gregorian month. |
| `func FixedFromCoptic(year, month, day int) int` | FixedFromCoptic converts a Coptic year, month and day to an RD fixed day. |
| `func FixedFromEthiopic(year, month, day int) int` | FixedFromEthiopic converts an Ethiopic year, month and day to an RD fixed day. |
| `func FixedFromGregorian(year, month, day int) int` | FixedFromGregorian converts a Gregorian year, month and day to an RD fixed day. |
| `func FixedFromHebrew(year, month, day int) int` | FixedFromHebrew converts a Hebrew year, month and day to an RD fixed day. |
| `func FixedFromISO(year, week, day int) int` | FixedFromISO converts an ISO week-date (year, week, weekday 1..7 with Monday = 1) to an RD fixed day. |
| `func FixedFromIslamic(year, month, day int) int` | FixedFromIslamic converts an Islamic year, month and day to an RD fixed day. |
| `func FixedFromJDN(jdn int) int` | FixedFromJDN converts a Julian Day Number to an RD fixed-day count. |
| `func FixedFromJulian(year, month, day int) int` | FixedFromJulian converts a Julian-calendar year, month and day to an RD fixed day. |
| `func FixedFromMJD(mjd int) int` | FixedFromMJD converts a Modified Julian Date to an RD fixed-day count. |
| `func FixedFromMoment(t float64) int` | FixedFromMoment returns the RD fixed day containing an RD moment (its floor). |
| `func FixedFromPersian(year, month, day int) int` | FixedFromPersian converts a Persian (Jalali) year, month and day to an RD fixed day. |
| `func FixedFromUnixDay(day int) int` | FixedFromUnixDay converts a day count since the Unix epoch to an RD fixed-day count. |
| `func FixedFromUnixSeconds(sec int64) int` | FixedFromUnixSeconds converts a count of Unix seconds to the RD fixed-day count of the day containing that instant. |
| `func FloorDiv(a, b int) int` | FloorDiv returns the floor of a divided by b (rounded toward negative infinity). |
| `func FloorMod(a, b int) int` | FloorMod returns the mathematical modulus of a by b, carrying the sign of b. |
| `func GoodFriday(year int) int` | GoodFriday returns the RD fixed day of Good Friday (the Friday before Western Easter) in the given Gregorian year. |
| `func GregorianDateDifference(y1, m1, d1, y2, m2, d2 int) int` | GregorianDateDifference returns the number of days from the first Gregorian date to the second (positive when the second date is later). |
| `func GregorianDayNumber(year, month, day int) int` | GregorianDayNumber returns the ordinal day-of-year (1 = 1 January) of the given Gregorian date. |
| `func GregorianDaysInMonth(year, month int) int` | GregorianDaysInMonth returns the number of days in the given Gregorian month of the given year. |
| `func GregorianDaysInYear(year int) int` | GregorianDaysInYear returns the number of days in the given Gregorian year (365 or 366). |
| `func GregorianDaysRemaining(year, month, day int) int` | GregorianDaysRemaining returns the number of days remaining in the year after the given Gregorian date. |
| `func GregorianEaster(year int) int` | GregorianEaster returns the RD fixed day of Western (Gregorian) Easter Sunday in the given Gregorian year, computed via the Gregorian computus. |
| `func GregorianLeapYear(year int) bool` | GregorianLeapYear reports whether the given Gregorian year is a leap year. |
| `func GregorianMonthName(month int) string` | GregorianMonthName returns the English name of a Gregorian month (1 = January). |
| `func GregorianNewYear(year int) int` | GregorianNewYear returns the RD fixed day of 1 January of the given year. |
| `func GregorianValid(year, month, day int) bool` | GregorianValid reports whether year, month and day form a valid Gregorian date. |
| `func GregorianYMD(fixed int) (year, month, day int)` | GregorianYMD returns the year, month and day components of an RD fixed day as a convenience over GregorianFromFixed. |
| `func GregorianYearEnd(year int) int` | GregorianYearEnd returns the RD fixed day of 31 December of the given year. |
| `func GregorianYearFromFixed(fixed int) int` | GregorianYearFromFixed returns the Gregorian year containing an RD fixed day. |
| `func GregorianYearRange(year int) (start, end int)` | GregorianYearRange returns the half-open RD range [start, end) spanned by the given Gregorian year, that is 1 January of the year up to but excluding… |
| `func HanukkahStart(gYear int) int` | HanukkahStart returns the RD fixed day of the first day of Hanukkah (25 Kislev) occurring in the given Gregorian year. |
| `func HebrewBirthday(birth HebrewDate, hYear int) int` | HebrewBirthday returns the RD fixed day on which the anniversary of the given Hebrew birth date falls in the given Hebrew year. |
| `func HebrewBirthdayInGregorian(birth HebrewDate, gYear int) []int` | HebrewBirthdayInGregorian returns every RD fixed day in the given Gregorian year on which the anniversary of a Hebrew birth date falls. |
| `func HebrewCalendarElapsedDays(year int) int` | HebrewCalendarElapsedDays returns the number of days elapsed from the Hebrew epoch to the mean conjunction (molad) of Tishri of the given year,… |
| `func HebrewDaysInMonth(year, month int) int` | HebrewDaysInMonth is an alias for LastDayOfHebrewMonth. |
| `func HebrewLeapYear(year int) bool` | HebrewLeapYear reports whether the given Hebrew year is embolismic (has a thirteenth month), under the 19-year Metonic cycle. |
| `func HebrewMolad(year, month int) float64` | HebrewMolad returns the RD moment (a real number of days since RD 0) of the mean lunar conjunction — the molad — that begins the given Hebrew… |
| `func HebrewMonthName(year, month int) string` | HebrewMonthName returns the English transliteration of a Hebrew month number. |
| `func HebrewMonthsElapsed(year, month int) int` | HebrewMonthsElapsed returns the number of lunar months from the Hebrew epoch to the start (molad) of the given month of the given Hebrew year. |
| `func HebrewNewYear(year int) int` | HebrewNewYear returns the RD fixed day of Rosh ha-Shanah (1 Tishri) of the given Hebrew year. |
| `func HebrewNewYearDelay(year int) int` | HebrewNewYearDelay returns the year-length correction (0, 1 or 2 days) applied to Rosh ha-Shanah so that the year has a permissible length. |
| `func HebrewValid(year, month, day int) bool` | HebrewValid reports whether year, month and day form a valid Hebrew date. |
| `func HebrewYearFromGregorian(gYear int) int` | HebrewYearFromGregorian returns the Hebrew year whose autumn new year (Rosh ha-Shanah) falls in the given Gregorian year. |
| `func ISOLongYear(year int) bool` | ISOLongYear reports whether the given ISO week-numbering year is a "long" year of 53 weeks (rather than 52). |
| `func ISOValid(year, week, day int) bool` | ISOValid reports whether year, week and day form a valid ISO week date. |
| `func ISOWeekOfYear(fixed int) int` | ISOWeekOfYear returns the ISO week number of an RD fixed day. |
| `func ISOWeeksInYear(year int) int` | ISOWeeksInYear returns the number of ISO weeks (52 or 53) in the given ISO week-numbering year. |
| `func ISOYearOf(fixed int) int` | ISOYearOf returns the ISO week-numbering year of an RD fixed day, which may differ from the Gregorian year near a year boundary. |
| `func IsWeekday(fixed int) bool` | IsWeekday reports whether an RD fixed day falls on a Monday through Friday. |
| `func IsWeekend(fixed int) bool` | IsWeekend reports whether an RD fixed day falls on a Saturday or Sunday. |
| `func IslamicDaysInMonth(year, month int) int` | IslamicDaysInMonth returns the number of days in the given month of the given Islamic year: odd months have 30 days, even months 29, and the twelfth… |
| `func IslamicDaysInYear(year int) int` | IslamicDaysInYear returns the number of days in the given Islamic year (354 or 355). |
| `func IslamicInGregorian(iMonth, iDay, gYear int) []int` | IslamicInGregorian returns every RD fixed day in the given Gregorian year on which the given Islamic month and day occur. |
| `func IslamicLeapYear(year int) bool` | IslamicLeapYear reports whether the given Islamic year is a leap year in the tabular calendar (11 leap years in each 30-year cycle). |
| `func IslamicMonthName(month int) string` | IslamicMonthName returns the English transliteration of an Islamic month number. |
| `func IslamicNewYear(year int) int` | IslamicNewYear returns the RD fixed day of 1 Muharram of the given Islamic year. |
| `func IslamicValid(year, month, day int) bool` | IslamicValid reports whether year, month and day form a valid tabular Islamic date. |
| `func JDFromMJD(mjd float64) float64` | JDFromMJD converts a Modified Julian Date to a Julian Date. |
| `func JDFromMoment(t float64) float64` | JDFromMoment converts an RD moment to a Julian Date. |
| `func JDNFromFixed(fixed int) int` | JDNFromFixed converts an RD fixed-day count to a Julian Day Number, the integer count of days at noon since the start of the Julian Period. |
| `func JulianDaysInMonth(year, month int) int` | JulianDaysInMonth returns the number of days in the given month of the given Julian-calendar year. |
| `func JulianDaysInYear(year int) int` | JulianDaysInYear returns the number of days in the given Julian-calendar year (365 or 366). |
| `func JulianEaster(year int) int` | JulianEaster returns the RD fixed day of Orthodox (Julian) Easter Sunday in the given year, computed via the Julian computus. |
| `func JulianInGregorian(jMonth, jDay, gYear int) []int` | JulianInGregorian returns every RD fixed day in the given Gregorian year on which the given Julian-calendar month and day occur. |
| `func JulianLeapYear(year int) bool` | JulianLeapYear reports whether the given Julian-calendar year is a leap year. |
| `func JulianValid(year, month, day int) bool` | JulianValid reports whether year, month and day form a valid Julian-calendar date. |
| `func JulianYearToAstronomical(year int) int` | JulianYearToAstronomical converts a Julian-calendar year that omits zero (…, -1 = 1 BCE, 1 = 1 CE, …) to a continuous astronomical year number. |
| `func KDayAfter(k Weekday, fixed int) int` | KDayAfter returns the RD of the first day that is weekday k strictly after the given fixed day. |
| `func KDayBefore(k Weekday, fixed int) int` | KDayBefore returns the RD of the last day that is weekday k strictly before the given fixed day. |
| `func KDayNearest(k Weekday, fixed int) int` | KDayNearest returns the RD of the day that is weekday k nearest to the given fixed day (ties resolve to the later day). |
| `func KDayOnOrAfter(k Weekday, fixed int) int` | KDayOnOrAfter returns the RD of the first day that is weekday k on or after the given fixed day. |
| `func KDayOnOrBefore(k Weekday, fixed int) int` | KDayOnOrBefore returns the RD of the last day that is weekday k on or before the given fixed day. |
| `func LastDayOfHebrewMonth(year, month int) int` | LastDayOfHebrewMonth returns the number of days (29 or 30) in the given month of the given Hebrew year. |
| `func LastKDay(k Weekday, fixed int) int` | LastKDay returns the RD of the last weekday k on or before the fixed anchor day. |
| `func LastMonthOfHebrewYear(year int) int` | LastMonthOfHebrewYear returns the number of the final month of the Hebrew year: 12 in a common year, 13 in a leap year. |
| `func LastWeekdayOfMonth(k Weekday, year, month int) int` | LastWeekdayOfMonth returns the RD fixed day of the last given weekday in the specified Gregorian month. |
| `func LongMarheshvan(year int) bool` | LongMarheshvan reports whether Marheshvan has 30 days in the given Hebrew year (a "complete" year). |
| `func LunationNumber(jd float64) int` | LunationNumber returns the integer lunation index k whose mean New Moon is nearest to the given Julian Ephemeris Day. |
| `func MJDFromFixed(fixed int) int` | MJDFromFixed converts an RD fixed-day count to a Modified Julian Date. |
| `func MJDFromJD(jd float64) float64` | MJDFromJD converts a Julian Date to a Modified Julian Date. |
| `func MeanNewMoonJD(k float64) float64` | MeanNewMoonJD returns the Julian Ephemeris Day of the mean New Moon for the given lunation index k, where k = 0 corresponds to the New Moon of 6… |
| `func MeanSeasonJD(year int, s Season) float64` | MeanSeasonJD returns the Julian Ephemeris Day of the mean (uncorrected) equinox or solstice for the given Gregorian year, using Meeus's polynomial… |
| `func MomentFromJD(jd float64) float64` | MomentFromJD converts a (possibly fractional) Julian Date to an RD moment, a real number of days since RD 0 at midnight. |
| `func MoonPhaseFixed(k int, phase MoonPhase) int` | MoonPhaseFixed returns the RD fixed day (in Terrestrial Time) of the requested principal phase for the lunation nearest index k. |
| `func MoonPhaseJD(k int, phase MoonPhase) float64` | MoonPhaseJD returns the Julian Ephemeris Day of the requested principal phase for the lunation nearest index k. |
| `func NextMoonPhaseJD(jd float64, phase MoonPhase) float64` | NextMoonPhaseJD returns the Julian Ephemeris Day of the first occurrence of the requested principal phase on or after the given Julian Ephemeris Day. |
| `func Nowruz(gYear int) int` | Nowruz returns the RD fixed day of Nowruz (the Persian new year, 1 Farvardin) occurring in the given Gregorian year. |
| `func NthKDay(n int, k Weekday, fixed int) int` | NthKDay returns the RD of the n-th weekday k relative to the given fixed anchor day. |
| `func NthWeekdayOfMonth(n int, k Weekday, year, month int) (int, bool)` | NthWeekdayOfMonth returns the RD fixed day of the n-th given weekday in the specified Gregorian month. |
| `func OmerCount(fixed int) (week, day int, ok bool)` | OmerCount returns the position of an RD fixed day in the counting of the Omer, the 49 days from the second night of Passover to the eve of Shavuot,… |
| `func PalmSunday(year int) int` | PalmSunday returns the RD fixed day of Palm Sunday (the Sunday before Western Easter) in the given Gregorian year. |
| `func Passover(gYear int) int` | Passover returns the RD fixed day of the first day of Passover (15 Nisan) occurring in the given Gregorian year. |
| `func Pentecost(year int) int` | Pentecost returns the RD fixed day of Pentecost (49 days after Western Easter) in the given Gregorian year. |
| `func PersianDaysInMonth(year, month int) int` | PersianDaysInMonth returns the number of days in the given month of the given Persian year: 31 for the first six months, 30 for months 7–11, and 29… |
| `func PersianDaysInYear(year int) int` | PersianDaysInYear returns the number of days in the given Persian year (365 or 366). |
| `func PersianLeapYear(year int) bool` | PersianLeapYear reports whether the given Persian year is a leap year in the arithmetical (2820-year) calendar. |
| `func PersianMonthName(month int) string` | PersianMonthName returns the transliterated name of a Persian month number. |
| `func PersianValid(year, month, day int) bool` | PersianValid reports whether year, month and day form a valid Persian date. |
| `func PersianYearFromFixed(fixed int) int` | PersianYearFromFixed returns the Persian year containing an RD fixed day. |
| `func PreviousMoonPhaseJD(jd float64, phase MoonPhase) float64` | PreviousMoonPhaseJD returns the Julian Ephemeris Day of the last occurrence of the requested principal phase on or before the given Julian Ephemeris… |
| `func Purim(gYear int) int` | Purim returns the RD fixed day of Purim occurring in the given Gregorian year. |
| `func RamadanStart(year int) int` | RamadanStart returns the RD fixed day of 1 Ramadan of the given Islamic year in the tabular calendar. |
| `func RoshHaShanah(gYear int) int` | RoshHaShanah returns the RD fixed day of Rosh ha-Shanah (1 Tishri) occurring in the given Gregorian year. |
| `func SeasonFixed(year int, s Season) int` | SeasonFixed returns the RD fixed day (in Terrestrial Time, at Greenwich) of the apparent equinox or solstice for the given Gregorian year. |
| `func SeasonJD(year int, s Season) float64` | SeasonJD returns the Julian Ephemeris Day of the apparent equinox or solstice for the given Gregorian year, applying Meeus's periodic corrections to… |
| `func Shavuot(gYear int) int` | Shavuot returns the RD fixed day of Shavuot (6 Sivan) occurring in the given Gregorian year. |
| `func ShortKislev(year int) bool` | ShortKislev reports whether Kislev has 29 days in the given Hebrew year (a "deficient" year). |
| `func SpringEquinoxJD(year int) float64` | SpringEquinoxJD returns the Julian Ephemeris Day of the March equinox of the given Gregorian year. |
| `func SubtractDays(fixed, n int) int` | SubtractDays returns the RD fixed day n days before the given fixed day. |
| `func Sukkot(gYear int) int` | Sukkot returns the RD fixed day of the first day of Sukkot (15 Tishri) occurring in the given Gregorian year. |
| `func SummerSolsticeJD(year int) float64` | SummerSolsticeJD returns the Julian Ephemeris Day of the June solstice of the given Gregorian year. |
| `func TimeFromMoment(t float64) float64` | TimeFromMoment returns the fractional time-of-day component of an RD moment, in the half-open interval [0, 1). |
| `func TishaBAv(gYear int) int` | TishaBAv returns the RD fixed day of Tisha be-Av (9 Av) occurring in the given Gregorian year, postponed to 10 Av when 9 Av falls on the Sabbath. |
| `func TuBShevat(gYear int) int` | TuBShevat returns the RD fixed day of Tu BiShvat (15 Shevat) occurring in the given Gregorian year. |
| `func USDaylightSavingEnd(year int) int` | USDaylightSavingEnd returns the RD fixed day on which United States daylight saving time ends (the first Sunday in November) in the given Gregorian… |
| `func USDaylightSavingStart(year int) int` | USDaylightSavingStart returns the RD fixed day on which United States daylight saving time begins (the second Sunday in March) in the given Gregorian… |
| `func USElectionDay(year int) int` | USElectionDay returns the RD fixed day of United States federal Election Day (the Tuesday after the first Monday in November) for the given Gregorian… |
| `func USIndependenceDay(year int) int` | USIndependenceDay returns the RD fixed day of United States Independence Day (4 July) for the given Gregorian year. |
| `func USLaborDay(year int) int` | USLaborDay returns the RD fixed day of United States Labor Day (the first Monday in September) for the given Gregorian year. |
| `func USMemorialDay(year int) int` | USMemorialDay returns the RD fixed day of United States Memorial Day (the last Monday in May) for the given Gregorian year. |
| `func USThanksgiving(year int) int` | USThanksgiving returns the RD fixed day of United States Thanksgiving (the fourth Thursday in November) for the given Gregorian year. |
| `func USThanksgivingCanada(year int) int` | USThanksgivingCanada returns the RD fixed day of Canadian Thanksgiving (the second Monday in October) for the given Gregorian year. |
| `func UnixDayFromFixed(fixed int) int` | UnixDayFromFixed converts an RD fixed-day count to the number of whole days since the Unix epoch (1 January 1970). |
| `func UnixSecondsFromFixed(fixed int) int64` | UnixSecondsFromFixed converts an RD fixed-day count (interpreted at midnight UTC) to the corresponding number of Unix seconds. |
| `func WeekDifference(fixed1, fixed2 int) int` | WeekDifference returns the signed number of whole weeks from fixed1 to fixed2 (truncated toward zero). |
| `func WeekOfYearGregorian(fixed int) int` | WeekOfYearGregorian returns the week number of an RD fixed day counting weeks from 1 January, where the first week is week 1 and weeks are… |
| `func WinterSolsticeJD(year int) float64` | WinterSolsticeJD returns the Julian Ephemeris Day of the December solstice of the given Gregorian year. |
| `func Yahrzeit(death HebrewDate, hYear int) int` | Yahrzeit returns the RD fixed day of the anniversary of a death in the given Hebrew year, applying the traditional rules for deaths on 30 Marheshvan,… |
| `func YahrzeitInGregorian(death HebrewDate, gYear int) []int` | YahrzeitInGregorian returns every RD fixed day in the given Gregorian year on which the yahrzeit of the given Hebrew death date falls (zero, one or… |
| `func YearMonthDayDifference(fixed1, fixed2 int) (years, months, days int)` | YearMonthDayDifference returns the difference between two RD fixed days broken into whole Gregorian years, months and days, with fixed1 <= fixed2… |
| `func YomKippur(gYear int) int` | YomKippur returns the RD fixed day of Yom Kippur (10 Tishri) occurring in the given Gregorian year. |

### Types

| Type | What it is |
| --- | --- |
| `CopticDate` | CopticDate is a date in the Coptic (Alexandrian) calendar. |
| `EthiopicDate` | EthiopicDate is a date in the Ethiopic calendar, which shares the Coptic structure but uses a different epoch. |
| `GregorianDate` | GregorianDate is a proleptic Gregorian date. |
| `HebrewDate` | HebrewDate is a date in the Hebrew (Jewish) calendar. |
| `ISODate` | ISODate is an ISO-8601 week date: a week-numbering year, an ISO week number (1–53) and a day of the week numbered 1 (Monday) through 7 (Sunday). |
| `IslamicDate` | IslamicDate is a date in the tabular (arithmetical) Islamic calendar. |
| `JulianDate` | JulianDate is a proleptic Julian-calendar date. |
| `MoonPhase` | MoonPhase identifies one of the four principal lunar phases. |
| `PersianDate` | PersianDate is a date in the arithmetical Persian (Jalali) calendar, which approximates the astronomical Persian calendar using a 2820-year… |
| `Season` | Season identifies one of the four astronomical season boundaries. |
| `Weekday` | Weekday identifies a day of the week. |

<details>
<summary><code>CopticDate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CopticFromFixed(fixed int) CopticDate` | CopticFromFixed converts an RD fixed day to a Coptic date. |
| `func NewCoptic(year, month, day int) (CopticDate, error)` | NewCoptic returns a validated CopticDate, or an error if the fields do not describe a real date. |
| `func (d CopticDate) Fixed() int` | Fixed returns the RD fixed day of the Coptic date. |
| `func (d CopticDate) JDN() int` | JDN returns the Julian Day Number of the Coptic date. |
| `func (d CopticDate) String() string` | String formats the Coptic date as "Day MonthName Year". |
| `func (d CopticDate) Valid() bool` | Valid reports whether the receiver is a valid Coptic date. |
| `func (d CopticDate) Weekday() Weekday` | Weekday returns the day of the week of the Coptic date. |

</details>

<details>
<summary><code>EthiopicDate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EthiopicFromFixed(fixed int) EthiopicDate` | EthiopicFromFixed converts an RD fixed day to an Ethiopic date. |
| `func NewEthiopic(year, month, day int) (EthiopicDate, error)` | NewEthiopic returns a validated EthiopicDate, or an error if the fields do not describe a real date. |
| `func (d EthiopicDate) Fixed() int` | Fixed returns the RD fixed day of the Ethiopic date. |
| `func (d EthiopicDate) JDN() int` | JDN returns the Julian Day Number of the Ethiopic date. |
| `func (d EthiopicDate) String() string` | String formats the Ethiopic date as "Day MonthName Year". |
| `func (d EthiopicDate) Valid() bool` | Valid reports whether the receiver is a valid Ethiopic date. |
| `func (d EthiopicDate) Weekday() Weekday` | Weekday returns the day of the week of the Ethiopic date. |

</details>

<details>
<summary><code>GregorianDate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GregorianEasterDate(year int) GregorianDate` | GregorianEasterDate returns Western Easter Sunday of the given year as a Gregorian date. |
| `func GregorianFromFixed(fixed int) GregorianDate` | GregorianFromFixed converts an RD fixed day to a Gregorian date. |
| `func JulianToGregorian(year, month, day int) GregorianDate` | JulianToGregorian converts a Julian-calendar date to the equivalent Gregorian date for the same instant. |
| `func NewGregorian(year, month, day int) (GregorianDate, error)` | NewGregorian returns a validated GregorianDate, or an error if the fields do not describe a real date. |
| `func OrthodoxEasterDate(year int) GregorianDate` | OrthodoxEasterDate returns Orthodox Easter Sunday of the given year as the equivalent Gregorian (civil) date. |
| `func (d GregorianDate) After(o GregorianDate) bool` | After reports whether the receiver follows the other Gregorian date. |
| `func (d GregorianDate) Before(o GregorianDate) bool` | Before reports whether the receiver precedes the other Gregorian date. |
| `func (d GregorianDate) Equal(o GregorianDate) bool` | Equal reports whether the receiver denotes the same day as the other date. |
| `func (d GregorianDate) Fixed() int` | Fixed returns the RD fixed day of the Gregorian date. |
| `func (d GregorianDate) JDN() int` | JDN returns the Julian Day Number of the Gregorian date. |
| `func (d GregorianDate) PlusDays(n int) GregorianDate` | PlusDays returns a new GregorianDate n days after the receiver. |
| `func (d GregorianDate) PlusMonths(n int) GregorianDate` | PlusMonths returns the GregorianDate n calendar months after the receiver. |
| `func (d GregorianDate) PlusYears(n int) GregorianDate` | PlusYears returns the GregorianDate n years after the receiver, clamping 29 February to 28 February in common years. |
| `func (d GregorianDate) String() string` | String formats the date in ISO-8601 extended form with an explicit year sign and at least four year digits (for example "+2024-03-20" or… |
| `func (d GregorianDate) Sub(o GregorianDate) int` | Sub returns the signed number of days from the other date to the receiver. |
| `func (d GregorianDate) Valid() bool` | Valid reports whether the receiver is a valid Gregorian date. |
| `func (d GregorianDate) Weekday() Weekday` | Weekday returns the day of the week of the Gregorian date. |

</details>

<details>
<summary><code>HebrewDate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HebrewFromFixed(fixed int) HebrewDate` | HebrewFromFixed converts an RD fixed day to a Hebrew date. |
| `func NewHebrew(year, month, day int) (HebrewDate, error)` | NewHebrew returns a validated HebrewDate, or an error if the fields do not describe a real Hebrew date. |
| `func (d HebrewDate) Fixed() int` | Fixed returns the RD fixed day of the Hebrew date. |
| `func (d HebrewDate) JDN() int` | JDN returns the Julian Day Number of the Hebrew date. |
| `func (d HebrewDate) String() string` | String formats the Hebrew date as "Day MonthName Year". |
| `func (d HebrewDate) Valid() bool` | Valid reports whether the receiver is a valid Hebrew date. |
| `func (d HebrewDate) Weekday() Weekday` | Weekday returns the day of the week of the Hebrew date. |

</details>

<details>
<summary><code>ISODate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ISOFromFixed(fixed int) ISODate` | ISOFromFixed converts an RD fixed day to an ISO week date. |
| `func NewISO(year, week, day int) (ISODate, error)` | NewISO returns a validated ISODate, or an error if the fields are not a real ISO week date. |
| `func (d ISODate) Fixed() int` | Fixed returns the RD fixed day of the ISO week date. |
| `func (d ISODate) JDN() int` | JDN returns the Julian Day Number of the ISO week date. |
| `func (d ISODate) String() string` | String formats the ISO week date as YYYY-Www-D. |
| `func (d ISODate) Valid() bool` | Valid reports whether the receiver is a valid ISO week date. |
| `func (d ISODate) Weekday() Weekday` | Weekday returns the day of the week of the ISO week date. |

</details>

<details>
<summary><code>IslamicDate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IslamicFromFixed(fixed int) IslamicDate` | IslamicFromFixed converts an RD fixed day to a tabular Islamic date. |
| `func NewIslamic(year, month, day int) (IslamicDate, error)` | NewIslamic returns a validated IslamicDate, or an error if the fields do not describe a real date. |
| `func (d IslamicDate) Fixed() int` | Fixed returns the RD fixed day of the Islamic date. |
| `func (d IslamicDate) JDN() int` | JDN returns the Julian Day Number of the Islamic date. |
| `func (d IslamicDate) String() string` | String formats the Islamic date as "Day MonthName Year AH". |
| `func (d IslamicDate) Valid() bool` | Valid reports whether the receiver is a valid Islamic date. |
| `func (d IslamicDate) Weekday() Weekday` | Weekday returns the day of the week of the Islamic date. |

</details>

<details>
<summary><code>JulianDate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GregorianToJulian(year, month, day int) JulianDate` | GregorianToJulian converts a Gregorian date to the equivalent Julian-calendar date for the same instant. |
| `func JulianFromFixed(fixed int) JulianDate` | JulianFromFixed converts an RD fixed day to a Julian-calendar date. |
| `func NewJulian(year, month, day int) (JulianDate, error)` | NewJulian returns a validated JulianDate, or an error if the fields do not describe a real date. |
| `func (d JulianDate) Fixed() int` | Fixed returns the RD fixed day of the Julian-calendar date. |
| `func (d JulianDate) JDN() int` | JDN returns the Julian Day Number of the Julian-calendar date. |
| `func (d JulianDate) String() string` | String formats the Julian date as YYYY-MM-DD. |
| `func (d JulianDate) Valid() bool` | Valid reports whether the receiver is a valid Julian-calendar date. |
| `func (d JulianDate) Weekday() Weekday` | Weekday returns the day of the week of the Julian-calendar date. |

</details>

<details>
<summary><code>PersianDate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPersian(year, month, day int) (PersianDate, error)` | NewPersian returns a validated PersianDate, or an error if the fields do not describe a real date. |
| `func PersianFromFixed(fixed int) PersianDate` | PersianFromFixed converts an RD fixed day to a Persian (Jalali) date. |
| `func (d PersianDate) Fixed() int` | Fixed returns the RD fixed day of the Persian date. |
| `func (d PersianDate) JDN() int` | JDN returns the Julian Day Number of the Persian date. |
| `func (d PersianDate) String() string` | String formats the Persian date as YYYY-MM-DD. |
| `func (d PersianDate) Valid() bool` | Valid reports whether the receiver is a valid Persian date. |
| `func (d PersianDate) Weekday() Weekday` | Weekday returns the day of the week of the Persian date. |

</details>

<details>
<summary><code>Weekday</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GregorianWeekday(year, month, day int) Weekday` | GregorianWeekday returns the day of the week of a Gregorian date. |
| `func WeekdayFromFixed(fixed int) Weekday` | WeekdayFromFixed returns the day of the week of an RD fixed day, with Sunday being zero. |
| `func (d Weekday) String() string` | String returns the English name of the weekday. |

</details>

### Constants

`Thoout`, `Paope`, `Athor`, `Koiak`, `Tobe`, `Meshir`, `Paremotep`, `Parmoute`, `Pashons`, `Paone`, `Epep`, `Mesore`, `Epagomenae`, `GregorianEpoch`, `JulianEpoch`, `ISOEpoch`, `HebrewEpoch`, `IslamicEpoch`, `PersianEpoch`, `CopticEpoch`, `EthiopicEpoch`, `MJDEpoch`, `UnixEpoch`, `January`, `February`, `March`, `April`, `May`, `June`, `July`, `August`, `September`, `October`, `November`, `December`, `Nisan`, `Iyyar`, `Sivan`, `Tammuz`, `Av`, `Elul`, `Tishri`, `Marheshvan`, `Kislev`, `Tevet`, `Shevat`, `Adar`, `AdarII`, `Muharram`, `Safar`, `RabiI`, `RabiII`, `JumadaI`, `JumadaII`, `Rajab`, `Shaban`, `Ramadan`, `Shawwal`, `DhuAlQadah`, `DhuAlHijjah`, `Farvardin`, `Ordibehesht`, `Khordad`, `Tir`, `Mordad`, `Shahrivar`, `Mehr`, `Aban`, `Azar`, `Dey`, `Bahman`, `Esfand`, `JDEpoch`, `MeanLunarMonth`, `MeanTropicalYear`

### Variables

`ErrRange`, `ErrInvalidDate`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/calendrical).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
