package calendrical

// Weekday identifies a day of the week. The zero value is Sunday, matching the
// convention used throughout this package where the week begins on Sunday.
type Weekday int

// Days of the week.
const (
	Sunday Weekday = iota
	Monday
	Tuesday
	Wednesday
	Thursday
	Friday
	Saturday
)

var weekdayNames = [...]string{
	"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday",
}

// String returns the English name of the weekday.
func (d Weekday) String() string {
	if d < Sunday || d > Saturday {
		return "Weekday(?)"
	}
	return weekdayNames[d]
}

// WeekdayFromFixed returns the day of the week of an RD fixed day, with Sunday
// being zero.
func WeekdayFromFixed(fixed int) Weekday {
	return Weekday(floorMod(fixed, 7))
}

// DayOfWeekFromFixed returns the day-of-week index (0 = Sunday) of an RD fixed
// day. It is the integer form of WeekdayFromFixed.
func DayOfWeekFromFixed(fixed int) int {
	return floorMod(fixed, 7)
}

// KDayOnOrBefore returns the RD of the last day that is weekday k on or before
// the given fixed day.
func KDayOnOrBefore(k Weekday, fixed int) int {
	return fixed - floorMod(fixed-int(k), 7)
}

// KDayOnOrAfter returns the RD of the first day that is weekday k on or after
// the given fixed day.
func KDayOnOrAfter(k Weekday, fixed int) int {
	return KDayOnOrBefore(k, fixed+6)
}

// KDayNearest returns the RD of the day that is weekday k nearest to the given
// fixed day (ties resolve to the later day).
func KDayNearest(k Weekday, fixed int) int {
	return KDayOnOrBefore(k, fixed+3)
}

// KDayBefore returns the RD of the last day that is weekday k strictly before
// the given fixed day.
func KDayBefore(k Weekday, fixed int) int {
	return KDayOnOrBefore(k, fixed-1)
}

// KDayAfter returns the RD of the first day that is weekday k strictly after
// the given fixed day.
func KDayAfter(k Weekday, fixed int) int {
	return KDayOnOrBefore(k, fixed+7)
}

// NthKDay returns the RD of the n-th weekday k relative to the given fixed
// anchor day. When n is positive it counts forward (n = 1 is the first k on or
// after the anchor); when n is negative it counts backward (n = -1 is the last
// k on or before the anchor). n must not be zero.
func NthKDay(n int, k Weekday, fixed int) int {
	if n > 0 {
		return 7*n + KDayBefore(k, fixed)
	}
	return 7*n + KDayAfter(k, fixed)
}

// FirstKDay returns the RD of the first weekday k on or after the fixed anchor
// day. It is a convenience wrapper for NthKDay(1, k, fixed).
func FirstKDay(k Weekday, fixed int) int {
	return NthKDay(1, k, fixed)
}

// LastKDay returns the RD of the last weekday k on or before the fixed anchor
// day. It is a convenience wrapper for NthKDay(-1, k, fixed).
func LastKDay(k Weekday, fixed int) int {
	return NthKDay(-1, k, fixed)
}
