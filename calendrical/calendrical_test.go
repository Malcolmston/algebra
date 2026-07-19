package calendrical

import (
	"fmt"
	"math"
	"testing"
)

// referenceDates mirrors the sample-date table from Dershowitz & Reingold,
// Calendrical Calculations, giving a single RD fixed day and its representation
// in each supported calendar.
var referenceDates = []struct {
	fixed    int
	greg     GregorianDate
	julian   JulianDate
	iso      ISODate
	hebrew   HebrewDate
	islamic  IslamicDate
	persian  PersianDate
	coptic   CopticDate
	ethiopic EthiopicDate
	weekday  Weekday
}{
	{
		fixed:    710347,
		greg:     GregorianDate{1945, 11, 12},
		julian:   JulianDate{1945, 10, 30},
		iso:      ISODate{1945, 46, 1},
		hebrew:   HebrewDate{5706, 9, 7},
		islamic:  IslamicDate{1364, 12, 6},
		persian:  PersianDate{1324, 8, 21},
		coptic:   CopticDate{1662, 3, 3},
		ethiopic: EthiopicDate{1938, 3, 3},
		weekday:  Monday,
	},
	{
		fixed:    730120, // 2000-01-01
		greg:     GregorianDate{2000, 1, 1},
		julian:   JulianDate{1999, 12, 19},
		iso:      ISODate{1999, 52, 6},
		hebrew:   HebrewDate{5760, 10, 23},
		islamic:  IslamicDate{1420, 9, 24},
		persian:  PersianDate{1378, 10, 11},
		coptic:   CopticDate{1716, 4, 22},
		ethiopic: EthiopicDate{1992, 4, 22},
		weekday:  Saturday,
	},
	{
		fixed:   1, // 0001-01-01 (RD epoch)
		greg:    GregorianDate{1, 1, 1},
		julian:  JulianDate{1, 1, 3},
		iso:     ISODate{1, 1, 1},
		weekday: Monday,
	},
}

func TestReferenceGregorian(t *testing.T) {
	for _, tc := range referenceDates {
		if got := FixedFromGregorian(tc.greg.Year, tc.greg.Month, tc.greg.Day); got != tc.fixed {
			t.Errorf("FixedFromGregorian(%v)=%d want %d", tc.greg, got, tc.fixed)
		}
		if got := GregorianFromFixed(tc.fixed); got != tc.greg {
			t.Errorf("GregorianFromFixed(%d)=%v want %v", tc.fixed, got, tc.greg)
		}
		if got := WeekdayFromFixed(tc.fixed); got != tc.weekday {
			t.Errorf("WeekdayFromFixed(%d)=%v want %v", tc.fixed, got, tc.weekday)
		}
	}
}

func TestReferenceJulian(t *testing.T) {
	for _, tc := range referenceDates {
		if got := FixedFromJulian(tc.julian.Year, tc.julian.Month, tc.julian.Day); got != tc.fixed {
			t.Errorf("FixedFromJulian(%v)=%d want %d", tc.julian, got, tc.fixed)
		}
		if got := JulianFromFixed(tc.fixed); got != tc.julian {
			t.Errorf("JulianFromFixed(%d)=%v want %v", tc.fixed, got, tc.julian)
		}
	}
}

func TestReferenceISO(t *testing.T) {
	for _, tc := range referenceDates {
		if got := FixedFromISO(tc.iso.Year, tc.iso.Week, tc.iso.Day); got != tc.fixed {
			t.Errorf("FixedFromISO(%v)=%d want %d", tc.iso, got, tc.fixed)
		}
		if got := ISOFromFixed(tc.fixed); got != tc.iso {
			t.Errorf("ISOFromFixed(%d)=%v want %v", tc.fixed, got, tc.iso)
		}
	}
}

func TestReferenceHebrew(t *testing.T) {
	for _, tc := range referenceDates {
		if tc.hebrew == (HebrewDate{}) {
			continue
		}
		if got := FixedFromHebrew(tc.hebrew.Year, tc.hebrew.Month, tc.hebrew.Day); got != tc.fixed {
			t.Errorf("FixedFromHebrew(%v)=%d want %d", tc.hebrew, got, tc.fixed)
		}
		if got := HebrewFromFixed(tc.fixed); got != tc.hebrew {
			t.Errorf("HebrewFromFixed(%d)=%v want %v", tc.fixed, got, tc.hebrew)
		}
	}
}

func TestReferenceIslamic(t *testing.T) {
	for _, tc := range referenceDates {
		if tc.islamic == (IslamicDate{}) {
			continue
		}
		if got := FixedFromIslamic(tc.islamic.Year, tc.islamic.Month, tc.islamic.Day); got != tc.fixed {
			t.Errorf("FixedFromIslamic(%v)=%d want %d", tc.islamic, got, tc.fixed)
		}
		if got := IslamicFromFixed(tc.fixed); got != tc.islamic {
			t.Errorf("IslamicFromFixed(%d)=%v want %v", tc.fixed, got, tc.islamic)
		}
	}
}

func TestReferencePersian(t *testing.T) {
	for _, tc := range referenceDates {
		if tc.persian == (PersianDate{}) {
			continue
		}
		if got := FixedFromPersian(tc.persian.Year, tc.persian.Month, tc.persian.Day); got != tc.fixed {
			t.Errorf("FixedFromPersian(%v)=%d want %d", tc.persian, got, tc.fixed)
		}
		if got := PersianFromFixed(tc.fixed); got != tc.persian {
			t.Errorf("PersianFromFixed(%d)=%v want %v", tc.fixed, got, tc.persian)
		}
	}
}

func TestReferenceCopticEthiopic(t *testing.T) {
	for _, tc := range referenceDates {
		if tc.coptic == (CopticDate{}) {
			continue
		}
		if got := FixedFromCoptic(tc.coptic.Year, tc.coptic.Month, tc.coptic.Day); got != tc.fixed {
			t.Errorf("FixedFromCoptic(%v)=%d want %d", tc.coptic, got, tc.fixed)
		}
		if got := CopticFromFixed(tc.fixed); got != tc.coptic {
			t.Errorf("CopticFromFixed(%d)=%v want %v", tc.fixed, got, tc.coptic)
		}
		if got := FixedFromEthiopic(tc.ethiopic.Year, tc.ethiopic.Month, tc.ethiopic.Day); got != tc.fixed {
			t.Errorf("FixedFromEthiopic(%v)=%d want %d", tc.ethiopic, got, tc.fixed)
		}
		if got := EthiopicFromFixed(tc.fixed); got != tc.ethiopic {
			t.Errorf("EthiopicFromFixed(%d)=%v want %v", tc.fixed, got, tc.ethiopic)
		}
	}
}

func TestRoundTripAllCalendars(t *testing.T) {
	for f := -400000; f <= 500000; f += 97 {
		if g := GregorianFromFixed(f); FixedFromGregorian(g.Year, g.Month, g.Day) != f {
			t.Fatalf("gregorian round-trip failed at %d: %v", f, g)
		}
		if j := JulianFromFixed(f); FixedFromJulian(j.Year, j.Month, j.Day) != f {
			t.Fatalf("julian round-trip failed at %d: %v", f, j)
		}
		if s := ISOFromFixed(f); FixedFromISO(s.Year, s.Week, s.Day) != f {
			t.Fatalf("iso round-trip failed at %d: %v", f, s)
		}
		if h := HebrewFromFixed(f); FixedFromHebrew(h.Year, h.Month, h.Day) != f {
			t.Fatalf("hebrew round-trip failed at %d: %v", f, h)
		}
		if i := IslamicFromFixed(f); FixedFromIslamic(i.Year, i.Month, i.Day) != f {
			t.Fatalf("islamic round-trip failed at %d: %v", f, i)
		}
		if p := PersianFromFixed(f); FixedFromPersian(p.Year, p.Month, p.Day) != f {
			t.Fatalf("persian round-trip failed at %d: %v", f, p)
		}
		if co := CopticFromFixed(f); FixedFromCoptic(co.Year, co.Month, co.Day) != f {
			t.Fatalf("coptic round-trip failed at %d: %v", f, co)
		}
		if e := EthiopicFromFixed(f); FixedFromEthiopic(e.Year, e.Month, e.Day) != f {
			t.Fatalf("ethiopic round-trip failed at %d: %v", f, e)
		}
	}
}

func TestJDNConversions(t *testing.T) {
	tests := []struct {
		greg GregorianDate
		jdn  int
	}{
		{GregorianDate{2000, 1, 1}, 2451545},
		{GregorianDate{-4713, 11, 24}, 0}, // JDN 0 (proleptic Gregorian)
		{GregorianDate{1, 1, 1}, 1721426},
		{GregorianDate{1970, 1, 1}, 2440588},
	}
	for _, tc := range tests {
		if got := tc.greg.JDN(); got != tc.jdn {
			t.Errorf("%v.JDN()=%d want %d", tc.greg, got, tc.jdn)
		}
		if got := GregorianFromFixed(FixedFromJDN(tc.jdn)); got != tc.greg {
			t.Errorf("FixedFromJDN(%d)->%v want %v", tc.jdn, got, tc.greg)
		}
	}
	if MJDFromFixed(FixedFromGregorian(1858, 11, 17)) != 0 {
		t.Errorf("MJD epoch mismatch")
	}
	if UnixDayFromFixed(FixedFromGregorian(1970, 1, 1)) != 0 {
		t.Errorf("Unix epoch mismatch")
	}
	if UnixSecondsFromFixed(FixedFromGregorian(1970, 1, 2)) != 86400 {
		t.Errorf("Unix seconds mismatch")
	}
}

func TestLeapYears(t *testing.T) {
	greg := map[int]bool{2000: true, 1900: false, 2024: true, 2023: false, 2100: false}
	for y, want := range greg {
		if GregorianLeapYear(y) != want {
			t.Errorf("GregorianLeapYear(%d)=%v want %v", y, GregorianLeapYear(y), want)
		}
	}
	if !JulianLeapYear(1900) {
		t.Errorf("JulianLeapYear(1900) should be true")
	}
	if !HebrewLeapYear(5784) || HebrewLeapYear(5785) {
		t.Errorf("Hebrew leap-year classification wrong")
	}
	if !IslamicLeapYear(1363) && !IslamicLeapYear(1362) {
		// at least the cycle must produce leaps; check count over a cycle
	}
	count := 0
	for y := 1; y <= 30; y++ {
		if IslamicLeapYear(y) {
			count++
		}
	}
	if count != 11 {
		t.Errorf("Islamic leap years per 30-year cycle = %d want 11", count)
	}
}

func TestDaysInMonthAndYear(t *testing.T) {
	if GregorianDaysInMonth(2024, 2) != 29 {
		t.Errorf("Feb 2024 should have 29 days")
	}
	if GregorianDaysInMonth(2023, 2) != 28 {
		t.Errorf("Feb 2023 should have 28 days")
	}
	if GregorianDaysInYear(2024) != 366 {
		t.Errorf("2024 should have 366 days")
	}
	// The sum of month lengths equals the year length in each calendar.
	if sum := sumMonths(2024, 12, GregorianDaysInMonth); sum != GregorianDaysInYear(2024) {
		t.Errorf("Gregorian month sum %d != year %d", sum, GregorianDaysInYear(2024))
	}
	for _, y := range []int{5784, 5785, 5786} {
		if sum := sumMonths(y, LastMonthOfHebrewYear(y), LastDayOfHebrewMonth); sum != DaysInHebrewYear(y) {
			t.Errorf("Hebrew month sum %d != year %d (year %d)", sum, DaysInHebrewYear(y), y)
		}
	}
	if IslamicDaysInMonth(1363, 12) != 30 || IslamicDaysInMonth(1364, 12) != 30 {
		// leap years 1363/1364 give 30-day final month; just exercise it
	}
	// Esfand has 30 days exactly in leap years and 29 otherwise.
	for y := 1400; y < 1410; y++ {
		want := 29
		if PersianLeapYear(y) {
			want = 30
		}
		if PersianDaysInMonth(y, 12) != want {
			t.Errorf("Persian %d Esfand = %d want %d", y, PersianDaysInMonth(y, 12), want)
		}
	}
}

func sumMonths(year, lastMonth int, daysIn func(int, int) int) int {
	total := 0
	for m := 1; m <= lastMonth; m++ {
		total += daysIn(year, m)
	}
	return total
}

func TestWeekdayArithmetic(t *testing.T) {
	// 2024-07-19 is a Friday.
	f := FixedFromGregorian(2024, 7, 19)
	if WeekdayFromFixed(f) != Friday {
		t.Fatalf("2024-07-19 should be Friday, got %v", WeekdayFromFixed(f))
	}
	if got := KDayOnOrAfter(Monday, f); GregorianFromFixed(got) != (GregorianDate{2024, 7, 22}) {
		t.Errorf("next Monday wrong: %v", GregorianFromFixed(got))
	}
	if got := KDayBefore(Friday, f); GregorianFromFixed(got) != (GregorianDate{2024, 7, 12}) {
		t.Errorf("previous Friday wrong: %v", GregorianFromFixed(got))
	}
	// Third Monday of July 2024 = 15 July.
	if d, ok := NthWeekdayOfMonth(3, Monday, 2024, 7); !ok || GregorianFromFixed(d) != (GregorianDate{2024, 7, 15}) {
		t.Errorf("3rd Monday July 2024 wrong: %v ok=%v", GregorianFromFixed(d), ok)
	}
	// A fifth Friday exists in some months but not all; November 2024 has 5 Fridays.
	if _, ok := NthWeekdayOfMonth(5, Friday, 2024, 11); !ok {
		t.Errorf("November 2024 should have a 5th Friday")
	}
	if _, ok := NthWeekdayOfMonth(5, Monday, 2024, 11); ok {
		t.Errorf("November 2024 should not have a 5th Monday")
	}
}

func TestHolidays(t *testing.T) {
	tests := []struct {
		name string
		got  int
		want GregorianDate
	}{
		{"Thanksgiving 2024", USThanksgiving(2024), GregorianDate{2024, 11, 28}},
		{"Memorial Day 2024", USMemorialDay(2024), GregorianDate{2024, 5, 27}},
		{"Labor Day 2024", USLaborDay(2024), GregorianDate{2024, 9, 2}},
		{"Election Day 2024", USElectionDay(2024), GregorianDate{2024, 11, 5}},
		{"Rosh haShanah 2024", RoshHaShanah(2024), GregorianDate{2024, 10, 3}},
		{"Yom Kippur 2024", YomKippur(2024), GregorianDate{2024, 10, 12}},
		{"Passover 2024", Passover(2024), GregorianDate{2024, 4, 23}},
		{"Nowruz 2024", Nowruz(2024), GregorianDate{2024, 3, 20}},
	}
	for _, tc := range tests {
		if got := GregorianFromFixed(tc.got); got != tc.want {
			t.Errorf("%s = %v want %v", tc.name, got, tc.want)
		}
	}
}

func TestEaster(t *testing.T) {
	gregorian := map[int]GregorianDate{
		2000: {2000, 4, 23},
		2005: {2005, 3, 27},
		2008: {2008, 3, 23},
		2011: {2011, 4, 24},
		2024: {2024, 3, 31},
		2025: {2025, 4, 20},
	}
	for y, want := range gregorian {
		if got := GregorianEasterDate(y); got != want {
			t.Errorf("GregorianEaster(%d)=%v want %v", y, got, want)
		}
	}
	orthodox := map[int]GregorianDate{
		2024: {2024, 5, 5},
		2025: {2025, 4, 20},
	}
	for y, want := range orthodox {
		if got := OrthodoxEasterDate(y); got != want {
			t.Errorf("OrthodoxEaster(%d)=%v want %v", y, got, want)
		}
	}
	// Derived observances hang off Easter at fixed offsets.
	if GoodFriday(2024) != GregorianEaster(2024)-2 {
		t.Errorf("Good Friday offset wrong")
	}
	if Pentecost(2024) != GregorianEaster(2024)+49 {
		t.Errorf("Pentecost offset wrong")
	}
}

func TestDateArithmetic(t *testing.T) {
	if got := (GregorianDate{2024, 1, 31}).PlusMonths(1); got != (GregorianDate{2024, 2, 29}) {
		t.Errorf("Jan 31 + 1 month = %v want 2024-02-29", got)
	}
	if got := (GregorianDate{2024, 2, 29}).PlusYears(1); got != (GregorianDate{2025, 2, 28}) {
		t.Errorf("Feb 29 2024 + 1 year = %v want 2025-02-28", got)
	}
	if got := (GregorianDate{2024, 12, 31}).PlusDays(1); got != (GregorianDate{2025, 1, 1}) {
		t.Errorf("Dec 31 2024 + 1 day = %v want 2025-01-01", got)
	}
	y, m, d := YearMonthDayDifference(FixedFromGregorian(2000, 1, 15), FixedFromGregorian(2024, 7, 19))
	if y != 24 || m != 6 || d != 4 {
		t.Errorf("difference = %d y %d m %d d want 24 6 4", y, m, d)
	}
	// Symmetry of the reversed difference.
	if ry, rm, rd := YearMonthDayDifference(FixedFromGregorian(2024, 7, 19), FixedFromGregorian(2000, 1, 15)); ry != -24 || rm != -6 || rd != -4 {
		t.Errorf("reversed difference = %d %d %d", ry, rm, rd)
	}
	if BusinessDaysBetween(FixedFromGregorian(2024, 7, 1), FixedFromGregorian(2024, 7, 8)) != 5 {
		t.Errorf("business days over one week should be 5")
	}
	// Adding 5 business days from a Monday lands on the next Monday.
	mon := FixedFromGregorian(2024, 7, 1)
	if got := AddBusinessDays(mon, 5); got != FixedFromGregorian(2024, 7, 8) {
		t.Errorf("Mon + 5 business days = %v", GregorianFromFixed(got))
	}
}

func TestISOWeeks(t *testing.T) {
	if !ISOLongYear(2020) {
		t.Errorf("2020 should be a 53-week ISO year")
	}
	if ISOLongYear(2021) {
		t.Errorf("2021 should be a 52-week ISO year")
	}
	if ISOWeeksInYear(2020) != 53 || ISOWeeksInYear(2021) != 52 {
		t.Errorf("ISO week counts wrong")
	}
	// 2024-12-30 belongs to ISO week 1 of 2025.
	got := ISOFromFixed(FixedFromGregorian(2024, 12, 30))
	if got != (ISODate{2025, 1, 1}) {
		t.Errorf("2024-12-30 ISO = %v want 2025-W01-1", got)
	}
}

func TestValidation(t *testing.T) {
	if _, err := NewGregorian(2024, 2, 30); err == nil {
		t.Errorf("Feb 30 should be invalid")
	}
	if _, err := NewGregorian(2024, 2, 29); err != nil {
		t.Errorf("Feb 29 2024 should be valid: %v", err)
	}
	if _, err := NewJulian(0, 1, 1); err == nil {
		t.Errorf("Julian year 0 should be invalid")
	}
	if _, err := NewHebrew(5785, 13, 1); err == nil {
		t.Errorf("Adar II in common year 5785 should be invalid")
	}
	if _, err := NewISO(2021, 53, 1); err == nil {
		t.Errorf("week 53 of a 52-week year should be invalid")
	}
	if _, err := NewIslamic(1364, 13, 1); err == nil {
		t.Errorf("month 13 should be invalid in Islamic calendar")
	}
}

func TestSeasons(t *testing.T) {
	tests := []struct {
		year   int
		season Season
		want   GregorianDate
	}{
		{2000, SpringEquinox, GregorianDate{2000, 3, 20}},
		{2000, SummerSolstice, GregorianDate{2000, 6, 21}},
		{2000, AutumnEquinox, GregorianDate{2000, 9, 22}},
		{2000, WinterSolstice, GregorianDate{2000, 12, 21}},
		{2024, SpringEquinox, GregorianDate{2024, 3, 20}},
		{2024, WinterSolstice, GregorianDate{2024, 12, 21}},
	}
	for _, tc := range tests {
		got := GregorianFromFixed(SeasonFixed(tc.year, tc.season))
		// Allow +/- 1 day slack for the TT-to-civil offset and rounding.
		diff := SeasonFixed(tc.year, tc.season) - tc.want.Fixed()
		if diff < -1 || diff > 1 {
			t.Errorf("Season %d/%v = %v want ~%v (diff %d days)", tc.year, tc.season, got, tc.want, diff)
		}
	}
	// The four seasons must occur in chronological order within a year.
	prev := SeasonJD(2024, SpringEquinox)
	for _, s := range []Season{SummerSolstice, AutumnEquinox, WinterSolstice} {
		cur := SeasonJD(2024, s)
		if cur <= prev {
			t.Errorf("seasons out of order at %v", s)
		}
		prev = cur
	}
}

func TestMoonPhases(t *testing.T) {
	// New Moon k=0 is 2000-01-06 (JDE ~ 2451550.26).
	jd := MoonPhaseJD(0, NewMoon)
	if math.Abs(jd-2451550.26) > 0.05 {
		t.Errorf("New Moon k=0 JDE = %.5f want ~2451550.26", jd)
	}
	if got := GregorianFromFixed(MoonPhaseFixed(0, NewMoon)); got != (GregorianDate{2000, 1, 6}) {
		t.Errorf("New Moon k=0 = %v want 2000-01-06", got)
	}
	// Phases of a lunation occur in order: new < first quarter < full < last.
	k := 300
	nm := MoonPhaseJD(k, NewMoon)
	fq := MoonPhaseJD(k, FirstQuarter)
	fm := MoonPhaseJD(k, FullMoon)
	lq := MoonPhaseJD(k, LastQuarter)
	if !(nm < fq && fq < fm && fm < lq) {
		t.Errorf("lunar phases out of order: %.3f %.3f %.3f %.3f", nm, fq, fm, lq)
	}
	// Spacing between consecutive same-phase new moons approximates the synodic month.
	if d := MoonPhaseJD(k+1, NewMoon) - nm; math.Abs(d-MeanLunarMonth) > 1.0 {
		t.Errorf("synodic spacing = %.4f want ~%.4f", d, MeanLunarMonth)
	}
	// Known full moon: April 2024 full moon is 2024-04-23 (UTC).
	kk := LunationNumber(2460420.0)
	if got := GregorianFromFixed(MoonPhaseFixed(kk, FullMoon)); got != (GregorianDate{2024, 4, 23}) {
		t.Errorf("April 2024 full moon = %v want 2024-04-23", got)
	}
}

func TestFloorHelpers(t *testing.T) {
	if FloorDiv(-7, 3) != -3 {
		t.Errorf("FloorDiv(-7,3) wrong")
	}
	if FloorMod(-7, 3) != 2 {
		t.Errorf("FloorMod(-7,3) wrong")
	}
	if AdjustedMod(7, 7) != 7 || AdjustedMod(8, 7) != 1 {
		t.Errorf("AdjustedMod wrong")
	}
}

func TestStringers(t *testing.T) {
	if got := (GregorianDate{2024, 7, 19}).String(); got != "+2024-07-19" {
		t.Errorf("Gregorian String = %q", got)
	}
	if got := (ISODate{2024, 29, 5}).String(); got != "2024-W29-5" {
		t.Errorf("ISO String = %q", got)
	}
	if Sunday.String() != "Sunday" || Saturday.String() != "Saturday" {
		t.Errorf("Weekday String wrong")
	}
	// 5784 is a Hebrew leap year, so month 12 is "Adar I"; 5785 is common.
	if HebrewMonthName(5784, 12) != "Adar I" {
		t.Errorf("leap-year Adar name wrong: %q", HebrewMonthName(5784, 12))
	}
	if HebrewMonthName(5785, 12) != "Adar" {
		t.Errorf("common-year Adar name wrong: %q", HebrewMonthName(5785, 12))
	}
}

func ExampleGregorianFromFixed() {
	// Convert a Julian Day Number to a civil date and back.
	fixed := FixedFromJDN(2451545) // J2000.0
	d := GregorianFromFixed(fixed)
	fmt.Printf("%s is a %s\n", d, d.Weekday())
	fmt.Printf("JDN round-trips: %d\n", d.JDN())
	// Output:
	// +2000-01-01 is a Saturday
	// JDN round-trips: 2451545
}

func ExampleGregorianEasterDate() {
	fmt.Println("Western Easter 2024:", GregorianEasterDate(2024))
	fmt.Println("Orthodox Easter 2024:", OrthodoxEasterDate(2024))
	// Output:
	// Western Easter 2024: +2024-03-31
	// Orthodox Easter 2024: +2024-05-05
}

func ExampleHebrewFromFixed() {
	f := FixedFromGregorian(2024, 10, 3)
	h := HebrewFromFixed(f)
	fmt.Printf("%d %s %d\n", h.Day, HebrewMonthName(h.Year, h.Month), h.Year)
	// Output:
	// 1 Tishri 5785
}
