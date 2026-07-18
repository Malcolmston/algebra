package physics

import "testing"

// TestExtendedConstantsValues pins the SI value carried by each exported
// extended constant, guarding against accidental edits to the literals.
func TestExtendedConstantsValues(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"MolarMassConstant", MolarMassConstant, 0.99999999965e-3},
		{"ClassicalElectronRadius", ClassicalElectronRadius, 2.8179403262e-15},
		{"ComptonWavelength", ComptonWavelength, 2.42631023867e-12},
		{"HartreeEnergy", HartreeEnergy, 4.3597447222071e-18},
		{"MagneticFluxQuantum", MagneticFluxQuantum, 2.067833848e-15},
		{"ConductanceQuantum", ConductanceQuantum, 7.748091729e-5},
		{"JosephsonConstant", JosephsonConstant, 483597.8484e9},
		{"VonKlitzingConstant", VonKlitzingConstant, 25812.80745},
		{"WienDisplacement", WienDisplacement, 2.897771955e-3},
		{"FirstRadiation", FirstRadiation, 3.741771852e-16},
		{"SecondRadiation", SecondRadiation, 1.438776877e-2},
		{"ThomsonCrossSection", ThomsonCrossSection, 6.6524587321e-29},
		{"NuclearMagneton", NuclearMagneton, 5.0507837461e-27},
		{"BohrMagneton", BohrMagneton, 9.2740100783e-24},
		{"MuonMass", MuonMass, 1.883531627e-28},
		{"TauMass", TauMass, 3.16754e-27},
		{"PlanckLength", PlanckLength, 1.616255e-35},
		{"PlanckTime", PlanckTime, 5.391247e-44},
		{"PlanckMass", PlanckMass, 2.176434e-8},
		{"PlanckTemperature", PlanckTemperature, 1.416784e32},
		{"StandardAtmosphere", StandardAtmosphere, 101325.0},
		{"MolarVolumeSTP", MolarVolumeSTP, 22.41396954e-3},
	}
	for _, tc := range cases {
		if tc.got != tc.want {
			t.Errorf("%s = %g, want %g", tc.name, tc.got, tc.want)
		}
	}
}

// TestExtendedConstantsDefensiveCopy verifies ExtendedConstants returns a
// stable-ordered, non-aliasing copy of the internal table.
func TestExtendedConstantsDefensiveCopy(t *testing.T) {
	a := ExtendedConstants()
	b := ExtendedConstants()
	if len(a) != len(extendedConstants) {
		t.Fatalf("length = %d, want %d", len(a), len(extendedConstants))
	}
	// Deterministic ordering: two calls agree element-by-element.
	for i := range a {
		if a[i] != b[i] {
			t.Errorf("index %d differs between calls: %+v vs %+v", i, a[i], b[i])
		}
		if a[i] != extendedConstants[i] {
			t.Errorf("index %d = %+v, want %+v", i, a[i], extendedConstants[i])
		}
	}
	// Mutating the returned slice must not affect the internal table.
	a[0].Value = -1
	if extendedConstants[0].Value == -1 {
		t.Error("ExtendedConstants aliases the internal table")
	}
}

// TestLookupByName exercises O(1) name lookups across the base and extended
// tables, including a base entry, an extended entry, and a miss.
func TestLookupByName(t *testing.T) {
	cases := []struct {
		name   string
		want   float64
		wantOK bool
	}{
		{"Speed of light", SpeedOfLight, true},
		{"Bohr magneton", BohrMagneton, true},
		{"Planck length", PlanckLength, true},
		{"Standard atmosphere", StandardAtmosphere, true},
		{"Nonexistent constant", 0, false},
		{"c", 0, false}, // symbol is not a name
	}
	for _, tc := range cases {
		got, ok := LookupByName(tc.name)
		if ok != tc.wantOK {
			t.Errorf("LookupByName(%q) ok = %v, want %v", tc.name, ok, tc.wantOK)
			continue
		}
		if ok && got.Value != tc.want {
			t.Errorf("LookupByName(%q).Value = %g, want %g", tc.name, got.Value, tc.want)
		}
	}
}

// TestLookupAny exercises O(1) symbol lookups spanning the base and extended
// tables, and confirms parity with the base-only Lookup for base symbols.
func TestLookupAny(t *testing.T) {
	cases := []struct {
		symbol string
		want   float64
		wantOK bool
	}{
		{"c", SpeedOfLight, true},
		{"μ_B", BohrMagneton, true},
		{"Φ0", MagneticFluxQuantum, true},
		{"atm", StandardAtmosphere, true},
		{"V_m", MolarVolumeSTP, true},
		{"nope", 0, false},
	}
	for _, tc := range cases {
		got, ok := LookupAny(tc.symbol)
		if ok != tc.wantOK {
			t.Errorf("LookupAny(%q) ok = %v, want %v", tc.symbol, ok, tc.wantOK)
			continue
		}
		if ok && got.Value != tc.want {
			t.Errorf("LookupAny(%q).Value = %g, want %g", tc.symbol, got.Value, tc.want)
		}
	}

	// A base symbol must resolve identically through LookupAny and Lookup.
	base, okBase := Lookup("c")
	any, okAny := LookupAny("c")
	if !okBase || !okAny || base != any {
		t.Errorf("LookupAny/Lookup parity for %q: %+v (%v) vs %+v (%v)", "c", any, okAny, base, okBase)
	}
}

// TestIndexesCoverEverything asserts every base and extended record is
// reachable through both O(1) indexes.
func TestIndexesCoverEverything(t *testing.T) {
	for _, c := range Constants() {
		if got, ok := LookupByName(c.Name); !ok || got != c {
			t.Errorf("base name %q not indexed correctly", c.Name)
		}
		if got, ok := LookupAny(c.Symbol); !ok || got != c {
			t.Errorf("base symbol %q not indexed correctly", c.Symbol)
		}
	}
	for _, c := range extendedConstants {
		if got, ok := LookupByName(c.Name); !ok || got != c {
			t.Errorf("extended name %q not indexed correctly", c.Name)
		}
		if got, ok := LookupAny(c.Symbol); !ok || got != c {
			t.Errorf("extended symbol %q not indexed correctly", c.Symbol)
		}
	}
}

func BenchmarkLookupByName(b *testing.B) {
	for i := 0; i < b.N; i++ {
		if _, ok := LookupByName("Planck length"); !ok {
			b.Fatal("missing")
		}
	}
}

func BenchmarkLookupAny(b *testing.B) {
	for i := 0; i < b.N; i++ {
		if _, ok := LookupAny("μ_B"); !ok {
			b.Fatal("missing")
		}
	}
}

// BenchmarkLookupBaseline measures the existing O(n) Lookup for comparison with
// the O(1) LookupAny above.
func BenchmarkLookupBaseline(b *testing.B) {
	for i := 0; i < b.N; i++ {
		if _, ok := Lookup("u"); !ok {
			b.Fatal("missing")
		}
	}
}
