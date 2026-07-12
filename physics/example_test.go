package physics_test

import (
	"fmt"

	"github.com/malcolmston/algebra/physics"
)

func ExampleConvert() {
	// One mile in metres (exact by the international definition of the mile).
	v, _ := physics.Convert(1, "mi", "m")
	fmt.Printf("%.3f m\n", v)
	// Output: 1609.344 m
}

func ExampleConvert_temperature() {
	// The boiling point of water at standard pressure.
	v, _ := physics.Convert(100, "C", "F")
	fmt.Printf("%.1f F\n", v)
	// Output: 212.0 F
}

func ExampleLookup() {
	c, _ := physics.Lookup("c")
	fmt.Printf("%s = %g %s\n", c.Name, c.Value, c.Unit)
	// Output: Speed of light = 2.99792458e+08 m/s
}

func ExampleCelsiusToKelvin() {
	fmt.Printf("%.2f K\n", physics.CelsiusToKelvin(0))
	// Output: 273.15 K
}

func ExampleDegToRad() {
	fmt.Printf("%.6f\n", physics.DegToRad(180))
	// Output: 3.141593
}

func ExampleMassEnergy() {
	// Rest energy of one kilogram, E = m·c².
	fmt.Printf("%.4e J\n", physics.MassEnergy(1))
	// Output: 8.9876e+16 J
}

func ExamplePhotonEnergy() {
	// Energy of a 5×10¹⁴ Hz (green-ish) photon.
	fmt.Printf("%.4e J\n", physics.PhotonEnergy(5e14))
	// Output: 3.3130e-19 J
}

func ExampleLorentzFactor() {
	// At 60% of the speed of light the Lorentz factor is exactly 1.25.
	fmt.Printf("%.4f\n", physics.LorentzFactor(0.6*physics.SpeedOfLight))
	// Output: 1.2500
}
