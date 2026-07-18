package logic

import "testing"

func TestBooleanOps(t *testing.T) {
	cases := []struct {
		name   string
		fn     func(a, b bool) bool
		ff, ft bool
		tf, tt bool
	}{
		{"And", And, false, false, false, true},
		{"Or", Or, false, true, true, true},
		{"Xor", Xor, false, true, true, false},
		{"Nand", Nand, true, true, true, false},
		{"Nor", Nor, true, false, false, false},
		{"Xnor", Xnor, true, false, false, true},
		{"Implies", Implies, true, true, false, true},
		{"Iff", Iff, true, false, false, true},
	}
	for _, c := range cases {
		got := [4]bool{c.fn(false, false), c.fn(false, true), c.fn(true, false), c.fn(true, true)}
		want := [4]bool{c.ff, c.ft, c.tf, c.tt}
		if got != want {
			t.Errorf("%s: got %v, want %v", c.name, got, want)
		}
	}
}

func TestNot(t *testing.T) {
	if Not(true) || !Not(false) {
		t.Fatalf("Not is wrong")
	}
}

func TestMajority(t *testing.T) {
	want := map[[3]bool]bool{
		{false, false, false}: false,
		{true, false, false}:  false,
		{true, true, false}:   true,
		{true, true, true}:    true,
		{false, true, true}:   true,
	}
	for in, w := range want {
		if got := Majority(in[0], in[1], in[2]); got != w {
			t.Errorf("Majority%v = %v, want %v", in, got, w)
		}
	}
}

func TestMux(t *testing.T) {
	if Mux(false, true, false) != true {
		t.Errorf("Mux(false) should select a")
	}
	if Mux(true, true, false) != false {
		t.Errorf("Mux(true) should select b")
	}
}

func TestReductions(t *testing.T) {
	if !AndAll() || !AndAll(true, true) || AndAll(true, false) {
		t.Errorf("AndAll wrong")
	}
	if OrAll() || !OrAll(false, true) || OrAll(false, false) {
		t.Errorf("OrAll wrong")
	}
	if XorAll() || !XorAll(true, false, false) || XorAll(true, true) {
		t.Errorf("XorAll wrong")
	}
	// Parity of three trues is odd -> true.
	if !XorAll(true, true, true) {
		t.Errorf("XorAll parity wrong")
	}
}
