package codingtheory

import "errors"

// ErrNotPrimitive is returned when a polynomial supplied to NewField does not
// generate the full multiplicative group of GF(2^m).
var ErrNotPrimitive = errors.New("codingtheory: polynomial is not primitive")

// ErrFieldParam is returned for out-of-range field parameters.
var ErrFieldParam = errors.New("codingtheory: invalid field parameter")

// Field is a finite field GF(2^m) represented by log/antilog tables generated
// from a primitive polynomial. Elements are ints in [0, 2^m); bit i is the
// coefficient of x^i in the polynomial-basis representation, and the generator
// alpha is the field element x (numeric value 2). Add and Sub are bitwise xor;
// Mul, Div, Inv and Pow use the tables for constant-time table lookups.
type Field struct {
	m    int   // extension degree
	poly int   // primitive polynomial with the x^m bit set
	size int   // 2^m
	exp  []int // exp[i] = alpha^i for i in [0, 2*(size-1))
	log  []int // log[x] = discrete log of x base alpha; log[0] is unused (-1)
}

// NewField constructs GF(2^m) from the primitive polynomial poly, which must
// be given with its x^m bit set (for example x^4+x+1 is 0b10011 = 19). It
// returns ErrFieldParam for m outside [1,24] and ErrNotPrimitive if poly does
// not generate every non-zero element.
func NewField(m, poly int) (*Field, error) {
	if m < 1 || m > 24 {
		return nil, ErrFieldParam
	}
	size := 1 << uint(m)
	if poly < size || poly >= 2*size {
		return nil, ErrFieldParam
	}
	f := &Field{m: m, poly: poly, size: size}
	f.exp = make([]int, 2*(size-1))
	f.log = make([]int, size)
	for i := range f.log {
		f.log[i] = -1
	}
	x := 1
	for i := 0; i < size-1; i++ {
		if f.log[x] != -1 {
			// x repeated before exhausting the group: poly is not primitive.
			return nil, ErrNotPrimitive
		}
		f.exp[i] = x
		f.log[x] = i
		x <<= 1
		if x&size != 0 {
			x ^= poly
		}
	}
	if x != 1 {
		return nil, ErrNotPrimitive
	}
	for i := size - 1; i < 2*(size-1); i++ {
		f.exp[i] = f.exp[i-(size-1)]
	}
	return f, nil
}

// NewGF2 returns GF(2) itself using the trivial primitive polynomial x.
func NewGF2() *Field { f, _ := NewField(1, 0b10); return f }

// NewGF4 returns GF(2^2) built from x^2+x+1.
func NewGF4() *Field { f, _ := NewField(2, 0b111); return f }

// NewGF8 returns GF(2^3) built from x^3+x+1.
func NewGF8() *Field { f, _ := NewField(3, 0b1011); return f }

// NewGF16 returns GF(2^4) built from x^4+x+1.
func NewGF16() *Field { f, _ := NewField(4, 0b10011); return f }

// NewGF32 returns GF(2^5) built from x^5+x^2+1.
func NewGF32() *Field { f, _ := NewField(5, 0b100101); return f }

// NewGF64 returns GF(2^6) built from x^6+x+1.
func NewGF64() *Field { f, _ := NewField(6, 0b1000011); return f }

// NewGF128 returns GF(2^7) built from x^7+x^3+1.
func NewGF128() *Field { f, _ := NewField(7, 0b10001001); return f }

// NewGF256 returns GF(2^8) built from x^8+x^4+x^3+x^2+1 (the polynomial 0x11d
// used by QR codes and many Reed-Solomon systems).
func NewGF256() *Field { f, _ := NewField(8, 0x11d); return f }

// M returns the extension degree m of GF(2^m).
func (f *Field) M() int { return f.m }

// Size returns the number of elements 2^m in the field.
func (f *Field) Size() int { return f.size }

// Order returns the multiplicative order 2^m-1 of the field's cyclic group.
func (f *Field) Order() int { return f.size - 1 }

// Poly returns the primitive polynomial (with the x^m bit set) defining the
// field.
func (f *Field) Poly() int { return f.poly }

// Add returns a+b in GF(2^m), which is the bitwise exclusive-or.
func (f *Field) Add(a, b int) int { return a ^ b }

// Sub returns a-b in GF(2^m); over characteristic two this equals Add.
func (f *Field) Sub(a, b int) int { return a ^ b }

// Mul returns a*b in GF(2^m) via log/antilog tables.
func (f *Field) Mul(a, b int) int {
	if a == 0 || b == 0 {
		return 0
	}
	return f.exp[f.log[a]+f.log[b]]
}

// Div returns a/b in GF(2^m). It panics if b is zero.
func (f *Field) Div(a, b int) int {
	if b == 0 {
		panic("codingtheory: division by zero in GF(2^m)")
	}
	if a == 0 {
		return 0
	}
	idx := f.log[a] - f.log[b]
	if idx < 0 {
		idx += f.size - 1
	}
	return f.exp[idx]
}

// Inv returns the multiplicative inverse of a in GF(2^m). It panics if a is
// zero.
func (f *Field) Inv(a int) int {
	if a == 0 {
		panic("codingtheory: inverse of zero in GF(2^m)")
	}
	return f.exp[(f.size-1)-f.log[a]]
}

// Pow returns a raised to the integer power e in GF(2^m); e may be negative.
// Pow(0, 0) is defined to be 1 and Pow(0, e) for e>0 is 0.
func (f *Field) Pow(a, e int) int {
	if a == 0 {
		if e == 0 {
			return 1
		}
		return 0
	}
	order := f.size - 1
	le := (f.log[a] * (e % order)) % order
	le %= order
	if le < 0 {
		le += order
	}
	return f.exp[le]
}

// Exp returns alpha^e, the e-th power of the field's primitive element. e may
// be any integer and is reduced modulo 2^m-1.
func (f *Field) Exp(e int) int {
	order := f.size - 1
	e %= order
	if e < 0 {
		e += order
	}
	return f.exp[e]
}

// Log returns the discrete logarithm of a non-zero element a to base alpha,
// i.e. the e in [0, 2^m-1) with alpha^e = a. It panics if a is zero.
func (f *Field) Log(a int) int {
	if a == 0 {
		panic("codingtheory: logarithm of zero in GF(2^m)")
	}
	return f.log[a]
}

// Alpha returns the primitive element alpha (numeric value 2) of the field.
func (f *Field) Alpha() int { return 2 }

// Elements returns a fresh slice of all 2^m field elements in ascending order.
func (f *Field) Elements() []int {
	out := make([]int, f.size)
	for i := range out {
		out[i] = i
	}
	return out
}

// Trace returns the field trace Tr(a) = a + a^2 + a^4 + ... + a^(2^(m-1)) of an
// element a, an element of GF(2) (0 or 1).
func (f *Field) Trace(a int) int {
	t := a
	x := a
	for i := 1; i < f.m; i++ {
		x = f.Mul(x, x)
		t ^= x
	}
	return t
}
