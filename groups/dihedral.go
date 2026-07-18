package groups

import "fmt"

// DihedralElement represents an element of a dihedral group D_n, the symmetry
// group of a regular n-gon with 2n elements. Every element is either a
// rotation r^k or a reflection s·r^k, encoded by the boolean Reflection and the
// rotation exponent Rotation taken modulo n. The identity is
// {Reflection: false, Rotation: 0}.
type DihedralElement struct {
	// Reflection is true for the n reflections and false for the n rotations.
	Reflection bool
	// Rotation is the rotation exponent k in [0, n).
	Rotation int
}

// DihedralOrder returns the order of the dihedral group D_n, which is 2n. The
// argument n must be positive.
func DihedralOrder(n int) int {
	if n <= 0 {
		panic("groups: DihedralOrder requires n > 0")
	}
	return 2 * n
}

// DihedralIdentity returns the identity element of every dihedral group, the
// rotation by zero.
func DihedralIdentity() DihedralElement {
	return DihedralElement{Reflection: false, Rotation: 0}
}

// DihedralRotation returns the rotation r^k in D_n with its exponent reduced
// modulo n. The argument n must be positive.
func DihedralRotation(n, k int) DihedralElement {
	if n <= 0 {
		panic("groups: DihedralRotation requires n > 0")
	}
	return DihedralElement{Reflection: false, Rotation: Mod(k, n)}
}

// DihedralReflection returns the reflection s·r^k in D_n with its exponent
// reduced modulo n. The argument n must be positive.
func DihedralReflection(n, k int) DihedralElement {
	if n <= 0 {
		panic("groups: DihedralReflection requires n > 0")
	}
	return DihedralElement{Reflection: true, Rotation: Mod(k, n)}
}

// DihedralCompose returns the product a·b in D_n using the defining relations
// r^n = s^2 = e and s·r = r^-1·s. Composition is applied left to right in the
// usual group sense (a followed by b under the convention a·b). The argument n
// must be positive.
func DihedralCompose(n int, a, b DihedralElement) DihedralElement {
	if n <= 0 {
		panic("groups: DihedralCompose requires n > 0")
	}
	// Write a = s^i r^j, b = s^k r^l with i,k in {0,1}.
	// s^i r^j · s^k r^l. Move r^j past s^k using r s = s r^-1, i.e.
	// r^j s^k = s^k r^{(-1)^k j}. So the product's reflection part is i xor k
	// and rotation part is (-1)^k * j + l.
	var rot int
	if b.Reflection {
		rot = Mod(-a.Rotation+b.Rotation, n)
	} else {
		rot = Mod(a.Rotation+b.Rotation, n)
	}
	return DihedralElement{
		Reflection: a.Reflection != b.Reflection,
		Rotation:   rot,
	}
}

// DihedralInverse returns the inverse of a in D_n. Reflections are their own
// inverse; the inverse of r^k is r^(n-k). The argument n must be positive.
func DihedralInverse(n int, a DihedralElement) DihedralElement {
	if n <= 0 {
		panic("groups: DihedralInverse requires n > 0")
	}
	if a.Reflection {
		return a
	}
	return DihedralElement{Reflection: false, Rotation: Mod(-a.Rotation, n)}
}

// DihedralElementOrder returns the order of a in D_n: 1 for the identity, 2 for
// every reflection, and n/gcd(k, n) for the rotation r^k. The argument n must
// be positive.
func DihedralElementOrder(n int, a DihedralElement) int {
	if n <= 0 {
		panic("groups: DihedralElementOrder requires n > 0")
	}
	if a.Reflection {
		return 2
	}
	if a.Rotation == 0 {
		return 1
	}
	return n / Gcd(a.Rotation, n)
}

// DihedralGroup returns all 2n elements of D_n: the n rotations r^0..r^(n-1)
// followed by the n reflections s·r^0..s·r^(n-1). The argument n must be
// positive.
func DihedralGroup(n int) []DihedralElement {
	if n <= 0 {
		panic("groups: DihedralGroup requires n > 0")
	}
	elems := make([]DihedralElement, 0, 2*n)
	for k := 0; k < n; k++ {
		elems = append(elems, DihedralElement{Reflection: false, Rotation: k})
	}
	for k := 0; k < n; k++ {
		elems = append(elems, DihedralElement{Reflection: true, Rotation: k})
	}
	return elems
}

// String renders a dihedral element as "r^k" for rotations and "s·r^k" for
// reflections, abbreviating the exponent-zero cases as "e", "r", "s".
func (a DihedralElement) String() string {
	if !a.Reflection {
		switch a.Rotation {
		case 0:
			return "e"
		case 1:
			return "r"
		default:
			return fmt.Sprintf("r^%d", a.Rotation)
		}
	}
	switch a.Rotation {
	case 0:
		return "s"
	case 1:
		return "s·r"
	default:
		return fmt.Sprintf("s·r^%d", a.Rotation)
	}
}
