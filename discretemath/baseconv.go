package discretemath

// discretemathDigits holds the canonical digit alphabet for radixes up to 36,
// using lowercase letters for digit values ten and above.
const discretemathDigits = "0123456789abcdefghijklmnopqrstuvwxyz"

// discretemathDigitValue returns the numeric value of an ASCII digit character
// for base conversion, accepting 0-9, a-z and A-Z. It returns -1 for any other
// byte.
func discretemathDigitValue(c byte) int {
	switch {
	case c >= '0' && c <= '9':
		return int(c - '0')
	case c >= 'a' && c <= 'z':
		return int(c-'a') + 10
	case c >= 'A' && c <= 'Z':
		return int(c-'A') + 10
	default:
		return -1
	}
}

// ToBaseUint returns the representation of n in the given base, which must be in
// the range 2..36. Digit values ten and above use lowercase letters. It returns
// an error for an out-of-range base.
func ToBaseUint(n uint64, base int) (string, error) {
	if base < 2 || base > 36 {
		return "", discretemathErrorf("ToBaseUint: base %d out of range [2,36]", base)
	}
	if n == 0 {
		return "0", nil
	}
	b := uint64(base)
	var buf [64]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = discretemathDigits[n%b]
		n /= b
	}
	return string(buf[i:]), nil
}

// ToBase returns the representation of the signed integer n in the given base,
// which must be in the range 2..36. Negative values are prefixed with a minus
// sign. It returns an error for an out-of-range base.
func ToBase(n int64, base int) (string, error) {
	if base < 2 || base > 36 {
		return "", discretemathErrorf("ToBase: base %d out of range [2,36]", base)
	}
	if n == 0 {
		return "0", nil
	}
	neg := n < 0
	var u uint64
	if neg {
		u = uint64(-(n + 1)) + 1
	} else {
		u = uint64(n)
	}
	s, _ := ToBaseUint(u, base)
	if neg {
		return "-" + s, nil
	}
	return s, nil
}

// FromBaseUint parses the string s as an unsigned integer in the given base,
// which must be in the range 2..36. An optional leading plus sign is accepted.
// It returns an error for an out-of-range base, an empty string, an invalid
// digit, or a value that overflows uint64.
func FromBaseUint(s string, base int) (uint64, error) {
	if base < 2 || base > 36 {
		return 0, discretemathErrorf("FromBaseUint: base %d out of range [2,36]", base)
	}
	i := 0
	if i < len(s) && s[i] == '+' {
		i++
	}
	if i >= len(s) {
		return 0, discretemathErrorf("FromBaseUint: no digits in %q", s)
	}
	b := uint64(base)
	var n uint64
	for ; i < len(s); i++ {
		d := discretemathDigitValue(s[i])
		if d < 0 || d >= base {
			return 0, discretemathErrorf("FromBaseUint: invalid digit %q for base %d", s[i], base)
		}
		hi := n
		n = n*b + uint64(d)
		if n < hi {
			return 0, discretemathErrorf("FromBaseUint: value %q overflows uint64", s)
		}
	}
	return n, nil
}

// FromBase parses the string s as a signed integer in the given base, which must
// be in the range 2..36. An optional leading plus or minus sign is accepted. It
// returns an error for an out-of-range base, an empty string, an invalid digit,
// or a value that does not fit in int64.
func FromBase(s string, base int) (int64, error) {
	if base < 2 || base > 36 {
		return 0, discretemathErrorf("FromBase: base %d out of range [2,36]", base)
	}
	neg := false
	rest := s
	if len(rest) > 0 && (rest[0] == '+' || rest[0] == '-') {
		neg = rest[0] == '-'
		rest = rest[1:]
	}
	u, err := FromBaseUint(rest, base)
	if err != nil {
		return 0, discretemathErrorf("FromBase: %v", err)
	}
	if neg {
		if u > 1<<63 {
			return 0, discretemathErrorf("FromBase: value %q underflows int64", s)
		}
		return -int64(u), nil
	}
	if u > 1<<63-1 {
		return 0, discretemathErrorf("FromBase: value %q overflows int64", s)
	}
	return int64(u), nil
}
