package codingtheory

// This file implements configurable cyclic redundancy checks following the
// standard parameterised CRC model (width, polynomial, initial value, input
// and output reflection, final xor).

// CRC describes a parameterised cyclic redundancy check. Width is the register
// width in bits (1..64); Poly is the generator polynomial with its implicit
// x^Width term omitted (normal, non-reflected form); Init is the initial
// register value; RefIn reflects each input byte; RefOut reflects the final
// register; XorOut is xored into the final value.
type CRC struct {
	Width  int
	Poly   uint64
	Init   uint64
	RefIn  bool
	RefOut bool
	XorOut uint64
}

// reflectBits reverses the low n bits of v.
func reflectBits(v uint64, n int) uint64 {
	var r uint64
	for i := 0; i < n; i++ {
		if v&(1<<uint(i)) != 0 {
			r |= 1 << uint(n-1-i)
		}
	}
	return r
}

// mask returns a bitmask of the low Width bits.
func (c CRC) mask() uint64 {
	if c.Width >= 64 {
		return ^uint64(0)
	}
	return (1 << uint(c.Width)) - 1
}

// Checksum computes the CRC of the byte slice data under the configured model.
func (c CRC) Checksum(data []byte) uint64 {
	mask := c.mask()
	crc := c.Init & mask
	topbit := uint64(1) << uint(c.Width-1)
	for _, b := range data {
		bb := uint64(b)
		if c.RefIn {
			bb = reflectBits(bb, 8)
		}
		crc ^= bb << uint(c.Width-8)
		for i := 0; i < 8; i++ {
			if crc&topbit != 0 {
				crc = (crc << 1) ^ c.Poly
			} else {
				crc <<= 1
			}
			crc &= mask
		}
	}
	if c.RefOut {
		crc = reflectBits(crc, c.Width)
	}
	return (crc ^ c.XorOut) & mask
}

// ChecksumBits computes the CRC of a big-endian bit slice (each entry 0 or 1),
// processing the message bit by bit. Reflection settings are ignored; the
// polynomial division is performed directly. This is useful for non-byte-
// aligned messages.
func (c CRC) ChecksumBits(bits []int) uint64 {
	mask := c.mask()
	crc := c.Init & mask
	topbit := uint64(1) << uint(c.Width-1)
	for _, bit := range bits {
		in := uint64(bit&1) << uint(c.Width-1)
		msb := crc & topbit
		crc = (crc << 1) & mask
		if (msb == 0) != (in == 0) {
			crc ^= c.Poly
		}
		crc &= mask
	}
	return (crc ^ c.XorOut) & mask
}

// Verify appends and checks: it reports whether the CRC of data equals the
// given expected checksum.
func (c CRC) Verify(data []byte, expected uint64) bool {
	return c.Checksum(data) == expected
}

// Append returns data with the Width/8 checksum bytes appended in big-endian
// order. Width must be a multiple of eight.
func (c CRC) Append(data []byte) []byte {
	sum := c.Checksum(data)
	nbytes := c.Width / 8
	out := append([]byte(nil), data...)
	for i := nbytes - 1; i >= 0; i-- {
		out = append(out, byte(sum>>uint(8*i)))
	}
	return out
}

// CRC8 returns the standard CRC-8 model (polynomial 0x07, zero init, no
// reflection).
func CRC8() CRC { return CRC{Width: 8, Poly: 0x07} }

// CRC8CCITT is not defined here; use CRC8 or a custom model.

// CRC16CCITTFalse returns the CRC-16/CCITT-FALSE model (poly 0x1021,
// init 0xFFFF, no reflection).
func CRC16CCITTFalse() CRC { return CRC{Width: 16, Poly: 0x1021, Init: 0xFFFF} }

// CRC16XModem returns the CRC-16/XMODEM model (poly 0x1021, zero init).
func CRC16XModem() CRC { return CRC{Width: 16, Poly: 0x1021} }

// CRC16IBM returns the CRC-16/IBM (ARC) model (poly 0x8005, reflected in and
// out).
func CRC16IBM() CRC {
	return CRC{Width: 16, Poly: 0x8005, RefIn: true, RefOut: true}
}

// CRC32 returns the standard CRC-32 (as used by zlib/PNG): poly 0x04C11DB7,
// init and xorout 0xFFFFFFFF, reflected in and out.
func CRC32() CRC {
	return CRC{Width: 32, Poly: 0x04C11DB7, Init: 0xFFFFFFFF, RefIn: true, RefOut: true, XorOut: 0xFFFFFFFF}
}

// CRC32C returns the Castagnoli CRC-32C model (poly 0x1EDC6F41, reflected).
func CRC32C() CRC {
	return CRC{Width: 32, Poly: 0x1EDC6F41, Init: 0xFFFFFFFF, RefIn: true, RefOut: true, XorOut: 0xFFFFFFFF}
}
