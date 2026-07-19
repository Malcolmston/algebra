package codingtheory

import "math"

// This file implements rate-1/n binary convolutional codes with zero-tail
// termination and both hard-decision and soft-decision Viterbi decoding.

// ConvCode is a rate-1/n binary convolutional code with constraint length K.
// The n generator polynomials are K-bit integers whose bit K-1 is the tap on
// the current input and bit 0 the tap on the oldest of the K-1 memory bits.
type ConvCode struct {
	n          int
	constraint int
	polys      []int
}

// NewConvCode returns a rate-1/n convolutional code with the given constraint
// length and generator polynomials (each a K-bit integer). It returns
// ErrFieldParam for inconsistent parameters.
func NewConvCode(constraint int, polys []int) (*ConvCode, error) {
	if constraint < 1 || len(polys) < 1 {
		return nil, ErrFieldParam
	}
	limit := 1 << uint(constraint)
	for _, p := range polys {
		if p <= 0 || p >= limit {
			return nil, ErrFieldParam
		}
	}
	return &ConvCode{n: len(polys), constraint: constraint, polys: append([]int(nil), polys...)}, nil
}

// NewConvCodeOctal is a convenience constructor taking generator polynomials in
// octal (for example []int{7, 5} for the classic K=3 rate-1/2 code).
func NewConvCodeOctal(constraint int, octals []int) (*ConvCode, error) {
	polys := make([]int, len(octals))
	for i, o := range octals {
		polys[i] = octalToInt(o)
	}
	return NewConvCode(constraint, polys)
}

// octalToInt converts a base-10 integer whose digits are octal (e.g. 7, 15)
// into its binary value.
func octalToInt(o int) int {
	val := 0
	place := 1
	for o > 0 {
		val += (o % 10) * place
		place *= 8
		o /= 10
	}
	return val
}

// N returns the number of output bits per input bit (the inverse code rate).
func (c *ConvCode) N() int { return c.n }

// Constraint returns the constraint length K.
func (c *ConvCode) Constraint() int { return c.constraint }

// Rate returns the nominal code rate 1/n.
func (c *ConvCode) Rate() float64 { return 1.0 / float64(c.n) }

// NumStates returns the number of trellis states 2^(K-1).
func (c *ConvCode) NumStates() int { return 1 << uint(c.constraint-1) }

// parityOf returns the parity of the low bits of x.
func parityOf(x int) int {
	p := 0
	for x != 0 {
		p ^= x & 1
		x >>= 1
	}
	return p
}

// outputs computes the n output bits for the register value reg = (u<<(K-1))|s.
func (c *ConvCode) outputs(reg int) []int {
	out := make([]int, c.n)
	for j := 0; j < c.n; j++ {
		out[j] = parityOf(reg & c.polys[j])
	}
	return out
}

// Encode encodes a message with zero-tail termination, flushing the encoder
// with K-1 trailing zero bits. The output has length n*(len(msg)+K-1).
func (c *ConvCode) Encode(msg []int) ([]int, error) {
	if !validBits(msg) {
		return nil, ErrBit
	}
	input := append(append([]int(nil), msg...), make([]int, c.constraint-1)...)
	state := 0
	out := make([]int, 0, len(input)*c.n)
	for _, u := range input {
		reg := (u << uint(c.constraint-1)) | state
		out = append(out, c.outputs(reg)...)
		state = reg >> 1
	}
	return out, nil
}

// EncodeUnterminated encodes a message without flushing, leaving the encoder in
// whatever state it ends in. The output has length n*len(msg).
func (c *ConvCode) EncodeUnterminated(msg []int) ([]int, error) {
	if !validBits(msg) {
		return nil, ErrBit
	}
	state := 0
	out := make([]int, 0, len(msg)*c.n)
	for _, u := range msg {
		reg := (u << uint(c.constraint-1)) | state
		out = append(out, c.outputs(reg)...)
		state = reg >> 1
	}
	return out, nil
}

// HardViterbiDecode decodes a zero-terminated received bit stream using the
// Viterbi algorithm with a Hamming-distance branch metric, returning the
// original message (without the K-1 flush bits). The stream length must be a
// multiple of n.
func (c *ConvCode) HardViterbiDecode(received []int) ([]int, error) {
	if len(received)%c.n != 0 {
		return nil, ErrLength
	}
	steps := len(received) / c.n
	if steps < c.constraint-1 {
		return nil, ErrLength
	}
	nStates := c.NumStates()
	const inf = int(1) << 30
	metric := make([]int, nStates)
	for i := range metric {
		metric[i] = inf
	}
	metric[0] = 0
	prev := make([][]int, steps)
	inBit := make([][]int, steps)
	for t := 0; t < steps; t++ {
		recv := received[t*c.n : (t+1)*c.n]
		nm := make([]int, nStates)
		pp := make([]int, nStates)
		ib := make([]int, nStates)
		for i := range nm {
			nm[i] = inf
			pp[i] = -1
		}
		for s := 0; s < nStates; s++ {
			if metric[s] == inf {
				continue
			}
			for u := 0; u < 2; u++ {
				reg := (u << uint(c.constraint-1)) | s
				out := c.outputs(reg)
				bm := 0
				for j := 0; j < c.n; j++ {
					if out[j] != recv[j] {
						bm++
					}
				}
				ns := reg >> 1
				cand := metric[s] + bm
				if cand < nm[ns] {
					nm[ns] = cand
					pp[ns] = s
					ib[ns] = u
				}
			}
		}
		metric, prev[t], inBit[t] = nm, pp, ib
	}
	return c.traceback(prev, inBit, steps)
}

// SoftViterbiDecode decodes a zero-terminated stream of soft values using the
// Viterbi algorithm with a squared-Euclidean branch metric. Each received value
// is a real number in which larger positive values favour a transmitted 0 and
// larger negative values a transmitted 1 (antipodal mapping 0->+1, 1->-1). The
// stream length must be a multiple of n.
func (c *ConvCode) SoftViterbiDecode(received []float64) ([]int, error) {
	if len(received)%c.n != 0 {
		return nil, ErrLength
	}
	steps := len(received) / c.n
	if steps < c.constraint-1 {
		return nil, ErrLength
	}
	nStates := c.NumStates()
	inf := math.Inf(1)
	metric := make([]float64, nStates)
	for i := range metric {
		metric[i] = inf
	}
	metric[0] = 0
	prev := make([][]int, steps)
	inBit := make([][]int, steps)
	for t := 0; t < steps; t++ {
		recv := received[t*c.n : (t+1)*c.n]
		nm := make([]float64, nStates)
		pp := make([]int, nStates)
		ib := make([]int, nStates)
		for i := range nm {
			nm[i] = inf
			pp[i] = -1
		}
		for s := 0; s < nStates; s++ {
			if math.IsInf(metric[s], 1) {
				continue
			}
			for u := 0; u < 2; u++ {
				reg := (u << uint(c.constraint-1)) | s
				out := c.outputs(reg)
				bm := 0.0
				for j := 0; j < c.n; j++ {
					expected := 1.0
					if out[j] == 1 {
						expected = -1.0
					}
					d := recv[j] - expected
					bm += d * d
				}
				ns := reg >> 1
				cand := metric[s] + bm
				if cand < nm[ns] {
					nm[ns] = cand
					pp[ns] = s
					ib[ns] = u
				}
			}
		}
		metric, prev[t], inBit[t] = nm, pp, ib
	}
	return c.traceback(prev, inBit, steps)
}

// traceback reconstructs the input sequence ending in state 0 and drops the
// K-1 termination bits.
func (c *ConvCode) traceback(prev, inBit [][]int, steps int) ([]int, error) {
	state := 0 // terminated code ends in the zero state
	bits := make([]int, steps)
	for t := steps - 1; t >= 0; t-- {
		if prev[t][state] < 0 {
			return nil, ErrTooManyErrors
		}
		bits[t] = inBit[t][state]
		state = prev[t][state]
	}
	msgLen := steps - (c.constraint - 1)
	return bits[:msgLen], nil
}
