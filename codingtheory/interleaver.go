package codingtheory

// This file implements block, row-column, and convolutional interleavers used
// to disperse burst errors before decoding, plus their de-interleavers.

// BlockInterleave writes the input row by row into a rows-by-cols matrix and
// reads it out column by column, returning the interleaved sequence. The input
// length must equal rows*cols.
func BlockInterleave(data []int, rows, cols int) ([]int, error) {
	if rows < 1 || cols < 1 || len(data) != rows*cols {
		return nil, ErrLength
	}
	out := make([]int, len(data))
	idx := 0
	for c := 0; c < cols; c++ {
		for r := 0; r < rows; r++ {
			out[idx] = data[r*cols+c]
			idx++
		}
	}
	return out, nil
}

// BlockDeinterleave inverts BlockInterleave, restoring the original row-major
// order. The input length must equal rows*cols.
func BlockDeinterleave(data []int, rows, cols int) ([]int, error) {
	if rows < 1 || cols < 1 || len(data) != rows*cols {
		return nil, ErrLength
	}
	out := make([]int, len(data))
	idx := 0
	for c := 0; c < cols; c++ {
		for r := 0; r < rows; r++ {
			out[r*cols+c] = data[idx]
			idx++
		}
	}
	return out, nil
}

// InterleavePermutation returns the permutation applied by a rows-by-cols block
// interleaver, i.e. out[i] = input[perm[i]].
func InterleavePermutation(rows, cols int) []int {
	perm := make([]int, rows*cols)
	idx := 0
	for c := 0; c < cols; c++ {
		for r := 0; r < rows; r++ {
			perm[idx] = r*cols + c
			idx++
		}
	}
	return perm
}

// PermuteInts applies a permutation to a slice, returning out[i]=data[perm[i]].
// It returns ErrLength if the lengths disagree.
func PermuteInts(data, perm []int) ([]int, error) {
	if len(data) != len(perm) {
		return nil, ErrLength
	}
	out := make([]int, len(data))
	for i, p := range perm {
		if p < 0 || p >= len(data) {
			return nil, ErrLength
		}
		out[i] = data[p]
	}
	return out, nil
}

// InvertPermutation returns the inverse of a permutation.
func InvertPermutation(perm []int) []int {
	inv := make([]int, len(perm))
	for i, p := range perm {
		inv[p] = i
	}
	return inv
}

// ConvolutionalInterleaver models a classic (branches, delay) convolutional
// interleaver: branch i delays its symbols by i*delay cells. It is stateful; a
// matching ConvolutionalDeinterleaver reverses it after an end-to-end latency.
type ConvolutionalInterleaver struct {
	branches int
	delay    int
	regs     [][]int // regs[i] holds branch i's delay line
	pos      []int   // write position per branch
	branch   int     // current branch (round-robin)
}

// NewConvolutionalInterleaver returns a convolutional interleaver with the given
// number of branches and per-branch delay increment.
func NewConvolutionalInterleaver(branches, delay int) (*ConvolutionalInterleaver, error) {
	if branches < 1 || delay < 0 {
		return nil, ErrFieldParam
	}
	ci := &ConvolutionalInterleaver{branches: branches, delay: delay}
	ci.regs = make([][]int, branches)
	ci.pos = make([]int, branches)
	for i := range ci.regs {
		ci.regs[i] = make([]int, i*delay)
	}
	return ci, nil
}

// Push feeds one symbol into the interleaver and returns the symbol emitted on
// the current branch (0 for symbols still filling an empty delay line).
func (ci *ConvolutionalInterleaver) Push(sym int) int {
	b := ci.branch
	ci.branch = (ci.branch + 1) % ci.branches
	if len(ci.regs[b]) == 0 {
		return sym
	}
	out := ci.regs[b][ci.pos[b]]
	ci.regs[b][ci.pos[b]] = sym
	ci.pos[b] = (ci.pos[b] + 1) % len(ci.regs[b])
	return out
}

// Interleave runs a fresh convolutional interleaver over an entire sequence and
// returns the emitted sequence of equal length.
func Interleave(data []int, branches, delay int) ([]int, error) {
	ci, err := NewConvolutionalInterleaver(branches, delay)
	if err != nil {
		return nil, err
	}
	out := make([]int, len(data))
	for i, s := range data {
		out[i] = ci.Push(s)
	}
	return out, nil
}
