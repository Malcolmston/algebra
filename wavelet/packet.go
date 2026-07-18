package wavelet

import (
	"errors"
	"math"
)

// PacketTree is a full (balanced) wavelet packet decomposition. Unlike the
// standard wavelet transform, which recurses only on the low-pass branch, the
// packet transform recurses on both branches, producing a complete binary tree
// of coefficient nodes. Nodes are held in natural (frequency-unordered) layout:
// at level L there are 2^L nodes, and the children of node i at level L are
// nodes 2i (low-pass) and 2i+1 (high-pass) at level L+1.
type PacketTree struct {
	w      Wavelet
	length int
	// nodes[level][index] holds the coefficients of one node.
	nodes [][][]float64
}

// WaveletPacket computes a full wavelet packet decomposition of signal with
// wavelet w down to the requested number of levels. It returns an error if
// levels is not positive, if it exceeds [MaxLevel] for the signal length, or if
// any intermediate node length is not even.
func WaveletPacket(signal []float64, w Wavelet, levels int) (*PacketTree, error) {
	if levels < 1 {
		return nil, errors.New("wavelet: WaveletPacket requires levels >= 1")
	}
	if levels > MaxLevel(len(signal)) {
		return nil, errors.New("wavelet: levels exceeds the maximum for this signal length")
	}
	nodes := make([][][]float64, levels+1)
	nodes[0] = [][]float64{append([]float64(nil), signal...)}
	for lvl := 0; lvl < levels; lvl++ {
		count := 1 << lvl
		nodes[lvl+1] = make([][]float64, 2*count)
		for i := 0; i < count; i++ {
			parent := nodes[lvl][i]
			if len(parent)%2 != 0 {
				return nil, errors.New("wavelet: intermediate node length is not even")
			}
			a, d := wavelet1DForward(parent, w.lo, w.hi)
			nodes[lvl+1][2*i] = a
			nodes[lvl+1][2*i+1] = d
		}
	}
	return &PacketTree{w: w, length: len(signal), nodes: nodes}, nil
}

// Levels returns the depth of the packet tree.
func (t *PacketTree) Levels() int { return len(t.nodes) - 1 }

// NodeCount returns the number of nodes at the given level, which is 2^level.
// It panics if level is out of range.
func (t *PacketTree) NodeCount(level int) int {
	if level < 0 || level >= len(t.nodes) {
		panic("wavelet: NodeCount level out of range")
	}
	return len(t.nodes[level])
}

// Node returns a copy of the coefficients stored at the given level and index.
// It panics if either coordinate is out of range.
func (t *PacketTree) Node(level, index int) []float64 {
	if level < 0 || level >= len(t.nodes) {
		panic("wavelet: Node level out of range")
	}
	if index < 0 || index >= len(t.nodes[level]) {
		panic("wavelet: Node index out of range")
	}
	return append([]float64(nil), t.nodes[level][index]...)
}

// Leaves returns copies of every node at the deepest level, left to right.
// Concatenated, they hold the full set of leaf coefficients of the decomposition.
func (t *PacketTree) Leaves() [][]float64 {
	deepest := t.nodes[len(t.nodes)-1]
	out := make([][]float64, len(deepest))
	for i, n := range deepest {
		out[i] = append([]float64(nil), n...)
	}
	return out
}

// Reconstruct rebuilds the original signal from the leaf nodes of the tree by
// inverting the packet transform level by level, achieving perfect
// reconstruction to within floating-point rounding.
func (t *PacketTree) Reconstruct() []float64 {
	// Work on a copy so the tree is left intact.
	cur := make([][]float64, len(t.nodes[len(t.nodes)-1]))
	for i, n := range t.nodes[len(t.nodes)-1] {
		cur[i] = append([]float64(nil), n...)
	}
	for lvl := len(t.nodes) - 1; lvl > 0; lvl-- {
		parentCount := len(cur) / 2
		next := make([][]float64, parentCount)
		for i := 0; i < parentCount; i++ {
			next[i] = wavelet1DInverse(cur[2*i], cur[2*i+1], t.w.lo, t.w.hi)
		}
		cur = next
	}
	return cur[0]
}

// ShannonEntropy returns the normalized Shannon entropy of x, a
// Coifman-Wickerhauser cost measure of coefficient concentration. Writing
// p_i = x_i^2 / sum(x^2), the entropy is -sum p_i * ln(p_i), with the
// convention 0*ln(0) = 0. It is 0 when the energy is concentrated in a single
// coefficient and ln(n) when it is spread uniformly. It returns 0 for a
// zero-energy input.
func ShannonEntropy(x []float64) float64 {
	e := Energy(x)
	if e == 0 {
		return 0
	}
	var h float64
	for _, v := range x {
		p := v * v / e
		if p > 0 {
			h -= p * math.Log(p)
		}
	}
	return h
}

// LogEnergy returns the log-energy cost measure sum_i ln(x_i^2) taken over the
// non-zero coefficients of x, another additive Coifman-Wickerhauser cost used
// for best-basis selection.
func LogEnergy(x []float64) float64 {
	var s float64
	for _, v := range x {
		if v != 0 {
			s += math.Log(v * v)
		}
	}
	return s
}
