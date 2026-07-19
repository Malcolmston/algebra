package codingtheory

import (
	"errors"
	"sort"
)

// This file implements Huffman source coding over integer symbol alphabets.

// ErrHuffman is returned for malformed Huffman inputs.
var ErrHuffman = errors.New("codingtheory: invalid Huffman input")

// HuffmanNode is a node of a Huffman code tree. Leaves carry a symbol and have
// nil children; internal nodes have two children and a Symbol of -1.
type HuffmanNode struct {
	Symbol int
	Weight int
	Left   *HuffmanNode
	Right  *HuffmanNode
}

// IsLeaf reports whether the node is a leaf.
func (n *HuffmanNode) IsLeaf() bool { return n.Left == nil && n.Right == nil }

// BuildHuffmanTree builds a Huffman tree for the given symbol weights (symbol
// id -> positive weight). Ties are broken deterministically by weight then
// symbol id so the tree is reproducible. It returns ErrHuffman if there are
// fewer than one symbol or a non-positive weight.
func BuildHuffmanTree(weights map[int]int) (*HuffmanNode, error) {
	if len(weights) == 0 {
		return nil, ErrHuffman
	}
	syms := make([]int, 0, len(weights))
	for s, w := range weights {
		if w <= 0 {
			return nil, ErrHuffman
		}
		syms = append(syms, s)
	}
	sort.Ints(syms)
	var nodes []*HuffmanNode
	for _, s := range syms {
		nodes = append(nodes, &HuffmanNode{Symbol: s, Weight: weights[s]})
	}
	if len(nodes) == 1 {
		// single symbol: give it a one-bit code via a padding parent.
		return &HuffmanNode{Symbol: -1, Weight: nodes[0].Weight, Left: nodes[0]}, nil
	}
	// order() picks the two lowest-weight nodes deterministically.
	extractMin := func() *HuffmanNode {
		best := 0
		for i := 1; i < len(nodes); i++ {
			if nodes[i].Weight < nodes[best].Weight ||
				(nodes[i].Weight == nodes[best].Weight && leafKey(nodes[i]) < leafKey(nodes[best])) {
				best = i
			}
		}
		n := nodes[best]
		nodes = append(nodes[:best], nodes[best+1:]...)
		return n
	}
	for len(nodes) > 1 {
		a := extractMin()
		b := extractMin()
		nodes = append(nodes, &HuffmanNode{Symbol: -1, Weight: a.Weight + b.Weight, Left: a, Right: b})
	}
	return nodes[0], nil
}

// leafKey returns a stable ordering key for tie-breaking (the smallest symbol
// id reachable from the node).
func leafKey(n *HuffmanNode) int {
	if n.IsLeaf() {
		return n.Symbol
	}
	l := leafKey(n.Left)
	if n.Right != nil {
		if r := leafKey(n.Right); r < l {
			l = r
		}
	}
	return l
}

// HuffmanCodes walks a Huffman tree and returns the code words as bit slices
// keyed by symbol id (left edge = 0, right edge = 1).
func HuffmanCodes(root *HuffmanNode) map[int][]int {
	codes := make(map[int][]int)
	var walk func(n *HuffmanNode, prefix []int)
	walk = func(n *HuffmanNode, prefix []int) {
		if n == nil {
			return
		}
		if n.IsLeaf() {
			codes[n.Symbol] = append([]int(nil), prefix...)
			return
		}
		walk(n.Left, append(prefix, 0))
		walk(n.Right, append(prefix, 1))
	}
	walk(root, nil)
	return codes
}

// BuildHuffmanCodes is a convenience wrapper returning the code map directly
// from symbol weights.
func BuildHuffmanCodes(weights map[int]int) (map[int][]int, error) {
	root, err := BuildHuffmanTree(weights)
	if err != nil {
		return nil, err
	}
	return HuffmanCodes(root), nil
}

// HuffmanEncode encodes a symbol sequence into a bit slice using the given code
// map. It returns ErrHuffman if a symbol has no code.
func HuffmanEncode(data []int, codes map[int][]int) ([]int, error) {
	var out []int
	for _, s := range data {
		cw, ok := codes[s]
		if !ok {
			return nil, ErrHuffman
		}
		out = append(out, cw...)
	}
	return out, nil
}

// HuffmanDecode decodes a bit slice back into symbols using a Huffman tree. It
// returns ErrHuffman if the bits do not correspond to a sequence of complete
// code words.
func HuffmanDecode(bits []int, root *HuffmanNode) ([]int, error) {
	if root == nil {
		return nil, ErrHuffman
	}
	var out []int
	node := root
	for _, b := range bits {
		if b == 0 {
			node = node.Left
		} else {
			node = node.Right
		}
		if node == nil {
			return nil, ErrHuffman
		}
		if node.IsLeaf() {
			out = append(out, node.Symbol)
			node = root
		}
	}
	if node != root {
		return nil, ErrHuffman
	}
	return out, nil
}

// HuffmanAverageLength returns the expected code-word length in bits per symbol
// for the given code map and symbol weights.
func HuffmanAverageLength(codes map[int][]int, weights map[int]int) float64 {
	total := 0
	sum := 0
	for s, w := range weights {
		total += w
		sum += w * len(codes[s])
	}
	if total == 0 {
		return 0
	}
	return float64(sum) / float64(total)
}
