package infotheory

import (
	"container/heap"
	"errors"
	"sort"
)

// ErrHuffman is returned by the Huffman routines when their input is empty or
// otherwise malformed (for example encoding a symbol with no codeword or
// decoding an incomplete bit string).
var ErrHuffman = errors.New("infotheory: invalid huffman input")

// HuffmanNode is a node of a Huffman code tree. Leaf nodes carry a Symbol and
// its Weight; internal nodes carry the summed Weight of their subtree and have
// non-nil Left and Right children. The empty Symbol is used for internal nodes.
type HuffmanNode struct {
	// Symbol is the source symbol carried by a leaf; empty for internal nodes.
	Symbol string
	// Weight is the probability or frequency of the leaf, or the subtree total.
	Weight float64
	// Left and Right are the child subtrees; both nil exactly at a leaf.
	Left, Right *HuffmanNode
}

// IsLeaf reports whether the node is a leaf (has no children).
func (n *HuffmanNode) IsLeaf() bool {
	return n.Left == nil && n.Right == nil
}

// infotheoryHuffItem is an entry in the priority queue used while building a
// Huffman tree. The order field records insertion sequence so that ties in
// weight are broken deterministically, giving reproducible codes.
type infotheoryHuffItem struct {
	node  *HuffmanNode
	order int
}

type infotheoryHuffHeap []infotheoryHuffItem

// Len reports the number of items in the heap, implementing heap.Interface.
func (h infotheoryHuffHeap) Len() int { return len(h) }

// Less reports whether item i orders before item j, implementing heap.Interface.
// Items are ordered by ascending node weight, breaking ties by ascending
// insertion order.
func (h infotheoryHuffHeap) Less(i, j int) bool {
	if h[i].node.Weight != h[j].node.Weight {
		return h[i].node.Weight < h[j].node.Weight
	}
	return h[i].order < h[j].order
}

// Swap exchanges items i and j, implementing heap.Interface.
func (h infotheoryHuffHeap) Swap(i, j int) { h[i], h[j] = h[j], h[i] }

// Push appends x to the heap, implementing heap.Interface.
func (h *infotheoryHuffHeap) Push(x any) {
	*h = append(*h, x.(infotheoryHuffItem))
}

// Pop removes and returns the last item of the heap, implementing heap.Interface.
func (h *infotheoryHuffHeap) Pop() any {
	old := *h
	n := len(old)
	item := old[n-1]
	*h = old[:n-1]
	return item
}

// BuildHuffmanTree builds an optimal prefix-code tree for the symbols in
// weights, whose keys are symbols and values are non-negative weights
// (probabilities or frequencies). Symbols are processed in ascending key order
// and ties in weight are broken deterministically, so the tree is reproducible.
// It returns the root of the tree, or ErrHuffman if weights is empty.
func BuildHuffmanTree(weights map[string]float64) (*HuffmanNode, error) {
	if len(weights) == 0 {
		return nil, ErrHuffman
	}
	symbols := make([]string, 0, len(weights))
	for s := range weights {
		symbols = append(symbols, s)
	}
	sort.Strings(symbols)

	h := &infotheoryHuffHeap{}
	order := 0
	for _, s := range symbols {
		*h = append(*h, infotheoryHuffItem{node: &HuffmanNode{Symbol: s, Weight: weights[s]}, order: order})
		order++
	}
	heap.Init(h)

	for h.Len() > 1 {
		a := heap.Pop(h).(infotheoryHuffItem)
		b := heap.Pop(h).(infotheoryHuffItem)
		parent := &HuffmanNode{
			Weight: a.node.Weight + b.node.Weight,
			Left:   a.node,
			Right:  b.node,
		}
		heap.Push(h, infotheoryHuffItem{node: parent, order: order})
		order++
	}
	return heap.Pop(h).(infotheoryHuffItem).node, nil
}

// Codewords returns the binary codeword for every symbol in the tree rooted at
// n, assigning "0" to each left branch and "1" to each right branch. A tree
// consisting of a single leaf assigns that symbol the codeword "0". The
// returned map is keyed by symbol.
func (n *HuffmanNode) Codewords() map[string]string {
	codes := make(map[string]string)
	if n == nil {
		return codes
	}
	if n.IsLeaf() {
		codes[n.Symbol] = "0"
		return codes
	}
	var walk func(node *HuffmanNode, prefix string)
	walk = func(node *HuffmanNode, prefix string) {
		if node.IsLeaf() {
			codes[node.Symbol] = prefix
			return
		}
		if node.Left != nil {
			walk(node.Left, prefix+"0")
		}
		if node.Right != nil {
			walk(node.Right, prefix+"1")
		}
	}
	walk(n, "")
	return codes
}

// HuffmanCodes builds a Huffman tree for weights and returns the map from each
// symbol to its binary codeword. It is a convenience wrapper combining
// BuildHuffmanTree and Codewords, and returns ErrHuffman for empty input.
func HuffmanCodes(weights map[string]float64) (map[string]string, error) {
	root, err := BuildHuffmanTree(weights)
	if err != nil {
		return nil, err
	}
	return root.Codewords(), nil
}

// AverageCodeLength returns the expected codeword length sum_s p_s len(code_s)
// for the given codewords under the probability distribution probs, both keyed
// by symbol. Symbols present in probs but absent from codes are ignored.
func AverageCodeLength(codes map[string]string, probs map[string]float64) float64 {
	var avg float64
	for s, p := range probs {
		if c, ok := codes[s]; ok {
			avg += p * float64(len(c))
		}
	}
	return avg
}

// HuffmanEncode concatenates the codewords for the given sequence of symbols
// into a single bit string of '0' and '1' characters. It returns ErrHuffman if
// any symbol has no codeword in codes.
func HuffmanEncode(codes map[string]string, symbols []string) (string, error) {
	var out []byte
	for _, s := range symbols {
		c, ok := codes[s]
		if !ok {
			return "", ErrHuffman
		}
		out = append(out, c...)
	}
	return string(out), nil
}

// HuffmanDecode decodes a bit string of '0' and '1' characters back into the
// sequence of symbols, walking the tree rooted at root one bit at a time ('0'
// goes left, '1' goes right). It returns ErrHuffman if the tree is nil, a bit
// other than '0' or '1' is encountered, or the bit string ends in the middle of
// a codeword.
func HuffmanDecode(root *HuffmanNode, bits string) ([]string, error) {
	if root == nil {
		return nil, ErrHuffman
	}
	var out []string
	// Single-leaf tree: every '0' bit decodes to the sole symbol.
	if root.IsLeaf() {
		for i := 0; i < len(bits); i++ {
			if bits[i] != '0' {
				return nil, ErrHuffman
			}
			out = append(out, root.Symbol)
		}
		return out, nil
	}
	node := root
	for i := 0; i < len(bits); i++ {
		switch bits[i] {
		case '0':
			node = node.Left
		case '1':
			node = node.Right
		default:
			return nil, ErrHuffman
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
		return nil, ErrHuffman // trailing partial codeword
	}
	return out, nil
}

// ShannonFanoCodes returns a Shannon-Fano prefix code for the symbols in
// weights, keyed by symbol. Symbols are sorted by descending weight (ties
// broken by symbol) and recursively split into two groups of as-equal-as-
// possible total weight, the upper group receiving a '0' bit and the lower a
// '1'. Shannon-Fano codes are prefix-free but, unlike Huffman codes, are not
// always optimal. It returns ErrHuffman for empty input.
func ShannonFanoCodes(weights map[string]float64) (map[string]string, error) {
	if len(weights) == 0 {
		return nil, ErrHuffman
	}
	type sym struct {
		name string
		w    float64
	}
	syms := make([]sym, 0, len(weights))
	for s, w := range weights {
		syms = append(syms, sym{s, w})
	}
	sort.Slice(syms, func(i, j int) bool {
		if syms[i].w != syms[j].w {
			return syms[i].w > syms[j].w
		}
		return syms[i].name < syms[j].name
	})

	codes := make(map[string]string, len(syms))
	if len(syms) == 1 {
		codes[syms[0].name] = "0"
		return codes, nil
	}
	for _, s := range syms {
		codes[s.name] = ""
	}

	var split func(lo, hi int)
	split = func(lo, hi int) {
		if hi-lo < 1 {
			return
		}
		// Total weight of syms[lo..hi].
		var total float64
		for i := lo; i <= hi; i++ {
			total += syms[i].w
		}
		// Find the split point minimising the imbalance between the halves.
		var acc float64
		best := lo
		bestDiff := total
		for i := lo; i < hi; i++ {
			acc += syms[i].w
			diff := acc - (total - acc)
			if diff < 0 {
				diff = -diff
			}
			if diff < bestDiff {
				bestDiff = diff
				best = i
			}
		}
		for i := lo; i <= best; i++ {
			codes[syms[i].name] += "0"
		}
		for i := best + 1; i <= hi; i++ {
			codes[syms[i].name] += "1"
		}
		split(lo, best)
		split(best+1, hi)
	}
	split(0, len(syms)-1)
	return codes, nil
}
