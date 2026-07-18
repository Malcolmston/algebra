package tensor

import (
	"fmt"
	"sort"
	"strings"
)

// Einsum evaluates an Einstein-summation expression over the given tensors.
//
// The spec is written in the familiar subscript notation: comma-separated input
// subscripts, optionally followed by "->" and the output subscripts, for example
// "ij,jk->ik" for a matrix product or "ii->" for a trace. Each subscript is a
// string of letters, one per axis of the corresponding tensor; a letter that
// repeats within a term selects that term's diagonal. Every letter shared across
// terms is contracted (summed) unless it appears in the output. Whitespace is
// ignored.
//
// If "->" is omitted, Einsum uses implicit mode: the output subscripts are the
// letters that appear exactly once across all inputs, in alphabetical order
// (so "ij,jk" yields "ik", and "ba" yields output "ab", the transpose).
//
// It returns [ErrSpec] for a malformed specification (wrong number of terms,
// a subscript length not matching a tensor's rank, an unknown or repeated output
// letter) and [ErrShape] when a shared letter is bound to axes of differing
// length. The result is a new tensor; contracting to no output letters yields a
// rank-0 scalar tensor.
func Einsum(spec string, tensors ...*Tensor) (*Tensor, error) {
	spec = strings.ReplaceAll(spec, " ", "")
	var lhs, rhs string
	hasArrow := strings.Contains(spec, "->")
	if hasArrow {
		parts := strings.SplitN(spec, "->", 2)
		lhs, rhs = parts[0], parts[1]
	} else {
		lhs = spec
	}
	if lhs == "" {
		return nil, fmt.Errorf("%w: empty subscript list", ErrSpec)
	}
	terms := strings.Split(lhs, ",")
	if len(terms) != len(tensors) {
		return nil, fmt.Errorf("%w: %d subscript terms for %d tensors", ErrSpec, len(terms), len(tensors))
	}

	sizes := make(map[rune]int)
	var order []rune // labels in first-appearance order
	for ti, term := range terms {
		labels := []rune(term)
		if len(labels) != tensors[ti].Rank() {
			return nil, fmt.Errorf("%w: term %q has %d labels for rank-%d tensor", ErrSpec, term, len(labels), tensors[ti].Rank())
		}
		for ax, lb := range labels {
			if !tensorIsLetter(lb) {
				return nil, fmt.Errorf("%w: subscript %q is not a letter", ErrSpec, string(lb))
			}
			d := tensors[ti].shape[ax]
			if s, ok := sizes[lb]; ok {
				if s != d {
					return nil, fmt.Errorf("%w: label %q bound to sizes %d and %d", ErrShape, string(lb), s, d)
				}
			} else {
				sizes[lb] = d
				order = append(order, lb)
			}
		}
	}

	var outLabels []rune
	if hasArrow {
		seen := make(map[rune]bool)
		for _, lb := range rhs {
			if _, ok := sizes[lb]; !ok {
				return nil, fmt.Errorf("%w: output label %q does not appear in inputs", ErrSpec, string(lb))
			}
			if seen[lb] {
				return nil, fmt.Errorf("%w: output label %q repeated", ErrSpec, string(lb))
			}
			seen[lb] = true
			outLabels = append(outLabels, lb)
		}
	} else {
		count := make(map[rune]int)
		for _, term := range terms {
			for _, lb := range term {
				count[lb]++
			}
		}
		for lb := range count {
			if count[lb] == 1 {
				outLabels = append(outLabels, lb)
			}
		}
		sort.Slice(outLabels, func(i, j int) bool { return outLabels[i] < outLabels[j] })
	}

	outShape := make([]int, len(outLabels))
	for i, lb := range outLabels {
		outShape[i] = sizes[lb]
	}
	var out *Tensor
	if len(outShape) == 0 {
		out = FromScalar(0)
	} else {
		out = New(outShape...)
	}

	// Precompute, for each term, the label rune at each of its axes.
	termLabels := make([][]rune, len(terms))
	for ti, term := range terms {
		termLabels[ti] = []rune(term)
	}

	assign := make(map[rune]int, len(order))
	oIdx := make([]int, len(outLabels))
	idxBuf := make([][]int, len(terms))
	for ti := range terms {
		idxBuf[ti] = make([]int, len(termLabels[ti]))
	}

	var rec func(k int)
	rec = func(k int) {
		if k == len(order) {
			prod := 1.0
			for ti := range terms {
				for ax, lb := range termLabels[ti] {
					idxBuf[ti][ax] = assign[lb]
				}
				prod *= tensors[ti].At(idxBuf[ti]...)
			}
			for i, lb := range outLabels {
				oIdx[i] = assign[lb]
			}
			out.addAt(prod, oIdx)
			return
		}
		lb := order[k]
		n := sizes[lb]
		for v := 0; v < n; v++ {
			assign[lb] = v
			rec(k + 1)
		}
	}
	rec(0)
	return out, nil
}

// tensorIsLetter reports whether r is an ASCII letter usable as an einsum
// subscript.
func tensorIsLetter(r rune) bool {
	return (r >= 'a' && r <= 'z') || (r >= 'A' && r <= 'Z')
}
