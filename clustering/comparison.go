package clustering

import (
	"math"
	"sort"
)

// ContingencyMatrix returns the contingency table between two labellings of the
// same samples, together with the sorted distinct labels of each. Entry [i][j]
// counts samples with the i-th true label and j-th predicted label.
func ContingencyMatrix(labelsTrue, labelsPred []int) (table [][]int, trueLabels, predLabels []int) {
	trueLabels = UniqueLabels(labelsTrue)
	predLabels = UniqueLabels(labelsPred)
	ti := make(map[int]int, len(trueLabels))
	for i, l := range trueLabels {
		ti[l] = i
	}
	pi := make(map[int]int, len(predLabels))
	for i, l := range predLabels {
		pi[l] = i
	}
	table = make([][]int, len(trueLabels))
	for i := range table {
		table[i] = make([]int, len(predLabels))
	}
	n := len(labelsTrue)
	if len(labelsPred) < n {
		n = len(labelsPred)
	}
	for s := 0; s < n; s++ {
		table[ti[labelsTrue[s]]][pi[labelsPred[s]]]++
	}
	return table, trueLabels, predLabels
}

// RandIndex returns the Rand index between two clusterings: the fraction of
// sample pairs that are either together in both clusterings or separate in both.
// It lies in [0, 1] with 1 for identical clusterings.
func RandIndex(labelsTrue, labelsPred []int) float64 {
	a, b, c, d := pairCounts(labelsTrue, labelsPred)
	total := a + b + c + d
	if total == 0 {
		return 1
	}
	return (a + d) / total
}

// AdjustedRandIndex returns the Rand index corrected for chance, in the range
// [-1, 1] (typically [0, 1]); 1 indicates identical clusterings and 0 indicates
// agreement no better than random labelling.
func AdjustedRandIndex(labelsTrue, labelsPred []int) float64 {
	table, _, _ := ContingencyMatrix(labelsTrue, labelsPred)
	var sumComb float64
	rowSums := make([]float64, len(table))
	colSums := make([]float64, 0)
	if len(table) > 0 {
		colSums = make([]float64, len(table[0]))
	}
	for i := range table {
		for j := range table[i] {
			sumComb += comb2(float64(table[i][j]))
			rowSums[i] += float64(table[i][j])
			colSums[j] += float64(table[i][j])
		}
	}
	var sumRow, sumCol float64
	for _, r := range rowSums {
		sumRow += comb2(r)
	}
	for _, c := range colSums {
		sumCol += comb2(c)
	}
	var n float64
	for _, r := range rowSums {
		n += r
	}
	totalComb := comb2(n)
	if totalComb == 0 {
		return 1
	}
	expected := sumRow * sumCol / totalComb
	maxIndex := 0.5 * (sumRow + sumCol)
	if maxIndex == expected {
		return 1
	}
	return (sumComb - expected) / (maxIndex - expected)
}

// pairCounts returns the four pair-agreement counts (a, b, c, d) used by the
// Rand and Fowlkes-Mallows indices: a = pairs together in both, b = together in
// true only, c = together in pred only, d = separate in both.
func pairCounts(labelsTrue, labelsPred []int) (a, b, c, d float64) {
	n := len(labelsTrue)
	if len(labelsPred) < n {
		n = len(labelsPred)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			sameTrue := labelsTrue[i] == labelsTrue[j]
			samePred := labelsPred[i] == labelsPred[j]
			switch {
			case sameTrue && samePred:
				a++
			case sameTrue && !samePred:
				b++
			case !sameTrue && samePred:
				c++
			default:
				d++
			}
		}
	}
	return a, b, c, d
}

func comb2(n float64) float64 {
	if n < 2 {
		return 0
	}
	return n * (n - 1) / 2
}

// FowlkesMallowsIndex returns the Fowlkes-Mallows index between two clusterings:
// the geometric mean of the pairwise precision and recall. It lies in [0, 1].
func FowlkesMallowsIndex(labelsTrue, labelsPred []int) float64 {
	a, b, c, _ := pairCounts(labelsTrue, labelsPred)
	if a == 0 {
		return 0
	}
	return a / math.Sqrt((a+b)*(a+c))
}

// MutualInformation returns the mutual information (in nats) between two
// clusterings viewed as discrete random variables.
func MutualInformation(labelsTrue, labelsPred []int) float64 {
	table, _, _ := ContingencyMatrix(labelsTrue, labelsPred)
	var n float64
	rowSums := make([]float64, len(table))
	var colSums []float64
	if len(table) > 0 {
		colSums = make([]float64, len(table[0]))
	}
	for i := range table {
		for j := range table[i] {
			v := float64(table[i][j])
			rowSums[i] += v
			colSums[j] += v
			n += v
		}
	}
	if n == 0 {
		return 0
	}
	var mi float64
	for i := range table {
		for j := range table[i] {
			nij := float64(table[i][j])
			if nij == 0 {
				continue
			}
			mi += (nij / n) * math.Log(nij*n/(rowSums[i]*colSums[j]))
		}
	}
	return mi
}

// entropyOfLabels returns the Shannon entropy (nats) of a labelling.
func entropyOfLabels(labels []int) float64 {
	counts := LabelCounts(labels)
	n := float64(len(labels))
	if n == 0 {
		return 0
	}
	var h float64
	for _, c := range counts {
		p := float64(c) / n
		if p > 0 {
			h -= p * math.Log(p)
		}
	}
	return h
}

// NormalizedMutualInformation returns the mutual information normalised by the
// arithmetic mean of the two cluster entropies, in the range [0, 1].
func NormalizedMutualInformation(labelsTrue, labelsPred []int) float64 {
	mi := MutualInformation(labelsTrue, labelsPred)
	ht := entropyOfLabels(labelsTrue)
	hp := entropyOfLabels(labelsPred)
	denom := 0.5 * (ht + hp)
	if denom == 0 {
		return 1
	}
	return mi / denom
}

// Homogeneity returns the homogeneity score: 1 when each cluster contains only
// members of a single true class. It lies in [0, 1].
func Homogeneity(labelsTrue, labelsPred []int) float64 {
	hc := conditionalEntropy(labelsTrue, labelsPred)
	ht := entropyOfLabels(labelsTrue)
	if ht == 0 {
		return 1
	}
	return 1 - hc/ht
}

// Completeness returns the completeness score: 1 when all members of a true
// class are assigned to the same cluster. It lies in [0, 1].
func Completeness(labelsTrue, labelsPred []int) float64 {
	hk := conditionalEntropy(labelsPred, labelsTrue)
	hp := entropyOfLabels(labelsPred)
	if hp == 0 {
		return 1
	}
	return 1 - hk/hp
}

// VMeasure returns the V-measure: the harmonic mean of homogeneity and
// completeness.
func VMeasure(labelsTrue, labelsPred []int) float64 {
	h := Homogeneity(labelsTrue, labelsPred)
	c := Completeness(labelsTrue, labelsPred)
	if h+c == 0 {
		return 0
	}
	return 2 * h * c / (h + c)
}

// conditionalEntropy returns H(A|B), the entropy of labelling a given labelling
// b, in nats.
func conditionalEntropy(a, b []int) float64 {
	table, _, _ := ContingencyMatrix(a, b)
	var n float64
	var colSums []float64
	if len(table) > 0 {
		colSums = make([]float64, len(table[0]))
	}
	for i := range table {
		for j := range table[i] {
			v := float64(table[i][j])
			colSums[j] += v
			n += v
		}
	}
	if n == 0 {
		return 0
	}
	var h float64
	for i := range table {
		for j := range table[i] {
			nij := float64(table[i][j])
			if nij == 0 {
				continue
			}
			h -= (nij / n) * math.Log(nij/colSums[j])
		}
	}
	return h
}

// Purity returns the clustering purity: for each predicted cluster the count of
// its most common true label is summed and divided by the total number of
// samples. It lies in [0, 1].
func Purity(labelsTrue, labelsPred []int) float64 {
	groups := ClusterIndices(labelsPred)
	n := len(labelsPred)
	if n == 0 {
		return 0
	}
	var correct int
	for _, members := range groups {
		counts := make(map[int]int)
		best := 0
		for _, m := range members {
			counts[labelsTrue[m]]++
			if counts[labelsTrue[m]] > best {
				best = counts[labelsTrue[m]]
			}
		}
		correct += best
	}
	return float64(correct) / float64(n)
}

// AdjustedMutualInformation returns the mutual information adjusted for chance
// using the expected mutual information under a hypergeometric model of
// randomness, normalised by the maximum of the two entropies. It lies in
// [-1, 1] with 1 for identical clusterings.
func AdjustedMutualInformation(labelsTrue, labelsPred []int) float64 {
	table, _, _ := ContingencyMatrix(labelsTrue, labelsPred)
	var n float64
	rowSums := make([]float64, len(table))
	var colSums []float64
	if len(table) > 0 {
		colSums = make([]float64, len(table[0]))
	}
	for i := range table {
		for j := range table[i] {
			v := float64(table[i][j])
			rowSums[i] += v
			colSums[j] += v
			n += v
		}
	}
	if n == 0 {
		return 1
	}
	mi := MutualInformation(labelsTrue, labelsPred)
	emi := expectedMutualInformation(rowSums, colSums, n)
	ht := entropyOfLabels(labelsTrue)
	hp := entropyOfLabels(labelsPred)
	maxH := math.Max(ht, hp)
	denom := maxH - emi
	if math.Abs(denom) < 1e-12 {
		return 1
	}
	return (mi - emi) / denom
}

// expectedMutualInformation computes the expected mutual information (nats)
// between two random labellings with the given marginal counts.
func expectedMutualInformation(a, b []float64, n float64) float64 {
	if n == 0 {
		return 0
	}
	logn := math.Log(n)
	var emi float64
	for i := range a {
		ai := a[i]
		for j := range b {
			bj := b[j]
			start := int(math.Max(1, ai+bj-n))
			end := int(math.Min(ai, bj))
			for nij := start; nij <= end; nij++ {
				fnij := float64(nij)
				term := (fnij / n) * (math.Log(fnij) + logn - math.Log(ai) - math.Log(bj))
				logProb := lgamma(ai+1) + lgamma(bj+1) + lgamma(n-ai+1) + lgamma(n-bj+1) -
					lgamma(n+1) - lgamma(fnij+1) - lgamma(ai-fnij+1) -
					lgamma(bj-fnij+1) - lgamma(n-ai-bj+fnij+1)
				emi += term * math.Exp(logProb)
			}
		}
	}
	return emi
}

func lgamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}

// BestLabelMatching returns a permutation that relabels labelsPred to best match
// labelsTrue by greedily pairing predicted clusters to true classes in order of
// overlap. It returns a map from predicted label to matched true label.
func BestLabelMatching(labelsTrue, labelsPred []int) map[int]int {
	table, trueLabels, predLabels := ContingencyMatrix(labelsTrue, labelsPred)
	type cell struct {
		i, j, v int
	}
	cells := make([]cell, 0)
	for i := range table {
		for j := range table[i] {
			cells = append(cells, cell{i, j, table[i][j]})
		}
	}
	sort.Slice(cells, func(a, b int) bool { return cells[a].v > cells[b].v })
	usedTrue := make(map[int]bool)
	usedPred := make(map[int]bool)
	mapping := make(map[int]int)
	for _, c := range cells {
		tl := trueLabels[c.i]
		pl := predLabels[c.j]
		if usedTrue[tl] || usedPred[pl] {
			continue
		}
		usedTrue[tl] = true
		usedPred[pl] = true
		mapping[pl] = tl
	}
	// Any unmatched predicted labels map to themselves.
	for _, pl := range predLabels {
		if _, ok := mapping[pl]; !ok {
			mapping[pl] = pl
		}
	}
	return mapping
}

// ClusteringAccuracy returns the fraction of samples correctly labelled after
// optimally matching predicted clusters to true classes with BestLabelMatching.
func ClusteringAccuracy(labelsTrue, labelsPred []int) float64 {
	mapping := BestLabelMatching(labelsTrue, labelsPred)
	n := len(labelsTrue)
	if n == 0 {
		return 0
	}
	var correct int
	for i := 0; i < n; i++ {
		if mapping[labelsPred[i]] == labelsTrue[i] {
			correct++
		}
	}
	return float64(correct) / float64(n)
}
