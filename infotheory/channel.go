package infotheory

import (
	"errors"
	"math"
)

// ErrChannelShape is returned when a channel transition matrix or an input
// distribution has an inconsistent or empty shape.
var ErrChannelShape = errors.New("infotheory: inconsistent channel shape")

// BSC models a binary symmetric channel: each transmitted bit is independently
// flipped with crossover probability Crossover, which must lie in [0,1].
type BSC struct {
	// Crossover is the probability that a transmitted bit is received flipped.
	Crossover float64
}

// Capacity returns the channel capacity of the binary symmetric channel in
// bits per channel use, C = 1 - H_b(p), where p is the crossover probability
// and H_b is the binary entropy function. It is achieved by a uniform input.
func (c BSC) Capacity() float64 {
	return 1 - BinaryEntropy(c.Crossover)
}

// Channel returns the discrete-memoryless Channel representation of the binary
// symmetric channel, a 2x2 transition matrix with rows [1-p, p] and [p, 1-p].
func (c BSC) Channel() Channel {
	p := c.Crossover
	return Channel{Transition: [][]float64{{1 - p, p}, {p, 1 - p}}}
}

// BSCCapacity returns the capacity in bits of a binary symmetric channel with
// crossover probability p. It is a convenience wrapper for BSC{p}.Capacity.
func BSCCapacity(p float64) float64 {
	return BSC{Crossover: p}.Capacity()
}

// BEC models a binary erasure channel: each transmitted bit is independently
// erased (replaced by an erasure symbol) with probability Erasure, which must
// lie in [0,1], and is otherwise received correctly.
type BEC struct {
	// Erasure is the probability that a transmitted bit is erased.
	Erasure float64
}

// Capacity returns the channel capacity of the binary erasure channel in bits
// per channel use, C = 1 - e, where e is the erasure probability. It is
// achieved by a uniform input.
func (c BEC) Capacity() float64 {
	return 1 - c.Erasure
}

// Channel returns the discrete-memoryless Channel representation of the binary
// erasure channel. The two inputs map to three outputs {0, 1, erasure} with
// transition rows [1-e, 0, e] and [0, 1-e, e].
func (c BEC) Channel() Channel {
	e := c.Erasure
	return Channel{Transition: [][]float64{{1 - e, 0, e}, {0, 1 - e, e}}}
}

// BECCapacity returns the capacity in bits of a binary erasure channel with
// erasure probability e. It is a convenience wrapper for BEC{e}.Capacity.
func BECCapacity(e float64) float64 {
	return BEC{Erasure: e}.Capacity()
}

// Channel is a discrete memoryless channel described by its transition
// probability matrix. Transition[i][j] is P(Y=j | X=i); each row must be a
// probability distribution over the output alphabet.
type Channel struct {
	// Transition holds the conditional output probabilities, indexed as
	// Transition[input][output].
	Transition [][]float64
}

// NumInputs returns the size of the channel's input alphabet.
func (ch Channel) NumInputs() int { return len(ch.Transition) }

// NumOutputs returns the size of the channel's output alphabet, or zero for a
// channel with no inputs.
func (ch Channel) NumOutputs() int {
	if len(ch.Transition) == 0 {
		return 0
	}
	return len(ch.Transition[0])
}

// OutputDistribution returns the distribution of the output Y induced by the
// input distribution input, that is P(Y=j) = sum_i input_i Transition[i][j]. It
// returns ErrChannelShape if input does not have one entry per channel input.
func (ch Channel) OutputDistribution(input []float64) ([]float64, error) {
	if len(input) != ch.NumInputs() || ch.NumInputs() == 0 {
		return nil, ErrChannelShape
	}
	py := make([]float64, ch.NumOutputs())
	for i, pi := range input {
		for j, t := range ch.Transition[i] {
			py[j] += pi * t
		}
	}
	return py, nil
}

// ConditionalEntropy returns the conditional entropy H(Y|X) in bits for the
// given input distribution, H(Y|X) = sum_i input_i H(Transition[i]). It returns
// ErrChannelShape for a mismatched input distribution.
func (ch Channel) ConditionalEntropy(input []float64) (float64, error) {
	if len(input) != ch.NumInputs() || ch.NumInputs() == 0 {
		return 0, ErrChannelShape
	}
	var h float64
	for i, pi := range input {
		h += pi * Entropy(ch.Transition[i])
	}
	return h, nil
}

// MutualInformation returns the mutual information I(X;Y) in bits between the
// channel input, distributed as input, and the channel output. It equals
// sum_i input_i sum_j T_ij log2(T_ij / py_j) with py the induced output
// distribution. It returns ErrChannelShape for a mismatched input.
func (ch Channel) MutualInformation(input []float64) (float64, error) {
	py, err := ch.OutputDistribution(input)
	if err != nil {
		return 0, err
	}
	var mi float64
	for i, pi := range input {
		if pi <= 0 {
			continue
		}
		for j, t := range ch.Transition[i] {
			if t <= 0 || py[j] <= 0 {
				continue
			}
			mi += pi * t * infotheoryLog2(t/py[j])
		}
	}
	if mi < 0 {
		return 0, nil
	}
	return mi, nil
}

// BlahutArimoto computes the capacity of the discrete memoryless channel ch by
// the Blahut-Arimoto algorithm. It returns the capacity in bits and the
// capacity-achieving input distribution. Iteration stops when the gap between
// the upper and lower capacity bounds falls below tol or after maxIter
// iterations, whichever comes first; a non-positive maxIter defaults to 1000
// and a non-positive tol to 1e-10. It returns ErrChannelShape for an empty
// channel.
func BlahutArimoto(ch Channel, tol float64, maxIter int) (capacity float64, input []float64, err error) {
	n := ch.NumInputs()
	if n == 0 || ch.NumOutputs() == 0 {
		return 0, nil, ErrChannelShape
	}
	if maxIter <= 0 {
		maxIter = 1000
	}
	if tol <= 0 {
		tol = 1e-10
	}
	px := UniformDistribution(n)
	d := make([]float64, n) // per-input KL divergence D(T_i || py) in bits
	for iter := 0; iter < maxIter; iter++ {
		py, _ := ch.OutputDistribution(px)
		var lower, upper float64
		for i := range px {
			d[i] = infotheoryChannelRowKL(ch.Transition[i], py)
			lower += px[i] * d[i]
			if d[i] > upper {
				upper = d[i]
			}
		}
		if upper-lower < tol {
			return lower, px, nil
		}
		// Multiplicative update: px_i <- px_i * 2^{d_i}, renormalised.
		var z float64
		next := make([]float64, n)
		for i := range px {
			next[i] = px[i] * math.Exp2(d[i])
			z += next[i]
		}
		if z <= 0 {
			return lower, px, nil
		}
		for i := range next {
			next[i] /= z
		}
		px = next
	}
	py, _ := ch.OutputDistribution(px)
	var lower float64
	for i := range px {
		lower += px[i] * infotheoryChannelRowKL(ch.Transition[i], py)
	}
	return lower, px, nil
}

// ChannelCapacity returns the capacity in bits of the discrete memoryless
// channel ch, computed via BlahutArimoto with default tolerance and iteration
// limit. The capacity-achieving input distribution is discarded.
func ChannelCapacity(ch Channel) float64 {
	c, _, _ := BlahutArimoto(ch, 0, 0)
	return c
}

// infotheoryChannelRowKL returns the KL divergence in bits of a transition row
// from the output distribution py, sum_j t_j log2(t_j/py_j), used as the update
// direction in the Blahut-Arimoto iteration.
func infotheoryChannelRowKL(t, py []float64) float64 {
	var d float64
	for j, tj := range t {
		if tj <= 0 || py[j] <= 0 {
			continue
		}
		d += tj * infotheoryLog2(tj/py[j])
	}
	return d
}
