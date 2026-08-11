# infotheory — implements classical information theory: entropy measures, divergences

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/infotheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/infotheory)

Package infotheory implements classical information theory: entropy measures,
divergences, discrete channel models and error-control codes.

The package is organised around discrete probability distributions represented
as plain []float64 slices whose entries are non-negative and sum to one, and
around joint distributions represented as row-major [][]float64 matrices whose
entries sum to one. Unless stated otherwise all entropic quantities are
expressed in bits (base-2 logarithms); companion routines suffixed with "Nat"
report the same quantity in nats (natural logarithms), and BitsToNats /
NatsToBits convert between the two.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/infotheory"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AverageCodeLength(codes map[string]string, probs map[string]float64) float64` | AverageCodeLength returns the expected codeword length sum_s p_s len(code_s) for the given codewords under the probability distribution probs, both… |
| `func BECCapacity(e float64) float64` | BECCapacity returns the capacity in bits of a binary erasure channel with erasure probability e. |
| `func BSCCapacity(p float64) float64` | BSCCapacity returns the capacity in bits of a binary symmetric channel with crossover probability p. |
| `func BhattacharyyaCoefficient(p, q []float64) float64` | BhattacharyyaCoefficient returns the Bhattacharyya coefficient BC(p,q) = sum_i sqrt(p_i q_i), the cosine of the angle between p and q under the… |
| `func BinaryEntropy(p float64) float64` | BinaryEntropy returns the binary entropy function H_b(p) = -p log2 p - (1-p) log2 (1-p), the entropy in bits of a Bernoulli random variable with… |
| `func BitsToNats(bits float64) float64` | BitsToNats converts an information quantity expressed in bits (base-2 units) to nats (natural-logarithm units) by multiplying by ln 2. |
| `func BlahutArimoto(ch Channel, tol float64, maxIter int) (capacity float64, input []float64, err error)` | BlahutArimoto computes the capacity of the discrete memoryless channel ch by the Blahut-Arimoto algorithm. |
| `func ChannelCapacity(ch Channel) float64` | ChannelCapacity returns the capacity in bits of the discrete memoryless channel ch, computed via BlahutArimoto with default tolerance and iteration… |
| `func CodeRate(k, n int) float64` | CodeRate returns the rate k/n of a block code that maps k information symbols to n coded symbols. |
| `func CodingEfficiency(entropyBits, avgLength float64) float64` | CodingEfficiency returns the efficiency of a code, the ratio of the source entropy in bits to the average codeword length: H / avgLength. |
| `func CollisionEntropy(p []float64) float64` | CollisionEntropy returns the Renyi entropy of order two, -log2(sum_i p_i^2), measured in bits. |
| `func ConditionalEntropyXgivenY(joint [][]float64) float64` | ConditionalEntropyXgivenY returns the conditional entropy H(X\|Y) = H(X,Y) - H(Y) of the joint distribution joint, measured in bits. |
| `func ConditionalEntropyYgivenX(joint [][]float64) float64` | ConditionalEntropyYgivenX returns the conditional entropy H(Y\|X) = H(X,Y) - H(X) of the joint distribution joint, measured in bits. |
| `func CrossEntropy(p, q []float64) float64` | CrossEntropy returns the cross entropy H(p,q) = -sum_i p_i log2 q_i between a true distribution p and a model distribution q, measured in bits. |
| `func CrossEntropyNat(p, q []float64) float64` | CrossEntropyNat returns the cross entropy H(p,q) measured in nats. |
| `func EmpiricalDistribution(samples []float64) (values []float64, dist []float64, err error)` | EmpiricalDistribution returns the empirical probability distribution of the samples together with the slice of distinct values, sorted ascending,… |
| `func Entropy(p []float64) float64` | Entropy returns the Shannon entropy H(p) = -sum_i p_i log2 p_i of the discrete probability distribution p, measured in bits. |
| `func EntropyBase(p []float64, base float64) float64` | EntropyBase returns the Shannon entropy of p measured in units of the given logarithm base (for example base 2 for bits, math.E for nats, base 10 for… |
| `func EntropyNat(p []float64) float64` | EntropyNat returns the Shannon entropy of p measured in nats (natural logarithm units): H(p) = -sum_i p_i ln p_i. |
| `func EntropyRate(p []float64, blockLength int) float64` | EntropyRate returns the entropy rate H(p)/L in bits per symbol of a source with per-block entropy H(p) emitting blocks of length L symbols. |
| `func FanoInequalityBound(errorProb float64, alphabet int) float64` | FanoInequalityBound returns the Fano inequality's upper bound on the conditional entropy H(X\|Y) in bits, given the probability of error errorProb of… |
| `func GiniImpurity(p []float64) float64` | GiniImpurity returns the Gini impurity 1 - sum_i p_i^2 of the distribution p, the probability that two independent samples fall in different… |
| `func Hamming74Decode(code byte) (data byte, errorPos int, corrected bool)` | Hamming74Decode decodes a 7-bit Hamming(7,4) codeword held in the low seven bits of code, correcting any single-bit error. |
| `func Hamming74Encode(data byte) byte` | Hamming74Encode encodes the four data bits in the low nibble of data (bit 0 is d1, bit 1 is d2, bit 2 is d3, bit 3 is d4) into a 7-bit Hamming(7,4)… |
| `func Hamming74Syndrome(code byte) int` | Hamming74Syndrome returns the error-position syndrome of a 7-bit Hamming(7,4) codeword held in the low seven bits of code. |
| `func HammingDistanceBits(a, b uint64) int` | HammingDistanceBits returns the Hamming distance between a and b, the number of bit positions in which their binary representations differ. |
| `func HammingDistanceBytes(a, b []byte) (int, error)` | HammingDistanceBytes returns the bit-level Hamming distance between two equal length byte slices, the total number of differing bits across all bytes. |
| `func HammingDistanceStrings(a, b string) (int, error)` | HammingDistanceStrings returns the symbol-level Hamming distance between two equal-length strings, the number of positions at which the corresponding… |
| `func HammingWeight(x uint64) int` | HammingWeight returns the Hamming weight of x, the number of 1-bits in its binary representation. |
| `func HartleyEntropy(p []float64) float64` | HartleyEntropy returns the Hartley (max-) entropy log2(k), where k is the number of outcomes with strictly positive probability in p, measured in… |
| `func HellingerDistance(p, q []float64) float64` | HellingerDistance returns the Hellinger distance sqrt( 1/2 sum_i (sqrt(p_i) - sqrt(q_i))^2 ) = sqrt(1 - BC(p,q)) between two distributions. |
| `func HuffmanCodes(weights map[string]float64) (map[string]string, error)` | HuffmanCodes builds a Huffman tree for weights and returns the map from each symbol to its binary codeword. |
| `func HuffmanDecode(root *HuffmanNode, bits string) ([]string, error)` | HuffmanDecode decodes a bit string of '0' and '1' characters back into the sequence of symbols, walking the tree rooted at root one bit at a time… |
| `func HuffmanEncode(codes map[string]string, symbols []string) (string, error)` | HuffmanEncode concatenates the codewords for the given sequence of symbols into a single bit string of '0' and '1' characters. |
| `func IsPrefixFree(codewords []string) bool` | IsPrefixFree reports whether the set of binary codewords is prefix-free: no codeword is a prefix of another. |
| `func IsProbabilityDistribution(p []float64, tol float64) bool` | IsProbabilityDistribution reports whether p consists of non-negative entries summing to one within the absolute tolerance tol. |
| `func JensenShannonDistance(p, q []float64) float64` | JensenShannonDistance returns the square root of the JensenShannonDivergence of p and q. |
| `func JensenShannonDivergence(p, q []float64) float64` | JensenShannonDivergence returns the Jensen-Shannon divergence between distributions p and q, measured in bits: JSD = 0.5 D(p\|\|m) + 0.5 D(q\|\|m)… |
| `func JointEntropy(joint [][]float64) float64` | JointEntropy returns the joint Shannon entropy H(X,Y) = -sum_ij p_ij log2 p_ij of the joint distribution joint, measured in bits. |
| `func KLDivergence(p, q []float64) float64` | KLDivergence returns the Kullback-Leibler divergence D(p\|\|q) = sum_i p_i log2(p_i/q_i) of distribution p from reference distribution q, measured in… |
| `func KLDivergenceNat(p, q []float64) float64` | KLDivergenceNat returns the Kullback-Leibler divergence D(p\|\|q) measured in nats (natural-logarithm units). |
| `func KraftEquality(lengths []int, base int) bool` | KraftEquality reports whether the codeword lengths make the Kraft sum equal to one within a small tolerance, the condition characterising a complete… |
| `func KraftInequality(lengths []int, base int) bool` | KraftInequality reports whether the codeword lengths satisfy the Kraft inequality sum_i base^{-l_i} <= 1, the necessary and sufficient condition for… |
| `func KraftSum(lengths []int, base int) float64` | KraftSum returns the Kraft sum sum_i base^{-l_i} of a code whose codeword lengths are l_i, over an alphabet of the given base (2 for binary). |
| `func MarginalX(joint [][]float64) []float64` | MarginalX returns the marginal distribution of the row variable X of the joint distribution given as a row-major matrix, where joint`i``j` is the… |
| `func MarginalY(joint [][]float64) []float64` | MarginalY returns the marginal distribution of the column variable Y of the joint distribution joint, where joint`i``j` is P(X=i, Y=j). |
| `func McMillanInequality(lengths []int, base int) bool` | McMillanInequality reports whether the codeword lengths satisfy the McMillan inequality sum_i base^{-l_i} <= 1. |
| `func MinEntropy(p []float64) float64` | MinEntropy returns the min-entropy -log2(max_i p_i) of the distribution p, measured in bits. |
| `func MutualInformation(joint [][]float64) float64` | MutualInformation returns the mutual information I(X;Y) = H(X) + H(Y) - H(X,Y) of the joint distribution joint, measured in bits. |
| `func NatsToBits(nats float64) float64` | NatsToBits converts an information quantity expressed in nats to bits by dividing by ln 2. |
| `func NormalizeCounts(counts []int) ([]float64, error)` | NormalizeCounts returns a probability distribution proportional to the supplied integer counts, dividing each count by their total. |
| `func NormalizeWeights(weights []float64) ([]float64, error)` | NormalizeWeights returns a probability distribution proportional to the supplied non-negative weights, that is each weight divided by their total. |
| `func NormalizedEntropy(p []float64) float64` | NormalizedEntropy returns the entropy efficiency H(p)/log2(n), where n is the number of outcomes len(p). |
| `func NormalizedMutualInformation(joint [][]float64) float64` | NormalizedMutualInformation returns the mutual information of joint scaled to [0,1] by the geometric mean of the marginal entropies: I(X;Y) /… |
| `func Parity(x uint64) int` | Parity returns the parity bit of x: 1 if x has an odd number of 1-bits and 0 otherwise (even parity). |
| `func Perplexity(p []float64) float64` | Perplexity returns the perplexity 2^H(p) of the distribution p, where H is the Shannon entropy in bits. |
| `func PinskerBound(tv float64) float64` | PinskerBound returns Pinsker's lower bound on the Kullback-Leibler divergence implied by a total variation distance tv, namely 2 tv^2 / ln 2 bits. |
| `func Redundancy(p []float64) float64` | Redundancy returns 1 - NormalizedEntropy(p), the fraction of the maximum entropy that is unused by the distribution p. |
| `func RenyiEntropy(p []float64, alpha float64) float64` | RenyiEntropy returns the Renyi entropy of order alpha of the distribution p, measured in bits. |
| `func ShannonFanoCodes(weights map[string]float64) (map[string]string, error)` | ShannonFanoCodes returns a Shannon-Fano prefix code for the symbols in weights, keyed by symbol. |
| `func Surprisal(p float64) float64` | Surprisal returns the self-information (surprisal) -log2 p of an outcome with probability p, measured in bits. |
| `func SurprisalBase(p, base float64) float64` | SurprisalBase returns the self-information of an outcome with probability p measured in units of the given logarithm base. |
| `func SurprisalNat(p float64) float64` | SurprisalNat returns the self-information of an outcome with probability p measured in nats: -ln p. |
| `func TotalVariationDistance(p, q []float64) float64` | TotalVariationDistance returns the total variation distance TV(p,q) = 1/2 sum_i \|p_i - q_i\| between two distributions, the largest difference the… |
| `func TsallisEntropy(p []float64, q float64) float64` | TsallisEntropy returns the Tsallis entropy of order q of the distribution p, S_q = (1 - sum_i p_i^q) / (q - 1). |
| `func UniformDistribution(n int) []float64` | UniformDistribution returns the uniform probability distribution over n outcomes, a slice of length n whose every entry is 1/n. |
| `func VariationOfInformation(joint [][]float64) float64` | VariationOfInformation returns the variation of information VI(X;Y) = H(X\|Y) + H(Y\|X) = 2 H(X,Y) - H(X) - H(Y) of the joint distribution joint,… |

### Types

| Type | What it is |
| --- | --- |
| `BEC` | BEC models a binary erasure channel: each transmitted bit is independently erased (replaced by an erasure symbol) with probability Erasure, which… |
| `BSC` | BSC models a binary symmetric channel: each transmitted bit is independently flipped with crossover probability Crossover, which must lie in [0,1]. |
| `Channel` | Channel is a discrete memoryless channel described by its transition probability matrix. |
| `HuffmanNode` | HuffmanNode is a node of a Huffman code tree. |

<details>
<summary><code>BEC</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c BEC) Capacity() float64` | Capacity returns the channel capacity of the binary erasure channel in bits per channel use, C = 1 - e, where e is the erasure probability. |
| `func (c BEC) Channel() Channel` | Channel returns the discrete-memoryless Channel representation of the binary erasure channel. |

</details>

<details>
<summary><code>BSC</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c BSC) Capacity() float64` | Capacity returns the channel capacity of the binary symmetric channel in bits per channel use, C = 1 - H_b(p), where p is the crossover probability… |
| `func (c BSC) Channel() Channel` | Channel returns the discrete-memoryless Channel representation of the binary symmetric channel, a 2x2 transition matrix with rows [1-p, p] and [p,… |

</details>

<details>
<summary><code>Channel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (ch Channel) ConditionalEntropy(input []float64) (float64, error)` | ConditionalEntropy returns the conditional entropy H(Y\|X) in bits for the given input distribution, H(Y\|X) = sum_i input_i H(Transition`i`). |
| `func (ch Channel) MutualInformation(input []float64) (float64, error)` | MutualInformation returns the mutual information I(X;Y) in bits between the channel input, distributed as input, and the channel output. |
| `func (ch Channel) NumInputs() int` | NumInputs returns the size of the channel's input alphabet. |
| `func (ch Channel) NumOutputs() int` | NumOutputs returns the size of the channel's output alphabet, or zero for a channel with no inputs. |
| `func (ch Channel) OutputDistribution(input []float64) ([]float64, error)` | OutputDistribution returns the distribution of the output Y induced by the input distribution input, that is P(Y=j) = sum_i input_i Transition`i``j`. |
| `func (ch Channel) Validate(tol float64) error` | Validate reports whether the transition matrix is a well-formed discrete memoryless channel. |

</details>

<details>
<summary><code>HuffmanNode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BuildHuffmanTree(weights map[string]float64) (*HuffmanNode, error)` | BuildHuffmanTree builds an optimal prefix-code tree for the symbols in weights, whose keys are symbols and values are non-negative weights… |
| `func (n *HuffmanNode) Codewords() map[string]string` | Codewords returns the binary codeword for every symbol in the tree rooted at n, assigning "0" to each left branch and "1" to each right branch. |
| `func (n *HuffmanNode) IsLeaf() bool` | IsLeaf reports whether the node is a leaf (has no children). |

</details>

### Variables

`ErrChannelRow`, `ErrChannelShape`, `ErrHamming`, `ErrHuffman`, `ErrNotNormalizable`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/infotheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
