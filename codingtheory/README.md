# codingtheory — implements classical error-correcting codes

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/codingtheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/codingtheory)

Package codingtheory implements classical error-correcting codes, source codes
and the finite-field arithmetic they are built on.

- Galois-field arithmetic: carry-less polynomial arithmetic over GF(2)
(GF2Poly* helpers) and a table-driven extension field GF(2^m) exposed as the
Field type, with the usual add/mul/div/inv/pow/log operations and ready-made
constructors NewGF4 … NewGF256. - Block codes over GF(2): general linear
codes described by generator and parity-check matrices (LinearCode), Hamming
and extended-Hamming codes, the perfect binary Golay(23,12) and extended
Golay(24,12) codes, repetition and parity-check codes, and Hadamard/Walsh
codes. - Cyclic and polynomial codes: generator-polynomial encoding, syndrome
computation, BCH codes with bounded-distance decoding, and Reed-Solomon codes
over GF(2^m) with Berlekamp-Massey / Chien / Forney decoding. - Stream codes:
configurable cyclic redundancy checks (CRC of arbitrary width and polynomial)
and rate-1/n convolutional codes with hard and soft Viterbi decoding. -
Modulation and framing helpers: binary and multi-symbol Gray codes, block and
convolutional interleavers, and LDPC parity-check utilities. - Source coding:
Huffman coding and a fixed-point arithmetic coder, plus Hamming-weight /
Hamming-distance metrics.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/codingtheory"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ArithmeticDecode(bits []int, freqs []int, count int) ([]int, error)` | ArithmeticDecode decodes exactly count symbols from a bit slice produced by ArithmeticEncode under the same frequency model. |
| `func ArithmeticEncode(symbols []int, freqs []int) ([]int, error)` | ArithmeticEncode encodes a sequence of symbols (each an index into freqs) under the static frequency model freqs, returning the compressed bit slice. |
| `func BinaryToGray(n uint64) uint64` | BinaryToGray converts a non-negative integer to its reflected binary Gray code, in which consecutive values differ in exactly one bit. |
| `func BitsToPoly(bits []int) int` | BitsToPoly converts a little-endian bit slice (index i = coefficient of x^i) into a GF(2) polynomial packed as an int. |
| `func BitsToUint(bits []int) uint64` | BitsToUint packs a bit slice (most-significant bit first) into a uint64. |
| `func BlockDeinterleave(data []int, rows, cols int) ([]int, error)` | BlockDeinterleave inverts BlockInterleave, restoring the original row-major order. |
| `func BlockInterleave(data []int, rows, cols int) ([]int, error)` | BlockInterleave writes the input row by row into a rows-by-cols matrix and reads it out column by column, returning the interleaved sequence. |
| `func BuildHuffmanCodes(weights map[int]int) (map[int][]int, error)` | BuildHuffmanCodes is a convenience wrapper returning the code map directly from symbol weights. |
| `func CumulativeFreqs(freqs []int) []int` | CumulativeFreqs returns the cumulative-frequency table of a frequency slice: a slice of length len(freqs)+1 where entry i is the sum of freqs[:i]. |
| `func DotGF2(a, b []int) int` | DotGF2 returns the GF(2) inner product (parity of the pairwise ands) of two equal-length bit vectors. |
| `func ExtendedGolayDecode(g *GolayCode, word []int) (msg []int, numErrors int, err error)` | ExtendedGolayDecode decodes a 24-bit extended Golay code word. |
| `func ExtendedGolayEncode(g *GolayCode, msg []int) ([]int, error)` | ExtendedGolayEncode encodes a 12-bit message into a 24-bit extended Golay code word by appending an overall even-parity bit to the Golay(23,12) code… |
| `func ExtendedHamming74Decode(word []int) (data []int, doubleError bool, err error)` | ExtendedHamming74Decode decodes an eight-bit extended Hamming code word. |
| `func ExtendedHamming74Encode(data []int) ([]int, error)` | ExtendedHamming74Encode encodes four data bits into an eight-bit extended Hamming(8,4) SECDED code word by appending an overall even-parity bit to… |
| `func GF2PolyAdd(a, b int) int` | GF2PolyAdd returns the sum a+b of two GF(2) polynomials, which equals their bitwise exclusive-or. |
| `func GF2PolyDegree(p int) int` | GF2PolyDegree returns the degree of the GF(2) polynomial p, i.e. |
| `func GF2PolyDiv(a, b int) int` | GF2PolyDiv returns the quotient of the GF(2) polynomial a divided by the non-zero polynomial b. |
| `func GF2PolyDivMod(a, b int) (quotient, remainder int)` | GF2PolyDivMod divides the GF(2) polynomial a by the non-zero polynomial b and returns the quotient and remainder, so that a = quotient*b + remainder… |
| `func GF2PolyEval(p, x int) int` | GF2PolyEval evaluates the GF(2) polynomial p at the bit value x (0 or 1), returning 0 or 1. |
| `func GF2PolyGCD(a, b int) int` | GF2PolyGCD returns a greatest common divisor of the GF(2) polynomials a and b using the Euclidean algorithm. |
| `func GF2PolyMod(a, b int) int` | GF2PolyMod returns the remainder of the GF(2) polynomial a modulo the non-zero polynomial b. |
| `func GF2PolyMul(a, b int) int` | GF2PolyMul returns the carry-less product a*b of two GF(2) polynomials. |
| `func GF2PolyMulMod(a, b, m int) int` | GF2PolyMulMod returns (a*b) mod m for GF(2) polynomials, with m non-zero. |
| `func GF2PolyPowMod(base, exp, m int) int` | GF2PolyPowMod returns (base^exp) mod m for GF(2) polynomials using square-and-multiply. |
| `func GF2PolyString(p int) string` | GF2PolyString renders a GF(2) polynomial in descending powers of x, e.g. |
| `func GF2PolySub(a, b int) int` | GF2PolySub returns the difference a-b of two GF(2) polynomials. |
| `func GallagerLDPC(n, wc, wr int, seed int64) ([][]int, error)` | GallagerLDPC constructs a regular (wc, wr) Gallager parity-check matrix of length n. |
| `func GrayBitsToBinary(bits []int) uint64` | GrayBitsToBinary converts a fixed-width Gray-code bit slice (most-significant bit first) back to its integer value. |
| `func GrayCodeBits(value uint64, nbits int) []int` | GrayCodeBits returns the Gray code of value with a fixed width of nbits, most-significant bit first. |
| `func GrayCodeSequence(n int) []uint64` | GrayCodeSequence returns the length-2^n reflected Gray code sequence as a slice of integers, where consecutive entries differ in exactly one bit. |
| `func GrayCodeTable(n int) [][]int` | GrayCodeTable returns the length-2^n reflected Gray code sequence as bit slices of width n (most-significant bit first). |
| `func GrayToBinary(g uint64) uint64` | GrayToBinary converts a reflected binary Gray code back to the ordinary binary integer it encodes. |
| `func GrayToNary(gray []int, radix int) []int` | GrayToNary inverts NaryToGray, converting modular n-ary Gray digits (most-significant first) back to ordinary base-radix digits. |
| `func HadamardDecode(word []int) ([]int, error)` | HadamardDecode maximum-likelihood decodes a length-2^m Hadamard code word by mapping bits to +/-1 and taking the fast Walsh-Hadamard transform: the… |
| `func HadamardEncode(msg []int) ([]int, error)` | HadamardEncode encodes an m-bit message into a length-2^m binary Hadamard code word: the bit at position x equals the GF(2) inner product of the… |
| `func HadamardMatrix(k int) [][]int` | HadamardMatrix returns the 2^k-by-2^k Sylvester-Hadamard matrix with entries +1 and -1, built recursively from H_1 = [[1]] and the Kronecker doubling… |
| `func Hamming74Decode(word []int) (data, corrected []int, errPos int, err error)` | Hamming74Decode corrects up to one bit error in a received seven-bit word and returns the four recovered data bits along with the corrected code word… |
| `func Hamming74Encode(data []int) ([]int, error)` | Hamming74Encode encodes four data bits into a seven-bit systematic Hamming(7,4) code word. |
| `func Hamming74Generator() [][]int` | Hamming74Generator returns the 4-by-7 generator matrix of the systematic Hamming(7,4) code with layout [data \| parity]. |
| `func Hamming74ParityCheck() [][]int` | Hamming74ParityCheck returns the 3-by-7 parity-check matrix matching Hamming74Generator. |
| `func Hamming74Syndrome(word []int) ([]int, error)` | Hamming74Syndrome returns the three-bit syndrome of a received Hamming(7,4) word. |
| `func HammingDistance(a, b []int) int` | HammingDistance returns the number of positions in which a and b differ. |
| `func HammingDistanceUint(a, b uint64) int` | HammingDistanceUint returns the number of differing bits between a and b. |
| `func HammingWeight(v []int) int` | HammingWeight returns the number of non-zero entries of v. |
| `func HammingWeightUint(x uint64) int` | HammingWeightUint returns the population count of an unsigned integer's bits. |
| `func HuffmanAverageLength(codes map[int][]int, weights map[int]int) float64` | HuffmanAverageLength returns the expected code-word length in bits per symbol for the given code map and symbol weights. |
| `func HuffmanCodes(root *HuffmanNode) map[int][]int` | HuffmanCodes walks a Huffman tree and returns the code words as bit slices keyed by symbol id (left edge = 0, right edge = 1). |
| `func HuffmanDecode(bits []int, root *HuffmanNode) ([]int, error)` | HuffmanDecode decodes a bit slice back into symbols using a Huffman tree. |
| `func HuffmanEncode(data []int, codes map[int][]int) ([]int, error)` | HuffmanEncode encodes a symbol sequence into a bit slice using the given code map. |
| `func IdentityGF2(n int) [][]int` | IdentityGF2 returns the n-by-n identity matrix over GF(2). |
| `func Interleave(data []int, branches, delay int) ([]int, error)` | Interleave runs a fresh convolutional interleaver over an entire sequence and returns the emitted sequence of equal length. |
| `func InterleavePermutation(rows, cols int) []int` | InterleavePermutation returns the permutation applied by a rows-by-cols block interleaver, i.e. |
| `func InverseWalshHadamardTransform(a []float64) []float64` | InverseWalshHadamardTransform returns the inverse Walsh-Hadamard transform, i.e. |
| `func InvertPermutation(perm []int) []int` | InvertPermutation returns the inverse of a permutation. |
| `func IsPowerOfTwo(n int) bool` | IsPowerOfTwo reports whether n is a positive power of two. |
| `func MajorityVote(bits []int) int` | MajorityVote returns the majority bit of a bit slice; ties resolve to 0. |
| `func MatMulGF2(a, b [][]int) [][]int` | MatMulGF2 returns the product a*b of two GF(2) matrices. |
| `func MatVecGF2(m [][]int, v []int) []int` | MatVecGF2 multiplies the GF(2) matrix m (rows of equal length) by the column vector v, returning m*v over GF(2). |
| `func NaryToGray(digits []int, radix int) []int` | NaryToGray converts a value expressed as base-radix digits (most-significant digit first) to its modular n-ary Gray code digits: consecutive values… |
| `func ParityBit(bits []int) int` | ParityBit returns the even-parity bit of a bit slice (the exclusive-or of all entries), so that appending it makes the total parity even. |
| `func PermuteInts(data, perm []int) ([]int, error)` | PermuteInts applies a permutation to a slice, returning out`i`=data[perm`i`]. |
| `func PolyToBits(p, n int) []int` | PolyToBits converts a GF(2) polynomial packed as an int into a length-n bit slice in little-endian order, where index i holds the coefficient of x^i. |
| `func RankGF2(m [][]int) int` | RankGF2 returns the rank of a GF(2) matrix computed by Gaussian elimination. |
| `func RepetitionDecode(word []int, n int) ([]int, error)` | RepetitionDecode decodes a repetition-coded word by majority vote over each block of n symbols, recovering the original message. |
| `func RepetitionEncode(msg []int, n int) ([]int, error)` | RepetitionEncode repeats each message bit n times, producing a code word of length n*len(msg). |
| `func ShannonCodeLength(p float64) float64` | ShannonCodeLength returns the ideal (fractional) number of bits to encode a symbol of probability p under an optimal arithmetic coder, -log2(p). |
| `func SingleParityCheck(word []int) bool` | SingleParityCheck reports whether a word has even parity (no detected single error). |
| `func SingleParityDecode(word []int) (msg []int, errorDetected bool, err error)` | SingleParityDecode strips the trailing parity bit and reports whether the received word passed the even-parity check (errorDetected is true when an… |
| `func SingleParityEncode(msg []int) ([]int, error)` | SingleParityEncode appends an even-parity bit to the message, producing a length-(k+1) code word of even weight. |
| `func TransposeGF2(m [][]int) [][]int` | TransposeGF2 returns the transpose of a rectangular GF(2) matrix. |
| `func UintToBits(x uint64, n int) []int` | UintToBits unpacks the low n bits of x into a bit slice, most-significant bit first. |
| `func VecMatGF2(v []int, m [][]int) []int` | VecMatGF2 multiplies the row vector v by the GF(2) matrix m, returning v*m over GF(2). |
| `func WalshCode(k, index int) []int` | WalshCode returns the index-th Walsh code of length 2^k as a bit slice, taken from the rows of the Sylvester-Hadamard matrix (a +1 entry maps to bit… |
| `func WalshHadamardTransform(a []float64) []float64` | WalshHadamardTransform returns the natural-ordered fast Walsh-Hadamard transform of a, whose length must be a power of two. |
| `func WalshHadamardTransformInt(a []int) []int` | WalshHadamardTransformInt is the integer-valued forward Walsh-Hadamard transform for a slice whose length is a power of two. |
| `func XORVectors(a, b []int) []int` | XORVectors returns the entrywise exclusive-or of two equal-length bit vectors. |

### Types

| Type | What it is |
| --- | --- |
| `BCHCode` | BCHCode is a binary BCH code of length n = 2^m-1 designed to correct t errors. |
| `CRC` | CRC describes a parameterised cyclic redundancy check. |
| `ConvCode` | ConvCode is a rate-1/n binary convolutional code with constraint length K. The n generator polynomials are K-bit integers whose bit K-1 is the tap on… |
| `ConvolutionalInterleaver` | ConvolutionalInterleaver models a classic (branches, delay) convolutional interleaver: branch i delays its symbols by i*delay cells. |
| `CyclicCode` | CyclicCode is a binary cyclic [n,k] code generated by a polynomial g(x) that divides x^n-1 over GF(2). |
| `Field` | Field is a finite field GF(2^m) represented by log/antilog tables generated from a primitive polynomial. |
| `GolayCode` | GolayCode is the perfect binary Golay code together with a complete syndrome-to-error table for exact correction of up to three bit errors. |
| `HammingCode` | HammingCode is a general binary Hamming code of length n = 2^r-1 with dimension k = n-r and minimum distance 3, correcting any single bit error. |
| `HuffmanNode` | HuffmanNode is a node of a Huffman code tree. |
| `LDPCCode` | LDPCCode wraps a binary parity-check matrix H (m rows by n columns) over GF(2). |
| `LinearCode` | LinearCode is a binary linear [n,k] block code over GF(2) described by a k-by-n generator matrix G and an (n-k)-by-n parity-check matrix H with… |
| `ReedSolomon` | ReedSolomon is a Reed-Solomon code over GF(2^m). |

<details>
<summary><code>BCHCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBCHCode(field *Field, t int) (*BCHCode, error)` | NewBCHCode returns the narrow-sense binary BCH code over the given field that corrects up to t errors. |
| `func NewBCHCodeB(field *Field, t, b int) (*BCHCode, error)` | NewBCHCodeB returns a binary BCH code over the given field correcting up to t errors with first root exponent b (b=1 gives the narrow-sense code). |
| `func (c *BCHCode) Cyclic() *CyclicCode` | Cyclic returns the underlying cyclic code. |
| `func (c *BCHCode) Decode(received []int) (corrected []int, numErrors int, err error)` | Decode corrects up to t bit errors in a received length-n word, returning the corrected code word and the number of corrected bits. |
| `func (c *BCHCode) DecodeMessage(received []int) (msg []int, numErrors int, err error)` | DecodeMessage decodes a received word and returns the k message bits. |
| `func (c *BCHCode) DesignedDistance() int` | DesignedDistance returns the designed minimum distance 2t+1. |
| `func (c *BCHCode) Encode(msg []int) ([]int, error)` | Encode systematically encodes a length-k message into a length-n code word. |
| `func (c *BCHCode) Field() *Field` | Field returns the underlying field. |
| `func (c *BCHCode) Generator() int` | Generator returns the generator polynomial packed as a GF(2) polynomial int. |
| `func (c *BCHCode) K() int` | K returns the message dimension. |
| `func (c *BCHCode) N() int` | N returns the code length 2^m-1. |
| `func (c *BCHCode) Syndromes(received []int) []int` | Syndromes returns the 2t syndrome values S_j = R(alpha^(b+j)) for j = 0 … 2t-1 of a received length-n binary word. |
| `func (c *BCHCode) T() int` | T returns the designed error-correction capability. |

</details>

<details>
<summary><code>CRC</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CRC16CCITTFalse() CRC` | CRC16CCITTFalse returns the CRC-16/CCITT-FALSE model (poly 0x1021, init 0xFFFF, no reflection). |
| `func CRC16IBM() CRC` | CRC16IBM returns the CRC-16/IBM (ARC) model (poly 0x8005, reflected in and out). |
| `func CRC16XModem() CRC` | CRC16XModem returns the CRC-16/XMODEM model (poly 0x1021, zero init). |
| `func CRC32() CRC` | CRC32 returns the standard CRC-32 (as used by zlib/PNG): poly 0x04C11DB7, init and xorout 0xFFFFFFFF, reflected in and out. |
| `func CRC32C() CRC` | CRC32C returns the Castagnoli CRC-32C model (poly 0x1EDC6F41, reflected). |
| `func CRC8() CRC` | CRC8 returns the standard CRC-8 model (polynomial 0x07, zero init, no reflection). |
| `func (c CRC) Append(data []byte) []byte` | Append returns data with the Width/8 checksum bytes appended in big-endian order. |
| `func (c CRC) Checksum(data []byte) uint64` | Checksum computes the CRC of the byte slice data under the configured model. |
| `func (c CRC) ChecksumBits(bits []int) uint64` | ChecksumBits computes the CRC of a big-endian bit slice (each entry 0 or 1), processing the message bit by bit. |
| `func (c CRC) Verify(data []byte, expected uint64) bool` | Verify appends and checks: it reports whether the CRC of data equals the given expected checksum. |

</details>

<details>
<summary><code>ConvCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewConvCode(constraint int, polys []int) (*ConvCode, error)` | NewConvCode returns a rate-1/n convolutional code with the given constraint length and generator polynomials (each a K-bit integer). |
| `func NewConvCodeOctal(constraint int, octals []int) (*ConvCode, error)` | NewConvCodeOctal is a convenience constructor taking generator polynomials in octal (for example []int{7, 5} for the classic K=3 rate-1/2 code). |
| `func (c *ConvCode) Constraint() int` | Constraint returns the constraint length K. |
| `func (c *ConvCode) Encode(msg []int) ([]int, error)` | Encode encodes a message with zero-tail termination, flushing the encoder with K-1 trailing zero bits. |
| `func (c *ConvCode) EncodeUnterminated(msg []int) ([]int, error)` | EncodeUnterminated encodes a message without flushing, leaving the encoder in whatever state it ends in. |
| `func (c *ConvCode) HardViterbiDecode(received []int) ([]int, error)` | HardViterbiDecode decodes a zero-terminated received bit stream using the Viterbi algorithm with a Hamming-distance branch metric, returning the… |
| `func (c *ConvCode) N() int` | N returns the number of output bits per input bit (the inverse code rate). |
| `func (c *ConvCode) NumStates() int` | NumStates returns the number of trellis states 2^(K-1). |
| `func (c *ConvCode) Rate() float64` | Rate returns the nominal code rate 1/n. |
| `func (c *ConvCode) SoftViterbiDecode(received []float64) ([]int, error)` | SoftViterbiDecode decodes a zero-terminated stream of soft values using the Viterbi algorithm with a squared-Euclidean branch metric. |

</details>

<details>
<summary><code>ConvolutionalInterleaver</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewConvolutionalInterleaver(branches, delay int) (*ConvolutionalInterleaver, error)` | NewConvolutionalInterleaver returns a convolutional interleaver with the given number of branches and per-branch delay increment. |
| `func (ci *ConvolutionalInterleaver) Push(sym int) int` | Push feeds one symbol into the interleaver and returns the symbol emitted on the current branch (0 for symbols still filling an empty delay line). |

</details>

<details>
<summary><code>CyclicCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCyclicCode(n, gen int) (*CyclicCode, error)` | NewCyclicCode builds the cyclic code of length n generated by the polynomial gen (packed as an int). |
| `func (c *CyclicCode) DecodeMessage(word []int) ([]int, error)` | DecodeMessage extracts the k message bits from a systematic code word, reading the high-order coefficients. |
| `func (c *CyclicCode) Encode(msg []int) ([]int, error)` | Encode systematically encodes a length-k message bit slice into a length-n code word. |
| `func (c *CyclicCode) EncodeNonSystematic(msg []int) ([]int, error)` | EncodeNonSystematic encodes a message by the direct product c(x)=m(x)*g(x), producing a length-n (non-systematic) code word. |
| `func (c *CyclicCode) EncodePoly(msg int) int` | EncodePoly systematically encodes a message polynomial (packed as an int, of degree < k) into a code-word polynomial of degree < n. |
| `func (c *CyclicCode) Generator() int` | Generator returns the generator polynomial packed as an int. |
| `func (c *CyclicCode) GeneratorBits() []int` | GeneratorBits returns the generator polynomial as a length-(n-k+1) bit slice. |
| `func (c *CyclicCode) GeneratorMatrix() [][]int` | GeneratorMatrix returns the k-by-n systematic generator matrix of the code with rows in little-endian coefficient order. |
| `func (c *CyclicCode) IsCodeword(word []int) bool` | IsCodeword reports whether a length-n word is a valid code word. |
| `func (c *CyclicCode) K() int` | K returns the message dimension. |
| `func (c *CyclicCode) N() int` | N returns the code length. |
| `func (c *CyclicCode) Syndrome(word []int) ([]int, error)` | Syndrome returns the syndrome of a received length-n word as a bit slice of length n-k. |
| `func (c *CyclicCode) SyndromePoly(word int) int` | SyndromePoly returns the syndrome polynomial r(x) mod g(x) of a received word polynomial. |

</details>

<details>
<summary><code>Field</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewField(m, poly int) (*Field, error)` | NewField constructs GF(2^m) from the primitive polynomial poly, which must be given with its x^m bit set (for example x^4+x+1 is 0b10011 = 19). |
| `func NewGF128() *Field` | NewGF128 returns GF(2^7) built from x^7+x^3+1. |
| `func NewGF16() *Field` | NewGF16 returns GF(2^4) built from x^4+x+1. |
| `func NewGF2() *Field` | NewGF2 returns GF(2) itself using the trivial primitive polynomial x. |
| `func NewGF256() *Field` | NewGF256 returns GF(2^8) built from x^8+x^4+x^3+x^2+1 (the polynomial 0x11d used by QR codes and many Reed-Solomon systems). |
| `func NewGF32() *Field` | NewGF32 returns GF(2^5) built from x^5+x^2+1. |
| `func NewGF4() *Field` | NewGF4 returns GF(2^2) built from x^2+x+1. |
| `func NewGF64() *Field` | NewGF64 returns GF(2^6) built from x^6+x+1. |
| `func NewGF8() *Field` | NewGF8 returns GF(2^3) built from x^3+x+1. |
| `func (f *Field) Add(a, b int) int` | Add returns a+b in GF(2^m), which is the bitwise exclusive-or. |
| `func (f *Field) Alpha() int` | Alpha returns the primitive element alpha (numeric value 2) of the field. |
| `func (f *Field) BerlekampMassey(s []int) []int` | BerlekampMassey returns the error-locator polynomial (little-endian, with constant term one) implied by the syndrome sequence s over GF(2^m) using… |
| `func (f *Field) ChienSearch(lambda []int, n int) []int` | ChienSearch returns the little-endian error positions (powers of x in [0,n)) for which the locator polynomial lambda evaluates to zero at alpha^{-i}. |
| `func (f *Field) Div(a, b int) int` | Div returns a/b in GF(2^m). |
| `func (f *Field) Elements() []int` | Elements returns a fresh slice of all 2^m field elements in ascending order. |
| `func (f *Field) Exp(e int) int` | Exp returns alpha^e, the e-th power of the field's primitive element. |
| `func (f *Field) Inv(a int) int` | Inv returns the multiplicative inverse of a in GF(2^m). |
| `func (f *Field) Log(a int) int` | Log returns the discrete logarithm of a non-zero element a to base alpha, i.e. |
| `func (f *Field) M() int` | M returns the extension degree m of GF(2^m). |
| `func (f *Field) MinimalPoly(e int) int` | MinimalPoly returns the minimal polynomial over GF(2) of the element alpha^e, packed as a GF(2) polynomial int (bit i is the coefficient of x^i). |
| `func (f *Field) Mul(a, b int) int` | Mul returns a*b in GF(2^m) via log/antilog tables. |
| `func (f *Field) Order() int` | Order returns the multiplicative order 2^m-1 of the field's cyclic group. |
| `func (f *Field) Poly() int` | Poly returns the primitive polynomial (with the x^m bit set) defining the field. |
| `func (f *Field) PolyAdd(a, b []int) []int` | PolyAdd returns the sum of two polynomials over GF(2^m). |
| `func (f *Field) PolyDegree(p []int) int` | PolyDegree returns the degree of a polynomial over GF(2^m). |
| `func (f *Field) PolyDeriv(p []int) []int` | PolyDeriv returns the formal derivative of p over GF(2^m). |
| `func (f *Field) PolyDivMod(a, b []int) (quotient, remainder []int)` | PolyDivMod divides a by the non-zero polynomial b over GF(2^m), returning the quotient and remainder with deg(remainder) < deg(b). |
| `func (f *Field) PolyEqual(a, b []int) bool` | PolyEqual reports whether a and b are the same polynomial, ignoring trailing zero coefficients. |
| `func (f *Field) PolyEval(p []int, x int) int` | PolyEval evaluates the polynomial p at the field element x using Horner's rule. |
| `func (f *Field) PolyMod(a, b []int) []int` | PolyMod returns a modulo the non-zero polynomial b over GF(2^m). |
| `func (f *Field) PolyMonic(p []int) []int` | PolyMonic returns p scaled so its leading coefficient is one. |
| `func (f *Field) PolyMul(a, b []int) []int` | PolyMul returns the product of two polynomials over GF(2^m). |
| `func (f *Field) PolyScale(p []int, s int) []int` | PolyScale multiplies every coefficient of p by the scalar s in GF(2^m). |
| `func (f *Field) PolySub(a, b []int) []int` | PolySub returns the difference of two polynomials; over GF(2^m) this equals PolyAdd. |
| `func (f *Field) Pow(a, e int) int` | Pow returns a raised to the integer power e in GF(2^m); e may be negative. |
| `func (f *Field) Size() int` | Size returns the number of elements 2^m in the field. |
| `func (f *Field) Sub(a, b int) int` | Sub returns a-b in GF(2^m); over characteristic two this equals Add. |
| `func (f *Field) Trace(a int) int` | Trace returns the field trace Tr(a) = a + a^2 + a^4 + ... |

</details>

<details>
<summary><code>GolayCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGolayCode() *GolayCode` | NewGolayCode returns the perfect binary Golay(23,12) code with its decoding table precomputed. |
| `func (g *GolayCode) Cyclic() *CyclicCode` | Cyclic returns the underlying cyclic code. |
| `func (g *GolayCode) Decode(word []int) (corrected []int, numErrors int, err error)` | Decode corrects up to three bit errors in a received 23-bit word and returns the corrected code word and the number of errors corrected. |
| `func (g *GolayCode) DecodeMessage(word []int) (msg []int, numErrors int, err error)` | DecodeMessage decodes a received 23-bit word and returns the 12 message bits. |
| `func (g *GolayCode) Encode(msg []int) ([]int, error)` | Encode systematically encodes a 12-bit message into a 23-bit Golay code word. |
| `func (g *GolayCode) K() int` | K returns the message dimension 12. |
| `func (g *GolayCode) MinimumDistance() int` | MinimumDistance returns the minimum distance 7 of the Golay(23,12) code. |
| `func (g *GolayCode) N() int` | N returns the code length 23. |
| `func (g *GolayCode) Syndrome(word []int) ([]int, error)` | Syndrome returns the 11-bit syndrome of a received 23-bit word. |

</details>

<details>
<summary><code>HammingCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewHammingCode(r int) (*HammingCode, error)` | NewHammingCode returns the Hamming code with r parity bits (r >= 2), of length 2^r-1 and dimension 2^r-1-r. |
| `func (c *HammingCode) Decode(word []int) (corrected []int, errPos int, err error)` | Decode corrects up to one bit error in a received length-n word and returns the corrected code word together with the error position (-1 if none). |
| `func (c *HammingCode) K() int` | K returns the message dimension 2^r-1-r. |
| `func (c *HammingCode) LinearCode() (*LinearCode, error)` | LinearCode converts the Hamming code into a general LinearCode built from its parity-check matrix. |
| `func (c *HammingCode) N() int` | N returns the code length 2^r-1. |
| `func (c *HammingCode) ParityCheck() [][]int` | ParityCheck returns a copy of the r-by-n parity-check matrix. |
| `func (c *HammingCode) R() int` | R returns the number of parity bits r. |
| `func (c *HammingCode) Syndrome(word []int) ([]int, error)` | Syndrome returns the r-bit syndrome H*wᵀ of a received length-n word. |
| `func (c *HammingCode) SyndromePosition(word []int) (int, error)` | SyndromePosition returns the 0-based position of a single bit error implied by a received word, or -1 if the syndrome is zero. |

</details>

<details>
<summary><code>HuffmanNode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BuildHuffmanTree(weights map[int]int) (*HuffmanNode, error)` | BuildHuffmanTree builds a Huffman tree for the given symbol weights (symbol id -> positive weight). |
| `func (n *HuffmanNode) IsLeaf() bool` | IsLeaf reports whether the node is a leaf. |

</details>

<details>
<summary><code>LDPCCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewLDPCCode(h [][]int) (*LDPCCode, error)` | NewLDPCCode builds an LDPC code from a rectangular binary parity-check matrix. |
| `func (c *LDPCCode) BitFlipDecode(received []int, maxIter int) (decoded []int, iters int, ok bool, err error)` | BitFlipDecode runs Gallager's hard-decision bit-flipping algorithm for up to maxIter iterations, flipping in each round every bit that fails a… |
| `func (c *LDPCCode) CheckDegrees() []int` | CheckDegrees returns the row weights (check-node degrees) of H. |
| `func (c *LDPCCode) Density() float64` | Density returns the fraction of ones in the parity-check matrix. |
| `func (c *LDPCCode) IsCodeword(word []int) bool` | IsCodeword reports whether a length-n word satisfies all parity checks. |
| `func (c *LDPCCode) IsRegular() bool` | IsRegular reports whether every variable node has the same degree and every check node has the same degree. |
| `func (c *LDPCCode) M() int` | M returns the number of parity checks (check nodes). |
| `func (c *LDPCCode) N() int` | N returns the code length (number of variable nodes). |
| `func (c *LDPCCode) ParityCheck() [][]int` | ParityCheck returns a copy of the parity-check matrix. |
| `func (c *LDPCCode) Syndrome(word []int) ([]int, error)` | Syndrome returns the length-m syndrome H*wᵀ of a received length-n word. |
| `func (c *LDPCCode) VariableDegrees() []int` | VariableDegrees returns the column weights (variable-node degrees) of H. |

</details>

<details>
<summary><code>LinearCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewLinearCode(g [][]int) (*LinearCode, error)` | NewLinearCode builds a linear code from a k-by-n generator matrix G. G must be non-empty, rectangular and full row rank; the parity-check matrix is… |
| `func NewLinearCodeGH(g, h [][]int) (*LinearCode, error)` | NewLinearCodeGH builds a linear code from explicit generator and parity-check matrices. |
| `func (c *LinearCode) Codewords() [][]int` | Codewords returns all 2^k code words in message-counting order. |
| `func (c *LinearCode) DecodeSyndrome(word []int, table map[uint64][]int) ([]int, error)` | DecodeSyndrome performs syndrome (coset-leader) decoding of a received word using a table built by SyndromeTable, returning the most likely code word. |
| `func (c *LinearCode) Encode(msg []int) ([]int, error)` | Encode encodes a length-k message into a length-n code word as m*G. It returns ErrLength if the message length is not k. |
| `func (c *LinearCode) Generator() [][]int` | Generator returns a copy of the generator matrix. |
| `func (c *LinearCode) IsCodeword(word []int) bool` | IsCodeword reports whether a length-n word is a valid code word (zero syndrome). |
| `func (c *LinearCode) K() int` | K returns the message dimension k. |
| `func (c *LinearCode) MinimumDistance() int` | MinimumDistance returns the minimum Hamming distance of the code, computed as the minimum weight over all non-zero code words. |
| `func (c *LinearCode) N() int` | N returns the code length n. |
| `func (c *LinearCode) ParityCheck() [][]int` | ParityCheck returns a copy of the parity-check matrix. |
| `func (c *LinearCode) Rate() float64` | Rate returns the code rate k/n. |
| `func (c *LinearCode) Redundancy() int` | Redundancy returns the number of parity checks n-k. |
| `func (c *LinearCode) Syndrome(word []int) ([]int, error)` | Syndrome returns the syndrome H*rᵀ of a received length-n word. |
| `func (c *LinearCode) SyndromeTable() map[uint64][]int` | SyndromeTable builds a standard-array coset-leader table mapping each syndrome (packed as a big-endian integer) to a minimum-weight error pattern. |

</details>

<details>
<summary><code>ReedSolomon</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewReedSolomon(field *Field, nsym int) (*ReedSolomon, error)` | NewReedSolomon returns a Reed-Solomon code over the given field with nsym parity symbols and first consecutive root exponent fcr=0. |
| `func NewReedSolomonFCR(field *Field, nsym, fcr int) (*ReedSolomon, error)` | NewReedSolomonFCR returns a Reed-Solomon code with nsym parity symbols and a configurable first-consecutive-root exponent fcr. |
| `func (rs *ReedSolomon) Correction() int` | Correction returns the maximum number of symbol errors the code can correct, nsym/2. |
| `func (rs *ReedSolomon) Decode(received []int) (corrected []int, numErrors int, err error)` | Decode corrects up to nsym/2 symbol errors in a received length-n word, returning the corrected code word and the number of corrected symbols. |
| `func (rs *ReedSolomon) DecodeMessage(received []int) (msg []int, numErrors int, err error)` | DecodeMessage decodes a received word and returns the k recovered message symbols. |
| `func (rs *ReedSolomon) Encode(msg []int) ([]int, error)` | Encode systematically encodes a length-k message into a length-n code word. |
| `func (rs *ReedSolomon) FCR() int` | FCR returns the first-consecutive-root exponent. |
| `func (rs *ReedSolomon) Field() *Field` | Field returns the underlying field. |
| `func (rs *ReedSolomon) GeneratorPoly() []int` | GeneratorPoly returns the generator polynomial as a big-endian symbol slice (index 0 is the highest-degree coefficient). |
| `func (rs *ReedSolomon) IsCodeword(word []int) bool` | IsCodeword reports whether a received length-n word has all-zero syndromes. |
| `func (rs *ReedSolomon) K() int` | K returns the message dimension n-nsym. |
| `func (rs *ReedSolomon) N() int` | N returns the code length 2^m-1. |
| `func (rs *ReedSolomon) NSym() int` | NSym returns the number of parity symbols. |
| `func (rs *ReedSolomon) Syndromes(received []int) []int` | Syndromes returns the nsym syndrome values S_j = R(alpha^(fcr+j)) of a received length-n word, in order j = 0 … nsym-1. |

</details>

### Constants

`Golay23Generator`

### Variables

`ErrBit`, `ErrFieldParam`, `ErrGenerator`, `ErrHuffman`, `ErrLength`, `ErrNotPrimitive`, `ErrRank`, `ErrShape`, `ErrTooManyErrors`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/codingtheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
