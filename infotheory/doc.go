// Package infotheory implements classical information theory: entropy
// measures, divergences, discrete channel models and error-control codes.
//
// The package is organised around discrete probability distributions
// represented as plain []float64 slices whose entries are non-negative and
// sum to one, and around joint distributions represented as row-major
// [][]float64 matrices whose entries sum to one. Unless stated otherwise all
// entropic quantities are expressed in bits (base-2 logarithms); companion
// routines suffixed with "Nat" report the same quantity in nats (natural
// logarithms), and BitsToNats / NatsToBits convert between the two.
//
// The following families are covered:
//
//   - entropy: Entropy (Shannon), EntropyBase, EntropyNat, BinaryEntropy,
//     Surprisal, Perplexity, NormalizedEntropy, Redundancy, GiniImpurity and
//     the Renyi family (RenyiEntropy, CollisionEntropy, MinEntropy,
//     HartleyEntropy, TsallisEntropy);
//   - joint and relative measures: JointEntropy, ConditionalEntropyYgivenX,
//     ConditionalEntropyXgivenY, MutualInformation, NormalizedMutualInformation,
//     VariationOfInformation, KLDivergence, CrossEntropy,
//     JensenShannonDivergence and JensenShannonDistance;
//   - channels: the binary symmetric channel (BSC), binary erasure channel
//     (BEC), the general discrete memoryless Channel, and the Blahut-Arimoto
//     algorithm for numerically computing channel capacity;
//   - source coding: Huffman coding (BuildHuffmanTree, HuffmanCodes,
//     HuffmanEncode, HuffmanDecode) and ShannonFanoCodes;
//   - error-control coding: Hamming(7,4) encode/decode with single-error
//     correction, Hamming weight and distance helpers;
//   - inequalities and rates: the Kraft-McMillan inequality (KraftSum,
//     KraftInequality, KraftEquality, McMillanInequality), the Fano
//     inequality bound and code-rate helpers.
//
// All routines are deterministic and depend only on the Go standard library.
package infotheory
