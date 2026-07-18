package crypto

import (
	"errors"
	"math/big"
	"math/rand"
)

// ElGamalPublicKey is an ElGamal public key over the multiplicative group
// modulo the prime P with generator G, where Y = G^x mod P for the secret x.
type ElGamalPublicKey struct {
	P *big.Int // prime modulus
	G *big.Int // generator
	Y *big.Int // public value G^x mod P
}

// ElGamalPrivateKey is an ElGamal private key: the secret exponent X together
// with the matching public key.
type ElGamalPrivateKey struct {
	Public ElGamalPublicKey
	X      *big.Int // secret exponent
}

// GenerateElGamalKey builds an ElGamal key pair over a fresh safe-prime group
// of the requested bit length, using generator G = 2 and a random secret
// exponent drawn from rng. bits must be at least 3.
func GenerateElGamalKey(bits int, rng *rand.Rand) *ElGamalPrivateKey {
	params := GenerateDHParams(bits, rng)
	x := params.DHGeneratePrivate(rng)
	y := ModExp(params.G, x, params.P)
	return &ElGamalPrivateKey{
		Public: ElGamalPublicKey{P: params.P, G: params.G, Y: y},
		X:      x,
	}
}

// ElGamalEncrypt encrypts the integer message m (0 <= m < P) under the public
// key, using the ephemeral randomness k drawn from rng. It returns the
// ciphertext pair (c1, c2) = (G^k, m·Y^k) mod P. Because k is random the same
// message encrypts to different ciphertexts on each call (unless rng is
// re-seeded identically). It returns an error if m is out of range.
func ElGamalEncrypt(pub ElGamalPublicKey, m *big.Int, rng *rand.Rand) (c1, c2 *big.Int, err error) {
	if m.Sign() < 0 || m.Cmp(pub.P) >= 0 {
		return nil, nil, errors.New("crypto: ElGamal message out of range [0, P)")
	}
	lo := big.NewInt(2)
	hi := new(big.Int).Sub(pub.P, cryptoTwo)
	k := cryptoRandRange(rng, lo, hi)
	c1 = ModExp(pub.G, k, pub.P)
	s := ModExp(pub.Y, k, pub.P)
	c2 = new(big.Int).Mul(m, s)
	c2.Mod(c2, pub.P)
	return c1, c2, nil
}

// ElGamalDecrypt recovers the message m from a ciphertext pair (c1, c2) using
// the private key, computing m = c2 · (c1^x)^-1 mod P. It returns an error if
// c1 is not invertible modulo P.
func ElGamalDecrypt(priv *ElGamalPrivateKey, c1, c2 *big.Int) (*big.Int, error) {
	p := priv.Public.P
	s := ModExp(c1, priv.X, p)
	sInv := new(big.Int).ModInverse(s, p)
	if sInv == nil {
		return nil, errors.New("crypto: ElGamal shared value not invertible")
	}
	m := new(big.Int).Mul(c2, sInv)
	m.Mod(m, p)
	return m, nil
}
