package crypto

import (
	"errors"
	"math/big"
	"math/rand"
)

// DHParams holds the public Diffie-Hellman domain parameters shared by both
// parties: a prime modulus P and a generator G of a (large) subgroup of the
// multiplicative group modulo P.
type DHParams struct {
	P *big.Int // prime modulus
	G *big.Int // generator
}

// GenerateDHParams produces Diffie-Hellman domain parameters using a freshly
// generated safe prime P of the requested bit length and the generator G = 2,
// which generates the prime-order subgroup of a safe-prime group. Randomness is
// drawn from rng. bits must be at least 3.
func GenerateDHParams(bits int, rng *rand.Rand) DHParams {
	p := GenerateSafePrime(bits, rng)
	return DHParams{P: p, G: big.NewInt(2)}
}

// DHGeneratePrivate returns a random Diffie-Hellman private exponent in the
// range [2, P-2] drawn from rng. Each party keeps this value secret.
func (params DHParams) DHGeneratePrivate(rng *rand.Rand) *big.Int {
	lo := big.NewInt(2)
	hi := new(big.Int).Sub(params.P, cryptoTwo)
	return cryptoRandRange(rng, lo, hi)
}

// DHPublicKey computes the public value G^priv mod P that a party transmits to
// its peer.
func (params DHParams) DHPublicKey(priv *big.Int) *big.Int {
	return ModExp(params.G, priv, params.P)
}

// DHSharedSecret computes the shared secret peerPub^priv mod P from this
// party's private exponent and the peer's public value. Both parties arrive at
// the same secret. It returns an error if the peer's public value is out of the
// range [1, P-1].
func (params DHParams) DHSharedSecret(priv, peerPub *big.Int) (*big.Int, error) {
	if peerPub.Sign() <= 0 || peerPub.Cmp(params.P) >= 0 {
		return nil, errors.New("crypto: Diffie-Hellman peer public value out of range")
	}
	return ModExp(peerPub, priv, params.P), nil
}
