package crypto

import (
	"errors"
	"math/big"
	"math/rand"
)

// RSAPublicKey is an RSA public key: the modulus N = p*q and the public
// exponent E. Encryption and signature verification use only these values.
type RSAPublicKey struct {
	N *big.Int // modulus, the product of two primes
	E *big.Int // public exponent, coprime to φ(N)
}

// RSAPrivateKey is an RSA private key. It embeds the matching public key and
// stores the private exponent D together with the secret primes P and Q and
// precomputed Chinese-Remainder-Theorem values (Dp, Dq, Qinv) that accelerate
// decryption via RSADecryptCRT.
type RSAPrivateKey struct {
	Public RSAPublicKey
	D      *big.Int // private exponent, the inverse of E modulo λ(N)
	P      *big.Int // first prime factor of N
	Q      *big.Int // second prime factor of N
	Dp     *big.Int // D mod (P-1)
	Dq     *big.Int // D mod (Q-1)
	Qinv   *big.Int // Q^-1 mod P
}

// GenerateRSAKey generates an RSA private key whose modulus has approximately
// the requested bit length, using the public exponent e (65537 is standard and
// used when e is nil) and drawing randomness from rng. bits must be at least 4
// so two distinct primes of bits/2 each can be produced. It returns an error if
// a suitable key cannot be formed with the chosen exponent (for example if e
// shares a factor with λ(N)); callers may simply retry. The determinism of rng
// makes generated keys reproducible.
func GenerateRSAKey(bits int, e *big.Int, rng *rand.Rand) (*RSAPrivateKey, error) {
	if bits < 4 {
		panic("crypto: GenerateRSAKey requires bits >= 4")
	}
	if e == nil {
		e = big.NewInt(65537)
	}
	if e.Sign() <= 0 || e.Bit(0) == 0 {
		return nil, errors.New("crypto: RSA public exponent must be a positive odd integer")
	}
	half := bits / 2
	for attempts := 0; attempts < 100; attempts++ {
		p := RandomPrime(half, rng)
		q := RandomPrime(bits-half, rng)
		if p.Cmp(q) == 0 {
			continue
		}
		n := new(big.Int).Mul(p, q)
		if n.BitLen() != bits {
			continue
		}
		// λ(N) = lcm(p-1, q-1).
		pm1 := new(big.Int).Sub(p, cryptoOne)
		qm1 := new(big.Int).Sub(q, cryptoOne)
		lambda := LCM(pm1, qm1)
		if GCD(e, lambda).Cmp(cryptoOne) != 0 {
			continue
		}
		d := new(big.Int).ModInverse(e, lambda)
		if d == nil {
			continue
		}
		key := &RSAPrivateKey{
			Public: RSAPublicKey{N: n, E: new(big.Int).Set(e)},
			D:      d,
			P:      p,
			Q:      q,
			Dp:     new(big.Int).Mod(d, pm1),
			Dq:     new(big.Int).Mod(d, qm1),
			Qinv:   new(big.Int).ModInverse(q, p),
		}
		return key, nil
	}
	return nil, errors.New("crypto: GenerateRSAKey exhausted attempts")
}

// PublicKey returns the public half of the private key.
func (k *RSAPrivateKey) PublicKey() RSAPublicKey {
	return k.Public
}

// RSAEncrypt performs the raw (textbook) RSA encryption primitive c = m^E mod
// N. The message m, interpreted as an integer, must satisfy 0 <= m < N;
// otherwise an error is returned. This is the pure mathematical operation with
// no padding and must not be used to encrypt real messages directly.
func RSAEncrypt(pub RSAPublicKey, m *big.Int) (*big.Int, error) {
	if m.Sign() < 0 || m.Cmp(pub.N) >= 0 {
		return nil, errors.New("crypto: RSA message out of range [0, N)")
	}
	return ModExp(m, pub.E, pub.N), nil
}

// RSADecrypt performs the raw RSA decryption primitive m = c^D mod N using the
// private exponent directly. The ciphertext c must satisfy 0 <= c < N.
func RSADecrypt(priv *RSAPrivateKey, c *big.Int) (*big.Int, error) {
	if c.Sign() < 0 || c.Cmp(priv.Public.N) >= 0 {
		return nil, errors.New("crypto: RSA ciphertext out of range [0, N)")
	}
	return ModExp(c, priv.D, priv.Public.N), nil
}

// RSADecryptCRT performs RSA decryption using the Chinese Remainder Theorem,
// exponentiating separately modulo P and Q and recombining. It is roughly four
// times faster than RSADecrypt for the same key and yields an identical result.
// The ciphertext c must satisfy 0 <= c < N.
func RSADecryptCRT(priv *RSAPrivateKey, c *big.Int) (*big.Int, error) {
	if c.Sign() < 0 || c.Cmp(priv.Public.N) >= 0 {
		return nil, errors.New("crypto: RSA ciphertext out of range [0, N)")
	}
	m1 := ModExp(c, priv.Dp, priv.P)
	m2 := ModExp(c, priv.Dq, priv.Q)
	// h = Qinv * (m1 - m2) mod P
	h := new(big.Int).Sub(m1, m2)
	h.Mul(h, priv.Qinv)
	h.Mod(h, priv.P)
	// m = m2 + h*Q
	m := new(big.Int).Mul(h, priv.Q)
	m.Add(m, m2)
	return m, nil
}

// RSASign produces a raw RSA signature s = m^D mod N over the integer message m
// (0 <= m < N). It is the textbook signing primitive without hashing or
// padding; RSAVerify inverts it.
func RSASign(priv *RSAPrivateKey, m *big.Int) (*big.Int, error) {
	if m.Sign() < 0 || m.Cmp(priv.Public.N) >= 0 {
		return nil, errors.New("crypto: RSA message out of range [0, N)")
	}
	return ModExp(m, priv.D, priv.Public.N), nil
}

// RSAVerify checks a raw RSA signature by computing s^E mod N and reporting
// whether it equals the original message m. It returns true exactly when the
// signature is valid for m under the given public key.
func RSAVerify(pub RSAPublicKey, m, s *big.Int) bool {
	if m.Sign() < 0 || m.Cmp(pub.N) >= 0 {
		return false
	}
	recovered := ModExp(s, pub.E, pub.N)
	return recovered.Cmp(m) == 0
}
