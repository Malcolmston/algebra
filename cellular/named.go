package cellular

// This file collects convenience constructors for the most-studied elementary
// rules so that callers can refer to them by name.

// Rule30 returns elementary rule 30, a chaotic, additive-looking rule used as a
// random-number source.
func Rule30() ElementaryRule { return ElementaryRule(30) }

// Rule54 returns elementary rule 54, a class-4 candidate with particle-like
// structures.
func Rule54() ElementaryRule { return ElementaryRule(54) }

// Rule60 returns elementary rule 60, an additive (Pascal-triangle mod 2) rule.
func Rule60() ElementaryRule { return ElementaryRule(60) }

// Rule90 returns elementary rule 90, the additive XOR rule that draws the
// Sierpinski triangle from a single seed.
func Rule90() ElementaryRule { return ElementaryRule(90) }

// Rule102 returns elementary rule 102, an additive rule mirror-equivalent to
// rule 60.
func Rule102() ElementaryRule { return ElementaryRule(102) }

// Rule110 returns elementary rule 110, proved capable of universal computation.
func Rule110() ElementaryRule { return ElementaryRule(110) }

// Rule150 returns elementary rule 150, the additive sum-mod-2 rule of all three
// neighbours.
func Rule150() ElementaryRule { return ElementaryRule(150) }

// Rule184 returns elementary rule 184, the traffic/ballistic-annihilation rule.
func Rule184() ElementaryRule { return ElementaryRule(184) }

// Rule18 returns elementary rule 18, a chaotic rule generating a Sierpinski-like
// fractal.
func Rule18() ElementaryRule { return ElementaryRule(18) }

// Rule22 returns elementary rule 22, a chaotic outer-totalistic rule.
func Rule22() ElementaryRule { return ElementaryRule(22) }

// Rule45 returns elementary rule 45, a chaotic rule often used for random
// generation.
func Rule45() ElementaryRule { return ElementaryRule(45) }

// Rule126 returns elementary rule 126, a chaotic rule.
func Rule126() ElementaryRule { return ElementaryRule(126) }

// Rule250 returns elementary rule 250, which fills space from a single seed.
func Rule250() ElementaryRule { return ElementaryRule(250) }

// Rule4 returns elementary rule 4, a class-2 rule whose isolated cells persist.
func Rule4() ElementaryRule { return ElementaryRule(4) }

// AdditiveElementaryRules returns the eight additive (GF(2)-linear) elementary
// rules in increasing order.
func AdditiveElementaryRules() []ElementaryRule {
	var out []ElementaryRule
	for i := 0; i < 256; i++ {
		r := ElementaryRule(i)
		if r.IsAdditive() {
			out = append(out, r)
		}
	}
	return out
}

// TotalisticElementaryRules returns the elementary rules that are totalistic in
// the strict sense (output depends only on the neighbourhood sum).
func TotalisticElementaryRules() []ElementaryRule {
	var out []ElementaryRule
	for i := 0; i < 256; i++ {
		r := ElementaryRule(i)
		if r.IsTotalistic() {
			out = append(out, r)
		}
	}
	return out
}

// SierpinskiRule is an alias for Rule90, whose single-seed evolution is the
// Sierpinski triangle.
func SierpinskiRule() ElementaryRule { return Rule90() }

// TrafficRule is an alias for Rule184, the elementary traffic model.
func TrafficRule() ElementaryRule { return Rule184() }
