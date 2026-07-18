package algebra

import (
	"math/big"
	"sort"
	"strconv"
	"strings"
)

// This file adds presentation-quality rendering for expression trees: LaTeX
// (LaTeX and LaTeXEq), a two-dimensional Unicode pretty printer (Pretty) and
// presentation MathML (MathML). All three are pure functions that only
// allocate, never mutate their input, and depend on the standard library only.
// Because they live in package algebra they type-switch directly on the
// unexported node types (*sum, *product, *power, *fn, *fn2, *integral, *bigOp)
// and the exported numeric/atomic types.
//
// The parenthesization and term-ordering logic mirrors print.go so that all
// renderers agree with String() on precedence.

// presHalf is the exact rational 1/2, recognised as a square root.
var presHalf = big.NewRat(1, 2)

// presIsHalf reports whether e is the exact rational 1/2.
func presIsHalf(e Expr) bool {
	r, ok := e.(*Rational)
	return ok && r.Val.Cmp(presHalf) == 0
}

// presSortedSumTerms returns the terms of a sum in the same descending
// polynomial-degree order used by (*sum).String, so every renderer agrees on
// term order.
func presSortedSumTerms(s *sum) []Expr {
	terms := append([]Expr(nil), s.args...)
	varSet := map[string]bool{}
	exps := make([]map[string]int, len(terms))
	for i, t := range terms {
		exps[i] = termExps(t)
		for v := range exps[i] {
			varSet[v] = true
		}
	}
	vars := make([]string, 0, len(varSet))
	for v := range varSet {
		vars = append(vars, v)
	}
	sort.Strings(vars)
	idx := map[Expr]int{}
	for i, t := range terms {
		idx[t] = i
	}
	sort.SliceStable(terms, func(i, j int) bool {
		di, dj := polyDegree(terms[i]), polyDegree(terms[j])
		if di != dj {
			return di > dj
		}
		ei, ej := exps[idx[terms[i]]], exps[idx[terms[j]]]
		for _, v := range vars {
			if ei[v] != ej[v] {
				return ei[v] > ej[v]
			}
		}
		return compareExpr(terms[i], terms[j]) < 0
	})
	return terms
}

// =========================================================================
// LaTeX
// =========================================================================

// LaTeX renders e as a LaTeX math fragment suitable for inclusion between
// dollar signs or in an equation environment.
//
// Rationals and negative or reciprocal powers become \frac{..}{..}; square
// roots and half powers become \sqrt{..}; every supported elementary function
// maps to its conventional LaTeX command (\sin, \cos, \Gamma, \psi, ...);
// the constants pi, E, I and oo render as \pi, e, i and \infty; a \cdot is
// inserted only where two factors would otherwise run together, and it is
// suppressed between a numeric coefficient and a symbol. Auto-sized
// \left(\right) delimiters wrap sums that appear under products or powers.
// The output is deterministic.
func LaTeX(e Expr) string { return ltxOf(e) }

// LaTeXEq renders the equation lhs = rhs as a LaTeX fragment, joining the two
// sides with " = ".
func LaTeXEq(lhs, rhs Expr) string { return ltxOf(lhs) + " = " + ltxOf(rhs) }

// ltxFnCmd maps single-argument function names to the LaTeX command that
// renders them as name\left(arg\right).
var ltxFnCmd = map[string]string{
	"sin": "\\sin", "cos": "\\cos", "tan": "\\tan",
	"cot": "\\cot", "sec": "\\sec", "csc": "\\csc",
	"sinh": "\\sinh", "cosh": "\\cosh", "tanh": "\\tanh", "coth": "\\coth",
	"asin": "\\arcsin", "acos": "\\arccos", "atan": "\\arctan",
	"exp": "\\exp", "log": "\\log", "ln": "\\ln", "arg": "\\arg",
	"gamma": "\\Gamma", "digamma": "\\psi",
}

// ltxOf renders e as LaTeX without any enclosing delimiters.
func ltxOf(e Expr) string {
	switch x := e.(type) {
	case *Integer:
		return x.Val.String()
	case *Rational:
		n, d := x.Val.Num(), x.Val.Denom()
		if n.Sign() < 0 {
			return "-\\frac{" + new(big.Int).Neg(n).String() + "}{" + d.String() + "}"
		}
		return "\\frac{" + n.String() + "}{" + d.String() + "}"
	case *Float:
		return strconv.FormatFloat(x.Val, 'g', -1, 64)
	case *Symbol:
		return x.Name
	case *Constant:
		return ltxConst(x)
	case *sum:
		return ltxSum(x)
	case *product:
		return ltxProduct(x)
	case *power:
		return ltxPow(x)
	case *fn:
		return ltxFn(x)
	case *fn2:
		return ltxFn2(x)
	case *integral:
		return "\\int " + ltxOf(x.arg) + " \\, d " + ltxOf(x.v)
	case *bigOp:
		return ltxBigOp(x)
	}
	return e.String()
}

// ltxConst renders a named constant in LaTeX.
func ltxConst(c *Constant) string {
	switch c.Name {
	case "pi":
		return "\\pi"
	case "E":
		return "e"
	case "I":
		return "i"
	case "oo":
		return "\\infty"
	case "-oo":
		return "-\\infty"
	}
	return c.Name
}

// ltxWrap parenthesizes e with auto-sized delimiters when it would need
// parentheses as the base or exponent of a power (mirroring needParen).
func ltxWrap(e Expr) string {
	if needParen(e) {
		return "\\left(" + ltxOf(e) + "\\right)"
	}
	return ltxOf(e)
}

// ltxMulFactor renders a product factor, wrapping only embedded sums.
func ltxMulFactor(f Expr) string {
	s := ltxOf(f)
	if _, ok := f.(*sum); ok {
		return "\\left(" + s + "\\right)"
	}
	return s
}

// ltxSum renders a sum, ordering terms like String and using + and -.
func ltxSum(s *sum) string {
	terms := presSortedSumTerms(s)
	var b strings.Builder
	for i, t := range terms {
		neg, mag := splitSign(t)
		ms := ltxOf(mag)
		if _, ok := mag.(*sum); ok {
			ms = "\\left(" + ms + "\\right)"
		}
		switch {
		case i == 0 && neg:
			b.WriteString("-" + ms)
		case i == 0:
			b.WriteString(ms)
		case neg:
			b.WriteString(" - " + ms)
		default:
			b.WriteString(" + " + ms)
		}
	}
	return b.String()
}

// ltxProduct renders a product, folding a leading -1 into a unary minus,
// collecting negative powers and rational coefficients into a \frac, and
// inserting \cdot only where adjacent factors would run together.
func ltxProduct(p *product) string {
	factors := p.factors
	neg := false
	if len(factors) > 1 && isMinusOne(factors[0]) {
		neg = true
		factors = factors[1:]
	}
	var num, den []string
	for _, f := range factors {
		switch x := f.(type) {
		case *Rational:
			n, d := x.Val.Num(), x.Val.Denom()
			ns := n.String()
			if n.Sign() < 0 {
				ns = new(big.Int).Neg(n).String()
				neg = !neg
			}
			num = append(num, ns)
			den = append(den, d.String())
		case *power:
			if isNum(x.exp) && numSign(x.exp) < 0 {
				den = append(den, ltxPowPositive(x.base, numNeg(x.exp)))
			} else {
				num = append(num, ltxMulFactor(f))
			}
		default:
			num = append(num, ltxMulFactor(f))
		}
	}
	numStr := ltxJoinFactors(num)
	if numStr == "" {
		numStr = "1"
	}
	var s string
	if len(den) > 0 {
		s = "\\frac{" + numStr + "}{" + ltxJoinFactors(den) + "}"
	} else {
		s = numStr
	}
	if neg {
		s = "-" + s
	}
	return s
}

// ltxJoinFactors joins rendered factor strings, dropping redundant unit
// numerators and inserting \cdot only between two numeric operands.
func ltxJoinFactors(tokens []string) string {
	kept := make([]string, 0, len(tokens))
	for _, t := range tokens {
		if t == "1" && len(tokens) > 1 {
			continue
		}
		kept = append(kept, t)
	}
	if len(kept) == 0 {
		return ""
	}
	var b strings.Builder
	b.WriteString(kept[0])
	for i := 1; i < len(kept); i++ {
		if ltxEndsDigit(kept[i-1]) && ltxStartsDigit(kept[i]) {
			b.WriteString(" \\cdot ")
		} else {
			b.WriteString(" ")
		}
		b.WriteString(kept[i])
	}
	return b.String()
}

// ltxEndsDigit reports whether s ends with an ASCII digit.
func ltxEndsDigit(s string) bool {
	if s == "" {
		return false
	}
	c := s[len(s)-1]
	return c >= '0' && c <= '9'
}

// ltxStartsDigit reports whether s begins with an ASCII digit.
func ltxStartsDigit(s string) bool {
	if s == "" {
		return false
	}
	c := s[0]
	return c >= '0' && c <= '9'
}

// ltxPow renders a power, using \sqrt for half powers and \frac for negative
// exponents.
func ltxPow(p *power) string {
	if presIsHalf(p.exp) {
		return "\\sqrt{" + ltxOf(p.base) + "}"
	}
	if isNum(p.exp) && numSign(p.exp) < 0 {
		return "\\frac{1}{" + ltxPowPositive(p.base, numNeg(p.exp)) + "}"
	}
	return ltxWrap(p.base) + "^{" + ltxOf(p.exp) + "}"
}

// ltxPowPositive renders base raised to the positive exponent pos, used to
// build the denominator of a reciprocal power.
func ltxPowPositive(base, pos Expr) string {
	if presIsHalf(pos) {
		return "\\sqrt{" + ltxOf(base) + "}"
	}
	if isOne(pos) {
		return ltxWrap(base)
	}
	return ltxWrap(base) + "^{" + ltxOf(pos) + "}"
}

// ltxFn renders a single-argument function application.
func ltxFn(f *fn) string {
	arg := ltxOf(f.arg)
	switch f.name {
	case "sqrt":
		return "\\sqrt{" + arg + "}"
	case "abs":
		return "\\left|" + arg + "\\right|"
	case "floor":
		return "\\left\\lfloor " + arg + " \\right\\rfloor"
	case "ceil":
		return "\\left\\lceil " + arg + " \\right\\rceil"
	case "conjugate":
		return "\\overline{" + arg + "}"
	case "re":
		return "\\operatorname{Re}\\left(" + arg + "\\right)"
	case "im":
		return "\\operatorname{Im}\\left(" + arg + "\\right)"
	case "factorial":
		if needParen(f.arg) || isFnCall(f.arg) {
			return "\\left(" + arg + "\\right)!"
		}
		return arg + "!"
	}
	if cmd, ok := ltxFnCmd[f.name]; ok {
		return cmd + "\\left(" + arg + "\\right)"
	}
	return "\\operatorname{" + f.name + "}\\left(" + arg + "\\right)"
}

// ltxFn2 renders a two-argument function application.
func ltxFn2(f *fn2) string {
	name := f.name
	if name == "beta" {
		name = "B"
	}
	return "\\operatorname{" + name + "}\\left(" + ltxOf(f.arg1) + ", " + ltxOf(f.arg2) + "\\right)"
}

// ltxBigOp renders a summation, product or limit.
func ltxBigOp(b *bigOp) string {
	body := ltxOf(b.body)
	if _, ok := b.body.(*sum); ok {
		body = "\\left(" + body + "\\right)"
	}
	switch b.kind {
	case "Limit":
		return "\\lim_{" + b.index + " \\to " + ltxOf(b.lo) + "} " + body
	case "Product":
		return "\\prod_{" + b.index + "=" + ltxOf(b.lo) + "}^{" + ltxOf(b.hi) + "} " + body
	default: // Sum
		return "\\sum_{" + b.index + "=" + ltxOf(b.lo) + "}^{" + ltxOf(b.hi) + "} " + body
	}
}

// =========================================================================
// Pretty (2D Unicode)
// =========================================================================

// prettyBlock is the internal two-dimensional layout primitive used by Pretty.
// lines holds the rows of the rendering, all padded to a common width, and
// baseline is the index of the row on which horizontally adjacent blocks are
// aligned (for example the fraction bar of a stacked fraction).
type prettyBlock struct {
	lines    []string
	baseline int
}

// Pretty renders e as a multi-line Unicode drawing in the style of a computer
// algebra system's pretty printer: fractions are stacked over a horizontal
// bar, exponents are raised onto the line above, square roots are drawn under
// a radical, and the terms of a sum are aligned on a common baseline. The
// result is the rows joined by newlines, with no trailing newline. The output
// is deterministic.
func Pretty(e Expr) string {
	b := ppOf(e)
	return strings.Join(b.lines, "\n")
}

// ppRuneLen returns the display width of s in runes (each supported glyph is a
// single fixed-width cell).
func ppRuneLen(s string) int { return len([]rune(s)) }

// ppPadRight right-pads s with spaces to a width of w runes.
func ppPadRight(s string, w int) string {
	if n := ppRuneLen(s); n < w {
		return s + strings.Repeat(" ", w-n)
	}
	return s
}

// ppCenter centers s within a field of w runes.
func ppCenter(s string, w int) string {
	pad := w - ppRuneLen(s)
	if pad <= 0 {
		return s
	}
	left := pad / 2
	return strings.Repeat(" ", left) + s + strings.Repeat(" ", pad-left)
}

// ppBlockWidth returns the common rune width of a block.
func ppBlockWidth(b prettyBlock) int {
	w := 0
	for _, l := range b.lines {
		if n := ppRuneLen(l); n > w {
			w = n
		}
	}
	return w
}

// ppNormalize pads every line to the widest line's width.
func ppNormalize(lines []string) []string {
	w := 0
	for _, l := range lines {
		if n := ppRuneLen(l); n > w {
			w = n
		}
	}
	out := make([]string, len(lines))
	for i, l := range lines {
		out[i] = ppPadRight(l, w)
	}
	return out
}

// ppAtom builds a single-line block whose baseline is that line.
func ppAtom(s string) prettyBlock { return prettyBlock{lines: []string{s}, baseline: 0} }

// ppHConcat concatenates blocks left to right, aligning them on their
// baselines and padding above and below as needed.
func ppHConcat(blocks ...prettyBlock) prettyBlock {
	if len(blocks) == 0 {
		return ppAtom("")
	}
	top, bottom := 0, 0
	for _, b := range blocks {
		if b.baseline > top {
			top = b.baseline
		}
		if d := len(b.lines) - b.baseline - 1; d > bottom {
			bottom = d
		}
	}
	height := top + bottom + 1
	rows := make([]string, height)
	for _, b := range blocks {
		w := ppBlockWidth(b)
		above := top - b.baseline
		blank := strings.Repeat(" ", w)
		for i := 0; i < height; i++ {
			switch {
			case i < above:
				rows[i] += blank
			case i-above < len(b.lines):
				rows[i] += ppPadRight(b.lines[i-above], w)
			default:
				rows[i] += blank
			}
		}
	}
	return prettyBlock{lines: rows, baseline: top}
}

// ppParen wraps a block in parentheses, using tall bracket glyphs when the
// block spans more than one line.
func ppParen(b prettyBlock) prettyBlock {
	h := len(b.lines)
	if h <= 1 {
		return ppHConcat(ppAtom("("), b, ppAtom(")"))
	}
	left := make([]string, h)
	right := make([]string, h)
	for i := 0; i < h; i++ {
		switch {
		case i == 0:
			left[i], right[i] = "⎛", "⎞"
		case i == h-1:
			left[i], right[i] = "⎝", "⎠"
		default:
			left[i], right[i] = "⎜", "⎟"
		}
	}
	lb := prettyBlock{lines: left, baseline: b.baseline}
	rb := prettyBlock{lines: right, baseline: b.baseline}
	return ppHConcat(lb, b, rb)
}

// ppFrac stacks num over den separated by a horizontal bar, aligning on the
// bar.
func ppFrac(num, den prettyBlock) prettyBlock {
	w := ppBlockWidth(num)
	if dw := ppBlockWidth(den); dw > w {
		w = dw
	}
	lines := make([]string, 0, len(num.lines)+1+len(den.lines))
	for _, l := range num.lines {
		lines = append(lines, ppCenter(l, w))
	}
	lines = append(lines, strings.Repeat("─", w))
	for _, l := range den.lines {
		lines = append(lines, ppCenter(l, w))
	}
	return prettyBlock{lines: lines, baseline: len(num.lines)}
}

// ppPow raises exp onto the lines above base, keeping base's baseline.
func ppPow(base, exp prettyBlock) prettyBlock {
	bw := ppBlockWidth(base)
	ew := ppBlockWidth(exp)
	lines := make([]string, 0, len(exp.lines)+len(base.lines))
	for _, l := range exp.lines {
		lines = append(lines, strings.Repeat(" ", bw)+ppPadRight(l, ew))
	}
	for _, l := range base.lines {
		lines = append(lines, ppPadRight(l, bw)+strings.Repeat(" ", ew))
	}
	return prettyBlock{lines: lines, baseline: len(exp.lines) + base.baseline}
}

// ppSqrt draws inner under a radical sign.
func ppSqrt(inner prettyBlock) prettyBlock {
	h := len(inner.lines)
	w := ppBlockWidth(inner)
	lines := make([]string, 0, h+1)
	lines = append(lines, strings.Repeat(" ", h+1)+strings.Repeat("_", w+1))
	for i := 0; i < h; i++ {
		slashCol := h - 1 - i
		stem := strings.Repeat(" ", slashCol) + "/" + strings.Repeat(" ", h-slashCol)
		lines = append(lines, stem+" "+ppPadRight(inner.lines[i], w))
	}
	return prettyBlock{lines: ppNormalize(lines), baseline: 1 + inner.baseline}
}

// ppOf renders e as a prettyBlock.
func ppOf(e Expr) prettyBlock {
	switch x := e.(type) {
	case *Integer:
		return ppAtom(x.Val.String())
	case *Rational:
		return ppFrac(ppAtom(x.Val.Num().String()), ppAtom(x.Val.Denom().String()))
	case *Float:
		return ppAtom(strconv.FormatFloat(x.Val, 'g', -1, 64))
	case *Symbol:
		return ppAtom(x.Name)
	case *Constant:
		return ppAtom(ppConst(x))
	case *sum:
		return ppSum(x)
	case *product:
		return ppProduct(x)
	case *power:
		return ppPowNode(x)
	case *fn:
		return ppFn(x)
	case *fn2:
		return ppHConcat(ppAtom(x.name), ppParen(ppHConcat(ppOf(x.arg1), ppAtom(", "), ppOf(x.arg2))))
	case *integral:
		return ppHConcat(ppAtom("∫ "), ppOf(x.arg), ppAtom(" d"), ppOf(x.v))
	case *bigOp:
		return ppBigOp(x)
	}
	return ppAtom(e.String())
}

// ppConst renders a named constant as a Unicode glyph.
func ppConst(c *Constant) string {
	switch c.Name {
	case "pi":
		return "π"
	case "E":
		return "e"
	case "I":
		return "i"
	case "oo":
		return "∞"
	case "-oo":
		return "-∞"
	}
	return c.Name
}

// ppBaseBlock renders a power base, wrapping it in parentheses when needed.
func ppBaseBlock(base Expr) prettyBlock {
	if needParen(base) {
		return ppParen(ppOf(base))
	}
	return ppOf(base)
}

// ppPowPositiveBlock renders base raised to the positive exponent pos.
func ppPowPositiveBlock(base, pos Expr) prettyBlock {
	if presIsHalf(pos) {
		return ppSqrt(ppOf(base))
	}
	if isOne(pos) {
		return ppBaseBlock(base)
	}
	return ppPow(ppBaseBlock(base), ppOf(pos))
}

// ppPowNode renders a power node, using a radical for half powers and a
// fraction for negative exponents.
func ppPowNode(p *power) prettyBlock {
	if presIsHalf(p.exp) {
		return ppSqrt(ppOf(p.base))
	}
	if isNum(p.exp) && numSign(p.exp) < 0 {
		return ppFrac(ppAtom("1"), ppPowPositiveBlock(p.base, numNeg(p.exp)))
	}
	return ppPow(ppBaseBlock(p.base), ppOf(p.exp))
}

// ppMulFactorBlock renders a product factor, wrapping embedded sums.
func ppMulFactorBlock(f Expr) prettyBlock {
	if _, ok := f.(*sum); ok {
		return ppParen(ppOf(f))
	}
	return ppOf(f)
}

// ppJoinMul joins factor blocks with the multiplication dot.
func ppJoinMul(blocks []prettyBlock) prettyBlock {
	parts := make([]prettyBlock, 0, len(blocks)*2)
	for i, b := range blocks {
		if i > 0 {
			parts = append(parts, ppAtom("⋅"))
		}
		parts = append(parts, b)
	}
	return ppHConcat(parts...)
}

// ppProduct renders a product, splitting negative powers and rational
// coefficients into a stacked fraction.
func ppProduct(p *product) prettyBlock {
	factors := p.factors
	neg := false
	if len(factors) > 1 && isMinusOne(factors[0]) {
		neg = true
		factors = factors[1:]
	}
	var num, den []prettyBlock
	for _, f := range factors {
		switch x := f.(type) {
		case *Rational:
			n, d := x.Val.Num(), x.Val.Denom()
			ns := n.String()
			if n.Sign() < 0 {
				ns = new(big.Int).Neg(n).String()
				neg = !neg
			}
			if !(ns == "1" && len(factors) > 1) {
				num = append(num, ppAtom(ns))
			}
			den = append(den, ppAtom(d.String()))
		case *power:
			if isNum(x.exp) && numSign(x.exp) < 0 {
				den = append(den, ppPowPositiveBlock(x.base, numNeg(x.exp)))
			} else {
				num = append(num, ppMulFactorBlock(f))
			}
		default:
			num = append(num, ppMulFactorBlock(f))
		}
	}
	if len(num) == 0 {
		num = append(num, ppAtom("1"))
	}
	res := ppJoinMul(num)
	if len(den) > 0 {
		res = ppFrac(res, ppJoinMul(den))
	}
	if neg {
		res = ppHConcat(ppAtom("-"), res)
	}
	return res
}

// ppSum renders a sum with terms aligned on a common baseline.
func ppSum(s *sum) prettyBlock {
	terms := presSortedSumTerms(s)
	parts := make([]prettyBlock, 0, len(terms)*2)
	for i, t := range terms {
		neg, mag := splitSign(t)
		mb := ppOf(mag)
		if _, ok := mag.(*sum); ok {
			mb = ppParen(mb)
		}
		switch {
		case i == 0 && neg:
			parts = append(parts, ppAtom("-"), mb)
		case i == 0:
			parts = append(parts, mb)
		case neg:
			parts = append(parts, ppAtom(" - "), mb)
		default:
			parts = append(parts, ppAtom(" + "), mb)
		}
	}
	return ppHConcat(parts...)
}

// ppFn renders a single-argument function application.
func ppFn(f *fn) prettyBlock {
	switch f.name {
	case "sqrt":
		return ppSqrt(ppOf(f.arg))
	case "abs":
		return ppHConcat(ppAtom("|"), ppOf(f.arg), ppAtom("|"))
	case "floor":
		return ppHConcat(ppAtom("⌊"), ppOf(f.arg), ppAtom("⌋"))
	case "ceil":
		return ppHConcat(ppAtom("⌈"), ppOf(f.arg), ppAtom("⌉"))
	case "factorial":
		arg := ppOf(f.arg)
		if needParen(f.arg) || isFnCall(f.arg) {
			arg = ppParen(arg)
		}
		return ppHConcat(arg, ppAtom("!"))
	}
	name := f.name
	switch f.name {
	case "gamma":
		name = "Γ"
	case "digamma":
		name = "ψ"
	}
	return ppHConcat(ppAtom(name), ppParen(ppOf(f.arg)))
}

// ppBigOp renders a summation, product or limit on a single baseline row with
// its bounds shown inline.
func ppBigOp(b *bigOp) prettyBlock {
	body := ppOf(b.body)
	if _, ok := b.body.(*sum); ok {
		body = ppParen(body)
	}
	switch b.kind {
	case "Limit":
		return ppHConcat(ppAtom("lim("+b.index+" → "+ppInline(b.lo)+") "), body)
	case "Product":
		return ppHConcat(ppAtom("∏("+b.index+"="+ppInline(b.lo)+".."+ppInline(b.hi)+") "), body)
	default:
		return ppHConcat(ppAtom("∑("+b.index+"="+ppInline(b.lo)+".."+ppInline(b.hi)+") "), body)
	}
}

// ppInline renders a bound expression as a single-line string for use inside a
// big-operator label.
func ppInline(e Expr) string { return e.String() }

// =========================================================================
// MathML
// =========================================================================

// MathML renders e as a presentation MathML fragment wrapped in a <math>
// element. Fractions become <mfrac>, powers <msup>, square roots <msqrt> and
// grouped subexpressions <mrow>, covering the same node set as LaTeX and
// Pretty. The output is deterministic.
func MathML(e Expr) string {
	return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">" + mmlOf(e) + "</math>"
}

// mmlInvTimes is the invisible-times operator inserted between juxtaposed
// factors.
const mmlInvTimes = "<mo>\u2062</mo>"

// mmlOf renders e as presentation MathML without the enclosing <math> element.
func mmlOf(e Expr) string {
	switch x := e.(type) {
	case *Integer:
		return "<mn>" + x.Val.String() + "</mn>"
	case *Rational:
		n, d := x.Val.Num(), x.Val.Denom()
		if n.Sign() < 0 {
			frac := "<mfrac><mn>" + new(big.Int).Neg(n).String() + "</mn><mn>" + d.String() + "</mn></mfrac>"
			return "<mrow><mo>-</mo>" + frac + "</mrow>"
		}
		return "<mfrac><mn>" + n.String() + "</mn><mn>" + d.String() + "</mn></mfrac>"
	case *Float:
		return "<mn>" + strconv.FormatFloat(x.Val, 'g', -1, 64) + "</mn>"
	case *Symbol:
		return "<mi>" + x.Name + "</mi>"
	case *Constant:
		return mmlConst(x)
	case *sum:
		return mmlSum(x)
	case *product:
		return mmlProduct(x)
	case *power:
		return mmlPow(x)
	case *fn:
		return mmlFn(x)
	case *fn2:
		return "<mrow><mi>" + x.name + "</mi><mo>\u2061</mo><mrow><mo>(</mo>" +
			mmlOf(x.arg1) + "<mo>,</mo>" + mmlOf(x.arg2) + "<mo>)</mo></mrow></mrow>"
	case *integral:
		return "<mrow><mo>∫</mo>" + mmlOf(x.arg) + "<mo>d</mo>" + mmlOf(x.v) + "</mrow>"
	case *bigOp:
		return mmlBigOp(x)
	}
	return "<mi>" + e.String() + "</mi>"
}

// mmlConst renders a named constant.
func mmlConst(c *Constant) string {
	switch c.Name {
	case "pi":
		return "<mi>π</mi>"
	case "E":
		return "<mi>e</mi>"
	case "I":
		return "<mi>i</mi>"
	case "oo":
		return "<mi>∞</mi>"
	case "-oo":
		return "<mrow><mo>-</mo><mi>∞</mi></mrow>"
	}
	return "<mi>" + c.Name + "</mi>"
}

// mmlBase renders a power base, wrapping it in a parenthesized <mrow> when
// needed.
func mmlBase(e Expr) string {
	if needParen(e) {
		return "<mrow><mo>(</mo>" + mmlOf(e) + "<mo>)</mo></mrow>"
	}
	return mmlOf(e)
}

// mmlMulFactor renders a product factor, wrapping embedded sums.
func mmlMulFactor(f Expr) string {
	if _, ok := f.(*sum); ok {
		return "<mrow><mo>(</mo>" + mmlOf(f) + "<mo>)</mo></mrow>"
	}
	return mmlOf(f)
}

// mmlSum renders a sum.
func mmlSum(s *sum) string {
	terms := presSortedSumTerms(s)
	var b strings.Builder
	b.WriteString("<mrow>")
	for i, t := range terms {
		neg, mag := splitSign(t)
		ms := mmlMulFactor(mag)
		switch {
		case i == 0 && neg:
			b.WriteString("<mo>-</mo>" + ms)
		case i == 0:
			b.WriteString(ms)
		case neg:
			b.WriteString("<mo>-</mo>" + ms)
		default:
			b.WriteString("<mo>+</mo>" + ms)
		}
	}
	b.WriteString("</mrow>")
	return b.String()
}

// mmlProduct renders a product, using <mfrac> for negative powers and rational
// coefficients.
func mmlProduct(p *product) string {
	factors := p.factors
	neg := false
	if len(factors) > 1 && isMinusOne(factors[0]) {
		neg = true
		factors = factors[1:]
	}
	var num, den []string
	for _, f := range factors {
		switch x := f.(type) {
		case *Rational:
			n, d := x.Val.Num(), x.Val.Denom()
			ns := n.String()
			if n.Sign() < 0 {
				ns = new(big.Int).Neg(n).String()
				neg = !neg
			}
			if !(ns == "1" && len(factors) > 1) {
				num = append(num, "<mn>"+ns+"</mn>")
			}
			den = append(den, "<mn>"+d.String()+"</mn>")
		case *power:
			if isNum(x.exp) && numSign(x.exp) < 0 {
				den = append(den, mmlPowPositive(x.base, numNeg(x.exp)))
			} else {
				num = append(num, mmlMulFactor(f))
			}
		default:
			num = append(num, mmlMulFactor(f))
		}
	}
	numStr := mmlJoin(num)
	if numStr == "" {
		numStr = "<mn>1</mn>"
	}
	var s string
	if len(den) > 0 {
		s = "<mfrac>" + mmlWrapRow(numStr) + mmlWrapRow(mmlJoin(den)) + "</mfrac>"
	} else {
		s = numStr
	}
	if neg {
		s = "<mrow><mo>-</mo>" + s + "</mrow>"
	}
	return s
}

// mmlWrapRow wraps s in an <mrow> so it forms a single fraction operand.
func mmlWrapRow(s string) string { return "<mrow>" + s + "</mrow>" }

// mmlJoin joins factor strings with the invisible-times operator.
func mmlJoin(tokens []string) string {
	return strings.Join(tokens, mmlInvTimes)
}

// mmlPow renders a power, using <msqrt> for half powers and <mfrac> for
// negative exponents.
func mmlPow(p *power) string {
	if presIsHalf(p.exp) {
		return "<msqrt>" + mmlOf(p.base) + "</msqrt>"
	}
	if isNum(p.exp) && numSign(p.exp) < 0 {
		return "<mfrac><mn>1</mn>" + mmlWrapRow(mmlPowPositive(p.base, numNeg(p.exp))) + "</mfrac>"
	}
	return "<msup>" + mmlBase(p.base) + mmlOf(p.exp) + "</msup>"
}

// mmlPowPositive renders base raised to the positive exponent pos.
func mmlPowPositive(base, pos Expr) string {
	if presIsHalf(pos) {
		return "<msqrt>" + mmlOf(base) + "</msqrt>"
	}
	if isOne(pos) {
		return mmlBase(base)
	}
	return "<msup>" + mmlBase(base) + mmlOf(pos) + "</msup>"
}

// mmlFn renders a single-argument function application.
func mmlFn(f *fn) string {
	arg := mmlOf(f.arg)
	switch f.name {
	case "sqrt":
		return "<msqrt>" + arg + "</msqrt>"
	case "abs":
		return "<mrow><mo>|</mo>" + arg + "<mo>|</mo></mrow>"
	case "floor":
		return "<mrow><mo>⌊</mo>" + arg + "<mo>⌋</mo></mrow>"
	case "ceil":
		return "<mrow><mo>⌈</mo>" + arg + "<mo>⌉</mo></mrow>"
	case "factorial":
		return "<mrow>" + mmlMulFactor(f.arg) + "<mo>!</mo></mrow>"
	}
	name := f.name
	switch f.name {
	case "gamma":
		name = "Γ"
	case "digamma":
		name = "ψ"
	}
	return "<mrow><mi>" + name + "</mi><mo>\u2061</mo><mrow><mo>(</mo>" + arg + "<mo>)</mo></mrow></mrow>"
}

// mmlBigOp renders a summation, product or limit.
func mmlBigOp(b *bigOp) string {
	body := mmlMulFactor(b.body)
	switch b.kind {
	case "Limit":
		under := "<mrow><mi>" + b.index + "</mi><mo>→</mo>" + mmlOf(b.lo) + "</mrow>"
		return "<mrow><munder><mo>lim</mo>" + under + "</munder>" + body + "</mrow>"
	case "Product":
		return "<mrow><munderover><mo>∏</mo>" + mmlBounds(b) + "</munderover>" + body + "</mrow>"
	default:
		return "<mrow><munderover><mo>∑</mo>" + mmlBounds(b) + "</munderover>" + body + "</mrow>"
	}
}

// mmlBounds renders the lower and upper bounds of a summation or product.
func mmlBounds(b *bigOp) string {
	lower := "<mrow><mi>" + b.index + "</mi><mo>=</mo>" + mmlOf(b.lo) + "</mrow>"
	upper := mmlOf(b.hi)
	return lower + upper
}
