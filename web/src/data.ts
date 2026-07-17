// Library content for the algebra documentation site. Mirrors the shape used by
// the malcolmston/go landing site's data.ts so the sibling sites stay in sync.
export interface Lib {
  id: string; name: string; icon: string; accent: string; pkg: string; node: string;
  repo: string; docs: string; tagline: string; blurb: string; tags: string[];
  features: string[]; node_code: string; go_code: string; integrate: string;
}

export const NODE_ACCENT = '#8cc84b';

export const ALGEBRA: Lib = {
  id:"algebra", name:"algebra", icon:'<i class="fa-solid fa-square-root-variable"></i>', accent:"#34c7a8",
  pkg:"github.com/malcolmston/algebra", node:"sympy/sympy",
  repo:"https://github.com/malcolmston/algebra", docs:"https://malcolmston.github.io/algebra/",
  tagline:"Symbolic algebra in Go.",
  blurb:"A from-scratch, standard-library-only computer-algebra system (CAS) for Go. It represents mathematical "+
    "expressions as immutable, value-based trees, parses ordinary infix strings, and manipulates them symbolically — "+
    "simplify, expand, differentiate, integrate, take limits and series, substitute, evaluate numerically or over the "+
    "complex plane, and solve polynomial equations of any degree. The canonicalizing Add / Mul / Pow constructors mean "+
    "mathematically equal expressions compare equal and print identically. Four companion subpackages add symbolic "+
    "matrix linear algebra, number theory, statistics and physics helpers. A faithful, idiomatic port of a subset of "+
    "Python's SymPy.",
  tags:["expression trees","infix parser","simplify","differentiation","integration","equation solving","complex numbers","linear algebra","number theory","math/big exact arithmetic"],
  features:[
    "<code>Expr</code> trees — immutable, value-based nodes: symbols, arbitrary-precision <code>Integer</code> / <code>Rational</code> (via <code>math/big</code>), floats, constants (<code>Pi</code>, <code>E</code>), sums, products, powers and elementary functions",
    "A precedence-climbing <code>Parse</code> for ordinary infix notation: <code>+ - * / ^</code>, unary signs, parentheses, functions (<code>sin</code>, <code>cos</code>, <code>tan</code>, <code>exp</code>, <code>log</code>, <code>sqrt</code>) and implicit multiplication (<code>2x</code>, <code>2(x+1)</code>)",
    "<code>Simplify</code> and <code>Expand</code> — canonicalization, like-term collection, algebraic identities and multinomial expansion",
    "<code>Diff</code> — symbolic differentiation with the sum, product, power, quotient and chain rules",
    "<code>Integrate</code> — symbolic antiderivatives of a documented subset, returning an unevaluated <code>Integral</code> rather than a wrong answer",
    "Full <b>trigonometry</b> &amp; <b>hyperbolic</b> functions (with inverses and exact special-angle values), <b>complex numbers</b> (<code>I</code>, <code>Conjugate</code>, <code>Re</code>/<code>Im</code>, <code>Arg</code>, Euler folds), and <b>special functions</b> (<code>Gamma</code>, <code>Beta</code>, <code>Erf</code>, <code>Factorial</code>)",
    "Advanced calculus — <code>Limit</code> (with L'Hôpital), <code>Series</code> (Taylor/Maclaurin), <code>Summation</code> / <code>Product</code>, and richer <code>Integrate</code> (by parts, partial fractions, arctan/asin forms)",
    "<code>Solve</code> — polynomial equations of any degree: exact rational, quadratic, cubic and quartic factors with complex conjugate roots and a numeric Durand–Kerner fallback for irreducible factors; <code>SolveSystem</code> solves linear systems by Gaussian elimination",
    "<b>Four subpackages</b>: <code>matrix</code> (symbolic linear algebra — <code>Det</code>, <code>Inverse</code>, <code>RREF</code>, <code>Eigenvalues</code>, <code>CharPoly</code>, solve <code>Ax=b</code>), <code>ntheory</code> (primes, factorization, modular arithmetic, combinatorics, integer sequences), <code>stats</code> (descriptive statistics, linear regression &amp; eight probability distributions) and <code>physics</code> (CODATA constants, unit conversion, kinematics/relativity/EM formulas)",
    "<code>Eval</code> / <code>Evalf</code> / <code>Evalc</code> — numeric &amp; complex evaluation. Zero dependencies — pure Go standard library"
  ],
  node_code:
`import sympy

x = sympy.Symbol("x")
e = x**2 + 2*x + 1

# f'(x) = 2*x + 2
print(sympy.simplify(sympy.diff(e, x)))

# solve x^2 - 5*x + 6 = 0  ->  [2, 3]
print(sympy.solve(x**2 - 5*x + 6, x))`,
  go_code:
`import "github.com/malcolmston/algebra"

x := algebra.Sym("x")

// f(x) = x^2 + 2*x + 1
e := algebra.MustParse("x^2 + 2*x + 1")

// f'(x) = 2*x + 2
fmt.Println(algebra.Simplify(algebra.Diff(e, x)))

// solve x^2 - 5*x + 6 = 0  ->  [2 3]
roots, _ := algebra.Solve(algebra.MustParse("x^2 - 5*x + 6"), x)
fmt.Println(roots)`,
  integrate:
`<span class="tok-c">// Expand a binomial power (multinomial expansion)</span>
algebra.Expand(algebra.MustParse("(x + 1)^3"))    <span class="tok-c">// x^3 + 3*x^2 + 3*x + 1</span>

<span class="tok-c">// Factor a quadratic into monic linear factors</span>
algebra.Factor(algebra.MustParse("x^2 - 1"), x)   <span class="tok-c">// (x - 1)*(x + 1)</span>

<span class="tok-c">// Symbolic integration of a documented subset</span>
algebra.Integrate(algebra.MustParse("cos(x)"), x) <span class="tok-c">// sin(x)</span>

<span class="tok-c">// Substitute y -> x + 1, then expand</span>
algebra.Subs(algebra.MustParse("y^2"), algebra.Sym("y"),
    algebra.MustParse("x + 1")).Expand()          <span class="tok-c">// x^2 + 2*x + 1</span>

<span class="tok-c">// Numeric evaluation of a bound expression</span>
algebra.Eval(algebra.MustParse("x^2 + 1"),
    map[string]float64{"x": 3})                   <span class="tok-c">// 10, nil</span>`
};
