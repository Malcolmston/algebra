// Package logic implements propositional logic and Boolean algebra using only
// the Go standard library.
//
// The package is organised into a few areas:
//
//   - Boolean primitives (ops.go): the elementary connectives [And], [Or],
//     [Not], [Xor], [Nand], [Nor], [Xnor], [Implies] and [Iff], together with
//     the variadic reductions [AndAll], [OrAll], [XorAll] and the derived
//     gates [Majority] and [Mux].
//   - Expression trees, parsing and evaluation (expr.go): the sealed [Expr]
//     interface implemented by [Var], [Const], [UnaryExpr] and [BinaryExpr];
//     the constructors [NewVar], [NewConst], [NewNot], [NewAnd], [NewOr],
//     [NewXor], [NewNand], [NewNor], [NewXnor], [NewImplies] and [NewIff]; the
//     [Parse]/[MustParse] recursive-descent parser; and the helpers [Vars] and
//     [EvalString].
//   - Truth tables (truthtable.go): the [TruthTable] and [TruthRow] types built
//     by [NewTruthTable], plus [Assignments], [Minterms], [Maxterms],
//     [IndexToAssignment] and [AssignmentToIndex].
//   - Semantic analysis (analyze.go): [IsTautology], [IsContradiction],
//     [IsSatisfiable], [IsContingency], [FindModel], [AllModels],
//     [CountModels], [Equivalent] and [Entails].
//   - Normal forms (normalform.go): [ToNNF], [ToCNF], [ToDNF], [Simplify],
//     [CNFString], [DNFString], [IsCNF] and [IsDNF].
//   - Two-level minimisation (minimize.go): the [Implicant] cube type, the
//     [QuineMcCluskey] engine with [PrimeImplicants] and
//     [EssentialPrimeImplicants], the cover routines [MinimizeSOP] and
//     [SOPString], the Gray-code helpers [GrayCode], [GrayEncode] and
//     [GrayDecode], the [KarnaughMap] grouping type built by [NewKarnaughMap],
//     and the bit utilities [PopCount] and [BitString].
//
// The syntax accepted by [Parse] uses C-like and word operators. Variables are
// identifiers; the constants are written 1/0, T/F or true/false. Connectives,
// from lowest to highest binding, are: <-> (iff), -> (implies), | (or), ^
// (xor), & (and) and the unary ! (not); the words and, or, not, xor, nand, nor,
// xnor, implies and iff are accepted as well, as are the aliases && * . for
// and, || + for or and ~ for not.
//
// Every routine is deterministic and depends only on the Go standard library.
package logic
