// Package auctions implements mechanism design and cooperative game theory in
// pure Go: sealed-bid auctions, the Vickrey-Clarke-Groves (VCG) mechanism with
// combinatorial winner determination, and the classical solution concepts of
// transferable-utility cooperative games together with two-sided matching.
//
// The package is organised around a few concrete types and a large collection
// of small, composable functions:
//
//   - sealed-bid auctions: first-price, second-price (Vickrey), third-price and
//     all-pay single-item formats with optional reserve prices, plus multi-unit
//     uniform-price, pay-as-bid and multi-unit Vickrey clearings. See Bid,
//     FirstPriceAuction, VickreyAuction and MultiUnitVickrey.
//
//   - the VCG mechanism: exact combinatorial winner determination over XOR and
//     OR bids by branch-and-bound, Clarke pivot payments and the resulting
//     bidder utilities. See CombinatorialBid, WinnerDetermination and
//     VCGMechanism.
//
//   - cooperative games: the CoopGame type carries a characteristic function on
//     2^n coalitions and computes the Shapley value, the Banzhaf value and power
//     index, marginal contributions, the core and its membership test, the
//     least-core value and the (pre)nucleolus via a self-contained two-phase
//     simplex, along with structural predicates (superadditivity, convexity,
//     monotonicity and more). Weighted voting games expose the Shapley-Shubik
//     and Banzhaf power indices, minimal winning coalitions, veto players and
//     dummies.
//
//   - bargaining: the Nash bargaining solution and the Kalai-Smorodinsky,
//     egalitarian and utilitarian solutions over a convex two-player feasible
//     set, with convex-hull and Pareto-frontier helpers.
//
//   - matching: the Gale-Shapley deferred-acceptance algorithm for stable
//     marriage (man- and woman-optimal), stability and blocking-pair tests,
//     Gale's top-trading-cycles algorithm for housing markets and serial
//     dictatorship.
//
// Everything uses only the Go standard library and is deterministic; the single
// Monte-Carlo estimator (ShapleyValueMonteCarlo) draws from a caller-supplied
// seed through math/rand so that results are reproducible. All numeric routines
// are validated in the accompanying tests against closed-form reference values.
package auctions
