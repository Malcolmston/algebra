// Package gametheory implements classical game theory: normal-form games,
// Nash equilibria (pure and mixed), zero-sum solutions via linear programming,
// dominated-strategy elimination, Pareto frontiers, cooperative solution
// concepts (the Shapley and Banzhaf values, the core), and a small toolkit of
// Prisoner's-dilemma and canonical-game helpers.
//
// The central type is Game, a two-player normal-form game holding separate
// payoff matrices for the row and column players. On top of it the package
// offers:
//
//   - equilibrium analysis: PureNashEquilibria and, via support enumeration,
//     the full set of mixed equilibria (MixedNashEquilibria);
//   - zero-sum solving: SolveZeroSum computes the value of the game and a pair
//     of optimal minimax strategies by linear programming (a self-contained
//     Bland's-rule simplex), with Maximin, Minimax and SaddlePoint helpers for
//     the pure case;
//   - dominance: strict and weak domination tests and iterated elimination;
//   - efficiency: Pareto frontiers and the social optimum;
//   - cooperation: CooperativeGame with the Shapley value, the Banzhaf value,
//     core membership tests and structural predicates (superadditivity,
//     monotonicity, convexity), plus WeightedVotingGame;
//   - the Prisoner's dilemma: PrisonersDilemma with iterated-play simulation
//     and stock strategies (TitForTat, GrimTrigger, WinStayLoseShift, ...),
//     and constructors for the canonical games MatchingPennies, StagHunt,
//     BattleOfTheSexes and GameOfChicken.
//
// All routines use only the Go standard library, are deterministic, and are
// validated in the accompanying tests against closed-form reference values.
package gametheory
