// Package socialchoice implements classical social-choice theory: ranked and
// cardinal voting rules, pairwise/tournament analysis, and seat apportionment.
//
// The central ranked-ballot type is Profile, a weighted collection of Ballot
// preference orders over a fixed set of candidates. Ballots may be truncated:
// a candidate that a ballot does not list is treated as ranked below every
// listed candidate and as tied with the other unlisted candidates. On top of
// Profile the package offers:
//
//   - positional rules: Plurality, AntiPlurality, Borda, Dowdall and the
//     generic PositionalScores, plus the iterative Borda methods Nanson and
//     Baldwin;
//   - pairwise analysis: the PairwiseMatrix of majority preferences, Condorcet
//     winners and losers, Copeland and MiniMax scores, and the Smith and
//     Schwartz sets computed from the majority tournament;
//   - Condorcet-consistent completions: the Schulze beatpath method and
//     Tideman's ranked pairs;
//   - elimination rules: instant-runoff voting, the single transferable vote
//     with Droop-quota Gregory surplus transfers, Coombs, Bucklin, the
//     contingent and supplementary votes, and the two-round runoff;
//   - cardinal rules: approval, range (score), STAR, cumulative and
//     majority-judgment voting.
//
// The apportionment layer allocates a fixed number of seats to parties from
// their vote totals: the highest-averages (divisor) methods D'Hondt/Jefferson,
// Sainte-Laguë/Webster, Adams, Dean and Huntington-Hill, and the
// largest-remainder (quota) methods Hamilton and Droop, together with quota-rule
// and Alabama-paradox diagnostics.
//
// All routines use only the Go standard library, are deterministic, resolve
// ties toward the lowest candidate or party index unless documented otherwise,
// and are validated in the accompanying tests against reference values such as
// the standard Tennessee-capital voting example.
package socialchoice
