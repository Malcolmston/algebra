// Package queueing provides queueing-theory models and formulas built entirely
// on the Go standard library.
//
// The package covers the classical single- and multi-server Markovian queues
// together with several non-Markovian and network extensions.
//
// Markovian queues: steady-state analysis of the M/M/1, M/M/c, M/M/1/K,
// M/M/c/K and M/M/∞ systems is provided through the [MM1], [MMc], [MM1K],
// [MMcK] and [MMInf] types. Each type exposes utilization, the stationary
// state distribution, mean number in system and in queue, mean sojourn and
// waiting times, blocking probabilities for the finite-capacity variants and
// waiting-time tail probabilities where a closed form exists.
//
// Non-Markovian queues: the [MG1] type implements the Pollaczek–Khinchine
// mean-value results for the M/G/1 queue given the mean and variance (or
// squared coefficient of variation) of a general service distribution, with
// [MD1] specializing to deterministic service. The [GG1] type collects the
// Kingman and Allen–Cunneen diffusion approximations for the general G/G/1 and
// G/G/c queues.
//
// Loss and delay systems: [ErlangB] and [ErlangC] evaluate the Erlang loss and
// delay formulas via numerically stable recurrences, with helpers for carried
// load, server sizing and the Erlang-C waiting-time distribution.
//
// Birth–death chains and networks: [BirthDeath] solves a general finite
// birth–death chain for its stationary distribution and moments, while
// [JacksonNetwork] solves the traffic equations of an open Jackson network and
// reports per-node utilization and occupancy together with network-wide
// throughput and sojourn time.
//
// Little's law and related identities are available as free functions such as
// [LittleL] and [LittleW].
//
// Unless stated otherwise, arrival and service rates must be strictly positive
// and probabilities lie in [0,1]; functions return NaN (or an error from a
// constructor) when arguments are out of range or when a stability condition
// is violated.
package queueing
