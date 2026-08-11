# queueing — provides queueing-theory models and formulas built entirely on the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/queueing.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/queueing)

Package queueing provides queueing-theory models and formulas built entirely
on the Go standard library.

The package covers the classical single- and multi-server Markovian queues
together with several non-Markovian and network extensions.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/queueing"
```

## Usage

This is the package's own `ExampleBirthDeath`, so it compiles and its output is
asserted on every `go test ./queueing/`.

```go
c, err := queueing.NewBirthDeath([]float64{1, 1, 1}, []float64{2, 2, 2})
	if err != nil {
		panic(err)
	}
	for n, p := range c.StationaryDistribution() {
		fmt.Printf("p%d = %.6f\n", n, p)
	}
	fmt.Printf("mean = %.6f\n", c.MeanNumber())
```

```
p0 = 0.533333
p1 = 0.266667
p2 = 0.133333
p3 = 0.066667
mean = 0.733333
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AllenCunneenLq(lambda, mu float64, c int, ca2, cs2 float64) float64` | AllenCunneenLq returns the Allen–Cunneen approximation of the mean queue length in a G/G/c queue, lambda times `AllenCunneenWq`. |
| `func AllenCunneenWq(lambda, mu float64, c int, ca2, cs2 float64) float64` | AllenCunneenWq returns the Allen–Cunneen approximation of the mean waiting time in a G/G/c queue with arrival rate lambda, per-server service rate… |
| `func BirthDeathStationary(birth, death []float64) []float64` | BirthDeathStationary returns the stationary distribution of a birth–death chain given standalone birth and death rate slices, or nil when the… |
| `func CarriedLoad(c int, a float64) float64` | CarriedLoad returns the load carried by an M/M/c/c loss system, a(1-B), where B is the Erlang-B blocking probability at offered load a. |
| `func ErlangB(c int, a float64) float64` | ErlangB returns the Erlang-B loss probability (blocking probability) of an M/M/c/c loss system with c servers and offered load a erlangs. |
| `func ErlangBSeries(c int, a float64) []float64` | ErlangBSeries returns the Erlang-B blocking probabilities for all server counts from 0 through c at offered load a, produced by a single sweep of the… |
| `func ErlangBServersFor(a, gos float64, maxServers int) int` | ErlangBServersFor returns the smallest number of servers whose Erlang-B blocking probability at offered load a is at or below the target grade of… |
| `func ErlangC(c int, a float64) float64` | ErlangC returns the Erlang-C delay probability, the probability that an arriving customer must wait in an M/M/c queue with c servers and offered load… |
| `func ErlangCDF(k int, lambda, x float64) float64` | ErlangCDF returns the cumulative distribution of an Erlang-k distribution at x>=0, 1 - e^{-lambda x} sum_{i=0}^{k-1} (lambda x)^i / i!. |
| `func ErlangCFromB(c int, a, b float64) float64` | ErlangCFromB returns the Erlang-C delay probability computed from a precomputed Erlang-B value b for c servers and offered load a. |
| `func ErlangCMeanWait(c int, lambda, mu float64) float64` | ErlangCMeanWait returns the mean queueing delay in an M/M/c queue with c servers, arrival rate lambda and per-server service rate mu, computed from… |
| `func ErlangCServersFor(a, p float64, maxServers int) int` | ErlangCServersFor returns the smallest number of servers whose Erlang-C delay probability at offered load a is at or below the target probability p. |
| `func ErlangCWaitTailProb(c int, lambda, mu, t float64) float64` | ErlangCWaitTailProb returns the probability that the queueing delay in an M/M/c queue exceeds t, C * exp(-(c*mu - lambda) t), where C is the Erlang-C… |
| `func ErlangMean(k int, lambda float64) float64` | ErlangMean returns the mean k/lambda of an Erlang-k distribution. |
| `func ErlangPDF(k int, lambda, x float64) float64` | ErlangPDF returns the density of an Erlang-k distribution (sum of k iid exponentials of rate lambda) at x>=0. |
| `func ExponentialCDF(rate, x float64) float64` | ExponentialCDF returns the exponential cumulative probability 1-exp(-rate*x) for x>=0. |
| `func ExponentialMean(rate float64) float64` | ExponentialMean returns the mean 1/rate of an exponential distribution. |
| `func ExponentialPDF(rate, x float64) float64` | ExponentialPDF returns the exponential probability density rate*exp(-rate*x) for x>=0. |
| `func ExponentialTail(rate, t float64) float64` | ExponentialTail returns the exponential survival probability P(X>t)=exp(-rate*t) for t>=0. |
| `func Factorial(n int) float64` | Factorial returns n! as a float64. |
| `func KingmanWq(rho, ca2, cs2, serviceMean float64) float64` | KingmanWq returns the Kingman heavy-traffic approximation of the mean G/G/1 waiting time given the utilization rho, the interarrival and service SCVs… |
| `func LittleL(lambda, w float64) float64` | LittleL returns the mean number in system L = lambda*W from the arrival rate and mean sojourn time via Little's law. |
| `func LittleLambda(l, w float64) float64` | LittleLambda returns the arrival rate lambda = L/W implied by Little's law from the mean number in system and mean sojourn time. |
| `func LittleLq(lambda, wq float64) float64` | LittleLq returns the mean queue length Lq = lambda*Wq from the arrival rate and mean waiting time via Little's law applied to the queue. |
| `func LittleW(l, lambda float64) float64` | LittleW returns the mean sojourn time W = L/lambda from the mean number in system and the arrival rate via Little's law. |
| `func LittleWq(lq, lambda float64) float64` | LittleWq returns the mean waiting time in queue Wq = Lq/lambda from the mean queue length and arrival rate via Little's law. |
| `func LogFactorial(n int) float64` | LogFactorial returns the natural logarithm of n! using the log-gamma function, remaining finite for large n. |
| `func MD1Lq(lambda, mu float64) float64` | MD1Lq returns the M/D/1 mean queue length rho^2/(2(1-rho)). |
| `func MD1Wq(lambda, mu float64) float64` | MD1Wq returns the M/D/1 mean waiting time. |
| `func MG1Lq(lambda, serviceMean, serviceVar float64) float64` | MG1Lq returns the M/G/1 Pollaczek–Khinchine mean queue length from the arrival rate, mean service time and service-time variance. |
| `func MG1Wq(lambda, serviceMean, serviceVar float64) float64` | MG1Wq returns the M/G/1 mean waiting time from the arrival rate, mean service time and service-time variance. |
| `func MM1KBlockingProb(lambda, mu float64, k int) float64` | MM1KBlockingProb returns the M/M/1/K loss probability P(N=K). |
| `func MM1KL(lambda, mu float64, k int) float64` | MM1KL returns the mean number in an M/M/1/K system. |
| `func MM1L(lambda, mu float64) float64` | MM1L returns the mean number in an M/M/1 system, rho/(1-rho). |
| `func MM1Lq(lambda, mu float64) float64` | MM1Lq returns the mean M/M/1 queue length rho^2/(1-rho). |
| `func MM1P0(lambda, mu float64) float64` | MM1P0 returns the M/M/1 empty-system probability 1-rho. |
| `func MM1Pn(lambda, mu float64, n int) float64` | MM1Pn returns the M/M/1 state probability (1-rho)rho^n. |
| `func MM1Rho(lambda, mu float64) float64` | MM1Rho returns the M/M/1 utilization lambda/mu. |
| `func MM1W(lambda, mu float64) float64` | MM1W returns the mean M/M/1 sojourn time 1/(mu-lambda). |
| `func MM1WaitTailProb(lambda, mu, t float64) float64` | MM1WaitTailProb returns P(W>t) for an M/M/1 queue. |
| `func MM1Wq(lambda, mu float64) float64` | MM1Wq returns the mean M/M/1 waiting time lambda/(mu(mu-lambda)). |
| `func MMInfL(lambda, mu float64) float64` | MMInfL returns the mean number in an M/M/infinity system, lambda/mu. |
| `func MMInfPn(lambda, mu float64, n int) float64` | MMInfPn returns the M/M/infinity state probability (Poisson mass). |
| `func MMcDelayProb(lambda, mu float64, c int) float64` | MMcDelayProb returns the Erlang-C delay probability for an M/M/c queue. |
| `func MMcKBlockingProb(lambda, mu float64, c, k int) float64` | MMcKBlockingProb returns the M/M/c/K loss probability P(N=K). |
| `func MMcKL(lambda, mu float64, c, k int) float64` | MMcKL returns the mean number in an M/M/c/K system. |
| `func MMcL(lambda, mu float64, c int) float64` | MMcL returns the mean number in an M/M/c system. |
| `func MMcLq(lambda, mu float64, c int) float64` | MMcLq returns the mean M/M/c queue length. |
| `func MMcP0(lambda, mu float64, c int) float64` | MMcP0 returns the M/M/c empty-system probability. |
| `func MMcPn(lambda, mu float64, c, n int) float64` | MMcPn returns the M/M/c state probability P(N=n). |
| `func MMcW(lambda, mu float64, c int) float64` | MMcW returns the mean M/M/c sojourn time. |
| `func MMcWq(lambda, mu float64, c int) float64` | MMcWq returns the mean M/M/c waiting time. |
| `func OfferedLoad(lambda, mu float64) float64` | OfferedLoad returns the offered load (traffic intensity in erlangs) a = lambda/mu, the mean number of busy servers in an unrestricted system. |
| `func OfferedLoadFromCarried(carried, blocking float64) float64` | OfferedLoadFromCarried returns the offered load implied by a measured carried load and blocking probability, carried/(1-blocking). |
| `func PoissonCDF(k int, lambda float64) float64` | PoissonCDF returns the Poisson cumulative probability P(N<=k) with mean lambda. |
| `func PoissonPMF(k int, lambda float64) float64` | PoissonPMF returns the Poisson probability mass P(N=k) with mean lambda, computed in log space for numerical stability. |
| `func PollaczekKhinchineLq(rho, scv float64) float64` | PollaczekKhinchineLq returns the Pollaczek–Khinchine mean queue length for an M/G/1 queue with utilization rho and service-time squared coefficient… |
| `func PollaczekKhinchineWq(rho, scv, serviceMean float64) float64` | PollaczekKhinchineWq returns the Pollaczek–Khinchine mean waiting time rho*(1+scv)/(2(1-rho)) * E`S`, where E`S` is the mean service time. |
| `func SojournFromWait(wq, mu float64) float64` | SojournFromWait returns the mean sojourn time W = Wq + 1/mu obtained by adding a mean service time to the mean waiting time. |
| `func SquaredCoefficientOfVariation(mean, variance float64) float64` | SquaredCoefficientOfVariation returns the squared coefficient of variation (SCV) variance/mean^2 of a distribution given its mean and variance. |
| `func TrafficFromMetrics(l, w float64) float64` | TrafficFromMetrics returns the arrival rate implied by a mean number in system L and a mean sojourn time W via Little's law, an alias for… |
| `func Utilization(lambda, mu float64, c int) float64` | Utilization returns the per-server utilization rho = lambda/(c*mu) of a c-server queue. |
| `func WaitFromSojourn(w, mu float64) float64` | WaitFromSojourn returns the mean waiting time Wq = W - 1/mu obtained by subtracting the mean service time from the mean sojourn time. |

### Types

| Type | What it is |
| --- | --- |
| `BirthDeath` | BirthDeath models a finite continuous-time birth–death chain on states 0,1,...,N. Birth`i` is the transition rate from state i to i+1 (for… |
| `GG1` | GG1 models a general single-server queue with renewal arrivals and general service (Kendall notation G/G/1). |
| `JacksonNetwork` | JacksonNetwork models an open Jackson network of M/M/c queues. |
| `MD1` | MD1 models a single-server queue with Poisson arrivals of rate Lambda and deterministic (constant) service time 1/Mu (Kendall notation M/D/1). |
| `MG1` | MG1 models a single-server queue with Poisson arrivals of rate Lambda and a general service-time distribution specified by its mean ServiceMean and… |
| `MM1` | MM1 models a single-server queue with Poisson arrivals of rate Lambda, exponential service of rate Mu, an infinite waiting room and FCFS discipline… |
| `MM1K` | MM1K models a single-server queue with Poisson arrivals of rate Lambda, exponential service of rate Mu and a finite system capacity K (including the… |
| `MMInf` | MMInf models an infinite-server queue with Poisson arrivals of rate Lambda and exponential service of rate Mu (Kendall notation M/M/∞). |
| `MMc` | MMc models a multi-server queue with Poisson arrivals of rate Lambda, C identical exponential servers each of rate Mu, an infinite waiting room and… |
| `MMcK` | MMcK models a multi-server queue with Poisson arrivals of rate Lambda, C exponential servers of rate Mu each and a finite system capacity K>=C.… |
| `PriorityMM1` | PriorityMM1 models a single-server queue with Poisson arrivals in several priority classes and exponential, class-dependent service. |

<details>
<summary><code>BirthDeath</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBirthDeath(birth, death []float64) (BirthDeath, error)` | NewBirthDeath constructs a `BirthDeath` chain from equal-length birth and death rate slices. |
| `func (c BirthDeath) DetailedBalance(p []float64, tol float64) bool` | DetailedBalance reports whether a candidate distribution p satisfies the detailed-balance equations p`i`*Birth`i` == p[i+1]*Death`i` within the given… |
| `func (c BirthDeath) MeanNumber() float64` | MeanNumber returns the stationary mean state value sum_n n*p_n, i.e. |
| `func (c BirthDeath) ProbN(n int) float64` | ProbN returns the stationary probability of state n, or 0 when n is outside the range 0..N. |
| `func (c BirthDeath) States() int` | States returns the number of states N+1 in the chain. |
| `func (c BirthDeath) StationaryDistribution() []float64` | StationaryDistribution returns the normalized stationary distribution over states 0..N obtained from the detailed-balance product p_n = p_0 *… |
| `func (c BirthDeath) Throughput() float64` | Throughput returns the stationary birth rate sum_n Birth`n`*p_n, which equals the long-run rate of upward transitions. |
| `func (c BirthDeath) VarianceNumber() float64` | VarianceNumber returns the variance of the stationary state value. |

</details>

<details>
<summary><code>GG1</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGG1(lambda, serviceMean, ca2, cs2 float64) (GG1, error)` | NewGG1 constructs a `GG1` queue from the arrival rate, mean service time and the squared coefficients of variation of the interarrival and service… |
| `func (q GG1) L() float64` | L returns the approximate mean number in system, Lq + rho. |
| `func (q GG1) Lq() float64` | Lq returns the approximate mean queue length, Lambda*Wq. |
| `func (q GG1) Rho() float64` | Rho returns the utilization rho = Lambda*E`S`. |
| `func (q GG1) W() float64` | W returns the approximate mean sojourn time, Wq + E`S`. |
| `func (q GG1) Wq() float64` | Wq returns the Kingman (VUT) approximation of the mean waiting time, (rho/(1-rho)) * ((Ca2+Cs2)/2) * E`S`. |

</details>

<details>
<summary><code>JacksonNetwork</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewJacksonNetwork(external []float64, routing [][]float64, serviceRate []float64, servers []int) (JacksonNetwork, error)` | NewJacksonNetwork constructs an open `JacksonNetwork`. |
| `func (net JacksonNetwork) ArrivalRates() []float64` | ArrivalRates solves the traffic equations lambda_j = gamma_j + sum_i lambda_i r_ij and returns the vector of per-node total arrival rates. |
| `func (net JacksonNetwork) MeanNumberAt(i int) float64` | MeanNumberAt returns the mean number of customers at node i, analyzing it as an independent M/M/c queue with the arrival rate from the traffic… |
| `func (net JacksonNetwork) MeanNumbers() []float64` | MeanNumbers returns the mean number of customers at every node. |
| `func (net JacksonNetwork) Nodes() int` | Nodes returns the number of nodes in the network. |
| `func (net JacksonNetwork) SojournTime() float64` | SojournTime returns the mean time a customer spends in the network, TotalMean/Throughput, by Little's law applied to the whole network. |
| `func (net JacksonNetwork) Stable() bool` | Stable reports whether every node has utilization strictly below one, the condition for the network to admit a product-form steady state. |
| `func (net JacksonNetwork) Throughput() float64` | Throughput returns the network throughput, the total external arrival rate (equal in steady state to the total departure rate). |
| `func (net JacksonNetwork) TotalMean() float64` | TotalMean returns the mean total number of customers in the network, summed over all nodes. |
| `func (net JacksonNetwork) Utilizations() []float64` | Utilizations returns the per-node utilization rho_i = lambda_i/(c_i mu_i). |

</details>

<details>
<summary><code>MD1</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMD1(lambda, mu float64) (MD1, error)` | NewMD1 constructs an `MD1` queue, returning an error for non-positive rates or when the stability condition Lambda < Mu is violated. |
| `func (q MD1) L() float64` | L returns the mean number in system, Lq + rho. |
| `func (q MD1) Lq() float64` | Lq returns the mean number waiting in the queue, rho^2/(2(1-rho)). |
| `func (q MD1) P0() float64` | P0 returns the probability the system is empty, 1-rho. |
| `func (q MD1) Rho() float64` | Rho returns the utilization rho = Lambda/Mu. |
| `func (q MD1) Stable() bool` | Stable reports whether Lambda < Mu. |
| `func (q MD1) W() float64` | W returns the mean sojourn time, Wq + 1/Mu. |
| `func (q MD1) Wq() float64` | Wq returns the mean waiting time in the queue, rho/(2 Mu (1-rho)). |

</details>

<details>
<summary><code>MG1</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMG1(lambda, serviceMean, serviceVar float64) (MG1, error)` | NewMG1 constructs an `MG1` queue from an arrival rate and the mean and variance of the service time. |
| `func NewMG1FromSCV(lambda, serviceMean, scv float64) (MG1, error)` | NewMG1FromSCV constructs an `MG1` queue from the arrival rate, mean service time and squared coefficient of variation scv = Var`S`/E`S`^2. |
| `func (q MG1) L() float64` | L returns the mean number in system, Lq + rho. |
| `func (q MG1) Lq() float64` | Lq returns the mean number waiting in the queue from the Pollaczek–Khinchine formula, (Lambda^2 Var`S` + rho^2) / (2(1-rho)). |
| `func (q MG1) MeanResidualService() float64` | MeanResidualService returns the mean residual (forward-recurrence) service time seen by an arrival at a busy instant, E[S^2]/(2 E`S`). |
| `func (q MG1) P0() float64` | P0 returns the probability that the system is empty, 1-rho, which holds for any M/G/1 queue. |
| `func (q MG1) Rho() float64` | Rho returns the utilization rho = Lambda*E`S`. |
| `func (q MG1) ServiceRate() float64` | ServiceRate returns the nominal service rate Mu = 1/E`S`. |
| `func (q MG1) ServiceSCV() float64` | ServiceSCV returns the squared coefficient of variation of the service time, Var`S`/E`S`^2. |
| `func (q MG1) Stable() bool` | Stable reports whether the utilization is below one. |
| `func (q MG1) W() float64` | W returns the mean sojourn time, Wq + E`S`. |
| `func (q MG1) Wq() float64` | Wq returns the mean waiting time in the queue, Lq/Lambda. |

</details>

<details>
<summary><code>MM1</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMM1(lambda, mu float64) (MM1, error)` | NewMM1 constructs an `MM1` queue, returning an error for non-positive rates or when the stability condition Lambda < Mu is violated. |
| `func (q MM1) IdleProb() float64` | IdleProb returns the probability the system is empty, P0 = 1 - rho. |
| `func (q MM1) L() float64` | L returns the mean number of customers in the system, rho/(1-rho). |
| `func (q MM1) Lq() float64` | Lq returns the mean number of customers waiting in the queue, rho^2/(1-rho). |
| `func (q MM1) Ls() float64` | Ls returns the mean number of customers in service, equal to rho. |
| `func (q MM1) MeanBusyCycle() float64` | MeanBusyCycle returns the mean length of a busy cycle (idle period plus busy period), 1/(Lambda(1-rho)). |
| `func (q MM1) MeanBusyPeriod() float64` | MeanBusyPeriod returns the mean length of a server busy period, 1/(Mu-Lambda). |
| `func (q MM1) P0() float64` | P0 returns the steady-state probability of an empty system, 1 - rho. |
| `func (q MM1) PAtLeast(n int) float64` | PAtLeast returns the probability of at least n customers in the system, P(N>=n) = rho^n. |
| `func (q MM1) Pn(n int) float64` | Pn returns the steady-state probability P(N=n) = (1-rho) rho^n of finding n customers in the system. |
| `func (q MM1) Rho() float64` | Rho returns the server utilization (traffic intensity) rho = Lambda/Mu, which also equals the long-run fraction of time the server is busy. |
| `func (q MM1) Stable() bool` | Stable reports whether the queue admits a steady state (Lambda < Mu). |
| `func (q MM1) Throughput() float64` | Throughput returns the departure rate of the system, which in steady state equals the arrival rate Lambda. |
| `func (q MM1) Utilization() float64` | Utilization returns the fraction of time the server is busy, equal to Rho. |
| `func (q MM1) VarianceN() float64` | VarianceN returns the variance of the number in system, rho/(1-rho)^2. |
| `func (q MM1) W() float64` | W returns the mean sojourn time (response time) in the system, 1/(Mu-Lambda). |
| `func (q MM1) WaitCDF(t float64) float64` | WaitCDF returns P(W<=t), the sojourn-time distribution 1-exp(-(Mu-Lambda)t). |
| `func (q MM1) WaitPercentile(p float64) float64` | WaitPercentile returns the sojourn-time quantile t such that P(W<=t)=p, i.e. |
| `func (q MM1) WaitTailProb(t float64) float64` | WaitTailProb returns P(W>t), the probability the sojourn time exceeds t, which is exp(-(Mu-Lambda)t) for the M/M/1 FCFS queue. |
| `func (q MM1) WaitqCDF(t float64) float64` | WaitqCDF returns P(Wq<=t) = 1 - rho*exp(-(Mu-Lambda)t), the queueing-delay distribution, which has an atom 1-rho at zero. |
| `func (q MM1) WaitqTailProb(t float64) float64` | WaitqTailProb returns P(Wq>t), the probability the queueing delay exceeds t, which is rho*exp(-(Mu-Lambda)t). |
| `func (q MM1) Wq() float64` | Wq returns the mean waiting time in the queue, Lambda/(Mu(Mu-Lambda)). |

</details>

<details>
<summary><code>MM1K</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMM1K(lambda, mu float64, k int) (MM1K, error)` | NewMM1K constructs an `MM1K` queue, returning an error for non-positive rates or a capacity below 1. |
| `func (q MM1K) BlockingProb() float64` | BlockingProb returns the probability P(N=K) that an arriving customer is blocked and lost. |
| `func (q MM1K) EffectiveArrivalRate() float64` | EffectiveArrivalRate returns the rate of accepted customers, Lambda*(1-BlockingProb). |
| `func (q MM1K) L() float64` | L returns the mean number of customers in the system. |
| `func (q MM1K) Lq() float64` | Lq returns the mean number of customers waiting in the queue, L-(1-P0). |
| `func (q MM1K) P0() float64` | P0 returns the steady-state probability that the system is empty. |
| `func (q MM1K) Pn(n int) float64` | Pn returns the steady-state probability P(N=n) for 0<=n<=K, and 0 outside that range. |
| `func (q MM1K) Rho() float64` | Rho returns the traffic intensity rho = Lambda/Mu, which may exceed one for a finite-capacity queue. |
| `func (q MM1K) Throughput() float64` | Throughput returns the departure rate, equal to the effective arrival rate. |
| `func (q MM1K) Utilization() float64` | Utilization returns the fraction of time the server is busy, 1-P0, which also equals EffectiveArrivalRate/Mu. |
| `func (q MM1K) W() float64` | W returns the mean sojourn time using the effective arrival rate, L/lambdaEff. |
| `func (q MM1K) Wq() float64` | Wq returns the mean waiting time in the queue, Lq/lambdaEff. |

</details>

<details>
<summary><code>MMInf</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMMInf(lambda, mu float64) (MMInf, error)` | NewMMInf constructs an `MMInf` queue, returning an error for non-positive rates. |
| `func (q MMInf) L() float64` | L returns the mean number of customers in the system, a = Lambda/Mu. |
| `func (q MMInf) Lq() float64` | Lq returns the mean number waiting in the queue, which is always 0. |
| `func (q MMInf) MeanBusyServers() float64` | MeanBusyServers returns the mean number of busy servers, equal to a. |
| `func (q MMInf) OfferedLoad() float64` | OfferedLoad returns the offered load a = Lambda/Mu, which equals the mean number in system. |
| `func (q MMInf) P0() float64` | P0 returns the probability that the system is empty, e^{-a}. |
| `func (q MMInf) Pn(n int) float64` | Pn returns the steady-state probability P(N=n) = e^{-a} a^n/n!, the Poisson mass with mean a = Lambda/Mu. |
| `func (q MMInf) VarianceN() float64` | VarianceN returns the variance of the number in system, equal to the mean a for the Poisson distribution. |
| `func (q MMInf) W() float64` | W returns the mean sojourn time, equal to a single service time 1/Mu. |
| `func (q MMInf) Wq() float64` | Wq returns the mean waiting time in the queue, which is always 0. |

</details>

<details>
<summary><code>MMc</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMMc(lambda, mu float64, c int) (MMc, error)` | NewMMc constructs an `MMc` queue, returning an error for non-positive rates, a non-positive server count, or when the stability condition Lambda <… |
| `func (q MMc) ErlangC() float64` | ErlangC returns the Erlang-C delay probability, the steady-state probability that an arriving customer must wait (finds all C servers busy). |
| `func (q MMc) L() float64` | L returns the mean number of customers in the system, Lq + a. |
| `func (q MMc) Lq() float64` | Lq returns the mean number of customers waiting in the queue, ErlangC * rho/(1-rho). |
| `func (q MMc) MeanBusyServers() float64` | MeanBusyServers returns the mean number of busy servers, equal to the offered load a = Lambda/Mu. |
| `func (q MMc) OfferedLoad() float64` | OfferedLoad returns the offered load a = Lambda/Mu in erlangs, the mean number of busy servers in steady state. |
| `func (q MMc) P0() float64` | P0 returns the steady-state probability that the system is empty. |
| `func (q MMc) Pn(n int) float64` | Pn returns the steady-state probability P(N=n) of finding n customers in the system. |
| `func (q MMc) ProbWait() float64` | ProbWait returns the probability that an arriving customer is delayed, equal to the Erlang-C value. |
| `func (q MMc) Rho() float64` | Rho returns the per-server utilization rho = Lambda/(C*Mu). |
| `func (q MMc) Stable() bool` | Stable reports whether the queue admits a steady state (Lambda < C*Mu). |
| `func (q MMc) Throughput() float64` | Throughput returns the departure rate, which equals Lambda in steady state. |
| `func (q MMc) Utilization() float64` | Utilization returns the per-server utilization, equal to Rho. |
| `func (q MMc) W() float64` | W returns the mean sojourn time in the system, Wq + 1/Mu. |
| `func (q MMc) WaitqCDF(t float64) float64` | WaitqCDF returns P(Wq<=t), the queueing-delay distribution. |
| `func (q MMc) WaitqPercentile(p float64) float64` | WaitqPercentile returns the queueing-delay quantile t such that P(Wq<=t)=p. |
| `func (q MMc) WaitqTailProb(t float64) float64` | WaitqTailProb returns P(Wq>t), the probability the queueing delay exceeds t, which is ErlangC * exp(-(C*Mu-Lambda) t). |
| `func (q MMc) Wq() float64` | Wq returns the mean waiting time in the queue, Lq/Lambda. |

</details>

<details>
<summary><code>MMcK</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMMcK(lambda, mu float64, c, k int) (MMcK, error)` | NewMMcK constructs an `MMcK` queue, returning an error for non-positive rates, a non-positive server count or a capacity smaller than the number of… |
| `func (q MMcK) BlockingProb() float64` | BlockingProb returns the loss probability P(N=K). |
| `func (q MMcK) EffectiveArrivalRate() float64` | EffectiveArrivalRate returns the accepted-customer rate Lambda*(1-BlockingProb). |
| `func (q MMcK) L() float64` | L returns the mean number of customers in the system. |
| `func (q MMcK) Lq() float64` | Lq returns the mean number of customers waiting in the queue. |
| `func (q MMcK) MeanBusyServers() float64` | MeanBusyServers returns the mean number of busy servers, which equals EffectiveArrivalRate/Mu. |
| `func (q MMcK) OfferedLoad() float64` | OfferedLoad returns the offered load a = Lambda/Mu in erlangs. |
| `func (q MMcK) P0() float64` | P0 returns the steady-state probability that the system is empty. |
| `func (q MMcK) Pn(n int) float64` | Pn returns the steady-state probability P(N=n) for 0<=n<=K, and 0 otherwise. |
| `func (q MMcK) ProbWait() float64` | ProbWait returns the probability that an accepted arrival must wait, i.e. |
| `func (q MMcK) Rho() float64` | Rho returns the nominal per-server traffic intensity Lambda/(C*Mu). |
| `func (q MMcK) Throughput() float64` | Throughput returns the departure rate, equal to the effective arrival rate. |
| `func (q MMcK) Utilization() float64` | Utilization returns the mean fraction of servers busy, MeanBusyServers / C. |
| `func (q MMcK) W() float64` | W returns the mean sojourn time L/lambdaEff. |
| `func (q MMcK) Wq() float64` | Wq returns the mean waiting time in the queue Lq/lambdaEff. |

</details>

<details>
<summary><code>PriorityMM1</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPriorityMM1(lambda, mu []float64) (PriorityMM1, error)` | NewPriorityMM1 constructs a `PriorityMM1` queue from per-class arrival and service rates. |
| `func (q PriorityMM1) ClassL(k int) float64` | ClassL returns the mean number of class-k customers in system, Lambda`k`*W_k. |
| `func (q PriorityMM1) ClassLq(k int) float64` | ClassLq returns the mean number of class-k customers waiting, Lambda`k`*Wq_k. |
| `func (q PriorityMM1) ClassRho(k int) float64` | ClassRho returns the per-class utilization rho_k = Lambda`k`/Mu`k`. |
| `func (q PriorityMM1) ClassW(k int) float64` | ClassW returns the mean sojourn time of class k, ClassWq(k)+1/Mu`k`. |
| `func (q PriorityMM1) ClassWq(k int) float64` | ClassWq returns the mean waiting time of class k under non-preemptive priority, R / ((1-sigma_{k-1})(1-sigma_k)), where sigma_j is the cumulative… |
| `func (q PriorityMM1) Classes() int` | Classes returns the number of priority classes. |
| `func (q PriorityMM1) MeanWait() float64` | MeanWait returns the arrival-weighted overall mean waiting time across all classes, (sum_k Lambda`k` Wq_k)/(sum_k Lambda`k`). |
| `func (q PriorityMM1) Rho() float64` | Rho returns the total utilization sum_k Lambda`k`/Mu`k`. |

</details>

### Variables

`ErrNonPositiveRate`, `ErrServers`, `ErrCapacity`, `ErrUnstable`, `ErrDimension`, `ErrNegative`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/queueing).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
