package auctions

import (
	"fmt"
	"math"
	"sort"
	"testing"
)

func approxScalar(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func approxVec(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

// ---------- sealed-bid auctions ----------

func TestSealedBidAuctions(t *testing.T) {
	bids := []Bid{{0, 5}, {1, 3}, {2, 8}, {3, 2}}
	tests := []struct {
		name       string
		run        func() (AuctionOutcome, error)
		wantWinner int
		wantPrice  float64
	}{
		{"first-price", func() (AuctionOutcome, error) { return FirstPriceAuction(bids, 0) }, 2, 8},
		{"second-price", func() (AuctionOutcome, error) { return SecondPriceAuction(bids, 0) }, 2, 5},
		{"vickrey", func() (AuctionOutcome, error) { return VickreyAuction(bids, 0) }, 2, 5},
		{"third-price", func() (AuctionOutcome, error) { return ThirdPriceAuction(bids, 0) }, 2, 3},
		{"second-price-reserve", func() (AuctionOutcome, error) { return SecondPriceAuction(bids, 6) }, 2, 6},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			out, err := tc.run()
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if out.Winner != tc.wantWinner {
				t.Fatalf("winner = %d, want %d", out.Winner, tc.wantWinner)
			}
			if !approxScalar(out.Price, tc.wantPrice, 1e-12) {
				t.Fatalf("price = %v, want %v", out.Price, tc.wantPrice)
			}
		})
	}
}

func TestReserveNotMet(t *testing.T) {
	bids := []Bid{{0, 5}, {1, 3}}
	out, err := FirstPriceAuction(bids, 10)
	if err != nil {
		t.Fatal(err)
	}
	if out.Winner != -1 {
		t.Fatalf("expected no sale, got winner %d", out.Winner)
	}
}

func TestAllPayAuction(t *testing.T) {
	bids := []Bid{{0, 5}, {1, 3}, {2, 8}}
	out, err := AllPayAuction(bids, 0)
	if err != nil {
		t.Fatal(err)
	}
	if out.Winner != 2 {
		t.Fatalf("winner = %d, want 2", out.Winner)
	}
	if !approxScalar(out.Revenue, 16, 1e-12) {
		t.Fatalf("revenue = %v, want 16", out.Revenue)
	}
}

func TestVickreyTruthfulness(t *testing.T) {
	// Bidder 2's utility (valuation 8) is 8 - price.
	bids := []Bid{{0, 5}, {1, 3}, {2, 8}}
	out, _ := SecondPriceAuction(bids, 0)
	u := BidderUtility(8, out, 2)
	if !approxScalar(u, 3, 1e-12) {
		t.Fatalf("utility = %v, want 3", u)
	}
}

// ---------- multi-unit ----------

func TestUniformAndVickreyMultiUnit(t *testing.T) {
	bids := []Bid{{0, 10}, {1, 8}, {2, 6}, {3, 4}}
	uni, err := UniformPriceAuction(bids, 2)
	if err != nil {
		t.Fatal(err)
	}
	// top 2 win, clearing price = 3rd highest = 6
	if len(uni.Winners) != 2 || uni.Winners[0] != 0 || uni.Winners[1] != 1 {
		t.Fatalf("winners = %v", uni.Winners)
	}
	for _, p := range uni.Payments {
		if !approxScalar(p, 6, 1e-12) {
			t.Fatalf("clearing price = %v, want 6", p)
		}
	}
	vick, _ := MultiUnitVickrey(bids, 2)
	for _, p := range vick.Payments {
		if !approxScalar(p, 6, 1e-12) {
			t.Fatalf("vickrey price = %v, want 6", p)
		}
	}
	// single unit reduces to second price
	one, _ := MultiUnitVickrey(bids, 1)
	if !approxScalar(one.Payments[0], 8, 1e-12) {
		t.Fatalf("single-unit vickrey price = %v, want 8", one.Payments[0])
	}
	pab, _ := PayAsBidAuction(bids, 2)
	if !approxScalar(pab.Revenue, 18, 1e-12) {
		t.Fatalf("pay-as-bid revenue = %v, want 18", pab.Revenue)
	}
}

// ---------- VCG / combinatorial ----------

func TestVCGSingleItemIsVickrey(t *testing.T) {
	bids := []Bid{{0, 5}, {1, 3}, {2, 8}}
	res, err := VCGMechanism(SingleItemToCombinatorial(bids), 1)
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(res.Welfare, 8, 1e-12) {
		t.Fatalf("welfare = %v, want 8", res.Welfare)
	}
	if !approxScalar(res.Payment[2], 5, 1e-12) {
		t.Fatalf("VCG payment = %v, want 5 (second price)", res.Payment[2])
	}
}

func TestCombinatorialVCG(t *testing.T) {
	// two items {0,1}. bidder0 wants {0,1}@10, bidder1 wants {0}@6, bidder2 wants {1}@5.
	bids := []CombinatorialBid{
		{Bidder: 0, Items: []int{0, 1}, Value: 10},
		{Bidder: 1, Items: []int{0}, Value: 6},
		{Bidder: 2, Items: []int{1}, Value: 5},
	}
	alloc, err := WinnerDetermination(bids, 2)
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(alloc.Value, 11, 1e-12) {
		t.Fatalf("welfare = %v, want 11", alloc.Value)
	}
	res, err := VCGMechanism(bids, 2)
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(res.Payment[1], 5, 1e-9) {
		t.Fatalf("payment[1] = %v, want 5", res.Payment[1])
	}
	if !approxScalar(res.Payment[2], 4, 1e-9) {
		t.Fatalf("payment[2] = %v, want 4", res.Payment[2])
	}
	if !approxScalar(res.Utility[1], 1, 1e-9) || !approxScalar(res.Utility[2], 1, 1e-9) {
		t.Fatalf("utilities = %v", res.Utility)
	}
}

// ---------- cooperative: Shapley / Banzhaf ----------

func gloveGame(t *testing.T) CoopGame {
	// v(12)=v(13)=v(123)=1? Use the game v(S)=1 if S⊇{1} and |S|>=2 style.
	// Simpler: v(1)=v(2)=v(3)=0, v(12)=1, v(13)=0, v(23)=0, v(123)=1.
	g, err := NewCoopGame(3, func(s Coalition) float64 {
		switch {
		case s == CoalitionFromMembers([]int{0, 1}):
			return 1
		case s == CoalitionFromMembers([]int{0, 1, 2}):
			return 1
		default:
			return 0
		}
	})
	if err != nil {
		t.Fatal(err)
	}
	return g
}

func TestShapleyClosedForm(t *testing.T) {
	// 2-player: phi1=(a+c-b)/2, phi2=(b+c-a)/2.
	a, b, c := 1.0, 2.0, 4.0
	g, _ := NewCoopGame(2, func(s Coalition) float64 {
		switch s {
		case 0b01:
			return a
		case 0b10:
			return b
		case 0b11:
			return c
		}
		return 0
	})
	phi := g.ShapleyValue()
	want := []float64{(a + c - b) / 2, (b + c - a) / 2}
	if !approxVec(phi, want, 1e-12) {
		t.Fatalf("shapley = %v, want %v", phi, want)
	}
}

func TestShapleyMonteCarloConverges(t *testing.T) {
	g := gloveGame(t)
	exact := g.ShapleyValue()
	approx := g.ShapleyValueMonteCarlo(200000, 42)
	if !approxVec(approx, exact, 5e-3) {
		t.Fatalf("monte carlo = %v, exact = %v", approx, exact)
	}
}

func TestWeightedVotingIndices(t *testing.T) {
	// [3; 2,1,1]
	w, err := NewWeightedVotingGame(3, []float64{2, 1, 1})
	if err != nil {
		t.Fatal(err)
	}
	ss := w.ShapleyShubikIndex()
	wantSS := []float64{2.0 / 3, 1.0 / 6, 1.0 / 6}
	if !approxVec(ss, wantSS, 1e-9) {
		t.Fatalf("shapley-shubik = %v, want %v", ss, wantSS)
	}
	bz := w.BanzhafIndex()
	wantBZ := []float64{3.0 / 5, 1.0 / 5, 1.0 / 5}
	if !approxVec(bz, wantBZ, 1e-9) {
		t.Fatalf("banzhaf = %v, want %v", bz, wantBZ)
	}
	if len(w.MinimalWinningCoalitions()) != 2 {
		t.Fatalf("expected 2 minimal winning coalitions")
	}
}

func TestDictatorAndDummy(t *testing.T) {
	// [3; 3,1,1]: player 0 is a dictator, players 1,2 dummies.
	w, _ := NewWeightedVotingGame(3, []float64{3, 1, 1})
	d, ok := w.Dictator()
	if !ok || d != 0 {
		t.Fatalf("dictator = %d,%v", d, ok)
	}
	if got := w.Dummies(); len(got) != 2 {
		t.Fatalf("dummies = %v", got)
	}
}

// ---------- core / least-core / nucleolus ----------

func TestCoreAndNucleolusGloveGame(t *testing.T) {
	g := gloveGame(t)
	// nucleolus = (0.5,0.5,0)
	nu, err := g.Nucleolus()
	if err != nil {
		t.Fatal(err)
	}
	if !approxVec(nu, []float64{0.5, 0.5, 0}, 1e-6) {
		t.Fatalf("nucleolus = %v, want [0.5 0.5 0]", nu)
	}
	if !g.InCore(nu, 1e-6) {
		t.Fatalf("nucleolus should be in the (non-empty) core")
	}
	if !g.CoreIsNonEmpty() {
		t.Fatalf("core should be non-empty")
	}
}

func TestLeastCoreEmptyCore(t *testing.T) {
	// simple majority: v(single)=0, v(pair)=2, v(N)=2 -> core empty, least core 2/3.
	g, _ := NewCoopGame(3, func(s Coalition) float64 {
		switch s.Size() {
		case 2:
			return 2
		case 3:
			return 2
		default:
			return 0
		}
	})
	if g.CoreIsNonEmpty() {
		t.Fatalf("core should be empty")
	}
	eps, x, ok := g.LeastCore()
	if !ok {
		t.Fatal("least core failed")
	}
	if !approxScalar(eps, 2.0/3, 1e-6) {
		t.Fatalf("least-core value = %v, want 2/3", eps)
	}
	nu, err := g.Nucleolus()
	if err != nil {
		t.Fatal(err)
	}
	if !approxVec(nu, []float64{2.0 / 3, 2.0 / 3, 2.0 / 3}, 1e-6) {
		t.Fatalf("nucleolus = %v, want [2/3 2/3 2/3]", nu)
	}
	_ = x
}

func TestConvexGameShapleyInCore(t *testing.T) {
	// convex game: v(S) = |S|^2.
	g, _ := NewCoopGame(3, func(s Coalition) float64 {
		n := float64(s.Size())
		return n * n
	})
	if !g.IsConvex() {
		t.Fatalf("game should be convex")
	}
	if !g.IsSuperadditive() {
		t.Fatalf("convex game should be superadditive")
	}
	phi := g.ShapleyValue()
	if !g.InCore(phi, 1e-9) {
		t.Fatalf("shapley value of a convex game must lie in the core, got %v", phi)
	}
	// symmetry: all players equal by symmetry -> each 3.
	if !approxVec(phi, []float64{3, 3, 3}, 1e-9) {
		t.Fatalf("shapley = %v, want [3 3 3]", phi)
	}
}

func TestTauValueConvexGame(t *testing.T) {
	g, _ := NewCoopGame(3, func(s Coalition) float64 {
		n := float64(s.Size())
		return n * n
	})
	tau, err := g.TauValue()
	if err != nil {
		t.Fatal(err)
	}
	if !approxVec(tau, []float64{3, 3, 3}, 1e-9) {
		t.Fatalf("tau = %v, want [3 3 3]", tau)
	}
	// utopia payoff M_i = 9 - 4 = 5.
	if !approxVec(g.UtopiaPayoff(), []float64{5, 5, 5}, 1e-9) {
		t.Fatalf("utopia = %v", g.UtopiaPayoff())
	}
}

func TestPrenucleolusMatchesNucleolus(t *testing.T) {
	g := gloveGame(t)
	pre, err := g.Prenucleolus()
	if err != nil {
		t.Fatal(err)
	}
	if !approxVec(pre, []float64{0.5, 0.5, 0}, 1e-6) {
		t.Fatalf("prenucleolus = %v", pre)
	}
}

func TestEpsilonCoreAndAdditive(t *testing.T) {
	// additive game v(S) = sum of members (members worth their index+1).
	w := []float64{1, 2, 3}
	g, _ := NewCoopGame(3, func(s Coalition) float64 {
		var sum float64
		for i := 0; i < 3; i++ {
			if s.Contains(i) {
				sum += w[i]
			}
		}
		return sum
	})
	if !g.IsAdditive() {
		t.Fatalf("game should be additive")
	}
	x := []float64{1, 2, 3}
	if !g.InCore(x, 1e-9) {
		t.Fatalf("marginal allocation should be in core of additive game")
	}
	if !g.InEpsilonCore(x, 0, 1e-9) {
		t.Fatalf("should be in 0-core")
	}
}

func TestSortedExcesses(t *testing.T) {
	g := gloveGame(t)
	ex := g.SortedExcesses([]float64{0.5, 0.5, 0})
	for i := 1; i < len(ex); i++ {
		if ex[i] > ex[i-1]+1e-12 {
			t.Fatalf("excesses not sorted descending: %v", ex)
		}
	}
}

// ---------- bargaining ----------

func TestNashBargainingSplitDollar(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {0, 1}}
	sol, err := NashBargainingSolution(pts, Point{0, 0})
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(sol.X, 0.5, 1e-9) || !approxScalar(sol.Y, 0.5, 1e-9) {
		t.Fatalf("nash = %v, want (0.5,0.5)", sol)
	}
}

func TestNashBargainingAsymmetric(t *testing.T) {
	pts := []Point{{0, 0}, {2, 0}, {0, 1}}
	sol, err := NashBargainingSolution(pts, Point{0, 0})
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(sol.X, 1, 1e-9) || !approxScalar(sol.Y, 0.5, 1e-9) {
		t.Fatalf("nash = %v, want (1,0.5)", sol)
	}
	ks, err := KalaiSmorodinskySolution(pts, Point{0, 0})
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(ks.X, 1, 1e-9) || !approxScalar(ks.Y, 0.5, 1e-9) {
		t.Fatalf("ks = %v, want (1,0.5)", ks)
	}
	util, err := UtilitarianSolution(pts)
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(util.X, 2, 1e-9) || !approxScalar(util.Y, 0, 1e-9) {
		t.Fatalf("utilitarian = %v, want (2,0)", util)
	}
}

func TestWeightedNashAndFeasibility(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {0, 1}}
	// symmetric weights recover the split-the-dollar midpoint.
	sol, err := WeightedNashBargainingSolution(pts, Point{0, 0}, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(sol.X, 0.5, 1e-6) || !approxScalar(sol.Y, 0.5, 1e-6) {
		t.Fatalf("weighted nash (0.5) = %v, want (0.5,0.5)", sol)
	}
	// heavier weight on player 0 shifts the split toward X.
	sol2, _ := WeightedNashBargainingSolution(pts, Point{0, 0}, 0.75)
	if !approxScalar(sol2.X, 0.75, 1e-6) || !approxScalar(sol2.Y, 0.25, 1e-6) {
		t.Fatalf("weighted nash (0.75) = %v, want (0.75,0.25)", sol2)
	}
	if !FeasibleContains(pts, Point{0.25, 0.25}) {
		t.Fatalf("interior point should be feasible")
	}
	if FeasibleContains(pts, Point{1, 1}) {
		t.Fatalf("(1,1) should be infeasible")
	}
}

func TestEgalitarianSolution(t *testing.T) {
	// unit square feasible set; egalitarian moves along (1,1) to (1,1).
	pts := []Point{{0, 0}, {1, 0}, {0, 1}, {1, 1}}
	sol, err := EgalitarianSolution(pts, Point{0, 0})
	if err != nil {
		t.Fatal(err)
	}
	if !approxScalar(sol.X, 1, 1e-9) || !approxScalar(sol.Y, 1, 1e-9) {
		t.Fatalf("egalitarian = %v, want (1,1)", sol)
	}
}

// ---------- matching ----------

func TestGaleShapley(t *testing.T) {
	// Men A(0),B(1); Women X(0),Y(1).
	menPrefs := [][]int{{0, 1}, {0, 1}}   // both prefer X then Y
	womenPrefs := [][]int{{1, 0}, {0, 1}} // X prefers B then A; Y prefers A then B
	m, err := GaleShapley(menPrefs, womenPrefs)
	if err != nil {
		t.Fatal(err)
	}
	// Expected man-optimal: A-Y, B-X.
	if m.ManToWoman[0] != 1 || m.ManToWoman[1] != 0 {
		t.Fatalf("matching = %v, want A-Y,B-X", m.ManToWoman)
	}
	if !IsStableMatching(menPrefs, womenPrefs, m) {
		t.Fatalf("matching should be stable")
	}
	// An unstable matching A-X, B-Y has a blocking pair.
	bad := Matching{ManToWoman: []int{0, 1}, WomanToMan: []int{0, 1}}
	if IsStableMatching(menPrefs, womenPrefs, bad) {
		t.Fatalf("A-X,B-Y should be unstable")
	}
	if CountBlockingPairs(menPrefs, womenPrefs, bad) == 0 {
		t.Fatalf("expected blocking pairs")
	}
}

func TestTopTradingCycles(t *testing.T) {
	// agents 0,1 want each other's house; agent 2 keeps its own.
	prefs := [][]int{
		{1, 2, 0},
		{0, 2, 1},
		{2, 0, 1},
	}
	got, err := TopTradingCycles(prefs)
	if err != nil {
		t.Fatal(err)
	}
	want := []int{1, 0, 2}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("ttc = %v, want %v", got, want)
		}
	}
}

func TestSerialDictatorship(t *testing.T) {
	prefs := [][]int{
		{0, 1, 2},
		{0, 1, 2},
		{0, 1, 2},
	}
	got, err := SerialDictatorship(prefs, []int{2, 0, 1})
	if err != nil {
		t.Fatal(err)
	}
	// order 2,0,1 -> agent2 takes 0, agent0 takes 1, agent1 takes 2.
	want := []int{1, 2, 0}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("serial dictatorship = %v, want %v", got, want)
		}
	}
}

// ---------- coalition helpers ----------

func TestCoalitionOps(t *testing.T) {
	c := CoalitionFromMembers([]int{0, 2, 3})
	if c.Size() != 3 {
		t.Fatalf("size = %d", c.Size())
	}
	if !c.Contains(2) || c.Contains(1) {
		t.Fatalf("membership wrong")
	}
	members := c.Members(4)
	sort.Ints(members)
	if fmt.Sprint(members) != "[0 2 3]" {
		t.Fatalf("members = %v", members)
	}
	if c.Complement(4) != SingletonCoalition(1) {
		t.Fatalf("complement wrong")
	}
}

// ---------- example ----------

func ExampleVickreyAuction() {
	bids := []Bid{{Bidder: 0, Value: 5}, {Bidder: 1, Value: 3}, {Bidder: 2, Value: 8}}
	out, _ := VickreyAuction(bids, 0)
	fmt.Printf("winner %d pays %.0f\n", out.Winner, out.Price)
	// Output: winner 2 pays 5
}

func ExampleCoopGame_ShapleyValue() {
	g, _ := NewCoopGame(3, func(s Coalition) float64 {
		n := float64(s.Size())
		return n * n
	})
	fmt.Println(g.ShapleyValue())
	// Output: [3 3 3]
}
