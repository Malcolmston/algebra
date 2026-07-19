package proofsystems

import "errors"

// ErrUnify is returned when two terms, literals, or lists have no unifier.
var ErrUnify = errors.New("proofsystems: terms are not unifiable")

// UnifyTerms computes a most-general unifier of two terms using Robinson's
// algorithm with the occurs-check. On success the returned substitution s
// satisfies s.ApplyTerm(a).Equal(s.ApplyTerm(b)); on failure it returns
// ErrUnify.
func UnifyTerms(a, b Term) (Substitution, error) {
	return unify(a, b, NewSubstitution())
}

// MostGeneralUnifier is a synonym for UnifyTerms returning the MGU of two
// terms.
func MostGeneralUnifier(a, b Term) (Substitution, error) {
	return UnifyTerms(a, b)
}

// IsUnifiable reports whether two terms have a unifier.
func IsUnifiable(a, b Term) bool {
	_, err := UnifyTerms(a, b)
	return err == nil
}

// UnifyTermLists unifies two equal-length lists of terms position by position,
// threading the accumulated substitution. Lists of differing length are not
// unifiable.
func UnifyTermLists(a, b []Term) (Substitution, error) {
	if len(a) != len(b) {
		return Substitution{}, ErrUnify
	}
	s := NewSubstitution()
	var err error
	for i := range a {
		s, err = unify(a[i], b[i], s)
		if err != nil {
			return Substitution{}, err
		}
	}
	return s, nil
}

func unify(a, b Term, s Substitution) (Substitution, error) {
	a = s.ApplyTerm(a)
	b = s.ApplyTerm(b)
	switch {
	case a.IsVar():
		return unifyVar(a.Name, b, s)
	case b.IsVar():
		return unifyVar(b.Name, a, s)
	default:
		if a.Name != b.Name || len(a.Args) != len(b.Args) {
			return Substitution{}, ErrUnify
		}
		var err error
		for i := range a.Args {
			s, err = unify(a.Args[i], b.Args[i], s)
			if err != nil {
				return Substitution{}, err
			}
		}
		return s, nil
	}
}

func unifyVar(v string, t Term, s Substitution) (Substitution, error) {
	if t.IsVar() && t.Name == v {
		return s, nil
	}
	if t.Occurs(v) {
		return Substitution{}, ErrUnify
	}
	// Compose the new binding into all existing images to keep s resolved.
	single := SingletonSubstitution(v, t)
	return s.Compose(single).Bind(v, t), nil
}

// MatchTerm computes a one-way match (pattern matching): a substitution s with
// domain among the variables of pattern such that s.ApplyTerm(pattern) equals
// subject. Variables occurring only in subject are treated as constants and are
// never bound. It returns ErrUnify when no such match exists.
func MatchTerm(pattern, subject Term) (Substitution, error) {
	s := NewSubstitution()
	s, ok := match(pattern, subject, s)
	if !ok {
		return Substitution{}, ErrUnify
	}
	return s, nil
}

func match(p, t Term, s Substitution) (Substitution, bool) {
	if p.IsVar() {
		if img, ok := s.Get(p.Name); ok {
			if img.Equal(t) {
				return s, true
			}
			return s, false
		}
		return s.Bind(p.Name, t), true
	}
	if p.Kind != t.Kind || p.Name != t.Name || len(p.Args) != len(t.Args) {
		return s, false
	}
	var ok bool
	for i := range p.Args {
		s, ok = match(p.Args[i], t.Args[i], s)
		if !ok {
			return s, false
		}
	}
	return s, true
}

// Matches reports whether pattern one-way matches subject.
func Matches(pattern, subject Term) bool {
	_, err := MatchTerm(pattern, subject)
	return err == nil
}

// UnifyLiterals unifies two first-order literals of the same predicate and
// arity (their signs are ignored, so a positive and a negative literal on the
// same predicate may unify — the caller decides how to use the result). It
// returns ErrUnify when the predicates, arities, or argument terms cannot be
// unified.
func UnifyLiterals(a, b FOLiteral) (Substitution, error) {
	if a.Pred != b.Pred || len(a.Args) != len(b.Args) {
		return Substitution{}, ErrUnify
	}
	return UnifyTermLists(a.Args, b.Args)
}
