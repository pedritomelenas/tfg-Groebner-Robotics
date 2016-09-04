from sympy import *

def s_polynomial(f, g):
	return expand(lcm(LM(f),LM(g))*(1/LT(f)*f-1 /LT(g)*g))

def buchberger(F):
	G = list(F)
	pairs = []
	for i, f1 in enumerate(F):
		for f2 in F[i + 1:]:
			pairs.append((f1, f2))
	while pairs:
		f1, f2 = pairs.pop(0)
		s = s_polynomial(f1, f2)
		h = reduced(s, G)[1]

		if h != 0:
			for g in G:
				pairs.append((g, h))
			G.append(h)
	return G
