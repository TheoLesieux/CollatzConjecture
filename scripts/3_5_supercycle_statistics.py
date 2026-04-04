"""
Supercycle statistics (Section 3.5).

Verifies supercycle decomposition and all numerical claims:
  - 849 total supercycles (376 A, 376 B, 97 C)
  - Asymptotic multipliers M_A^inf ~ 5.025, M_B^inf ~ 3.608, M_C^inf ~ 3.625
  - P_max < 0.76, M_min ~ 3.59, M_max ~ 5.03
  - Bracket bounds at N0 = 5000
"""
import mpmath

from _collatz_common import load_data, subband_name

# Load data
D = load_data()
ns  = D['ns']
eps = D['eps']
an  = D['an']
N   = D['N']

# Assign sub-bands
subs = [subband_name(e) for e in eps]

# Find Ie indices (supercycle hubs)
ie_idx = [i for i in range(N) if subs[i] == 'Ie']
print(f"Number of Ie indices: {len(ie_idx)}")

# Classify each Ie-to-Ie supercycle
sc_types = []
sc_data = []  # (type, n_start, P, M)

for k in range(len(ie_idx) - 1):
    i0, i1 = ie_idx[k], ie_idx[k + 1]
    L = i1 - i0
    n0 = ns[i0]
    P = an[i1] / an[i0]                               # return product
    M = mpmath.fsum(an[i0:i1]) / an[i0]                    # multiplier
    route = set(subs[i] for i in range(i0, i1 + 1))

    if L == 7 and 'Id' in route:
        sc_type = 'A'
    elif L == 5 and 'Ib' in route:
        sc_type = 'B'
    elif L == 5 and 'Ic' in route:
        sc_type = 'C'
    else:
        sc_type = '?'

    sc_types.append(sc_type)
    sc_data.append((sc_type, n0, P, M))

count_A = sum(1 for t in sc_types if t == 'A')
count_B = sum(1 for t in sc_types if t == 'B')
count_C = sum(1 for t in sc_types if t == 'C')
count_X = sum(1 for t in sc_types if t == '?')
total   = len(sc_types)

print(f"\n=== SUPERCYCLE COUNTS ===")
print(f"Total supercycles: {total}")
print(f"  A (L=7): {count_A}")
print(f"  B (L=5): {count_B}")
print(f"  C (L=5): {count_C}")
if count_X > 0:
    print(f"  Unknown: {count_X}")

# Asymptotic multipliers (n >= 2000)
def asymptotic_stats(label):
    vals_P = [P for (t, n, P, M) in sc_data if t == label and n >= 2000]
    vals_M = [M for (t, n, P, M) in sc_data if t == label and n >= 2000]
    P_mean = mpmath.fsum(vals_P) / len(vals_P)
    M_mean = mpmath.fsum(vals_M) / len(vals_M)
    return P_mean, M_mean, vals_P, vals_M

print(f"\n=== ASYMPTOTIC PARAMETERS (n >= 2000) ===")
for label in ['A', 'B', 'C']:
    P_mean, M_mean, Ps, Ms = asymptotic_stats(label)
    print(f"  {label}: P^inf = {mpmath.nstr(P_mean, 10)}, M^inf = {mpmath.nstr(M_mean, 10)} "
          f"(N={len(Ps)})")

# Observed extrema (n >= 200)
all_P = [P for (t, n, P, M) in sc_data if n >= 200]
all_M = [M for (t, n, P, M) in sc_data if n >= 200]
P_max = max(all_P)
P_min = min(all_P)
M_max = max(all_M)
M_min = min(all_M)

print(f"\n=== EXTREMA (n >= 200) ===")
print(f"  P_max = {mpmath.nstr(P_max, 10)}  (article claims P_max < 0.76)")
print(f"  P_min = {mpmath.nstr(P_min, 10)}")
print(f"  M_max = {mpmath.nstr(M_max, 10)}  (article claims ~5.03)")
print(f"  M_min = {mpmath.nstr(M_min, 10)}  (article claims ~3.59)")
print(f"  All P < 1? {all(p < 1 for p in all_P)}")

# Bracket bounds
u_N0 = mpmath.fsum(an)                   # u_{5000}
last_ie = ie_idx[-1]
a_nJ = an[last_ie]
n_J  = ns[last_ie]

# Full geometric tail from n_J: a_{n_J} * M / (1-P)
geom_lower = a_nJ * M_min / (1 - P_max)
geom_upper = a_nJ * M_max / (1 - P_min)

# Already counted from n_J to N0=5000
already_counted = mpmath.fsum(an[last_ie:])

# Tail beyond N0 = 5000
tail_lower = geom_lower - already_counted
tail_upper = geom_upper - already_counted

print(f"\n=== BRACKET AT N0=5000 ===")
print(f"  u_N0 = {mpmath.nstr(u_N0, 30)}")
print(f"  |u_N0 - 1/2| = {mpmath.nstr(abs(u_N0 - mpmath.mpf('0.5')), 10)}")
print(f"  a_{{n_J}} = {mpmath.nstr(a_nJ, 10)}  (n_J = {n_J})")
print(f"  Geometric series [n_J..inf]: [{mpmath.nstr(geom_lower, 6)}, {mpmath.nstr(geom_upper, 6)}]")
print(f"  Already counted [n_J..5000]: {mpmath.nstr(already_counted, 6)}")
print(f"  Tail beyond 5000: [{mpmath.nstr(tail_lower, 6)}, {mpmath.nstr(tail_upper, 6)}]")
print(f"  Article claims: 3.0e-125 <= tail <= 4.6e-124")

# Supercycle word structure
print(f"\n=== SUPERCYCLE WORD ===")
word = ''.join(sc_types[:30])
print(f"  First 30: {word}")
# Check AB alternation
ab_word = [t for t in sc_types if t in ('A', 'B')]
alternating = all(ab_word[i] != ab_word[i + 1] for i in range(len(ab_word) - 1))
print(f"  A-B perfectly alternating? {alternating}")

# C insertion gaps
c_positions = [i for i, t in enumerate(sc_types) if t == 'C']
c_gaps = [c_positions[i + 1] - c_positions[i] for i in range(len(c_positions) - 1)]
unique_gaps = sorted(set(c_gaps))
print(f"  Gaps between C insertions: {unique_gaps} (expect [7, 9])")
print(f"  First 15 gaps: {c_gaps[:15]}")
