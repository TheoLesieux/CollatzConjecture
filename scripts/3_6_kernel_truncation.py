"""
High-precision kernel-depth truncation test using mpmath.
Verifies that u^{[D]} > 1/2 for all D, with 50-digit precision.
"""
import mpmath

from _collatz_common import s, eps, load_data

# ── Load b_n at full 50-digit precision from data file ───────
_D = load_data()
bn = [mpmath.mpf(0)] + _D['bn']   # 1-indexed: bn[n] for n=1..N

def c_plus(d):
    """c_+(d) = C(s(d), d) / 2^s(d)"""
    sd = s(d)
    return mpmath.binomial(sd, d) / mpmath.power(2, sd)

def c_minus(d):
    """c_-(d) = C(s(d)-1, d) / 2^(s(d)-1)"""
    sd = s(d)
    return mpmath.binomial(sd - 1, d) / mpmath.power(2, sd - 1)

N = 2000  # total terms

print(f"Computing with {mpmath.mp.dps}-digit precision, N={N}")
print()

# Precompute backbone (already loaded)
b_arr = bn[:N + 1]

# Precompute c_+, c_-
cp_arr = [mpmath.mpf(0)] + [c_plus(d) for d in range(1, N + 1)]
cm_arr = [mpmath.mpf(0)] + [c_minus(d) for d in range(1, N + 1)]

# Precompute epsilon
eps_arr = [mpmath.mpf(0)] + [eps(n) for n in range(1, N + 1)]

def kernel(n, d):
    """c(n,d) using precomputed arrays"""
    if eps_arr[n - d] + eps_arr[d] < 1:
        return cp_arr[d]
    else:
        return cm_arr[d]

# Compute exact a_n (full kernel) for reference
print("Computing exact a_n (full kernel)...")
a_exact = [mpmath.mpf(0)] * (N + 1)
for n in range(1, N + 1):
    val = b_arr[n]
    for d in range(1, n):
        val -= kernel(n, d) * a_exact[n - d]
    a_exact[n] = val
u_exact = sum(a_exact[1:])
print(f"  u_exact = {mpmath.nstr(u_exact, 30)}")
print(f"  u_exact - 1/2 = {mpmath.nstr(u_exact - mpmath.mpf('0.5'), 15)}")
print()

# Test depths
D_values = [1, 2, 3, 5, 7, 10, 15, 20, 30, 41, 53, 100, 200, 306, 500, 665, 1000, 2000]

print("=== EXACT KERNEL TRUNCATION (mpmath) ===")
print(f"{'D':>6} | {'u^{[D]}':>35} | {'u^{[D]} - 1/2':>20} | above 1/2?")
print("-" * 95)

prev_u = None
all_above = True
all_decreasing = True

for D in D_values:
    if D > N:
        break
    # Solve truncated renewal: a_n^{[D]} = b_n - sum_{d=1}^{min(D,n-1)} c(n,d) * a_{n-d}^{[D]}
    a_trunc = [mpmath.mpf(0)] * (N + 1)
    for n in range(1, N + 1):
        val = b_arr[n]
        for d in range(1, min(D, n - 1) + 1):
            val -= kernel(n, d) * a_trunc[n - d]
        a_trunc[n] = val
    u_trunc = sum(a_trunc[1:])
    diff = u_trunc - mpmath.mpf('0.5')
    above = diff > 0
    
    mono_str = ""
    if prev_u is not None:
        if u_trunc < prev_u:
            mono_str = "  \\ (decreasing)"
        elif u_trunc > prev_u:
            mono_str = "  / (INCREASING!)"
            all_decreasing = False
        else:
            mono_str = "  = (same)"
    
    if not above:
        all_above = False
    
    print(f"D={D:>5} | {mpmath.nstr(u_trunc, 35):>35} | {mpmath.nstr(diff, 15):>20} | {'YES' if above else 'NO!'}{mono_str}")
    prev_u = u_trunc

print()
print(f"All u^{{[D]}} > 1/2? {'YES' if all_above else 'NO'}")
print(f"Monotonically decreasing? {'YES' if all_decreasing else 'NO'}")
print()

# Also check: B, P+, P-, sum K
B_val = sum(b_arr[1:])
Pp = sum(cp_arr[1:])
Pm = sum(cm_arr[1:])
print(f"B  = {mpmath.nstr(B_val, 25)}")
print(f"P+ = {mpmath.nstr(Pp, 25)}")
print(f"P- = {mpmath.nstr(Pm, 25)}")
print(f"B/(1+P+) = {mpmath.nstr(B_val / (1 + Pp), 25)}  (lower sandwich)")
print(f"B/(1+P-) = {mpmath.nstr(B_val / (1 + Pm), 25)}  (upper sandwich)")
print(f"B - P+/2 + P-/2 check: does B = P+ - P-/2?  B={mpmath.nstr(B_val, 15)}, P+ - P-/2 = {mpmath.nstr(Pp - Pm/2, 15)}")
