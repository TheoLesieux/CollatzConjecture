"""
Shared constants, helpers, and data loading for all Collatz analysis scripts.
"""
import os
import mpmath

# ── mpmath precision (matches _generate_data.py) ───────────
mpmath.mp.dps = 50

# ── Core constants (mpmath — full precision) ────────────────
gamma = mpmath.log(3) / mpmath.log(2)   # log_2(3) ≈ 1.58496…
alpha = gamma - 1                        # ≈ 0.58496…
beta  = 2 - gamma                        # ≈ 0.41504…
lam0  = gamma**gamma / (3 * alpha**alpha)

# ── Paths ───────────────────────────────────────────────────
REPO_ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FILE      = os.path.join(REPO_ROOT, 'data', 'collatz_data.txt')
FIG_DIR        = os.path.join(REPO_ROOT, 'paper', 'figures')

# ── Elementary functions ────────────────────────────────────
def s(n):
    """Stopping time s(n) = ceil(gamma * n), with s(0) = 1."""
    if n == 0:
        return 1
    return int(mpmath.ceil(gamma * n))

def eps(n):
    """Sturmian phase: eps(n) = ceil(gamma*n) - gamma*n  (mpmath)."""
    return mpmath.mpf(s(n)) - gamma * n

def jbeta(j):
    """Fractional part {j * beta}  (float64 for threshold tables)."""
    return float(j * beta) % 1

# ── Sub-band classification (11 bands for rho) ─────────────
thresholds_11 = sorted([
    0, jbeta(5), jbeta(10), float(2*gamma - 3), jbeta(3),
    float(beta), float(alpha), jbeta(4), jbeta(2), jbeta(7), jbeta(12), 1.0
])
sub_names = [
    'Ia', 'Ib', 'Ic', 'Id', 'Ie', 'II',
    'IIIa', 'IIIb', 'IIIc', 'IIId', 'IIIe',
]
sub_colors = [
    '#08306b', '#2166ac', '#4393c3', '#6baed6', '#92c5de',
    '#33a02c',
    '#e31a1c', '#fb6a4a', '#f46d43', '#fdae61', '#fee090',
]

def subband(e):
    for j in range(len(thresholds_11) - 1):
        if thresholds_11[j] <= e < thresholds_11[j + 1]:
            return j
    return len(thresholds_11) - 2

def subband_name(e):
    return sub_names[subband(e)]

# ── 3-band classification ──────────────────────────────────
def band3(e):
    if e < float(beta):   return 0   # Band I
    elif e < float(alpha): return 1  # Band II
    else:                  return 2  # Band III

# ── Data loading ────────────────────────────────────────────
def load_data():
    """Load collatz_data.txt and return a dict of mpmath / int values.

    Columns: n  epsilon(n)  s(n)  a_n  log2(a_n)  R(n)  b_n  v_n  B_n

    Float columns are parsed as mpmath.mpf (50-digit precision).
    ns and sn are plain Python int lists.
    vn and Bn are lists of exact Python ints.
    """
    return _load_full(DATA_FILE)

def _load_full(path):
    """Parse collatz_data.txt; all float columns as mpmath.mpf, integers exact."""
    ns_l, sn_l = [], []
    eps_l, an_l, log2_l, Rn_l, bn_l = [], [], [], [], []
    vn_l, Bn_l = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            ns_l.append(int(parts[0]))
            eps_l.append(mpmath.mpf(parts[1].strip()))
            sn_l.append(int(parts[2]))
            an_l.append(mpmath.mpf(parts[3].strip()))
            log2_l.append(mpmath.mpf(parts[4].strip()))
            Rn_l.append(mpmath.mpf(parts[5].strip()))
            bn_l.append(mpmath.mpf(parts[6].strip()))
            vn_l.append(int(parts[7]))
            Bn_l.append(int(parts[8]))
    return dict(
        ns      = ns_l,
        eps = eps_l,
        sn  = sn_l,
        an  = an_l,
        log2_an = log2_l,
        Rn  = Rn_l,
        bn  = bn_l,
        vn  = vn_l,
        Bn  = Bn_l,
        N       = len(ns_l),
    )
