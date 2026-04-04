"""
Base data generator — single source of truth.

Computes for n = 1..N_MAX all fundamental quantities using mpmath for
maximum floating-point precision (DPS decimal digits):

  v_n   (exact integer) — prime-combination count (Corollary 2.8 recurrence)
  B_n   (exact integer) — backbone integer  B_n = C(s(n)-1, n-1)
  a_n   = v_n / 2^{s(n)}          (mpmath)
  b_n   = B_n / 2^{s(n)}          (mpmath)
  R(n)  = n · v_n / B_n           (mpmath)
  epsilon(n) = ceil(γn) – γn      (mpmath, high-precision γ)
  log2(a_n)                        (mpmath)

Output: data/collatz_data.txt
Tab-separated columns:
  n  epsilon(n)  s(n)  a_n  log2(a_n)  R(n)  b_n  v_n  B_n

v_n and B_n are exact decimal integers; all float columns are written
with DPS significant digits.
"""
import sys
import os
from math import comb

import mpmath

# ── Locate repo root and adjust sys.path ────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from _collatz_common import gamma, s, eps

# ── v_n recurrence (Corollary 2.8) ──────────────────────────
def compute_v(n_max):
    """v_n = C(s(n)-1,n-1) - sum_{i<n} C(s(n)-s(i),n-i)*v_i"""
    v = [0] * (n_max + 1)
    v[0] = 1
    for n in range(1, n_max + 1):
        sn = s(n)
        v[n] = comb(sn - 1, n - 1) - sum(
            comb(sn - s(i), n - i) * v[i] for i in range(1, n)
        )
        if n % 500 == 0:
            print(f"  recurrence n={n}", file=sys.stderr)
    return v

# ── Precision (inherited from _collatz_common via mpmath.mp.dps) ─────
DPS   = mpmath.mp.dps
_log2 = mpmath.log(2)

_REPO_ROOT  = os.path.dirname(_SCRIPT_DIR)
N_MAX       = 5000
OUTPUT_PATH = os.path.join(_REPO_ROOT, "data", "collatz_data.txt")

def _fmt(x):
    return mpmath.nstr(x, DPS, strip_zeros=False)


def main():
    print(f"Computing v_n for n = 1..{N_MAX} (this may take several minutes)…",
          file=sys.stderr)
    v = compute_v(N_MAX)
    print("Recurrence done.  Writing output file…", file=sys.stderr)

    with open(OUTPUT_PATH, "w") as f:
        f.write(f"# Collatz data: n=1..{N_MAX},  gamma={mpmath.nstr(gamma, DPS)}\n")
        f.write(
            "# Columns: n  epsilon(n)  s(n)  a_n  log2(a_n)  R(n)  b_n  v_n  B_n\n"
        )
        f.write(
            "# Definitions: a_n = v_n/2^s(n),  b_n = B_n/2^s(n),  R(n) = n*v_n/B_n\n"
        )
        f.write(
            f"# Float columns computed with mpmath at {DPS} decimal digits of precision\n"
        )

        for n in range(1, N_MAX + 1):
            sn    = s(n)
            eps_n = eps(n)
            vn    = v[n]                      # exact Python integer
            Bn    = comb(sn - 1, n - 1)       # exact Python integer

            pow2    = mpmath.power(2, sn)
            a_n     = mpmath.mpf(vn) / pow2
            b_n     = mpmath.mpf(Bn) / pow2
            log2_an = mpmath.log(a_n) / _log2 if vn > 0 else mpmath.mpf('-inf')
            Rn      = mpmath.mpf(n * vn) / mpmath.mpf(Bn) if Bn != 0 else mpmath.mpf(0)

            f.write(
                f"{n}\t{_fmt(eps_n)}\t{sn}\t"
                f"{_fmt(a_n)}\t{_fmt(log2_an)}\t"
                f"{_fmt(Rn)}\t{_fmt(b_n)}\t"
                f"{vn}\t{Bn}\n"
            )

            if n % 500 == 0:
                print(f"  n={n} written", file=sys.stderr)

    print(f"Wrote {OUTPUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
