"""Collatz Conjecture — Residue-Class Tree Explorer

Explores the Collatz conjecture by tracking residue-class pairs (α, β)
through the reduced Collatz map T(n) = n/2 if even, (3n+1)/2 if odd.

Every positive integer belongs to a residue class (a, b) meaning x ≡ b (mod a).
The algorithm starts with (2,0) [even] and (2,1) [odd], applies the Collatz map,
and tracks when each class first decreases. The cumulative probability of decrease
should converge to 1 if the conjecture is true.

Usage:
    python main.py [--steps 100000]

Structure:
    main.py              — Entry point: runs the pair exploration
    classes/pair.py      — Pair data class (seed + transformed tuples)
    classes/pair_handler.py — Queue-based BFS through the residue-class tree
    classes/analytics.py — Progress display and density table
"""
