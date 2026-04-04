# Residue-Class Tree Explorer

## Introduction

This module provides an independent constructive exploration of the Collatz conjecture at the level of residue classes.  Rather than counting first-descent classes via the analytic recurrence (Corollary 2.8 in the paper), it discovers them directly by applying the reduced Collatz map to residue-class pairs $(M, r)$:

$$T(n) = \begin{cases} n/2 & \text{if } n \text{ is even} \\ (3n+1)/2 & \text{if } n \text{ is odd} \end{cases}$$

Every positive integer belongs to a residue class $(M, r)$ with modulus $M = 2^n$ and remainder $0 \leq r < M$, meaning $x \equiv r \pmod{M}$. The algorithm seeds with $(M, r) = (2, 0)$ (even numbers) and $(2, 1)$ (odd numbers), applies the reduced Collatz map $T$, and detects when each class achieves its **first decrease**. The cumulative first-decrease density should converge to $1$ if the conjecture is true.

## The concept

A residue class is represented by a pair $(M, r)$. A positive integer $x$ belongs to it if $x = M k + r$ for some integer $k \geq 0$.

Two pairs are tracked per node in the exploration tree:
- **Seed pair**: the original residue class $(M_0, r_0)$, never modified
- **Transformed pair**: the class $(M^{(m)}, r^{(m)})$ after $m$ applications of $T$

A **first decrease** is detected when the transformed pair represents values strictly smaller than the seed pair — that is, when the modulus condition $3^p < 2^n$ (where $p$ is the rank, i.e. the number of $\uparrow$ steps) guarantees $x^{(n)} < x$ for all $k \geq 1$.

### The four cases

Each application of $T$ to the current pair $(M, r)$ falls into one of four cases depending on the parities of $M$ and $r$:

| Case | Condition | Step | Transformation |
|------|-----------|------|----------------|
| Even-even | $M$ even, $r$ even | $\downarrow$ | $(M, r) \to (M/2,\; r/2)$ |
| Even-odd | $M$ even, $r$ odd | $\uparrow$ | $(M, r) \to (3M/2,\; (3r+1)/2)$ |
| Odd-odd | $M$ odd, $r$ odd | — | Split (see below) |
| Odd-even | $M$ odd, $r$ even | — | Split (see below) |

**Even modulus (deterministic).** When $M$ is even, the parity of $x = M k + r$ is determined by $r$ alone (since $M k$ is always even). The step type ($\downarrow$ or $\uparrow$) is therefore known and the transformation can be applied directly.

**Odd modulus (split).** When $M$ is odd, the parity of $x = M k + r$ depends on whether $k$ is even or odd. Since both cases occur, the pair must be split into two children with doubled modulus $2M$, each of which has even modulus and can then be transformed:

**Odd-odd** ($M$ odd, $r$ odd):
- Even $k$: $x = 2M k' + r$ has $r$ odd → apply $\uparrow$: $(M, r) \to (3M,\; (3r+1)/2)$
- Odd $k$: $x = 2M k' + (M + r)$ has $M + r$ even → apply $\downarrow$: $(M, r) \to (M,\; (M+r)/2)$

**Odd-even** ($M$ odd, $r$ even):
- Even $k$: $x = 2M k' + r$ has $r$ even → apply $\downarrow$: $(M, r) \to (M,\; r/2)$
- Odd $k$: $x = 2M k' + (M + r)$ has $M + r$ odd → apply $\uparrow$: $(M, r) \to (3M,\; (3(M+r)+1)/2)$

In both split cases, the seed modulus doubles to $2M$ (halving the density), and the seed remainder is $r$ for the even-$k$ child and $M + r$ for the odd-$k$ child. This ensures every integer is covered exactly once across all leaves of the exploration tree.

## Usage

```bash
cd residue_tree
python main.py --steps 100000
```

## Structure

```
main.py                    — Entry point
classes/pair.py            — Pair data class (seed + transformed pairs)
classes/pair_handler.py    — Queue-based BFS through the residue-class tree
classes/analytics.py       — Progress display
```

## Relation to the paper

The paper derives the prime combination counts $v_p$ analytically via the recurrence in Corollary 2.8, using only the stopping time $s(p) = \lceil p\gamma \rceil$ and binomial coefficients. This module arrives at the same $v_p$ values through a completely independent route: direct constructive enumeration of residue-class pairs by BFS.

The two approaches agree exactly on every $v_p$ reached (Table A.1 in the paper), which serves as mutual validation. Crucially, the BFS also provides the **only available analysis of the $k = 0$ outcast condition**: the paper's Theorem 2.5 proves first decrease for $k \geq 1$, but the single representative $x = r < 2^{s(p)}$ in each residue class requires separate treatment (Remark 2.6). The BFS tracks all integers including $k = 0$, and confirms that no outcasts exist at any completed rank. The safety margins $r^{(0)}/\tau$ correlate tightly with the Sturmian phase $\varepsilon(p)$, reflecting the continued-fraction structure of $\gamma = \log_2 3$ (Figure A.1 in the paper).
