class Pair:
    """A node in the residue-class exploration tree.

    Represents the residue class {x = M·k + r : k ≥ 0} under the reduced
    Collatz map T(n) = n/2 if even, (3n+1)/2 if odd.

    Attributes
    ----------
    seed : tuple[int, int]
        (M^(0), r^(0)) — the original residue class, fixed at creation.
        M^(0) is always a power of 2; the natural density of the class is 1/M^(0).
    transformed : tuple[int, int]
        (M^(m), r^(m)) — the residue class after m applications of T.
        Updated in place as the pair is processed.
    """

    def __init__(self, seed_M, seed_r,
                 trans_M=None, trans_r=None):
        self.seed = (seed_M, seed_r)
        if trans_M is None:
            self.transformed = (seed_M, seed_r)
        else:
            self.transformed = (trans_M, trans_r)

    def __repr__(self):
        return f"Pair(seed={self.seed}, transformed={self.transformed})"
