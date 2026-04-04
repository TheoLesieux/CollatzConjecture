from collections import deque
from classes.pair import Pair


class PairHandler:
    """Breadth-first exploration of the residue-class tree.

    Starting from the two root classes (2, 0) and (2, 1), every pair is
    either a *leaf* (first decrease detected) or an *internal node* that
    spawns one or two children via the reduced Collatz map T.

    The reduced Collatz map is T(x) = x/2 if x even (↓ step),
    T(x) = (3x+1)/2 if x odd (↑ step).

    Attributes
    ----------
    queue : deque[Pair]
        FIFO queue of pairs awaiting processing.
    densities : dict[int, float]
        Maps step count n = log₂(M^(0)) to the accumulated first-decrease
        density at that depth.  The step count n equals s(p) = ⌈γp⌉ in
        the paper's notation, where p is the rank (number of ↑ steps).
    total_density : float
        Running sum of all recorded densities, i.e. Σ a_p so far.
    """

    def __init__(self):
        self.queue = deque()
        self.densities = {}
        self.total_density = 0.0

    # ------------------------------------------------------------------
    # Decrease detection
    # ------------------------------------------------------------------

    @staticmethod
    def _has_decreased(pair):
        """Return True if the transformed class is strictly below the seed.

        A decrease means that for every k ≥ 1 in the seed class, the
        transformed value x^(m) = M^(m)k + r^(m) satisfies x^(m) < x = M^(0)k + r^(0).
        This holds when M^(m) < M^(0) and r^(m) ≤ r^(0), or M^(m) = M^(0) and r^(m) < r^(0).
        """
        M_s, r_s = pair.seed
        M_t, r_t = pair.transformed
        if M_t == M_s:
            return r_t < r_s
        if M_t < M_s:
            return r_t <= r_s
        return False

    # ------------------------------------------------------------------
    # Deterministic transformations (α even)
    # ------------------------------------------------------------------

    @staticmethod
    def _apply_down_step(pair):
        """↓ step: M even, r even → (M/2, r/2)."""
        M, r = pair.transformed
        pair.transformed = (M // 2, r // 2)
        return pair

    @staticmethod
    def _apply_up_step(pair):
        """↑ step: M even, r odd → (3M/2, (3r+1)/2)."""
        M, r = pair.transformed
        pair.transformed = (3 * M // 2, (3 * r + 1) // 2)
        return pair

    # ------------------------------------------------------------------
    # Split transformations (α odd)
    # ------------------------------------------------------------------

    @staticmethod
    def _split_odd_odd(pair):
        """Split when M odd, r odd: parity of x depends on k.

        The class {Mk + r} is refined to modulus 2M.  Each child then
        has even modulus, so one T-step can be applied immediately:

        Even-k child: {2Mk' + r}, r odd → ↑ step: (M, r) → (3M, (3r+1)/2)
        Odd-k child:  {2Mk' + (M+r)}, M+r even → ↓ step: (M, r) → (M, (M+r)/2)
        """
        M, r = pair.transformed
        a_s, b_s = pair.seed
        even_child = Pair(2 * a_s, b_s, 3 * M, (3 * r + 1) // 2)
        odd_child = Pair(2 * a_s, a_s + b_s, M, (M + r) // 2)
        return even_child, odd_child

    @staticmethod
    def _split_odd_even(pair):
        """Split when M odd, r even: parity of x depends on k.

        Even-k child: {2Mk' + r}, r even → ↓ step: (M, r) → (M, r/2)
        Odd-k child:  {2Mk' + (M+r)}, M+r odd → ↑ step: (M, r) → (3M, (3(M+r)+1)/2)
        """
        M, r = pair.transformed
        a_s, b_s = pair.seed
        even_child = Pair(2 * a_s, b_s, M, r // 2)
        odd_child = Pair(2 * a_s, a_s + b_s, 3 * M, (3 * (M + r) + 1) // 2)
        return even_child, odd_child

    # ------------------------------------------------------------------
    # BFS engine
    # ------------------------------------------------------------------

    def seed(self):
        """Initialise the tree with the two base residue classes mod 2."""
        self.queue.append(Pair(2, 0))  # even integers
        self.queue.append(Pair(2, 1))  # odd integers

    def step(self):
        """Process one pair from the queue.

        Returns False when the queue is empty (exploration complete).
        """
        if not self.queue:
            return False

        pair = self.queue.popleft()

        # --- Leaf: first decrease detected ---
        if self._has_decreased(pair):
            M_seed = pair.seed[0]
            density = 1 / M_seed          # natural density 1/M^(0)
            n = M_seed.bit_length() - 1   # n = log₂(M^(0)) = s(p)
            self.densities[n] = self.densities.get(n, 0) + density
            self.total_density += density
            return True

        # --- Internal node: apply T or split ---
        M, r = pair.transformed

        if M % 2 == 0:
            # M even → parity determined by r
            if r % 2 == 0:
                self.queue.append(self._apply_down_step(pair))
            else:
                self.queue.append(self._apply_up_step(pair))
        else:
            # M odd → must split on parity of k
            if r % 2 == 1:
                c1, c2 = self._split_odd_odd(pair)
            else:
                c1, c2 = self._split_odd_even(pair)
            self.queue.append(c1)
            self.queue.append(c2)

        return True

    def run(self, max_steps=100_000):
        """Run up to *max_steps* BFS iterations.

        Returns the cumulative first-decrease density reached.
        """
        self.seed()
        for i in range(max_steps):
            if not self.step():
                break
            if (i + 1) % 10_000 == 0:
                print(f"Step {i+1:>8,}: queue={len(self.queue):>7,}  "
                      f"P={self.total_density:.6f}")
        return self.total_density
