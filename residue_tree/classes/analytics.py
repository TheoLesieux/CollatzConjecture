class Analytics:
    """Pretty-prints the density table from a completed PairHandler run."""

    def __init__(self, handler):
        self.handler = handler

    def print_summary(self):
        """Print first-decrease density by step count n = log₂(M^(0)) = s(p)."""
        steps = sorted(self.handler.densities.keys())
        cumulative = 0.0
        print(f"\n{'n':<6} {'a(n)':<22} {'Cumulative':<22}")
        print("-" * 50)
        for n in steps:
            cumulative += self.handler.densities[n]
            print(f"{n:<6} {self.handler.densities[n]:<22.15f} {cumulative:<22.15f}")
        print(f"\nTotal density (excl. 1/2): {self.handler.total_density:.15f}")
        print(f"Total density (incl. 1/2): {self.handler.total_density + 0.5:.15f}")
        print(f"Queue size: {len(self.handler.queue)}")
