"""Collatz pair exploration — entry point.

Explores the residue-class tree via BFS under the reduced Collatz map
T(n) = n/2 if even, (3n+1)/2 if odd, and prints the cumulative
first-decrease density by step count.

Usage:
    cd residue_tree
    python main.py [--steps 100000]
"""
import argparse
import sys
sys.path.insert(0, '.')

from classes.pair_handler import PairHandler
from classes.analytics import Analytics


def main():
    parser = argparse.ArgumentParser(description="Collatz pair exploration")
    parser.add_argument("--steps", type=int, default=100_000,
                        help="Maximum BFS iterations (default: 100000)")
    args = parser.parse_args()

    handler = PairHandler()
    density = handler.run(max_steps=args.steps)

    analytics = Analytics(handler)
    analytics.print_summary()

    print(f"\nFinal cumulative density: {density:.10f}")
    print(f"Expected if conjecture true: 1.0")


if __name__ == "__main__":
    main()
