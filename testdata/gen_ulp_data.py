#!/usr/bin/env python3

# SPDX-License-Identifier: MIT

"""
ULP stress test reference data generator.

Generates high-precision reference values for trigonometric and hyperbolic
functions using gmpy2, for ULP error measurement in Go tests.

Usage:
    python testdata/gen_ulp_data.py
    python testdata/gen_ulp_data.py --prec 128
    python testdata/gen_ulp_data.py --output-dir testdata
"""

import sys
import json
import random
import argparse
import math
import gmpy2

PRECISIONS = [64, 128, 256, 1024]
REF_MARGIN = 64  # extra bits for high-precision reference

# Go's math/big exponent limits
GO_EMAX = 2147483647
GO_EMIN = -2147483648


def configure_context(prec):
    """Set gmpy2 context to the given precision with Go-compatible exponent limits."""
    ctx = gmpy2.get_context()
    ctx.precision = prec
    ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
    ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())


def generate_test_points(prec, base_seed=42):
    """Generate deterministic test points for the given precision.

    Args:
        prec: Target precision in bits.
        base_seed: Base random seed; combined with prec for per-precision seed.

    Returns a list of dicts, each containing:
        x (str): hex float representation of the input
        sin_ref, cos_ref, sin_3x_ref, sinh_ref, cosh_ref, cosh_2x_ref (str): hex floats
    """
    ref_prec = prec + REF_MARGIN
    configure_context(ref_prec)

    seed = base_seed + prec  # different seed per precision for varied random points
    rng = random.Random(seed)

    points = []

    # Helper to add a point
    def add_point(x_val, name=""):
        """Compute all reference values for a given x and append the test point."""
        # Ensure x_val is mpfr for format('a') compatibility
        if not isinstance(x_val, gmpy2.mpfr):
            x_val = gmpy2.mpfr(str(x_val))
        # Compute with gmpy2 at reference precision
        sin_ref = gmpy2.sin(x_val)
        cos_ref = gmpy2.cos(x_val)
        sin_3x_ref = gmpy2.sin(3 * x_val)
        sinh_ref = gmpy2.sinh(x_val)
        cosh_ref = gmpy2.cosh(x_val)
        cosh_2x_ref = gmpy2.cosh(2 * x_val)

        points.append({
            "name": name,
            "x": format(x_val, 'a'),
            "sin_ref": format(sin_ref, 'a'),
            "cos_ref": format(cos_ref, 'a'),
            "sin_3x_ref": format(sin_3x_ref, 'a'),
            "sinh_ref": format(sinh_ref, 'a'),
            "cosh_ref": format(cosh_ref, 'a'),
            "cosh_2x_ref": format(cosh_2x_ref, 'a'),
        })

    # --- Trigonometric test points ---
    # Unit circle key points (as multiples of π)
    pi_val = gmpy2.const_pi()
    trig_key_points = [
        (0, "0"),
        (pi_val / 6, "π/6"),
        (pi_val / 4, "π/4"),
        (pi_val / 3, "π/3"),
        (pi_val / 2, "π/2"),
        (2 * pi_val / 3, "2π/3"),
        (3 * pi_val / 4, "3π/4"),
        (5 * pi_val / 6, "5π/6"),
        (pi_val, "π"),
        (3 * pi_val / 2, "3π/2"),
        (2 * pi_val, "2π"),
        (-pi_val / 2, "-π/2"),
        (-pi_val, "-π"),
        (-2 * pi_val, "-2π"),
    ]
    for x, name in trig_key_points:
        add_point(x, name)

    # Small-to-moderate values
    for x in ["1e-10", "1e-5", "0.5", "1.0", "2.0", "10.0"]:
        add_point(gmpy2.mpfr(x), x)

    # Large values (stress reducePi2)
    for x in ["100", "1e5", "1e10", "1e20"]:
        add_point(gmpy2.mpfr(x), x)

    # Negative small-to-moderate
    for x in ["-0.5", "-1.0", "-2.0", "-10.0", "-100", "-1e5", "-1e10", "-1e20"]:
        add_point(gmpy2.mpfr(x), x)

    # Random trigonometric test points in [-2π, 2π]
    for i in range(10):
        x = rng.uniform(-2 * math.pi, 2 * math.pi)
        add_point(gmpy2.mpfr(str(x)), f"rand_trig_{i}")

    # --- Hyperbolic test points ---
    # Systematic hyperbolic test points
    for x in ["0", "0.5", "1.0", "2.0", "5.0", "10.0"]:
        add_point(gmpy2.mpfr(x), x)
    for x in ["-0.5", "-1.0", "-2.0", "-5.0", "-10.0"]:
        add_point(gmpy2.mpfr(x), x)

    # Small hyperbolic test points
    for x in ["1e-10", "1e-5", "-1e-10", "-1e-5"]:
        add_point(gmpy2.mpfr(x), x)

    # Random hyperbolic test points in [-10, 10]
    for i in range(10):
        x = rng.uniform(-10, 10)
        add_point(gmpy2.mpfr(str(x)), f"rand_hyp_{i}")

    # Remove duplicate x values (same hex representation)
    seen = set()
    deduped = []
    for p in points:
        if p["x"] not in seen:
            seen.add(p["x"])
            deduped.append(p)

    return deduped


def generate_json(prec, base_seed=42):
    """Generate JSON-serializable data for a given precision."""
    test_points = generate_test_points(prec, base_seed=base_seed)
    return {
        "prec": prec,
        "ref_prec": prec + REF_MARGIN,
        "points": test_points,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate ULP stress test reference data for bigmath."
    )
    parser.add_argument(
        "--prec",
        type=int,
        default=None,
        help="Single precision to generate (default: all precisions)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Base random seed (default: 42)",
    )
    parser.add_argument(
        "--output-dir",
        default="testdata",
        help="Output directory (default: testdata)",
    )
    args = parser.parse_args()

    precs = [args.prec] if args.prec is not None else PRECISIONS

    for prec in precs:
        print(f"Generating reference data for prec={prec}... ", end="", flush=True)
        data = generate_json(prec, base_seed=args.seed)
        filename = f"{args.output_dir}/ulp_data_{prec}.json"
        with open(filename, "w") as f:
            json.dump(data, f, indent=2)
        print(f"wrote {len(data['points'])} points to {filename}")


if __name__ == "__main__":
    main()
