#!/usr/bin/env python3

# SPDX-License-Identifier: MIT

import sys
import json
import argparse
import gmpy2
import builtins


def format_or_nan(val):
    """Format an mpfr value as hex float string.
    Returns (hex_string, is_nan) tuple.
    is_nan is True if the value is NaN."""
    if gmpy2.is_nan(val):
        return None, True
    return format(val, 'a'), False


def round_input(val, prec):
    """Round an mpfr value to the given precision, matching Go's SetPrec(prec).Parse(s,0)
    behavior: the decimal string is parsed at the target precision, so the result is
    already rounded to prec bits. Re-rounding via mpfr creation at prec ensures the
    value is at exactly prec bits regardless of any context precision changes."""
    saved = gmpy2.get_context().precision
    ctx = gmpy2.get_context()
    ctx.precision = prec
    r = gmpy2.mpfr(val)
    ctx.precision = saved
    return r


def generate_json_tests(input_file, output_file, precision):
    ctx = gmpy2.get_context()
    ctx.precision = precision
    # Attempt to match Go's math/big exponent limits as closely as MPFR allows
    # Go: -2^31 to 2^31-1. MPFR is often limited to 2^30-1 or 2^62-1.
    GO_EMAX = 2147483647
    GO_EMIN = -2147483648
    ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
    ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())

    cases = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            func_name = parts[0]
            raw_args = parts[1:]

            # Detect !panic annotation
            is_panic = "!panic" in raw_args
            str_args = [a for a in raw_args if a != "!panic"]

            # For constant functions, ignore passed args (they are placeholders)
            if func_name.startswith("const_"):
                str_args = []

            # Parse arguments at target precision to match Go's
            # SetPrec(prec).Parse(s,0) — the input decimal is rounded to
            # prec bits at parse time.
            mpfr_args = [round_input(gmpy2.mpfr(arg), precision) for arg in str_args]

            # Lookup function: gmpy2 -> globals -> builtins
            if func_name.startswith("const_"):
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    raise AttributeError(f"constant function '{func_name}' not found in gmpy2")
                result = func()
            elif func_name in ("floor", "ceil"):
                if gmpy2.is_infinite(mpfr_args[0]):
                    result = mpfr_args[0]  # ±Inf stays ±Inf
                else:
                    func = gmpy2.floor if func_name == "floor" else gmpy2.ceil
                    result = gmpy2.mpfr(func(mpfr_args[0]))
            elif func_name == "acot":
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = gmpy2.atan2(1, mpfr_args[0])
                ctx.precision = saved_prec
            elif func_name == "asec":
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = gmpy2.acos(1 / mpfr_args[0])
                ctx.precision = saved_prec
            elif func_name == "acsc":
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = gmpy2.asin(1 / mpfr_args[0])
                ctx.precision = saved_prec
            elif func_name == "acoth":
                # Compute at 2x precision to avoid double-rounding from 1/x step.
                # Use direct formula ½·ln((x+1)/(x-1)) matching Go's implementation.
                # Handle Inf explicitly since (x+1)/(x-1) = ∓1 for ±Inf leading to NaN.
                if gmpy2.is_infinite(mpfr_args[0]):
                    result = gmpy2.mpfr(0)
                    if mpfr_args[0] < 0:
                        result = -result
                else:
                    saved_prec = gmpy2.get_context().precision
                    ctx = gmpy2.get_context()
                    ctx.precision = precision * 2
                    t = mpfr_args[0]
                    result = gmpy2.log((t + 1) / (t - 1)) / 2
                    ctx.precision = saved_prec
            elif func_name == "asech":
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = gmpy2.acosh(1 / mpfr_args[0])
                ctx.precision = saved_prec
            elif func_name == "acsch":
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = gmpy2.asinh(1 / mpfr_args[0])
                ctx.precision = saved_prec
            elif func_name == "lgamma":
                lg, sign = gmpy2.lgamma(mpfr_args[0])
                # Wrap sign int as mpfr for the tuple handler
                result = (lg, gmpy2.mpfr(sign))
            else:
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    func = globals().get(func_name)
                if func is None:
                    func = getattr(builtins, func_name, None)
                if func is None:
                    raise AttributeError(f"function '{func_name}' not found")
                result = func(*mpfr_args)

            # Handle tuple results (sin_cos, sinh_cosh)
            if isinstance(result, tuple):
                entry = {"fn": func_name, "args": str_args}
                hex_val0, nan0 = format_or_nan(result[0])
                hex_val1, nan1 = format_or_nan(result[1])
                if nan0 or nan1 or is_panic:
                    entry["panics"] = True
                    if nan0 or nan1:
                        entry["reason"] = "gmpy2 returned NaN"
                    if is_panic and not (nan0 or nan1):
                        print(
                            f"Warning: {func_name} {str_args} tagged !panic "
                            f"but gmpy2 did not produce NaN",
                            file=sys.stderr,
                        )
                else:
                    entry["res"] = hex_val0
                    entry["res2"] = hex_val1
                cases.append(entry)
                continue

            # Single result
            hex_val, is_nan = format_or_nan(result)

            entry = {"fn": func_name, "args": str_args}
            if is_nan or is_panic:
                entry["panics"] = True
                if is_nan:
                    entry["reason"] = "gmpy2 returned NaN"
                if is_panic and not is_nan:
                    print(
                        f"Warning: {func_name} {str_args} tagged !panic "
                        f"but gmpy2 did not produce NaN",
                        file=sys.stderr,
                    )
                if is_nan and not is_panic:
                    print(
                        f"Warning: {func_name} {str_args} NOT tagged !panic "
                        f"but gmpy2 did produce NaN",
                        file=sys.stderr,
                    )
            else:
                entry["res"] = hex_val

            cases.append(entry)

    output = {"prec": precision, "cases": cases}

    out = open(output_file, 'w') if output_file else sys.stdout
    try:
        json.dump(output, out, indent=2)
        out.write('\n')
    finally:
        if output_file:
            out.close()


def generate_cplx_tests(input_file, output_file, precision):
    """Generate complex golden test data from cplx_data.txt format.

    Each input line:
        <func_name> <re1> <im1> [<re2> <im2> ...] [!panic]

    Parses real/imag pairs into gmpy2.mpc values, dispatches via gmpy2,
    and writes res_re/res_im for non-panic cases.
    """
    ctx = gmpy2.get_context()
    ctx.precision = precision
    GO_EMAX = 2147483647
    GO_EMIN = -2147483648
    ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
    ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())

    cases = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            func_name = parts[0]
            raw_pairs = parts[1:]

            # Detect !panic annotation
            is_panic = "!panic" in raw_pairs
            tokens = [t for t in raw_pairs if t != "!panic"]

            # Pair-wise parse into mpc values, rounding each component to
            # target precision to match Go's SetPrec(prec).Parse(s,0).
            mpc_args = []
            for i in range(0, len(tokens), 2):
                re = round_input(gmpy2.mpfr(tokens[i]), precision)
                im = round_input(gmpy2.mpfr(tokens[i + 1]), precision)
                mpc_args.append(gmpy2.mpc(re, im))

            # Build args pairs as strings for JSON output
            str_args = [[tokens[i], tokens[i + 1]] for i in range(0, len(tokens), 2)]

            # Lookup function and call with mpc_args
            # Some Go function names differ from gmpy2 names, handled here
            if func_name == "cbrt":
                # gmpy2.cbrt does not support mpc. Compute via polar
                # decomposition at 2x precision to match Go's sequential
                # Float.Cbrt + Sincos computation path.
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                z = mpc_args[0]
                rho = gmpy2.cbrt(abs(z))
                theta = gmpy2.phase(z) / 3
                result = gmpy2.mpc(rho * gmpy2.cos(theta), rho * gmpy2.sin(theta))
                ctx.precision = saved_prec
            elif func_name == "cot":
                # gmpy2 does not provide cot for mpc, compute as cos/sin
                # Compute at 2x precision then round down to match Go's
                # sequential Cos/Sin → Quo at target precision.
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                result = gmpy2.cos(mpc_args[0]) / gmpy2.sin(mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "sec":
                # gmpy2 does not provide sec for mpc, compute as 1/cos
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = 1 / gmpy2.cos(mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "csc":
                # gmpy2 does not provide csc for mpc, compute as 1/sin
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = 1 / gmpy2.sin(mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "coth":
                # gmpy2 does not provide coth for mpc, compute as 1/tanh
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = 1 / gmpy2.tanh(mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "sech":
                # gmpy2 does not provide sech for mpc, compute as 1/cosh
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = 1 / gmpy2.cosh(mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "csch":
                # gmpy2 does not provide csch for mpc, compute as 1/sinh
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                result = 1 / gmpy2.sinh(mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "acot":
                # gmpy2 does not provide acot for mpc, compute as π/2 - atan(z)
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                result = gmpy2.const_pi() / 2 - gmpy2.atan(mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "asec":
                # gmpy2 does not provide asec for mpc, compute as acos(1/z)
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                result = gmpy2.acos(1 / mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "acsc":
                # gmpy2 does not provide acsc for mpc, compute as asin(1/z)
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                result = gmpy2.asin(1 / mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "acoth":
                # gmpy2 does not provide acoth for mpc, compute as atanh(1/z)
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                result = gmpy2.atanh(1 / mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "asech":
                # gmpy2 does not provide asech for mpc, compute as acosh(1/z)
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                result = gmpy2.acosh(1 / mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "acsch":
                # gmpy2 does not provide acsch for mpc, compute as asinh(1/z)
                saved_prec = gmpy2.get_context().precision
                ctx = gmpy2.get_context()
                ctx.precision = precision * 2
                ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
                ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())
                result = gmpy2.asinh(1 / mpc_args[0])
                ctx.precision = saved_prec
            elif func_name == "quo":
                result = gmpy2.div(*mpc_args)
            elif func_name == "neg":
                result = -mpc_args[0]
            elif func_name == "conj":
                result = mpc_args[0].conjugate()
            elif func_name == "pow":
                result = builtins.pow(*mpc_args)
            elif func_name.startswith("const_"):
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    raise AttributeError(f"constant function '{func_name}' not found in gmpy2")
                result = func()
            else:
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    raise AttributeError(f"complex function '{func_name}' not found in gmpy2")
                result = func(*mpc_args)

            # Extract real and imaginary parts
            hex_re, nan_re = format_or_nan(result.real)
            hex_im, nan_im = format_or_nan(result.imag)

            entry = {"fn": func_name, "args": str_args}

            if nan_re or nan_im or is_panic:
                entry["panics"] = True
                if nan_re or nan_im:
                    entry["reason"] = "gmpy2 returned NaN"
                if is_panic and not (nan_re or nan_im):
                    print(
                        f"Warning: {func_name} {tokens} tagged !panic "
                        f"but gmpy2 did not produce NaN",
                        file=sys.stderr,
                    )
                if (nan_re or nan_im) and not is_panic:
                    print(
                        f"Warning: complex {func_name} {tokens} NOT tagged !panic "
                        f"but reference produced NaN",
                        file=sys.stderr,
                    )
            else:
                entry["res_re"] = hex_re
                entry["res_im"] = hex_im

            cases.append(entry)

    output = {"prec": precision, "mode": "cplx", "cases": cases}

    out = open(output_file, 'w') if output_file else sys.stdout
    try:
        json.dump(output, out, indent=2)
        out.write('\n')
    finally:
        if output_file:
            out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate test data in JSON format.")
    parser.add_argument("input", help="Input text file")
    parser.add_argument("-o", "--output", help="Output JSON file")
    parser.add_argument(
        "-p", "--precision", type=int, default=256,
        help="Precision in bits (default: 256)",
    )
    parser.add_argument(
        "-m", "--mode", choices=["float", "cplx"], default="float",
        help="Test data mode (default: float)",
    )

    args = parser.parse_args()
    if args.mode == "cplx":
        generate_cplx_tests(args.input, args.output, args.precision)
    else:
        generate_json_tests(args.input, args.output, args.precision)
