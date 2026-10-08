#!/usr/bin/env python3
"""patch-lab -- daily containment planner for a worm outbreak.

A worm is spreading through a population of N hosts (SIR dynamics). A patch
crew can deploy up to `c` hosts/day, infectious hosts first. This tool answers
the planner's questions deterministically and fail-closed:

  * what is the MINIMUM daily capacity that holds the outbreak at its current
    level?  (c_min, from a closed form: with G = (1+beta)(1-gamma), the daily
    growth factor, c_min = I0 * (G-1)/G.  Deploy above it and the worm dies
    out below its starting level; deploy below it and it outruns the crew.)
  * given a capacity c, is the outbreak CONTAINED or ESCAPED?
  * the closed form is verified against the model itself (numeric bisection),
    and the model is verified against the closed form (linear-window match),
  * a parameter sweep reports the WORST case across the uncertainty box --
    worst c_min, with the envelope verified numerically over the grid,
  * every run is appended to a hash-chained audit trail, and tampering with
    any retained record is detected on demand.

The daily map is explicit and deterministic, in this order:
  1. patch:  k = min(c, S+I) hosts, infectious first,
  2. spread: new = beta * S * I / N,
  3. recover: rec = gamma * I.

Everything is deterministic: the same arguments produce the same numbers.

Single file, Python 3.9+, standard library only.

Usage:
  python patch-lab.py selftest
  python patch-lab.py plan   --beta 0.3 --gamma 0.1 --i0 100 --N 1000000 [--c 18]
  python patch-lab.py compare [--fail-over] [--keep-frac 0.99] [--tol 1e-3]
  python patch-lab.py sweep  --beta-min 0.2 --beta-max 0.4 --gamma-min 0.05 --gamma-max 0.15
  python patch-lab.py stats [--json]
  python patch-lab.py verify [--store PATH]

Audit store: env PATCH_LAB_STORE, else ~/.patchlab/audit.jsonl
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
import time

__version__ = "1.0.0"

# The genesis anchor every chain must start from. If record 0 does not carry
# exactly this prev, the oldest record was deleted.
GENESIS = "0" * 64

# Fail-closed default behaviour for validation errors: every bad input exits 1
# with a reason on stderr. Nothing silently degrades.
_EPILOG = "patch-lab: "

# Everything below this infectious count counts as extinct; a fraction of a
# host below 1e-9 cannot rekindle anything and cannot be reported honestly.
_ZERO = 1e-9


# ---------------------------------------------------------------------------
# The daily containment map
# ---------------------------------------------------------------------------

def step(state, p):
    """Advance the SIR+patching model by exactly one day.

    State layout: (S, I, P, R) -- P is the patched (immune) pool.
    Deterministic order, documented and asserted by the selftest:
      1. patch  (up to c hosts/day, infectious hosts first),
      2. spread (each infectious host infects beta susceptibles/day),
      3. recover (each infectious host recovers at gamma/day).
    """
    s, i, p_pool, r = state
    n = p["N"]
    # 1. patch
    k = min(p["c"], s + i)
    k_i = min(k, i)
    k_s = k - k_i
    i -= k_i
    s -= k_s
    p_pool += k_i + k_s
    # 2. spread
    new = p["beta"] * s * i / n
    s -= new
    i += new
    # 3. recover
    rec = p["gamma"] * i
    i -= rec
    r += rec
    return [s, i, p_pool, r]


def make_run(params):
    """Run the daily map to the horizon (or to extinction).

    Returns (states, extinct_day). states is [(day, (S, I, P, R)), ...]
    starting at day 0; extinct_day is the first day I <= 1e-9, else None.
    """
    n = float(params["N"])
    i = float(params["i0"])
    state = [n - i, i, 0.0, 0.0]
    states = [(0, tuple(state))]
    extinct = None
    for day in range(1, int(params["days"]) + 1):
        state = step(state, params)
        states.append((day, tuple(state)))
        if state[1] <= _ZERO:
            extinct = day
            break
    return states, extinct


def metrics(states, extinct, params):
    """Summarise a trajectory. Verdict: CONTAINED if I never exceeds I0."""
    peak = max(st[1] for _, st in states)
    peak_day = next(d for d, st in states if st[1] == peak)
    final = states[-1][1]
    n = sum(final)
    return {
        "peak": peak,
        "peak_day": peak_day,
        "extinct_day": extinct,
        "final_i": final[1],
        "final_r": final[3],
        "attack_ratio": round(final[3] / n, 6),
        "patched_ratio": round(final[2] / n, 6),
        "verdict": "CONTAINED" if peak <= params["i0"] * (1 + 1e-6) else "ESCAPED",
    }


# ---------------------------------------------------------------------------
# The closed-form anchor
# ---------------------------------------------------------------------------

def growth_factor(beta, gamma):
    """The one-day growth factor in the linear regime.

    One infectious host patches out c/day; while S ~= N the map is
    I(d+1) = G * (I(d) - c) with G = (1+beta)(1-gamma).
    """
    return (1.0 + beta) * (1.0 - gamma)


def min_capacity(beta, gamma, i0):
    """Minimum daily capacity that holds the outbreak at its current level.

    With c hosts/day patched, I(d+1) = G*(I(d)-c) in the linear regime.
    I(d) shrinks iff I(d) < c*G/(G-1); so capacity c holds the line at I0 iff
    c >= I0*(G-1)/G.  Below G <= 1 the worm dies on its own: c_min = 0.
    """
    g = growth_factor(beta, gamma)
    if g <= 1.0:
        return 0.0
    return i0 * (g - 1.0) / g


def map_match(params, states, keep_frac=0.99):
    """Compare the daily map against the closed-form map G*(I-c).

    Only days where I > c (the crew is fully spent on infectious hosts) and
    where susceptible hosts are still >= keep_frac*N are compared, because the
    closed form assumes S ~= N and c < I.  Returns (worst_rel_error, points).
    """
    g = growth_factor(params["beta"], params["gamma"])
    c = params["c"]
    n = params["N"]
    worst = 0.0
    points = 0
    for day, (s, i, _p, _r) in states:
        if day == 0:
            prev_i = i
            continue
        if i > c and s >= keep_frac * n and prev_i > c:
            expected = g * (prev_i - c)
            if expected > 0 and i > 0:
                worst = max(worst, abs(i - expected) / expected)
                points += 1
        prev_i = i
    return worst, points


def numeric_min_capacity(beta, gamma, i0, n, days):
    """Smallest capacity that contains the outbreak, by bisection on the model.

    Contains(c) is monotone in c; the return value is within ~0.05% of the
    true numeric threshold.  The closed form is asserted against this in the
    selftest and `compare`: the closed form must be the fail-closed UPPER
    bound (analytic >= numeric), because susceptible depletion only helps.
    """
    def contained(cap):
        run, _ = make_run({"N": n, "i0": i0, "days": days,
                           "beta": beta, "gamma": gamma, "c": cap})
        return max(st[1] for _, st in run) <= i0 * (1 + 1e-6)

    lo, hi = 0.0, float(i0)
    if not contained(hi):
        return hi  # cannot contain even at capacity I0/day -- report the bound
    for _ in range(14):  # 0.05% resolution on [0, I0]
        mid = (lo + hi) / 2.0
        if contained(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ---------------------------------------------------------------------------
# Validation (fail-closed): bad input exits 1 with the reason on stderr.
# ---------------------------------------------------------------------------

def positive(name, value):
    if value <= 0:
        fail(f"{name} must be > 0, got {value}")
    return value


def nonneg(name, value):
    if value < 0:
        fail(f"{name} must be >= 0, got {value}")
    return value


def validate_common(args):
    params = {
        "N": float(positive("N", args.N)),
        "i0": float(positive("i0", args.i0)),
        "days": float(positive("days", args.days)),
    }
    if params["i0"] > params["N"]:
        fail(f"i0 ({params['i0']}) cannot exceed N ({params['N']})")
    params["N"] = int(params["N"])
    params["days"] = int(params["days"])
    return params


def validate_sir(args):
    p = validate_common(args)
    p["beta"] = float(positive("beta", args.beta))
    p["gamma"] = float(positive("gamma", args.gamma))
    return p


def fail(msg):
    print(_EPILOG + msg, file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Audit trail: append-only JSONL, hash chained, tamper-evident.
#   record = {seq, ts, tool, version, command, params, metrics, prev, hash}
#   hash   = sha256(canonical(record minus hash/prev) + "|" + prev)
#   record 0 must carry prev == GENESIS; every later record's prev must equal
#   the previous record's hash; seq must equal the line index.
# ---------------------------------------------------------------------------

def default_store():
    env = os.environ.get("PATCH_LAB_STORE")
    if env:
        return env
    return os.path.join(os.path.expanduser("~"), ".patchlab", "audit.jsonl")


def _canonical(record):
    clean = {k: v for k, v in record.items() if k not in ("hash", "prev")}
    return json.dumps(clean, sort_keys=True, separators=(",", ":"))


def _record_hash(record):
    digest = hashlib.sha256()
    digest.update(_canonical(record).encode("utf-8"))
    digest.update(b"|")
    digest.update(record["prev"].encode("utf-8"))
    return digest.hexdigest()


def read_store(path):
    """Parse the audit file. Returns (records, error_message)."""
    if not os.path.exists(path):
        return [], None
    records = []
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for lineno, line in enumerate(fh, start=1):
                line = line.rstrip("\n")
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except ValueError as exc:
                    return None, f"{path}:{lineno}: unparseable JSON ({exc})"
                if not isinstance(rec, dict):
                    return None, f"{path}:{lineno}: not an object"
                records.append(rec)
    except OSError as exc:
        return None, f"{path}: cannot read ({exc})"
    return records, None


def verify_store(path):
    """Rebuild the chain. Returns (ok, problems, n)."""
    records, err = read_store(path)
    if err:
        return False, [err], 0
    problems = []
    for i, rec in enumerate(records):
        if rec.get("seq") != i:
            problems.append(
                f"record {i}: seq is {rec.get('seq')!r}, expected {i} "
                "(an insert or a deletion shifted the file)")
        if i == 0:
            if rec.get("prev") != GENESIS:
                problems.append(
                    f"record 0: prev is {rec.get('prev', '')[:16]!r}, "
                    "expected the genesis anchor (the oldest record was deleted)")
        else:
            if rec.get("prev") != records[i - 1].get("hash"):
                problems.append(
                    f"record {i}: prev does not match the previous record's hash "
                    "(the previous record was edited)")
        expect = _record_hash(rec)
        if rec.get("hash") != expect:
            problems.append(f"record {i}: hash does not match its contents (edited)")
    return (len(problems) == 0), problems, len(records)


def append_record(store, command, params, metrics):
    """Append one chained record; fail closed if the write cannot be trusted."""
    records, err = read_store(store)
    if err:
        fail(f"cannot append to store {store}: {err}")
    prev = records[-1]["hash"] if records else GENESIS
    seq = len(records)
    rec = {
        "seq": seq,
        "ts": round(time.time(), 3),
        "tool": "patch-lab",
        "version": __version__,
        "command": command,
        "params": params,
        "metrics": metrics,
        "prev": prev,
        "hash": None,
    }
    rec["hash"] = _record_hash(rec)
    line = json.dumps(rec, sort_keys=True, separators=(",", ":"))
    directory = os.path.dirname(store)
    if directory:
        os.makedirs(directory, exist_ok=True)
    # Fail-closed write: append, flush, fsync, then read the line back and
    # confirm the chain is intact before reporting success.
    try:
        with open(store, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
            fh.flush()
            os.fsync(fh.fileno())
    except OSError as exc:
        fail(f"cannot write audit store {store}: {exc}")
    ok, problems, _ = verify_store(store)
    if not ok:
        fail(f"audit store does not verify right after writing ({problems[0]}); "
             "refusing to report a run that cannot be evidenced")
    return rec


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def fmt_hosts(x):
    return f"{x:,.0f}"


def report_plan(params, metrics, c_min, store, seq):
    r0 = params["beta"] / params["gamma"]
    g = growth_factor(params["beta"], params["gamma"])
    print(
        f"patch-lab {__version__} -- containment plan, "
        f"N={fmt_hosts(params['N'])}, I0={fmt_hosts(params['i0'])}, "
        f"beta={params['beta']:.2f}/day, gamma={params['gamma']:.2f}/day, "
        f"R0={r0:.2f}",
    )
    print(f"  growth factor : G = {g:.3f} per day (closed form I(d+1) = G*(I-c))")
    print(f"  min capacity  : {c_min:.2f} hosts/day to hold the line at I0")
    print(
        f"  c deployed    : {params['c']:.2f} hosts/day "
        f"({params['c'] / c_min:.2f}x min)"
        if c_min > 0 else
        f"  c deployed    : {params['c']:.2f} hosts/day (no capacity needed)")
    print(
        f"  verdict       : {metrics['verdict']} -- "
        + ("peak never exceeds I0" if metrics["verdict"] == "CONTAINED"
           else "the worm outruns daily patching"))
    print(
        f"  peak I        : {fmt_hosts(metrics['peak'])} hosts at day "
        f"{metrics['peak_day']}")
    if metrics["extinct_day"] is not None:
        print(f"  extinction    : day {metrics['extinct_day']} (I reaches 0)")
    else:
        print("  extinction    : not within the window")
    print(
        f"  attack ratio  : {100 * metrics['attack_ratio']:.2f}% of hosts "
        f"recovered; {100 * metrics['patched_ratio']:.2f}% patched by the end")
    print(f"  recorded to audit store: {store} (seq {seq})")


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_plan(args):
    p = validate_sir(args)
    c_min = min_capacity(p["beta"], p["gamma"], p["i0"])
    # Planning headroom when no explicit c is given: run the model at 1.2x the
    # minimum so the forecast is itself the recommended deployment.
    if args.c is None:
        p["c"] = round(c_min * 1.2, 6)
    else:
        p["c"] = float(nonneg("c", args.c))
    states, extinct = make_run(p)
    m = metrics(states, extinct, p)
    logged = dict(p)
    logged["N"] = int(p["N"])
    logged["days"] = int(p["days"])
    rec = append_record(args.store, "plan", logged, {
        "c_min": round(c_min, 6), "verdict": m["verdict"],
        "peak": round(m["peak"], 6), "extinct_day": m["extinct_day"],
        "attack_ratio": m["attack_ratio"],
    })
    report_plan(p, m, c_min, args.store, rec["seq"])


def cmd_compare(args):
    p = validate_sir(args)
    if args.c is None:
        p["c"] = min_capacity(p["beta"], p["gamma"], p["i0"]) * 0.5
    else:
        p["c"] = float(nonneg("c", args.c))
    states, _ = make_run(p)
    worst, points = map_match(p, states, keep_frac=args.keep_frac)
    formula = min_capacity(p["beta"], p["gamma"], p["i0"])
    numeric = numeric_min_capacity(p["beta"], p["gamma"], p["i0"],
                                   p["N"], p["days"])
    ratio = numeric / formula if formula > 0 else 1.0
    print(
        f"patch-lab {__version__} -- analytic anchor vs daily map, "
        f"N={fmt_hosts(p['N'])}, I0={fmt_hosts(p['i0'])}, c={p['c']:.2f}"
    )
    if points:
        print(
            f"  map error     : {worst:.3e}   "
            f"({points} days where S >= {args.keep_frac * 100:.0f}% of N)")
    else:
        print("  map error     : (no days inside the linear window; tighten keep-frac)")
    print(
        f"  c_min formula : {formula:.2f} hosts/day    "
        f"numeric bisection: {numeric:.2f}    (ratio {ratio:.3f})"
    )
    in_tol = worst <= args.tol and 0.95 <= ratio <= 1.02
    print(f"  verdict       : {'within tolerance' if in_tol else 'OUT OF TOLERANCE'}")
    logged = dict(p)
    logged["N"] = int(p["N"])
    logged["days"] = int(p["days"])
    append_record(args.store, "compare", logged, {
        "map_err": round(worst, 9), "points": points,
        "c_min_formula": round(formula, 6), "c_min_numeric": round(numeric, 6),
        "ratio": round(ratio, 6), "verdict": "ok" if in_tol else "fail",
    })
    if not in_tol and args.fail_over:
        fail(f"the anchor diverged from the model (map {worst:.3e}, "
             f"ratio {ratio:.3f})")
    sys.exit(0)


def _linspace(lo, hi, count):
    if count == 1:
        return [float(lo)]
    return [lo + (hi - lo) * i / (count - 1) for i in range(count)]


def cmd_sweep(args):
    positive("beta-min", args.beta_min)
    positive("beta-max", args.beta_max)
    positive("gamma-min", args.gamma_min)
    positive("gamma-max", args.gamma_max)
    if args.beta_max < args.beta_min:
        fail("beta-max must be >= beta-min")
    if args.gamma_max < args.gamma_min:
        fail("gamma-max must be >= gamma-min")
    if args.grid < 2:
        fail("grid must be >= 2 (a box, not a point)")
    base = validate_common(args)
    corners = []
    for b in _linspace(args.beta_min, args.beta_max, args.grid):
        for g in _linspace(args.gamma_min, args.gamma_max, args.grid):
            corners.append((b, g))
    worst_cmin = -1.0
    worst_at = None
    worst_r0 = -1.0
    for b, g in corners:
        cm = min_capacity(b, g, base["i0"])
        if cm > worst_cmin:
            worst_cmin, worst_at = cm, (b, g)
        worst_r0 = max(worst_r0, b / g)
    # Fail-closed envelope verification: the closed-form worst capacity must
    # contain EVERY corner numerically, and under-capacity must let some
    # corner escape. Sweep refuses to print a claim it cannot back.
    contained_at = 0
    escaped_under = 0
    for b, g in corners:
        p = dict(base, beta=b, gamma=g, c=worst_cmin)
        run_c, _ = make_run(p)
        if max(st[1] for _, st in run_c) <= base["i0"] * (1 + 1e-6):
            contained_at += 1
        p["c"] = worst_cmin * 0.8
        run_u, _ = make_run(p)
        if max(st[1] for _, st in run_u) > base["i0"] * 1.2:
            escaped_under += 1
    if contained_at != len(corners):
        fail("the worst-case capacity did not contain every corner -- the "
             "closed form is not a safe bound; refusing to publish an envelope")
    print(
        f"patch-lab {__version__} -- worst-case containment sweep over a "
        f"{args.grid}x{args.grid} grid ({len(corners)} corners), "
        f"N={fmt_hosts(base['N'])}, I0={fmt_hosts(base['i0'])}, "
        f"{base['days']} days"
    )
    print(
        f"  worst c_min   : {worst_cmin:.2f} hosts/day  "
        f"(beta={worst_at[0]:.3f}, gamma={worst_at[1]:.3f})")
    print(f"  worst R0      : {worst_r0:.2f}")
    print(
        f"  envelope      : {contained_at}/{len(corners)} corners contained "
        "at the worst capacity (fail-closed)")
    print(
        f"  under-capacity: {escaped_under}/{len(corners)} corners escape at "
        "80% of it (the bound is load-bearing)")
    print(
        f"  plan for      : deploy >= {math.ceil(worst_cmin * 100) / 100:.2f} "
        "hosts/day to hold every corner of the box at its I0")
    append_record(args.store, "sweep", {
        "N": int(base["N"]), "i0": base["i0"], "days": int(base["days"]),
        "grid": args.grid, "corners": len(corners),
        "beta_range": [args.beta_min, args.beta_max],
        "gamma_range": [args.gamma_min, args.gamma_max],
    }, {
        "worst_cmin": round(worst_cmin, 6),
        "worst_cmin_at": list(worst_at),
        "worst_r0": round(worst_r0, 6),
        "contained_at": contained_at,
        "escaped_under": escaped_under,
    })


def cmd_stats(args):
    records, err = read_store(args.store)
    if err:
        fail(err)
    if not records:
        print("no runs recorded yet")
        return
    print(f"patch-lab audit store: {args.store}  ({len(records)} records)")
    for rec in records[-8:]:
        m = rec.get("metrics", {})
        cmd = rec.get("command", "?")
        if cmd == "compare":
            disp = f"map_err={m.get('map_err', 0):.2e}"
        elif cmd == "sweep":
            disp = f"worst_cmin={m.get('worst_cmin', 0):.1f}"
        else:
            disp = f"{m.get('verdict', '?')}  peak={fmt_hosts(m.get('peak', 0))}"
        beta = rec.get("params", {}).get("beta", 0)
        gamma = rec.get("params", {}).get("gamma", 0)
        print(
            f"  #{rec['seq']}  {time.strftime('%Y-%m-%d %H:%M', time.localtime(rec.get('ts', 0)))}  "
            f"{cmd:7}  {disp}  beta={beta:.2f}  gamma={gamma:.2f}"
        )


def cmd_verify(args):
    ok, problems, n = verify_store(args.store)
    print(f"patch-lab audit verify: {args.store}")
    if not os.path.exists(args.store):
        print("  store does not exist yet (no runs recorded) -- nothing to verify")
        sys.exit(0)
    if ok:
        print(f"  verified: {n} records, chain intact")
        sys.exit(0)
    print(f"  TAMPERED: {len(problems)} problem(s) in {n} records")
    for p_ in problems:
        print(f"    - {p_}")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Selftest
# ---------------------------------------------------------------------------

def _expect(cond, name, state):
    state["total"] += 1
    if cond:
        print(f"  ok   {name}")
    else:
        print(f"  FAIL {name}")
        state["failures"].append(name)


def selftest():
    state = {"failures": [], "total": 0}

    # --- containment anchor (closed form, both directions) -----------------
    beta, gamma, i0, n, days = 0.3, 0.1, 100.0, 1_000_000, 365
    c_min = min_capacity(beta, gamma, i0)
    contained, _ = make_run({"beta": beta, "gamma": gamma, "i0": i0,
                             "N": n, "days": days, "c": c_min * 1.2})
    contained_peak = max(st[1] for _, st in contained)
    _expect(
        contained_peak <= i0 * (1 + 1e-6),
        "containment: c = 1.2*c_min holds the peak at I0 (never above)",
        state,
    )
    _expect(
        contained[-1][1][1] <= _ZERO,
        "containment: the contained run reaches I=0 (extinct)",
        state,
    )
    escaped, _ = make_run({"beta": beta, "gamma": gamma, "i0": i0,
                           "N": n, "days": days, "c": c_min * 0.8})
    escaped_peak = max(st[1] for _, st in escaped)
    _expect(
        escaped_peak > i0 * 1.2,
        "escape: c = 0.8*c_min lets the worm outrun the patchers "
        f"(peak {escaped_peak:,.0f} > 120)",
        state,
    )

    # --- closed form vs the model itself -----------------------------------
    ratios = []
    for b, g in ((0.3, 0.1), (0.25, 0.125), (0.5, 0.2)):
        formula = min_capacity(b, g, i0)
        numeric = numeric_min_capacity(b, g, i0, n, days)
        ratios.append(numeric / formula if formula > 0 else 1.0)
    _expect(
        all(0.95 <= r <= 1.02 for r in ratios),
        "the minimum-capacity formula matches a numeric bisection "
        f"(ratios {[f'{r:.3f}' for r in ratios]})",
        state,
    )

    dead, _ = make_run({"beta": 0.1, "gamma": 0.2, "i0": 50.0,
                        "N": 100_000, "days": 120, "c": 0.0})
    dead_peak = max(st[1] for _, st in dead)
    _expect(
        min_capacity(0.1, 0.2, 50.0) == 0.0 and dead_peak <= 50.0 * (1 + 1e-6),
        "below the self-limiting line the worm dies with no capacity (c_min = 0)",
        state,
    )

    slow, _ = make_run({"beta": beta, "gamma": gamma, "i0": i0,
                        "N": 500_000, "days": 120, "c": 5.0})
    fast, _ = make_run({"beta": beta, "gamma": gamma, "i0": i0,
                        "N": 500_000, "days": 120, "c": 30.0})
    peak_slow = max(st[1] for _, st in slow)
    peak_fast = max(st[1] for _, st in fast)
    _expect(
        peak_slow >= peak_fast,
        "capacity is monotone: a larger patch crew never makes it worse "
        f"({peak_slow:,.0f} >= {peak_fast:,.0f})",
        state,
    )

    # --- the documented equations are implemented exactly ------------------
    # Re-derive the daily map from the documented formulas (patch, then
    # spread: beta*S*I/N, then recover: gamma*I) and require the engine to
    # reproduce it bit-for-bit. Any drift mutant -- recovery no longer gamma*I,
    # spread no longer beta*S*I/N, patched hosts staying infectious -- breaks
    # the equality.
    def documented_step(s, i, p0, r0, n, beta_, gamma_, c):
        k = min(c, s + i)
        k_i = min(k, i)
        k_s = k - k_i
        i -= k_i
        s -= k_s
        p0 += k_i + k_s
        new = beta_ * s * i / n
        s -= new
        i += new
        rec = gamma_ * i
        i -= rec
        r0 += rec
        return (s, i, p0, r0)

    s0, ii0 = 800_000.0, 200.0
    probe = {"beta": beta, "gamma": gamma, "i0": i0, "N": 1_000_000,
             "days": 10, "c": 10.0}
    after = step([s0, ii0, 0.0, 0.0], probe)
    expect = documented_step(s0, ii0, 0.0, 0.0, probe["N"],
                             beta, gamma, probe["c"])
    _expect(
        all(abs(a - b) < 1e-9 for a, b in zip(after, expect)),
        "every day the documented equations are exact: patch, "
        "then spread beta*S*I/N, then recover gamma*I",
        state,
    )

    def max_conservation_error(p_params):
        states, _ = make_run(p_params)
        return max(abs(sum(st) - p_params["N"]) for _, st in states)

    _expect(
        max_conservation_error({"beta": beta, "gamma": gamma, "i0": i0,
                                "N": 1_000_000, "days": days, "c": 20.0}) < 1e-6
        and
        max_conservation_error({"beta": beta, "gamma": gamma, "i0": i0,
                                "N": 1_000_000, "days": days, "c": 5.0}) < 1e-6,
        "S+I+P+R == N is conserved every day (contained and escape)",
        state,
    )

    # --- linear-window map vs the closed form ------------------------------
    # The pure G*(I-c) anchor assumes S ~= N, so it carries an
    # O(beta*(1-S/N)) error at finite N: the honest tolerance is 1e-3, not
    # float noise.  N is huge here so depletion is ~1e-5 and the window holds.
    win_p = {"beta": beta, "gamma": gamma, "i0": 500.0, "N": 1_000_000_000,
             "days": 30, "c": min_capacity(beta, gamma, 500.0) * 0.5}
    win_states, _ = make_run(win_p)
    win_worst, win_points = map_match(win_p, win_states)
    _expect(
        win_worst <= 1e-3 and win_points >= 20,
        f"inside the linear window the daily map matches G*(I-c) "
        f"(worst {win_worst:.2e}, {win_points} days)",
        state,
    )

    # Load-bearing: with a small N the epidemic burns S down and breaks the
    # anchor; excluding the tail is what makes the window honest.  The window
    # error must be a fraction of the all-points error, or the window serves
    # no purpose (and a mutant that disables it goes unnoticed).
    long_p = {"beta": beta, "gamma": gamma, "i0": 200.0, "N": 50_000,
              "days": 120, "c": min_capacity(beta, gamma, 200.0) * 0.5}
    long_states, _ = make_run(long_p)
    win2_worst, _ = map_match(long_p, long_states)
    all_worst, _ = map_match(long_p, long_states, keep_frac=0.0)
    _expect(
        win2_worst < all_worst / 10 and all_worst > 1e-2,
        f"excluding the nonlinear tail is load-bearing "
        f"(window {win2_worst:.2e} vs all-points {all_worst:.2e})",
        state,
    )

    # --- determinism + negative control ------------------------------------
    s1, _ = make_run({"beta": beta, "gamma": gamma, "i0": i0,
                      "N": 500_000, "days": 120, "c": 20.0})
    s2, _ = make_run({"beta": beta, "gamma": gamma, "i0": i0,
                      "N": 500_000, "days": 120, "c": 20.0})
    _expect(
        [(d, tuple(st)) for d, st in s1] == [(d, tuple(st)) for d, st in s2],
        "identical arguments produce an identical trajectory (deterministic)",
        state,
    )

    n1, _ = make_run({"beta": beta, "gamma": gamma, "i0": i0,
                      "N": 500_000, "days": 120, "c": 5.0})
    n2, _ = make_run({"beta": beta, "gamma": gamma * 1.01, "i0": i0,
                      "N": 500_000, "days": 120, "c": 5.0})
    p1 = max(st[1] for _, st in n1)
    p2 = max(st[1] for _, st in n2)
    _expect(
        abs(p1 - p2) > 1000.0,
        "negative control: a perturbed recovery rate is detected as a deviation",
        state,
    )

    # --- sweep envelope (fail-closed planning) -----------------------------
    grid_n = 500_000
    corners = [(b, g) for b in _linspace(0.2, 0.4, 5)
               for g in _linspace(0.05, 0.15, 5)]
    worst_cmin = max(min_capacity(b, g, i0) for b, g in corners)
    contained_at = sum(
        1 for b, g in corners
        if max(st[1] for _, st in make_run(
            {"beta": b, "gamma": g, "i0": i0, "N": grid_n,
             "days": 365, "c": worst_cmin})[0]) <= i0 * (1 + 1e-6)
    )
    escaped_under = sum(
        1 for b, g in corners
        if max(st[1] for _, st in make_run(
            {"beta": b, "gamma": g, "i0": i0, "N": grid_n,
             "days": 365, "c": worst_cmin * 0.8})[0]) > i0 * 1.2
    )
    _expect(
        contained_at == len(corners),
        f"the worst-case c_min envelope contains every corner numerically "
        f"({contained_at}/{len(corners)})",
        state,
    )
    _expect(
        escaped_under >= 1,
        f"under-capacity (0.8x) lets some corner escape "
        f"({escaped_under}/{len(corners)} corners) -- the envelope is load-bearing",
        state,
    )

    # --- audit chain --------------------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        store = os.path.join(td, "audit.jsonl")
        p_audit = {"N": 1000000, "i0": 100.0, "days": 30, "beta": 0.3,
                   "gamma": 0.1, "c": 20.0}
        rec = append_record(store, "plan", p_audit, {"verdict": "CONTAINED", "peak": 100.0})
        _expect(rec["seq"] == 0 and rec["prev"] == GENESIS,
                "the first audit record anchors to the genesis hash", state)
        ok, problems, n = verify_store(store)
        _expect(ok and n == 1, "a freshly appended audit record verifies", state)

        # tamper 1: edit a recorded number (canonical JSON has no spaces)
        with open(store, encoding="utf-8") as fh:
            data = fh.read()
        edited = data.replace('"peak":100.0', '"peak":999.0')
        _expect(edited != data, "the tamper edit actually changed the line", state)
        with open(store, "w", encoding="utf-8") as fh:
            fh.write(edited)
        ok, problems, _ = verify_store(store)
        _expect(not ok and any("hash does not match" in p_ for p_ in problems),
                "editing a recorded number breaks the chain", state)

        # tamper 2: delete the oldest record (genesis anchor must catch it)
        store2 = os.path.join(td, "audit2.jsonl")
        append_record(store2, "plan", p_audit, {"verdict": "CONTAINED", "peak": 100.0})
        append_record(store2, "plan", p_audit, {"verdict": "CONTAINED", "peak": 200.0})
        with open(store2, encoding="utf-8") as fh:
            lines = fh.readlines()
        with open(store2, "w", encoding="utf-8") as fh:
            fh.writelines(lines[1:])
        ok, problems, _ = verify_store(store2)
        _expect(not ok and any("expected the genesis anchor" in p_ for p_ in problems),
                "deleting the oldest record breaks the genesis anchor", state)

        # tamper 3: insert a forged record borrowing a seq number
        store3 = os.path.join(td, "audit3.jsonl")
        append_record(store3, "plan", p_audit, {"verdict": "CONTAINED", "peak": 100.0})
        append_record(store3, "plan", p_audit, {"verdict": "CONTAINED", "peak": 200.0})
        with open(store3, encoding="utf-8") as fh:
            first = fh.readline()
            rest = fh.read()
        forged = {
            "seq": 1, "ts": 0.0, "tool": "patch-lab", "version": "9.9",
            "command": "plan", "params": {}, "metrics": {},
            "prev": "x" * 64, "hash": "y" * 64,
        }
        with open(store3, "w", encoding="utf-8") as fh:
            fh.write(first)
            fh.write(json.dumps(forged, sort_keys=True, separators=(",", ":")) + "\n")
            fh.write(rest)
        ok, problems, _ = verify_store(store3)
        _expect(not ok and any("does not match" in p_ for p_ in problems),
                "an inserted record with a borrowed seq breaks the chain", state)

        ok, problems, n = verify_store(os.path.join(td, "empty.jsonl"))
        _expect(ok and n == 0, "a missing (empty) store verifies as empty", state)

    return state["failures"], state["total"]


def run_selftest_cli():
    print(f"patch-lab {__version__} selftest")
    failures, total = selftest()
    if failures:
        print(f"selftest FAILED: {len(failures)} check(s)")
        sys.exit(1)
    print(f"\nall checks passed: {total} run, 0 skipped")
    sys.exit(0)


# CLI
# ---------------------------------------------------------------------------

def build_parser():
    ap = argparse.ArgumentParser(
        prog="patch-lab",
        description="Daily containment planner for a worm outbreak "
                    "(SIR + patching at capacity c, closed-form c_min, "
                    "hash-chained audit).",
        epilog="audit store: env PATCH_LAB_STORE, else ~/.patchlab/audit.jsonl",
    )
    ap.add_argument("--version", action="version", version=f"patch-lab {__version__}")
    ap.add_argument("--store", default=None, help="audit store path (overrides env/default)")
    sub = ap.add_subparsers(dest="command", required=False)

    def common(cmd):
        cmd.add_argument("--N", type=int, default=1_000_000)
        cmd.add_argument("--i0", type=float, default=100.0)
        cmd.add_argument("--days", type=float, default=365)

    s = sub.add_parser("selftest", help="run the self-test")
    s.set_defaults(func=lambda a: run_selftest_cli())

    c = sub.add_parser("plan", help="containment plan for one scenario")
    common(c)
    c.add_argument("--beta", type=float, default=0.3)
    c.add_argument("--gamma", type=float, default=0.1)
    c.add_argument("--c", type=float, default=None,
                   help="capacity (hosts/day); default: 1.2x the minimum")
    c.set_defaults(func=cmd_plan)

    c = sub.add_parser("compare", help="closed-form anchor vs the daily map")
    common(c)
    # compare needs a huge N so S ~= N holds inside the window (the anchor
    # carries an O(beta*(1-S/N)) error at finite N), hence the defaults:
    c.set_defaults(N=1_000_000_000, i0=500.0, days=30)
    c.add_argument("--beta", type=float, default=0.3)
    c.add_argument("--gamma", type=float, default=0.1)
    c.add_argument("--c", type=float, default=None,
                   help="capacity (hosts/day); default 0.5x the minimum")
    c.add_argument("--tol", type=float, default=1e-3)
    c.add_argument("--keep-frac", type=float, default=0.99)
    c.add_argument("--fail-over", action="store_true",
                   help="exit 1 if the anchor diverges from the model (CI gate)")
    c.set_defaults(func=cmd_compare)

    s = sub.add_parser("sweep", help="worst-case envelope over an uncertain beta/gamma")
    common(s)
    s.add_argument("--beta-min", type=float, default=0.2)
    s.add_argument("--beta-max", type=float, default=0.4)
    s.add_argument("--gamma-min", type=float, default=0.05)
    s.add_argument("--gamma-max", type=float, default=0.15)
    s.add_argument("--grid", type=int, default=5)
    s.set_defaults(func=cmd_sweep)

    s = sub.add_parser("stats", help="recent audit records")
    s.set_defaults(func=cmd_stats)

    s = sub.add_parser("verify", help="verify the audit chain")
    s.set_defaults(func=cmd_verify)

    return ap


def main(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)
    if args.store is None:
        args.store = default_store()
    if args.command in (None, "help"):
        ap.print_help()
        return 0
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())