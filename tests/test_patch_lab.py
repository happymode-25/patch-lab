"""patch-lab end-to-end suite.

Every test drives the real CLI as a subprocess (the tool is a single
hyphenated file, so importing is deliberately not how we test it). The store
for each test lives in a fresh temp directory, so runs are isolated.

The suite asserts outcomes the user can observe: exit codes, report fields,
fail-closed behaviour on bad input, the closed-form c_min anchor being
reported and load-bearing, and the hash-chained audit trail detecting every
class of tamper we know about.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

TOOL = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "..", "patch-lab.py"))


def run_cli(args, store=None, extra_env=None):
    env = dict(os.environ)
    if store is not None:
        env["PATCH_LAB_STORE"] = store
    if extra_env:
        env.update(extra_env)
    return subprocess.run(
        [sys.executable, TOOL] + [str(a) for a in args],
        capture_output=True, text=True, env=env,
        cwd=os.path.dirname(TOOL), timeout=120,
    )


def read_store(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return [json.loads(l) for l in fh if l.strip()]


def field(pattern, text):
    """Extract one `label : value` style field, commas stripped."""
    return re.search(pattern, text).group(1).replace(",", "")


class PatchLabE2E(unittest.TestCase):

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.store = os.path.join(self._td.name, "audit.jsonl")

    # --- core commands ------------------------------------------------------

    def test_selftest_green(self):
        r = run_cli(["selftest"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("all checks passed", r.stdout)
        self.assertNotIn("FAIL", r.stdout)

    def test_plan_contained_with_sufficient_capacity(self):
        # beta=0.3, gamma=0.1, I0=100 -> c_min ~14.53; c=18 must hold the line.
        r = run_cli(["plan", "--beta", "0.3", "--gamma", "0.1", "--i0", "100",
                     "--N", "1000000", "--c", "18"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("containment plan", r.stdout)
        self.assertIn("verdict       : CONTAINED", r.stdout)
        self.assertIn("peak I        : 100 hosts at day 0", r.stdout)
        self.assertIn("extinction    :", r.stdout)
        self.assertIn("recorded to audit store", r.stdout)

    def test_plan_reports_min_capacity(self):
        # Without --c the plan must surface the analytic c_min for the scenario.
        r = run_cli(["plan", "--beta", "0.3", "--gamma", "0.1", "--i0", "100"],
                    store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        c_min = float(field(r"min capacity\s+: ([\d.]+)", r.stdout))
        self.assertTrue(14.0 < c_min < 15.0, c_min)  # 100*(1.17-1)/1.17 = 14.53
        self.assertIn("CONTAINED", r.stdout)  # planning headroom 1.2x holds

    def test_plan_under_capacity_escapes(self):
        # c=10 < c_min: the worm must outrun the crew and the report must say so.
        r = run_cli(["plan", "--beta", "0.3", "--gamma", "0.1", "--i0", "100",
                     "--N", "1000000", "--c", "10"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("verdict       : ESCAPED", r.stdout)
        peak = float(field(r"peak I\s+: ([\d,]+)", r.stdout))
        self.assertGreater(peak, 120.0)  # > 1.2*I0 as the selftest requires

    def test_plan_zero_capacity_is_valid(self):
        # c=0 is a legitimate scenario (no patching at all), not an error.
        r = run_cli(["plan", "--c", "0"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("verdict       : ESCAPED", r.stdout)

    def test_plan_dead_worm_cmin_zero(self):
        # (1+beta)(1-gamma) = 0.88 <= 1: the worm dies on its own, no capacity.
        r = run_cli(["plan", "--beta", "0.1", "--gamma", "0.2", "--i0", "50",
                     "--N", "100000"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("min capacity  : 0.00", r.stdout)
        self.assertIn("no capacity needed", r.stdout)
        self.assertIn("verdict       : CONTAINED", r.stdout)

    def test_plan_deterministic_two_runs(self):
        r1 = run_cli(["plan", "--c", "18"], store=self.store)
        r2 = run_cli(["plan", "--c", "18"], store=self.store)
        self.assertEqual(r1.returncode, 0)
        self.assertEqual(r2.returncode, 0)
        # The audit line legitimately differs (seq increments); the *plan
        # report* must be byte-identical (trajectory is deterministic).
        body = lambda r: [l for l in r.stdout.splitlines()
                          if "recorded to audit store" not in l]
        self.assertEqual(body(r1), body(r2))

    # --- fail-closed input handling ----------------------------------------

    def test_invalid_i0_fails_closed(self):
        r = run_cli(["plan", "--i0", "5000000", "--N", "1000000"],
                    store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("cannot exceed N", r.stderr)

    def test_negative_beta_fails_closed(self):
        r = run_cli(["plan", "--beta", "-0.5"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("beta must be > 0", r.stderr)

    def test_zero_beta_fails_closed(self):
        r = run_cli(["plan", "--beta", "0"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("beta must be > 0", r.stderr)

    def test_zero_gamma_fails_closed(self):
        # No recovery at all is division-by-nonsense for the model; refused.
        r = run_cli(["plan", "--gamma", "0"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("gamma must be > 0", r.stderr)

    def test_negative_c_fails_closed(self):
        r = run_cli(["plan", "--c", "-1"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("c must be >= 0", r.stderr)

    def test_unknown_command_exits_2(self):
        r = run_cli(["frobnicate"], store=self.store)
        self.assertEqual(r.returncode, 2)

    def test_version_flag(self):
        r = run_cli(["--version"], store=self.store)
        self.assertEqual(r.returncode, 0)
        self.assertIn("patch-lab 1.0.0", r.stdout)

    # --- analytic anchor ----------------------------------------------------

    def test_compare_green_defaults(self):
        r = run_cli(["compare"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("verdict       : within tolerance", r.stdout)
        self.assertIn("map error     :", r.stdout)
        self.assertIn("numeric bisection", r.stdout)

    def test_compare_fail_over_really_fails(self):
        # Negative control through the real CLI: small N burns S down and the
        # pure G*(I-c) anchor breaks; --fail-over MUST go red.
        r = run_cli(["compare", "--N", "50000", "--i0", "200", "--days", "120",
                     "--keep-frac", "0.0", "--fail-over"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("diverged", r.stderr)
        self.assertIn("OUT OF TOLERANCE", r.stdout)

    # --- worst-case sweep ---------------------------------------------------

    def test_sweep_green(self):
        r = run_cli(["sweep", "--grid", "2"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("worst c_min   :", r.stdout)
        self.assertIn("plan for      :", r.stdout)

    def test_sweep_worst_bounds_every_corner(self):
        # Fail-closed planning: the envelope (worst c_min over the box) must be
        # at least as pessimistic as any single corner, AND it must actually
        # contain every corner in the numeric model (9/9 on a 3x3 grid).
        r = run_cli(["sweep", "--grid", "3", "--beta-min", "0.2", "--beta-max",
                     "0.4", "--gamma-min", "0.05", "--gamma-max", "0.15"],
                    store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        worst = float(field(r"worst c_min\s+: ([\d.]+)", r.stdout))
        self.assertIn("envelope      : 9/9 corners contained", r.stdout)
        for b, g in ((0.2, 0.15), (0.4, 0.05), (0.3, 0.1)):
            s = run_cli(["plan", "--beta", str(b), "--gamma", str(g)],
                        store=self.store)
            self.assertEqual(s.returncode, 0, s.stderr)
            single = float(field(r"min capacity\s+: ([\d.]+)", s.stdout))
            self.assertGreaterEqual(worst, single,
                                    f"envelope {worst} below corner {single}")

    # --- audit trail --------------------------------------------------------

    def test_run_creates_audit_record(self):
        run_cli(["plan"], store=self.store)
        records = read_store(self.store)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["command"], "plan")
        self.assertEqual(records[0]["seq"], 0)
        self.assertEqual(records[0]["tool"], "patch-lab")
        self.assertEqual(len(records[0]["hash"]), 64)

    def test_verify_clean_after_runs(self):
        run_cli(["plan"], store=self.store)
        run_cli(["compare"], store=self.store)
        r = run_cli(["verify"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("verified: 2 records, chain intact", r.stdout)

    def test_verify_empty_missing_store(self):
        r = run_cli(["verify"], store=self.store)
        self.assertEqual(r.returncode, 0)
        self.assertIn("nothing to verify", r.stdout)

    def test_verify_catches_middle_edit(self):
        run_cli(["plan"], store=self.store)
        run_cli(["compare"], store=self.store)
        run_cli(["plan"], store=self.store)
        with open(self.store, encoding="utf-8") as fh:
            data = fh.read()
        edited = data.replace('"version":"1.0.0"', '"version":"9.9.9"', 1)
        self.assertNotEqual(edited, data)
        with open(self.store, "w", encoding="utf-8") as fh:
            fh.write(edited)
        r = run_cli(["verify"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("TAMPERED", r.stdout)
        self.assertIn("hash does not match", r.stdout)

    def test_verify_catches_prefix_delete(self):
        run_cli(["plan"], store=self.store)
        run_cli(["plan"], store=self.store)
        with open(self.store, encoding="utf-8") as fh:
            lines = fh.readlines()
        with open(self.store, "w", encoding="utf-8") as fh:
            fh.writelines(lines[1:])
        r = run_cli(["verify"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("genesis anchor", r.stdout)

    def test_verify_catches_forged_insert(self):
        run_cli(["plan"], store=self.store)
        run_cli(["plan"], store=self.store)
        with open(self.store, encoding="utf-8") as fh:
            first = fh.readline()
            rest = fh.read()
        forged = {
            "seq": 1, "ts": 0.0, "tool": "patch-lab", "version": "9.9",
            "command": "plan", "params": {}, "metrics": {},
            "prev": "x" * 64, "hash": "y" * 64,
        }
        with open(self.store, "w", encoding="utf-8") as fh:
            fh.write(first)
            fh.write(json.dumps(forged, sort_keys=True, separators=(",", ":")) + "\n")
            fh.write(rest)
        r = run_cli(["verify"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("TAMPERED", r.stdout)

    def test_run_refused_on_tampered_store(self):
        # Fail-closed: once the store no longer verifies, new runs must be
        # refused rather than silently appending to a broken chain.
        run_cli(["plan"], store=self.store)
        with open(self.store, encoding="utf-8") as fh:
            data = fh.read()
        with open(self.store, "w", encoding="utf-8") as fh:
            fh.write(data.replace('"version":"1.0.0"', '"version":"0.0.1"', 1))
        r = run_cli(["plan"], store=self.store)
        self.assertEqual(r.returncode, 1)
        self.assertIn("audit store does not verify", r.stderr)

    def test_env_store_location_respected(self):
        alt = os.path.join(self._td.name, "elsewhere", "run.jsonl")
        r = run_cli(["plan"], store=alt)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.path.exists(alt))
        self.assertEqual(len(read_store(alt)), 1)

    def test_stats_lists_runs(self):
        run_cli(["plan"], store=self.store)
        run_cli(["compare"], store=self.store)
        r = run_cli(["stats"], store=self.store)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("#0", r.stdout)
        self.assertIn("#1", r.stdout)
        self.assertIn("plan", r.stdout)
        self.assertIn("compare", r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)