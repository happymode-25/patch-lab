# patch-lab — build report

**Version:** 1.0.0 · **Repo:** https://github.com/happymode-25/patch-lab (private) ·
**Built:** Oct 2026 · **Basis:** a containment planner for the worm-lab
outbreak scenarios, to the egress-broker/worm-lab bar (single file, stdlib
only, fail-closed, tamper-evident, mutation-swept, CI-gated).

## 1. What was built

`patch-lab.py` is a one-file, zero-dependency (Python 3.9+ stdlib) daily
containment planner. A worm spreads through `N` hosts with SIR dynamics while
a patch crew deploys up to `c` hosts/day, infectious first. The tool answers
the planner's question — *what is the minimum capacity that holds the
outbreak?* — with a closed form, verifies that closed form against the model,
and sweeps the uncertainty box for the worst case. Every run is appended to a
hash-chained JSONL audit store.

## 2. Security-engineering stance

Same three fail-closed guarantees as egress-broker and worm-lab:

- **Fail-closed input.** `i0 > N`, `β ≤ 0`, `γ ≤ 0`, `days ≤ 0`, `c < 0`,
  unknown commands: refused with an exit code and a reason. Nothing is
  silently clamped. (`c = 0` stays legal: "no patching capacity" is a real
  scenario, not an error.)
- **Fail-closed planning.** `sweep` reports the *worst* `c_min` across an
  uncertainty box — and verifies numerically, before printing anything, that
  deploying it contains every corner and that 80% of it lets some corner
  escape. If the envelope does not hold, sweep exits 1: a planner gets no
  number the model cannot back.
- **Fail-closed evidence.** If the audit chain does not verify *immediately
  after a write*, the tool refuses to report the run.

## 3. The model, the closed form, and how both are proved

```
patch    k = min(c, S + I)   infectious first (k_I = min(k, I))
spread   new = β·S·I/N
recover  rec = γ·I
G = (1+β)(1−γ)                          (linear-regime daily growth factor)
I(d+1) = G·(I(d) − c)                   (closed-form map while S ≈ N)
c_min = I₀·(G−1)/G,  or 0 when G ≤ 1    (hold the line at I₀)
```

The closed form is the load-bearing analytic result: give the crew `c ≥ c_min`
hosts/day and the outbreak shrinks from day one and never exceeds `I₀`; give
them less and it grows without bound until susceptibles run out.

Selftest checks the properties that matter, not the code path:

| Property | How it is checked |
|---|---|
| Containment anchor | `c = 1.2·c_min` keeps peak ≤ I₀ and goes extinct; `c = 0.8·c_min` peaks > 1.2·I₀ |
| Formula vs model | closed-form `c_min` vs numeric bisection on the model, ratio in [0.95, 1.02] |
| Exact equations | one day re-derived independently from the documented map, bit-for-bit equal — recovery is exactly γ·I, spread exactly β·S·I/N |
| Conservation | S+I+P+R == N every day (contained and escape) |
| Map vs closed form | inside the linear window the daily map matches `G·(I−c)` to 1e-3, and the window is *load-bearing*: on a nonlinear run, excluding late points changes the error by an order of magnitude |
| Determinism | identical arguments ⇒ byte-identical trajectories; a 1% perturbed recovery rate is detected |
| Worst-case sweep | the envelope contains every grid corner; 80% of it lets corners escape |

## 4. Audit trail

`{seq, ts, tool, version, command, params, metrics, prev, hash}` with
`hash = sha256(canonical(record minus hash/prev) + "|" + prev)`; record 0 is
anchored to a genesis sentinel. `verify` checks seq contiguity, prev linkage,
and content hash *independently* — all problems are reported, not just the
first. Editing a number, deleting the oldest record, and inserting a forged
record are all caught and named.

Honest limitation: trimming the tail of the file cannot be detected from
inside the file. External anchoring (log shipper, SIEM, append-only media) is
the documented answer if that matters.

## 5. Defects found and fixed during bring-up

| Defect | Symptom | Fix |
|---|---|---|
| exact-equation check derived the wrong expected spread | it used I₀ instead of post-patch I, so the check could not pass | the check re-derives the full documented day independently and compares bit-for-bit |
| pure `G·(I−c)` anchor assumed S ≈ N too strongly | windowed map error ~1e-3 at finite N, so the float-noise tolerance failed | honest 1e-3 tolerance at N=10⁹, and a *relative* load-bearing check (window ≪ all-points) |
| the load-bearing check compared against an absolute 1e-6 | a small-N window legitimately carries ~1e-3 of anchor drift | `win < all/10 and all > 1e-2` — the window must change the verdict by an order of magnitude |
| selftest passed iterables instead of states | `make_run(...)` returns `(states, extinct)`; the envelope checks unpacked the tuple and crashed after 12 checks | index `[0]` before max |
| e2e `invalid_i0` used `--i0 999999` with default `N=1e6` | i0 < N, so the "violation" was legal and the test failed | force `--i0 5000000 --N 1000000` |
| mutation canary anchored to a line that does not exist | `growth_factor` returns the expression, never binds `g` | canary targets `g = growth_factor(beta, gamma)` instead |

## 6. Verification — what the green is, and that it can fail

| Gate | Result |
|---|---|
| `ci/import_gate.py` | 0 non-stdlib imports across the repo |
| `selftest` | 21/21 checks |
| `tests/test_patch_lab.py` | 27/27 tests, real subprocess runs of the CLI |
| `ci/prove_suites_can_fail.py` | recovery-drift mutation and injected-check mutation both turn the suites red; restore is green again |
| `ci/mutation.py` | 10 real mutants killed (map, formula, window, envelope, audit, cli), semantically-identical canary survives, live file hash unchanged (scratch-copy only) |
| CI matrix | Python 3.9 / 3.11 / 3.13 × linux / macos / windows, count floors: 21 selftest checks, 27 e2e tests |

The suite is not trusted on a green run alone: `prove_suites_can_fail` breaks
the code on purpose and requires the failures to be *named checks*, not
crashes, and the mutation sweep requires the canary to survive — otherwise
every verdict in the run is void.

## 7. Reproduce

```powershell
git clone https://github.com/happymode-25/patch-lab   # private; needs access
cd patch-lab
python ci/import_gate.py
python patch-lab.py selftest
python tests/test_patch_lab.py
python ci/prove_suites_can_fail.py
python ci/mutation.py
```

## 8. Relationship to the previous builds

The trilogy is one discipline across three claims: egress-broker proves that
egress decisions and audit rows are trustworthy (prevent), worm-lab proves the
same for propagation numbers and their scenario record (model), and patch-lab
proves the containment answer — *this many hosts/day holds the line, and here
is the bound to plan for when you do not know the parameters* — against the
model it sits on. The patterns carried over deliberately: a separate pure
selftest plus a subprocess e2e suite, count floors in CI, an import gate, a
prove-suites-can-fail harness, a canary-anchored mutation sweep on a scratch
copy with the live file hashed before and after, a `.gitattributes` with
`eol=lf` from the first commit, and a README that states guarantees *and*
their limits.