# patch-lab

A daily containment planner for a worm outbreak, built to the same bar as
[egress-broker](https://github.com/happymode-25/egress-broker) and
[worm-lab](https://github.com/happymode-25/worm-lab):

- **One file, zero dependencies.** `patch-lab.py` runs on Python 3.9+ with the
  standard library alone. The daily model, the closed form, the envelope sweep
  — no numpy, no scipy, nothing to pip-install.
- **Fail-closed by default.** Invalid input is refused with an exit code and a
  reason; `sweep` refuses to publish an envelope it cannot verify; the audit
  store refuses to accept a run it cannot evidence.
- **Tamper-evident audit trail.** Every run is appended to a hash-chained JSONL
  store. `verify` rebuilds the chain and names exactly what is wrong if
  anything was edited, deleted, or forged.
- **Provably testable.** A 21-check selftest, a 27-test end-to-end suite, a
  `prove_suites_can_fail` harness, and a 10-mutant + canary mutation sweep —
  all enforced in CI with count floors.

## What it does

A worm is spreading through `N` hosts (SIR dynamics). A patch crew can deploy
up to `c` hosts/day, **infectious hosts first**. The daily map is explicit and
deterministic — patch, then spread, then recover:

```
patch    k = min(c, S + I)   infectious first
spread   new = β·S·I/N
recover  rec = γ·I
```

The question a planner actually asks is: **what is the minimum daily capacity
that holds the outbreak at its current level?** That has a closed form. With
`G = (1+β)(1−γ)` (the one-day growth factor net of recovery), the linear-regime
map is `I(d+1) = G·(I(d) − c)`, which shrinks iff `I < c·G/(G−1)`. So:

```
c_min = I₀·(G − 1)/G        (G > 1;  c_min = 0 when G ≤ 1, the worm dies alone)
```

Deploy at or above `c_min` and the outbreak shrinks from its first day and
never exceeds `I₀`. Deploy below it and it outruns the crew. The comparison is
not believed on faith: `compare` verifies the closed form against the model
itself (numeric bisection) and the model against the closed form
(linear-window match), and `sweep` verifies its own worst-case envelope
numerically before it prints a single number.

## When to use it

- As a **containment planner**: given `β`, `γ`, `I₀`, `N`, what capacity holds
  the line, and what happens at your actual crew size?
- As a **fail-closed estimator** when the parameters are *not* known: `sweep`
  runs the whole uncertainty box and reports the worst `c_min` to plan for,
  and refuses to print the envelope unless it contains every corner.
- As an **auditable record**: every run lands in the hash-chained store, so
  the capacity you committed to can be proved later to be the capacity you
  committed to.

## Quick start

```powershell
python patch-lab.py selftest                 # 21 checks, proves the model + anchor

python patch-lab.py plan --beta 0.3 --gamma 0.1 --i0 100 --N 1000000 --c 18
python patch-lab.py plan --beta 0.3 --gamma 0.1 --i0 100   # reports c_min, plans at 1.2x

python patch-lab.py compare                  # closed form vs daily map
python patch-lab.py compare --fail-over      # CI-style gate: exit 1 out of tol

python patch-lab.py sweep --beta-min 0.2 --beta-max 0.4 `
                        --gamma-min 0.05 --gamma-max 0.15 --grid 5

python patch-lab.py stats                    # recent runs from the audit store
python patch-lab.py verify                   # prove the chain is untampered
```

The audit store lives at `~/.patchlab/audit.jsonl`, or wherever
`PATCH_LAB_STORE` points, or where `--store` says.

## Commands

| Command | What it gives you |
|---|---|
| `selftest` | 21 checks: containment anchor both directions, formula vs bisection, exact equations, conservation, linear-window match, determinism, sweep envelope, audit tamper detection |
| `plan [--c X]` | Minimum capacity `c_min`, then a run at your `c` (or 1.2× `c_min`): verdict CONTAINED/ESCAPED, peak, extinction day, attack ratio |
| `compare [--fail-over]` | Daily map vs `G·(I−c)`, and the closed form vs numeric bisection; `--fail-over` gates on tolerance |
| `sweep` | Worst-case `c_min` over an uncertain β/γ box, with the envelope verified numerically — plan for the bound, not the guess |
| `stats` | Recent runs from the audit store |
| `verify` | Rebuild the hash chain; exit 1 and name the tamper if anything is wrong |

Everything is deterministic: the same arguments produce the same numbers, and
`--version` is available everywhere.

## Fail-closed behaviour

- `i0 > N`, `β ≤ 0`, `γ ≤ 0`, `days ≤ 0`, `c < 0`, unknown commands → exit 1/2
  with a reason, no partial run, no clamped silence. (`c = 0` is *valid*: no
  patching at all is a legitimate scenario.)
- `sweep` computes the worst `c_min` from the closed form, then verifies
  numerically that deploying it contains *every* corner of the grid; if it
  does not, sweep exits 1 rather than print a bound it cannot back.
- If the audit chain does not verify *immediately after a write*, the tool
  refuses to report the run rather than record one it cannot evidence.

## Audit trail

Each record is `{seq, ts, tool, version, command, params, metrics, prev, hash}`
where `hash = sha256(canonical(record minus hash/prev) + "|" + prev)`. Record 0
must carry the genesis anchor. `verify` checks three independent things per
record — seq contiguity, prev linkage, and content hash — so it reports *all*
problems, not just the first. Editing a number, deleting the oldest record,
and inserting a forged record are all caught and named.

What it cannot do, stated honestly: trimming the *tail* of the file is not
detectable from inside the file. If that matters, anchor the last hash
externally (log shipper, SIEM, append-only media).

## Verification

| Gate | What it proves | Floor |
|---|---|---|
| `ci/import_gate.py` | zero third-party imports anywhere | — |
| `selftest` | model, closed form, envelope, audit, determinism | 21 checks |
| `tests/test_patch_lab.py` | real CLI behaviour over subprocesses | 27 tests |
| `ci/prove_suites_can_fail.py` | the suites fail when the code breaks and pass when restored | — |
| `ci/mutation.py` | 10 real mutants killed, semantically-identical canary survives | 11/11 |

Run them locally from the repo root:

```powershell
python ci/import_gate.py
python patch-lab.py selftest
python tests/test_patch_lab.py
python ci/prove_suites_can_fail.py
python ci/mutation.py
```

The mutation sweep works on a scratch copy — the live `patch-lab.py` is hashed
before and after, so a sweep can never leave a mutant behind in the repo.