"""Interleaved A/B micro-benchmark harness.

This machine is noisy — identical code has varied 40% between runs on file-based
benchmarks.  Three rules make a measurement here mean something, and this module
exists so they are not re-litigated per benchmark:

1.  **In-memory only.**  ``BytesIO``, never a file, so the page cache and the
    filesystem are out of the loop.
2.  **Min-of-N, N >= 7.**  The minimum is the least contaminated sample; the mean
    measures the noise, not the code.
3.  **Interleaved in one process.**  Baseline and candidate alternate, and the
    order flips each round, so a machine that drifts during the run drifts
    through both arms equally.

Anything under ~5% is noise, and :func:`ab` says so rather than reporting a
number that reads as a result.

Usage::

    from bench_ab import ab
    ab("N-drop filter", baseline_callable, candidate_callable)
"""
from __future__ import annotations

import gc
import time

#: Below this, a difference is not distinguishable from machine noise.
NOISE_FLOOR_PCT = 5.0


def _time_once(fn, reps: int) -> float:
    gc.collect()
    gc.disable()
    try:
        t0 = time.perf_counter()
        for _ in range(reps):
            fn()
        return time.perf_counter() - t0
    finally:
        gc.enable()


def ab(name: str, baseline, candidate, *, rounds: int = 9, reps: int = 1,
       quiet: bool = False) -> dict:
    """Measure *candidate* against *baseline*, interleaved, min-of-*rounds*.

    Returns a dict with the two minima, the speedup as a multiplier, the change
    as a percentage, and whether the result clears the noise floor.
    """
    if rounds < 7:
        raise ValueError("rounds must be >= 7; fewer does not survive this machine")

    base_times: list[float] = []
    cand_times: list[float] = []
    for i in range(rounds):
        # Flip the order every round so systematic drift cancels rather than
        # accruing to whichever arm happens to run second.
        if i % 2 == 0:
            b = _time_once(baseline, reps)
            c = _time_once(candidate, reps)
        else:
            c = _time_once(candidate, reps)
            b = _time_once(baseline, reps)
        base_times.append(b)
        cand_times.append(c)

    base = min(base_times)
    cand = min(cand_times)
    speedup = base / cand if cand else float("inf")
    pct = (speedup - 1.0) * 100.0
    significant = abs(pct) >= NOISE_FLOOR_PCT

    result = {
        "name": name,
        "baseline_s": base,
        "candidate_s": cand,
        "speedup": speedup,
        "pct": pct,
        "significant": significant,
    }
    if not quiet:
        verdict = f"{pct:+.1f}%" if significant else f"{pct:+.1f}% (NOISE)"
        print(f"  {name:<52s} {base * 1e3:8.2f} ms -> {cand * 1e3:8.2f} ms  "
              f"{speedup:5.2f}x  {verdict}")
    return result


def header(title: str) -> None:
    print(f"\n{title}")
    print("  " + "-" * 92)
    print(f"  {'benchmark':<52s} {'baseline':>11s}    {'candidate':>8s}  "
          f"{'ratio':>5s}  change")
    print("  " + "-" * 92)
