"""Regenerate every table and figure in the manuscript.

    python studies/run_all.py            # everything
    python studies/run_all.py s3 s6      # selected studies only
    python studies/run_all.py figures    # figures only, from cached JSON

All outputs land in results/tables, results/figures and results/raw, each
carrying a provenance header that records the git commit, the seed, and the
statement that the contents are simulated.
"""
from __future__ import annotations

import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from opcem import report                                        # noqa: E402

STUDIES = {
    "s0": ("front-end validation", "s0_frontend"),
    "s1": ("transport calibration", "s1_transport"),
    "s2": ("descriptor fidelity and forecast skill", "s2_descriptors"),
    "s3": ("event-aligned estimator and power", "s3_events"),
    "s4": ("transport-correction level", "s4_transport_effect"),
    "s5": ("projection bias in mobility", "s5_projection"),
    "s6": ("beam perturbation index", "s6_beam"),
    "s7": ("representativeness and sensitivity", "s7_representativeness"),
    "s8": ("closed-loop control", "s8_control"),
}


def run_study(key: str) -> float:
    label, mod_name = STUDIES[key]
    print(f"\n=== {key.upper()}  {label} " + "=" * (48 - len(label)))
    t0 = time.time()
    mod = __import__(mod_name)
    if key == "s4":
        out = {"s4": mod.run_s4(), "s4b": mod.run_s4b()}
    else:
        out = mod.run()
    dt = time.time() - t0
    _print(out)
    print(f"--- {key} done in {dt:.0f} s")
    return dt


def _print(out, indent: int = 0) -> None:
    import pandas as pd
    pad = " " * indent
    if isinstance(out, pd.DataFrame):
        print(pad + out.to_string(index=False).replace("\n", "\n" + pad))
    elif isinstance(out, dict):
        for k, v in out.items():
            print(f"{pad}[{k}]")
            _print(v, indent + 2)
    else:
        print(pad + str(out))


def main(argv) -> int:
    report.ensure_dirs()
    keys = [a for a in argv if a in STUDIES]
    want_figs = ("figures" in argv) or not argv
    if not keys and not argv:
        keys = list(STUDIES)
    total = 0.0
    failures = []
    for k in keys:
        try:
            total += run_study(k)
        except Exception:
            failures.append(k)
            traceback.print_exc()
    if want_figs:
        print("\n=== FIGURES " + "=" * 48)
        t0 = time.time()
        try:
            import figures
            figures.make_all()
            print(f"--- figures done in {time.time() - t0:.0f} s")
        except Exception:
            failures.append("figures")
            traceback.print_exc()
    print(f"\nTotal study wall time: {total:.0f} s  "
          f"(git {report.git_sha()})")
    if failures:
        print(f"FAILED: {', '.join(failures)}")
        return 1
    print("All requested outputs regenerated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
