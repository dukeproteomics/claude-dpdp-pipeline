#!/usr/bin/env python3
"""Wait until the analyst sends a decision from a review page, then print it and exit.

Claude runs this in the background after opening a step's page; it is notified when it exits.
Usage: python3 pipeline/wait_decision.py --run=<run_dir> --step=04_noise [--timeout=7200]
Exit 0 with the decision JSON on stdout; exit 3 on timeout; exit 4 if a pending decision is for a different step.
"""
import argparse, json, os, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
ap.add_argument("--step", required=True)
ap.add_argument("--timeout", type=float, default=7200)
a = ap.parse_args()
path = os.path.join(a.run, "review", "pending_decision.json")
end = time.time() + a.timeout
while time.time() < end:
    if os.path.exists(path):
        try:
            rec = json.load(open(path))
        except Exception:
            time.sleep(0.5)
            continue
        print(json.dumps(rec, indent=1))
        sys.exit(0 if rec.get("step") == a.step else 4)
    time.sleep(1)
print(json.dumps({"timeout": True, "step": a.step}))
sys.exit(3)
