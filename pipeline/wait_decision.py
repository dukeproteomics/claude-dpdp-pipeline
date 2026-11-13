#!/usr/bin/env python3
"""Wait until the analyst sends a decision from a review page, then print it and exit.

Claude runs this in the background after opening a step's page; it is notified when it exits.
Usage: python3 pipeline/wait_decision.py --run=<run_dir> --step=04_noise [--timeout=7200]
Exit 0 with the decision JSON on stdout; exit 3 on timeout; exit 4 if a pending decision is for a different step.
Requests no step command records ("other" free text, "plots", "fix", "stop") are moved to review/decisions_done/ once read
(and noted in decisions.jsonl as step "request"); the rest stay until the step command that applies them archives them.
"""
import argparse, json, os, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
ap.add_argument("--step", required=True)
ap.add_argument("--timeout", type=float, default=7200)
a = ap.parse_args()
TAKEN_ON_READ = {"other", "plots", "fix", "stop"}
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
        if rec.get("step") != a.step:
            sys.exit(4)
        # Requests that no step command records (free text, plots, fix, stop) are set aside now that Claude has read them,
        # so the next wait doesn't return the same one again (log_decision / page_args archive the others when they run).
        if rec.get("choice") in TAKEN_ON_READ:
            done = os.path.join(a.run, "review", "decisions_done")
            os.makedirs(done, exist_ok=True)
            os.replace(path, os.path.join(done, time.strftime("%Y%m%d_%H%M%S") + "_" + a.step + "_" + rec["choice"] + ".json"))
            line = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "step": "request", "analyst": rec.get("analyst") or os.environ.get("USER", ""),
                    "settings": {"step": a.step, "choice": rec.get("choice"), "fields": rec.get("fields") or {}},
                    "note": rec.get("note") or "", "source": "review page"}
            with open(os.path.join(a.run, "decisions.jsonl"), "a") as f:
                f.write(json.dumps(line) + "\n")
        sys.exit(0)
    time.sleep(1)
print(json.dumps({"timeout": True, "step": a.step}))
sys.exit(3)
