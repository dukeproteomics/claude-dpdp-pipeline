#!/usr/bin/env python3
"""PreToolUse hook (Bash): a pipeline command that applies an analyst decision (--decided_by=user) may only run if the
analyst has answered since the last such command. One decision per answer: Claude cannot chain several applies, or
apply before the analyst has said anything.

"Answered" = a user message typed by the analyst, or an AskUserQuestion answer, in the session transcript, or a decision for
this step sent from the review page (<run>/review/pending_decision.json, written by pipeline/review_server.py).
Exit 2 blocks the tool call; the message on stderr goes back to Claude.
"""
import json, os, re, sys


def page_decision_for(cmd):
    """True if the command's run has a pending review-page decision for the command's step."""
    m = re.search(r"--db=(\S+)", cmd)
    if not m:
        return False
    run_dir = os.path.dirname(os.path.abspath(m.group(1).strip("'\"")))
    s = re.search(r"steps/(\d\d)_\w+\.R", cmd) or re.search(r"--step=(\d\d)_\w+", cmd)
    try:
        rec = json.load(open(os.path.join(run_dir, "review", "pending_decision.json")))
    except Exception:
        return False
    return bool(s) and str(rec.get("step", "")).startswith(s.group(1) + "_")

def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        return 0
    if data.get("tool_name") != "Bash":
        return 0
    cmd = (data.get("tool_input") or {}).get("command", "")
    if "--decided_by=user" not in cmd or "pipeline/" not in cmd:
        return 0
    if page_decision_for(cmd):
        return 0
    path = data.get("transcript_path")
    if not path:
        return 0
    ask_ids, events = set(), []   # events: ("answer"|"apply", order)
    try:
        with open(path) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                msg = rec.get("message") or {}
                content = msg.get("content")
                if rec.get("type") == "assistant" and isinstance(content, list):
                    for c in content:
                        if c.get("type") == "tool_use":
                            if c.get("name") == "AskUserQuestion":
                                ask_ids.add(c.get("id"))
                            elif c.get("name") == "Bash" and "--decided_by=user" in (c.get("input") or {}).get("command", ""):
                                events.append("apply")
                elif rec.get("type") == "user":
                    if rec.get("isMeta") or rec.get("isCompactSummary"):
                        continue
                    if isinstance(content, str) and content.strip():
                        events.append("answer")
                    elif isinstance(content, list):
                        for c in content:
                            if c.get("type") == "text" and c.get("text", "").strip() and not c.get("text", "").startswith("<"):
                                events.append("answer")
                            elif c.get("type") == "tool_result" and c.get("tool_use_id") in ask_ids:
                                events.append("answer")
    except OSError:
        return 0
    # The current command is not in the transcript yet. Block if an apply already happened after the last answer.
    last_answer = max((i for i, e in enumerate(events) if e == "answer"), default=-1)
    if any(e == "apply" for e in events[last_answer + 1:]) or last_answer == -1:
        sys.stderr.write("DPDP: this command applies an analyst decision (--decided_by=user), but the analyst has not answered since "
                         "the last decision was applied. Show the step's review page, ask the analyst (AskUserQuestion), and run it "
                         "only after they choose.\n")
        return 2
    return 0

if __name__ == "__main__":
    sys.exit(main())
