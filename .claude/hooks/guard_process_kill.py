#!/usr/bin/env python3
"""PreToolUse hook: refuse pattern-based process kills, warn about pattern greps.

`pkill -f planemo` matches the shell that issues it, because that shell's own
command line contains the pattern. It has killed the issuing session repeatedly
in this project - once in the middle of clearing a full disk. `pgrep -f` has the
same self-match but is merely misleading: it reports its own pipeline as a hit,
which sends the reader chasing a PID that no longer exists.

Kill by PID instead, resolved from something that is not a string you just
typed - the port, or the parent of the tree:

    P=$(ss -ltnpH 'sport = :9090' | grep -oE 'pid=[0-9]+' | cut -d= -f2 | sort -u)
    kill $P
"""
import json
import re
import sys

DENY = (
    (re.compile(r"\bpkill\b[^|;&\n]*(-f|--full)\b"),
     "pkill -f matches the shell that runs it, so this can kill the session"),
    (re.compile(r"\bpkill\b\s+(-\w+\s+)*[A-Za-z]"),
     "pkill matches by name and can hit unrelated processes"),
    (re.compile(r"\bkillall\b"),
     "killall matches by name and can hit unrelated processes"),
)
WARN = re.compile(r"\bpgrep\b[^|;&\n]*(-f|--full)\b")

SAFE_FORM = (
    "Resolve the PID from something other than a string you just typed, then "
    "kill that:\n"
    "  by port:   P=$(ss -ltnpH 'sport = :9090' | grep -oE 'pid=[0-9]+' | "
    "cut -d= -f2 | sort -u); kill $P\n"
    "  by argv:   ps -u \"$(id -u)\" -o pid=,args= | PAT=<pattern> awk "
    "'ENVIRON[\"PAT\"] != \"\" && index($0, ENVIRON[\"PAT\"]) {print $1}'\n"
    "  (pattern via the environment, so it is not in the command line being "
    "searched; set PAT on awk, not on ps - an empty pattern matches every "
    "process. Print the PIDs and count them before killing.)"
)


def main():
    try:
        payload = json.load(sys.stdin)
    except Exception:
        return 0
    command = (payload.get("tool_input") or {}).get("command") or ""
    for pattern, why in DENY:
        if pattern.search(command):
            print(json.dumps({
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": f"{why}. {SAFE_FORM}",
                },
            }))
            return 0
    if WARN.search(command):
        note = ("pgrep -f matches its own pipeline, so one of the PIDs it prints "
                "is the search itself and changes on every run. Exclude it, or "
                "pass the pattern through the environment "
                "(ps -u \"$(id -u)\" -o pid=,args= | PAT=<pattern> awk "
                "'ENVIRON[\"PAT\"] != \"\" && index($0, ENVIRON[\"PAT\"])'; "
                "PAT goes on awk, not on ps).")
        print(json.dumps({
            "systemMessage": note,
            "hookSpecificOutput": {
                "hookEventName": "PreToolUse",
                "additionalContext": note,
            },
        }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
