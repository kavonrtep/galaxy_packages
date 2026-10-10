#!/usr/bin/env python3
"""PostToolUse hook: parse any tool XML just written, and say so when it broke.

Malformed tool XML is cheap to make and expensive to notice: planemo lint
reports only "Error parsing file <path>" with no line, and a Galaxy server skips
the tool. Three causes have accounted for every occurrence in this repository:

  * '--' inside an XML comment, which is illegal. It happens because the natural
    way to document a wrapper is to name the flags it passes.
  * a stray character left after an element, e.g. '>' or '"' on its own line.
  * an unescaped '<' or '&' in generated text.

Reads the PostToolUse payload on stdin, and is silent unless the file is a .xml
inside the repository and fails to parse.
"""
import json
import os
import sys
import xml.dom.minidom

REPO = "/mnt/ssd/galaxy_packages"

HINTS = (
    ("invalid token",
     "the usual cause here is '--' inside an XML comment, which is illegal; "
     "a stray '>' or '\"' left after an element is the other"),
    ("mismatched tag", "an element is not closed, or is closed in the wrong order"),
    ("not well-formed", "check for an unescaped '<' or '&' in text content"),
)


def main():
    try:
        payload = json.load(sys.stdin)
    except Exception:
        return 0
    path = ((payload.get("tool_response") or {}).get("filePath")
            or (payload.get("tool_input") or {}).get("file_path") or "")
    if not path.endswith(".xml") or not os.path.realpath(path).startswith(REPO):
        return 0
    if not os.path.exists(path):
        return 0
    try:
        xml.dom.minidom.parse(path)
    except Exception as exc:
        rel = os.path.relpath(path, REPO)
        text = str(exc)
        message = f"{rel} is not well-formed XML: {text}"
        for needle, hint in HINTS:
            if needle in text:
                message += f" -- {hint}"
                break
        print(json.dumps({
            "systemMessage": message,
            "hookSpecificOutput": {
                "hookEventName": "PostToolUse",
                "additionalContext": message + ". Fix it before running planemo.",
            },
        }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
