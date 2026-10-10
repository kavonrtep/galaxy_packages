#!/usr/bin/env python3
"""Regenerate the classification-vocabulary help blocks in carp/macros.xml.

Both CARP tools document which classes may appear in a library header. That list
is CARP's, not ours - `classification_vocabulary.yaml` in the image - so it is
generated from the pinned @CONTAINER_TAG@ instead of being retyped, and it lives
in macros.xml as tokens so the two tools cannot drift apart.

Run it whenever @CONTAINER_TAG@ in carp/macros.xml changes:

    python3 scripts/update_carp_vocabulary_help.py            # rewrite the tokens
    python3 scripts/update_carp_vocabulary_help.py --check     # exit 1 if stale

--check is the guard against the help quietly describing an older vocabulary
than the image the tool actually runs. It reads the vocabulary of the pinned tag
and compares; it does not go looking for newer CARP releases.

The vocabulary is read from the matching git tag on GitHub, because the container
tag and the release tag are the same string and the image is built from that tag.
That keeps this a text download rather than a 4 GB image pull; pass
--from-container to read the copy inside the image instead, which is the stronger
check if you suspect the image and the tag disagree (set SINGULARITY_CACHEDIR
first, or the pull lands in the default cache).
"""
import argparse
import re
import subprocess
import sys
import urllib.request
from pathlib import Path
from xml.sax.saxutils import escape

import yaml

REPO = Path(__file__).resolve().parent.parent
MACROS = REPO / "carp" / "macros.xml"
IMAGE = "oras://ghcr.io/kavonrtep/carp/sif:{tag}"
VOCAB_IN_IMAGE = "/opt/pipeline/classification_vocabulary.yaml"
RAW = ("https://raw.githubusercontent.com/kavonrtep/CARP/{tag}/"
       "classification_vocabulary.yaml")
INDENT = "  "
# Every generated line carries this prefix, so a token placed at column 0 in a
# <help> block lands as a correctly indented RST literal block. Indenting in the
# help instead would only indent the token's first line.
BASE = "    "


def container_tag(text):
    m = re.search(r'<token name="@CONTAINER_TAG@">([^<]+)</token>', text)
    if not m:
        sys.exit("no @CONTAINER_TAG@ in carp/macros.xml")
    return m.group(1).strip()


def read_vocabulary(tag, from_container=False, local=None, attempts=3):
    if local:
        body = Path(local).read_bytes()
        return yaml.safe_load(body), f"{local}"
    if from_container:
        out = subprocess.run(
            ["singularity", "exec", IMAGE.format(tag=tag), "cat", VOCAB_IN_IMAGE],
            capture_output=True, timeout=1800, check=True).stdout
        if not out.strip():
            sys.exit(f"{VOCAB_IN_IMAGE} is empty in the {tag} image")
        return yaml.safe_load(out), f"{tag} image"
    last = None
    for i in range(attempts):
        try:
            with urllib.request.urlopen(RAW.format(tag=tag), timeout=60) as fh:
                body = fh.read()
            break
        except Exception as exc:               # DNS and TLS failures included
            last = exc
            print(f"attempt {i + 1}/{attempts} failed: {exc}", file=sys.stderr)
    else:
        sys.exit(f"could not read the vocabulary from GitHub: {last}\n"
                 f"pass --vocabulary with a local copy, or --from-container")
    if len(body) < 1000:            # raw.githubusercontent's 404 page is 14 bytes
        sys.exit(f"no vocabulary at tag {tag} on GitHub ({len(body)} bytes)")
    return yaml.safe_load(body), f"github tag {tag}"


def path_list(paths):
    """One class per line, written in full, parents before children.

    Full paths rather than an indented tree: the string on the line is exactly
    what goes after the '#' in a FASTA header, so it can be copied as it stands,
    and the hierarchy is still legible from the shared prefixes.
    """
    return [BASE + p for p in paths]


def render(vocab):
    classes = list(vocab["classifications"])
    special = vocab["special_classes"]

    rdna = [c for c in classes if c.split("/")[0] == "rDNA"]
    blocks = {}
    blocks["CLASS_LIST"] = path_list(classes)

    tandem = ["Satellite", "Satellite/<name>"] + rdna
    lines = path_list(tandem)
    lines[1] += "      any name you choose, e.g. Satellite/FabTR_PisTR-B"
    blocks["CLASS_LIST_TANDEM"] = lines

    blocks["CLASS_SPECIAL"] = [
        BASE + "%-16s %s%s" % (name, body.get("description", ""),
                        "" if not body.get("accepts_subpath") else
                        " (a sub-name may be appended, e.g. %s/<name>)" % name)
        for name, body in special.items()]

    blocks["VOCABULARY_SOURCE"] = [
        "CARP %%s (vocabulary version %s)" % vocab.get("version")]
    return blocks


def splice(text, name, lines):
    """Replace one generated token body, keeping the markers."""
    begin = f"<!-- BEGIN GENERATED {name} -->"
    end = f"<!-- END GENERATED {name} -->"
    if begin not in text or end not in text:
        sys.exit(f"markers for {name} not found in carp/macros.xml")
    head, rest = text.split(begin, 1)
    _, tail = rest.split(end, 1)
    # The rendered text goes into an XML text node, and it contains '<' - the
    # <name> placeholder for a satellite name. Unescaped, that is a start tag and
    # macros.xml stops parsing.
    body = escape("\n".join(lines))
    return f"{head}{begin}\n    <token name=\"@{name}@\">{body}</token>\n    {end}{tail}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if the tokens are out of date")
    ap.add_argument("--from-container", action="store_true",
                    help="read the vocabulary inside the image rather than from "
                         "the git tag (pulls the image if it is not cached)")
    ap.add_argument("--vocabulary", metavar="FILE",
                    help="read a local copy of classification_vocabulary.yaml "
                         "instead of fetching it; for an offline run, state in "
                         "the commit where the copy came from")
    a = ap.parse_args()

    text = MACROS.read_text()
    tag = container_tag(text)
    vocab, origin = read_vocabulary(tag, a.from_container, a.vocabulary)
    blocks = render(vocab)
    blocks["VOCABULARY_SOURCE"] = [blocks["VOCABULARY_SOURCE"][0] % tag]

    new = text
    for name, lines in blocks.items():
        new = splice(new, name, lines)

    counts = ", ".join(f"{n}={len(v)} line(s)" for n, v in blocks.items())
    if a.check:
        if new == text:
            print(f"up to date with CARP {tag} (read from {origin}): {counts}")
            return 0
        print(f"OUT OF DATE: carp/macros.xml does not match CARP {tag} "
              f"(read from {origin}). Run without --check to regenerate.",
              file=sys.stderr)
        return 1
    MACROS.write_text(new)
    print(f"carp/macros.xml updated from CARP {tag} ({origin}): {counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
