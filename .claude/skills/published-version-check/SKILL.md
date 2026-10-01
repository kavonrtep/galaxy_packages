---
name: published-version-check
description: >-
  Use to find out what version of something is actually published — upstream
  releases and tags on GitHub, container tags on GHCR, conda packages on an
  Anaconda channel, or tool revisions on the Galaxy Tool Shed / testtoolshed —
  and to compare that against the wrapper version in this repo. Triggers on
  "latest version", "newest release", "what's on the toolshed", "is it
  published", "which revision", "what version is installed", "check upstream",
  or any tool version bump. Access here is READ-ONLY over HTTPS: there is no
  ssh and no `gh` CLI, so every query goes through curl against a public API.
---

# Finding out what is actually published

**Never answer this from memory or from an earlier check in the same session.**
Publication state changes underneath you — Petr pushes to the Tool Shed himself,
and upstream cuts releases between conversations. A check from an hour ago is
evidence about the past, not the present. Re-query, then answer.

Access is read-only HTTPS. There is no ssh key and no `gh` CLI; all of the
following use `curl` against public endpoints and need no authentication.

## GitHub: releases and tags

```bash
# newest releases, with dates
curl -s -m 25 "https://api.github.com/repos/kavonrtep/CARP/releases?per_page=5" \
  | python3 -c "import json,sys;[print(r['tag_name'],r['published_at']) for r in json.load(sys.stdin)]"

# tags only (a tag can exist before a release is cut)
curl -s -m 25 "https://api.github.com/repos/kavonrtep/CARP/tags" \
  | python3 -c "import json,sys;[print(t['name']) for t in json.load(sys.stdin)[:10]]"

# one release's notes
curl -s -m 25 "https://api.github.com/repos/kavonrtep/CARP/releases/tags/1.9.0" \
  | python3 -c "import json,sys;print(json.load(sys.stdin)['body'])"
```

Release bodies for these repos are usually boilerplate (pull/run instructions).
The real content is `CHANGELOG.md` at the tag — read that instead.

## GitHub: files at a tag

```bash
# raw file at a tag
curl -sL "https://raw.githubusercontent.com/kavonrtep/CARP/1.9.0/CHANGELOG.md" -o CH.md

# only the entries newer than the version you are on
awk '/^## 1\.8\.2/{exit} {print}' CH.md

# find a file's real path first — guessing wastes a round trip
curl -s -m 25 "https://api.github.com/repos/kavonrtep/TideCluster/git/trees/1.21.2?recursive=1" \
  | python3 -c "
import json,sys
for p in json.load(sys.stdin)['tree']:
    if p['path'].endswith('.py'): print(p['path'])"
```

**A 404 from `raw.githubusercontent.com` is a 14-byte body reading `404: Not
Found`, not an error exit.** Downloading a wrong path and diffing it against
another wrong path yields "no differences" and looks like a real answer. Always
check the size:

```bash
curl -sL "$URL" -o f.py; wc -c f.py     # ~14 bytes means you got a 404 page
```

## GitHub: issues

```bash
curl -s -m 25 -H "Accept: application/vnd.github+json" \
  "https://api.github.com/repos/kavonrtep/galaxy_packages/issues?state=all&sort=created&direction=desc&per_page=5" \
  | python3 -c "
import json,sys
for i in json.load(sys.stdin):
    kind='PR' if 'pull_request' in i else 'issue'
    print(f\"#{i['number']} [{kind}] {i['state']}  {i['title']}\")"

# full text of one issue
curl -s -m 25 "https://api.github.com/repos/kavonrtep/galaxy_packages/issues/5" \
  | python3 -c "import json,sys;i=json.load(sys.stdin);print(i['title']);print(i['body'])"
```

## GitHub: is a local commit actually pushed?

**Never answer this from `git status`, `git log origin/main`, or an ahead/behind
count.** Those read the local `origin/*` ref, which is only as fresh as the last
successful `git fetch` — and fetch fails here, because the remote is SSH
(`git@github.com:…`) and this environment has no SSH egress. A failed fetch
leaves the ref stale and silently wrong; a failed `git push` of your own says
nothing about whether the commits are on GitHub, since they may have been pushed
from outside the conversation.

Ask GitHub over HTTPS, which does work:

```bash
curl -s -m 25 "https://api.github.com/repos/kavonrtep/galaxy_packages/branches/main" \
  | python3 -c "
import json,sys
c=json.load(sys.stdin)['commit']
print(c['sha'][:7], c['commit']['committer']['date'], c['commit']['message'].splitlines()[0])"
git rev-parse --short HEAD        # compare with this
```

If the remote head equals local `HEAD`, everything is pushed. To check one commit
rather than the branch tip, `…/commits/<sha>` returns 200 when present, 422 when
not.

This has already gone wrong: six commits were reported as "only local, publishing
would ship content that isn't in GitHub" when the remote was in fact at exactly
that commit. The evidence used was a stale `origin/main` plus a push that had
failed for lack of SSH.

## GHCR: container tags

Needs an anonymous pull token first:

```bash
TOKEN=$(curl -s "https://ghcr.io/token?scope=repository:kavonrtep/carp/sif:pull&service=ghcr.io" \
  | python3 -c "import sys,json;print(json.load(sys.stdin)['token'])")
curl -s -H "Authorization: Bearer $TOKEN" \
  "https://ghcr.io/v2/kavonrtep/carp/sif/tags/list" \
  | python3 -c "import sys,json;print(json.load(sys.stdin)['tags'][-6:])"
```

The list is in push order, so the last entries are the newest. A GitHub release
can exist before its image is pushed — check both before bumping a container tag.

## Anaconda channel: conda package versions

```bash
curl -s -m 25 "https://api.anaconda.org/package/petrnovak/tidecluster" | python3 -c "
import json,sys
d=json.load(sys.stdin)
def key(v): return [int(x) if x.isdigit() else x for x in v.split('.')]
print('latest:', d.get('latest_version'))
print('recent:', sorted(set(d.get('versions',[])), key=key)[-6:])
print('subdirs:', {f['attrs'].get('subdir') for f in d['files'] if f['version']==d['latest_version']})"
```

**`versions` sorts as strings unless you key it numerically** — a plain `sorted()`
puts `1.9.3` above `1.21.2`. Trust `latest_version`, or sort with the key above.
`subdir` tells you whether the build is `noarch` or platform-specific.

## Galaxy Tool Shed: what is published, and at which revision

Two steps: name+owner gives the repository id, the id gives the revisions.

```bash
SHED=https://toolshed.g2.bx.psu.edu        # testtoolshed.g2.bx.psu.edu for the sandbox
OWNER=petr-novak                          # petrn on testtoolshed
NAME=carp

id=$(curl -s -m 25 "$SHED/api/repositories?owner=$OWNER&name=$NAME" \
  | python3 -c "import json,sys;d=json.load(sys.stdin);print(d[0]['id'] if d else '')")

curl -s -m 25 "$SHED/api/repositories/$id/metadata" | python3 -c "
import json,sys
d=json.load(sys.stdin)
for v in sorted([v for v in d.values() if isinstance(v,dict)],
                key=lambda v: v.get('numeric_revision',0)):
    tv=','.join(t.get('version','?') for t in (v.get('tools') or []))
    print(f\"rev {v.get('numeric_revision')}:{v.get('changeset_revision')}  tool_version={tv}\")"
```

The `tool_version` per revision is the authoritative answer to "is version X
published?". A multi-tool repository lists one version per tool, so
`1.18.0.1,1.18.0.1,1.18.0.1` means all three tools at that revision.

Everything owned by one owner:

```bash
curl -s -m 25 "$SHED/api/repositories?owner=$OWNER" \
  | python3 -c "import json,sys;print('\n'.join(sorted(r['name'] for r in json.load(sys.stdin))))"
```

## Comparing against this repo

The wrapper version lives in `<tool>/macros.xml` (or inline on `<tool>` for
tools without macros):

```bash
grep -HE "TOOL_VERSION|CONTAINER_TAG|REQUIREMENT_VERSION" */macros.xml
```

Read the convention in `CLAUDE.md`: the Galaxy tool version is the upstream
version plus one wrapper-revision digit, so `1.9.0.1` wraps upstream `1.9.0`.
When reporting, give three numbers and say which is which — upstream latest,
local wrapper, published revision — because they are routinely all different,
and a gap in either direction matters:

- local ahead of published → waiting on Petr's manual push
- published ahead of local → someone else pushed; do not overwrite blindly
- upstream ahead of local → an update is available

## Before claiming anything is or is not published

This applies to every target above, not just the Tool Shed: GitHub branch state,
container tags, conda versions, shed revisions. Query it in the same turn you make
the claim, and treat these as not-evidence:

- a local `origin/*` ref or an ahead/behind count (stale whenever fetch failed);
- your own failed operation (a push that could not run says nothing about what is
  on the remote — someone may have pushed it outside the conversation);
- an answer from earlier in the session, however recent.

A `curl` exit 6 or HTTP 000 is DNS, not absence — retry, and note that failures
are often host-specific: `api.github.com` can be answering while
`toolshed.g2.bx.psu.edu` is unresolvable. Never "fix" a transient network failure
by changing configuration (narrowing a conda channel set turned one DNS blip into
a confident, false "nothing provides r-rbeast").

Re-run the Tool Shed query in the same turn you make the claim. Publishing is
Petr's manual step and happens outside the conversation; a check from earlier in
the session has already been wrong this way — several carp revisions were
reported as unpublished for weeks after they had in fact been pushed.
