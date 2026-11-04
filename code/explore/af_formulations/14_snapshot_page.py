#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# code/explore/af_formulations/14_snapshot_page.py
#
# Freeze earlier versions of the AF formulations page so readers can open them
# from the page's Version menu, and write CHANGELOG.md.
#
# docs/af_formulations/versions.json is the source of truth (id, date, ref =
# git tag or commit, title, plain-language changes). For every version other than `latest`, this
# copies the page as it was at that commit (index.html, data.js and methods.md
# when present) from git into docs/af_formulations/versions/<id>/ and:
#   * adds the shared version menu and an "earlier version" banner (versions.js);
#   * closes comments (hides the buttons and blocks the giscus loader), so new
#     comments always land on the live page's threads;
#   * marks the copy noindex.
# Existing snapshots are left alone (frozen); --force rebuilds them. The live
# page is never touched. Idempotent.
#
# To release a new version: tag the commit (git tag v0.2), add an entry at the
# top of versions.json with "ref": "v0.2", move "latest" to it, re-run this
# script (it freezes the previous version from its tag), commit, push with tags.
#
# Run from the repo root:  python3 code/explore/af_formulations/14_snapshot_page.py [--force]
# ------------------------------------------------------------------------------
import json, subprocess, sys
from pathlib import Path

DOCS = Path("docs/af_formulations")
FILES = ["index.html", "data.js", "methods.md"]
GISCUS = "https://giscus.app/client.js"


def git_show(commit, path):
    r = subprocess.run(["git", "show", f"{commit}:{path}"], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def snapshot(v, force):
    out = DOCS / "versions" / v["id"]
    if (out / "index.html").exists() and not force:
        return "kept"
    ref = v.get("ref") or v["commit"]
    html = git_show(ref, DOCS / "index.html")
    if html is None:
        sys.exit(f"{v['id']}: no docs/af_formulations/index.html at {ref}")
    out.mkdir(parents=True, exist_ok=True)
    for f in FILES[1:]:
        b = git_show(ref, DOCS / f)
        if b is not None:
            (out / f).write_bytes(b)
    s = html.decode("utf-8")
    s = s.replace(GISCUS, "about:blank")  # comments closed on archived copies
    head = '<meta name="robots" content="noindex">\n<style>.cbtn, #cdrawer { display: none !important; }</style>\n'
    tail = f'<script src="../../versions.js" data-version="{v["id"]}" data-root="../../"></script>\n'
    assert s.count("</head>") == 1 and s.count("</body>") == 1, f"{v['id']}: unexpected page structure"
    s = s.replace("</head>", head + "</head>", 1).replace("</body>", tail + "</body>", 1)
    (out / "index.html").write_text(s, encoding="utf-8")
    return "written"


def changelog(V):
    L = ["# AF formulations page: change history", "",
         "Generated from `docs/af_formulations/versions.json` by `code/explore/af_formulations/14_snapshot_page.py`; edit the JSON, not this file.",
         "Each version links to its frozen copy on the page and to the git tag it was released from.", ""]
    for v in V["versions"]:
        latest = v["id"] == V["latest"]
        page = "./" if latest else f"versions/{v['id']}/"
        L.append(f"## {v['id']} · {v['title']} · {v['date']}{' (latest)' if latest else ''}")
        L.append("")
        ref = v.get("ref") or v["commit"]
        code = f"[`{ref}`](https://github.com/CC4A-SSA/cc4adapt/tree/{ref})"
        L.append(f"Code and page source: {code} · [page]({page}){' · **results changed**' if v.get('results_changed') else ''}")
        L.append("")
        L.append(v["summary"])
        L.append("")
        L += [f"- {c}" for c in v["changes"]]
        L.append("")
    (DOCS / "CHANGELOG.md").write_text("\n".join(L), encoding="utf-8")


def main():
    if not DOCS.is_dir():
        sys.exit("Run from the cc4adapt repo root")
    force = "--force" in sys.argv
    V = json.loads((DOCS / "versions.json").read_text(encoding="utf-8"))
    ids = [v["id"] for v in V["versions"]]
    assert V["latest"] in ids and len(set(ids)) == len(ids), "versions.json: latest missing or duplicate ids"
    for v in V["versions"]:
        if v["id"] != V["latest"]:
            print(f"{v['id']:5} {v.get('ref') or v['commit']}  {snapshot(v, force)}")
    changelog(V)
    print("CHANGELOG.md written")


if __name__ == "__main__":
    main()
