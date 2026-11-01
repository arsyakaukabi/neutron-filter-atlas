"""Check the documentation (README.md and wiki/*.md): links, anchors, file paths and page form.

- every relative link resolves to a file in the repository, and every #anchor to a heading of the target page
  (GitHub's rule: lower case, punctuation removed, spaces to hyphens);
- every repository path written in backticks (configs/..., src/..., pipelines/..., scripts/..., tools/..., tests/...,
  wiki/..., data/...) exists in the repository or, for data files, in data/MANIFEST.csv;
- no absolute local path (drive letters, home directories) and no marker of unfinished work;
- every wiki page starts with a title and, within its first ten lines, "What this page is for:" and
  "Who should read it:", and has at most MAX_LINES lines (README.md at most README_LINES).
Usage: python tools/docs/check_docs.py      (exit 1 on any problem)
"""
import csv
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MAX_LINES, README_LINES = 300, 25
LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+)\)")
CODE_PATH = re.compile(r"`((?:configs|src|pipelines|scripts|tools|tests|wiki|data)/[^`\s*<>]+)`")
ABSOLUTE = re.compile(r"(?<![\w/])[A-Za-z]:[\\/]|/home/\w|/Users?/\w|~/\w")
MARKER = re.compile(r"\bTODO\b|\bFIXME\b|\bXXX\b")


def slug(heading):
    s = heading.strip().lower()
    s = re.sub(r"[^\w\- ]", "", s)
    return s.replace(" ", "-")


def anchors(path):
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    out, seen = set(), {}
    for h in re.findall(r"^#{1,6}\s+(.+?)\s*$", text, re.M):
        s = slug(h.replace("`", ""))
        n = seen.get(s, 0)
        out.add(s if n == 0 else f"{s}-{n}")
        seen[s] = n + 1
    return out


def manifest_paths():
    m = REPO / "data" / "MANIFEST.csv"
    if not m.exists():
        return set()
    paths = {r["path"] for r in csv.DictReader(open(m, encoding="utf-8"))}
    dirs = {"/".join(p.split("/")[:i]) for p in paths for i in range(1, p.count("/") + 1)}
    return paths | dirs


def strip_code_blocks(text):
    return re.sub(r"```.*?```", "", text, flags=re.S)


def main():
    pages = [REPO / "README.md"] + sorted((REPO / "wiki").glob("*.md"))
    known = manifest_paths()
    problems = []
    for page in pages:
        rel = page.relative_to(REPO).as_posix()
        text = page.read_text(encoding="utf-8")
        lines = text.splitlines()
        limit = README_LINES if page.name == "README.md" and page.parent == REPO else MAX_LINES
        if len(lines) > limit:
            problems.append(f"{rel}: {len(lines)} lines (limit {limit})")
        if page.parent.name == "wiki":
            head = [l for l in lines[:10] if l.strip()]
            if not (head and head[0].startswith("# ") and any(l.startswith("What this page is for:") for l in head)
                    and any(l.startswith("Who should read it:") for l in head)):
                problems.append(f"{rel}: missing title or the two orientation lines")
        prose = strip_code_blocks(text)
        for target in LINK.findall(prose):
            if re.match(r"[a-z]+:", target):
                continue
            file_part, _, anchor = target.partition("#")
            dest = (page.parent / file_part).resolve() if file_part else page
            if not dest.exists():
                problems.append(f"{rel}: broken link {target}")
                continue
            if anchor and dest.suffix == ".md" and anchor not in anchors(dest):
                problems.append(f"{rel}: no heading for #{anchor} in {dest.relative_to(REPO).as_posix()}")
        for p in CODE_PATH.findall(text):
            p = p.rstrip(".,;:)")
            if "NN" in p:   # placeholder such as pipelines/NN_name.py
                continue
            if p.startswith("data/") and p.rstrip("/") in known:
                continue
            if not (REPO / p).exists():
                problems.append(f"{rel}: path `{p}` does not exist")
        for n, line in enumerate(lines, 1):
            if ABSOLUTE.search(line):
                problems.append(f"{rel}:{n}: absolute local path")
            if MARKER.search(line):
                problems.append(f"{rel}:{n}: unfinished-work marker")
    print(f"checked {len(pages)} pages: {len(problems)} problems")
    for p in problems:
        print("  " + p)
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
