"""EXFOR: fetch X4 entries from the IAEA exfor_master GitHub mirror and parse (N,TOT) data into tidy rows."""
import re
import urllib.request
from pathlib import Path

import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.paths import D

XC = D("exfor")
URL = cfg("benchmark")["exfor_url"]
E_UNITS = {"EV": 1.0, "KEV": 1e3, "MEV": 1e6, "MILLI-EV": 1e-3, "MEV-": 1e6}
S_UNITS = {"B": 1.0, "MB": 1e-3, "B/SR": 1.0, "KB": 1e3}
NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")


def fetch(entry):
    e = str(entry).zfill(5)
    p = XC / f"{e}.x4"
    if not p.exists():
        with urllib.request.urlopen(URL.format(d=e[:3], e=e), timeout=60) as r:
            p.write_bytes(r.read())
    return p


_FORT = re.compile(r"^([-+]?\d*\.?\d*)([-+]\d+)$")


def x4float(s):
    """EXFOR numbers: '2.69   -09' and '1.23-4' mean 2.69e-9 and 1.23e-4 (Fortran style); also E/D forms."""
    s = (s or "").strip().replace(" ", "")
    if not s:
        return float("nan")
    m = _FORT.match(s)
    if m and m.group(1) not in ("", "+", "-", "."):
        return float(m.group(1)) * 10.0 ** int(m.group(2))
    try:
        return float(s.replace("D", "E").replace("d", "e"))
    except ValueError:
        return float("nan")


def _split_fields(line):
    """X4 fixed-width 11-character fields (6 per data line)."""
    line = line.rstrip("\n")
    return [line[i:i + 11].strip() for i in range(0, 66, 11)]


def _bib(lines):
    bib, key = {}, None
    for ln in lines:
        k = ln[:10].strip()
        if k in ("ENDBIB", "NOCOMMON", "COMMON", "DATA", "ENDSUBENT"):
            break
        if k:
            key = k
            bib[key] = ln[11:].rstrip()
        elif key:
            bib[key] += " " + ln[11:].strip()
    return bib


def _table(lines, start):
    """Parse a COMMON or DATA block starting at index `start` (the COMMON/DATA line)."""
    nf = int(lines[start][11:22])
    nrow = int(lines[start][22:33])
    hlines = (nf + 5) // 6
    heads, units = [], []
    for j in range(hlines):
        heads += _split_fields(lines[start + 1 + j])
        units += _split_fields(lines[start + 1 + hlines + j])
    heads, units = heads[:nf], units[:nf]
    rows, i = [], start + 1 + 2 * hlines
    while len(rows) < nrow:
        vals = []
        for j in range(hlines):
            vals += _split_fields(lines[i + j])
        i += hlines
        rows.append([x4float(v) for v in vals[:nf]])
    return heads, units, rows


def parse(entry):
    """All (N,TOT)-type subentries of an entry as a tidy DataFrame (E in eV, sigma/err in barn)."""
    lines = fetch(entry).read_text(errors="replace").splitlines()
    main_bib, out = {}, []
    sub_starts = [i for i, ln in enumerate(lines) if ln.startswith("SUBENT")]
    for si, s in enumerate(sub_starts):
        end = sub_starts[si + 1] if si + 1 < len(sub_starts) else len(lines)
        block = lines[s:end]
        subent = block[0].split()[1]
        bib = _bib(block[2:]) if len(block) > 2 and block[1].startswith("BIB") else {}
        if subent.endswith("001"):
            main_bib = bib
            continue
        common = {}
        for i, ln in enumerate(block):
            if ln.startswith("COMMON") and ln[11:22].strip().isdigit():
                h, u, r = _table(block, i)
                common = {hh: (vv, uu) for hh, uu, vv in zip(h, u, r[0])} if r else {}
            if ln.startswith("DATA ") and ln[11:22].strip().isdigit():
                h, u, rows = _table(block, i)
                reaction = bib.get("REACTION", main_bib.get("REACTION", ""))
                if "TOT" not in reaction.upper():
                    continue
                df = pd.DataFrame(rows, columns=[f"{a}|{b}" for a, b in zip(h, u)])
                ecol = next((c for c in df.columns if c.split("|")[0] in ("EN", "EN-CM", "E", "WVE-LN")), None)
                scol = next((c for c in df.columns if c.split("|")[0] == "DATA"), None)
                if ecol is None or scol is None:
                    continue
                dcol = next((c for c in df.columns if c.split("|")[0] in ("ERR-T", "DATA-ERR", "ERR-S")), None)
                eu, su = ecol.split("|")[1], scol.split("|")[1]
                ef, sf = E_UNITS.get(eu), S_UNITS.get(su)
                if ecol.startswith("WVE-LN") and eu in ("ANGSTROM", "A"):
                    ef = "wavelength"
                if ef is None or sf is None:
                    continue
                temp = None
                for key in ("TEMP", "TEMP-ERR"):
                    if key == "TEMP" and key in common:
                        temp = common[key][0]
                tcol = next((c for c in df.columns if c.split("|")[0] == "TEMP"), None)
                energy = 0.0818042 / df[ecol] ** 2 if ef == "wavelength" else df[ecol] * ef  # E[eV]=81.8042meV*A^2/lambda^2
                t = pd.DataFrame({"E_eV": energy, "sigma_b": df[scol] * sf,
                                  "err_b": (df[dcol] * sf if dcol and "PER-CENT" not in dcol
                                            else df[dcol] / 100 * df[scol] * sf if dcol else float("nan")),
                                  "temp_K": df[tcol] if tcol else temp})
                t["entry"], t["subentry"], t["reaction"] = str(entry).zfill(5), subent, reaction
                t["title"] = main_bib.get("TITLE", "")
                t["reference"] = main_bib.get("REFERENCE", "")
                t["author"] = main_bib.get("AUTHOR", "")
                t["sample"] = (main_bib.get("SAMPLE", "") + " " + bib.get("SAMPLE", "")).strip()
                t["status"] = bib.get("STATUS", main_bib.get("STATUS", ""))
                t["sub_bib"] = " | ".join(f"{k}: {v}" for k, v in bib.items())[:500]
                out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()
