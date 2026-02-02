"""Parse and compare VCF panels from validator logs or truth files."""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

VariantKey = Tuple[int, str, str]
VariantRow = Tuple[int, str, str, str]  # pos, ref, alt, gt


@dataclass(frozen=True)
class Panel:
    rows: Tuple[VariantRow, ...]
    source: str = "unknown"  # oracle | native | other | truth | empty

    @property
    def n_sites(self) -> int:
        return len(self.rows)

    @property
    def keys(self) -> Set[VariantKey]:
        return {(p, r, a) for p, r, a, _ in self.rows}

    @property
    def positions(self) -> Set[int]:
        return {p for p, _, _, _ in self.rows}

    def fingerprint(self) -> str:
        key = "|".join(f"{p}:{r}>{a}:{g}" for p, r, a, g in sorted(self.rows))
        return hashlib.md5(key.encode()).hexdigest()[:16]

    def classify_log(self, log: str) -> str:
        if "CLNDN" in log or "ONCDN" in log:
            return "oracle"
        if "niome-native" in log or "niome_miner" in log:
            return "native"
        return "other"


def parse_vcf_text(text: str, source: str = "unknown") -> Panel:
    rows: List[VariantRow] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5 or not parts[0].startswith("chr7"):
            continue
        pos = int(parts[1])
        ref, alt = parts[3], parts[4]
        gt = "."
        if len(parts) > 9:
            sample = parts[9]
            gt = sample.split(":")[0] if ":" in sample else sample
        rows.append((pos, ref, alt, gt))
    return Panel(tuple(rows), source=source)


def parse_log_vcf(log: str) -> Panel:
    if "Miner VCF\n" not in log:
        return Panel((), source="empty")
    body = log.split("Miner VCF\n", 1)[1]
    panel = parse_vcf_text(body, source="unknown")
    return Panel(panel.rows, source=panel.classify_log(log) if panel.rows else "empty")


def load_truth_vcf(path: Path) -> Panel:
    if not path.is_file():
        return Panel((), source="empty")
    return parse_vcf_text(
        path.read_text(encoding="utf-8", errors="replace"), source="truth"
    )


def band_from_count(n: int) -> str:
    if n <= 14:
        return "low"
    if n <= 22:
        return "mid"
    if n <= 29:
        return "high"
    return "ultra"


def jaccard(a: Set[VariantKey], b: Set[VariantKey]) -> float:
    if not a and not b:
        return 1.0
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def compare_panels(truth: Panel, miner: Panel) -> Dict:
    t_keys, m_keys = truth.keys, miner.keys
    t_pos, m_pos = truth.positions, miner.positions
    tp = t_keys & m_keys
    n_t, n_m = truth.n_sites, miner.n_sites
    count_penalty = min(n_m, n_t) / max(n_m, n_t) if n_t and n_m else 0.0
    prec = len(tp) / n_m if n_m else 0.0
    rec = len(tp) / n_t if n_t else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0

    gt_match = 0
    gt_total = 0
    t_gt = {(p, r, a): g for p, r, a, g in truth.rows}
    m_gt = {(p, r, a): g for p, r, a, g in miner.rows}
    for k in tp:
        gt_total += 1
        if t_gt[k] == m_gt[k]:
            gt_match += 1

    return {
        "truth_n": n_t,
        "miner_n": n_m,
        "tp": len(tp),
        "fp": len(m_keys - t_keys),
        "fn": len(t_keys - m_keys),
        "pos_tp": len(t_pos & m_pos),
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "count_penalty": count_penalty,
        "gt_match": gt_match,
        "gt_total": gt_total,
        "missed_positions": sorted(t_pos - m_pos),
        "extra_positions": sorted(m_pos - t_pos),
        "jaccard": jaccard(t_keys, m_keys),
    }


def read_key_from_urls(read1: str, read2: str) -> str:
    parts = []
    for url in (read1 or "", read2 or ""):
        base = url.split("?")[0] if url else ""
        parts.append(base.rsplit("/", 1)[-1] if base else "")
    raw = "|".join(parts)
    return hashlib.md5(raw.encode()).hexdigest()[:12]


def region_length(region: str) -> int:
    _, rest = region.split(":")
    start, end = rest.split("-")
    return int(end) - int(start)


def extract_source_rev(log: str) -> Optional[str]:
    m = re.search(r"##source=niome_miner[_-](\S+)", log)
    return m.group(1) if m else None
