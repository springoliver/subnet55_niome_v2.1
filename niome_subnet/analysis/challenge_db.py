"""
Build and query a cumulative NIOME challenge database from all Results/.

Output: Results/niome_challenge_db/
  manifest.json       — build metadata + round index
  rounds/<key>.json   — per-round summary (no full validator logs)
  training/
    band_stats.json
    truth_position_hotspots.json
    strategy_calibration.json
    fleet_timeline.jsonl
    top_miner_patterns.json

Rebuild after adding new round folders:
  python scripts/build_challenge_db.py
  python scripts/build_challenge_db.py --rounds 5.24.01 5.24.02
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from .vcf_panel import (
    Panel,
    band_from_count,
    compare_panels,
    extract_source_rev,
    jaccard,
    load_truth_vcf,
    parse_log_vcf,
    read_key_from_urls,
    region_length,
)

DEFAULT_RESULTS_ROOT = Path(__file__).resolve().parents[2] / "Results"
DB_DIR_NAME = "niome_challenge_db"
VALIDATORS = {119, 154, 58}
DEFAULT_FLEET_UIDS = {
    141, 139, 71, 50, 99, 209, 235, 226, 124, 9, 97, 217, 36, 92, 225, 38, 155,
}
FLEET_STRATEGY = {
    50: "high_recall",
    99: "high_recall",
    209: "high_recall",
    235: "high_recall",
    226: "win",
    124: "v5_style",
    9: "v10",
    97: "v5_style",
    217: "v10",
    141: "auto",
    139: "auto",  # legacy UID if still registered on-chain
    71: "v5_style",
}


@dataclass
class RoundPaths:
    key: str
    dir: Path
    rel: str


def discover_round_dirs(results_root: Path) -> List[RoundPaths]:
    """All folders with miner.json + task.json or problem.json."""
    found: List[RoundPaths] = []
    seen_task_ids: Set[str] = set()

    def scan(base: Path, prefix: str = "") -> None:
        if not base.is_dir():
            return
        for child in sorted(base.iterdir()):
            if not child.is_dir():
                continue
            rel = f"{prefix}/{child.name}" if prefix else child.name
            has_miner = (child / "miner.json").is_file()
            has_task = (child / "task.json").is_file() or (child / "problem.json").is_file()
            if has_miner and has_task:
                key = rel.replace("/", "__")
                found.append(RoundPaths(key=key, dir=child, rel=rel))
                continue
            # recurse one level for old_version_result only
            if child.name == "old_version_result" or prefix == "":
                scan(child, rel)

    scan(results_root)
    return found


def discover_truth_archives(results_root: Path) -> Dict[str, Path]:
    """task_id -> truth.vcf path (correct_answer_* and real_correct_result)."""
    out: Dict[str, Path] = {}
    for task_json in results_root.rglob("task.json"):
        try:
            data = json.loads(task_json.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        tid = data.get("task_id")
        if not tid:
            continue
        parent = task_json.parent
        for sub in ("real_correct_result", "."):
            vcf = parent / sub / "truth.vcf" if sub != "." else parent / "truth.vcf"
            if vcf.is_file():
                out[tid] = vcf
                break
    for task_json in results_root.rglob("problem.json"):
        try:
            data = json.loads(task_json.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        tid = data.get("task_id") or data.get("id")
        if not tid:
            continue
        vcf = task_json.parent / "truth.vcf"
        if vcf.is_file():
            out[tid] = vcf
    return out


def load_task(round_dir: Path) -> Dict[str, Any]:
    for name in ("task.json", "problem.json"):
        p = round_dir / name
        if p.is_file():
            return json.loads(p.read_text(encoding="utf-8"))
    return {}


def load_miner_items(round_dir: Path) -> List[Dict]:
    p = round_dir / "miner.json"
    if not p.is_file():
        return []
    data = json.loads(p.read_text(encoding="utf-8", errors="replace"))
    return data if isinstance(data, list) else data.get("items", [])


def best_entry_per_uid(
    items: List[Dict], task_id: str, validators: Set[int]
) -> Dict[int, Dict]:
    pool = [
        x
        for x in items
        if x.get("task_id") == task_id and x.get("validator_uid") in validators
    ]
    by_uid: Dict[int, Dict] = {}
    for x in pool:
        uid = x["miner_uid"]
        sc = float(x.get("final_score") or 0)
        if uid not in by_uid or sc > float(by_uid[uid].get("final_score") or 0):
            by_uid[uid] = x
    return by_uid


def summarize_uid_entry(entry: Dict, truth: Optional[Panel]) -> Dict[str, Any]:
    log = entry.get("log", "")
    panel = parse_log_vcf(log)
    out: Dict[str, Any] = {
        "miner_uid": entry["miner_uid"],
        "validator_uid": entry.get("validator_uid"),
        "final_score": float(entry.get("final_score") or 0),
        "vcf_score": float(entry.get("vcf_score") or 0),
        "annotation_score": float(entry.get("annotation_score") or 0),
        "n_sites": panel.n_sites,
        "panel_class": panel.source,
        "fingerprint": panel.fingerprint() if panel.rows else "",
        "source_rev": extract_source_rev(log),
    }
    if truth and truth.rows:
        cmp = compare_panels(truth, panel)
        out["vs_truth"] = {
            k: cmp[k]
            for k in (
                "tp",
                "fp",
                "fn",
                "precision",
                "recall",
                "f1",
                "count_penalty",
                "gt_match",
                "gt_total",
                "jaccard",
            )
        }
        out["vs_truth"]["missed_positions"] = cmp["missed_positions"][:40]
        out["vs_truth"]["extra_positions"] = cmp["extra_positions"][:40]
    return out


def ingest_round(
    rp: RoundPaths,
    truth_by_task: Dict[str, Path],
    fleet_uids: Set[int],
) -> Optional[Dict[str, Any]]:
    task = load_task(rp.dir)
    task_id = task.get("task_id") or task.get("id")
    if not task_id:
        return None

    gc = task.get("genome_context") or {}
    region = gc.get("region", "")
    rlen = region_length(region) if region and ":" in region else 0
    inp = task.get("input") or {}
    r1 = inp.get("read1_fastq") or inp.get("reads_1") or ""
    r2 = inp.get("read2_fastq") or inp.get("reads_2") or ""
    rk = read_key_from_urls(r1, r2)

    truth_path = rp.dir / "real_correct_result" / "truth.vcf"
    if not truth_path.is_file():
        truth_path = truth_by_task.get(task_id)
    truth_panel = load_truth_vcf(truth_path) if truth_path else Panel((), "empty")
    has_truth = bool(truth_panel.rows)

    items = load_miner_items(rp.dir)
    by_uid = best_entry_per_uid(items, task_id, VALIDATORS)
    if not by_uid:
        return None

    top_uid = max(by_uid, key=lambda u: float(by_uid[u].get("final_score") or 0))
    top = summarize_uid_entry(by_uid[top_uid], truth_panel if has_truth else None)

    oracle_best = None
    native_best = None
    for uid, entry in by_uid.items():
        s = summarize_uid_entry(entry, truth_panel if has_truth else None)
        if s["panel_class"] == "oracle":
            if oracle_best is None or s["final_score"] > oracle_best["final_score"]:
                oracle_best = s
        if s["panel_class"] == "native":
            if native_best is None or s["final_score"] > native_best["final_score"]:
                native_best = s

    ref_panel = truth_panel if has_truth else Panel(
        tuple(
            (p, r, a, g)
            for p, r, a, g in parse_log_vcf(
                by_uid[top_uid].get("log", "")
            ).rows
        ),
        top.get("panel_class", "other"),
    )
    band = band_from_count(ref_panel.n_sites if ref_panel.rows else top["n_sites"])

    fleet: Dict[str, Any] = {}
    fp_groups: Dict[str, List[int]] = defaultdict(list)
    for uid in fleet_uids:
        if uid not in by_uid:
            continue
        s = summarize_uid_entry(by_uid[uid], truth_panel if has_truth else None)
        s["configured_strategy"] = FLEET_STRATEGY.get(uid, "?")
        fleet[str(uid)] = s
        if s["fingerprint"]:
            fp_groups[s["fingerprint"]].append(uid)

    duplicate_panels = [
        {"fingerprint": fp, "uids": uids, "n_sites": fleet[str(uids[0])]["n_sites"]}
        for fp, uids in fp_groups.items()
        if len(uids) > 1
    ]

    return {
        "round_key": rp.key,
        "round_path": rp.rel,
        "task_id": task_id,
        "region": region,
        "region_len": rlen,
        "read_key": rk,
        "band": band,
        "has_truth": has_truth,
        "truth_path": str(truth_path) if has_truth and truth_path else None,
        "truth_n": truth_panel.n_sites,
        "top_miner": top,
        "oracle_best": oracle_best,
        "native_best": native_best,
        "fleet": fleet,
        "fleet_duplicate_panels": duplicate_panels,
        "n_miners_scored": len(by_uid),
        "recommended_strategy": _strategy_for_band(band),
    }


def _strategy_for_band(b: str) -> str:
    """Match task_strategy._BAND_TO_STRATEGY (native winners, not recall merge)."""
    return {
        "low": "v5_style",
        "mid": "v10",
        "high": "v5_style",
        "ultra": "v5_style",
    }.get(b, "v10")


def build_training_aggregates(rounds: List[Dict[str, Any]]) -> Dict[str, Any]:
    band_truth_counts: Dict[str, List[int]] = defaultdict(list)
    band_top_counts: Dict[str, List[int]] = defaultdict(list)
    position_freq: Counter = Counter()
    missed_freq: Counter = Counter()
    fleet_scores: Dict[str, List[float]] = defaultdict(list)
    native_gaps: List[Dict[str, Any]] = []

    for rd in rounds:
        band = rd["band"]
        if rd["has_truth"]:
            band_truth_counts[band].append(rd["truth_n"])
            truth_path = rd.get("truth_path")
            if truth_path:
                tp = load_truth_vcf(Path(truth_path))
                for pos in tp.positions:
                    position_freq[pos] += 1
        band_top_counts[band].append(rd["top_miner"]["n_sites"])

        for uid, rec in rd.get("fleet", {}).items():
            strat = rec.get("configured_strategy", "?")
            fleet_scores[strat].append(rec["final_score"])
            vs = rec.get("vs_truth")
            if vs:
                for pos in vs.get("missed_positions", []):
                    missed_freq[pos] += 1

        nb = rd.get("native_best")
        if nb:
            best_fleet = max(
                (float(rec["final_score"]) for rec in rd.get("fleet", {}).values()),
                default=0.0,
            )
            native_gaps.append(
                {
                    "round": rd["round_key"],
                    "band": band,
                    "native_score": nb["final_score"],
                    "fleet_best": best_fleet,
                    "gap": nb["final_score"] - best_fleet,
                }
            )

    def _stats(vals: List[int]) -> Dict[str, Any]:
        if not vals:
            return {}
        s = sorted(vals)
        return {
            "n": len(s),
            "min": s[0],
            "max": s[-1],
            "median": s[len(s) // 2],
            "mean": round(sum(s) / len(s), 2),
        }

    calibration = {}
    for b in ("low", "mid", "high", "ultra"):
        calibration[b] = {
            "strategy": _strategy_for_band(b),
            "truth_site_counts": _stats(band_truth_counts[b]),
            "top_site_counts": _stats(band_top_counts[b]),
            "curriculum_target_suggested": _suggest_curriculum(b, band_truth_counts[b], band_top_counts[b]),
        }

    top_patterns = []
    for rd in sorted(rounds, key=lambda x: x["round_key"]):
        top_patterns.append(
            {
                "round": rd["round_key"],
                "task_id": rd["task_id"][:8],
                "band": rd["band"],
                "top_uid": rd["top_miner"]["miner_uid"],
                "top_n": rd["top_miner"]["n_sites"],
                "top_class": rd["top_miner"]["panel_class"],
                "top_score": rd["top_miner"]["final_score"],
                "truth_n": rd.get("truth_n"),
            }
        )

    return {
        "band_stats": {
            b: calibration[b]
            for b in ("low", "mid", "high", "ultra")
        },
        "strategy_calibration": calibration,
        "truth_position_hotspots": [
            {"pos": p, "count": c}
            for p, c in position_freq.most_common(80)
        ],
        "commonly_missed_positions": [
            {"pos": p, "count": c}
            for p, c in missed_freq.most_common(60)
        ],
        "fleet_score_by_strategy": {
            k: {
                "n": len(v),
                "mean": round(sum(v) / len(v), 4) if v else 0,
                "max": round(max(v), 4) if v else 0,
            }
            for k, v in fleet_scores.items()
        },
        "native_gap_rounds": native_gaps,
        "top_miner_patterns": top_patterns,
    }


def _suggest_curriculum(
    band: str,
    truth_counts: List[int],
    top_counts: List[int],
) -> int:
    pool = truth_counts or top_counts
    if not pool:
        return {"low": 14, "mid": 20, "high": 28, "ultra": 32}.get(band, 30)
    s = sorted(pool)
    median = s[len(s) // 2]
    return max(10, min(36, median))


def build_database(
    results_root: Optional[Path] = None,
    fleet_uids: Optional[Set[int]] = None,
    only_rounds: Optional[Iterable[str]] = None,
) -> Path:
    results_root = Path(results_root or DEFAULT_RESULTS_ROOT)
    db_dir = results_root / DB_DIR_NAME
    rounds_dir = db_dir / "rounds"
    training_dir = db_dir / "training"
    rounds_dir.mkdir(parents=True, exist_ok=True)
    training_dir.mkdir(parents=True, exist_ok=True)

    fleet_uids = fleet_uids or DEFAULT_FLEET_UIDS
    truth_by_task = discover_truth_archives(results_root)
    all_round_paths = discover_round_dirs(results_root)

    if only_rounds:
        want = {r.replace("/", "__") for r in only_rounds}
        all_round_paths = [rp for rp in all_round_paths if rp.key in want or rp.rel in only_rounds]

    ingested: List[Dict[str, Any]] = []
    for rp in all_round_paths:
        rec = ingest_round(rp, truth_by_task, fleet_uids)
        if not rec:
            continue
        out_path = rounds_dir / f"{rp.key}.json"
        out_path.write_text(json.dumps(rec, indent=2), encoding="utf-8")
        ingested.append(rec)

    # merge with existing round json if incremental partial build
    existing_keys = {r["round_key"] for r in ingested}
    for p in rounds_dir.glob("*.json"):
        if p.stem not in existing_keys:
            try:
                ingested.append(json.loads(p.read_text(encoding="utf-8")))
            except json.JSONDecodeError:
                pass

    ingested.sort(key=lambda x: x["round_key"])
    training = build_training_aggregates(ingested)
    (training_dir / "band_stats.json").write_text(
        json.dumps(training["band_stats"], indent=2), encoding="utf-8"
    )
    (training_dir / "strategy_calibration.json").write_text(
        json.dumps(training["strategy_calibration"], indent=2), encoding="utf-8"
    )
    (training_dir / "truth_position_hotspots.json").write_text(
        json.dumps(training["truth_position_hotspots"], indent=2), encoding="utf-8"
    )
    (training_dir / "commonly_missed_positions.json").write_text(
        json.dumps(training["commonly_missed_positions"], indent=2), encoding="utf-8"
    )
    (training_dir / "fleet_score_by_strategy.json").write_text(
        json.dumps(training["fleet_score_by_strategy"], indent=2), encoding="utf-8"
    )
    (training_dir / "top_miner_patterns.json").write_text(
        json.dumps(training["top_miner_patterns"], indent=2), encoding="utf-8"
    )
    with (training_dir / "fleet_timeline.jsonl").open("w", encoding="utf-8") as fh:
        for rd in ingested:
            for uid, rec in rd.get("fleet", {}).items():
                line = {
                    "round": rd["round_key"],
                    "task_id": rd["task_id"][:8],
                    "band": rd["band"],
                    "uid": int(uid),
                    "strategy": rec.get("configured_strategy"),
                    "score": rec["final_score"],
                    "n_sites": rec["n_sites"],
                    "fingerprint": rec.get("fingerprint"),
                    "has_truth": rd["has_truth"],
                }
                if rec.get("vs_truth"):
                    line["jaccard_truth"] = rec["vs_truth"].get("jaccard")
                fh.write(json.dumps(line) + "\n")

    manifest = {
        "version": 1,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "results_root": str(results_root),
        "n_rounds": len(ingested),
        "n_with_truth": sum(1 for r in ingested if r["has_truth"]),
        "truth_tasks_indexed": len(truth_by_task),
        "rounds": [
            {
                "key": r["round_key"],
                "path": r["round_path"],
                "task_id": r["task_id"],
                "band": r["band"],
                "has_truth": r["has_truth"],
                "truth_n": r.get("truth_n"),
                "top_n": r["top_miner"]["n_sites"],
            }
            for r in ingested
        ],
    }
    (db_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    return db_dir


def load_manifest(results_root: Optional[Path] = None) -> Dict[str, Any]:
    root = Path(results_root or DEFAULT_RESULTS_ROOT)
    path = root / DB_DIR_NAME / "manifest.json"
    if not path.is_file():
        raise FileNotFoundError(
            f"Challenge DB not built. Run: python scripts/build_challenge_db.py\n"
            f"Expected: {path}"
        )
    return json.loads(path.read_text(encoding="utf-8"))


def load_round(key: str, results_root: Optional[Path] = None) -> Dict[str, Any]:
    root = Path(results_root or DEFAULT_RESULTS_ROOT)
    path = root / DB_DIR_NAME / "rounds" / f"{key}.json"
    return json.loads(path.read_text(encoding="utf-8"))
