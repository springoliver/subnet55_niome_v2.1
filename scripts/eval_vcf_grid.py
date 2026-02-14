#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

GENOTYPE_HET = {"0/1", "1/0", "0|1", "1|0"}
GENOTYPE_HOM_ALT = {"1/1", "1|1"}


@dataclass(frozen=True)
class Panel:
    rows: Tuple[Tuple[int, str, str, str], ...]
    source: str = "unknown"

    @property
    def n_sites(self) -> int:
        return len(self.rows)

    @property
    def keys(self) -> set:
        return {(pos, ref, alt) for pos, ref, alt, _ in self.rows}

    @property
    def positions(self) -> set:
        return {pos for pos, _, _, _ in self.rows}


def compare_panels(truth: Panel, miner: Panel) -> Dict[str, object]:
    t_keys = truth.keys
    m_keys = miner.keys
    t_pos = truth.positions
    m_pos = miner.positions
    tp_set = t_keys & m_keys
    n_t = truth.n_sites
    n_m = miner.n_sites
    prec = len(tp_set) / n_m if n_m else 0.0
    rec = len(tp_set) / n_t if n_t else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    t_gt = {(pos, ref, alt): gt for pos, ref, alt, gt in truth.rows}
    m_gt = {(pos, ref, alt): gt for pos, ref, alt, gt in miner.rows}
    gt_match = 0
    gt_total = 0
    for key in tp_set:
        gt_total += 1
        if t_gt.get(key) == m_gt.get(key):
            gt_match += 1
    return {
        "truth_n": n_t,
        "miner_n": n_m,
        "tp": len(tp_set),
        "fp": len(m_keys - t_keys),
        "fn": len(t_keys - m_keys),
        "pos_tp": len(t_pos & m_pos),
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "gt_match": gt_match,
        "gt_total": gt_total,
        "jaccard": len(tp_set) / len(t_keys | m_keys) if (t_keys | m_keys) else 1.0,
    }


def load_truth_vcf(path: Path) -> Panel:
    if not path.is_file():
        return Panel((), source="empty")
    text = path.read_text(encoding="utf-8", errors="replace")
    variants = parse_vcf_text(text)
    return panel_from_variants(variants, source="truth")


@dataclass(frozen=True)
class ParsedVariant:
    pos: int
    ref: str
    alt: str
    gt: str
    qual: float
    info: Dict[str, str]


def normalize_variant(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    if alt == "." or ref == ".":
        return pos, ref, alt
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return pos, ref, alt


def parse_info(info_text: str) -> Dict[str, str]:
    items = [token for token in info_text.split(";") if token]
    info: Dict[str, str] = {}
    for item in items:
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = ""
    return info


def extract_af(info: Dict[str, str]) -> Optional[float]:
    for key in ("AF", "AF_ESP", "AF1", "AF_INDEL"):
        if key in info:
            try:
                return float(info[key])
            except ValueError:
                return None
    return None


def parse_vcf_text(text: str) -> List[ParsedVariant]:
    variants: List[ParsedVariant] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5 or not parts[0].startswith("chr7"):
            continue
        pos = int(parts[1])
        ref = parts[3]
        alt = parts[4]
        if alt == ".":
            continue
        qual = 0.0
        if len(parts) > 5 and parts[5] not in {".", ""}:
            try:
                qual = float(parts[5])
            except ValueError:
                qual = 0.0
        info_text = parts[7] if len(parts) > 7 else ""
        info = parse_info(info_text)
        gt = "."
        if len(parts) > 9:
            sample = parts[9]
            format_keys = parts[8].split(":")
            sample_values = sample.split(":")
            sample_dict = {k: v for k, v in zip(format_keys, sample_values) if k}
            gt = sample_dict.get("GT", gt)
            if "DP" in sample_dict and sample_dict["DP"]:
                info.setdefault("DP", sample_dict["DP"])
        pos, ref, alt = normalize_variant(pos, ref, alt)
        variants.append(ParsedVariant(pos=pos, ref=ref, alt=alt, gt=gt, qual=qual, info=info))
    return variants


def panel_from_variants(variants: Iterable[ParsedVariant], source: str = "miner") -> Panel:
    rows = tuple((v.pos, v.ref, v.alt, v.gt) for v in variants)
    return Panel(rows, source=source)


def filter_variants(
    variants: Iterable[ParsedVariant],
    min_af: float,
    min_dp: int,
    min_qual: float,
    require_gt: str,
) -> List[ParsedVariant]:
    result: List[ParsedVariant] = []
    for v in variants:
        af = extract_af(v.info)
        if min_af > 0.0:
            if af is None or af < min_af:
                continue
        dp = 0
        if "DP" in v.info and v.info["DP"]:
            try:
                dp = int(v.info["DP"])
            except ValueError:
                dp = 0
        if dp < min_dp:
            continue
        if v.qual < min_qual:
            continue
        if require_gt == "het" and v.gt not in GENOTYPE_HET:
            continue
        if require_gt == "hom-alt" and v.gt not in GENOTYPE_HOM_ALT:
            continue
        result.append(v)
    return result


def choose_miner_entry(items: List[dict], miner_uid: Optional[int], best_by: str) -> dict:
    rows = [item for item in items if isinstance(item.get("log"), str) and "Miner VCF" in item["log"]]
    if not rows:
        raise ValueError("No miner entries with Miner VCF found in the JSON file.")
    if miner_uid is not None:
        for item in rows:
            if item.get("miner_uid") == miner_uid:
                return item
        raise ValueError(f"Miner UID {miner_uid} not found in miner.json entries.")
    if best_by in {"f1", "f1_score"}:
        key = "f1_score"
    else:
        key = best_by
    return max(rows, key=lambda item: float(item.get(key, 0.0)))


def load_miner_variants_from_json(json_path: Path, miner_uid: Optional[int], best_by: str) -> Tuple[List[ParsedVariant], dict]:
    data = json.loads(json_path.read_text(encoding="utf-8"))
    items = data.get("items") if isinstance(data, dict) else []
    if not isinstance(items, list):
        raise ValueError("miner.json does not contain an items list.")
    entry = choose_miner_entry(items, miner_uid, best_by)
    log = entry["log"]
    _, _, vcf_text = log.partition("Miner VCF\n")
    variants = parse_vcf_text(vcf_text)
    if not variants:
        raise ValueError("Selected miner entry contains no parsed variants.")
    return variants, entry


def format_float(value: float) -> str:
    return f"{value:.4f}"


def print_top_results(results: List[dict], title: str, top_n: int = 10) -> None:
    print(f"\n{title}")
    print("min_af\tmin_dp\tmin_qual\trequire_gt\ttp\tfp\tfn\tprec\trec\tf1\tcount\tgt_match/gt_total")
    for row in sorted(results, key=lambda m: m["f1"], reverse=True)[:top_n]:
        print(
            "\t".join(
                [
                    format_float(row["min_af"]),
                    str(row["min_dp"]),
                    format_float(row["min_qual"]),
                    row["require_gt"],
                    str(row["tp"]),
                    str(row["fp"]),
                    str(row["fn"]),
                    format_float(row["precision"]),
                    format_float(row["recall"]),
                    format_float(row["f1"]),
                    str(row["miner_n"]),
                    f"{row['gt_match']}/{row['gt_total']}",
                ]
            )
        )


def save_csv(results: List[dict], csv_path: Path) -> None:
    fieldnames = [
        "min_af",
        "min_dp",
        "min_qual",
        "require_gt",
        "miner_n",
        "tp",
        "fp",
        "fn",
        "precision",
        "recall",
        "f1",
        "gt_match",
        "gt_total",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow({k: row[k] for k in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate AF/DP/QUAL/GT filter grids against same-round truth VCF."
    )
    parser.add_argument("--round", required=True, help="Round folder under Results, e.g. 5.21.01")
    parser.add_argument("--miner-json", default=None, help="Path to miner.json for the chosen round.")
    parser.add_argument("--miner-uid", type=int, default=None, help="Miner UID to evaluate.")
    parser.add_argument(
        "--best-by",
        choices=("f1", "precision", "recall", "f1_score"),
        default="f1",
        help="Choose the miner entry by best metric when miner_uid is not provided.",
    )
    parser.add_argument(
        "--truth-file",
        default=None,
        help="Explicit truth VCF path. Defaults to Results/<round>/real_correct_result/truth.vcf.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="Optional CSV path to record every grid row.",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=10,
        help="How many top rows to print for each summary.",
    )
    args = parser.parse_args()

    round_path = Path("Results") / args.round
    miner_json_path = Path(args.miner_json) if args.miner_json else round_path / "miner.json"
    truth_path = Path(args.truth_file) if args.truth_file else round_path / "real_correct_result" / "truth.vcf"

    if not miner_json_path.is_file():
        print(f"ERROR: miner.json not found at {miner_json_path}", file=sys.stderr)
        return 1
    if not truth_path.is_file():
        print(f"ERROR: truth VCF not found at {truth_path}", file=sys.stderr)
        return 1

    truth_panel = load_truth_vcf(truth_path)
    if truth_panel.n_sites == 0:
        print(f"ERROR: truth VCF at {truth_path} parsed zero variants.", file=sys.stderr)
        return 1

    miner_variants, miner_entry = load_miner_variants_from_json(miner_json_path, args.miner_uid, args.best_by)
    baseline_panel = panel_from_variants(miner_variants, source="miner")
    baseline = compare_panels(truth_panel, baseline_panel)
    print("Round:", args.round)
    print("Truth:", truth_path)
    print(
        "Selected miner entry:",
        f"miner_uid={miner_entry.get('miner_uid')} f1_score={miner_entry.get('f1_score')} precision={miner_entry.get('precision')} recall={miner_entry.get('recall')}",
    )
    print(
        "Baseline: n_miner=",
        baseline["miner_n"],
        "tp=",
        baseline["tp"],
        "fp=",
        baseline["fp"],
        "fn=",
        baseline["fn"],
        "prec=",
        format_float(baseline["precision"]),
        "rec=",
        format_float(baseline["recall"]),
        "f1=",
        format_float(baseline["f1"]),
        "gt_match=",
        f"{baseline['gt_match']}/{baseline['gt_total']}",
    )

    af_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25]
    dp_values = [0, 10, 15, 20, 25, 30]
    qual_values = [0.0, 20.0, 40.0, 60.0, 80.0]
    require_gt_values = ["any", "het", "hom-alt"]

    results: List[dict] = []
    for min_af in af_values:
        for min_dp in dp_values:
            for min_qual in qual_values:
                for require_gt in require_gt_values:
                    filtered = filter_variants(
                        miner_variants,
                        min_af=min_af,
                        min_dp=min_dp,
                        min_qual=min_qual,
                        require_gt=require_gt,
                    )
                    filtered_panel = panel_from_variants(filtered, source="miner")
                    metrics = compare_panels(truth_panel, filtered_panel)
                    metrics.update(
                        {
                            "min_af": min_af,
                            "min_dp": min_dp,
                            "min_qual": min_qual,
                            "require_gt": require_gt,
                            "miner_n": filtered_panel.n_sites,
                        }
                    )
                    results.append(metrics)

    if args.output_csv:
        save_csv(results, Path(args.output_csv))
        print(f"Saved grid results to {args.output_csv}")

    print_top_results(sorted(results, key=lambda r: r["f1"], reverse=True), "Top combos by F1", top_n=args.max_results)

    thresholds = [0.75, 0.70]
    for threshold in thresholds:
        filtered_results = [r for r in results if r["recall"] >= threshold]
        if not filtered_results:
            continue
        top_precision = max(filtered_results, key=lambda r: r["precision"])
        print(f"\nBest precision with recall >= {threshold:.2f}:")
        print(
            "min_af\tmin_dp\tmin_qual\trequire_gt\ttp\tfp\tfn\tprec\trec\tf1\tcount\tgt_match/gt_total"
        )
        print(
            "\t".join(
                [
                    format_float(top_precision["min_af"]),
                    str(top_precision["min_dp"]),
                    format_float(top_precision["min_qual"]),
                    top_precision["require_gt"],
                    str(top_precision["tp"]),
                    str(top_precision["fp"]),
                    str(top_precision["fn"]),
                    format_float(top_precision["precision"]),
                    format_float(top_precision["recall"]),
                    format_float(top_precision["f1"]),
                    str(top_precision["miner_n"]),
                    f"{top_precision['gt_match']}/{top_precision['gt_total']}",
                ]
            )
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
