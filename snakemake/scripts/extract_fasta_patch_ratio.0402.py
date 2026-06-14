#!/usr/bin/env python3

import argparse
import io
import gzip
import os
import re
from typing import Dict, Iterator, List, Tuple, Optional, Literal, TextIO


ComponentKind = Literal["impuT2T", "segment"]

IMPUT2T_REGEX = re.compile(r".*_impuT2T_\d+$")

# Temporary: print first/last endpoint flank ratio and pass/fail (set False to silence)
DEBUG_ENDPOINT_FLANK_CHECKS = True


def open_maybe_gzip(path: str) -> TextIO:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, mode="rb"))
    return open(path, "r")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter AGP patch blocks by segment ratio, optionally disassemble failing patches "
            "by removing the largest 'segment' component, and emit a filtered AGP and a "
            "corresponding FASTA built from the input contig FASTA."
        )
    )
    parser.add_argument("-i", "--agp", required=True, help="Input AGP file (.agp or .agp.gz)")
    parser.add_argument("-f", "--fasta", required=True, help="Input contig FASTA (bgzip/plain)")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix for filtered files")
    parser.add_argument(
        "--ratio-threshold",
        type=float,
        default=0.10,
        help="Maximum allowed segment ratio in a (sub)patch [default: 0.10]",
    )
    parser.add_argument(
        "--segment-flank-ratio",
        type=float,
        default=1,
        help="Max allowed ratio of (largest) segment length to shortest flank; exceed triggers disassembly [default: 1]",
    )
    parser.add_argument(
        "--endpoint-segment-flank-ratio",
        type=float,
        default=3.0,
        help=(
            "Relaxed max ratio of first (resp. last) segment length to total W-span before (resp. after) "
            "that segment (all consecutive patch components, e.g. multiple impuT2T); exceed triggers "
            "disassembly unless edge-protected [default: 2]"
        ),
    )
    parser.add_argument(
        "--disable-endpoint-segment-flank-filter",
        action="store_true",
        help=(
            "Disable first/last segment vs prefix/suffix W-span endpoint filtering "
            "(legacy behavior)."
        ),
    )
    parser.add_argument(
        "--min-contig-frac",
        type=float,
        default=0.10,
        help="Keep individual contigs if length >= frac * longest kept patch length [default: 0.10]",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write header/comment lines in filtered AGP",
    )
    parser.add_argument(
        "--edge-info",
        type=str,
        default=None,
        help=(
            "Optional edge log (.edge.log). If provided, a filtered edge log "
            "will be written as <output_prefix>.filtered.edge.log where any "
            "edge with last-column True and probability < 0.5 is relabeled "
            "to False."
        ),
    )
    parser.add_argument(
        "--long-segment-bp",
        type=int,
        default=2_000_000,
        help=(
            "Any segment longer than this (bp) is removed (disassembled) unless the "
            "impuT2T-impuT2T edge across that segment has probability >= 0.5 in "
            "--edge-info; without an edge log, long segments are not protected "
            "[default: 2000000]"
        ),
    )
    return parser.parse_args()


def classify_component(component_id: str) -> ComponentKind:
    """impuT2T if ID matches *_impuT2T_<digits>; otherwise treat as segment fill (incl. *_segment_* and plain donor IDs)."""
    if IMPUT2T_REGEX.match(component_id):
        return "impuT2T"
    return "segment"


class AgpEntry:
    def __init__(self, object_name: str, object_beg: int, object_end: int, part_num: int,
                 component_type: str, component_id: str, component_beg: Optional[int],
                 component_end: Optional[int], orientation: Optional[str]):
        self.object_name = object_name
        self.object_beg = object_beg
        self.object_end = object_end
        self.part_num = part_num
        self.component_type = component_type
        self.component_id = component_id
        self.component_beg = component_beg
        self.component_end = component_end
        self.orientation = orientation

    @property
    def length(self) -> int:
        return self.object_end - self.object_beg + 1


def parse_agp(agp_path: str) -> Dict[str, List[AgpEntry]]:
    per_object: Dict[str, List[AgpEntry]] = {}
    with open_maybe_gzip(agp_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            obj = fields[0]
            try:
                obj_beg = int(fields[1])
                obj_end = int(fields[2])
            except ValueError:
                continue
            try:
                part_num = int(fields[3])
            except ValueError:
                part_num = 0
            comp_type = fields[4]
            comp_id = fields[5]
            comp_beg: Optional[int] = None
            comp_end: Optional[int] = None
            orient: Optional[str] = None
            if comp_type == "W":
                try:
                    comp_beg = int(fields[6])
                    comp_end = int(fields[7])
                except ValueError:
                    continue
                orient = fields[8]
            entry = AgpEntry(obj, obj_beg, obj_end, part_num, comp_type, comp_id, comp_beg, comp_end, orient)
            per_object.setdefault(obj, []).append(entry)
    # Keep original order (by object_beg)
    for obj, entries in per_object.items():
        entries.sort(key=lambda e: (e.object_beg, e.part_num))
    return per_object


def compute_segment_ratio(entries: List[AgpEntry]) -> Tuple[int, int, float]:
    total_len = 0
    segment_len = 0
    for e in entries:
        if e.component_type != "W":
            continue
        total_len += e.length
        if classify_component(e.component_id) == "segment":
            segment_len += e.length
    ratio = (segment_len / total_len) if total_len > 0 else 0.0
    #print(f"RATIO: total_len: {total_len}, segment_len: {segment_len}, ratio: {ratio}")
    return total_len, segment_len, ratio


def find_largest_segment_index(entries: List[AgpEntry]) -> Optional[int]:
    max_len = -1
    max_idx: Optional[int] = None
    for i, e in enumerate(entries):
        if e.component_type != "W":
            continue
        if classify_component(e.component_id) == "segment":
            if e.length > max_len:
                max_len = e.length
                max_idx = i
    return max_idx


def segment_flank_ratio_exceeds(entries: List[AgpEntry], threshold: float) -> bool:
    """
    True if the largest segment's length / shortest_flank_length > threshold.
    The largest segment divides the patch into 'before' and 'after' fragments;
    shortest fragment = min(before_len, after_len). If no segment or shortest
    flank is 0, returns True (treat as exceeds / trigger disassembly).
    """
    idx = find_largest_segment_index(entries)
    if idx is None:
        return False
    segment_len = entries[idx].length
    before_len = sum(e.length for e in entries[:idx] if e.component_type == "W")
    after_len = sum(e.length for e in entries[idx + 1 :] if e.component_type == "W")
    shortest_flank = min(before_len, after_len)
    #print(f"FLANK: segment_len: {segment_len}, before_len: {before_len}, after_len: {after_len}, shortest_flank: {shortest_flank}, threshold: {threshold}, ratio: {segment_len / shortest_flank}, threshold: {threshold}, {segment_len / shortest_flank > threshold}")
    if shortest_flank <= 0:
        return True
    return (segment_len / shortest_flank) > threshold


def find_first_segment_index(entries: List[AgpEntry]) -> Optional[int]:
    for i, e in enumerate(entries):
        if e.component_type == "W" and classify_component(e.component_id) == "segment":
            return i
    return None


def find_last_segment_index(entries: List[AgpEntry]) -> Optional[int]:
    for i in range(len(entries) - 1, -1, -1):
        e = entries[i]
        if e.component_type == "W" and classify_component(e.component_id) == "segment":
            return i
    return None


def _w_span_before_index(entries: List[AgpEntry], idx: int) -> int:
    """Total bp of W lines strictly before index idx (patch prefix)."""
    return sum(e.length for e in entries[:idx] if e.component_type == "W")


def _w_span_after_index(entries: List[AgpEntry], idx: int) -> int:
    """Total bp of W lines strictly after index idx (patch suffix)."""
    return sum(e.length for e in entries[idx + 1 :] if e.component_type == "W")


def endpoint_segment_flank_exceeds(
    entries: List[AgpEntry],
    endpoint_threshold: float,
    edge_probs: Dict[Tuple[str, str], float],
    patch_label: str,
    prob_threshold: float = 0.5,
) -> bool:
    """
    First segment vs total W-span before it, last segment vs total W-span after it
    (multiple impuT2T blocks before/after count as one combined reference span).

    Fail if segment_len / ref_span > endpoint_threshold, unless that segment is
    edge-log protected (then print and do not fail for that endpoint).
    """
    fi = find_first_segment_index(entries)
    li = find_last_segment_index(entries)
    fails = False

    # --- First segment vs cumulative prefix W-span ---
    if fi is None:
        if DEBUG_ENDPOINT_FLANK_CHECKS:
            print(f"{patch_label} [endpoint-check] first: no segment W-line — skipped")
    else:
        prefix_len = _w_span_before_index(entries, fi)
        if prefix_len <= 0:
            if DEBUG_ENDPOINT_FLANK_CHECKS:
                print(
                    f"{patch_label} [endpoint-check] first: no W-lines before first segment — skipped"
                )
        else:
            seg_len = entries[fi].length
            r = seg_len / prefix_len
            passes = r <= endpoint_threshold
            if DEBUG_ENDPOINT_FLANK_CHECKS:
                print(
                    f"{patch_label} [endpoint-check] first: segment '{entries[fi].component_id}' "
                    f"len={seg_len} vs prefix W-span={prefix_len} "
                    f"ratio={r:.4f} threshold={endpoint_threshold} → "
                    f"{'PASS' if passes else 'FAIL threshold'}"
                )
            if not passes:
                prot = segment_at_index_protected_by_confident_edge(
                    entries, fi, edge_probs, prob_threshold
                )
                edge = get_impuT2T_edge_across_segment_at_index(entries, fi)
                if prot:
                    ep = _get_edge_prob(edge_probs, edge[0], edge[1]) if edge else None
                    print(
                        f"{patch_label}: first segment vs prefix W-span {prefix_len} "
                        f"ratio={r:.4f} (>{endpoint_threshold}) — would split, but edge-protected"
                        + (
                            f" ({edge[0]} -> {edge[1]}, p={ep:.3f})"
                            if edge and ep is not None
                            else (" (" + edge[0] + " -> " + edge[1] + ")" if edge else "")
                        )
                    )
                else:
                    print(
                        f"{patch_label}: first segment vs prefix W-span {prefix_len} "
                        f"ratio={r:.4f} (>{endpoint_threshold}) — endpoint flank fail"
                    )
                    fails = True

    # --- Last segment vs cumulative suffix W-span ---
    if li is None:
        if DEBUG_ENDPOINT_FLANK_CHECKS:
            print(f"{patch_label} [endpoint-check] last: no segment W-line — skipped")
    else:
        suffix_len = _w_span_after_index(entries, li)
        if suffix_len <= 0:
            if DEBUG_ENDPOINT_FLANK_CHECKS:
                print(
                    f"{patch_label} [endpoint-check] last: no W-lines after last segment — skipped"
                )
        else:
            seg_len = entries[li].length
            r = seg_len / suffix_len
            passes = r <= endpoint_threshold
            if DEBUG_ENDPOINT_FLANK_CHECKS:
                print(
                    f"{patch_label} [endpoint-check] last: segment '{entries[li].component_id}' "
                    f"len={seg_len} vs suffix W-span={suffix_len} "
                    f"ratio={r:.4f} threshold={endpoint_threshold} → "
                    f"{'PASS' if passes else 'FAIL threshold'}"
                )
            if not passes:
                prot = segment_at_index_protected_by_confident_edge(
                    entries, li, edge_probs, prob_threshold
                )
                edge = get_impuT2T_edge_across_segment_at_index(entries, li)
                if prot:
                    ep = _get_edge_prob(edge_probs, edge[0], edge[1]) if edge else None
                    print(
                        f"{patch_label}: last segment vs suffix W-span {suffix_len} "
                        f"ratio={r:.4f} (>{endpoint_threshold}) — would split, but edge-protected"
                        + (
                            f" ({edge[0]} -> {edge[1]}, p={ep:.3f})"
                            if edge and ep is not None
                            else (" (" + edge[0] + " -> " + edge[1] + ")" if edge else "")
                        )
                    )
                else:
                    print(
                        f"{patch_label}: last segment vs suffix W-span {suffix_len} "
                        f"ratio={r:.4f} (>{endpoint_threshold}) — endpoint flank fail"
                    )
                    fails = True

    return fails


def patch_segment_flank_checks_fail(
    entries: List[AgpEntry],
    segment_flank_ratio: float,
    endpoint_segment_flank_ratio: float,
    edge_probs: Dict[Tuple[str, str], float],
    patch_label: str,
    use_endpoint_filter: bool = True,
) -> bool:
    """True if longest-segment flank fails OR any unprotected endpoint flank fails."""
    # Run both checks always so first/last endpoint diagnostics still print even when
    # the longest-segment flank rule fails first (short-circuit hid endpoint messages).
    longest_fails = segment_flank_ratio_exceeds(entries, segment_flank_ratio)
    endpoint_fails = False
    if use_endpoint_filter:
        endpoint_fails = endpoint_segment_flank_exceeds(
            entries, endpoint_segment_flank_ratio, edge_probs, patch_label
        )
    return bool(longest_fails or endpoint_fails)


def _impuT2T_edge_nodes(prev: AgpEntry, curr: AgpEntry) -> Tuple[str, str]:
    """
    Return (u, v) for the edge from prev's right end to curr's left end.
    Convention: '+' => right=id_e, left=id_b; '-' => right=id_b, left=id_e.
    """
    u = f"{prev.component_id}_e" if prev.orientation == "+" else f"{prev.component_id}_b"
    v = f"{curr.component_id}_b" if curr.orientation == "+" else f"{curr.component_id}_e"
    return (u, v)


def iter_impuT2T_edges(entries: List[AgpEntry]) -> Iterator[Tuple[str, str]]:
    """Yield (u, v) for each consecutive impuT2T-impuT2T edge in entries (skipping segments in between)."""
    prev_imp: Optional[AgpEntry] = None
    for e in entries:
        if e.component_type != "W" or classify_component(e.component_id) != "impuT2T":
            continue
        if prev_imp is not None:
            yield _impuT2T_edge_nodes(prev_imp, e)
        prev_imp = e


def _get_edge_prob(edge_probs: Dict[Tuple[str, str], float], u: str, v: str) -> Optional[float]:
    """Look up probability for edge (u, v) or (v, u)."""
    if not edge_probs:
        return None
    prob = edge_probs.get((u, v))
    return prob if prob is not None else edge_probs.get((v, u))


def load_edge_log(path: str) -> Tuple[Dict[Tuple[str, str], float], List[Tuple[str, List[str]]]]:
    """
    Load edge log file. Returns (edge_probs, lines) where edge_probs maps (u,v)/(v,u) to
    probability and lines is a list of (raw_line, fields) for writing filtered.edge.log.
    """
    edge_probs: Dict[Tuple[str, str], float] = {}
    lines: List[Tuple[str, List[str]]] = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                lines.append((line, []))
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                lines.append((line, fields))
                continue
            u, v = fields[0], fields[1]
            try:
                prob = float(fields[-2])
            except ValueError:
                lines.append((line, fields))
                continue
            edge_probs[(u, v)] = prob
            edge_probs[(v, u)] = prob
            lines.append((line, fields))
    return edge_probs, lines


def get_impuT2T_edge_across_segment_at_index(
    entries: List[AgpEntry],
    idx: int,
) -> Optional[Tuple[str, str]]:
    """
    Return (u, v) for the impuT2T-impuT2T adjacency that would be broken if the
    segment W-line at entries[idx] were removed (nearest impuT2T on each side).
    """
    if idx < 0 or idx >= len(entries):
        return None
    prev_imp: Optional[AgpEntry] = None
    for i in range(idx - 1, -1, -1):
        e = entries[i]
        if e.component_type == "W" and classify_component(e.component_id) == "impuT2T":
            prev_imp = e
            break
    next_imp: Optional[AgpEntry] = None
    for i in range(idx + 1, len(entries)):
        e = entries[i]
        if e.component_type == "W" and classify_component(e.component_id) == "impuT2T":
            next_imp = e
            break
    if prev_imp is None or next_imp is None:
        return None
    return _impuT2T_edge_nodes(prev_imp, next_imp)


def segment_at_index_protected_by_confident_edge(
    entries: List[AgpEntry],
    idx: int,
    edge_probs: Dict[Tuple[str, str], float],
    prob_threshold: float = 0.5,
) -> bool:
    """True iff edge across segment at idx exists in edge log with prob >= threshold."""
    edge = get_impuT2T_edge_across_segment_at_index(entries, idx)
    if edge is None:
        return False
    if not edge_probs:
        return False
    u, v = edge
    prob = _get_edge_prob(edge_probs, u, v)
    return prob is not None and prob >= prob_threshold


def find_longest_unprotected_segment_over_threshold(
    entries: List[AgpEntry],
    edge_probs: Dict[Tuple[str, str], float],
    threshold_bp: int,
    prob_threshold: float = 0.5,
) -> Optional[int]:
    """
    Index of the longest 'segment' W-line longer than threshold_bp that is not
    protected by a confident (>= prob_threshold) edge in edge_probs. Without an
    edge log, no long segment is protected.
    """
    best_idx: Optional[int] = None
    best_len = -1
    for i, e in enumerate(entries):
        if e.component_type != "W":
            continue
        if classify_component(e.component_id) != "segment":
            continue
        if e.length <= threshold_bp:
            continue
        if segment_at_index_protected_by_confident_edge(
            entries, i, edge_probs, prob_threshold
        ):
            continue
        if e.length > best_len:
            best_len = e.length
            best_idx = i
    return best_idx


def print_long_segment_threshold_report(
    entries: List[AgpEntry],
    edge_probs: Dict[Tuple[str, str], float],
    threshold_bp: int,
    context: str,
    prob_threshold: float = 0.5,
) -> None:
    """Print each segment W-line longer than threshold_bp and whether it is edge-protected."""
    for i, e in enumerate(entries):
        if e.component_type != "W":
            continue
        if classify_component(e.component_id) != "segment":
            continue
        if e.length <= threshold_bp:
            continue
        protected = segment_at_index_protected_by_confident_edge(
            entries, i, edge_probs, prob_threshold
        )
        edge_pair = get_impuT2T_edge_across_segment_at_index(entries, i)
        if edge_pair is not None:
            u, v = edge_pair
            p = _get_edge_prob(edge_probs, u, v)
            if p is None:
                prob_note = "no entry in edge log"
            else:
                prob_note = f"p={p:.3f}"
            bridge = f"{u} -> {v} ({prob_note})"
        else:
            bridge = "no flanking impuT2T pair (cannot define bridging edge)"
        print(
            f"{context}: long segment '{e.component_id}' "
            f"object {e.object_beg}-{e.object_end} len={e.length} "
            f"(threshold {threshold_bp} bp) protected={protected} | {bridge}"
        )


def get_broken_patch_edge_with_prob(
    entries: List[AgpEntry],
    edge_probs: Dict[Tuple[str, str], float],
) -> Optional[Tuple[str, str, Optional[float]]]:
    """
    Return the single impuT2T-impuT2T edge that would be broken if this patch
    were disassembled by removing its largest segment, along with its
    probability from the edge log (or None if missing).
    """
    idx = find_largest_segment_index(entries)
    if idx is None:
        return None

    prev_imp: Optional[AgpEntry] = None
    for i in range(idx - 1, -1, -1):
        e = entries[i]
        if e.component_type != "W" or classify_component(e.component_id) != "impuT2T":
            continue
        prev_imp = e
        break

    next_imp: Optional[AgpEntry] = None
    for i in range(idx + 1, len(entries)):
        e = entries[i]
        if e.component_type != "W" or classify_component(e.component_id) != "impuT2T":
            continue
        next_imp = e
        break

    if prev_imp is None or next_imp is None:
        return None

    u, v = _impuT2T_edge_nodes(prev_imp, next_imp)
    return (u, v, _get_edge_prob(edge_probs, u, v))


def patch_has_high_prob_edge(
    entries: List[AgpEntry],
    edge_probs: Dict[Tuple[str, str], float],
    prob_threshold: float = 0.5,
) -> bool:
    """
    Return True if the specific impuT2T-impuT2T edge that would be broken by
    removing the largest segment has probability >= prob_threshold.
    """
    broken = get_broken_patch_edge_with_prob(entries, edge_probs)
    if broken is None:
        return False
    _, _, prob = broken
    return prob is not None and prob >= prob_threshold


def count_w_components(entries: List[AgpEntry]) -> int:
    return sum(1 for e in entries if e.component_type == "W")


def _should_skip_single_segment_part(
    part: List[AgpEntry], total_patch_len: int, min_contig_frac: float
) -> bool:
    """Skip appending this part if it has exactly one segment and that segment's length ratio < min_contig_frac."""
    segment_entries = [e for e in part if e.component_type == "W"]
    if len(segment_entries) != 1:
        return False
    segment_len = segment_entries[0].length
    ratio = (segment_len / total_patch_len) if total_patch_len > 0 else 0.0
    #print(f"segment_len: {segment_len}, total_patch_len: {total_patch_len}, ratio: {ratio}, {min_contig_frac}, {ratio < min_contig_frac}")
    return ratio < min_contig_frac


def disassemble_patch_at_index(
    entries: List[AgpEntry], idx: int, min_contig_frac: float
) -> List[List[AgpEntry]]:
    """Split patch by dropping the segment W-line at idx (left and right sub-patches)."""
    if idx < 0 or idx >= len(entries):
        return []
    mid = entries[idx]
    if mid.component_type != "W" or classify_component(mid.component_id) != "segment":
        return []
    total_patch_len = sum(e.length for e in entries if e.component_type == "W")

    left = [e for e in entries[:idx] if e.component_type in ("W", "N")]
    right = [e for e in entries[idx + 1 :] if e.component_type in ("W", "N")]
    parts: List[List[AgpEntry]] = []
    if count_w_components(left) >= 2 or not _should_skip_single_segment_part(
        left, total_patch_len, min_contig_frac
    ):
        parts.append(left)
    if count_w_components(right) >= 2 or not _should_skip_single_segment_part(
        right, total_patch_len, min_contig_frac
    ):
        parts.append(right)
    return parts


def disassemble_patch(entries: List[AgpEntry], min_contig_frac: float) -> List[List[AgpEntry]]:
    idx = find_largest_segment_index(entries)
    if idx is None:
        return []
    return disassemble_patch_at_index(entries, idx, min_contig_frac)


def choose_disassembly_split(
    entries: List[AgpEntry],
    edge_probs: Dict[Tuple[str, str], float],
    min_contig_frac: float,
    long_segment_bp: int,
    prob_threshold: float = 0.5,
) -> List[List[AgpEntry]]:
    """
    Prefer removing the longest unprotected segment longer than long_segment_bp;
    otherwise remove the largest segment (existing behavior).
    """
    long_idx = find_longest_unprotected_segment_over_threshold(
        entries, edge_probs, long_segment_bp, prob_threshold
    )
    if long_idx is not None:
        return disassemble_patch_at_index(entries, long_idx, min_contig_frac)
    return disassemble_patch(entries, min_contig_frac)


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def read_fasta_sequences(fasta_path: str) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    current: Optional[str] = None
    with open_maybe_gzip(fasta_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            else:
                if current is not None:
                    seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def extract_patch_sequence(patch_name: str, entries: List[AgpEntry], patch_sequences: Dict[str, str]) -> str:
    """Extract sequence for a patch block from the patch sequences FASTA."""
    # Get the full patch sequence
    full_seq = patch_sequences.get(patch_name)
    if full_seq is None:
        return ""
    
    # Calculate the start and end positions for this patch block
    start_pos = entries[0].object_beg - 1  # Convert to 0-based
    end_pos = entries[-1].object_end  # Already 1-based, so use as-is
    
    # Extract the sequence for this specific range
    return full_seq[start_pos:end_pos]


def _part_w_span_total(part: List[AgpEntry]) -> int:
    return sum(e.length for e in part if e.component_type == "W")


def _nonempty_part_indices(parts: List[List[AgpEntry]]) -> List[int]:
    return [i for i, p in enumerate(parts) if _part_w_span_total(p) > 0]


def _stem_for_part_sibling(
    base_name: str, obj: str, part_idx: int, num_nonempty_parts: int
) -> str:
    """If only one nonempty sibling, keep base_name; else append _partNN under base_name/obj."""
    if num_nonempty_parts <= 1:
        return base_name
    if base_name != obj:
        return f"{base_name}_part{part_idx:02d}"
    return f"{obj}_part{part_idx:02d}"


def _kept_part_object_name(
    part_entries: List[AgpEntry],
    part_idx: int,
    num_nonempty_parts: int,
    obj: str,
    base_name: str,
) -> str:
    """
    - Single W-line that is impuT2T -> use that component_id.
    - Only one surviving nonempty part from this split -> use base_name (no _partNN).
    - Else -> base_name_partNN / obj_partNN (legacy).
    """
    w_only = [e for e in part_entries if e.component_type == "W"]
    if (
        len(w_only) == 1
        and classify_component(w_only[0].component_id) == "impuT2T"
    ):
        return w_only[0].component_id
    return _stem_for_part_sibling(base_name, obj, part_idx, num_nonempty_parts)


def reindex_agp_entries(entries: List[AgpEntry], object_name: str) -> List[AgpEntry]:
    new_entries: List[AgpEntry] = []
    cursor = 1
    part = 1
    for e in entries:
        length = e.length
        if e.component_type not in ("W", "N"):
            continue
        new_entries.append(
            AgpEntry(
                object_name,
                cursor,
                cursor + length - 1,
                part,
                e.component_type,
                e.component_id,
                e.component_beg,
                e.component_end,
                e.orientation,
            )
        )
        cursor += length
        part += 1
    return new_entries


def write_filtered_agp(objects: List[Tuple[str, List[AgpEntry]]], output_path: str, write_header: bool) -> None:
    with open(output_path, "w") as out:
        if write_header:
            out.write("## agp-version 2.0\n")
            out.write("# Filtered by extract_fasta_patch_ratio.py\n")
        for obj_name, entries in objects:
            for e in entries:
                if e.component_type == "W":
                    out.write(
                        f"{obj_name}\t{e.object_beg}\t{e.object_end}\t{e.part_num}\tW\t{e.component_id}\t{e.component_beg}\t{e.component_end}\t{e.orientation}\n"
                    )
                elif e.component_type == "N":
                    # Write a simple gap line; column semantics simplified
                    out.write(
                        f"{obj_name}\t{e.object_beg}\t{e.object_end}\t{e.part_num}\tN\t{e.length}\tscaffold\tyes\tpaired-ends\n"
                    )


def _process_parts_recursive(
    parts: List[List[AgpEntry]],
    obj: str,
    base_name: str,
    kept_objects: List[Dict[str, object]],
    longest_kept_patch_len: int,
    ratio_threshold: float,
    segment_flank_ratio: float,
    endpoint_segment_flank_ratio: float,
    use_endpoint_filter: bool,
    min_contig_frac: float,
    edge_probs: Dict[Tuple[str, str], float],
    long_segment_bp: int,
) -> int:
    """Recursively keep or disassemble parts: keep if ratio and segment-flank ratio OK, else split again."""
    nonempty_idx = _nonempty_part_indices(parts)
    num_nonempty = len(nonempty_idx)
    for part_idx, part_entries in enumerate(parts):
        p_total, _p_seg, p_ratio = compute_segment_ratio(part_entries)
        if p_total <= 0:
            continue
        long_idx = find_longest_unprotected_segment_over_threshold(
            part_entries, edge_probs, long_segment_bp, prob_threshold=0.5
        )
        # Avoid `a and not f()` short-circuit: if segment ratio fails, we must still run flank/endpoint checks.
        ratio_ok = p_ratio <= ratio_threshold
        flank_ok = not patch_segment_flank_checks_fail(
            part_entries,
            segment_flank_ratio,
            endpoint_segment_flank_ratio,
            edge_probs,
            base_name,
            use_endpoint_filter,
        )
        passes_inner = ratio_ok and flank_ok
        if passes_inner and long_idx is None:
            new_name = _kept_part_object_name(
                part_entries, part_idx, num_nonempty, obj, base_name
            )
            orig_start = part_entries[0].object_beg
            orig_end = part_entries[-1].object_end
            new_entries = reindex_agp_entries(
                [e for e in part_entries if e.component_type in ("W", "N")], new_name
            )
            kept_objects.append({
                "name": new_name,
                "entries": new_entries,
                "source": obj,
                "orig_start": orig_start,
                "orig_end": orig_end,
            })
            # Use the shortest kept patch length as reference.
            if longest_kept_patch_len == 0 or p_total < longest_kept_patch_len:
                longest_kept_patch_len = p_total
        else:
            sub_parts = choose_disassembly_split(
                part_entries,
                edge_probs,
                min_contig_frac,
                long_segment_bp,
                prob_threshold=0.5,
            )
            if sub_parts:
                child_base = _stem_for_part_sibling(
                    base_name, obj, part_idx, num_nonempty
                )
                longest_kept_patch_len = _process_parts_recursive(
                    sub_parts,
                    obj,
                    child_base,
                    kept_objects,
                    longest_kept_patch_len,
                    ratio_threshold,
                    segment_flank_ratio,
                    endpoint_segment_flank_ratio,
                    use_endpoint_filter,
                    min_contig_frac,
                    edge_probs,
                    long_segment_bp,
                )
            # else: cannot split further, part is discarded
    return longest_kept_patch_len


def main() -> None:
    args = parse_args()

    per_object = parse_agp(args.agp)

    # Each kept object: {name, entries, source, orig_start, orig_end}
    kept_objects: List[Dict[str, object]] = []
    # Track a reference patch length. Historically this tracked the longest kept
    # patch; now we intentionally track the *shortest* kept patch length so that
    # the min_contig_len threshold is set relative to the smallest accepted
    # patch. A value of 0 still means "no patch kept yet".
    longest_kept_patch_len = 0

    # Load edge log once if provided (used for protection during patch filtering and for filtered.edge.log).
    edge_probs: Dict[Tuple[str, str], float] = {}
    edge_log_lines: List[Tuple[str, List[str]]] = []
    if args.edge_info:
        edge_probs, edge_log_lines = load_edge_log(args.edge_info)

    # Process patch objects first
    for obj, entries in per_object.items():
        if not obj.startswith("patch"):
            continue
        total_len, seg_len, ratio = compute_segment_ratio(entries)
        # Protect patches that contain at least one high-probability edge:
        # such patches should not be disassembled or dropped based on segment
        # ratio / flank ratio alone.
        protect_patch = False
        if edge_probs:
            if patch_has_high_prob_edge(entries, edge_probs, prob_threshold=0.5):
                protect_patch = True

        # Avoid short-circuit and keep split triggers explicit:
        # - endpoint fails are hard failures (must split), even if patch is edge-protected.
        # - longest-flank fails can be overridden by patch-level edge protection.
        ratio_ok = ratio <= args.ratio_threshold
        longest_flank_fails = segment_flank_ratio_exceeds(entries, args.segment_flank_ratio)
        endpoint_fails = False
        if not args.disable_endpoint_segment_flank_filter:
            endpoint_fails = endpoint_segment_flank_exceeds(
                entries,
                args.endpoint_segment_flank_ratio,
                edge_probs,
                obj,
            )
        flank_ok = not (longest_flank_fails or endpoint_fails)
        passes_ratio = ratio_ok and flank_ok

        long_unprot_idx = find_longest_unprotected_segment_over_threshold(
            entries, edge_probs, args.long_segment_bp, prob_threshold=0.5
        )
        #print_long_segment_threshold_report(
        #    entries,
        #    edge_probs,
        #    args.long_segment_bp,
        #    context=f"Patch {obj}",
        #    prob_threshold=0.5,
        #)

        # Print hints when edge information actually changes the outcome:
        # - If the patch would fail ratio/flank filters but is kept due to a
        #   high-probability edge, report that specific edge.
        # - If the patch is filtered (fails ratio/flank and is not protected),
        #   report which edge across the largest segment is being broken
        #   (without probability if no edge log).
        if not passes_ratio:
            broken = get_broken_patch_edge_with_prob(entries, edge_probs)
            if protect_patch and broken is not None:
                u, v, prob = broken
                print(
                    f"Patch {obj} would fail ratio/flank filters (ratio={ratio:.4f}) "
                    f"but is protected by edge {u} -> {v} with probability "
                    f"{0.0 if prob is None else prob:.3f}"
                )
            elif (not protect_patch) and broken is not None:
                u, v, prob = broken
                if prob is not None:
                    print(
                        f"Filtering patch {obj} by ratio/flank (ratio={ratio:.4f}); "
                        f"breaking edge {u} -> {v} with probability {prob:.3f}"
                    )
                else:
                    print(
                        f"Filtering patch {obj} by ratio/flank (ratio={ratio:.4f}); "
                        f"breaking edge {u} -> {v}"
                    )

        # Keep full patch if:
        # - edge-protected and endpoint checks pass (protection overrides ratio/longest-flank/long-segment split), OR
        # - unprotected but fully passes ratio+flank and has no long unprotected segment over threshold.
        if (protect_patch and (not endpoint_fails)) or (
            (not protect_patch) and flank_ok and ratio_ok and long_unprot_idx is None
        ):
            # Keep full patch (no unprotected segment longer than --long-segment-bp)
            new_name = obj
            # Original coordinate span in the full patch
            orig_start = entries[0].object_beg
            orig_end = entries[-1].object_end
            new_entries = reindex_agp_entries([e for e in entries if e.component_type in ("W", "N")], new_name)
            kept_objects.append({
                "name": new_name,
                "entries": new_entries,
                "source": obj,
                "orig_start": orig_start,
                "orig_end": orig_end,
            })
            # Use the shortest kept patch length as reference.
            if longest_kept_patch_len == 0 or total_len < longest_kept_patch_len:
                longest_kept_patch_len = total_len
            continue

        # Disassemble: prefer longest unprotected segment > --long-segment-bp, else largest segment
        parts = choose_disassembly_split(
            entries,
            edge_probs,
            args.min_contig_frac,
            args.long_segment_bp,
            prob_threshold=0.5,
        )
        longest_kept_patch_len = _process_parts_recursive(
            parts,
            obj,
            obj,
            kept_objects,
            longest_kept_patch_len,
            args.ratio_threshold,
            args.segment_flank_ratio,
            args.endpoint_segment_flank_ratio,
            not args.disable_endpoint_segment_flank_filter,
            args.min_contig_frac,
            edge_probs,
            args.long_segment_bp,
        )

    # Process individual contig objects
    min_contig_len = int(longest_kept_patch_len * args.min_contig_frac)
    for obj, entries in per_object.items():
        if obj.startswith("patch"):
            continue
        # Keep only unit objects that are a single W spanning
        w_entries = [e for e in entries if e.component_type == "W"]
        if len(w_entries) != 1:
            continue
        e = w_entries[0]
        if e.length >= min_contig_len:
            new_name = obj
            new_entries = reindex_agp_entries([e], new_name)
            kept_objects.append({
                "name": new_name,
                "entries": new_entries,
                "source": obj,
                "orig_start": e.object_beg,
                "orig_end": e.object_end,
            })

    # Write filtered AGP
    agp_out = args.output_prefix + ".filtered.agp"
    # Adapt writer to new kept_objects structure
    def _write_filtered_agp(objects: List[Dict[str, object]], output_path: str, write_header: bool) -> None:
        with open(output_path, "w") as out:
            if write_header:
                out.write("## agp-version 2.0\n")
                out.write("# Filtered by extract_fasta_patch_ratio.py\n")
            for obj in objects:
                obj_name = obj["name"]  # type: ignore[index]
                for e in obj["entries"]:  # type: ignore[index]
                    if e.component_type == "W":
                        out.write(
                            f"{obj_name}\t{e.object_beg}\t{e.object_end}\t{e.part_num}\tW\t{e.component_id}\t{e.component_beg}\t{e.component_end}\t{e.orientation}\n"
                        )
                    elif e.component_type == "N":
                        out.write(
                            f"{obj_name}\t{e.object_beg}\t{e.object_end}\t{e.part_num}\tN\t{e.length}\tscaffold\tyes\tpaired-ends\n"
                        )

    _write_filtered_agp(kept_objects, agp_out, write_header=not args.no_header)

    # Build filtered FASTA
    patch_sequences = read_fasta_sequences(args.fasta)
    fasta_out = args.output_prefix + ".filtered.fasta"
    
    with open(fasta_out, "w") as out_fa:
        for obj in kept_objects:
            print(obj["name"], obj["source"], obj["orig_start"], obj["orig_end"])
            obj_name = obj["name"]  # type: ignore[index]
            # Use the original source patch name and original coordinate span
            source = obj["source"]  # type: ignore[index]
            orig_start = int(obj["orig_start"])  # type: ignore[index]
            orig_end = int(obj["orig_end"])  # type: ignore[index]
            full_seq = patch_sequences.get(str(source), "")
            seq = full_seq[orig_start - 1:orig_end] if full_seq else ""
            if not seq:
                print(f"No sequence found for {obj_name}")
                continue
            out_fa.write(f">{obj_name}\n")
            out_fa.write(seq + "\n")

    # Write filtered edge log when edge info was loaded (single source: edge_log_lines).
    if args.edge_info and edge_log_lines:
        used_edges_after: set[tuple[str, str]] = set()
        for obj in kept_objects:
            obj_name = obj["name"]  # type: ignore[index]
            if not str(obj_name).startswith("patch"):
                continue
            for u, v in iter_impuT2T_edges(obj["entries"]):  # type: ignore[index]
                if (u, v) not in edge_probs:
                    raise ValueError(
                        f"Edge ({u}, {v}) used in filtered AGP has no counterpart in edge-info log {args.edge_info}"
                    )
                used_edges_after.add((u, v))

        edge_out = args.output_prefix + ".filtered.edge.log"
        with open(edge_out, "w") as out_edge:
            for line, fields in edge_log_lines:
                if len(fields) < 3:
                    out_edge.write(line + "\n" if line else "\n")
                    continue
                u, v = fields[0], fields[1]
                used_flag = fields[-1]
                try:
                    prob = float(fields[-2])
                except ValueError:
                    out_edge.write(line + "\n")
                    continue
                still_used = (u, v) in used_edges_after or (v, u) in used_edges_after
                if used_flag == "True" and (not still_used) and prob < 0.5:
                    fields[-1] = "False"
                out_edge.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()


