#!/usr/bin/env python3

import argparse
import io
import gzip
import os
import re
from typing import Dict, List, Tuple, Optional, Literal, TextIO


ComponentKind = Literal["impuT2T", "segment", "other"]

IMPUT2T_REGEX = re.compile(r".*_impuT2T_\d+$")
SEGMENT_REGEX = re.compile(r".*_segment_\d+$")


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
    return parser.parse_args()


def classify_component(component_id: str) -> ComponentKind:
    if IMPUT2T_REGEX.match(component_id):
        return "impuT2T"
    if SEGMENT_REGEX.match(component_id):
        return "segment"
    return "other"


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


def disassemble_patch(entries: List[AgpEntry], min_contig_frac: float) -> List[List[AgpEntry]]:
    idx = find_largest_segment_index(entries)
    if idx is None:
        return []
    total_patch_len = sum(e.length for e in entries if e.component_type == "W")

    left = [e for e in entries[:idx] if e.component_type in ("W", "N")]
    right = [e for e in entries[idx + 1 :] if e.component_type in ("W", "N")]
    parts: List[List[AgpEntry]] = []
    #print("Split!", count_w_components(left), count_w_components(right))
    if count_w_components(left) >= 2 or not _should_skip_single_segment_part(left, total_patch_len, min_contig_frac):
        parts.append(left)
    if count_w_components(right) >= 2 or not _should_skip_single_segment_part(right, total_patch_len, min_contig_frac):
        parts.append(right)
    return parts


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


def main() -> None:
    args = parse_args()

    per_object = parse_agp(args.agp)

    # Each kept object: {name, entries, source, orig_start, orig_end}
    kept_objects: List[Dict[str, object]] = []
    longest_kept_patch_len = 0

    # Process patch objects first
    for obj, entries in per_object.items():
        if not obj.startswith("patch"):
            continue
        total_len, seg_len, ratio = compute_segment_ratio(entries)
        if ratio <= args.ratio_threshold:
            # Keep full patch
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
            if total_len > longest_kept_patch_len:
                longest_kept_patch_len = total_len
            continue

        # Disassemble by removing largest segment
        parts = disassemble_patch(entries, args.min_contig_frac)
        for part_idx, part_entries in enumerate(parts):
            p_total, p_seg, p_ratio = compute_segment_ratio(part_entries)
            if p_total <= 0:
                continue
            if p_ratio <= args.ratio_threshold:
                new_name = f"{obj}_part{part_idx:02d}"
                # Original coordinate span for this subpatch in the full patch object
                orig_start = part_entries[0].object_beg
                orig_end = part_entries[-1].object_end
                new_entries = reindex_agp_entries([e for e in part_entries if e.component_type in ("W", "N")], new_name)
                kept_objects.append({
                    "name": new_name,
                    "entries": new_entries,
                    "source": obj,
                    "orig_start": orig_start,
                    "orig_end": orig_end,
                })
                if p_total > longest_kept_patch_len:
                    longest_kept_patch_len = p_total

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
            obj_name = obj["name"]  # type: ignore[index]
            if str(obj_name).startswith("patch"):
                # Use the original source patch name and original coordinate span
                source = obj["source"]  # type: ignore[index]
                orig_start = int(obj["orig_start"])  # type: ignore[index]
                orig_end = int(obj["orig_end"])  # type: ignore[index]
                full_seq = patch_sequences.get(str(source), "")
                seq = full_seq[orig_start - 1:orig_end] if full_seq else ""
            else:
                # For individual contigs, get the full sequence by name
                seq = patch_sequences.get(str(obj_name), "")
            
            if not seq:
                continue
            out_fa.write(f">{obj_name}\n")
            out_fa.write(seq + "\n")


if __name__ == "__main__":
    main()


