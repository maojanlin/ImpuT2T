#!/home/mlin77/miniconda3/envs/conda_ragoo/bin/python

"""
MIT License

Copyright (c) 2021 Michael Alonge <malonge11@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Modified by Mao-Jan Lin in 2025 <mj.maojanlin@gmail.com> for the impuT2T patching project
"""

import os
import sys
import argparse
import subprocess
import time

import pysam

from ragtag_utilities.utilities import log, run_oae, get_ragtag_version
from ragtag_utilities.AlignmentReader import PAFReader
from ragtag_utilities.ContigAlignment import ContigAlignment
from ragtag_utilities.Aligner import Minimap2Aligner
from ragtag_utilities.Aligner import UnimapAligner
from ragtag_utilities.Aligner import NucmerAligner
from ragtag_utilities.ScaffoldGraph import AGPMultiScaffoldGraph
from ragtag_utilities.ScaffoldGraph import PatchScaffoldGraph
from ragtag_utilities.ScaffoldGraph import Alignment


def read_genome_alignments(aln_file, flag_primary=False, query_blacklist=[], ref_blacklist=[]):
    tmp_ctg_alns = dict()
    aln_reader = PAFReader(aln_file)
    for aln_line in aln_reader.parse_alignments():
        if flag_primary == True and aln_line.flag_primary == False:
            # filter the secondary alignments if the tp:A:P/tp:A:S is specified
            continue
        # Check that the contig and reference in this alignment are allowed.
        if aln_line.query_header not in query_blacklist and aln_line.ref_header not in ref_blacklist:
            if aln_line.query_header not in tmp_ctg_alns:
                tmp_ctg_alns[aln_line.query_header] = [aln_line.query_header, aln_line.query_len,
                                                       [aln_line.query_start], [aln_line.query_end], [aln_line.strand],
                                                       [aln_line.ref_header], [aln_line.ref_len],
                                                       [aln_line.ref_start], [aln_line.ref_end],
                                                       [aln_line.num_match], [aln_line.aln_len],
                                                       [aln_line.mapq]]
            else:
                tmp_ctg_alns[aln_line.query_header][2].append(aln_line.query_start)
                tmp_ctg_alns[aln_line.query_header][3].append(aln_line.query_end)
                tmp_ctg_alns[aln_line.query_header][4].append(aln_line.strand)
                tmp_ctg_alns[aln_line.query_header][5].append(aln_line.ref_header)
                tmp_ctg_alns[aln_line.query_header][6].append(aln_line.ref_len)
                tmp_ctg_alns[aln_line.query_header][7].append(aln_line.ref_start)
                tmp_ctg_alns[aln_line.query_header][8].append(aln_line.ref_end)
                tmp_ctg_alns[aln_line.query_header][9].append(aln_line.num_match)
                tmp_ctg_alns[aln_line.query_header][10].append(aln_line.aln_len)
                tmp_ctg_alns[aln_line.query_header][11].append(aln_line.mapq)

    ctg_alns = dict()
    for i in tmp_ctg_alns:
        ctg_alns[i] = ContigAlignment(
            tmp_ctg_alns[i][0],
            tmp_ctg_alns[i][1],
            tmp_ctg_alns[i][2],
            tmp_ctg_alns[i][3],
            tmp_ctg_alns[i][4],
            tmp_ctg_alns[i][5],
            tmp_ctg_alns[i][6],
            tmp_ctg_alns[i][7],
            tmp_ctg_alns[i][8],
            tmp_ctg_alns[i][9],
            tmp_ctg_alns[i][10],
            tmp_ctg_alns[i][11]
        )
    return ctg_alns


def build_aln_scaffold_graph(ctg_alns, components_fn, max_term_dist, terminal_matches, debug=False):
    """
    Build a directed scaffold graph from filtered alignments
    :param ctg_alns: query sequence -> ContigAlignment object
    :param components_fn: name of FASTA file with all relevant sequences
    :param max_term_dist: maximum alignment distance from a sequence terminus
    :return: PatchScaffoldGraph
    """
    sg = PatchScaffoldGraph(components_fn)
    if debug:
        print(sg, components_fn, sg.edges, sg.nodes)
        for key, item in ctg_alns.items():
            print(key, item)

    for query_seq in ctg_alns:
        als = ctg_alns[query_seq]
        als.sort_by_query()
        last = None
        last_reversed = False

        # dict_ref_end_connectable:
        #   L_end: list of booleans, True if the reference end is connectable to the left end of the reference
        #   R_end: list of booleans, True if the reference end is connectable to the right end of the reference
        #   checking if a contig is split into multiple segments, if so, only boundary segments are connectable for the corresponding ends
        dict_ref_end_connectable = {'L_end': [], 'R_end': []}
        for i in range(als.num_alns):
            if debug:
                print(i, als.ref_headers[i], '(', als.ref_starts[i], als.ref_ends[i],')')
            duplicate_segments = []
            for j in range(als.num_alns):
                if i != j and als.ref_headers[i] == als.ref_headers[j]:
                    duplicate_segments.append(j)
            if len(duplicate_segments) > 0:
                # multiple segments, can only be connected to one end
                dict_ref_end_connectable['L_end'].append(False)
                dict_ref_end_connectable['R_end'].append(False)

                if als.ref_starts[i] <= min(als.ref_starts[ele] for ele in duplicate_segments):
                    if als.strands[i] == '+':
                        dict_ref_end_connectable['L_end'][i] = True
                    else:
                        dict_ref_end_connectable['R_end'][i] = True
                if als.ref_ends[i] >= max(als.ref_ends[ele] for ele in duplicate_segments):
                    if als.strands[i] == '+':
                        dict_ref_end_connectable['R_end'][i] = True
                    else:
                        dict_ref_end_connectable['L_end'][i] = True
            else:
                # single segment, can be connected to both ends
                dict_ref_end_connectable['L_end'].append(True)
                dict_ref_end_connectable['R_end'].append(True)
        if debug:
            print(dict_ref_end_connectable['L_end'])
            print(dict_ref_end_connectable['R_end'])
            for i in range(als.num_alns):
                print(i, als.ref_headers[i], dict_ref_end_connectable['L_end'][i], dict_ref_end_connectable['R_end'][i], als.ref_starts[i], als.ref_ends[i], '\t|\t', als.strands[i], als.query_starts[i], als.query_ends[i])

        # Iterate over each alignment for this query sequence
        for i in range(als.num_alns):
            # Don't connect irrelevant middle alignments
            if not dict_ref_end_connectable['L_end'][i] and not dict_ref_end_connectable['R_end'][i]:
                continue

            cur_reversed = False

            # For each reference/query alignment terminus, determine if it is close to the sequence terminus
            ref_left_end, ref_right_end = als.ref_start_end(i, max_term_dist)
            query_left_end, query_right_end = als.query_start_end(i, max_term_dist)

            #print("++++++++++++++++++++", i, als.ref_headers[i], dict_ref_end_connectable['L_end'][i], dict_ref_end_connectable['R_end'][i], '(', ref_left_end, ref_right_end,')')
            #print("++++++++++++++++++++", i, als.ref_headers[i], dict_ref_end_connectable['L_end'][i], dict_ref_end_connectable['R_end'][i], '(', query_left_end, query_right_end,')')

            if als.strands[i] == '-':
                cur_reversed = True
            """
            # Determine if we are reversing the reference sequence of this alignment
            if ref_left_end and ref_right_end:
                # Entire reference aligns - have to look at strand in alignment
                if als.strands[i] == '-':
                    cur_reversed = True
            elif ref_left_end:
                # Beginning of reference - we would expect suffix of read if same strand
                if query_left_end and query_right_end:
                    if als.strands[i] == '-':
                        cur_reversed = True
                elif query_left_end:
                    cur_reversed = True
            else:
                # End of contig - we would expect prefix of read if same strand
                if query_left_end and query_right_end:
                    if als.strands[i] == '-':
                        cur_reversed = True
                elif query_right_end:
                    cur_reversed = True
            """

            if last is not None:
                if debug:
                    print("====================", i, last, als.ref_headers[last], '(', als.strands[last], als.ref_starts[last], als.ref_ends[last],')',  als.ref_headers[i], '(', als.strands[i], als.ref_starts[i], als.ref_ends[i],')')
                if als.ref_headers[last] != als.ref_headers[i]:
                    if debug:
                        print("====================", i, last, als.ref_headers[last], '(', als.strands[last], als.ref_starts[last], als.ref_ends[last],')',  als.ref_headers[i], '(', als.strands[i], als.ref_starts[i], als.ref_ends[i],')')
                    if dict_ref_end_connectable['R_end'][last] == False or dict_ref_end_connectable['L_end'][i] == False:
                        if debug:
                            print("Not connectable: my_right_end", dict_ref_end_connectable['R_end'][last], "their_left_end", dict_ref_end_connectable['L_end'][i])
                    else:
                        my_query_end_offset = als.ref_lens[last] - als.ref_ends[last]
                        if last_reversed:
                            my_query_end_offset = als.ref_starts[last]

                        their_query_start_offset = als.ref_starts[i]
                        if cur_reversed:
                            their_query_start_offset = als.ref_lens[i] - als.ref_ends[i]

                        my_query_end = als.query_ends[last] + my_query_end_offset
                        their_query_start = als.query_starts[i] - their_query_start_offset
                        overlap = my_query_end - their_query_start

                        #print("....................", last_reversed, cur_reversed, my_query_end_offset, als.query_ends[last], their_query_start_offset, als.query_starts[i])
                        #print("....................", overlap, my_query_end, their_query_start, als.ref_lens[last], als.ref_lens[i])
                        if overlap <= als.ref_lens[last] and overlap <= als.ref_lens[i]:
                            if debug:
                                print("Process the connection !!!!!!!!!!!!!!!!!!!!!!", overlap, als.ref_lens[last], als.ref_lens[i])
                                print(terminal_matches.keys())
                            #print(last_reversed, cur_reversed)
                            #print(terminal_matches[als.ref_headers[last]][0].residue_matches, terminal_matches[als.ref_headers[i]][0].residue_matches)
                            #print(terminal_matches[als.ref_headers[last]][1].residue_matches, terminal_matches[als.ref_headers[i]][1].residue_matches)

                            # Determine the scaffold graph nodes
                            u = als.ref_headers[last] + "_e"
                            v = als.ref_headers[i] + "_b"
                            if last_reversed:
                                u = als.ref_headers[last] + "_b"
                            if cur_reversed:
                                v = als.ref_headers[i] + "_e"

                            pair_score = 0
                            if last_reversed:
                                pair_score += terminal_matches[als.ref_headers[last]][0]
                                print("last_reversed", pair_score)
                            else:
                                pair_score += terminal_matches[als.ref_headers[last]][1]
                                print("last_forward", pair_score)
                            if cur_reversed:
                                pair_score += terminal_matches[als.ref_headers[i]][1]
                                print("cur_reversed", terminal_matches[als.ref_headers[i]][1])
                            else:
                                pair_score += terminal_matches[als.ref_headers[i]][0]
                                print("cur_forward", terminal_matches[als.ref_headers[i]][0])

                            alignment = Alignment(
                                u,
                                v,
                                query_seq,
                                als.query_len,
                                my_query_end,
                                their_query_start,
                                0,  # Always on the query's forward strand
                                is_gap=False
                            )
                            sg.add_edge(u, v, alignment, pair_score)
            if debug:
                print("!!!!!!!!!!!!!!!!!!!!!! RESET!")
            last = i
            last_reversed = cur_reversed

    sg.remove_heavier_than(2)
    return sg


def filter_high_quality_alns(aln_full_len, minimum_length=10000, length_ratio=0.05, max_min_id_ratio=0.7):
    """filter the full length then map the terminal alignments.
    
    Args:
        aln_full_len (dict): Dictionary of ContigAlignment objects for full length alignments
        
    Returns:
        tuple: (selected_full) where:
            - selected_full is a dict of filtered ContigAlignment objects
    """
    selected_full = {}
    
    # First filter full length alignments to get the longest alignment covering terminal regions
    for query_name, full_aln in aln_full_len.items():
        # filter the full length alignments to get the longest alignment covering terminal regions
        max_identity = 0
        for i in range(full_aln.num_alns):
            if full_aln.residue_matches[i] > minimum_length and full_aln.residue_matches[i]/full_aln.aln_lens[i] > max_identity:
                max_identity = full_aln.residue_matches[i]/full_aln.aln_lens[i]

        selected_indices = []
        print(query_name, max_identity)
        for i in range(full_aln.num_alns):
            if full_aln.residue_matches[i]/full_aln.aln_lens[i] < max_identity*max_min_id_ratio: # 0.7 of max_identity
                continue
            if full_aln.query_len*length_ratio < full_aln.query_ends[i]-full_aln.query_starts[i]: # 0.05 of query length
                selected_indices.append(i)
            elif minimum_length < full_aln.query_ends[i]-full_aln.query_starts[i]: # query length > 10000
                selected_indices.append(i)
            else:
                pass

        if len(selected_indices) > 0:
            selected_full[query_name] = ContigAlignment(
                full_aln.query_header,
                full_aln.query_len,
                [full_aln.query_starts[i] for i in selected_indices],
                [full_aln.query_ends[i] for i in selected_indices],
                [full_aln.strands[i] for i in selected_indices],
                [full_aln.ref_headers[i] for i in selected_indices],
                [full_aln.ref_lens[i] for i in selected_indices],
                [full_aln.ref_starts[i] for i in selected_indices],
                [full_aln.ref_ends[i] for i in selected_indices],
                [full_aln.residue_matches[i] for i in selected_indices],
                [full_aln.aln_lens[i] for i in selected_indices],
                [full_aln.mapqs[i] for i in selected_indices]
            )
    return selected_full



def reverse_alignments(alns):
    """Reverse alignments by swapping query and reference information.
    
    Args:
        alns (dict): Dictionary of ContigAlignment objects keyed by query headers
        
    Returns:
        dict: Dictionary of ContigAlignment objects keyed by reference headers,
             with query and reference information swapped
    """
    # First collect all alignments by reference header
    tmp_ref_alns = {}
    for query_header, aln in alns.items():
        for i in range(aln.num_alns):
            ref_header = aln.ref_headers[i]
            if ref_header not in tmp_ref_alns:
                tmp_ref_alns[ref_header] = [
                    ref_header,                    # Now this is query header
                    aln.ref_lens[i],              # Now this is query len
                    [aln.ref_starts[i]],          # Now these are query starts
                    [aln.ref_ends[i]],            # Now these are query ends
                    [aln.strands[i]],             # Strands stay the same
                    [aln.query_header],           # Now these are ref headers
                    [aln.query_len],              # Now these are ref lens
                    [aln.query_starts[i]],        # Now these are ref starts
                    [aln.query_ends[i]],          # Now these are ref ends
                    [aln.residue_matches[i]],     # These metrics stay the same
                    [aln.aln_lens[i]],           
                    [aln.mapqs[i]]
                ]
            else:
                tmp_ref_alns[ref_header][2].append(aln.ref_starts[i])
                tmp_ref_alns[ref_header][3].append(aln.ref_ends[i])
                tmp_ref_alns[ref_header][4].append(aln.strands[i])
                tmp_ref_alns[ref_header][5].append(aln.query_header)
                tmp_ref_alns[ref_header][6].append(aln.query_len)
                tmp_ref_alns[ref_header][7].append(aln.query_starts[i])
                tmp_ref_alns[ref_header][8].append(aln.query_ends[i])
                tmp_ref_alns[ref_header][9].append(aln.residue_matches[i])
                tmp_ref_alns[ref_header][10].append(aln.aln_lens[i])
                tmp_ref_alns[ref_header][11].append(aln.mapqs[i])
    
    # Convert to ContigAlignment objects
    ref_alns = {}
    for ref_header in tmp_ref_alns:
        ref_alns[ref_header] = ContigAlignment(
            tmp_ref_alns[ref_header][0],  # query_header (former ref_header)
            tmp_ref_alns[ref_header][1],  # query_len (former ref_len)
            tmp_ref_alns[ref_header][2],  # query_starts (former ref_starts)
            tmp_ref_alns[ref_header][3],  # query_ends (former ref_ends)
            tmp_ref_alns[ref_header][4],  # strands (unchanged)
            tmp_ref_alns[ref_header][5],  # ref_headers (former query_headers)
            tmp_ref_alns[ref_header][6],  # ref_lens (former query_lens)
            tmp_ref_alns[ref_header][7],  # ref_starts (former query_starts)
            tmp_ref_alns[ref_header][8],  # ref_ends (former query_ends)
            tmp_ref_alns[ref_header][9],  # residue_matches (unchanged)
            tmp_ref_alns[ref_header][10], # aln_lens (unchanged)
            tmp_ref_alns[ref_header][11]  # mapqs (unchanged)
        )
    
    return ref_alns

def output_graph_info(graph, output_path):
    with open(output_path, "w") as f:
        for u, v in graph.edges:
            f.write("\t".join([u,v, str(graph[u][v]["score"]), "\n"]))


def select_terminals(alns_selected_full, terminal_ratio=0.05, map_lenth=5000):
    # Only select if the sequence contains the first or last 0.05 part of the query
    # default length is 5000

    dict_head_end_range = {}
    dict_tail_end_range = {}
    for query_name, full_aln in alns_selected_full.items():
        head_best_idx = None
        head_best_start = full_aln.query_len
        tail_best_idx = None
        tail_best_end = 0

        for i in range(full_aln.num_alns):
            #0~full_aln.query_len*0.05, and full_aln.query_len*0.95~full_aln.query_len
            if full_aln.query_starts[i] < full_aln.query_len*terminal_ratio:
                if full_aln.query_starts[i] < head_best_start:
                    head_best_start = full_aln.query_starts[i]
                    head_best_idx = i
            if full_aln.query_ends[i] > full_aln.query_len*(1-terminal_ratio):
                if full_aln.query_ends[i] > tail_best_end:
                    tail_best_end = full_aln.query_ends[i]
                    tail_best_idx = i
        extend_length = min(full_aln.query_len, map_lenth)
        if head_best_idx is not None:
            dict_head_end_range[query_name] = (full_aln.query_starts[head_best_idx], full_aln.query_starts[head_best_idx]+extend_length)
        if tail_best_idx is not None:
            dict_tail_end_range[query_name] = (full_aln.query_ends[tail_best_idx]-extend_length, full_aln.query_ends[tail_best_idx])
    return dict_head_end_range, dict_tail_end_range


def write_terminal_fasta(input_fasta, output_fasta,  dict_head_end_range, dict_tail_end_range):
    """
    Extracts head and tail sequences from each record in input_fasta and writes to output_fasta.
    Does not use BioPython.
    """
    def write_record(out_f, header, seq, suffix, start, end):
        subseq = seq[start:end]
        out_f.write(f">{header}_{suffix} {start}-{end}\n{subseq}\n")

    with open(input_fasta, "r") as in_f, open(output_fasta, "w") as out_f:
        header = None
        seq_lines = []
        for line in in_f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    if header in dict_head_end_range:
                        seq = ''.join(seq_lines)
                        write_record(out_f, header, seq, 'head', dict_head_end_range[header][0], dict_head_end_range[header][1])
                    if header in dict_tail_end_range:
                        seq = ''.join(seq_lines)
                        write_record(out_f, header, seq, 'tail', dict_tail_end_range[header][0], dict_tail_end_range[header][1])
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        # Write last record
        if header in dict_head_end_range:
            seq = ''.join(seq_lines)
            write_record(out_f, header, seq, 'head', dict_head_end_range[header][0], dict_head_end_range[header][1])
        if header in dict_tail_end_range:
            seq = ''.join(seq_lines)
            write_record(out_f, header, seq, 'tail', dict_tail_end_range[header][0], dict_tail_end_range[header][1])


def align_terminals_to_reference(ends_fasta, reference_fasta, output_paf, threads, minimap2_path="minimap2"):
    """
    Runs minimap2 -x asm5 to align ends_fasta to reference_fasta, writes output to output_paf.
    """
    cmd = [minimap2_path, "-x", "asm5", "-t", str(threads), reference_fasta, ends_fasta]
    with open(output_paf, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


def parse_ends_paf(paf_file):
    """
    Reads the PAF file and extracts relevant alignments using PAFReader.
    Returns a dict of ContigAlignment objects.
    """
    return read_genome_alignments(paf_file)


def process_terminals(input_fasta, reference_fasta, full_paf_path, dict_head_end_range, dict_tail_end_range, minimap2_path="minimap2", threads=1):
    """
    Orchestrates the extraction, alignment, and parsing of terminal sequences.
    """
    ends_fasta = input_fasta + ".ends.fa"
    ends_paf = full_paf_path + ".ends.paf"
    write_terminal_fasta(input_fasta, ends_fasta, dict_head_end_range, dict_tail_end_range)
    align_terminals_to_reference(ends_fasta, reference_fasta, ends_paf, threads, minimap2_path)
    return parse_ends_paf(ends_paf)

def map_terminal_to_full(alns_selected_full, alns_terminal):
    """
    For each query in alns_selected_full, find the closest contig_name_head and contig_name_tail terminal alignments
    to the start and end of the full alignment (on reference coordinates), and report their residue_matches.
    Returns a dict mapping query_name to (head_residue_matches, tail_residue_matches).
    """
    terminal_matches = {}
    for query_name, full_aln in alns_selected_full.items():
        # For each full alignment, find the closest head and tail terminal alignments
        best_head_dist = float('inf')
        best_head_residue = 0
        best_tail_dist = float('inf')
        best_tail_residue = 0
        # For each alignment in the full alignment
        for i in range(full_aln.num_alns):
            full_ref_start = full_aln.ref_starts[i]
            full_ref_end = full_aln.ref_ends[i]
            # Find all head and tail terminal alignments for this query
            head_key = f"{query_name}_head"
            tail_key = f"{query_name}_tail"
            # Check head
            if head_key in alns_terminal:
                term_aln = alns_terminal[head_key]
                # Assume only one alignment per terminal
                term_ref_start = term_aln.ref_starts[0]
                dist = abs(full_ref_start - term_ref_start)
                if dist < best_head_dist:
                    best_head_dist = dist
                    best_head_residue = term_aln.residue_matches[0]
            # Check tail
            if tail_key in alns_terminal:
                term_aln = alns_terminal[tail_key]
                term_ref_end = term_aln.ref_ends[0]
                dist = abs(full_ref_end - term_ref_end)
                if dist < best_tail_dist:
                    best_tail_dist = dist
                    best_tail_residue = term_aln.residue_matches[0]
        terminal_matches[query_name] = (best_head_residue, best_tail_residue)
    return terminal_matches


def filter_boundary_alignments(selected_full, boundary_ratio=0.05):
    """
    Returns a dict with only those alignments from selected_full that overlap the left or right boundary
    of the query sequence (within boundary_ratio of the length).
    """
    selected_boundary = {}
    for query_name, full_aln in selected_full.items():
        # Find indices that overlap the boundary
        boundary_indices = [
            i for i in range(full_aln.num_alns)
            if full_aln.query_starts[i] < full_aln.query_len * boundary_ratio
            or full_aln.query_ends[i] > full_aln.query_len * (1 - boundary_ratio)
        ]
        if boundary_indices:  # Only assign if there is at least one
            selected_boundary[query_name] = ContigAlignment(
                full_aln.query_header,
                full_aln.query_len,
                [full_aln.query_starts[i] for i in boundary_indices],
                [full_aln.query_ends[i] for i in boundary_indices],
                [full_aln.strands[i] for i in boundary_indices],
                [full_aln.ref_headers[i] for i in boundary_indices],
                [full_aln.ref_lens[i] for i in boundary_indices],
                [full_aln.ref_starts[i] for i in boundary_indices],
                [full_aln.ref_ends[i] for i in boundary_indices],
                [full_aln.residue_matches[i] for i in boundary_indices],
                [full_aln.aln_lens[i] for i in boundary_indices],
                [full_aln.mapqs[i] for i in boundary_indices]
            )
    return selected_boundary


def main():
    parser = argparse.ArgumentParser(description='Combine, filter full-length and terminal PAF alignments and build scaffold graph')
    parser.add_argument('-fl', '--full_length', required=True, help='Full length alignment PAF file')
    parser.add_argument('-f', '--fasta', required=True, help='Components FASTA file')
    parser.add_argument('-r', '--reference', required=True, help='Reference FASTA file for terminal alignment')
    parser.add_argument('-o', '--output', help='Output scaffold graph report file')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads [1]')
    parser.add_argument('--terminal_ratio', type=float, default=0.05, help='Terminal ratio [0.05]')
    parser.add_argument('--map_length', type=int, default=5000, help='Map length [5000]')
    parser.add_argument('--min_contig_length', type=int, default=10000, help='Minimum contig length [10000]')
    parser.add_argument('--max_min_id_ratio', type=float, default=0.5, help='Maximum minimum identity ratio [0.7]')
    parser.add_argument('--debug', action='store_true', help='Debug mode')
    parser.add_argument('--info_only', action='store_true', help='Output only the agp and info files without the patched fasta')
    args = parser.parse_args()

    full_len_paf = args.full_length
    output_path = args.output
    components_fn = args.fasta
    reference_fn = args.reference
    threads = args.threads
    TERMINAL_RATIO = args.terminal_ratio
    MAP_LENGTH = args.map_length
    MIN_CONTIG_LENGTH = args.min_contig_length
    MAX_MIN_ID_RATIO = args.max_min_id_ratio
    DEBUG = args.debug
    INFO_ONLY = args.info_only
    
    alns_full_len = read_genome_alignments(full_len_paf, True)

    alns_selected_full = filter_high_quality_alns(alns_full_len, minimum_length=MIN_CONTIG_LENGTH, length_ratio=TERMINAL_RATIO, max_min_id_ratio=MAX_MIN_ID_RATIO)
    dict_head_end_range, dict_tail_end_range = select_terminals(alns_selected_full, terminal_ratio=TERMINAL_RATIO, map_lenth=MAP_LENGTH)
    # Generate, align, and parse terminal sequences
    alns_terminal = process_terminals(components_fn, reference_fn, full_len_paf, dict_head_end_range, dict_tail_end_range, threads=threads)
    terminal_matches = map_terminal_to_full(alns_selected_full, alns_terminal)

    if DEBUG:
        print("selected_full", "--------------------------------")
        for key, item in alns_selected_full.items():
            print(key, item.query_header, item.query_len, item.num_alns, item.ref_starts, item.ref_ends, item.strands)
        print("terminal_matches", "--------------------------------")
        for key, item in terminal_matches.items():
            print(key, item)
    
    ref_alns = reverse_alignments(alns_selected_full)
    aln_psg = build_aln_scaffold_graph(ref_alns, components_fn, TERMINAL_RATIO, terminal_matches, debug=DEBUG)
    match_psg = aln_psg.max_weight_matching()
    match_psg.write_agp(output_path+'.agp', components_fn, add_suffix_to_unplaced=False)
    if not INFO_ONLY:
        cmd = [
            "./ragtag_agp2fa.py",
            output_path + ".agp",
            components_fn
        ]
        run_oae(cmd, output_path + ".fasta", "test.log")

    output_graph_info(match_psg, output_path + ".info")





if __name__ == "__main__":
    main()
