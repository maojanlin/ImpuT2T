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


def get_connectable_fragments(ctg_alns, components_fn, terminal_matches, TERMINAL_RATIO=0.05, debug=False, blacklist_pairs=None):
    """
    Build a directed scaffold graph from filtered alignments
    """
    def check_connectable(idx_l, idx_r, als, dict_ref_end_connectable):
        """
        idx_l start should be on the left of idx_r start
        return (True, type)
        - type 0: ends not connectable
        - type 1: two ends contain
        - type 2: one end contain
        - type 3: reverse one end contain
        - type 4: overlap
        - type 5: disjoint
        5, 4 are good, 2, 3 is acceptable, 1 and 0 are negligible
        """
        # Check bounds to avoid index errors
        if idx_l >= len(dict_ref_end_connectable['L_end']) or idx_r >= len(dict_ref_end_connectable['R_end']):
            return False, 0
            
        if idx_l < idx_r:
            #if dict_ref_end_connectable['R_end'][idx_l] and dict_ref_end_connectable['L_end'][idx_r]:
            if als.query_ends[idx_l] <= als.query_starts[idx_r]:
                return True, 5
            elif als.query_ends[idx_l] <= als.query_ends[idx_r]:
                return True, 4
            else: # contain
                if dict_ref_end_connectable['R_end'][idx_r]:
                    return False, 1
                else:
                    return True, 2
        if idx_l > idx_r:
            if als.query_ends[idx_l] <= als.query_starts[idx_r]:
                if dict_ref_end_connectable['L_end'][idx_l]:
                    return False, 1
                else:
                    return True, 2
        return False, 0

    # notice that the query_seq is the reference name, and the ref_seq is the query contig name
    dict_all_connections = dict() # key: query_seq, value: list of connections
    for query_seq in ctg_alns:
        als = ctg_alns[query_seq]
        als.sort_by_query()
        last = None
        last_reversed = False

        # dict_ref_end_connectable:
        #   L_end: list of booleans, True if the reference end is connectable to the left end of the reference
        #   R_end: list of booleans, True if the reference end is connectable to the right end of the reference
        #   checking if a contig is split into multiple segments, if so, only boundary segments are connectable for the corresponding ends
        #   NEW: segments must also overlap the terminal regions (TERMINAL_RATIO) of the contig
        dict_ref_end_connectable = {'L_end': [], 'R_end': []}
        for i in range(als.num_alns):
            if debug:
                print(i, als.ref_headers[i], '(', als.ref_starts[i], als.ref_ends[i],')')
            duplicate_segments = []
            for j in range(als.num_alns):
                if i != j and als.ref_headers[i] == als.ref_headers[j]:
                    duplicate_segments.append(j)
            
            # Check if segment overlaps terminal regions
            if als.strands[i] == '+':

                overlaps_left_terminal = als.ref_starts[i] < als.ref_lens[i] * TERMINAL_RATIO
                if terminal_matches[als.ref_headers[i]]['head'].get((query_seq, als.query_starts[i])) == None:
                    terminal_matches[als.ref_headers[i]]['head'][(query_seq, als.query_starts[i])] = 0
                overlaps_right_terminal = als.ref_ends[i] > als.ref_lens[i] * (1 - TERMINAL_RATIO)
                if terminal_matches[als.ref_headers[i]]['tail'].get((query_seq, als.query_ends[i])) == None:
                    terminal_matches[als.ref_headers[i]]['tail'][(query_seq, als.query_ends[i])] = 0
            else:
                overlaps_right_terminal = als.ref_starts[i] < als.ref_lens[i] * (TERMINAL_RATIO)
                if terminal_matches[als.ref_headers[i]]['head'].get((query_seq, als.query_ends[i])) == None:
                    terminal_matches[als.ref_headers[i]]['head'][(query_seq, als.query_ends[i])] = 0
                overlaps_left_terminal = als.ref_ends[i] > als.ref_lens[i] * (1 - TERMINAL_RATIO)
                if terminal_matches[als.ref_headers[i]]['tail'].get((query_seq, als.query_starts[i])) == None:
                    terminal_matches[als.ref_headers[i]]['tail'][(query_seq, als.query_starts[i])] = 0
                    overlaps_left_terminal = als.ref_ends[i] > als.ref_lens[i] * (1 - TERMINAL_RATIO)
            
            if len(duplicate_segments) > 0:
                # multiple segments, can only be connected to one end
                dict_ref_end_connectable['L_end'].append(False)
                dict_ref_end_connectable['R_end'].append(False)

                if als.ref_starts[i] <= min(als.ref_starts[ele] for ele in duplicate_segments):
                    if als.strands[i] == '+':
                        # Only mark as connectable if it overlaps the left terminal region
                        dict_ref_end_connectable['L_end'][i] = overlaps_left_terminal
                    else:
                        # Only mark as connectable if it overlaps the right terminal region
                        dict_ref_end_connectable['R_end'][i] = overlaps_right_terminal
                if als.ref_ends[i] >= max(als.ref_ends[ele] for ele in duplicate_segments):
                    if als.strands[i] == '+':
                        # Only mark as connectable if it overlaps the right terminal region
                        dict_ref_end_connectable['R_end'][i] = overlaps_right_terminal
                    else:
                        # Only mark as connectable if it overlaps the left terminal region
                        dict_ref_end_connectable['L_end'][i] = overlaps_left_terminal
            else:
                # single segment, can be connected to both ends only if they overlap terminal regions
                dict_ref_end_connectable['L_end'].append(overlaps_left_terminal)
                dict_ref_end_connectable['R_end'].append(overlaps_right_terminal)
        if debug:
            print("L_end:", dict_ref_end_connectable['L_end'])
            print("R_end:", dict_ref_end_connectable['R_end'])
            for i in range(als.num_alns):
                print(i, als.ref_headers[i], dict_ref_end_connectable['L_end'][i], dict_ref_end_connectable['R_end'][i], als.ref_starts[i], als.ref_ends[i], '\t|\t', als.strands[i], als.query_starts[i], als.query_ends[i])
        
        list_connectable = []
        dict_headers_to_idx = {}
        for i in range(als.num_alns):
            if dict_ref_end_connectable['L_end'][i] or dict_ref_end_connectable['R_end'][i]:
                list_connectable.append(i)
                if als.ref_headers[i] not in dict_headers_to_idx:
                    dict_headers_to_idx[als.ref_headers[i]] = [i]
                else:
                    dict_headers_to_idx[als.ref_headers[i]].append(i)

        ##################################################################################################################
        # merge overlapping splitted segments
        ##################################################################################################################
        merged_segments = []  # New list to store merged segments
        segments_to_remove = set()  # Track which segments should be removed
        
        for header, list_idx in dict_headers_to_idx.items():
            if len(list_idx) == 1:
                merged_segments.append(list_idx[0])
                continue
            
            # Ensure we only have 2 segments as expected
            if len(list_idx) > 2:
                if debug:
                    print(f"Warning: More than 2 segments for {header}, skipping merge")
                merged_segments.extend(list_idx)
                continue
            
            # Check if segments can be merged
            idx_0, idx_1 = list_idx[0], list_idx[1]
            
            # Only merge if strands are the same
            if als.strands[idx_0] != als.strands[idx_1]:
                if debug:
                    print(f"Cannot merge segments for {header}: different strands ({als.strands[idx_0]} vs {als.strands[idx_1]})")
                merged_segments.extend(list_idx)
                continue
            
            # Check if segments overlap, query_ends and query_starts are the real position on the reference/pangenome
            if als.query_ends[idx_0] >= als.query_starts[idx_1]:
                # Segments overlap, can merge
                if debug:
                    print(f"Merging overlapping segments for {header}: {idx_0} and {idx_1}")
                    print(f"  Ref positions: ({als.ref_starts[idx_0]}, {als.ref_ends[idx_0]}) and ({als.ref_starts[idx_1]}, {als.ref_ends[idx_1]})")
                
                # Find which segment has L_end and which has R_end
                l_end_idx = None
                r_end_idx = None
                
                if dict_ref_end_connectable['L_end'][idx_0]:
                    l_end_idx = idx_0
                    r_end_idx = idx_1
                    #elif dict_ref_end_connectable['L_end'][idx_1]:
                    #    l_end_idx = idx_1
                    #    r_end_idx = idx_0
                else:
                    # Neither has L_end, skip merge
                    if debug:
                        print(f"Cannot merge segments for {header}: neither segment has L_end")
                        print(f"left side segment not having L_end avialable: {dict_ref_end_connectable['L_end']}")
                    merged_segments.extend(list_idx)
                    continue
                
                # Update the L_end segment to have both L_end and R_end = True
                dict_ref_end_connectable['L_end'][l_end_idx] = True
                dict_ref_end_connectable['R_end'][l_end_idx] = True
                
                # Update the alignment information for the merged segment
                # Use start positions from L_end segment and end positions from R_end segment
                als.ref_starts[l_end_idx] = min(als.ref_starts[idx_0], als.ref_starts[idx_1])
                als.ref_ends[l_end_idx] = max(als.ref_ends[idx_0], als.ref_ends[idx_1])
                als.query_starts[l_end_idx] = als.query_starts[l_end_idx] #min(als.query_starts[idx_0], als.query_starts[idx_1])
                als.query_ends[l_end_idx] = als.query_ends[r_end_idx] #max(als.query_ends[idx_0], als.query_ends[idx_1])
                
                # Add the merged segment to our new list
                merged_segments.append(l_end_idx)
                
                # Mark the R_end segment for removal
                segments_to_remove.add(r_end_idx)
                
                if debug:
                    print(f"  Merged segment {l_end_idx}: ref({als.ref_starts[l_end_idx]}, {als.ref_ends[l_end_idx]}), query({als.query_starts[l_end_idx]}, {als.query_ends[l_end_idx]})")
                    print(f"  Updated connectability: L_end[{l_end_idx}]={dict_ref_end_connectable['L_end'][l_end_idx]}, R_end[{l_end_idx}]={dict_ref_end_connectable['R_end'][l_end_idx]}")
            else:
                # No overlap, no need to merge - add both segments as-is without any modification
                if debug:
                    print(f"No overlap for segments in {header}: {idx_0} and {idx_1}")
                merged_segments.extend(list_idx)
        
        # Update list_connectable to use the merged segments
        list_connectable = sorted([idx for idx in merged_segments if idx not in segments_to_remove])
        
        if debug:
            #print(f"Original segments: {len(als.num_alns)}")
            print(f"Merged segments: {len(list_connectable)}")
            print(f"Segments removed: {len(segments_to_remove)}")
            print(f"Final connectable list: {list_connectable}")
            for idx in list_connectable:
                print(idx, als.ref_headers[idx], dict_ref_end_connectable['L_end'][idx], dict_ref_end_connectable['R_end'][idx], als.ref_starts[idx], als.ref_ends[idx], '\t|\t', als.strands[idx], als.query_starts[idx], als.query_ends[idx])
            print("--------------------------------"*3)
            
        ##################################################################################################################
        # Filter complete contained segments
        ##################################################################################################################
        right_most_position = 0
        left_of_right_most_position = 0
        contained_segments = []
        for idx in list_connectable:
            if als.query_ends[idx] > right_most_position:
                right_most_position = als.query_ends[idx]
                left_of_right_most_position = als.query_starts[idx]
            elif dict_ref_end_connectable['L_end'][idx] and dict_ref_end_connectable['R_end'][idx]:
                #print("contained", idx, als.ref_headers[idx], als.ref_starts[idx], right_most_position)
                contained_segments.append(idx)
            else:
                # split case, check if the other segment is also contained
                list_index = dict_headers_to_idx[als.ref_headers[idx]]
                if len(list_index) == 2:
                    if list_index[0] == idx:  # First segment
                        if als.query_ends[list_index[1]] < right_most_position:
                            contained_segments.append(idx)
                            if debug:
                                print("contained1", idx, als.ref_headers[idx], als.ref_starts[idx], right_most_position)
                    elif list_index[1] == idx:  # Second segment
                        if als.query_starts[list_index[0]] < right_most_position and als.query_starts[list_index[0]] > left_of_right_most_position:
                            contained_segments.append(idx)
                            if debug:
                                print("contained2", idx, als.ref_headers[idx], als.ref_starts[idx], right_most_position)
                else:
                    contained_segments.append(idx)
        if debug:
            print("List connectable:", list_connectable)
            for idx in list_connectable:
                if idx not in contained_segments:
                    print(idx, als.ref_headers[idx], dict_ref_end_connectable['L_end'][idx], dict_ref_end_connectable['R_end'][idx], als.ref_starts[idx], als.ref_ends[idx], '\t|\t', als.strands[idx], als.query_starts[idx], als.query_ends[idx])


        ##################################################################################################################
        # Build candidate connections
        ##################################################################################################################
        list_connectable_left = []
        list_connectable_right = []
        for idx in list_connectable:
            if idx not in contained_segments:
                if dict_ref_end_connectable['L_end'][idx]:
                    list_connectable_left.append(idx)
                if dict_ref_end_connectable['R_end'][idx]:
                    list_connectable_right.append(idx)

        list_candidate_connections = []
        for idx_l in list_connectable_right:
            break_mapping_position = -1
            for idx_r in list_connectable_left:
                #if idx_l >= idx_r:
                #    continue
                if break_mapping_position != -1 and als.query_starts[idx_r] > break_mapping_position:
                    break
                if als.ref_headers[idx_l] != als.ref_headers[idx_r]:
                    # Check for overlap before determining connectability
                    
                    is_connectable, connection_type = check_connectable(idx_l, idx_r, als, dict_ref_end_connectable)
                    if is_connectable:
                        my_query_end_offset = als.ref_lens[idx_l] - als.ref_ends[idx_l]
                        if als.strands[idx_l] == '-':
                            my_query_end_offset = als.ref_starts[idx_l]

                        their_query_start_offset = als.ref_starts[idx_r]
                        if als.strands[idx_r] == '-':
                            their_query_start_offset = als.ref_lens[idx_r] - als.ref_ends[idx_r]

                        my_query_end = als.query_ends[idx_l] + my_query_end_offset
                        their_query_start = als.query_starts[idx_r] - their_query_start_offset
                        overlap = my_query_end - their_query_start

                        # Check if overlap is acceptable
                        if overlap > als.ref_lens[idx_l] or overlap > als.ref_lens[idx_r]:
                            if debug:
                                print(f"Overlap too large for connection {idx_l} -> {idx_r}: {overlap} > min({als.ref_lens[idx_l]}, {als.ref_lens[idx_r]})")
                                print(f"  Query positions: {my_query_end} - {their_query_start} = {overlap}")
                                print(f"  Ref lengths: {als.ref_lens[idx_l]}, {als.ref_lens[idx_r]}")
                            continue  # Skip this connection due to excessive overlap
                        # Skip blacklisted contig pairs and continue searching
                        if blacklist_pairs is not None:
                            left_ctg = als.ref_headers[idx_l]
                            right_ctg = als.ref_headers[idx_r]
                            if (left_ctg, right_ctg) in blacklist_pairs:
                                if debug:
                                    print(f"Skipping blacklisted pair: {left_ctg} -> {right_ctg}")
                                continue
                            # also respect unordered set style if provided
                            if frozenset((left_ctg, right_ctg)) in blacklist_pairs:
                                if debug:
                                    print(f"Skipping blacklisted pair: {{ {left_ctg}, {right_ctg} }}")
                                continue

                        list_candidate_connections.append((idx_l, idx_r, connection_type))
                        if connection_type >= 4:
                            # Only early-break if we accepted a non-blacklisted strong connection
                            break_mapping_position = als.query_starts[idx_r] + 100000
                        #    break
        
        if debug:
            print("================================", als.query_header, "================================")
            for idx_l, idx_r, connection_type in list_candidate_connections:
                if connection_type >= 4:
                    print("Connectable", idx_l, idx_r,  als.ref_headers[idx_l], '->', als.ref_headers[idx_r])
                elif connection_type == 2:
                    print("Right Contain", idx_l, idx_r, als.ref_headers[idx_l], '->', als.ref_headers[idx_r])
                elif connection_type == 3:
                    print("Left Contain", idx_l, idx_r, als.ref_headers[idx_r], '->', als.ref_headers[idx_l])
                else:
                    print("Not Connectable", idx_l, idx_r, als.ref_headers[idx_l], '???', als.ref_headers[idx_r])

        dict_all_connections[query_seq] = list_candidate_connections

    return dict_all_connections


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

        pair_left_pos = [(idx, 0, als.query_starts[idx]) for idx in range(als.num_alns)]
        pair_right_pos = [(idx, 1, als.query_ends[idx]) for idx in range(als.num_alns) ]
        all_pairs = pair_left_pos + pair_right_pos

        all_pairs.sort(key=lambda x: x[2])

        set_contained = set()
        largest_right_pos = 0
        if debug:
            print("--------------------------------"*3)
        for i in range(als.num_alns):
            
            if als.query_ends[i] < largest_right_pos:
                if debug:
                    print("Find contained", i, als.ref_headers[i], als.query_ends[i], largest_right_pos)
                set_contained.add(i)
            if dict_ref_end_connectable['L_end'][i] and dict_ref_end_connectable['R_end'][i]:
                largest_right_pos = max(largest_right_pos, als.query_ends[i])
        if debug:
            print("--------------------------------"*3)

        # Iterate over each alignment for this query sequence
        #for i in list_traverse:
        for i in range(als.num_alns):
            # Don't connect irrelevant middle alignments
            if not dict_ref_end_connectable['L_end'][i] and not dict_ref_end_connectable['R_end'][i]:
                continue
            if i in set_contained:
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
                                if debug:
                                    print("last_reversed", pair_score)
                            else:
                                pair_score += terminal_matches[als.ref_headers[last]][1]
                                if debug:
                                    print("last_forward", pair_score)
                            if cur_reversed:
                                pair_score += terminal_matches[als.ref_headers[i]][1]
                                if debug:
                                    print("cur_reversed", terminal_matches[als.ref_headers[i]][1])
                            else:
                                pair_score += terminal_matches[als.ref_headers[i]][0]
                                if debug:
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
                            sg.add_edge(u, v, alignment, pair_score, -overlap)
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
        for i in range(full_aln.num_alns):
            segment_length = full_aln.query_ends[i]-full_aln.query_starts[i]
            if (full_aln.query_len*0.2 > segment_length or segment_length < 500000) and \
               full_aln.residue_matches[i]/full_aln.aln_lens[i] < max_identity*max_min_id_ratio: # 0.7 of max_identity
                continue
            if full_aln.query_len*length_ratio < segment_length: # 0.05 of query length
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
            f.write("\t".join([u,v, str(graph[u][v]["score"]), str(graph[u][v]["dist"]), "\n"]))


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
        extend_length = min(full_aln.query_len, int(map_lenth*0.5))
        quarter_length = min(full_aln.query_len, int(map_lenth*0.25))
        if head_best_idx is not None:
            #dict_head_end_range[query_name] = (full_aln.query_starts[head_best_idx], full_aln.query_starts[head_best_idx]+extend_length)
            #dict_head_end_range[query_name] = (full_aln.query_starts[head_best_idx], max(map_lenth, full_aln.query_starts[head_best_idx]+extend_length))
            dict_head_end_range[query_name] = (max(0, full_aln.query_starts[head_best_idx]-quarter_length), max(map_lenth, full_aln.query_starts[head_best_idx]+quarter_length))
        if tail_best_idx is not None:
            #dict_tail_end_range[query_name] = (full_aln.query_ends[tail_best_idx]-extend_length, full_aln.query_ends[tail_best_idx])
            #dict_tail_end_range[query_name] = (min(full_aln.query_len-map_lenth, full_aln.query_ends[tail_best_idx]-extend_length), full_aln.query_ends[tail_best_idx])
            dict_tail_end_range[query_name] = (min(full_aln.query_len-map_lenth, full_aln.query_ends[tail_best_idx]-quarter_length), min(full_aln.query_len, full_aln.query_ends[tail_best_idx]+quarter_length))
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
    #cmd = [minimap2_path, "-x", "asm5", "-t", str(threads), reference_fasta, ends_fasta]
    cmd = [minimap2_path, "-cx", "asm5", "-t", str(threads), reference_fasta, ends_fasta]
    with open(output_paf, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


def parse_ends_paf(paf_file):
    """
    Reads the PAF file and extracts relevant alignments using PAFReader.
    Returns a dict of ContigAlignment objects.
    """
    return read_genome_alignments(paf_file)


def process_terminals(input_fasta, reference_fasta, dict_head_end_range, dict_tail_end_range, output_path, minimap2_path="minimap2", threads=1):
    """
    Orchestrates the extraction, alignment, and parsing of terminal sequences.
    """
    ends_fasta = output_path + ".ends.fa"
    ends_paf = output_path + ".ends.paf"
    write_terminal_fasta(input_fasta, ends_fasta, dict_head_end_range, dict_tail_end_range)
    align_terminals_to_reference(ends_fasta, reference_fasta, ends_paf, threads, minimap2_path)
    return parse_ends_paf(ends_paf)

def map_terminal_to_full(alns_selected_full, alns_terminal):
    """
    For each query in alns_selected_full, find the best matching terminal alignments (head and tail)
    against the selected full-length alignment segments, constrained to the same reference (ref_header).
    
    Data structure:
    - Layer 1: contig_name -> 'head', 'tail'
    - Layer 2: (ref_name, ref_start, ref_end, strand) -> best_residue_matches
    
    Matching logic:
    - Forward strand (+): head covers full_start±5000, tail covers full_end±5000
    - Reverse strand (-): head covers full_end±5000, tail covers full_start±5000
    - Only consider terminals with matching ref_header
    - Keep best residue_matches per (ref_name, ref_start, ref_end, strand) combination
    
    Returns: dict with two-layer structure for terminal matches
    """
    terminal_matches = {}
    
    for query_name, full_aln in alns_selected_full.items():
        terminal_matches[query_name] = {"head": {}, "tail": {}}
        
        head_key = f"{query_name}_head"
        tail_key = f"{query_name}_tail"
        
        # Iterate over all full segments
        for i in range(full_aln.num_alns):
            full_ref_start = full_aln.ref_starts[i]
            full_ref_end = full_aln.ref_ends[i]
            full_ref_header = full_aln.ref_headers[i]
            full_strand = full_aln.strands[i]
            
            # For forward strand: head->start, tail->end
            # For reverse strand: head->end, tail->start
            if full_strand == '+':
                head_target = full_ref_start
                tail_target = full_ref_end
            else:  # reverse strand
                head_target = full_ref_end
                tail_target = full_ref_start
            
            # Check head terminal alignments
            if head_key in alns_terminal:
                term_aln = alns_terminal[head_key]
                current_best_dist = 10001
                for j in range(term_aln.num_alns):
                    if term_aln.ref_headers[j] != full_ref_header:
                        continue
                    
                    # Check if terminal covers target position ±5000
                    term_start = term_aln.ref_starts[j]
                    term_end = term_aln.ref_ends[j]
                    min_dist = min(abs(term_start - head_target), abs(term_end - head_target))
                    if min_dist <= 10000:
                        key = (full_ref_header, head_target)
                        if min_dist < current_best_dist:
                            current_best_dist = min_dist
                            terminal_matches[query_name]["head"][key] = term_aln.residue_matches[j]
            
            # Check tail terminal alignments
            if tail_key in alns_terminal:
                term_aln = alns_terminal[tail_key]
                current_best_dist = 10001
                for j in range(term_aln.num_alns):
                    if term_aln.ref_headers[j] != full_ref_header:
                        continue
                    
                    # Check if terminal covers target position ±5000
                    term_start = term_aln.ref_starts[j]
                    term_end = term_aln.ref_ends[j]
                    min_dist = min(abs(term_start - tail_target), abs(term_end - tail_target))
                    if min_dist <= 10000:
                        key = (full_ref_header, tail_target)
                        if min_dist < current_best_dist:
                            current_best_dist = min_dist
                            terminal_matches[query_name]["tail"][key] = term_aln.residue_matches[j]
    
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

def generate_relevant_fasta(components_fn, reference_fn, output_path):
    """
    Generate the relevant fasta file for patching.
    """
    cmd = f"cat {components_fn} {reference_fn} > {output_path}.relevant.fasta"
    subprocess.run(cmd, shell=True, check=True)


def main():
    parser = argparse.ArgumentParser(description='Combine, filter full-length and terminal PAF alignments and build scaffold graph')
    parser.add_argument('-fl', '--full_length', required=True, help='Full length alignment PAF file')
    parser.add_argument('-q', '--query', required=True, help='Query FASTA file')
    parser.add_argument('-r', '--reference', required=True, help='Reference FASTA file for terminal alignment')
    parser.add_argument('-o', '--output', help='Output scaffold graph report file')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads [1]')
    parser.add_argument('--terminal_ratio', type=float, default=0.05, help='Terminal ratio [0.05]')
    parser.add_argument('--map_length', type=int, default=5000, help='Map length [5000]')
    parser.add_argument('--min_contig_length', type=int, default=10000, help='Minimum contig length [10000]')
    parser.add_argument('--max_min_id_ratio', type=float, default=0.5, help='Maximum minimum identity ratio [0.7]')
    parser.add_argument('--debug', action='store_true', help='Debug mode')
    parser.add_argument('--info_only', action='store_true', help='Output only the agp and info files without the patched fasta')
    parser.add_argument('--connect_mode', choices=['closest', 'longest'], default='closest', help='Connection mode [closest]')
    parser.add_argument('--blacklist', help='Blacklist file')
    args = parser.parse_args()

    full_len_paf = args.full_length
    output_path = args.output
    components_fn = args.query
    reference_fn = args.reference
    threads = args.threads
    TERMINAL_RATIO = args.terminal_ratio
    MAP_LENGTH = args.map_length
    MIN_CONTIG_LENGTH = args.min_contig_length
    MAX_MIN_ID_RATIO = args.max_min_id_ratio
    DEBUG = args.debug
    INFO_ONLY = args.info_only
    CONNECT_MODE = args.connect_mode
    BLACKLIST = args.blacklist
    
    set_blacklist = set()
    if BLACKLIST:
        with open(BLACKLIST, "r") as f:
            for line in f:
                line = line.rstrip().split(',')
                set_blacklist.add((line[0], line[1]))
                set_blacklist.add((line[1], line[0]))
    else:
        set_blacklist = set()
        
    alns_full_len = read_genome_alignments(full_len_paf, True)

    alns_selected_full = filter_high_quality_alns(alns_full_len, minimum_length=MIN_CONTIG_LENGTH, length_ratio=TERMINAL_RATIO, max_min_id_ratio=MAX_MIN_ID_RATIO)
    dict_head_end_range, dict_tail_end_range = select_terminals(alns_selected_full, terminal_ratio=TERMINAL_RATIO, map_lenth=MAP_LENGTH)
    # Generate, align, and parse terminal sequences

    print("Processing terminals with minimap2 realignment...")
    alns_terminal = process_terminals(components_fn, reference_fn, dict_head_end_range, dict_tail_end_range, output_path, threads=threads)
    terminal_matches = map_terminal_to_full(alns_selected_full, alns_terminal)
    
    ref_alns = reverse_alignments(alns_selected_full)
    #aln_psg = build_aln_scaffold_graph(ref_alns, components_fn, TERMINAL_RATIO, terminal_matches, debug=DEBUG)
    if DEBUG:
        print("ref_alns", "---------------------------------------------------------------------------------------")
    dict_all_connections = get_connectable_fragments(ref_alns, components_fn, terminal_matches, TERMINAL_RATIO*2, debug=DEBUG, blacklist_pairs=set_blacklist)

    dict_node_connections = {}
    if CONNECT_MODE == 'closest':
        # sorted the connections by the two nodes
        for key, item in dict_all_connections.items():
            for connection in item:
                idx_l, idx_r, connection_type = connection

                als = ref_alns[key]
                
                my_query_end_offset = als.ref_lens[idx_l] - als.ref_ends[idx_l]
                if als.strands[idx_l] == '-':
                    my_query_end_offset = als.ref_starts[idx_l]

                their_query_start_offset = als.ref_starts[idx_r]
                if als.strands[idx_r] == '-':
                    their_query_start_offset = als.ref_lens[idx_r] - als.ref_ends[idx_r]

                my_query_end = als.query_ends[idx_l] + my_query_end_offset
                their_query_start = als.query_starts[idx_r] - their_query_start_offset
                overlap = my_query_end - their_query_start

                u = als.ref_headers[idx_l] + "_e"
                v = als.ref_headers[idx_r] + "_b"
                if als.strands[idx_l] == "-":
                    u = als.ref_headers[idx_l] + "_b"
                if als.strands[idx_r] == "-":
                    v = als.ref_headers[idx_r] + "_e"

                if u in dict_node_connections:
                    dict_node_connections[u].append((abs(overlap), key, idx_l, idx_r, connection_type))
                else:
                    dict_node_connections[u] = [(abs(overlap), key, idx_l, idx_r, connection_type)]
                if v in dict_node_connections:
                    dict_node_connections[v].append((abs(overlap), key, idx_l, idx_r, connection_type))
                else:
                    dict_node_connections[v] = [(abs(overlap), key, idx_l, idx_r, connection_type)]
        # filter the connections and only keep the one with the cloest distance
        removed_dict = {}
        for node, connections in dict_node_connections.items():
            connections.sort(key=lambda x: x[0])
            if connections[0][-1] >= 4:
                # remove others if the length is greater than 100,000 bp:
                base_length = connections[0][0]
                for connection in connections[1:]:
                    key, idx_l, idx_r, connection_type = connection[1:]
                    if connection[-1] >= 4 and connection[0] < base_length + 100000:
                        continue
                    else:
                        if key not in removed_dict:
                            removed_dict[key] = set()
                        removed_dict[key].add((idx_l, idx_r, connection_type))
                        if DEBUG:
                            print("REMOVED!", key, idx_l, idx_r, ref_alns[key].ref_headers[idx_l], ref_alns[key].ref_headers[idx_r], connection_type)
        for key, item in dict_all_connections.items():
            new_item = []
            if key not in removed_dict:
                continue
            for connection in item:
                if connection not in removed_dict[key]:
                    new_item.append(connection)
            dict_all_connections[key] = new_item

    sg = PatchScaffoldGraph(components_fn)
    if DEBUG:
        print("---------------------------------------------------------------------------------------")
    for key, item in dict_all_connections.items():
        if DEBUG:
            print("================================", key, "================================")
        for connection in item:
            idx_l, idx_r, connection_type = connection
            als = ref_alns[key]
            
            # Print connection information similar to what's in get_connectable_fragments()
            if DEBUG:
                if connection_type >= 4:
                    print("Connectable", idx_l, idx_r, als.ref_headers[idx_l], '->', als.ref_headers[idx_r])
                elif connection_type == 2:
                    print("Right Contain", idx_l, idx_r, als.ref_headers[idx_l], '->', als.ref_headers[idx_r])
                elif connection_type == 3:
                    print("Left Contain", idx_l, idx_r, als.ref_headers[idx_r], '->', als.ref_headers[idx_l])
                else:
                    print("Not Connectable", idx_l, idx_r, als.ref_headers[idx_l], '???', als.ref_headers[idx_r])
            
            # Print detailed alignment information
            #print(f"  Left segment {idx_l}: {als.ref_headers[idx_l]} ({als.ref_starts[idx_l]}-{als.ref_ends[idx_l]}) strand:{als.strands[idx_l]}")
            #print(f"  Right segment {idx_r}: {als.ref_headers[idx_r]} ({als.ref_starts[idx_r]}-{als.ref_ends[idx_r]}) strand:{als.strands[idx_r]}")
            #print(f"  Query positions: {als.query_starts[idx_l]}-{als.query_ends[idx_l]} -> {als.query_starts[idx_r]}-{als.query_ends[idx_r]}")
            #print(f"  Connection type: {connection_type}")
            #print()
            # Determine the scaffold graph nodes
            u = als.ref_headers[idx_l] + "_e"
            l_match = 'tail', (key, als.query_ends[idx_l])
            v = als.ref_headers[idx_r] + "_b"
            r_match = 'head', (key, als.query_starts[idx_r])
            if als.strands[idx_l] == "-":
                u = als.ref_headers[idx_l] + "_b"
                l_match = 'head', (key, als.query_ends[idx_l])
            if als.strands[idx_r] == "-":
                v = als.ref_headers[idx_r] + "_e"
                r_match = 'tail', (key, als.query_starts[idx_r])
            
            if DEBUG:
                print("key", key)
                print((als.ref_headers[idx_l], als.query_starts[idx_l], als.query_ends[idx_l]))
                print(l_match, terminal_matches[als.ref_headers[idx_l]][l_match[0]])
                print(r_match, terminal_matches[als.ref_headers[idx_r]][r_match[0]])
            left_pair_score = terminal_matches[als.ref_headers[idx_l]][l_match[0]][l_match[1]]
            right_pair_score = terminal_matches[als.ref_headers[idx_r]][r_match[0]][r_match[1]]
            pair_score = left_pair_score + right_pair_score
            if DEBUG:
                print("u, v, pair_score", u, v, pair_score)
            """
            pair_score = 0
            if als.strands[idx_l] == "-":
                pair_score += terminal_matches[als.ref_headers[idx_l]][0]
            else:
                pair_score += terminal_matches[als.ref_headers[idx_l]][1]
            if als.strands[idx_r] == "-":
                pair_score += terminal_matches[als.ref_headers[idx_r]][1]
            else:
                pair_score += terminal_matches[als.ref_headers[idx_r]][0]
            """

            # Calculate overlap (similar to the logic in get_connectable_fragments)
            my_query_end_offset = als.ref_lens[idx_l] - als.ref_ends[idx_l]
            if als.strands[idx_l] == '-':
                my_query_end_offset = als.ref_starts[idx_l]

            their_query_start_offset = als.ref_starts[idx_r]
            if als.strands[idx_r] == '-':
                their_query_start_offset = als.ref_lens[idx_r] - als.ref_ends[idx_r]

            my_query_end = als.query_ends[idx_l] + my_query_end_offset
            their_query_start = als.query_starts[idx_r] - their_query_start_offset
            overlap = my_query_end - their_query_start
            if DEBUG:
                print("my_query_end, their_query_start, distance", my_query_end, their_query_start, -overlap)

            # Now create the alignment with the correct variables
            alignment = Alignment(
                u,
                v,
                key,  # Use 'key' instead of 'query_seq' since that's what you have
                als.query_len,
                my_query_end,
                their_query_start,
                0,  # Always on the query's forward strand
                is_gap=False
            )
            sg.add_edge(u, v, alignment, pair_score, -overlap)
    if DEBUG:
        print("sg.edges", "---------------------------------------------------------------------------------------")
        for u, v in sg.edges:
            print(u, v, sg[u][v]["weight"], sg[u][v]["score"], sg[u][v]["dist"])
    
    output_graph_info(sg, output_path + ".info")
    sg.straight_to_agp(output_path+'.agp', components_fn, add_suffix_to_unplaced=False)
    
    match_psg = sg.max_weight_matching()
    if DEBUG:
        print("match_psg.edges", "---------------------------------------------------------------------------------------")
        for u, v in match_psg.edges:
            print(u, v, match_psg[u][v]["weight"], match_psg[u][v]["score"], match_psg[u][v]["dist"])
    match_psg.write_agp(output_path+'.linear.agp', components_fn, add_suffix_to_unplaced=False)
    #match_psg.write_agp('0825.test.agp', components_fn, add_suffix_to_unplaced=False)
    #output_graph_info(match_psg, output_path + ".info")
    
    print(f"Complete individual patching for {full_len_paf}!")
    exit()





if __name__ == "__main__":
    main()
