#!/usr/bin/env python3
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

Modified by Mao-Jan Lin in 2025 <mj.maojanlin@gmail.com>
"""

from ragtag_utilities.utilities import log, run_oae, get_ragtag_version


import csv
import os
import sys
import argparse
from typing import List, Dict, Tuple

def read_agp_file(file_path: str) -> List[Dict]:
    """Read an AGP file and return a list of dictionaries containing the AGP entries."""
    agp_entries = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                entry = {
                    'object': fields[0],
                    'object_beg': int(fields[1]),
                    'object_end': int(fields[2]),
                    'part_num': int(fields[3]),
                    'component_type': fields[4],
                    'component_id': fields[5],
                    'component_beg': int(fields[6]),
                    'component_end': int(fields[7]),
                    'orientation': fields[8]
                }
                agp_entries.append(entry)
    return agp_entries

def reverse_segments(segments: List[Dict]) -> List[Dict]:
    """Reverse the segments."""
    reversed_segments = []
    max_position = max(segment['object_end'] for segment in segments)
    for segment in segments[::-1]:
        reversed_segments.append({
            'object': segment['object'],
            'object_beg': max_position - segment['object_end'] + 1,
            'object_end': max_position - segment['object_beg'] + 1,
            'orientation': '+' if segment['orientation'] == '-' else '-',
            'part_num': segment['part_num'],
            'component_type': segment['component_type'],
            'component_id': segment['component_id'],
            'component_beg': segment['component_beg'],
            'component_end': segment['component_end']
        })
    return reversed_segments

def find_connection_segments(agp_entries: List[Dict], start_node: str, end_node: str) -> List[Dict]:
    """Find all segments between start_node and end_node in the AGP entries."""
    segments = []
    found_start = False
    found_end = False
    
    for entry in agp_entries:
        if entry['component_id'] == start_node:
            found_start = True
            segments.append(entry)
            if found_end:
                flag_forward = False
                break
        elif entry['component_id'] == end_node:
            found_end = True
            segments.append(entry)
            if found_start:
                flag_forward = True
                break
        elif found_start or found_end:
            segments.append(entry)
    
    if flag_forward:
        return segments
    else:
        return reverse_segments(segments)

def adjust_coordinates(segments: List[Dict], start_pos: int = 1, start_num: int = 1) -> List[Dict]:
    """Adjust the coordinates of segments to be continuous."""
    current_pos = start_pos
    adjusted_segments = []
    
    for idx, segment in enumerate(segments):
        segment_length = segment['object_end'] - segment['object_beg'] + 1
        adjusted_segment = segment.copy()
        adjusted_segment['part_num'] = start_num + idx
        adjusted_segment['object_beg'] = current_pos
        adjusted_segment['object_end'] = current_pos + segment_length - 1
        adjusted_segments.append(adjusted_segment)
        current_pos += segment_length
    
    return adjusted_segments

def merge_segments(core_segments: List[Dict], connection_segments: List[Dict]) -> List[Dict]:
    """Merge core segments with connection segments."""
    if not core_segments:
        offset = connection_segments[0]['object_beg'] - 1
        offset_segments = []
        for idx, segment in enumerate(connection_segments):
            offset_segments.append({
            'object': segment['object'],
            'object_beg': segment['object_beg'] - offset,
            'object_end': segment['object_end'] - offset,
            'part_num': idx + 1,
            'component_type': segment['component_type'],
            'component_id': segment['component_id'],
            'component_beg': segment['component_beg'],
            'component_end': segment['component_end'],
            'orientation': segment['orientation']
        })
        return offset_segments

    #print("Core segment:\n", core_segments[-1])
    
    start_node = connection_segments[0]['component_id']
    if core_segments[-1]['component_id'] == start_node:
        if core_segments[-1]['orientation'] != connection_segments[0]['orientation']:
            print("WARNING: Core segment and connection segment have different orientations")
            return False
        elif core_segments[-1]['component_beg'] > connection_segments[0]['component_end']:
            print(f"WARNING: Core segment begins after connection segment, core_beg, connection_end: {core_segments[-1]['component_beg'], connection_segments[0]['component_end']}")
            if core_segments[-1]['component_beg'] - connection_segments[0]['component_end'] > connection_segments[1]['component_end'] - connection_segments[0]['component_beg']:
                return False
            else:
                if core_segments[-2]['orientation'] == '+':
                    element_beg = core_segments[-2]['component_beg']
                    element_end = core_segments[-2]['component_end'] - (core_segments[-1]['component_beg'] - connection_segments[0]['component_end'])
                else:
                    element_beg = core_segments[-2]['component_beg'] + (core_segments[-1]['component_beg'] - connection_segments[0]['component_end'])
                    element_end = core_segments[-2]['component_end']
                # try to connect to the connection_segments[1:]
                connected_segment = {
                    'object': core_segments[-2]['object'],
                    'object_beg': core_segments[-2]['object_beg'],
                    'object_end': core_segments[-2]['object_end'] + element_end - element_beg,
                    'part_num': core_segments[-2]['part_num'],
                    'component_type': core_segments[-2]['component_type'],
                    'component_id': core_segments[-2]['component_id'],
                    'component_beg': element_beg,
                    'component_end': element_end,
                    'orientation': core_segments[-2]['orientation']
                }
                return core_segments[:-2] + [connected_segment] + adjust_coordinates(connection_segments[1:], connected_segment['object_end'] + 1, connected_segment['part_num'] + 1)
        else:
            if core_segments[-1]['orientation'] == '+':
                element_beg = core_segments[-1]['component_beg']
                element_end = connection_segments[0]['component_end']
            else:
                element_beg = connection_segments[0]['component_beg']
                element_end = core_segments[-1]['component_end']
            connected_segment = {
                'object': core_segments[-1]['object'],
                'object_beg': core_segments[-1]['object_beg'],
                'object_end': core_segments[-1]['object_beg'] + element_end - element_beg,
                'part_num': core_segments[-1]['part_num'],
                'component_type': core_segments[-1]['component_type'],
                'component_id': connection_segments[0]['component_id'],
                'component_beg': element_beg,
                'component_end': element_end,
                'orientation': core_segments[-1]['orientation']
            }
            return core_segments[:-1] + [connected_segment] + adjust_coordinates(connection_segments[1:], connected_segment['object_end'] + 1, connected_segment['part_num'] + 1)
    else:
        print(f"WARNING: Core segment does not end with start node, last_core_node, start_node: {core_segments[-1]['component_id']}, {start_node}")
        return False


def process_connections(input_path: str, output_path: str, connections_file: str):
    """Process all connections and generate combined AGP output."""
    all_segments = []
    current_pos = 1
    
    with open(connections_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) != 3:
                continue
                
            start_node, end_node, agp_file = row
            agp_path = os.path.join(input_path, f"{agp_file}.agp")
            print(f"Processing {agp_path}")
            
            if not os.path.exists(agp_path):
                print(f"Warning: AGP file not found: {agp_path}")
                continue
            
            agp_entries = read_agp_file(agp_path)
            segments = find_connection_segments(agp_entries, start_node[:-2], end_node[:-2])
            #print('--------------------------------'*3)
            #print(segments)
            
            if segments:
                #adjusted_segments = adjust_coordinates(segments, current_pos)
                #all_segments.extend(adjusted_segments)
                #current_pos = adjusted_segments[-1]['object_end'] + 1
                updated_segments = merge_segments(all_segments, segments)
                if updated_segments:
                    all_segments = updated_segments
                else:
                    print("Failed to merge segments")
                    break
    
    # Write the combined AGP output
    with open(output_path, 'w') as f:
        for segment in all_segments:
            f.write(f"{segment['object']}\t{segment['object_beg']}\t{segment['object_end']}\t"
                   f"{segment['part_num']}\t{segment['component_type']}\t{segment['component_id']}\t"
                   f"{segment['component_beg']}\t{segment['component_end']}\t{segment['orientation']}\n")
    


def main():
    parser = argparse.ArgumentParser(description='Combine AGP files based on connection information.')
    parser.add_argument('--input_path', help='Directory containing the input AGP files')
    parser.add_argument('--output_path', help='Path for the output combined AGP file')
    parser.add_argument('--connections_file', help='Path to the CSV file containing the connections')
    parser.add_argument('--fasta', help='Path to the FASTA file')
    args = parser.parse_args()
    
    
    process_connections(args.input_path, args.output_path+'.agp', args.connections_file)

    if args.fasta:
        cmd = [
            "./ragtag_agp2fa.py",
            args.output_path + ".agp",
            args.fasta
        ]
        run_oae(cmd, args.output_path + ".fasta", "test.log")

if __name__ == "__main__":
    main()
