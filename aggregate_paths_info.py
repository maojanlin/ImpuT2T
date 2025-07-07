import argparse
import os
from collections import defaultdict
from typing import List, Dict
from combine_connections import read_agp_file, merge_segments
import subprocess
from ragtag_utilities.utilities import run_oae

def parse_sample_list(list_sample: str) -> list[str]:
    with open(list_sample) as f:
        return [line.strip() for line in f if line.strip()]

class EdgeAggregator:
    def __init__(self, input_path: str, samples: list[str]):
        self.input_path = input_path
        self.samples = samples
        self.edges = {}

    def aggregate(self):
        """
        edges:
        - key: (u, v) # sorted
        - value: (weight, max_score, max_sample_id)
        """
        for sample in self.samples:
            if os.path.exists(self.input_path + sample + ".info"):
                with open(self.input_path + sample + ".info", "r") as f:
                    for line in f.readlines():
                        u, v, score = line.strip().split("\t")
                        pair = tuple(sorted([u, v]))
                        if pair not in self.edges:
                            self.edges[pair] = [1, int(score), sample]
                        else:
                            if int(score) > self.edges[pair][1]:
                                self.edges[pair] = [self.edges[pair][0] + 1, int(score), sample]
                            else:
                                self.edges[pair][0] += 1
            else:
                print(f"Warning: {sample} not found in {self.input_path}")
        return self.edges

class GraphBuilder:
    def __init__(self, edges, weight_type, weight_ratio):
        self.edges = edges
        self.weight_type = weight_type
        self.weight_ratio = weight_ratio
        self.graph = defaultdict(list)
        self.weights = {}
        self.nodes = set()
        self.dict_node_score = defaultdict(list)
        self.dict_edge_info = {}

    def build(self):
        edges_sorted = sorted(self.edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
        for pair, values in edges_sorted:
            u, v = pair
            w = values[self.weight_type]
            self.dict_node_score[u].append((w, v))
            self.dict_node_score[v].append((w, u))
            self.dict_edge_info[(u, v)] = values
            self.dict_edge_info[(v, u)] = values

        filtered_edges = set()
        for node, scores in self.dict_node_score.items():
            max_score = max(scores, key=lambda x: x[0])[0]
            for score, neighbor in scores:
                if score < max_score * self.weight_ratio:
                    filtered_edges.add((node, neighbor))
                    filtered_edges.add((neighbor, node))

        for pair, values in edges_sorted:
            u, v = pair
            w = values[self.weight_type]
            if w < 6000:
                continue
            if (u, v) in filtered_edges or (v, u) in filtered_edges:
                continue
            # undirected graph assignment
            self.graph[u].append(v)
            self.graph[v].append(u)
            self.weights[(u, v)] = w
            self.weights[(v, u)] = w
            self.nodes.update([u, v])
            self.nodes.update([v, u])

        for node in self.nodes:
            if node.endswith("_b"):
                contig = node[:-2]
                e_node = f"{contig}_e"
                self.graph[node].append(e_node)
                self.weights[(node, e_node)] = 0
        return self.graph, self.weights, self.nodes, self.dict_edge_info

class PathFinder:
    def __init__(self, graph, weights, nodes):
        self.graph = graph
        self.weights = weights
        self.nodes = nodes

    def dfs_edge(self, results, node, visited_contigs, path, score, visited_nodes, depth=0):
        path.append(node)
        contig = node[:-2]
        visited_contigs.add(contig)
        if node.endswith("_b"):
            opposite = f"{contig}_e"
        else:
            opposite = f"{contig}_b"
        visited_nodes.add(node)
        if opposite not in visited_nodes:
            self.dfs_node(results, opposite, visited_contigs, path, score, visited_nodes, depth + 1)
        else:
            results.append((score, list(path)))
        visited_nodes.remove(node)
        path.pop()

    def dfs_node(self, results, node, visited_contigs, path, score, visited_nodes, depth=0):
        path.append(node)
        contig = node[:-2]
        visited_contigs.add(contig)
        visited_nodes.add(node)
        results.append((score, list(path)))
        if self.graph.get(node):
            for neighbor in self.graph[node]:
                if neighbor not in visited_nodes:
                    edge_score = self.weights[(node, neighbor)]
                    self.dfs_edge(results, neighbor, visited_contigs, path, score + edge_score, visited_nodes, depth + 1)
        else:
            print(f"{'  ' * depth}dfs_node: node {node} not found in graph, it's an end node")
        visited_nodes.remove(node)
        path.pop()

    def find_paths(self, edges_sorted):
        final_results = []
        nodes = self.nodes
        graph_to_process = self.graph.copy()
        weights_to_process = self.weights.copy()
        edges_to_process = edges_sorted.copy()
        visited_nodes = set()
        visited_contigs = set()
        while len(nodes - visited_nodes) > 0:
            print("Remaining contigs number: ", len(nodes - visited_nodes))
            start_node = '_'.join(edges_to_process[0][0][0].split("_")[:-1]) + '_b'
            print("================================================"*2)
            print(f"Starting DFS from {start_node}")
            results = []
            self.dfs_edge(results, start_node, visited_contigs, [], 0, visited_nodes)
            final_results.append(results)
            for contig in visited_contigs:
                visited_nodes.add(contig+'_b')
                visited_nodes.add(contig+'_e')
            edges_to_process = [edge for edge in edges_to_process if edge[0][0] not in visited_nodes]
            graph_to_process = {node: [neighbor for neighbor in graph_to_process[node] if neighbor not in visited_nodes] for node in graph_to_process}
            weights_to_process = {pair: weights_to_process[pair] for pair in weights_to_process if pair[0] not in visited_nodes and pair[1] not in visited_nodes}
        return final_results

class AGPAssembler:
    def __init__(self, agp_dir):
        self.agp_dir = agp_dir

    def find_connection_segments_in_agp(self, agp_path, start_node, end_node) -> tuple[list[dict], bool]:
        def parse_node(node):
            if node.endswith('_b'):
                return node[:-2], 'b'
            elif node.endswith('_e'):
                return node[:-2], 'e'
            else:
                return node, None
            
        def check_orientation(strand_s, strand_e, start_end, end_end, flag_forward):
            start_match = end_match = False
            if flag_forward == False:
                strand_s *= -1
                strand_e *= -1
            if (start_end == 'e' and strand_s == 1) or (start_end == 'b' and strand_s == -1):
                start_match = True
            if (end_end == 'b' and strand_e == 1) or (end_end == 'e' and strand_e == -1):
                end_match = True
            return start_match and end_match
        
        start_contig, start_end = parse_node(start_node)
        end_contig, end_end = parse_node(end_node)
        agp_entries = read_agp_file(agp_path)
        start_indices = [i for i, l in enumerate(agp_entries) if l['component_id'] == start_contig]
        end_indices = [i for i, l in enumerate(agp_entries) if l['component_id'] == end_contig]
        if not start_indices or not end_indices:
            print("No valid node found!!")
            return [], False
        candidate_pairs = []
        for s in start_indices:
            for e in end_indices:
                if abs(e - s) <= 2:
                    candidate_pairs.append((s, e))
        if len(candidate_pairs) == 0:
            print("No valid candidate pairs found!!")
            return [], False
        for s, e in candidate_pairs:
            strand_s = 1 if agp_entries[s]['orientation'] == '+' else -1
            strand_e = 1 if agp_entries[e]['orientation'] == '+' else -1
            if s <= e:
                if check_orientation(strand_s, strand_e, start_end, end_end, True):
                    return agp_entries[s:e+1], True
            else:
                if check_orientation(strand_s, strand_e, start_end, end_end, False):
                    return agp_entries[s:e+1], False
        print("No valid connection found!!")
        return [], False

    def extract_and_merge(self, path, dict_edge_info, agp_path):
        all_segments = []
        for idx in range(1, len(path)-1, 2):
            u, v = path[idx], path[idx+1]
            edge_info = dict_edge_info[(u, v)]
            sample = edge_info[2]
            agp_file = agp_path + f"{sample}.agp"
            print(f"AGP connection for {u} -> {v} in sample {sample}:")
            segments, flag_forward = self.find_connection_segments_in_agp(agp_file, u, v)
            if segments:
                updated_segments = merge_segments(all_segments, segments)
                if updated_segments:
                    all_segments = updated_segments
                else:
                    print("Fail to merge segments")
                    continue
        return all_segments
    
class SequenceCollector:
    def __init__(self, ref_path):
        self.ref_path = ref_path
    
    def collect(self, relevant_ref_ids, contig_fasta_path, output_path):
        list_files = [self.ref_path + id + ".fa" for id in relevant_ref_ids]
        with open(output_path + ".relevant_seq.fasta", "w") as out_f:
            for fasta in [contig_fasta_path] + list_files:
                with open(fasta, "r") as in_f:
                    content = in_f.read()
                    out_f.write(content)
                    if not content.endswith('\n'):
                        out_f.write('\n')
        return output_path + ".relevant_seq.fasta"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", required=True, help="Input path")
    parser.add_argument("-l", "--list_sample", required=True, help="file containing list of samples")
    parser.add_argument("-w", "--weight", help="weight type: allele_freq or max_alignment_score [AF|AS]", default="AS")
    parser.add_argument("-wr", "--weight_ratio", help="weight ratio, discard the edge less than ratio of the max weight the same node [0.7]", type=float, default=0.7)
    parser.add_argument("--debug", help="debug mode", action="store_true")
    parser.add_argument("--contig_fasta", help="contig fasta file, only specify if you want the patch fasta output")
    parser.add_argument("--ref_path", help="reference path, only specify if you want the patch fasta output")
    parser.add_argument("-o", "--output", help="output path", default="patch.result")
    args = parser.parse_args()

    # Collect edges from all samples
    if args.weight == "AF":
        weight_type = 0
    elif args.weight == "AS":
        weight_type = 1
    else:
        raise ValueError(f"Invalid weight type: {args.weight}")

    samples = parse_sample_list(args.list_sample)

    edge_aggregator = EdgeAggregator(args.path, samples)
    edges = edge_aggregator.aggregate()
    if args.debug:
        for pair, values in sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2])):
            print(pair, values)

    graph_builder = GraphBuilder(edges, weight_type, args.weight_ratio)
    graph, weights, nodes, dict_edge_info = graph_builder.build()
    edges_sorted = sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))

    path_finder = PathFinder(graph, weights, nodes)
    final_results = path_finder.find_paths(edges_sorted)
    agp_assembler = AGPAssembler(args.path)
    components = []

    for results in final_results:
        top_paths = sorted(results, key=lambda x: -x[0])
        len_threshold = len(top_paths[0][1])
        for score, path in top_paths:
            if score < top_paths[0][0] * 1:
                break
            print(f"Score: {score}, Path: {path}")
            print(f"Path length: {len(path)}")
            all_segments = agp_assembler.extract_and_merge(path, dict_edge_info, args.path)
            print("***********", len(all_segments), len(path), path)
            if args.debug:
                print('--------------------------------'*3)
                for seg in all_segments:
                    print('\t'.join([str(v) for v in seg.values()]))
            components.append(all_segments)
            print("================================"*3)

    # Assemble components into a single AGP file
    with open(args.output + '.agp', 'w') as f:
        for idx, component in enumerate(components):
            object = 'patch' + "{0:08}".format(idx)
            print(f"idx: {idx}, Object: {object}, component length: {len(component)}")
            for segment in component:
                f.write(f"{object}\t{segment['object_beg']}\t{segment['object_end']}\t"
                   f"{segment['part_num']}\t{segment['component_type']}\t{segment['component_id']}\t"
                   f"{segment['component_beg']}\t{segment['component_end']}\t{segment['orientation']}\n")

    if args.contig_fasta and args.ref_path:
        relevant_samples = set()
        for component in components:
            for segment in component:
                if segment['component_id'].split("_")[0] in samples:
                    relevant_samples.add(segment['component_id'].split("_")[0])
        
        print(f"Relevant samples: {relevant_samples}")
        print(f"Relevant sequences saved to {args.output}.relevant_seq.fasta...")
        sequence_collector = SequenceCollector(args.ref_path)
        sequence_collector.collect(relevant_samples, args.contig_fasta, args.output)
        cmd = [
            "./ragtag_agp2fa.py",
            args.output + ".agp",
            args.output + ".relevant_seq.fasta"
        ]
        run_oae(cmd, args.output + ".patch.fasta", "test.log")


if __name__ == "__main__":
    main()
