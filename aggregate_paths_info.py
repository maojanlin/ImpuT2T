import argparse
import os
from collections import defaultdict
from typing import List, Dict
from combine_connections import read_agp_file, merge_segments, reverse_segments
import subprocess
from ragtag_utilities.utilities import run_oae
import pysam

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
    def __init__(self, edges, weight_type, weight_ratio, AS_threshold, AF_threshold, sample_size):
        self.edges = edges
        self.weight_type = weight_type
        self.weight_ratio = weight_ratio
        self.AS_threshold = AS_threshold
        self.AF_threshold = AF_threshold
        self.sample_size = sample_size
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
            AF, AS, _ = values
            w = values[self.weight_type]
            if AS < self.AS_threshold: # filtering through AS threshold
                continue
            if AF <= self.AF_threshold * self.sample_size * 2: # filtering through AF threshold
                continue
            if (u, v) in filtered_edges or (v, u) in filtered_edges:
                continue
            print("\t", (u, v), values)
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
            elif node.endswith("_e"):
                contig = node[:-2]
                b_node = f"{contig}_b"
                self.graph[node].append(b_node)
                self.weights[(node, b_node)] = 0
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

    def find_end_nodes(self):
        set_contigs = set([node[:-2] for node in self.nodes])
        end_nodes = set()
        for contig in set_contigs:
            s_node = contig + '_b'
            e_node = contig + '_e'
            if s_node not in self.graph:
                end_nodes.add(s_node)
            if e_node not in self.graph:
                end_nodes.add(e_node)
        return end_nodes
    
    # ========================= sub-functions for find_paths =========================
    def give_start_nodes(self, edges_to_process, end_nodes, flag_exhaustive):
        if len(end_nodes) == 0:
            start_node = '_'.join(edges_to_process[0][0][0].split("_")[:-1]) + '_b'
            start_e_node = '_'.join(edges_to_process[0][0][0].split("_")[:-1]) + '_e'
            start_nodes = [start_node, start_e_node]
        else:
            if flag_exhaustive == False:
                start_nodes = sorted(end_nodes, key=lambda x: int(x.split("_")[-2]))[:1]
            else:
                start_nodes = sorted(end_nodes, key=lambda x: int(x.split("_")[-2]))
        return start_nodes
    
    def find_candidate_paths(self, results, end_nodes):
        dict_candidate_paths = defaultdict(list)
        top_results = sorted(results, key=lambda x: -x[0])
        max_score = top_results[0][0]
        for score, path in top_results:
            terminal_node = path[-1]
            if terminal_node in end_nodes:
                dict_candidate_paths[terminal_node].append((score, path))
            elif score >= max_score:
                dict_candidate_paths["longest"].append((score, path))
        dict_trimmed_paths = defaultdict(list)
        for terminal_node, paths in dict_candidate_paths.items():
            if terminal_node == "longest":
                dict_trimmed_paths[terminal_node] = paths
            else:
                sorted_paths = sorted(paths, key=lambda x: -x[0])
                for path in sorted_paths:
                    if path[0] >= sorted_paths[0][0]:
                        dict_trimmed_paths[terminal_node].append(path)
                    else:
                        break
        return dict_trimmed_paths
    # =================================================================================

    def search_paths(self, final_results, current_path, visited_contigs, visited_nodes, edges_sorted, end_nodes, flag_exhaustive):
        new_edges_sorted = [edge for edge in edges_sorted if edge[0][0] not in visited_nodes and edge[0][1] not in visited_nodes]
        #if len(new_edges_sorted) == 0:
        #    final_results.append(current_path)
        #    print("......"*20, "Empty edges found!")
        #    return final_results
        start_nodes = self.give_start_nodes(new_edges_sorted, end_nodes, flag_exhaustive)
        #print(f"visited_nodes: {visited_nodes}")
        #print(f"[0]: {start_nodes[0]}")
        if start_nodes[0] in visited_nodes:
            final_results.append(current_path)
            print("......"*20, "Duplicate start node found!")
            return final_results
        
        for idx, start_node in enumerate(start_nodes):
            results = []
            nodes = self.nodes
            #print("Remaining contigs number: ", len(nodes - visited_nodes))
            print(f"Starting DFS from {start_node}")
            self.dfs_edge(results, start_node, visited_contigs.copy(), [], 0, visited_nodes)
            # choose the longest path with the highest score or the path stop with end-nodes
            dict_candidate_paths = self.find_candidate_paths(results, end_nodes)
            for terminal_node, paths in dict_candidate_paths.items():
                print(f"terminal_node: {terminal_node}")
                for path in paths:
                    print(f"path: {path}")
            print("--------------------------------------------------"*2)
            for idy, (terminal_node, paths) in enumerate(dict_candidate_paths.items()):
                for idz, path in enumerate(sorted(paths, key=lambda x: -x[0])):
                    print(idx, idy, idz, terminal_node, path)
                    contig_in_path = set([node[:-2] for node in path[1]])
                    set_visited_contigs = visited_contigs.union(contig_in_path)
                    set_visited_nodes = set()
                    for visited_contig in set_visited_contigs:
                        set_visited_nodes.update([visited_contig + '_b', visited_contig + '_e'])

                    if len(nodes - set_visited_nodes) == 0:
                        final_results.append(current_path + [path])
                        print("......"*20, "Final result found!")
                        #return final_results
                    else:
                        new_end_nodes = end_nodes - set_visited_nodes
                        #print(f"new_end_nodes: {new_end_nodes}")
                        #print(f"path: {path}")
                        self.search_paths(final_results, current_path + [path], set_visited_contigs, set_visited_nodes, edges_sorted, new_end_nodes, flag_exhaustive)

        return final_results
    
    def standardize_traversals(self, paths, end_nodes):
        new_paths = []
        total_score = 0
        non_end_node_count = 0
        dangling_contig_count = 0
        for score, path in sorted(paths, key=lambda x: (-x[0], x[1])):
            total_score += score
            if len(path) <= 2:
                dangling_contig_count += 1
            else:
                if path[-1] not in end_nodes:
                    non_end_node_count += 1
                if path[0] not in end_nodes:
                    non_end_node_count += 1
            # check the orientation of the path
            if int(path[-1].split("_")[-2]) > int(path[0].split("_")[-2]):
                new_paths.append((score, tuple(path)))
            else:
                new_paths.append((score, tuple(path[::-1])))
        new_paths = sorted(new_paths, key=lambda x: int(x[1][0].split("_")[-2]))
        return (total_score, non_end_node_count, dangling_contig_count, tuple(new_paths))                
                    

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
            start_e_node = '_'.join(edges_to_process[0][0][0].split("_")[:-1]) + '_e'
            if start_e_node not in graph_to_process:
                print("Reverse strand for start node", start_node, "->", start_e_node)
                start_node = start_e_node
            #print(start_node, graph_to_process[start_node])
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
                    return agp_entries[e:s+1], False
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
                if not flag_forward:
                    segments = reverse_segments(segments)
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
    
def output_agp(components, contig_fasta_path, output_path):
    used_contigs = set()
    contig_lengths = {}
    contig_order = []
    fai = pysam.FastaFile(contig_fasta_path)
    contig_order = list(fai.references)
    contig_lengths = {name: fai.get_reference_length(name) for name in contig_order}
    fai.close()
    
    # Create a mapping from first component to component index
    # Create a list of all items to output, sorted by their first contig's position in FASTA
    output_items = []
    component_first_contig = {}
    component_idx = 0
    for component in components:
        if component:
            first_contig = component[0]['component_id']
            output_items.append(('patch', component_idx, component, first_contig))
            component_first_contig[first_contig] = component_idx
            for segment in component:
                used_contigs.add(segment['component_id'])
            component_idx += 1
    
    # Add unplaced contigs
    for contig in contig_order:
        if contig not in used_contigs:
            output_items.append(('unplaced', contig, None, contig))
    
    # Sort by position in FASTA order
    contig_to_position = {contig: idx for idx, contig in enumerate(contig_order)}
    output_items.sort(key=lambda x: contig_to_position.get(x[3], float('inf')))
    
    with open(output_path + '.agp', 'w') as f:
        for item_type, item_id, component, first_contig in output_items:
            if item_type == 'patch':
                object = 'patch' + "{0:08}".format(item_id)
                print(f"idx: {item_id}, Object: {object}, component length: {len(component)}")
                for segment in component:
                    f.write(f"{object}\t{segment['object_beg']}\t{segment['object_end']}\t"
                       f"{segment['part_num']}\t{segment['component_type']}\t{segment['component_id']}\t"
                       f"{segment['component_beg']}\t{segment['component_end']}\t{segment['orientation']}\n")
            else:  # unplaced
                f.write(f"{item_id}\t1\t{contig_lengths[item_id]}\t1\tW\t{item_id}\t1\t{contig_lengths[item_id]}\t+\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", required=True, help="Input path")
    parser.add_argument("-l", "--list_sample", required=True, help="file containing list of samples")
    parser.add_argument("-w", "--weight", help="weight type: allele_freq or max_alignment_score [AF|AS]", default="AS")
    parser.add_argument("-wr", "--weight_ratio", help="weight ratio, discard the edge less than ratio of the max weight the same node [0.7]", type=float, default=0.7)
    parser.add_argument("-ast", "--AS_threshold", help="AS threshold, discard the edge less than threshold [6000]", type=int, default=6000)
    parser.add_argument("-aft", "--AF_threshold", help="AF threshold, discard the edge less than threshold [0.005]", type=float, default=0.005)
    parser.add_argument("--debug", help="debug mode", action="store_true")
    parser.add_argument("--contig_fasta", required=True, help="contig fasta file, needed for unplaced contigs")
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

    graph_builder = GraphBuilder(edges, weight_type, args.weight_ratio, args.AS_threshold, args.AF_threshold, len(samples))
    graph, weights, nodes, dict_edge_info = graph_builder.build()
    edges_sorted = sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
    print("graph", "--------------------------------"*2)
    for ele in graph:
        print(ele, graph[ele])
    print("--------------------------------"*2)

    path_finder = PathFinder(graph, weights, nodes)
    end_nodes = path_finder.find_end_nodes()
    print("================================================"*2)
    print("End nodes: ", end_nodes)
    print("================================================"*2)
    test_results = path_finder.search_paths([], [], set(), set(), edges_sorted, end_nodes, True)
    print("TEST_RESULTS", "########################"*2)
    print("TEST_RESULTS", "########################"*2)
    print("TEST_RESULTS", "########################"*2)
    print("TEST_RESULTS", "########################"*2)
    standardized_result = set()
    for ele in test_results:
        standardized_result.add(path_finder.standardize_traversals(ele, end_nodes))
    sorted_standardized_result = sorted(list(standardized_result), key=lambda x: -x[0])
    for ele in sorted_standardized_result:
        print(ele)
    
    
    components = []
    agp_assembler = AGPAssembler(args.path)
    total_score_1st, non_end_node_count_1st, dangling_contig_count_1st, paths = sorted_standardized_result[0]
    for score, path in paths:
        print(f"Total Score: {total_score_1st}, Score: {score}, Path: {path}")
        print(f"Path length: {len(path)}")
        all_segments = agp_assembler.extract_and_merge(path, dict_edge_info, args.path)
        print("***********", len(all_segments), len(path), path)
        if args.debug:
            print('--------------------------------'*3)
            for seg in all_segments:
                print('\t'.join([str(v) for v in seg.values()]))
        components.append(all_segments)
    # get the alternative components if other components has higher score, lower number of non-end nodes or dangling contigs
    alt_components = []
    for total_score, non_end_node_count, dangling_contig_count, paths in sorted_standardized_result[1:]:
        alt_component = []
        if total_score >= total_score_1st or non_end_node_count < non_end_node_count_1st or dangling_contig_count < dangling_contig_count_1st:
            print("Alternative component found!")
            for score, path in paths:
                print(f"Total Score: {total_score}, Score: {score}, Path: {path}")
                print(f"Path length: {len(path)}")
                all_segments = agp_assembler.extract_and_merge(path, dict_edge_info, args.path)
                print("***********", len(all_segments), len(path), path)
                alt_component.append(all_segments)
        if alt_component:
            alt_components.append(alt_component)

    # Assemble components into a single AGP file
    output_agp(components, args.contig_fasta, args.output)
    for idx, alt_component in enumerate(alt_components):
        output_agp(alt_component, args.contig_fasta, args.output + f".alt_{idx}")

    if args.ref_path:
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
