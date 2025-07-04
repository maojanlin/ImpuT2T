import argparse
import os
from collections import defaultdict

def aggregate_edges(input_path, samples):
    """
    edges:
    - key: (u, v) # sorted
    - value: (weight, max_score, max_sample_id)
    """
    edges = {}
    for sample in samples:
        if os.path.exists(input_path + sample + ".info"):
            with open(input_path + sample + ".info", "r") as f:
                for line in f.readlines():
                    u, v, score = line.strip().split("\t")
                    pair = tuple(sorted([u, v]))
                    if pair not in edges:
                        edges[pair] = [1, int(score), sample]
                    else:
                        if int(score) > edges[pair][1]:
                            edges[pair] = [edges[pair][0] + 1, int(score), sample]
                        else:
                            edges[pair][0] += 1
        else:
            print(f"Warning: {sample} not found in {input_path}")
    return edges

def parse_sample_list(list_sample: str) -> list[str]:
    with open(list_sample) as f:
        return [line.strip() for line in f if line.strip()]

def dfs_edge(results, graph, weights, node, visited_contigs, path, score, visited_nodes, depth=0):
    """
    DFS function to explore edges between contig ends (beginning/end).
    """
    indent = "  " * depth
    #print(f"{indent}dfs_edge: entering {node}, path: {path}")
    
    # Add current node to the path
    path.append(node)
    contig = node[:-2]
    visited_contigs.add(contig)

    # Determine the opposite end of the current contig
    if node.endswith("_b"):
        opposite = f"{contig}_e"
    else:
        opposite = f"{contig}_b"

    #print(f"{indent}dfs_edge: opposite of {node} is {opposite}")
    #print(f"{indent}dfs_edge: visited_nodes = {visited_nodes}")

    # Mark current node as visited
    visited_nodes.add(node)
    
    # ALWAYS explore the opposite end first
    if opposite not in visited_nodes:
        #print(f"{indent}dfs_edge: calling dfs_node({opposite})")
        dfs_node(results, graph, weights, opposite, visited_contigs, path, score, visited_nodes, depth + 1)
    else:
        #print(f"{indent}dfs_edge: opposite {opposite} already visited, saving path")
        results.append((score, list(path)))
    
    # Backtrack
    visited_nodes.remove(node)
    path.pop()
    #print(f"{indent}dfs_edge: backtracking from {node}")

def dfs_node(results, graph, weights, node, visited_contigs, path, score, visited_nodes, depth=0):
    """
    DFS function to explore nodes and their neighbors in the graph.
    """
    indent = "  " * depth
    #print(f"{indent}dfs_node: entering {node}, path: {path}")
    
    # Add current node to the path
    path.append(node)
    contig = node[:-2]
    visited_contigs.add(contig)
    visited_nodes.add(node)
    
    #print(f"{indent}dfs_node: saving path: {path}")
    results.append((score, list(path)))
    
    # Explore all neighbors of the current node
    #print(f"{indent}dfs_node: neighbors of {node}: {graph[node]}")
    if graph.get(node):
        for neighbor in graph[node]:
            #print(f"{indent}dfs_node: checking neighbor {neighbor}, visited: {neighbor in visited_nodes}")
            if neighbor not in visited_nodes:
                edge_score = weights[(node, neighbor)]
                #print(f"{indent}dfs_node: calling dfs_edge({neighbor}) with score {edge_score}")
                dfs_edge(results, graph, weights, neighbor, visited_contigs, path, score + edge_score, visited_nodes, depth + 1)
            else:
                pass
                #print(f"{indent}dfs_node: neighbor {neighbor} already visited, skipping")
    else:
        print(f"{indent}dfs_node: node {node} not found in graph, it's an end node")
    
    # Backtrack
    visited_nodes.remove(node)
    path.pop()

def find_connection_segments_in_agp(agp_path, start_node, end_node) -> tuple[list[str], bool]:
    """
    Extracts AGP lines between start_node and end_node (inclusive) from the AGP file,
    considering the _b/_e suffix for directionality and both forward and reverse order.
    """
    def parse_node(node):
        if node.endswith('_b'):
            return node[:-2], 'b'
        elif node.endswith('_e'):
            return node[:-2], 'e'
        else:
            return node, None

    def check_orientation(strand_s, strand_e, start_end, end_end, flag_forward):
        print(strand_s, strand_e, start_end, end_end, flag_forward)
        if flag_forward == False:
            strand_s *= -1
            strand_e *= -1
        
        if start_end == 'e' and strand_s == 1 or start_end == 'b' and strand_s == -1:
            start_match = True
        if end_end == 'b' and strand_e == 1 or end_end == 'e' and strand_e == -1:
            end_match = True
        return start_match and end_match

    start_contig, start_end = parse_node(start_node)
    end_contig, end_end = parse_node(end_node)

    lines = []
    with open(agp_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue
            lines.append(line.strip())

    # Find all indices for start and end contigs
    start_indices = [i for i, l in enumerate(lines) if l.split("\t")[5] == start_contig]
    end_indices = [i for i, l in enumerate(lines) if l.split("\t")[5] == end_contig]

    if not start_indices or not end_indices:
        print("No valid node found!!")
        return [], False

    # Find the closest pair (forward or reverse)
    candidate_pairs = []
    for s in start_indices:
        for e in end_indices:
            if abs(e - s) <= 2:
                candidate_pairs.append((s, e))

    if len(candidate_pairs) == 0:
        print("No valid candidate pairs found!!")
        return [], False

    for s, e in candidate_pairs:
        # check if the contigs agree with start_end and end_end
        strand_s = 1 if lines[s].split("\t")[8] == '+' else -1
        strand_e = 1 if lines[e].split("\t")[8] == '+' else -1
        if s <= e: # Forward order
            if check_orientation(strand_s, strand_e, start_end, end_end, True):
                return lines[s:e+1], True
        else: # Reverse order
            if check_orientation(strand_s, strand_e, start_end, end_end, False):
                return lines[s:e+1], False
    print("No valid connection found!!")
    return [], False

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", required=True, help="Input path")
    parser.add_argument("-l", "--list_sample", required=True, help="file containing list of samples")
    parser.add_argument("-w", "--weight", help="weight type: allele_freq or max_alignment_score [AF|AS]", default="AS")
    #parser.add_argument("-wt", "--weight_threshold", help="weight threshold [AF:5|AS:6000]", type=int)
    parser.add_argument("-wr", "--weight_ratio", help="weight ratio, discard the edge less than ratio of the max weight the same node [0.7]", type=float, default=0.7)
    args = parser.parse_args()

    # Collect edges from all samples
    if args.weight == "AF":
        weight_type = 0
    elif args.weight == "AS":
        weight_type = 1
    else:
        raise ValueError(f"Invalid weight type: {args.weight}")

    samples = parse_sample_list(args.list_sample)
    edges = aggregate_edges(args.path, samples)
    for pair, values in sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2])):
        print(pair, values)

    # Build the graph
    edges_sorted = sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
    graph = defaultdict(list)
    weights = {}
    nodes = set()
    dict_node_score = defaultdict(list)
    dict_edge_info = {}

    for pair, values in edges_sorted:
        u, v = pair
        w = values[weight_type] # 0 for allele frequency, 1 for max alignment score
        dict_node_score[u].append((w, v))
        dict_node_score[v].append((w, u))
        dict_edge_info[(u, v)] = values
        dict_edge_info[(v, u)] = values

    filtered_edges = set()
    for node, scores in dict_node_score.items():
        max_score = max(scores, key=lambda x: x[0])[0]
        for score, neighbor in scores:
            if score < max_score * args.weight_ratio:
                filtered_edges.add((node, neighbor))
                filtered_edges.add((neighbor, node))

    for pair, values in edges_sorted:
        u, v = pair
        w = values[weight_type] # 0 for allele frequency, 1 for max alignment score
        if (u, v) in filtered_edges or (v, u) in filtered_edges:
            continue

        # bidirectional edges
        graph[u].append(v)
        graph[v].append(u)
        weights[(u, v)] = w
        weights[(v, u)] = w
        nodes.update([u, v])
        nodes.update([v, u])

    # Add implicit internal b -> e links with weight = 0
    for node in nodes:
        if node.endswith("_b"):
            contig = node[:-2]
            e_node = f"{contig}_e"
            graph[node].append(e_node)
            weights[(node, e_node)] = 0

    
    # Find the best path and iterate until all contigs are visited
    
    final_results = []
    contigs = set(contig[:-2] for contig in nodes)
    graph_to_process = graph
    weights_to_process = weights
    edges_to_process = edges_sorted
    visited_nodes = set()
    visited_contigs = set()
    while len(nodes - visited_nodes) > 0:
        print("Remaining contigs number: ", len(nodes - visited_nodes))
        start_node = '_'.join(edges_to_process[0][0][0].split("_")[:-1]) + '_b'
        print("================================================"*2)
        print(f"Starting DFS from {start_node}")
        results = []
        dfs_edge(results, graph_to_process, weights_to_process, start_node, visited_contigs, [], 0, visited_nodes)
        final_results.append(results)

        for contig in visited_contigs:
            visited_nodes.add(contig+'_b')
            visited_nodes.add(contig+'_e')

        edges_to_process = [edge for edge in edges_to_process if edge[0][0] not in visited_nodes]
        graph_to_process = {node: [neighbor for neighbor in graph_to_process[node] if neighbor not in visited_nodes] for node in graph_to_process}
        weights_to_process = {pair: weights_to_process[pair] for pair in weights_to_process if pair[0] not in visited_nodes and pair[1] not in visited_nodes}
        """
        print("--------------------------------"*3)
        print(len(visited_contigs))
        print(len(edges_to_process))
        print(len(graph_to_process))
        print(len(weights_to_process))
        """
    # Show top 3 paths
    for results in final_results:
        top_paths = sorted(results, key=lambda x: -x[0])#[:3]
        len_threshold = len(top_paths[0][1])
        for score, path in top_paths:
            #if len(path) < len_threshold:
            #    break
            if score < top_paths[0][0] * 1:
                break
            print(f"Score: {score}, Path: {path}")
            print(f"Path length: {len(path)}")
            for idx in range(1, len(path)-1, 2):
                u, v = path[idx], path[idx+1]
                edge_info = dict_edge_info[(u, v)]
                sample = edge_info[2]  # Assuming sample is stored as the third element
                agp_file = args.path + f"{sample}.agp"
                print(f"AGP connection for {u} -> {v} in sample {sample}:")
                segments, flag_forward = find_connection_segments_in_agp(agp_file, u, v)  # Remove _b/_e if needed
                if flag_forward == False:
                    segments = list(reversed(segments))
                for seg in segments:
                    print(seg)
                print("---" * 10)
            print("---"*20)
