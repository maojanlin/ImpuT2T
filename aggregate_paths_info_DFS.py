import argparse
import os
from collections import defaultdict
from combine_connections import read_agp_file, merge_segments, reverse_segments
from ragtag_utilities.utilities import run_oae, get_ragtag_version, get_impuT2T_version

os.environ["MKL_THREADING_LAYER"] = "GNU"
import pysam
import numpy as np
from sklearn.mixture import GaussianMixture


def parse_sample_list(list_sample: str) -> list[str]:
    with open(list_sample) as f:
        return [line.strip() for line in f if line.strip()]

class EdgeAggregator:
    def __init__(self, input_path: str, samples: list[str], gmm: bool):
        self.input_path = input_path
        self.samples = samples
        self.gmm = gmm
        self.raw_edges = {}
        self.edges = {}

    def fit_gmm(self, x, max_components=10, cutoff=0.1):
        bics = []
        models = []
        x = x.reshape(-1, 1)
        for k in range(1, max_components + 1):
            gmm = GaussianMixture(n_components=k, covariance_type='full').fit(x)
            weights = gmm.weights_
            significant_components = weights > cutoff
            if significant_components.sum() < len(significant_components):
                break
            bics.append(gmm.bic(x))
            models.append(gmm)
        best_k = np.argmin(bics)
        best_model = models[best_k]
        return best_model, best_k+1
    
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
                        u, v, score, dist = line.strip().split("\t")
                        pair = tuple(sorted([u, v]))
                        if pair not in self.raw_edges:
                            self.raw_edges[pair] = [[sample], [(int(score), int(dist))]]
                        else:
                            self.raw_edges[pair][0].append(sample)
                            self.raw_edges[pair][1].append((int(score), int(dist)))
            else:
                print(f"Warning: {sample} not found in {self.input_path}")
        if self.gmm:
            for pair, info in self.raw_edges.items():
                scores = [s for s, d in info[1]]
                dists = [d for s, d in info[1]]
                relevant_samples = info[0]
                if len(relevant_samples) > 10:
                    gmm, best_k = self.fit_gmm(np.array(dists))
                    means = gmm.means_.flatten()
                    stds = np.sqrt(gmm.covariances_.flatten())  # for 1D case
                    print(f"====================================== Best k: {best_k} ======================================")
                    print(f"======== pair: {pair} ========")
                    print(f"GMM means: {means}, stds: {stds}")
                    gmm_scores = []
                    gmm_samples = []
                    for mean, std in zip(means, stds):
                        candidate_dists = []
                        candidate_scores = []
                        candidate_samples = []
                        for idx in range(len(dists)):
                            upper_bound = mean + 1.5*std
                            lower_bound = mean - 1.5*std if mean - 1.5*std > 0 and mean > 0 else mean - std
                            if dists[idx] <= upper_bound and dists[idx] >= lower_bound:
                                candidate_dists.append(dists[idx])
                                candidate_scores.append(scores[idx])
                                candidate_samples.append(relevant_samples[idx])
                        if len(candidate_dists) > 0:
                            best_close_idx = np.argmax(candidate_scores)
                            gmm_scores.append(candidate_scores[best_close_idx])
                            gmm_samples.append(candidate_samples[best_close_idx])
                    final_selected_idx = np.argmax(gmm_scores)
                    self.edges[pair] = [len(relevant_samples), gmm_scores[final_selected_idx], gmm_samples[final_selected_idx]]
                else:
                    idx = np.argmax(scores)
                    self.edges[pair] = [len(relevant_samples), scores[idx], relevant_samples[idx]]
        else:
            for pair, info in self.raw_edges.items():
                idx = np.argmax(info[1], key=lambda x: x[0])
                self.edges[pair] = [len(info[0]), info[1][idx][0], info[0][idx]]

        return self.edges

class GraphBuilder:
    def __init__(self, edges, weight_type, weight_ratio, AS_threshold, AF_threshold, contig_fasta, sample_size):
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
        self.filtered_edges = set()
        fai = pysam.FastaFile(contig_fasta)
        self.contig_lengths = {contig: fai.get_reference_length(contig) for contig in fai.references}

    def build(self):
        edges_sorted = sorted(self.edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
        for pair, values in edges_sorted:
            u, v = pair
            w = values[self.weight_type]
            self.dict_node_score[u].append((w, v))
            self.dict_node_score[v].append((w, u))
            self.dict_edge_info[(u, v)] = values
            self.dict_edge_info[(v, u)] = values

        for node, scores in self.dict_node_score.items():
            max_score = max(scores, key=lambda x: x[0])[0]
            for score, neighbor in scores:
                if score < max_score * self.weight_ratio:
                    self.filtered_edges.add((node, neighbor))
                    self.filtered_edges.add((neighbor, node))

        dict_node_AF_count = defaultdict(int)
        for pair, values in edges_sorted:
            AF, AS, _ = values
            dict_node_AF_count[pair[0]] += AF
            dict_node_AF_count[pair[1]] += AF

        for pair, values in edges_sorted:
            u, v = pair
            AF, AS, _ = values
            #w = values[self.weight_type]
            w = AS*AF*2/(dict_node_AF_count[u]+dict_node_AF_count[v])
            if AS < self.AS_threshold: # filtering through AS threshold
                self.filtered_edges.add((u, v))
                self.filtered_edges.add((v, u))
                continue
            if AF <= self.AF_threshold * self.sample_size * 2: # filtering through AF threshold
                self.filtered_edges.add((u, v))
                self.filtered_edges.add((v, u))
                continue
            if (u, v) in self.filtered_edges or (v, u) in self.filtered_edges:
                self.filtered_edges.add((u, v))
                self.filtered_edges.add((v, u))
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
            elif node.endswith("_e"):
                contig = node[:-2]
                b_node = f"{contig}_b"
                self.graph[node].append(b_node)
        
        for contig in self.contig_lengths:
            b_node = f"{contig}_b"
            e_node = f"{contig}_e"
            self.weights[(b_node, e_node)] = self.contig_lengths[contig]/10
            self.weights[(e_node, b_node)] = self.contig_lengths[contig]/10
        return self.graph, self.weights, self.nodes, self.dict_edge_info

class GraphBuilder_2nd_pass:
    def __init__(self, edges, weight_type, weight_ratio, AS_threshold, AF_threshold, contig_fasta, sample_size):
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
        self.filtered_edges = set()
        fai = pysam.FastaFile(contig_fasta)
        self.contig_lengths = {contig: fai.get_reference_length(contig) for contig in fai.references}
        self.dict_compressed_2_component = {}

    def convert_first_stage_results(self, first_stage_results):
        target_path = first_stage_results[0]

        dict_node_2_compressed = {}
        self.dict_compressed_2_component = {}
        set_used_contigs = set()

        _, _, _, components = target_path
        for component in components:
            local_score, local_path = component
            average_order = sum([int(node.split("_")[-2]) for node in local_path]) / len(local_path)
            compressed_name = f"impuT2T_{average_order:.2f}"
            u = compressed_name + '_b'
            v = compressed_name + '_e'
            self.dict_compressed_2_component[compressed_name] = component
            dict_node_2_compressed[local_path[0]] = u
            dict_node_2_compressed[local_path[-1]] = v
            set_used_contigs.update([node[:-2] for node in local_path])
            #self.nodes.update([u, v])
            #self.nodes.update([v, u])
            #self.graph[u].append(v)
            #self.graph[v].append(u)
            self.weights[(u, v)] = local_score
            self.weights[(v, u)] = local_score

            print(f"u: {u}, v: {v}, component: {component}, weight: {component[0]}")
        
        edges_sorted = sorted(self.edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
        compressed_edges = []
        print(f"edges_sorted: {edges_sorted[0]}")
        for pair, values in edges_sorted:
            u, v = pair
            if u in dict_node_2_compressed and v in dict_node_2_compressed:
                new_u = dict_node_2_compressed[u]
                new_v = dict_node_2_compressed[v]
                compressed_edges.append(((new_u, new_v), values))
            elif u in dict_node_2_compressed and v[:-2] not in set_used_contigs:
                new_u = dict_node_2_compressed[u]
                compressed_edges.append(((new_u, v), values))
            elif v in dict_node_2_compressed and u[:-2] not in set_used_contigs:
                new_v = dict_node_2_compressed[v]
                compressed_edges.append(((u, new_v), values))
            else:
                if u[:-2] in set_used_contigs or v[:-2] in set_used_contigs:
                    continue
                else:
                    compressed_edges.append(((u, v), values))

        #print(f"compressed_edges: {compressed_edges[0]}")

        for pair, values in compressed_edges:
            w = values[self.weight_type]
            self.dict_node_score[u].append((w, v))
            self.dict_node_score[v].append((w, u))
            self.dict_edge_info[(u, v)] = values
            self.dict_edge_info[(v, u)] = values

        dict_node_AF_count = defaultdict(int)
        for pair, values in compressed_edges:
            AF, AS, _ = values
            dict_node_AF_count[pair[0]] += AF
            dict_node_AF_count[pair[1]] += AF

        for pair, values in compressed_edges:
            u, v = pair
            AF, AS, _ = values
            #w = values[self.weight_type]
            w = AS*AF*2/(dict_node_AF_count[u]+dict_node_AF_count[v])
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
            elif node.endswith("_e"):
                contig = node[:-2]
                b_node = f"{contig}_b"
                self.graph[node].append(b_node)
        
        for contig in self.contig_lengths:
            if contig in set_used_contigs:
                print(f"contig: {contig} is used in the first stage")
                continue
            b_node = f"{contig}_b"
            e_node = f"{contig}_e"
            self.weights[(b_node, e_node)] = self.contig_lengths[contig]/10
            self.weights[(e_node, b_node)] = self.contig_lengths[contig]/10
        return self.graph, self.weights, self.nodes, self.dict_edge_info

    def recover_unused_components(self, results):
        updated_results = [] # counterpart of traversals
        print(results)
        for traversal in results:
            path_info = traversal[3]
            set_used_contigs = set()
            updated_components = []
            for component in path_info:
                set_used_contigs.update(node[:-2] for node in component[1])
                updated_components.append(component)
            for component_contig in self.dict_compressed_2_component.keys():
                if component_contig not in set_used_contigs:
                    print(f"component_contig: {component_contig} {self.weights[(f'{component_contig}_b', f'{component_contig}_e')]}")
                    edge_info = (f'{component_contig}_b', f'{component_contig}_e')
                    updated_components.append((self.weights[edge_info], edge_info))
            updated_components = sorted(updated_components, key=lambda x: float(x[1][0].split("_")[-2]))
            updated_results.append(tuple([traversal[0], traversal[1], traversal[2], tuple(updated_components)]))
        return updated_results

    def decompress_path(self, path):
        new_path = []
        idx = 0
        print(f"self.dict_compressed_2_component: {self.dict_compressed_2_component}")
        while idx < len(path):
            node = path[idx]
            print(f"node: {node}, {node[:-2] in self.dict_compressed_2_component}")
            if node[:-2] in self.dict_compressed_2_component:
                if node.endswith("_b"):
                    new_path += list(self.dict_compressed_2_component[node[:-2]][1])
                    #print(f"find_path: {self.dict_compressed_2_component[node[:-2]]}")
                else:
                    new_path += list(self.dict_compressed_2_component[node[:-2]][1][::-1])
                    #print(f"find_path: {self.dict_compressed_2_component[node[:-2]][::-1]}")
                idx += 2
            else:
                new_path.append(node)
                idx += 1
        return tuple(new_path)


class PathFinder:
    def __init__(self, graph, weights, nodes):
        self.graph = graph
        self.weights = weights
        self.nodes = nodes

    def dfs_edge(self, results, node, visited_contigs, path, score, visited_nodes, depth=0):
        path.append(node)
        print(f"{'-|' * depth}dfs_edge: node {node}")
        contig = node[:-2]
        visited_contigs.add(contig)
        if node.endswith("_b"):
            opposite = f"{contig}_e"
        else:
            opposite = f"{contig}_b"
        node_score = self.weights[(node, opposite)]
        visited_nodes.add(node)
        if opposite not in visited_nodes:
            self.dfs_node(results, opposite, visited_contigs, path, score + node_score, visited_nodes, depth + 1)
        else:
            results.append((score + node_score, list(path)))
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
        print(f"set_contigs: {set_contigs}")
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
                start_nodes = sorted(end_nodes, key=lambda x: float(x.split("_")[-2]))[:1]
            else:
                start_nodes = sorted(end_nodes, key=lambda x: float(x.split("_")[-2]))
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
        #print(f"end_nodes: {end_nodes}")
        #print(f"new_edges_sorted: {new_edges_sorted}")
        #print(f"visited_nodes: {visited_nodes}")
        print(f"[0]: {start_nodes[0]}")
        if start_nodes[0] in visited_nodes:
            final_results.append(current_path)
            #print("......"*20, "Duplicate start node found!")
            return final_results
        
        for idx, start_node in enumerate(start_nodes):
            results = []
            nodes = self.nodes
            print("Remaining contigs number: ", len(nodes - visited_nodes))
            #print("Remaining contigs number: ", nodes - visited_nodes)
            #print(f"Starting DFS from {start_node}")
            self.dfs_edge(results, start_node, visited_contigs.copy(), [], 0, visited_nodes)
            # choose the longest path with the highest score or the path stop with end-nodes
            dict_candidate_paths = self.find_candidate_paths(results, end_nodes)
            #for terminal_node, paths in dict_candidate_paths.items():
            #    print(f"terminal_node: {terminal_node}")
            #    for path in paths:
            #        print(f"path: {path}")
            print("--------------------------------------------------"*2)
            for idy, (terminal_node, paths) in enumerate(dict_candidate_paths.items()):
                for idz, path in enumerate(sorted(paths, key=lambda x: -x[0])):
                    print(idx, idy, idz, terminal_node, path)
                    contig_in_path = set([node[:-2] for node in path[1]])
                    set_visited_contigs = visited_contigs.union(contig_in_path)
                    set_visited_nodes = set()
                    for visited_contig in set_visited_contigs:
                        set_visited_nodes.update([visited_contig + '_b', visited_contig + '_e'])

                    remaining_edges = [edge for edge in edges_sorted if edge[0][0] not in set_visited_nodes and edge[0][1] not in set_visited_nodes]
                    if len(nodes - set_visited_nodes) == 0 or len(remaining_edges) == 0:
                        final_results.append(current_path + [path])
                        #print("......"*20, "Final result found!")
                        #return final_results
                    else:
                        new_end_nodes = end_nodes - set_visited_nodes
                        #print(f"nodes - set_visited_nodes: {nodes - set_visited_nodes}")
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
                # check the orientation of the path
                if path[0].endswith("_b"):
                    new_paths.append((int(score), tuple(path)))
                else:
                    new_paths.append((int(score), tuple(path[::-1])))
            else:
                if path[-1] not in end_nodes:
                    non_end_node_count += 1
                if path[0] not in end_nodes:
                    non_end_node_count += 1
                # check the orientation of the path
                if float(path[-1].split("_")[-2]) > float(path[0].split("_")[-2]):
                    new_paths.append((int(score), tuple(path)))
                else:
                    new_paths.append((int(score), tuple(path[::-1])))
        new_paths = sorted(new_paths, key=lambda x: float(x[1][0].split("_")[-2]))
        return (int(total_score), non_end_node_count, dangling_contig_count, tuple(new_paths))

    def second_pass_merge(self, paths, possible_connections):
        # Build mapping from end nodes to component indices
        print('--------------------------------'*2, "second_pass_merge")
        end_to_idx = {}
        for idx, path in enumerate(paths):
            end_to_idx[path[1][0]] = idx  # start
            end_to_idx[path[1][-1]] = idx  # end

        # Build adjacency list for connections
        adj = defaultdict(list)
        in_degree = defaultdict(int)
        print("possible_connections: ", sorted(possible_connections, key=lambda x: (x[0], x[1])))
        for idx, idj, (from_node, to_node) in sorted(possible_connections, key=lambda x: (x[0], x[1])):
            adj[idx].append((idj, from_node, to_node))
            in_degree[idj] += 1

        print("adj: ", adj)
        print("in_degree: ", in_degree)

        # Find starting components (those not pointed to by any connection)
        start_indices = [i for i in range(len(paths)) if in_degree[i] == 0]

        merged_paths = []
        merged_indices = []
        visited = set()

        # Track available ends for each component
        available_ends = {}
        for idx, path in enumerate(paths):
            available_ends[idx] = {path[1][0], path[1][-1]}  # Both ends start available
        
        for start in start_indices:
            if start in visited:
                continue
            current = start
            merged = list(paths[current][1])
            component_indices = [current]
            visited.add(current)
            
            # Mark the start end as used
            available_ends[current].discard(paths[current][1][0])
            
            while adj[current]:
                next_idx, from_node, to_node = adj[current][0]
                
                # Validate: from_node must be available in current component
                if from_node not in available_ends[current]:
                    print(f"Warning: {from_node} not available in component {current}")
                    # Remove this connection and try the next one
                    adj[current].pop(0)
                    continue
                    
                # Validate: to_node must be available in next component
                if to_node not in available_ends[next_idx]:
                    print(f"Warning: {to_node} not available in component {next_idx}")
                    adj[current].pop(0)
                    continue
                
                # Mark both ends as used
                available_ends[current].discard(from_node)
                available_ends[next_idx].discard(to_node)
                
                # Now do the merge with proper orientation
                if to_node == paths[next_idx][1][0]:
                    merged.extend(paths[next_idx][1])
                elif to_node == paths[next_idx][1][-1]:
                    merged.extend(paths[next_idx][1][::-1])
                else:
                    print("******error: "*100, to_node, "not found in ", paths[next_idx][1])
                    exit()
                    
                component_indices.append(next_idx)
                visited.add(next_idx)
                current = next_idx
                
            merged_paths.append(tuple(merged))
            merged_indices.append(tuple(component_indices))

        # Add any unvisited components (not connected to anything)
        for i in range(len(paths)):
            if i not in visited:
                merged_paths.append(paths[i][1])
                merged_indices.append((i,))

        return merged_indices, merged_paths


    def second_pass_connection(self, sorted_standardized_result, filtered_edges):
        print("filtered_edges: ", sorted(filtered_edges))
        second_pass_result = []
        print("@@@@@@@@@@@@@@@@@@@@@@"*47)
        print(self.nodes)
        print(self.graph)
        print(self.weights)
        print("@@@@@@@@@@@@@@@@@@@@@@"*47)

        for score, non_end_node_count, dangling_contig_count, paths in sorted_standardized_result:
            end_nodes = []
            for path_info in paths:
                _, path = path_info
                end_nodes.append((path[0], path[-1]))

            possible_connections = set()
            for idx, (end_node_b, end_node_e) in enumerate(end_nodes):
                same_contig_flag_i = False
                if end_node_b[:-2] == end_node_e[:-2]:
                    same_contig_flag_i = True
                for idj in range(idx+1, len(end_nodes)):
                    same_contig_flag_j = False
                    (end_node_b_2, end_node_e_2) = end_nodes[idj]
                    if end_node_b_2[:-2] == end_node_e_2[:-2]:
                        same_contig_flag_j = True
                    if end_node_b not in (end_node_b_2, end_node_e_2) and end_node_e not in (end_node_b_2, end_node_e_2):
                        if same_contig_flag_i == True:
                            combination_i = [end_node_e, end_node_b]
                        else:
                            combination_i = [end_node_e]
                        if same_contig_flag_j == True:
                            combination_j = [end_node_b_2, end_node_e_2]
                        else:
                            combination_j = [end_node_b_2]
                        #possible_connections.add((end_node_b, end_node_b_2))
                        #possible_connections.add((end_node_b, end_node_e_2))
                        #possible_connections.add((end_node_e, end_node_e_2))
                        for comb_i in combination_i:
                            for comb_j in combination_j:
                                if (comb_i, comb_j) in filtered_edges:
                                    possible_connections.add((idx, idj, (comb_i, comb_j)))
                        #if (end_node_e, end_node_b_2) in filtered_edges:
                        #    possible_connections.add((idx, idj, (end_node_e, end_node_b_2)))
            
            merged_indices, merged_paths = self.second_pass_merge(paths, possible_connections)
            print("Merged paths: ", merged_paths)
            print("Candidate connections: ", possible_connections)
            print("--------------------------------"*2)
            new_dangling_contig_count = 0
            for new_path in merged_paths:
                if len(new_path) <= 2:
                    new_dangling_contig_count += 1
            new_paths = []
            for idx, path in enumerate(merged_paths):
                new_score = sum(paths[idy][0] for idy in merged_indices[idx])
                new_paths.append((new_score, tuple(path)))
            second_pass_result.append((score, non_end_node_count, new_dangling_contig_count, tuple(new_paths)))
        return second_pass_result

                    

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
        list_files = [self.ref_path + id + ".fasta" for id in relevant_ref_ids]
        with open(output_path + ".relevant_seq.fasta", "w") as out_f:
            for fasta in [contig_fasta_path] + list_files:
                with open(fasta, "r") as in_f:
                    content = in_f.read()
                    out_f.write(content)
                    if not content.endswith('\n'):
                        out_f.write('\n')
        return output_path + ".relevant_seq.fasta"
    
def output_agp(components, contig_fasta_path, output_path):
    pragma_line = f"## agp-version {get_ragtag_version()}"
    comment_line = f"# AGP created by impuT2T {get_impuT2T_version()}"

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
        f.write(pragma_line + '\n')
        f.write(comment_line + '\n')
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
    parser = argparse.ArgumentParser("impuT2T patching version: " + get_impuT2T_version())
    parser.add_argument("-p", "--path", required=True, help="Input path")
    parser.add_argument("-l", "--list_sample", required=True, help="file containing list of samples")
    parser.add_argument("-w", "--weight", help="weight type: allele_freq or max_alignment_score [AF|AS]", default="AS")
    parser.add_argument("-wr", "--weight_ratio", help="weight ratio, discard the edge less than ratio of the max weight the same node [0.7]", type=float, default=0.7)
    parser.add_argument("-ast", "--AS_threshold", help="AS threshold, discard the edge less than threshold [6000]", type=int, default=6000)
    parser.add_argument("-aft", "--AF_threshold", help="AF threshold, discard the edge less than threshold [0.005]", type=float, default=0.005)
    parser.add_argument("--debug", help="debug mode", action="store_true")
    parser.add_argument("--contig_fasta", required=True, help="contig fasta file, needed for unplaced contigs")
    parser.add_argument("--ref_path", help="reference path, only specify if you want the patch fasta output")
    parser.add_argument("--second_pass", help="second pass connection", action="store_true")
    parser.add_argument("--gmm", help="use gmm for edges selection", action="store_true")
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

    edge_aggregator = EdgeAggregator(args.path, samples, args.gmm)
    edges = edge_aggregator.aggregate()
    if args.debug:
        for pair, values in sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2])):
            print(pair, values)

    graph_builder = GraphBuilder(edges, weight_type, args.weight_ratio, args.AS_threshold, args.AF_threshold, args.contig_fasta, len(samples))
    graph, weights, nodes, dict_edge_info = graph_builder.build()
    edges_sorted = sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
    print("graph", "--------------------------------"*2)
    for ele in graph:
        print(ele, graph[ele])
    print("--------------------------------"*2)
    for ele in weights:
        print("weights", ele, weights[ele])
    print("--------------------------------"*2)
    for ele in edges_sorted:
        print(ele)
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
    print("############", "########################"*2)
    
    if args.second_pass:
        print("Try connecting with filtered connections in the first pass")
        graph_builder_2nd = GraphBuilder_2nd_pass(edges, weight_type, args.weight_ratio, args.AS_threshold, args.AF_threshold, args.contig_fasta, len(samples))
        graph_2nd, weights_2nd, nodes_2nd, dict_edge_info_2nd = graph_builder_2nd.convert_first_stage_results(sorted_standardized_result)
        print("graph_2nd", "--------------------------------"*2)
        for ele in graph_2nd:
            print(ele, graph_2nd[ele])
        print("--------------------------------"*2)
        for ele in weights_2nd:
            print("weights_2nd", ele, weights_2nd[ele])
        print("--------------------------------"*2)
        for ele in dict_edge_info_2nd:
            print("dict_edge_info_2nd", ele, dict_edge_info_2nd[ele])
        print("--------------------------------"*2)
        for ele in nodes_2nd:
            print("nodes_2nd", nodes_2nd)
        print("--------------------------------"*2)

        path_finder_2nd = PathFinder(graph_2nd, weights_2nd, nodes_2nd)
        end_nodes_2nd = path_finder_2nd.find_end_nodes()
        print("================================================"*2)
        print("End nodes: ", end_nodes_2nd)
        print("================================================"*2)
        if len(end_nodes_2nd) == 0:
            sorted_decompressed_result_2nd = sorted_standardized_result
            print("No end nodes found, use the first pass results")
        else:
            second_pass_results = path_finder_2nd.search_paths([], [], set(), set(), edges_sorted, end_nodes_2nd, True)
            print("TEST_RESULTS_2ND", "########################"*2)
            #second_pass_results = path_finder.second_pass_connection(sorted_standardized_result, graph_2nd, weights_2nd, nodes_2nd, dict_edge_info_2nd)
            print("SECOND_PASS_RESULTS", "########################"*2)
            print("SECOND_PASS_RESULTS", "########################"*2)

            standardized_result_2nd = set()
            for ele in second_pass_results:
                standardized_result_2nd.add(path_finder_2nd.standardize_traversals(ele, end_nodes_2nd))
            sorted_standardized_result_2nd = sorted(list(standardized_result_2nd), key=lambda x: -x[0])

            recover_components = graph_builder_2nd.recover_unused_components(sorted_standardized_result_2nd)
            recover_components = sorted(list(recover_components), key=lambda x: -x[0])
            print("Recover components: ", recover_components)
            decompressed_result_2nd = []
            for ele in recover_components:
                total_score, non_end_node_count, dangling_contig_count, paths = ele
                print(type(paths))
                new_paths = []
                for path_info in paths:
                    score, path = path_info
                    print(f"path: {path}")
                    print(f"decompressed path: {graph_builder_2nd.decompress_path(path)}")
                    new_paths.append((score, graph_builder_2nd.decompress_path(path)))
                decompressed_result_2nd.append((total_score, non_end_node_count, dangling_contig_count, tuple(new_paths)))
            sorted_decompressed_result_2nd = sorted(decompressed_result_2nd, key=lambda x: -x[0])
            for ele in sorted_decompressed_result_2nd:
                print(ele)
            print("############", "########################"*2)
        
    
    components = []
    agp_assembler = AGPAssembler(args.path)
    total_score_1st, non_end_node_count_1st, dangling_contig_count_1st, paths = sorted_decompressed_result_2nd[0]
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
        relevant_samples = sorted(list(relevant_samples))
        
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
