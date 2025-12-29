import argparse
import os
import copy
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Set, Tuple, Dict
from combine_connections import read_agp_file, merge_segments, reverse_segments
from ragtag_utilities.utilities import get_ragtag_version, get_impuT2T_version
import ragtag_agp2fa 

os.environ["MKL_THREADING_LAYER"] = "GNU"
import pysam
import numpy as np
from sklearn.mixture import GaussianMixture
import sys

import time

# ==================== Constants ====================

# Simple, frequently-used constants at module level
DEFAULT_TIMEOUT = 300
RETRY_TIMEOUT = 300
CONNECTION_SCORE_MULTIPLIER = 10

# Grouped constants in classes (for organization)
class GMMConfig:
    """GMM-related constants for edge selection"""
    MAX_COMPONENTS = 10
    WEIGHT_CUTOFF = 0.1
    MIN_SAMPLES = 10
    STD_MULTIPLIER = 1.5
    LOW_CONFIDENCE_SCORE = 3000
    LARGE_STD_THRESHOLD = 200000  # bp
    CLOSE_TO_ZERO_DISTANCE = 100000  # bp
    SCORE_RATIO_THRESHOLD = 0.5
    TOP_SCORES_COUNT = 10

class GraphConfig:
    """Graph building constants"""
    DEFAULT_FILTER_NUM = 200
    RETRY_FILTER_NUM_FIRST_PASS = 20
    RETRY_FILTER_NUM_SECOND_PASS = 40

class TimeoutException(Exception):
    pass



def parse_sample_list(list_sample: str) -> list[str]:
    with open(list_sample) as f:
        return [line.strip() for line in f if line.strip()]


# ==================== Node Utility Functions ====================
class NodeUtils:
    """Utility functions for working with graph nodes (contig_b, contig_e format)."""
    
    @staticmethod
    def get_contig(node: str) -> str:
        """Extract contig name from node (removes _b or _e suffix)."""
        return node[:-2]
    
    @staticmethod
    def is_begin(node: str) -> bool:
        """Check if node is a begin node (ends with _b)."""
        return node.endswith("_b")
    
    @staticmethod
    def get_opposite(node: str) -> str:
        """Get the opposite end of a contig node."""
        contig = NodeUtils.get_contig(node)
        return f"{contig}_e" if NodeUtils.is_begin(node) else f"{contig}_b"
    

# ==================== Graph Builder ====================
@dataclass
class GraphBuildConfig:
    """Configuration for building graphs"""
    edges: dict
    weight_type: int
    weight_ratio: float
    contig_fasta: str
    sample_size: int
    filter_num: int = 20
    first_pass_results: list[list[str]] = None
    first_pass_scores: list[float] = None
    
    def build_graph(self):
        """Build a graph with these parameters"""
        graph_builder = GraphBuilder(self.edges, self.weight_type, self.weight_ratio, self.contig_fasta, self.sample_size)
        if self.first_pass_results is not None and self.first_pass_scores is not None:
            return graph_builder.build_with_results(self.first_pass_results, self.first_pass_scores, filter_num=self.filter_num)
        else:
            return graph_builder.build(filter_num=self.filter_num)

# ==================== Edge Aggregator ====================
class EdgeAggregator:
    def __init__(self, input_path: str, samples: list[str], gmm: bool):
        self.input_path = input_path
        self.samples = samples
        self.gmm = gmm
        self.raw_edges = {}
        self.edges = {}

    def fit_gmm(self, x, max_components=GMMConfig.MAX_COMPONENTS, cutoff=GMMConfig.WEIGHT_CUTOFF):
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
    
    def _extract_gmm_candidates(self, means, stds, dists, scores, relevant_samples):
        """Extract best candidate for each GMM component within bounds."""
        gmm_scores, gmm_samples, gmm_dists = [], [], []
        
        for mean, std in zip(means, stds):
            upper_bound = mean + GMMConfig.STD_MULTIPLIER * std
            lower_bound = mean - GMMConfig.STD_MULTIPLIER * std
            if lower_bound < 0 and mean > 0:
                lower_bound = mean - std
            
            # Find candidates within bounds
            candidate_dists = []
            candidate_scores = []
            candidate_samples = []
            for idx in range(len(dists)):
                if lower_bound <= dists[idx] <= upper_bound:
                    candidate_dists.append(dists[idx])
                    candidate_scores.append(scores[idx])
                    candidate_samples.append(relevant_samples[idx])
            
            if len(candidate_dists) > 0:
                best_close_idx = np.argmax(candidate_scores)
                gmm_scores.append(candidate_scores[best_close_idx])
                gmm_samples.append(candidate_samples[best_close_idx])
                gmm_dists.append(candidate_dists[best_close_idx])
        
        return gmm_scores, gmm_samples, gmm_dists
    
    def _select_with_fallback(self, gmm_scores, gmm_samples, gmm_dists, stds,
                              scores, dists, relevant_samples):
        """
        Select best candidate, applying fallback strategies if low confidence.
        Returns (selected_idx, selected_score, selected_sample, selected_distance)
        """
        # Initial selection from GMM results
        final_selected_idx = np.argmax(gmm_scores)
        max_score = max(gmm_scores)
        
        # High confidence - return immediately
        if (max_score >= GMMConfig.LOW_CONFIDENCE_SCORE or 
            stds[final_selected_idx] <= GMMConfig.LARGE_STD_THRESHOLD):
            selected_sample = gmm_samples[final_selected_idx]
            selected_distance = gmm_dists[final_selected_idx]
            print(f"Selected sample: {selected_sample}, distance: {selected_distance}")
            return (final_selected_idx, gmm_scores[final_selected_idx], 
                   selected_sample, selected_distance)
        
        # Low confidence - try fallback strategies
        print(f"std: {stds[final_selected_idx]} is too large, try others")
        
        # Strategy 1: Smaller std components
        smaller_scores = []
        smaller_idx = []
        for gmm_idx, std in enumerate(stds):
            if std < GMMConfig.LARGE_STD_THRESHOLD:
                smaller_scores.append(gmm_scores[gmm_idx])
                smaller_idx.append(gmm_idx)
        
        if len(smaller_scores) > 0:
            final_selected_idx = smaller_idx[np.argmax(smaller_scores)]
            selected_sample = gmm_samples[final_selected_idx]
            selected_distance = gmm_dists[final_selected_idx]
            print(f"Selected sample: {selected_sample}, distance: {selected_distance}")
            return (final_selected_idx, gmm_scores[final_selected_idx], 
                   selected_sample, selected_distance)
        
        # Strategy 2: Close to zero distance
        print(f"NO SMALLER STD FOUND, FORCE TO SELECT THE CLOSEST TO 0!!!!")
        small_max_score = 0
        small_idx = -1
        for idx_dist, dist in enumerate(dists):
            if abs(dist) < GMMConfig.CLOSE_TO_ZERO_DISTANCE:
                if scores[idx_dist] > small_max_score:
                    small_max_score = scores[idx_dist]
                    small_idx = idx_dist
        
        if small_idx != -1 and small_max_score > max_score * GMMConfig.SCORE_RATIO_THRESHOLD:
            final_selected_idx = small_idx
            max_score = small_max_score
            print(f"Selected sample: {relevant_samples[final_selected_idx]}, distance: {dists[final_selected_idx]}")
            return (final_selected_idx, max_score, 
                   relevant_samples[final_selected_idx], dists[final_selected_idx])
        
        # Strategy 3: Top scores with smallest distance
        print(f"NO SMALLER SCORE FOUND, FORCE TO SELECT THE CLOSEST TO 0!!!!")
        top_10_score = sorted(scores, reverse=True)[:GMMConfig.TOP_SCORES_COUNT]
        print(f"Top 10 score: {top_10_score})")
        small_max_dist = gmm_dists[final_selected_idx]
        small_idx = -1
        for idx_score, score in enumerate(scores):
            if score in top_10_score:
                if dists[idx_score] < small_max_dist:
                    small_max_dist = dists[idx_score]
                    small_idx = idx_score
        
        if small_idx != -1:
            final_selected_idx = small_idx
            max_score = scores[final_selected_idx]
            print(f"Selected sample: {relevant_samples[final_selected_idx]}, distance: {dists[final_selected_idx]}")
            return (final_selected_idx, max_score, 
                   relevant_samples[final_selected_idx], dists[final_selected_idx])
        
        # Final fallback: Use original GMM selection
        selected_sample = gmm_samples[final_selected_idx]
        selected_distance = gmm_dists[final_selected_idx]
        print(f"Selected sample: {selected_sample}, distance: {selected_distance}")
        return (final_selected_idx, gmm_scores[final_selected_idx], 
               selected_sample, selected_distance)
    
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
                if len(relevant_samples) > GMMConfig.MIN_SAMPLES:
                    # Fit GMM model
                    gmm, best_k = self.fit_gmm(np.array(dists))
                    means = gmm.means_.flatten()
                    stds = np.sqrt(gmm.covariances_.flatten())  # for 1D case
                    print(f"====================================== Best k: {best_k} ======================================")
                    print(f"======== pair: {pair} ========")
                    print(f"GMM means: {means}, stds: {stds}")
                    
                    # Extract candidates from GMM components
                    gmm_scores, gmm_samples, gmm_dists = self._extract_gmm_candidates(
                        means, stds, dists, scores, relevant_samples
                    )
                    
                    # Select best candidate with fallback if needed
                    selected_idx, selected_score, selected_sample, selected_dist = \
                        self._select_with_fallback(
                            gmm_scores, gmm_samples, gmm_dists, stds,
                            scores, dists, relevant_samples
                        )
                    
                    self.edges[pair] = [
                        len(relevant_samples), selected_score, 
                        selected_sample, selected_dist
                    ]
                else:
                    idx = np.argmax(scores)
                    self.edges[pair] = [len(relevant_samples), scores[idx], relevant_samples[idx], dists[idx]]
        else:
            for pair, info in self.raw_edges.items():
                idx = np.argmax(info[1], key=lambda x: x[0])
                self.edges[pair] = [len(info[0]), info[1][idx][0], info[0][idx]]

        return self.edges

class GraphBuilder:
    def __init__(self, edges, weight_type, weight_ratio, contig_fasta, sample_size):
        self.edges = edges
        self.weight_type = weight_type
        self.weight_ratio = weight_ratio
        self.sample_size = sample_size
        self.graph = defaultdict(list)
        self.weights = {}
        self.nodes = set()
        self.dict_node_score = defaultdict(list)
        self.dict_edge_info = {}
        self.dict_node_outgoing_edges = defaultdict(list)
        self.filtered_edges = set()
        fai = pysam.FastaFile(contig_fasta)
        self.contig_lengths = {contig: fai.get_reference_length(contig) for contig in fai.references}
        self.contig_avg_length = sum(self.contig_lengths.values()) / len(self.contig_lengths)
        self.dict_compressed_2_component = {}

    def _get_filter_function(self, filter_type: str):
        """Get filter function for sorting edges."""
        if filter_type == "AF":
            return lambda x: x[1][0]
        elif filter_type == "AS":
            return lambda x: x[1][1]
        else:
            raise ValueError(f"Invalid filter type: {filter_type}")
    
    def _build_node_score_and_edge_info(self, edges_to_process):
        """
        Build dict_node_score and dict_edge_info from edges.
        edges_to_process: list of ((u, v), values) tuples
        """
        for pair, values in edges_to_process:
            u, v = pair
            w = values[self.weight_type]
            self.dict_node_score[u].append((w, v))
            self.dict_node_score[v].append((w, u))
            self.dict_edge_info[(u, v)] = values
            self.dict_edge_info[(v, u)] = values
    
    def _add_edge_to_graph(self, u, v, w):
        """
        Add an edge to the graph with all necessary updates.
        """
        # Undirected graph assignment
        self.graph[u].append(v)
        self.graph[v].append(u)
        self.weights[(u, v)] = w
        self.weights[(v, u)] = w
        self.nodes.update([u, v, NodeUtils.get_opposite(u), NodeUtils.get_opposite(v)])
        self.dict_node_outgoing_edges[u].append((v, w))
        self.dict_node_outgoing_edges[v].append((u, w))
    
    def _add_internal_edges(self):
        """
        Add edges connecting begin and end nodes of the same contig.
        """
        for node in self.nodes:
            if node.endswith("_b"):
                contig = node[:-2]
                e_node = f"{contig}_e"
                self.graph[node].append(e_node)
            elif node.endswith("_e"):
                contig = node[:-2]
                b_node = f"{contig}_b"
                self.graph[node].append(b_node)
    
    def _add_contig_weights(self, used_contigs=None):
        """
        Add weights for traversing within contigs (b->e, e->b).
        used_contigs: set of contigs to skip (for second pass)
        """
        for contig in self.contig_lengths:
            if used_contigs and contig in used_contigs:
                print(f"contig: {contig} is used in the first stage")
                continue
            b_node = f"{contig}_b"
            e_node = f"{contig}_e"
            if used_contigs:
                # Second pass: set to 0 for used contigs
                self.weights[(b_node, e_node)] = 0
                self.weights[(e_node, b_node)] = 0
            else:
                # First pass: use normalized length
                self.weights[(b_node, e_node)] = self.contig_lengths[contig] / self.contig_avg_length
                self.weights[(e_node, b_node)] = self.contig_lengths[contig] / self.contig_avg_length
    
    def _apply_local_edge_filtering(self):
        """
        Apply local edge filtering based on weight ratio.
        """
        for node, scores in self.dict_node_score.items():
            max_score = max(scores, key=lambda x: x[0])[0]
            for score, neighbor in scores:
                if score < max_score * self.weight_ratio:
                    self.filtered_edges.add((node, neighbor))
                    self.filtered_edges.add((neighbor, node))

    def _build_compressed_nodes(self, connected_components, connected_components_scores):
        """Build compressed nodes from connected components (second pass only)."""
        dict_node_2_compressed = {}
        self.dict_compressed_2_component = {}
        set_used_nodes = set()
        set_used_contigs = set()
        
        for idx, component in enumerate(connected_components):
            if len(component) == 1:
                continue
            average_order = sum([int(pair[0].split("_")[-2]) for pair in component]) / len(component)
            average_order += 0.01 * int(component[0][0].split("_")[-2])
            compressed_name = f"impuT2T_{average_order:.2f}"
            self.dict_compressed_2_component[compressed_name] = component
            u = compressed_name + '_b'
            v = compressed_name + '_e'
            dict_node_2_compressed[component[0][0]] = u
            dict_node_2_compressed[component[-1][1]] = v
            for idy in range(len(component)-1):
                set_used_nodes.update([component[idy][1], component[idy+1][0]])
            set_used_contigs.update([pair[0][:-2] for pair in component])
            
            self.nodes.update([u, v, NodeUtils.get_opposite(u), NodeUtils.get_opposite(v)])
            self.weights[(u, v)] = connected_components_scores[idx]
            self.weights[(v, u)] = connected_components_scores[idx]
        
        print("dict_node_2_compressed", dict_node_2_compressed)
        print("set_used_contigs", set_used_contigs)
        print("set_used_nodes", set_used_nodes)
        print("--------------------------------"*2)
        
        return dict_node_2_compressed, set_used_nodes, set_used_contigs
    
    def _compress_edges(self, edges_sorted, dict_node_2_compressed, set_used_nodes, filter_num):
        """Compress edges based on compressed nodes (second pass only)."""
        compressed_edges = []
        add_counter = 0
        
        for pair, values in edges_sorted:
            if add_counter > filter_num:
                self.filtered_edges.add(pair)
                self.filtered_edges.add(pair[::-1])
                break
            
            u, v = pair
            # Skip if one of the nodes is internal compressed node
            if u in set_used_nodes or v in set_used_nodes:
                continue
            elif u in dict_node_2_compressed and v in dict_node_2_compressed:
                new_u = dict_node_2_compressed[u]
                new_v = dict_node_2_compressed[v]
                compressed_edges.append(((new_u, new_v), values))
            elif u in dict_node_2_compressed:
                new_u = dict_node_2_compressed[u]
                compressed_edges.append(((new_u, v), values))
            elif v in dict_node_2_compressed:
                new_v = dict_node_2_compressed[v]
                compressed_edges.append(((u, new_v), values))
            else:
                compressed_edges.append(((u, v), values))
            add_counter += 1
        
        return compressed_edges

    def build(self, filter_num: int = GraphConfig.DEFAULT_FILTER_NUM, filter_type: str = "AF"):
        filter_func = self._get_filter_function(filter_type)
        edges_sorted = sorted(self.edges.items(), key=filter_func, reverse=True)
        
        # Build node scores and edge info
        self._build_node_score_and_edge_info(edges_sorted)
        
        # Apply local edge filtering
        self._apply_local_edge_filtering()
        
        # Process edges and add to graph
        add_counter = 0 # counter for filtering edges
        for pair, values in edges_sorted:
            u, v = pair
            AF, AS, _, edge_dist = values
            w = AF / (2 * self.sample_size)
            
            if (u, v) in self.filtered_edges or (v, u) in self.filtered_edges or add_counter > filter_num:
                self.filtered_edges.add((u, v))
                self.filtered_edges.add((v, u))
                continue

            add_counter += 1
            print("\t", (u, v), values)
            self._add_edge_to_graph(u, v, w)
        
        # Add contig edges and weights
        self._add_internal_edges()
        self._add_contig_weights()
        
        return self.graph, self.weights, self.nodes, self.dict_edge_info, self.dict_node_outgoing_edges
        
    def build_with_results(self, connected_components, connected_components_scores, filter_num: int = GraphConfig.DEFAULT_FILTER_NUM, filter_type: str = "AF"):
        filter_func = self._get_filter_function(filter_type)
        
        # Build compressed nodes from connected components (unique to second pass)
        dict_node_2_compressed, set_used_nodes, set_used_contigs = \
            self._build_compressed_nodes(connected_components, connected_components_scores)
        
        # Process and compress edges (unique to second pass)
        edges_sorted = sorted(self.edges.items(), key=filter_func, reverse=True)
        compressed_edges = self._compress_edges(edges_sorted, dict_node_2_compressed, 
                                                 set_used_nodes, filter_num)
        
        # Build node scores and edge info from compressed edges
        self._build_node_score_and_edge_info(compressed_edges)
        
        # Process compressed edges and add to graph
        for pair, values in compressed_edges:
            u, v = pair
            AF, AS, *_ = values
            w = AF / (2 * self.sample_size)
            
            print("\t", (u, v), values)
            self._add_edge_to_graph(u, v, w)
        
        # Add contig edges and weights (with used_contigs)
        self._add_internal_edges()
        self._add_contig_weights(used_contigs=set_used_contigs)
        
        return self.graph, self.weights, self.nodes, self.dict_edge_info, self.dict_node_outgoing_edges

    def decompress_component(self, component):
        """
        Decompress the component by replacing the compressed nodes with the original nodes
        """
        def invert_path(path: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
            new_path = []
            for edge in path[::-1]:
                new_path.append((edge[1], edge[0]))
            return new_path

        new_component = []
        for path in component:
            new_path = []
            for pair in path:
                if pair[0][:-2] in self.dict_compressed_2_component:
                    if pair[0].endswith("_b"):
                        new_path += self.dict_compressed_2_component[pair[0][:-2]]
                    else:
                        new_path += invert_path( self.dict_compressed_2_component[pair[0][:-2]] )
                else:
                    new_path.append(pair)
            new_component.append(new_path)
        return new_component



class PathFinder:
    def __init__(self, graph, weights, nodes, dict_node_outgoing_edges):
        self.graph = graph
        self.weights = weights
        self.nodes = nodes
        self.dict_node_outgoing_edges = dict_node_outgoing_edges

        # variables for search_paths
        self.num_nodes = len(nodes)

        #self.searched_paths = set()
        self.best_score_component = {'score':0, 'length':self.num_nodes, 'component':[]}
        self.best_connected_component = {'score':0, 'length':self.num_nodes, 'component':[]}
        
        # Timeout management
        self.search_start_time = None
        self.search_timeout_seconds = None

    def find_end_nodes_unvisited(self, visited_nodes: Set[str]) -> Set[str]:
        remaining_nodes = self.nodes - visited_nodes

        end_nodes = set()
        for node in remaining_nodes:
            list_outgoing_edges = self.dict_node_outgoing_edges[node]
            flag_end_node = True
            for edge in list_outgoing_edges:
                if edge[0] not in visited_nodes and edge[1] not in visited_nodes:
                    flag_end_node = False
                    break
            if flag_end_node:
                end_nodes.add(node)
        
        if end_nodes:
            return end_nodes
        else:
            return remaining_nodes

    def search_paths_DFS(self, current_path: List[str], current_node: str, connected_components: List[List[str]],
                    visited_nodes: Set[str], visited_contigs: Set[str], depth: int) -> None:
        """
        DFS search for all possible paths.
        1. Check if all the nodes are visited
            1.1 all visited -> calculate and record the best scored path
            1.2 keep visited until all the nodes are visited
        2. Tunneling to the other end of current node, then traverse all outgoing edges for the other end
            2.1 if no outgoing edges -> find new starting node in the remaining unvisited nodes
        """
        
        # Check timeout at the start of each recursive call
        if self.search_timeout_seconds is not None and self.search_start_time is not None:
            elapsed_time = time.time() - self.search_start_time
            if elapsed_time > self.search_timeout_seconds:
                raise TimeoutException(f"Search timeout after {elapsed_time:.2f} seconds")

        # All nodes are visited -> calculate and record only the best scored path
        if self.num_nodes == len(visited_nodes):
            # CRITICAL: Append a copy of current_path, not a reference
            # Shallow copy is sufficient since tuples are immutable and current_path is never modified in place
            connected_components.append(list(current_path))
            #print("Traversal completed!")
            score, length, _ = self.component_score(connected_components)
            if score > self.best_score_component['score']:
                # CRITICAL: Store a deep copy, not a reference (connected_components will continue to be modified)
                self.best_score_component = {
                    'score': score, 
                    'length': length, 
                    'component': copy.deepcopy(connected_components)
                }
                print("best_score_component updated: ", self.best_score_component)
                print("--------------------------------"*2)
            if length <= self.best_connected_component['length']:
                if score > self.best_connected_component['score']:
                    # CRITICAL: Store a deep copy, not a reference
                    self.best_connected_component = {
                        'score': score, 
                        'length': length, 
                        'component': copy.deepcopy(connected_components)
                    }
                    print("best_connected_component updated: ", self.best_connected_component)
                    print("--------------------------------"*2)
            return

        # Go to the other end of the same contig / connected contigs
        opposite_node = NodeUtils.get_opposite(current_node)
        outgoing_edges = self.dict_node_outgoing_edges[opposite_node]
        unvisited_edges = [edge for edge in outgoing_edges if edge[0] not in visited_nodes and edge[1] not in visited_nodes]

        
        if len(unvisited_edges) == 0:
            #print("Reach End!") #, current_path)
            # Reaching the end of this path, start a new path from the remaining end nodes
            # CRITICAL: Append a copy of current_path, not a reference
            # Shallow copy is sufficient since tuples are immutable and current_path is never modified in place
            connected_components.append(list(current_path))
            remaining_end_nodes = self.find_end_nodes_unvisited(visited_nodes)
            for end_node in remaining_end_nodes:
                #print("---" * depth + '||', len(connected_components), "| end_node: ", end_node)
                new_path = [(end_node, NodeUtils.get_opposite(end_node))]
                new_visited_nodes = visited_nodes.union([end_node, NodeUtils.get_opposite(end_node)])
                new_visited_contigs = visited_contigs.union([NodeUtils.get_contig(end_node)])
                # CRITICAL: Use deep copy - shallow copy only copies the outer list, inner lists are still references
                self.search_paths_DFS(new_path, end_node, copy.deepcopy(connected_components), new_visited_nodes, new_visited_contigs, depth + 1)
            #print("---" * depth + '||', current_node)
            return

        # Keep traversing the outgoing edges
        for outgoing_edge in unvisited_edges:
            new_path = current_path + [(outgoing_edge[0], NodeUtils.get_opposite(outgoing_edge[0]))]
            new_visited_nodes = visited_nodes.union([outgoing_edge[0], NodeUtils.get_opposite(outgoing_edge[0])])
            new_visited_contigs = visited_contigs.union([NodeUtils.get_contig(outgoing_edge[0])])
            # CRITICAL: Must copy connected_components here! All recursive calls share the same list otherwise
            # When one call appends to connected_components (line 739), all other calls see that change
            self.search_paths_DFS(new_path, outgoing_edge[0], copy.deepcopy(connected_components), new_visited_nodes, new_visited_contigs, depth + 1)        
        #print("---" * depth + '|', current_node)
        return

    def component_score(self, component: List[List[str]]) -> float:
        """
        Calculate individual score of each path, use HHI with power 1.5 to sum up the total score
        """
        scores = []
        for path in component:
            contig_score = 0.0
            #print("path: ", path)
            for edge in path:
                #print("edge: ", edge, "weight: ", self.weights.get((edge[0], edge[1]), 0.0))
                contig_score += self.weights.get((edge[0], edge[1]), 0.0)
            connection_score = 0.0
            for idx in range(len(path) - 1):
                #print("connection: ", path[idx][1], path[idx+1][0], "weight: ", self.weights.get((path[idx][1], path[idx+1][0]), 0.0))
                connection_score += self.weights.get((path[idx][1], path[idx+1][0]), 0.0) * CONNECTION_SCORE_MULTIPLIER
                
            path_score = contig_score + connection_score
            scores.append(path_score)
            
            #scores.append(connection_score)
        HHI_score = sum([score**1.5 for score in scores])
        return (int(HHI_score), len(component), scores)

    def invert_path(self, path: List[str]) -> List[str]:
        """
        Invert the path by reversing the order of the each edge and internal nodes
        """
        new_path = []
        for edge in path[::-1]:
            new_path.append((edge[1], edge[0]))
        return new_path

    def smaller_than(self, node_1: str, node_2: str) -> bool:
        return float(node_1.split("_")[-2]) < float(node_2.split("_")[-2])

    def sort_component(self, component: List[List[str]]) -> List[List[str]]:
        sorted_component = []
        for path in component:
            if self.smaller_than(path[0][0], path[-1][0]):
                sorted_component.append(path)
            else:
                sorted_component.append(self.invert_path(path))
        sorted_component = sorted(sorted_component, key=lambda x: float(x[0][0].split("_")[-2]))
        return sorted_component
        

class AGPAssembler:
    def __init__(self, agp_dir, contig_fasta):
        self.agp_dir = agp_dir
        fai = pysam.FastaFile(contig_fasta)
        self.contig_lengths = {contig: fai.get_reference_length(contig) for contig in fai.references}
        fai.close()

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
            
            # Check if all entries belong to the same AGP block
            if s <= e:
                candidate_entries = agp_entries[s:e+1]
            else:
                candidate_entries = agp_entries[e:s+1]
            
            # Verify all entries have the same object ID (first column)
            object_ids = [entry['object'] for entry in candidate_entries]
            if len(set(object_ids)) != 1:
                continue  # Skip this pair if entries don't belong to the same block
            
            if s <= e:
                if check_orientation(strand_s, strand_e, start_end, end_end, True):
                    return candidate_entries, True
            else:
                if check_orientation(strand_s, strand_e, start_end, end_end, False):
                    return candidate_entries, False
        print("No valid connection found!!")
        return [], False

    def _calculate_trim(self, segment, is_left, other_orientation):
        """Calculate trim value for a segment based on orientation, position, and other segment's orientation."""
        contig_len = self.contig_lengths[segment['component_id']]
        if segment['orientation'] == '+':
            # For + orientation: trim is overhang at the right end (left segment) or left end (right segment)
            if is_left:
                return contig_len - segment['component_end']
            else:
                return segment['component_beg']
        else:
            if is_left:
                return segment['component_beg']
            else:
                return contig_len - segment['component_end']
    
    def _apply_adjustment(self, segment, offset, is_left, left_orientation, right_orientation):
        """Apply offset adjustment to segment. Direction depends on both segments' orientations."""
        # If orientations differ, both segments adjust in the direction of the left segment's orientation
        if left_orientation != right_orientation:
            # Mixed orientations: both extend in the direction of the left segment's orientation
            if left_orientation == '+':
                segment['component_end'] += offset
                if is_left:
                    segment['object_end'] += offset
                else:
                    segment['object_beg'] -= offset
            else:  # left is '-'
                segment['component_beg'] -= offset
                if is_left:
                    segment['object_end'] += offset
                else:
                    segment['object_beg'] -= offset
        else:
            # Same orientations: left extends in its orientation direction, right extends opposite
            if is_left:
                # Left segment: extend in direction of its orientation
                if segment['orientation'] == '+':
                    segment['component_end'] += offset
                    segment['object_end'] += offset
                else:  # orientation == '-'
                    segment['component_beg'] -= offset
                    segment['object_end'] += offset
            else:
                # Right segment: extend opposite to its orientation (toward the overlap)
                if segment['orientation'] == '+':
                    segment['component_beg'] -= offset
                    segment['object_beg'] -= offset
                else:  # orientation == '-'
                    segment['component_end'] += offset
                    segment['object_beg'] -= offset
    
    def middle_connection_for_overlap_segments(self, segments):
        if len(segments) == 3: # if it's not overlapping, there are patch segments in the middle
            return segments
        elif len(segments) == 2:
            print("--------------------------------"*2)
            print("Overlapping segments found!!")
            for idx, segment in enumerate(segments):
                print(f"Segment {idx}: {segment}")
            print("--------------------------------"*2)
            segment_left, segment_right = segments
            
            # Calculate trims based on orientations
            left_trim = self._calculate_trim(segment_left, is_left=True, other_orientation=segment_right['orientation'])
            right_trim = self._calculate_trim(segment_right, is_left=False, other_orientation=segment_left['orientation'])
            
            # Calculate middle point and offsets
            middle_trim = int((left_trim + right_trim) / 2)
            offset_left = left_trim - middle_trim
            offset_right = right_trim - middle_trim

            # Apply adjustments
            left_orient = segment_left['orientation']
            right_orient = segment_right['orientation']
            self._apply_adjustment(segment_left, offset_left, is_left=True, 
                                 left_orientation=left_orient, right_orientation=right_orient)
            self._apply_adjustment(segment_right, -offset_left, is_left=False, 
                                 left_orientation=left_orient, right_orientation=right_orient)
            
            return [segment_left, segment_right]
        print("Invalid number of segments!!")
        return []

    def extract_and_merge(self, path, dict_edge_info, agp_path):
        all_segments = []
        for idx in range(len(path)-1):
            u, v = path[idx][1], path[idx+1][0]
            edge_info = dict_edge_info[(u, v)]
            sample = edge_info[2]
            agp_file = agp_path + f"{sample}.agp"
            print(f"AGP connection for {u} -> {v} in sample {sample}:")
            segments, flag_forward = self.find_connection_segments_in_agp(agp_file, u, v)
            segments = self.middle_connection_for_overlap_segments(segments)
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

def run_path_search_with_retry(
    path_finder: PathFinder,
    graph_config: GraphBuildConfig,
    initial_timeout: int = 300,
    retry_timeout: int = 300,
    ) -> PathFinder:
    """
    Run path search with timeout handling and automatic retry.
    Args:
        path_finder: PathFinder instance to use
        initial_timeout: Initial timeout in seconds (default 300)
        retry_timeout: Timeout for retry attempt in seconds (default 300)

    Returns:
        path_finder if search completed successfully, partially searched path_finder if timed out
    """
    end_nodes = path_finder.find_end_nodes_unvisited(visited_nodes=set())

    # Initial attempt
    path_finder.search_start_time = time.time()
    path_finder.search_timeout_seconds = initial_timeout
    for end_node in end_nodes:
        try:
            path_finder.search_paths_DFS(
                current_path=[(end_node, NodeUtils.get_opposite(end_node))],
                current_node=end_node,
                connected_components=[],
                visited_nodes={end_node, NodeUtils.get_opposite(end_node)},
                visited_contigs={NodeUtils.get_opposite(end_node)},
                depth=0
            )
            return path_finder
        except TimeoutException as e:
            elapsed = time.time() - path_finder.search_start_time
            print(f"Timeout occurred after {elapsed:.2f} seconds: {e}")
            # Retry with filtered graph
            print("Restarting with a filtered graph...")
            graph, weights, nodes, _, dict_node_outgoing_edges = graph_config.build_graph()
            path_finder = PathFinder(graph, weights, nodes, dict_node_outgoing_edges)
            path_finder.search_start_time = time.time()
            path_finder.search_timeout_seconds = retry_timeout
            #end_nodes = path_finder.find_end_nodes_unvisited(visited_nodes=set())
            for end_node in end_nodes:
                try:
                    path_finder.search_paths_DFS(
                        current_path=[(end_node, NodeUtils.get_opposite(end_node))],
                        current_node=end_node,
                        connected_components=[],
                        visited_nodes={end_node, NodeUtils.get_opposite(end_node)},
                        visited_contigs={NodeUtils.get_contig(end_node)},
                        depth=0
                    )
                    return path_finder
                except TimeoutException as e2:
                    elapsed2 = time.time() - path_finder.search_start_time
                    print(f"Retry timeout occurred after {elapsed2:.2f} seconds: {e2}")
                    return path_finder


def main():                        
    parser = argparse.ArgumentParser("impuT2T patching version: " + get_impuT2T_version())
    parser.add_argument("-p", "--path", required=True, help="Input path")
    parser.add_argument("-l", "--list_sample", required=True, help="file containing list of samples")
    parser.add_argument("-w", "--weight", help="weight type: allele_freq or max_alignment_score [AF|AS]", default="AS")
    parser.add_argument("-wr", "--weight_ratio", help="weight ratio, discard the edge less than ratio of the max weight the same node [0.7]", type=float, default=0.7)
    #parser.add_argument("-ast", "--AS_threshold", help="AS threshold, discard the edge less than threshold [6000]", type=int, default=6000)
    #parser.add_argument("-aft", "--AF_threshold", help="AF threshold, discard the edge less than threshold [0.005]", type=float, default=0.005)
    parser.add_argument("--timeout", help=f"timeout for full path traversal [{DEFAULT_TIMEOUT}(s)]", type=int, default=DEFAULT_TIMEOUT)
    parser.add_argument("--retry_timeout", help=f"timeout for retry attempt [{RETRY_TIMEOUT}(s)]", type=int, default=RETRY_TIMEOUT)
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

    graph_builder = GraphBuilder(edges, weight_type, args.weight_ratio, args.contig_fasta, len(samples))
    graph, weights, nodes, dict_edge_info, dict_node_outgoing_edges = graph_builder.build()
    edges_sorted = sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
    """
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
    """

    path_finder = PathFinder(graph, weights, nodes, dict_node_outgoing_edges)
    graph_config = GraphBuildConfig(
        edges=edges,
        weight_type=weight_type,
        weight_ratio=args.weight_ratio,
        contig_fasta=args.contig_fasta,
        sample_size=len(samples),
        filter_num=GraphConfig.RETRY_FILTER_NUM_FIRST_PASS
    )
    path_finder = run_path_search_with_retry(path_finder, graph_config, args.timeout, args.retry_timeout)

    
    first_pass_results = path_finder.best_score_component['component']
    first_pass_results = path_finder.sort_component(first_pass_results)
    _, _, first_pass_scores = path_finder.component_score(first_pass_results)
    #first_pass_results_connected = path_finder.sort_component(first_pass_results_connected)
    print(first_pass_results)
    print(first_pass_scores)

    #print(first_pass_results_connected)
    print("##########################################"*2)
    
    if args.second_pass and len(first_pass_results) > 1:
        print("Try connect with previously filtered edges")
        graph_builder_2nd_pass = GraphBuilder(edges, weight_type, args.weight_ratio, args.contig_fasta, len(samples))
        graph_2nd_pass, weights_2nd_pass, nodes_2nd_pass, _, dict_node_outgoing_edges_2nd_pass = graph_builder_2nd_pass.build_with_results(first_pass_results, first_pass_scores)
        """
        print("graph_2nd_pass", "--------------------------------"*2)
        for ele in graph_2nd_pass:
            print(ele, graph_2nd_pass[ele])
        print("--------------------------------"*2)
        for key in weights_2nd_pass:
            print("weights_2nd_pass", key, weights_2nd_pass[key])
        print("--------------------------------"*2)
        for key in dict_node_outgoing_edges_2nd_pass:
            print("dict_node_outgoing_edges_2nd_pass", key, dict_node_outgoing_edges_2nd_pass[key])
        print("--------------------------------"*2)
        for ele in nodes_2nd_pass:
            print("nodes_2nd_pass", ele)
        print("--------------------------------"*2)
        for ele in graph_builder_2nd_pass.dict_compressed_2_component:
            print("dict_compressed_2_component", ele, graph_builder_2nd_pass.dict_compressed_2_component[ele])
        print("--------------------------------"*2)
        """

        path_finder_2nd = PathFinder(graph_2nd_pass, weights_2nd_pass, nodes_2nd_pass, dict_node_outgoing_edges_2nd_pass)
        graph_config_2nd = GraphBuildConfig(
            edges=edges,
            weight_type=weight_type,
            weight_ratio=args.weight_ratio,
            contig_fasta=args.contig_fasta,
            sample_size=len(samples),
            filter_num=GraphConfig.RETRY_FILTER_NUM_SECOND_PASS,
            first_pass_results=first_pass_results,
            first_pass_scores=first_pass_scores
        )
        path_finder_2nd = run_path_search_with_retry(path_finder_2nd, graph_config_2nd, args.timeout, args.retry_timeout)

        print("--------------------------------"*2)
        print("dict_compressed_2_component", graph_builder_2nd_pass.dict_compressed_2_component)
        print("--------------------------------"*2)
        second_pass_results = path_finder_2nd.best_score_component['component']
        _, _, second_pass_scores = path_finder_2nd.component_score(second_pass_results)
        second_pass_results = graph_builder_2nd_pass.decompress_component(second_pass_results)
        second_pass_results = path_finder_2nd.sort_component(second_pass_results)
        print(second_pass_results)
        print(second_pass_scores)

        components = []
        agp_assembler = AGPAssembler(args.path, args.contig_fasta)
        for path in second_pass_results:
            all_segments = agp_assembler.extract_and_merge(path, dict_edge_info, args.path)
            components.append(all_segments)

    else:
        # Process the first pass results without second pass
        components = []
        agp_assembler = AGPAssembler(args.path, args.contig_fasta)
        for path in first_pass_results:
            all_segments = agp_assembler.extract_and_merge(path, dict_edge_info, args.path)
            components.append(all_segments)

    # Assemble components into a single AGP file
    output_agp(components, args.contig_fasta, args.output)

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

        old_stdout = sys.stdout
        output_path = args.output + ".patch.fasta"
        with open(output_path, "w") as out_f:
            sys.stdout = out_f
            input_args = argparse.Namespace(agp=args.output + ".agp", components=args.output + ".relevant_seq.fasta")
            ragtag_agp2fa.main(input_args)
        sys.stdout = old_stdout
    

if __name__ == "__main__":
    main()
