import argparse
import os
import copy
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Set, Tuple, Dict, Optional, Any
from combine_connections import read_agp_file, merge_segments, reverse_segments
from ragtag_utilities.utilities import get_ragtag_version, get_impuT2T_version
import ragtag_agp2fa 

os.environ["MKL_THREADING_LAYER"] = "GNU"
import pysam
import math
import numpy as np
from sklearn.mixture import GaussianMixture
import sys

import time

# ==================== Constants ====================

# Simple, frequently-used constants at module level
DEFAULT_TIMEOUT = 300
RETRY_TIMEOUT = 300
CONNECTION_SCORE_MULTIPLIER = 10

# Edges with distance above this and probability below threshold are excluded from
# graph building / path finding but still written to .edge.log from gaussian_info.
LONG_LOW_CONF_EDGE_DIST_BP = 2_000_000
LONG_LOW_CONF_EDGE_PROB_MAX = 0.5


def edge_forbidden_long_low_conf(values) -> bool:
    """
    True iff this edge must not participate in the graph: dist strictly greater than
    LONG_LOW_CONF_EDGE_DIST_BP and probability strictly less than LONG_LOW_CONF_EDGE_PROB_MAX.
    values: [AF, AS, sample, dist, prob, ...] as in EdgeAggregator.edges.
    """
    if values is None or len(values) < 5:
        return False
    try:
        edge_dist = int(values[3])
        prob = float(values[4])
    except (TypeError, ValueError):
        return False
    if math.isnan(prob):
        return False
    return edge_dist > LONG_LOW_CONF_EDGE_DIST_BP and prob < LONG_LOW_CONF_EDGE_PROB_MAX


class OldSimplifiedConfig:
    """Old-simplified GMM champion selection (matches edge_distribution_FP.ipynb)."""
    MAX_COMPONENTS = 10
    WEIGHT_CUTOFF = 0.1
    MIN_SAMPLES = 10
    GMM_N_INIT = 5
    GMM_RANDOM_STATE = 0
    CHAMPION_BAND_STD_MULTIPLIER = 1.5
    CHAMPION_BAND_EXTRA_BP = 10.0
    LEGIT_SIGMA_MAX = 200000.0
    LEGIT_SCORE_MIN = 3000.0
    FALLBACK_TOP_N = 10
    FALLBACK_SCORE_FRAC = 0.5
    ABSOLUTE_SCORE_MIN = 6000
    NEAR_TOP_DELTA = 3
    TIE_SIGMA_MULTIPLIER = 100.0
    TIE_BIG_COUNT = 50

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


def info_prefix_to_tiebreak_dir(input_path: str) -> str:
    """
    input_path is the aggregate -p prefix, e.g. .../CN1.asm1. ending with a dot.
    Tiebreak artifacts go under a tiebreak/ folder next to the .info files.
    """
    base = input_path[:-1] if input_path.endswith(".") else input_path
    info_dir = os.path.dirname(base) or "."
    return os.path.join(info_dir, "tiebreak")


def read_info_edges(path: str) -> Dict[Tuple[str, str], Tuple[int, int]]:
    """Parse one *.info file: sorted (u,v) -> (score, dist)."""
    out: Dict[Tuple[str, str], Tuple[int, int]] = {}
    if not os.path.isfile(path):
        return out
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            u, v, score_s, dist_s = parts[0], parts[1], parts[2], parts[3]
            pair = tuple(sorted([u, v]))
            try:
                out[pair] = (int(score_s), int(dist_s))
            except ValueError:
                continue
    return out


def parse_tiebreak_manifest(
    manifest_path: str,
    default_ref: Optional[str] = None,
) -> Dict[str, Tuple[str, str, str]]:
    """
    TSV columns: sample_id, full_length.paf, query.fa, and optional reference.fa (omit 4th if --tiebreak_ref).
    If only 3 columns, default_ref must be set.
    """
    mapping: Dict[str, Tuple[str, str, str]] = {}
    with open(manifest_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 3:
                if default_ref is None:
                    raise ValueError(
                        f"Tiebreak manifest line has 3 columns but --tiebreak_ref is not set: {line!r}"
                    )
                sample, fl, q = parts
                mapping[sample.strip()] = (fl.strip(), q.strip(), default_ref.strip())
            elif len(parts) >= 4:
                sample, fl, q, ref = parts[0], parts[1], parts[2], parts[3]
                mapping[sample.strip()] = (fl.strip(), q.strip(), ref.strip())
            else:
                raise ValueError(f"Invalid tiebreak manifest line (need 3 or 4 columns): {line!r}")
    return mapping


@dataclass
class TiebreakConfig:
    tiebreak_dir: str
    map_length: int
    threads: int
    jobs: int
    from_paf_script: str
    sample_inputs: Dict[str, Tuple[str, str, str]]
    reuse: bool


def _run_tiebreak_sample_worker(args: Tuple[TiebreakConfig, str]) -> Tuple[str, Optional[str]]:
    """Top-level worker for ProcessPoolExecutor (picklable). Returns (sample, error_msg)."""
    cfg, sample = args
    try:
        run_tiebreak_from_paf(cfg, sample)
        return sample, None
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        return sample, str(e)


def run_tiebreak_from_paf(
    cfg: TiebreakConfig,
    sample: str,
) -> None:
    fl, q, r = cfg.sample_inputs[sample]
    out_prefix = os.path.join(cfg.tiebreak_dir, sample)
    info_path = out_prefix + ".info"
    if cfg.reuse and os.path.isfile(info_path):
        print(f"Tiebreak: reuse existing {info_path}")
        return
    for label, pth in (("-fl", fl), ("-q", q), ("-r", r)):
        if not os.path.isfile(pth):
            raise FileNotFoundError(f"Tiebreak {label} path for sample {sample} not found: {pth}")
    os.makedirs(cfg.tiebreak_dir, exist_ok=True)
    cmd = [
        sys.executable,
        cfg.from_paf_script,
        "-fl",
        fl,
        "-q",
        q,
        "-r",
        r,
        "-o",
        out_prefix,
        "-t",
        str(cfg.threads),
        "--map_length",
        str(cfg.map_length),
        "--info_only",
    ]
    print(f"Tiebreak: running from_paf for sample {sample} -> {info_path}")
    subprocess.run(cmd, check=True)


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
    def __init__(
        self,
        input_path: str,
        samples: list[str],
        sample_num: int,
        max_hap_total: int,
        max_hap_x: int,
        max_hap_y: int,
        choose_max_edge: bool = False,
        score_norm_max: float = 10000.0,
        tiebreak_cfg: Optional[TiebreakConfig] = None,
        enable_optional_tie_skips: bool = True,
    ):
        self.input_path = input_path
        self.samples = samples
        self.sample_num = sample_num
        self.max_hap_total = int(max_hap_total)
        self.max_hap_x = int(max_hap_x)
        self.max_hap_y = int(max_hap_y)
        self.choose_max_edge = bool(choose_max_edge)
        self.score_norm_max = float(score_norm_max)
        self.tiebreak_cfg = tiebreak_cfg
        self.enable_optional_tie_skips = bool(enable_optional_tie_skips)
        self.raw_edges = {}
        self.edges = {} # edge_occurence, edge_score, edge_sample, edge_distance, edge_probability
        self.gaussian_info = {}
        self.pending_tiebreak = {}
        # Aggregate-time tiebreak skips (sigma/distance/tie-count): log only if pair is used in final traversal.
        self.pending_skip_tiebreak_log: Dict[Tuple[str, str], Dict[str, Any]] = {}
        self._tiebreak_lookup = {}
        self._tiebreak_attempted_samples = set()
        self.tiebreak_resolution_lines = []

    def _max_hap_for_pair(self, pair: tuple[str, str]) -> int:
        a, b = str(pair[0]), str(pair[1])
        if ("_chrX_" in a) and ("_chrX_" in b):
            return int(self.max_hap_x)
        if ("_chrY_" in a) and ("_chrY_" in b):
            return int(self.max_hap_y)
        return int(self.max_hap_total)

    @classmethod
    def _fit_gmm_old_1d(cls, dists: np.ndarray):
        """Sklearn GMM + BIC, same as edge_distribution_FP old/ old_simplified."""
        cfg = OldSimplifiedConfig
        x = np.asarray(dists, dtype=float).reshape(-1, 1)
        bics, models = [], []
        for k in range(1, int(cfg.MAX_COMPONENTS) + 1):
            gmm = GaussianMixture(
                n_components=k,
                covariance_type="full",
                init_params="kmeans",
                n_init=int(cfg.GMM_N_INIT),
                random_state=int(cfg.GMM_RANDOM_STATE),
            ).fit(x)
            weights = gmm.weights_
            significant = weights > float(cfg.WEIGHT_CUTOFF)
            if significant.sum() < len(significant):
                break
            bics.append(gmm.bic(x))
            models.append(gmm)
        if not models:
            return None, None
        best_k_idx = int(np.argmin(bics))
        return models[best_k_idx], best_k_idx + 1

    @classmethod
    def _old_simplified_champions_per_component(cls, means, stds, dists, scores, samples):
        """One champion per GMM component; band 1.5σ+extra; tighten lo to μ−σ when lo<0 and μ>0."""
        cfg = OldSimplifiedConfig
        out = []
        mult = float(cfg.CHAMPION_BAND_STD_MULTIPLIER)
        extra = float(cfg.CHAMPION_BAND_EXTRA_BP)
        for mu, sig in zip(means, stds):
            half = mult * float(sig) + extra
            lo, hi = float(mu) - half, float(mu) + half
            if lo < 0 and float(mu) > 0:
                lo = float(mu) - float(sig)
            best_i = None
            for i in range(len(dists)):
                di, si = float(dists[i]), float(scores[i])
                if lo <= di <= hi:
                    if best_i is None or si > float(scores[best_i]):
                        best_i = i
            if best_i is not None:
                out.append(
                    (
                        float(mu),
                        float(sig),
                        float(dists[best_i]),
                        float(scores[best_i]),
                        str(samples[best_i]),
                    )
                )
        return out

    @classmethod
    def _select_winner_old_simplified(cls, scores, dists, samples):
        """
        Return (sel_dist, sel_score, sel_sample, sel_sigma, n_components).
        sel_sigma is the GMM component σ for legit_champion picks; -1.0 for greedy/union-style picks.
        n_components is len(GMM means) when a model was fit, else -1.
        """
        cfg = OldSimplifiedConfig
        s_list = [float(x) for x in scores]
        d_list = [float(x) for x in dists]
        n = len(d_list)
        idx_max = int(np.argmax(s_list))

        def greedy_return(nc: int = -1):
            return (
                float(d_list[idx_max]),
                float(s_list[idx_max]),
                str(samples[idx_max]),
                -1.0,
                int(nc),
            )

        if n < int(cfg.MIN_SAMPLES):
            return greedy_return()
        model, _bk = cls._fit_gmm_old_1d(np.asarray(d_list, dtype=float))
        if model is None:
            return greedy_return()
        means = model.means_.flatten()
        stds = np.sqrt(model.covariances_.flatten())
        n_comp = len(means)
        champs = cls._old_simplified_champions_per_component(means, stds, d_list, s_list, samples)
        if not champs:
            return greedy_return(n_comp)
        legit_max_sig = float(cfg.LEGIT_SIGMA_MAX)
        legit_min_score = float(cfg.LEGIT_SCORE_MIN)
        legit = [
            (d, s, sig, sm)
            for _mu, sig, d, s, sm in champs
            if float(sig) < legit_max_sig or float(s) >= legit_min_score
        ]
        if legit:
            legit.sort(key=lambda t: (-t[1], abs(t[0])))
            d, s, sig_w, sm = legit[0]
            return float(d), float(s), str(sm), float(sig_w), int(n_comp)
        max_s = max(s_list)
        thr = float(cfg.FALLBACK_SCORE_FRAC) * max_s
        order = sorted(range(n), key=lambda i: s_list[i], reverse=True)
        topn = int(min(int(cfg.FALLBACK_TOP_N), n))
        cand = set(order[:topn]) | {i for i in range(n) if s_list[i] > thr}
        if not cand:
            return greedy_return(n_comp)
        best_i = min(cand, key=lambda i: (abs(d_list[i]), -s_list[i]))
        return (
            float(d_list[best_i]),
            float(s_list[best_i]),
            str(samples[best_i]),
            -1.0,
            int(n_comp),
        )

    @classmethod
    def _gaussian_champion_tie_indices(
        cls,
        scores: List,
        dists: List,
        samples: List[str],
    ) -> List[int]:
        """
        Sample indices that share the same ambiguity as the GMM "champion" path:
        - Legit path: same integer flank score as the winning champion, distance within
          that winner's GMM-component band (1.5σ+10bp, with the same lo tightening as champions).
        - Fallback path: same minimization key (abs(dist), -score) among fallback candidates.
        - Greedy / no-GMM: all indices tied on the maximum flank score.
        """
        cfg = OldSimplifiedConfig
        s_list = [float(x) for x in scores]
        d_list = [float(x) for x in dists]
        n = len(d_list)
        if n == 0:
            return []
        idx_max = int(np.argmax(s_list))
        max_s = max(s_list)

        def _max_score_tie_indices() -> List[int]:
            return [i for i in range(n) if s_list[i] == max_s]

        if n < int(cfg.MIN_SAMPLES):
            return _max_score_tie_indices()
        model, _bk = cls._fit_gmm_old_1d(np.asarray(d_list, dtype=float))
        if model is None:
            return _max_score_tie_indices()
        means = model.means_.flatten()
        stds = np.sqrt(model.covariances_.flatten())
        n_comp = len(means)
        champs = cls._old_simplified_champions_per_component(means, stds, d_list, s_list, samples)
        if not champs:
            return _max_score_tie_indices()
        legit_max_sig = float(cfg.LEGIT_SIGMA_MAX)
        legit_min_score = float(cfg.LEGIT_SCORE_MIN)
        legit = [
            (d, s, sig, sm)
            for _mu, sig, d, s, sm in champs
            if float(sig) < legit_max_sig or float(s) >= legit_min_score
        ]
        if legit:
            legit.sort(key=lambda t: (-t[1], abs(t[0])))
            win_d, win_s, _win_sig, win_sm = legit[0]
            win_si = int(round(float(win_s)))
            lo = hi = None
            for mu, sig, d, s, sm in champs:
                if str(sm) != str(win_sm):
                    continue
                if abs(float(d) - float(win_d)) > 1e-3 or abs(float(s) - float(win_s)) > 1e-3:
                    continue
                mult = float(cfg.CHAMPION_BAND_STD_MULTIPLIER)
                extra = float(cfg.CHAMPION_BAND_EXTRA_BP)
                half = mult * float(sig) + extra
                lo, hi = float(mu) - half, float(mu) + half
                if lo < 0 and float(mu) > 0:
                    lo = float(mu) - float(sig)
                break
            if lo is None or hi is None:
                return _max_score_tie_indices()
            tie_idx = [
                i
                for i in range(n)
                if lo <= float(d_list[i]) <= hi and int(round(float(s_list[i]))) == win_si
            ]
            return tie_idx if tie_idx else _max_score_tie_indices()

        thr = float(cfg.FALLBACK_SCORE_FRAC) * max_s
        order = sorted(range(n), key=lambda i: s_list[i], reverse=True)
        topn = int(min(int(cfg.FALLBACK_TOP_N), n))
        cand = set(order[:topn]) | {i for i in range(n) if s_list[i] > thr}
        if not cand:
            return _max_score_tie_indices()
        best_i = min(cand, key=lambda i: (abs(d_list[i]), -s_list[i]))
        best_key = (abs(d_list[best_i]), -s_list[best_i])
        return [i for i in cand if (abs(d_list[i]), -s_list[i]) == best_key]

    @classmethod
    def _pick_within_indices(
        cls,
        scores: List,
        dists: List,
        samples: List[str],
        idxs: List[int],
        choose_max_edge: bool,
    ) -> Tuple[float, float, str]:
        """Deterministic pick among idxs when tiebreak subprocess is unavailable."""
        if not idxs:
            raise ValueError("_pick_within_indices: empty idxs")
        if choose_max_edge:
            max_sc = max(int(scores[i]) for i in idxs)
            pool = [i for i in idxs if int(scores[i]) == max_sc]
            best_i = min(pool)
        else:
            pool_sorted = sorted(
                idxs,
                key=lambda i: (-float(scores[i]), abs(float(dists[i])), str(samples[i])),
            )
            best_i = pool_sorted[0]
        return float(dists[best_i]), float(scores[best_i]), str(samples[best_i])

    @classmethod
    def _band_bounds(cls, mu: float, sig: float) -> Tuple[float, float]:
        cfg = OldSimplifiedConfig
        half = float(cfg.CHAMPION_BAND_STD_MULTIPLIER) * float(sig) + float(cfg.CHAMPION_BAND_EXTRA_BP)
        lo, hi = float(mu) - half, float(mu) + half
        if lo < 0 and float(mu) > 0:
            lo = float(mu) - float(sig)
        return lo, hi

    @classmethod
    def _pick_highest_score_closest_to_mean_idx(
        cls,
        idxs: List[int],
        scores: List,
        dists: List,
        samples: List[str],
        mu: float,
    ) -> int:
        if not idxs:
            raise ValueError("_pick_highest_score_closest_to_mean_idx: empty idxs")
        return sorted(
            idxs,
            key=lambda i: (-float(scores[i]), abs(float(dists[i]) - float(mu)), str(samples[i])),
        )[0]

    @classmethod
    def _pick_closest_to_mean_idx(
        cls,
        idxs: List[int],
        dists: List,
        samples: List[str],
        mu: float,
    ) -> int:
        if not idxs:
            raise ValueError("_pick_closest_to_mean_idx: empty idxs")
        return sorted(
            idxs,
            key=lambda i: (
                abs(float(dists[i]) - float(mu)),
                abs(float(dists[i])),
                str(samples[i]),
            ),
        )[0]

    @classmethod
    def _pick_fallback_low_score_idx(
        cls,
        scores: List,
        dists: List,
        samples: List[str],
    ) -> int:
        cfg = OldSimplifiedConfig
        if not scores:
            raise ValueError("_pick_fallback_low_score_idx: empty scores")
        s_list = [int(round(float(x))) for x in scores]
        max_s = max(s_list)
        thr = float(cfg.FALLBACK_SCORE_FRAC) * float(max_s)
        order = sorted(range(len(s_list)), key=lambda i: s_list[i], reverse=True)
        topn = int(min(int(cfg.FALLBACK_TOP_N), len(s_list)))
        cand = set(order[:topn]) | {i for i in range(len(s_list)) if float(s_list[i]) > thr}
        if not cand:
            cand = set(range(len(s_list)))
        return sorted(
            cand,
            key=lambda i: (
                abs(float(dists[i])),
                -int(s_list[i]),
                str(samples[i]),
            ),
        )[0]

    @classmethod
    def _summarize_tie_distribution(
        cls,
        dists: List,
        tie_idx: List[int],
    ) -> Tuple[float, float, bool, bool]:
        """
        Return (mu, sigma, all_negative, within_100_sigma) for tie candidates.
        Sigma is the standard deviation of tie distances.
        """
        if not tie_idx:
            return 0.0, 0.0, False, True
        vals = np.asarray([float(dists[i]) for i in tie_idx], dtype=float)
        mu = float(np.mean(vals))
        sigma = float(np.std(vals))
        all_negative = bool(np.all(vals < 0))
        if sigma <= 0:
            within_100_sigma = True
        else:
            within_100_sigma = bool(np.all(np.abs(vals - mu) <= float(OldSimplifiedConfig.TIE_SIGMA_MULTIPLIER) * sigma))
        return mu, sigma, all_negative, within_100_sigma

    def _tiebreak_info_path(self, sample: str) -> str:
        return os.path.join(self.tiebreak_cfg.tiebreak_dir, sample + ".info")

    def _load_tiebreak_lookup_for_sample(self, sample: str) -> None:
        if sample in self._tiebreak_attempted_samples:
            return
        self._tiebreak_attempted_samples.add(sample)
        tb_path = self._tiebreak_info_path(sample)
        if not os.path.isfile(tb_path):
            print(f"Tiebreak: no info file for sample {sample} ({tb_path})")
            return
        for p, sd in read_info_edges(tb_path).items():
            self._tiebreak_lookup[(p, sample)] = sd

    def _collect_needed_tiebreak_samples(
        self, used_pairs_sorted: Set[Tuple[str, str]]
    ) -> Set[str]:
        needed: Set[str] = set()
        for pair in used_pairs_sorted:
            if pair not in self.pending_tiebreak:
                continue
            for sm in self.pending_tiebreak[pair]["unique_tie_samples"]:
                if sm in self.tiebreak_cfg.sample_inputs:
                    needed.add(str(sm))
        return needed

    def _samples_missing_tiebreak_info(self, samples: Set[str]) -> List[str]:
        missing: List[str] = []
        for sample in sorted(samples):
            if self.tiebreak_cfg.reuse and os.path.isfile(self._tiebreak_info_path(sample)):
                continue
            missing.append(sample)
        return missing

    def _prefetch_tiebreak_parallel(self, samples: List[str]) -> None:
        if not samples:
            return
        jobs = max(1, int(self.tiebreak_cfg.jobs))
        print(f"Tiebreak: prefetch {len(samples)} sample(s) with {jobs} worker(s)")
        if jobs == 1 or len(samples) == 1:
            for sample in samples:
                try:
                    run_tiebreak_from_paf(self.tiebreak_cfg, sample)
                except (FileNotFoundError, subprocess.CalledProcessError) as e:
                    print(f"Tiebreak: skipping sample {sample} ({e})")
            return
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(_run_tiebreak_sample_worker, (self.tiebreak_cfg, sample)): sample
                for sample in samples
            }
            for future in as_completed(futures):
                sample = futures[future]
                try:
                    _, err = future.result()
                except Exception as e:
                    err = str(e)
                if err:
                    print(f"Tiebreak: skipping sample {sample} ({err})")

    def _ensure_tiebreak_sample_loaded(self, sample: str) -> None:
        """Load tiebreak .info for one sample (subprocess runs happen in prefetch)."""
        self._load_tiebreak_lookup_for_sample(sample)

    def _resolve_tie_for_pair(self, pair: Tuple[str, str], pending_info: dict) -> tuple[int, int, str, int]:
        """
        Resolve one deferred tie for a used edge.
        Returns (selected_score, selected_dist, selected_sample, tiebreak_used_flag).
        """
        scores = pending_info["scores"]
        dists = pending_info["dists"]
        relevant_samples = pending_info["relevant_samples"]
        tie_idx = pending_info["tie_idx"]
        unique_tie_samples = pending_info["unique_tie_samples"]
        tie_mu = float(pending_info.get("tie_mu", 0.0))
        original_selected_score = int(pending_info.get("original_selected_score", 0))
        force_no_subprocess = bool(pending_info.get("force_no_subprocess", False))
        force_reason = str(pending_info.get("force_reason", "fallback_within_ties"))

        if not force_no_subprocess:
            for sm in unique_tie_samples:
                if sm in self.tiebreak_cfg.sample_inputs:
                    self._ensure_tiebreak_sample_loaded(sm)

        tb_score_by_idx: Dict[int, int] = {}
        for i in tie_idx:
            sm = str(relevant_samples[i])
            tb_row = self._tiebreak_lookup.get((pair, sm))
            if tb_row is not None and not force_no_subprocess:
                tb_score_by_idx[i] = int(tb_row[0])
            else:
                tb_score_by_idx[i] = int(round(float(scores[i])))

        def pick_idx_with_mu(idxs: List[int], mu: float) -> int:
            return sorted(
                idxs,
                key=lambda i: (
                    -int(tb_score_by_idx.get(i, int(round(float(scores[i]))))),
                    -float(scores[i]),
                    abs(float(dists[i]) - float(mu)),
                    abs(float(dists[i])),
                    str(relevant_samples[i]),
                ),
            )[0]

        win_i = pick_idx_with_mu(tie_idx, tie_mu)

        selected_score = int(original_selected_score if original_selected_score > 0 else int(round(float(scores[win_i]))))
        selected_dist = int(round(float(dists[win_i])))
        selected_sample = str(relevant_samples[win_i])
        tiebreak_used = 0 if force_no_subprocess else int(any(self._tiebreak_lookup.get((pair, str(relevant_samples[i]))) is not None for i in tie_idx))
        resolution = force_reason if force_no_subprocess else ("tiebreak_subprocess" if tiebreak_used else "fallback_within_ties")

        self.tiebreak_resolution_lines.append(
            "\t".join(
                [
                    pair[0],
                    pair[1],
                    resolution,
                    str(len(tie_idx)),
                    ",".join(unique_tie_samples),
                    selected_sample,
                    str(selected_score),
                    str(selected_dist),
                ]
            )
        )
        return selected_score, selected_dist, selected_sample, tiebreak_used

    def resolve_tiebreaks_for_used_pairs(self, used_pairs_sorted: Set[Tuple[str, str]], dict_edge_info: dict) -> None:
        """
        Resolve deferred ties only for edges that are actually used in final paths.
        Append tiebreak_resolutions rows for deferred subprocess ties and for aggregate-time
        skipped tiebreaks (only when that edge is used).
        Updates self.edges/self.gaussian_info in place and syncs dict_edge_info.
        """
        if self.tiebreak_cfg is None:
            return

        needed_samples = self._collect_needed_tiebreak_samples(used_pairs_sorted)
        missing_samples = self._samples_missing_tiebreak_info(needed_samples)
        self._prefetch_tiebreak_parallel(missing_samples)
        for sample in sorted(needed_samples):
            self._load_tiebreak_lookup_for_sample(sample)

        for pair in sorted(used_pairs_sorted):
            if pair in self.pending_skip_tiebreak_log:
                log = self.pending_skip_tiebreak_log[pair]
                self.tiebreak_resolution_lines.append(
                    "\t".join(
                        [
                            pair[0],
                            pair[1],
                            "skipped:" + "|".join(log["skip_reasons"]),
                            str(log["n_tie"]),
                            ",".join(log["unique_tie_samples"]),
                            log["winner_sample"],
                            str(log["winner_score"]),
                            str(log["winner_dist"]),
                        ]
                    )
                )
            if pair not in self.pending_tiebreak:
                continue
            pending_info = self.pending_tiebreak[pair]
            selected_score, selected_dist, selected_sample, tiebreak_used = self._resolve_tie_for_pair(pair, pending_info)

            # Update edge row in place so any shared references remain valid.
            self.edges[pair][1] = int(selected_score)
            self.edges[pair][2] = str(selected_sample)
            self.edges[pair][3] = int(selected_dist)

            max_hap = self._max_hap_for_pair(pair)
            probability = predict_edge_groundtruth_prob(
                num_haplotypes=float(self.edges[pair][0]),
                chosen_score=float(self.edges[pair][1]),
                max_hap=float(max_hap),
                max_score=self.score_norm_max,
            )
            self.edges[pair][4] = probability

            # gaussian_info = [std, dist, selected_score, max_score, selected_sample, n_haps, n_comp, prob, tiebreak_used]
            self.gaussian_info[pair][1] = int(selected_dist)
            self.gaussian_info[pair][2] = int(selected_score)
            self.gaussian_info[pair][4] = str(selected_sample)
            self.gaussian_info[pair][7] = probability
            self.gaussian_info[pair][8] = int(tiebreak_used)

            # Keep dict_edge_info aligned (for AGP sample lookup).
            for k in list(dict_edge_info.keys()):
                if tuple(sorted([k[0], k[1]])) == pair:
                    dict_edge_info[k] = self.edges[pair]

        if self.tiebreak_resolution_lines:
            res_path = os.path.join(self.tiebreak_cfg.tiebreak_dir, "tiebreak_resolutions.tsv")
            os.makedirs(self.tiebreak_cfg.tiebreak_dir, exist_ok=True)
            with open(res_path, "w") as rf:
                rf.write(
                    "\t".join(
                        [
                            "u",
                            "v",
                            "resolution",
                            "n_gaussian_tie",
                            "tied_samples",
                            "winner_sample",
                            "winner_score",
                            "winner_dist",
                        ]
                    )
                    + "\n"
                )
                rf.write("\n".join(self.tiebreak_resolution_lines) + "\n")
            print(f"Tiebreak: wrote {res_path}")

    def aggregate(self):
        """
        edges:
        - key: (u, v) # sorted
        - value: (weight, max_score, max_sample_id)
        """
        for sample in self.samples:
            if os.path.exists(self.input_path + sample + ".info"):
                # .info lists each undirected edge twice (u,v) and (v,u); one row per pair per sample.
                seen_pair_this_sample: Set[Tuple[str, str]] = set()
                with open(self.input_path + sample + ".info", "r") as f:
                    for line in f.readlines():
                        u, v, score, dist = line.strip().split("\t")
                        pair = tuple(sorted([u, v]))
                        if pair in seen_pair_this_sample:
                            continue
                        seen_pair_this_sample.add(pair)
                        if pair not in self.raw_edges:
                            self.raw_edges[pair] = [[sample], [(int(score), int(dist))]]
                        else:
                            self.raw_edges[pair][0].append(sample)
                            self.raw_edges[pair][1].append((int(score), int(dist)))
            else:
                print(f"Warning: {sample} not found in {self.input_path}")

        for pair, info in self.raw_edges.items():
            relevant_samples = info[0]
            scores = [s for s, _d in info[1]]
            dists = [d for _s, d in info[1]]
            # Ensure aligned (score,dist,sample) lists.
            if len(relevant_samples) != len(scores) or len(scores) != len(dists):
                m = min(len(relevant_samples), len(scores), len(dists))
                relevant_samples = relevant_samples[:m]
                scores = scores[:m]
                dists = dists[:m]
            if len(scores) == 0:
                continue

            # Global max-score baseline (for logging / probability features).
            idx_max = int(np.argmax(scores))
            max_score = int(scores[idx_max])
            max_dist = int(dists[idx_max])
            max_sample = str(relevant_samples[idx_max])

            tiebreak_used = 0
            selected_score: int
            selected_dist: int
            selected_sample: str
            selected_std: float
            n_components: int
            cfg = OldSimplifiedConfig
            s_list = [int(round(float(s))) for s in scores]
            d_list = [int(round(float(d))) for d in dists]
            sorted_idx = sorted(range(len(s_list)), key=lambda i: s_list[i], reverse=True)
            s1 = int(s_list[sorted_idx[0]])
            s2 = int(s_list[sorted_idx[1]]) if len(sorted_idx) > 1 else -10**9
            top_idx = [i for i in range(len(s_list)) if s_list[i] == s1]
            unique_top = len(top_idx) == 1

            # 1) Absolute winner
            if s1 > int(cfg.ABSOLUTE_SCORE_MIN) and unique_top and s2 < (s1 - int(cfg.NEAR_TOP_DELTA)):
                win_i = top_idx[0]
                selected_score = int(s_list[win_i])
                selected_dist = int(d_list[win_i])
                selected_sample = str(relevant_samples[win_i])
                selected_std = -1.0
                n_components = -1
            # 2) Low-score fallback (legacy): no sample > 6000
            elif s1 <= int(cfg.ABSOLUTE_SCORE_MIN):
                win_i = self._pick_fallback_low_score_idx(s_list, d_list, relevant_samples)
                selected_score = int(s_list[win_i])
                selected_dist = int(d_list[win_i])
                selected_sample = str(relevant_samples[win_i])
                selected_std = -1.0
                n_components = -1
            else:
                print(f"High-score ambiguous set: {pair}!!!!!!!!!!!!!!!!!!")
                # High-score ambiguous set:
                # - Multiple samples tied at top s1: only those ties (ignore lower scores).
                # - Unique top: only s1 row(s) plus all samples at the second-highest distinct score (not a ±3 band).
                runner_idx: List[int] = []
                if not unique_top:
                    tie_idx = list(top_idx)
                else:
                    below_scores = [s_list[i] for i in range(len(s_list)) if s_list[i] < s1]
                    if not below_scores:
                        tie_idx = list(top_idx)
                    else:
                        s_second = max(below_scores)
                        runner_idx = [i for i in range(len(s_list)) if s_list[i] == s_second]
                        tie_idx = list(top_idx) + runner_idx
                        print(tie_idx)
                        print(runner_idx)
                        for i in tie_idx:
                            print(f"tie_idx: {i}, s_list: {s_list[i]}, d_list: {d_list[i]}, relevant_samples: {relevant_samples[i]}")


                unique_tie_samples = sorted({str(relevant_samples[i]) for i in tie_idx})
                tie_mu, tie_sigma, all_negative, within_100_sigma = self._summarize_tie_distribution(d_list, tie_idx)
                resolve_mu = float(tie_mu)

                # 3) 100σ skip case
                if self.enable_optional_tie_skips and within_100_sigma:
                    if unique_top:
                        win_i = top_idx[0]
                    else:
                        win_i = self._pick_closest_to_mean_idx(top_idx, d_list, relevant_samples, tie_mu)
                    selected_score = int(s_list[win_i])
                    selected_dist = int(d_list[win_i])
                    selected_sample = str(relevant_samples[win_i])
                    selected_std = float(tie_sigma)
                    n_components = -1
                    self.pending_skip_tiebreak_log[pair] = {
                        "skip_reasons": ["within_100sigma"],
                        "n_tie": len(tie_idx),
                        "unique_tie_samples": unique_tie_samples,
                        "winner_sample": selected_sample,
                        "winner_score": selected_score,
                        "winner_dist": selected_dist,
                    }
                # 4) all-negative skip case
                elif all_negative:
                    if unique_top:
                        win_i = top_idx[0]
                    else:
                        win_i = self._pick_closest_to_mean_idx(top_idx, d_list, relevant_samples, tie_mu)
                    selected_score = int(s_list[win_i])
                    selected_dist = int(d_list[win_i])
                    selected_sample = str(relevant_samples[win_i])
                    selected_std = float(tie_sigma)
                    n_components = -1
                    self.pending_skip_tiebreak_log[pair] = {
                        "skip_reasons": ["all_negative_ties"],
                        "n_tie": len(tie_idx),
                        "unique_tie_samples": unique_tie_samples,
                        "winner_sample": selected_sample,
                        "winner_score": selected_score,
                        "winner_dist": selected_dist,
                    }
                else:
                    # 5) >50 tie candidates with positive distances present.
                    if self.enable_optional_tie_skips and len(unique_tie_samples) > int(cfg.TIE_BIG_COUNT) and any(float(d_list[i]) > 0.0 for i in tie_idx):
                        if all(int(s_list[i]) == s1 for i in tie_idx):
                            # Case 1: only top-tier (multiple tied at s1) -> closest to mean among top.
                            win_i = self._pick_closest_to_mean_idx(top_idx, d_list, relevant_samples, tie_mu)
                            selected_score = int(s_list[win_i])
                            selected_dist = int(d_list[win_i])
                            selected_sample = str(relevant_samples[win_i])
                            selected_std = float(tie_sigma)
                            n_components = -1
                            self.pending_skip_tiebreak_log[pair] = {
                                "skip_reasons": ["tied_samples_gt_50_all_max"],
                                "n_tie": len(tie_idx),
                                "unique_tie_samples": unique_tie_samples,
                                "winner_sample": selected_sample,
                                "winner_score": selected_score,
                                "winner_dist": selected_dist,
                            }
                        elif unique_top and len(runner_idx) > int(cfg.TIE_BIG_COUNT):
                            # Case 2: unique top + many at second score: μ from runner tier only; one s2 champion vs top.
                            runner_mu = float(np.mean([float(d_list[i]) for i in runner_idx]))
                            runner_champion = self._pick_closest_to_mean_idx(
                                runner_idx, d_list, relevant_samples, runner_mu
                            )
                            top_champion = top_idx[0]
                            tie_idx = sorted({top_champion, runner_champion})
                            unique_tie_samples = sorted({str(relevant_samples[i]) for i in tie_idx})
                            resolve_mu = runner_mu
                            selected_score = int(s1)
                            selected_dist = int(d_list[top_champion])
                            selected_sample = f"__TIEBREAK__:{'|'.join(unique_tie_samples)}"
                            selected_std = float(tie_sigma)
                            n_components = -1
                            if self.tiebreak_cfg is not None and len(unique_tie_samples) > 1:
                                self.pending_tiebreak[pair] = {
                                    "scores": list(s_list),
                                    "dists": list(d_list),
                                    "relevant_samples": list(relevant_samples),
                                    "tie_idx": list(tie_idx),
                                    "unique_tie_samples": list(unique_tie_samples),
                                    "tie_mu": float(resolve_mu),
                                    "original_selected_score": int(s1),
                                    "force_no_subprocess": False,
                                    "force_reason": "",
                                }
                            else:
                                selected_sample = str(relevant_samples[top_champion])
                        else:
                            # Many tied samples but no runner-tier-only collapse: full tiebreak on full tie_idx.
                            top_champion = self._pick_closest_to_mean_idx(top_idx, d_list, relevant_samples, tie_mu)
                            selected_score = int(s1)
                            selected_dist = int(d_list[top_champion])
                            selected_sample = f"__TIEBREAK__:{'|'.join(unique_tie_samples)}"
                            selected_std = float(tie_sigma)
                            n_components = -1
                            if self.tiebreak_cfg is not None and len(unique_tie_samples) > 1:
                                self.pending_tiebreak[pair] = {
                                    "scores": list(s_list),
                                    "dists": list(d_list),
                                    "relevant_samples": list(relevant_samples),
                                    "tie_idx": list(tie_idx),
                                    "unique_tie_samples": list(unique_tie_samples),
                                    "tie_mu": float(resolve_mu),
                                    "original_selected_score": int(s1),
                                    "force_no_subprocess": False,
                                    "force_reason": "",
                                }
                            else:
                                selected_sample = str(relevant_samples[top_champion])
                    # 6) Full tiebreak on current ambiguous set
                    else:
                        top_champion = self._pick_closest_to_mean_idx(top_idx, d_list, relevant_samples, tie_mu)
                        selected_score = int(s1)
                        selected_dist = int(d_list[top_champion])
                        selected_sample = f"__TIEBREAK__:{'|'.join(unique_tie_samples)}"
                        selected_std = float(tie_sigma)
                        n_components = -1
                        if self.tiebreak_cfg is not None and len(unique_tie_samples) > 1:
                            self.pending_tiebreak[pair] = {
                                "scores": list(s_list),
                                "dists": list(d_list),
                                "relevant_samples": list(relevant_samples),
                                "tie_idx": list(tie_idx),
                                "unique_tie_samples": list(unique_tie_samples),
                                "tie_mu": float(resolve_mu),
                                "original_selected_score": int(s1),
                                "force_no_subprocess": False,
                                "force_reason": "",
                            }
                        else:
                            selected_sample = str(relevant_samples[top_champion])

            num_haplotypes = int(len(relevant_samples))
            max_hap = self._max_hap_for_pair(pair)
            probability = predict_edge_groundtruth_prob(
                num_haplotypes=float(num_haplotypes),
                chosen_score=float(selected_score),
                max_hap=float(max_hap),
                max_score=self.score_norm_max,
            )

            self.gaussian_info[pair] = [
                int(selected_std) if selected_std != -1 else -1,
                int(selected_dist),
                int(selected_score),
                int(max_score),
                str(selected_sample),
                int(len(relevant_samples)),
                int(n_components) if n_components is not None else -1,
                probability,
                int(tiebreak_used),
            ]
            self.edges[pair] = [
                int(len(relevant_samples)),
                int(selected_score),
                str(selected_sample),
                int(selected_dist),
                probability,
            ]

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
        if filter_type == "Conf":
            return lambda x: x[1][4]
        elif filter_type == "AF":
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

    def _apply_aggressive_filtering(self):
        """
        Apply aggressive edge filtering based on weight ratio.
        """
        for node, scores in self.dict_node_score.items():
            #print("!"*10, node, "!"*10, scores)
            for score, neighbor in scores:
                if score < 0.6:
                    self.filtered_edges.add((node, neighbor))
                    self.filtered_edges.add((neighbor, node))

    def _apply_white_list_edge_filtering(self):
        """
        Apply white list edge filtering based on white list file.
        """
        set_white_list = set()
        for node, scores in self.dict_node_score.items():
            max_score = max(scores, key=lambda x: x[0])[0]
            for score, neighbor in scores:
                if score >= max_score * 1:
                    set_white_list.add((node, neighbor))
                    set_white_list.add((neighbor, node))
        for node, scores in self.dict_node_score.items():
            for score, neighbor in scores:
                if (node, neighbor) not in set_white_list:
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
            if edge_forbidden_long_low_conf(values):
                continue
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

    def build(self, filter_num: int = GraphConfig.DEFAULT_FILTER_NUM, flag_filter: bool = False, filter_type: str = "Conf"):
        filter_func = self._get_filter_function(filter_type)
        edges_sorted = sorted(self.edges.items(), key=filter_func, reverse=True)
        # Node scores / edge info for pathfinding: omit long + low-confidence edges
        edges_for_graph = [
            (p, v) for p, v in edges_sorted if not edge_forbidden_long_low_conf(v)
        ]
        self._build_node_score_and_edge_info(edges_for_graph)

        # Apply local edge filtering
        if flag_filter:
            self._apply_local_edge_filtering()
            self._apply_aggressive_filtering()
        #self._apply_white_list_edge_filtering()

        # Process edges and add to graph
        add_counter = 0 # counter for filtering edges
        for pair, values in edges_sorted:
            u, v = pair
            AF, AS, _, edge_dist, prob = values
            #w = AF / (2 * self.sample_size) * (AS / 10000)
            w = prob

            if edge_forbidden_long_low_conf(values):
                self.filtered_edges.add((u, v))
                self.filtered_edges.add((v, u))
                continue

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
            if edge_forbidden_long_low_conf(values):
                continue
            u, v = pair
            AF, AS, *_ = values
            w = AF / (2 * self.sample_size) * (AS / 10000)
            
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



class PathFinderGreedy:
    def __init__(self, weights, base_components):
        self.weights = weights
        
        self.used_nodes = set()
        if base_components:
            self.connected_components = base_components
            for component in self.connected_components:
                for node in component:
                    self.used_nodes.add(node)
        else:
            self.connected_components= []
    
    def find_paths(self):
        for pair, confidence in self.weights.items():
            head_node, tail_node = pair
            if confidence < 0.05:
                continue
            if head_node == NodeUtils.get_opposite(tail_node):
                continue
            if head_node in self.used_nodes or tail_node in self.used_nodes:
                continue
        for ele in sorted(self.weights.items(), key=lambda x: x[1], reverse=True):
            pair, confidence = ele
            head_node, tail_node = ele[0]
            #print("head_node: ", head_node, "tail_node: ", tail_node)
            if confidence < 0.05:
                continue
            if head_node == NodeUtils.get_opposite(tail_node):
                continue
            if head_node in self.used_nodes or tail_node in self.used_nodes:
                continue
            overlap_components = []
            for component in self.connected_components:
                if NodeUtils.get_opposite(head_node) == component[0]:
                    overlap_components.append(('0', component))
                elif NodeUtils.get_opposite(head_node) == component[-1]:
                    overlap_components.append(('1', component))
                if NodeUtils.get_opposite(tail_node) == component[0]:
                    overlap_components.append(('2', component))
                elif NodeUtils.get_opposite(tail_node) == component[-1]:
                    overlap_components.append(('3', component))
            if len(overlap_components) == 0:
                print("\tA!")
                self.connected_components.append([head_node, tail_node])
                self.used_nodes.add(head_node)
                self.used_nodes.add(tail_node)
            elif len(overlap_components) == 1:
                print("\tB!")
                direction, component = overlap_components[0]
                if direction == '0':
                    self.connected_components.append([tail_node, head_node] + component)
                elif direction == '1':
                    self.connected_components.append(component + [head_node, tail_node])
                elif direction == '2':
                    self.connected_components.append([head_node, tail_node] + component)
                elif direction == '3':
                    self.connected_components.append(component + [tail_node, head_node])
                self.used_nodes.add(head_node)
                self.used_nodes.add(tail_node)
                self.connected_components.remove(component)
            elif len(overlap_components) == 2:
                print("\tC!")
                overlap_components = sorted(overlap_components, key=lambda x: x[0])
                direction_1, component_1 = overlap_components[0]
                direction_2, component_2 = overlap_components[1]
                if component_1 == component_2:
                    print("Skip loop connection!")
                    continue
                if direction_1 == '0' and direction_2 == '2':
                    new_component = component_1[::-1] + [head_node, tail_node] + component_2
                elif direction_1 == '0' and direction_2 == '3':
                    new_component = component_1[::-1] + [head_node, tail_node] + component_2[::-1]
                elif direction_1 == '1' and direction_2 == '2':
                    new_component = component_1 + [head_node, tail_node] + component_2
                elif direction_1 == '1' and direction_2 == '3':
                    new_component = component_1 + [head_node, tail_node] + component_2[::-1]
                # sort the component by the node position
                if new_component[-1].split("_")[-2] < new_component[0].split("_")[-2]:
                    new_component = new_component[::-1]
                self.connected_components.append(new_component)
                self.connected_components.remove(component_1)
                self.connected_components.remove(component_2)
                self.used_nodes.add(head_node)
                self.used_nodes.add(tail_node)
            else:
                print("Error: more than 2 overlap components")
            #print(head_node, tail_node)
            print(self.connected_components)
        return self.connected_components




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
        used_pairs = set()
        for idx in range(len(path)-1):
            u, v = path[idx][1], path[idx+1][0]
            used_pairs.add((u, v))
            used_pairs.add((v, u))
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
        return all_segments, used_pairs
    
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
    timeout_threshold: int = 3,
    ) -> PathFinder:
    """
    Run path search with timeout handling. No retry with filtered graph;
    on timeout keep best path so far and continue to next end_node. Break
    the end_node loop after timeout_threshold timeouts.
    Args:
        path_finder: PathFinder instance to use
        initial_timeout: Timeout in seconds per DFS attempt (default 300)
        retry_timeout: Unused (kept for API compatibility)
        timeout_threshold: Stop trying more end_nodes after this many timeouts (default 3)

    Returns:
        path_finder with best path(s) found so far (complete or partial).
    """
    end_nodes = path_finder.find_end_nodes_unvisited(visited_nodes=set())
    print("end_nodes: ", end_nodes)
    print("================================"*2)

    path_finder.search_start_time = time.time()
    path_finder.search_timeout_seconds = initial_timeout
    timeout_count = 0
    for end_node in end_nodes:
        try:
            print("end_node: ", end_node)
            print("--------------------------------"*2)
            path_finder.search_paths_DFS(
                current_path=[(end_node, NodeUtils.get_opposite(end_node))],
                current_node=end_node,
                connected_components=[],
                visited_nodes={end_node, NodeUtils.get_opposite(end_node)},
                visited_contigs={NodeUtils.get_opposite(end_node)},
                depth=0
            )
            path_finder.search_start_time = time.time()
        except TimeoutException as e:
            elapsed = time.time() - path_finder.search_start_time
            print(f"Timeout occurred after {elapsed:.2f} seconds: {e}")
            print("Keeping best path so far and continuing to next end_node.")
            timeout_count += 1
            if timeout_count >= timeout_threshold:
                print(f"Reached timeout threshold ({timeout_threshold}); stopping end_node loop.")
                break
            path_finder.search_start_time = time.time()
    return path_finder

def _capacity_fractions(
    num_haplotypes: float,
    chosen_score: float,
    max_hap: float,
    max_score: float,
) -> tuple[float, float]:
    """
    Same normalization as build_frac_xy in the training cell.

    max_hap: maximum possible num_haplotypes for this panel (POP-470 -> 470, POP-300 -> 300).
    num_haplotypes is clamped to [0, max_hap] before dividing.
    """
    nh = min(max(float(num_haplotypes), 0.0), float(max_hap))
    hap_frac = nh / float(max_hap)
    score_frac = float(chosen_score) / float(max_score)
    return hap_frac, score_frac


def predict_edge_groundtruth_prob(
    num_haplotypes: float,
    chosen_score: float,
    max_hap: float,
    max_score: float = 10000.0,
) -> float:
    """
    P(edge is in benchmark groundtruth) from HG002 POP-470 hap_frac/score_frac LR.

    max_hap: panel maximum num_haplotypes (same as edge.log cap for that cohort).
    hap_frac = min(num_haplotypes, max_hap) / max_hap; score_frac = chosen_score / max_score.
    """
    if max_hap <= 0 or max_score <= 0:
        raise ValueError("max_hap and max_score must be > 0.")
    if num_haplotypes is None or chosen_score is None:
        return float("nan")
    if num_haplotypes == -1 or chosen_score == -1:
        return float("nan")

    # Hard-coded from edge_analysis.0516 hap_frac/score_frac train cell (HG002-470)
    scaler_center_hap_frac = 0.03191489
    scaler_center_score_frac = 0.9367
    scaler_scale_hap_frac = 0.30638298
    scaler_scale_score_frac = 0.5004
    lr_coef_hap_frac = 2.064144
    lr_coef_score_frac = 1.653317
    lr_intercept = -2.430398471792

    hap_frac, score_frac = _capacity_fractions(num_haplotypes, chosen_score, max_hap, max_score)
    x0 = (hap_frac - scaler_center_hap_frac) / scaler_scale_hap_frac
    x1 = (score_frac - scaler_center_score_frac) / scaler_scale_score_frac
    logit = lr_intercept + lr_coef_hap_frac * x0 + lr_coef_score_frac * x1
    return 1.0 / (1.0 + math.exp(-logit))




def main():                        
    parser = argparse.ArgumentParser("impuT2T patching version: " + get_impuT2T_version())
    parser.add_argument("-p", "--path", required=True, help="Input path")
    parser.add_argument("-l", "--list_sample", required=True, help="file containing list of samples")
    parser.add_argument("-n", "--sample_num", help="sample number, used for chrY, chrX normalization")
    parser.add_argument("-w", "--weight", help="weight type: allele_freq or max_alignment_score [Conf|AF|AS]", default="Conf")
    parser.add_argument("-wr", "--weight_ratio", help="weight ratio, discard the edge less than ratio of the max weight the same node [0.5]", type=float, default=0.5)
    #parser.add_argument("-ast", "--AS_threshold", help="AS threshold, discard the edge less than threshold [6000]", type=int, default=6000)
    #parser.add_argument("-aft", "--AF_threshold", help="AF threshold, discard the edge less than threshold [0.005]", type=float, default=0.005)
    parser.add_argument("--timeout", help=f"timeout for full path traversal [{DEFAULT_TIMEOUT}(s)]", type=int, default=DEFAULT_TIMEOUT)
    parser.add_argument("--retry_timeout", help=f"timeout for retry attempt [{RETRY_TIMEOUT}(s)]", type=int, default=RETRY_TIMEOUT)
    parser.add_argument("--timeout_threshold", help="stop trying more end_nodes after this many timeouts [3]", type=int, default=3)
    parser.add_argument("--debug", help="debug mode", action="store_true")
    parser.add_argument("--contig_fasta", required=True, help="contig fasta file, needed for unplaced contigs")
    parser.add_argument("--ref_path", help="reference path, only specify if you want the patch fasta output")
    parser.add_argument("--second_pass", help="second pass connection", action="store_true")
    parser.add_argument(
        "--choose_max_edge",
        action="store_true",
        help="Per edge, choose the haplotype with maximum alignment score (greedy). "
        "Default is old_simplified GMM + champion bands (1.5σ+10bp, tighten lower to μ−σ when lo<0 and μ>0).",
    )
    parser.add_argument(
        "--max_hap_total",
        type=int,
        help="Total haplotypes available (default: number of samples in list_sample).",
    )
    parser.add_argument(
        "--max_hap_x",
        type=int,
        help="Haplotype count to use for chrX edges (default: round(3/4 * max_hap_total)).",
    )
    parser.add_argument(
        "--max_hap_y",
        type=int,
        help="Haplotype count to use for chrY edges (default: round(1/4 * max_hap_total)).",
    )
    parser.add_argument(
        "--score_norm_max",
        type=float,
        default=10000.0,
        help="Normalization max score used in edge probability model [10000].",
    )
    parser.add_argument(
        "--tiebreak",
        action="store_true",
        help="After the usual GMM champion logic, if several haplotypes tie (same flank score within the "
        "winner's GMM distance band, or same fallback key), rerun from_paf_to_multi_connections with "
        "--tiebreak_map_length when >=2 tied samples have tiebreak/*.info rows; otherwise pick within the "
        "tie set (argmax score if --choose_max_edge, else same ordering as GMM fallback).",
    )
    parser.add_argument(
        "--tiebreak_manifest",
        help="TSV per line: sample_id, full_length.paf, query.fa, optional reference.fa (omit 4th if --tiebreak_ref).",
    )
    parser.add_argument(
        "--tiebreak_ref",
        help="Reference FASTA used for all samples when manifest lines have only 3 columns.",
    )
    parser.add_argument(
        "--tiebreak_map_length",
        type=int,
        default=25000,
        help="Terminal flank length (--map_length) for tiebreak reruns [25000].",
    )
    parser.add_argument(
        "--tiebreak_threads",
        type=int,
        default=1,
        help="Threads per from_paf_to_multi_connections (-t) tiebreak subprocess [1]. "
        "Use with --tiebreak_jobs for parallel samples (jobs × threads ≈ total cores).",
    )
    parser.add_argument(
        "--tiebreak_jobs",
        type=int,
        default=1,
        help="Parallel tiebreak subprocesses after greedy traversal resolves which samples tied [1].",
    )
    parser.add_argument(
        "--tiebreak_dir",
        help="Directory for tiebreak outputs (default: <info_dir>/tiebreak next to -p prefix).",
    )
    parser.add_argument(
        "--tiebreak_from_paf",
        help="Path to from_paf_to_multi_connections.py (default: alongside this script).",
    )
    parser.add_argument(
        "--tiebreak_no_reuse",
        action="store_true",
        help="Always rerun tiebreak from_paf even if tiebreak/<sample>.info already exists.",
    )
    parser.add_argument(
        "--no_optional_tie_skips",
        action="store_true",
        help="Disable optional tie skip shortcuts (within_100sigma and >50 tied samples). "
        "All-negative tie skip remains enabled.",
    )
    parser.add_argument("-o", "--output", help="output path", default="patch.result")
    args = parser.parse_args()

    if args.tiebreak and not args.tiebreak_manifest:
        parser.error("--tiebreak requires --tiebreak_manifest")

    # Collect edges from all samples
    if args.weight == "Conf":
        weight_type = 4
    elif args.weight == "AF":
        weight_type = 0
    elif args.weight == "AS":
        weight_type = 1
    else:
        raise ValueError(f"Invalid weight type: {args.weight}")

    samples = parse_sample_list(args.list_sample)
    if args.sample_num:
        sample_num = int(args.sample_num)
    else:
        sample_num = len(samples)
    max_hap_total = int(args.max_hap_total) if args.max_hap_total is not None else int(len(samples))
    if args.max_hap_x is None:
        # max_hap_x = int(round(0.75 * max_hap_total))
        max_hap_x = max_hap_total
    else:
        max_hap_x = int(args.max_hap_x)
    if args.max_hap_y is None:
        #max_hap_y = int(round(0.25 * max_hap_total))
        max_hap_y = max_hap_total * 0.5
    else:
        max_hap_y = int(args.max_hap_y)
    # basic guardrails
    max_hap_total = max(1, max_hap_total)
    max_hap_x = max(1, max_hap_x)
    max_hap_y = max(1, max_hap_y)

    tiebreak_cfg: Optional[TiebreakConfig] = None
    if args.tiebreak:
        manifest = parse_tiebreak_manifest(args.tiebreak_manifest, args.tiebreak_ref)
        tb_dir = args.tiebreak_dir or info_prefix_to_tiebreak_dir(args.path)
        fp_script = args.tiebreak_from_paf or os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "from_paf_to_multi_connections.py",
        )
        tiebreak_cfg = TiebreakConfig(
            tiebreak_dir=tb_dir,
            map_length=int(args.tiebreak_map_length),
            threads=int(args.tiebreak_threads),
            jobs=max(1, int(args.tiebreak_jobs)),
            from_paf_script=fp_script,
            sample_inputs=manifest,
            reuse=not bool(args.tiebreak_no_reuse),
        )

    edge_aggregator = EdgeAggregator(
        args.path,
        samples,
        sample_num,
        max_hap_total,
        max_hap_x,
        max_hap_y,
        choose_max_edge=bool(args.choose_max_edge),
        score_norm_max=float(args.score_norm_max),
        tiebreak_cfg=tiebreak_cfg,
        enable_optional_tie_skips=not bool(args.no_optional_tie_skips),
    )
    edges = edge_aggregator.aggregate()
    if args.debug:
        for pair, values in sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2])):
            print(pair, values)

    graph_builder = GraphBuilder(edges, weight_type, args.weight_ratio, args.contig_fasta, len(samples))
    graph, weights, nodes, dict_edge_info, dict_node_outgoing_edges = graph_builder.build(flag_filter=True, filter_type=args.weight)
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
    
    path_finder = PathFinder(graph, weights, nodes, dict_node_outgoing_edges)
    graph_config = GraphBuildConfig(
        edges=edges,
        weight_type=weight_type,
        weight_ratio=args.weight_ratio,
        contig_fasta=args.contig_fasta,
        sample_size=len(samples),
        filter_num=GraphConfig.RETRY_FILTER_NUM_FIRST_PASS
    )
    path_finder = run_path_search_with_retry(path_finder, graph_config, args.timeout, args.retry_timeout, args.timeout_threshold)
    
    first_pass_results = path_finder.best_score_component['component']
    first_pass_results = path_finder.sort_component(first_pass_results)
    _, _, first_pass_scores = path_finder.component_score(first_pass_results)
    #first_pass_results_connected = path_finder.sort_component(first_pass_results_connected)
    print(first_pass_results)
    print("--------------------------------"*2)

    base_components = []
    for component in first_pass_results:
        new_component = []
        for pair in component:
            new_component.append(pair[0])
            new_component.append(pair[1])
        if len(new_component) > 2:
            base_components.append(new_component[1:-1])
        else:
            print("Warning: component is too short, skipping", new_component)
            continue
    print(base_components)

    graph_builder = GraphBuilder(edges, weight_type, args.weight_ratio, args.contig_fasta, len(samples))
    graph, weights, nodes, dict_edge_info, dict_node_outgoing_edges = graph_builder.build(flag_filter=False, filter_type=args.weight)
    edges_sorted = sorted(edges.items(), key=lambda x: int(x[0][0].split("_")[-2]))
 
    print("weights", "--------------------------------"*2)
    for ele in weights:
        print("weights", ele, weights[ele])
    print("--------------------------------"*2)

    path_finder_greedy = PathFinderGreedy(weights, base_components)
    connected_components = path_finder_greedy.find_paths()
    greedy_results = []
    for component in connected_components:
        greedy_component = [NodeUtils.get_opposite(component[0])] + component + [NodeUtils.get_opposite(component[-1])]
        greedy_component = [tuple(greedy_component[i:i+2]) for i in range(0, len(greedy_component), 2)]
        greedy_results.append(greedy_component)
    
    print("greedy_results", "--------------------------------"*2)
    for greedy_component in greedy_results:
        print(greedy_component)

    # Resolve deferred tiebreaks only for edges actually used by final greedy traversal.
    used_pairs_sorted = set()
    for path in greedy_results:
        for idx in range(len(path) - 1):
            u, v = path[idx][1], path[idx + 1][0]
            used_pairs_sorted.add(tuple(sorted([u, v])))
    edge_aggregator.resolve_tiebreaks_for_used_pairs(used_pairs_sorted, dict_edge_info)

    components = []
    total_used_pairs = set()
    agp_assembler = AGPAssembler(args.path, args.contig_fasta)
    for path in greedy_results:
        all_segments, used_pairs = agp_assembler.extract_and_merge(path, dict_edge_info, args.path)
        total_used_pairs.update(used_pairs)
        components.append(all_segments)

    # Assemble components into a single AGP file
    output_agp(components, args.contig_fasta, args.output)

    # Output edge info in a edge_log file
    with open(args.output + ".edge.log", "w") as f:
        for pair, info in sorted(edge_aggregator.gaussian_info.items(), key=lambda x: int(x[0][0].split("_")[-2])):
            probability = info[7]
            if pair in total_used_pairs:
                f.write("\t".join(pair) + "\t")
                f.write("\t".join([str(x) for x in info]) + "\t")
                f.write(f"{probability:.3f}\t".replace("nan", "0.000"))
                f.write("True\n")
            else:
                f.write("\t".join(pair) + "\t")
                f.write("\t".join([str(x) for x in info]) + "\t")
                f.write(f"{probability:.3f}\t".replace("nan", "0.000"))
                f.write("False\n")


    if args.ref_path:
        # Collect relevant samples based on the edges actually used in the final patch graph.
        # This is more robust than parsing component_id prefixes and works even when
        # reference FASTA sequence names do not encode the sample ID (e.g. YAO#1#...).
        relevant_samples = set()
        for pair in total_used_pairs:
            if pair in dict_edge_info:
                sample_id = dict_edge_info[pair][2]
                if sample_id in samples:
                    relevant_samples.add(sample_id)
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
