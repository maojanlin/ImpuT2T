"""Integration tests require pysam, networkx, sklearn, numpy and minimap2 on PATH.

IMPUT2T_BASELINE_SCRIPT optionally runs the unmodified scorer as an independent oracle.
"""
import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "snakemake/scripts"
sys.path.insert(0, str(SCRIPTS))
import from_paf_to_multi_connections as scorer
import targeted_tiebreak as targeted
import pysam

spec = importlib.util.spec_from_file_location("aggregate", SCRIPTS / "aggregate_paths_info_semi-greedy_greedy-tiebreak_0517.py")
aggregate = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = aggregate
spec.loader.exec_module(aggregate)


def fixture(root):
    rng = random.Random(93)
    sequence = ''.join(rng.choices('ACGT', k=260000))
    ref = root / 'ref.fa'
    ref.write_text('>chrFixture\n' + sequence + '\n')
    query = root / 'query.fa'
    paf = root / 'full.paf'
    fasta, rows = [], []
    # Both orientations, underscores in names, a short contig, and overlapping fragments.
    for i, (start, length, strand) in enumerate([(0, 35000, '+'), (36000, 34000, '-'), (71000, 8000, '+'), (80000, 40000, '+'), (115000, 40000, '-'), (170000, 45000, '+')]):
        name = f'fragment_{i}'
        seq = sequence[start:start + length]
        if strand == '-':
            seq = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        fasta.append(f'>{name}\n{seq}\n')
        rows.append(f'{name}\t{length}\t0\t{length}\t{strand}\tchrFixture\t{len(sequence)}\t{start}\t{start + length}\t{length}\t{length}\t60\ttp:A:P\n')
    query.write_text(''.join(fasta))
    paf.write_text(''.join(rows))
    pysam.faidx(str(query))
    return str(paf), str(query), str(ref)


class TargetedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.inputs = fixture(self.root)
        self.cfg = aggregate.TiebreakConfig(str(self.root / 'tb'), 25000, 1, 1, str(SCRIPTS / 'from_paf_to_multi_connections.py'), {'donorA': self.inputs, 'donorB': self.inputs}, True)

    def tearDown(self):
        self.temp.cleanup()

    def full_rows(self):
        graph, _, _ = scorer.build_connection_evidence(*self.inputs, str(self.root / 'full'), MAP_LENGTH=25000)
        return {tuple(sorted((u, v))): (graph[u][v]['score'], graph[u][v]['dist']) for u, v in graph.edges}

    def test_parity_and_only_requested_terminals(self):
        full = self.full_rows()
        self.assertGreaterEqual(len(full), 3)
        baseline = os.environ.get('IMPUT2T_BASELINE_SCRIPT')
        if baseline:
            out = str(self.root / 'oracle')
            subprocess.run([sys.executable, baseline, '-fl', self.inputs[0], '-q', self.inputs[1], '-r', self.inputs[2], '-o', out, '--map_length', '25000', '--info_only'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            self.assertEqual(full, aggregate.read_info_edges(out + '.info'))
        # Verify each requested pair independently, including short and reversed terminals.
        for n, pair in enumerate(sorted(full)):
            _, rows, windows = scorer.build_connection_evidence(*self.inputs, str(self.root / f'target{n}'), MAP_LENGTH=25000, requested_pairs={pair})
            self.assertEqual(rows, {pair: full[pair]})
            self.assertLessEqual(set(windows), set(pair))
            headers = [line.split()[0][1:] for line in (self.root / f'target{n}.ends.fa').read_text().splitlines() if line.startswith('>')]
            expected = {node[:-2] + ('_head' if node.endswith('_b') else '_tail') for node in windows}
            self.assertEqual(set(headers), expected)

    def test_parity_with_split_segments_and_duplicate_connections(self):
        paf = Path(self.inputs[0])
        lines = paf.read_text().splitlines()
        first = lines[0].split('\t')
        left, right = first.copy(), first.copy()
        left[3] = left[8] = left[9] = left[10] = '22000'
        right[2] = right[7] = '18000'
        right[9] = right[10] = '17000'
        # Overlapping primary segments of one fragment merge in the legacy geometry.
        lines = ['\t'.join(left), '\t'.join(right), *lines[1:]]
        ref = Path(self.inputs[2])
        seq = ref.read_text().splitlines()[1]
        ref.write_text(ref.read_text() + '>chrDuplicate\n' + seq + '\n')
        duplicates = []
        for line in lines:
            fields = line.split('\t'); fields[5] = 'chrDuplicate'
            duplicates.append('\t'.join(fields))
        paf.write_text('\n'.join(lines + duplicates) + '\n')
        self.test_parity_and_only_requested_terminals()

    def test_cache_batch_rollback(self):
        cache = targeted.EvidenceCache(self.root / 'rollback.sqlite3')
        try:
            good = ('fragment_0_e', 'fragment_1_e')
            bad = ('fragment_2_e', 'fragment_3_b')
            with self.assertRaises(ValueError):
                cache.put_batch('fingerprint', 'donorA', {good: {'status': 'scored', 'score': 1}, bad: {'status': 'error'}})
            self.assertIsNone(cache.get('fingerprint', 'donorA', good))
        finally:
            cache.close()

    def test_reuse_export_refresh_stale_and_missing(self):
        pair = ('fragment_0_e', 'fragment_1_e')
        absent = ('missing_b', 'missing_e')
        requests = {'donorA': {pair, absent}}
        cold = targeted.resolve_requests(self.cfg, requests)
        self.assertEqual(cold[(pair, 'donorA')]['status'], 'scored')
        self.assertEqual(cold[(absent, 'donorA')]['status'], 'no_connection')
        self.assertFalse((self.root / 'tb/targeted/info').exists())
        self.cfg.cache_mode = 'only'
        self.cfg.dump_info = True
        with patch.object(targeted, '_worker', side_effect=AssertionError('warm cache aligned')):
            warm = targeted.resolve_requests(self.cfg, requests)
        self.assertTrue(all(r['origin'] == 'cache' for r in warm.values()))
        self.assertEqual(aggregate.read_info_edges(str(self.root / 'tb/targeted/info/donorA.info')), {pair: (cold[(pair, 'donorA')]['score'], cold[(pair, 'donorA')]['distance'])})
        self.assertTrue(json.loads((self.root / 'tb/targeted/info.metadata.json').read_text())['partial'])
        self.cfg.map_length = 24000
        with self.assertRaisesRegex(RuntimeError, 'cache-only'):
            targeted.resolve_requests(self.cfg, requests)
        self.cfg.map_length = 25000
        # Same-size content change must invalidate provenance.
        paf = Path(self.inputs[0]); paf.write_text(paf.read_text().replace('\t60\t', '\t59\t'))
        with self.assertRaisesRegex(RuntimeError, 'cache-only'):
            targeted.resolve_requests(self.cfg, requests)
        self.cfg.cache_mode = 'refresh'
        refreshed = targeted.resolve_requests(self.cfg, requests)
        self.assertTrue(all(r['origin'] == 'computed' for r in refreshed.values()))
        self.cfg.cache_mode = 'only'
        with self.assertRaisesRegex(RuntimeError, 'cache-only'):
            targeted.resolve_requests(self.cfg, {'unknown': {pair}})

    def test_donor_batches_and_request_collection(self):
        a = ('fragment_0_e', 'fragment_1_e')
        b = ('fragment_2_e', 'fragment_3_b')
        c = ('unused_b', 'unused_e')
        pending = {a: {'unique_tie_samples': ['donorA']}, b: {'unique_tie_samples': ['donorB']}, c: {'unique_tie_samples': ['donorA', 'donorB']}}
        requests = targeted.collect_requests({a, b}, pending)
        self.assertEqual(requests, {'donorA': {a}, 'donorB': {b}})
        self.cfg.jobs = 2
        records = targeted.resolve_requests(self.cfg, requests)
        self.assertEqual(set(records), {(a, 'donorA'), (b, 'donorB')})
        pending[a]['force_no_subprocess'] = True
        self.assertEqual(targeted.collect_requests({a, b}, pending), {'donorB': {b}})
        with patch.object(targeted, '_worker', side_effect=AssertionError('empty aligned')):
            self.assertEqual(targeted.resolve_requests(self.cfg, {}), {})

    def test_failed_worker_not_cached(self):
        pair = ('fragment_0_e', 'fragment_1_e')
        with patch.object(targeted, '_worker', side_effect=subprocess.CalledProcessError(1, 'minimap2')):
            result = targeted.resolve_requests(self.cfg, {'donorA': {pair}})
        self.assertEqual(result[(pair, 'donorA')]['status'], 'error')
        self.cfg.cache_mode = 'only'
        with self.assertRaisesRegex(RuntimeError, 'cache-only'):
            targeted.resolve_requests(self.cfg, {'donorA': {pair}})

    def test_indexed_extraction_matches_python_slices(self):
        heads = {'fragment_2': (-17000, 25000)}
        tails = {'fragment_2': (-17000, 8000)}
        indexed = self.root / 'indexed.fa'
        scan = self.root / 'scan.fa'
        scorer.write_terminal_fasta(self.inputs[1], str(indexed), heads, tails)
        os.unlink(self.inputs[1] + '.fai')
        scorer.write_terminal_fasta(self.inputs[1], str(scan), heads, tails)
        self.assertEqual(indexed.read_bytes(), scan.read_bytes())

    def test_cli_agp_and_edge_decisions_match_legacy(self):
        prefix = self.root / 'initial.'
        for donor in ('donorA', 'donorB'):
            out = str(prefix) + donor
            subprocess.run([sys.executable, str(SCRIPTS / 'from_paf_to_multi_connections.py'), '-fl', self.inputs[0], '-q', self.inputs[1], '-r', self.inputs[2], '-o', out, '--info_only'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        donors = self.root / 'donors.txt'
        donors.write_text('donorA\ndonorB\n')
        manifest = self.root / 'manifest.tsv'
        manifest.write_text(''.join('\t'.join([donor, *self.inputs]) + '\n' for donor in ('donorA', 'donorB')))
        for backend in ('legacy', 'targeted'):
            cmd = [sys.executable, str(SCRIPTS / 'aggregate_paths_info_semi-greedy_greedy-tiebreak_0517.py'), '-p', str(prefix), '-l', str(donors), '-o', str(self.root / backend), '--contig_fasta', self.inputs[1], '--tiebreak', '--tiebreak_manifest', str(manifest), '--tiebreak_dir', str(self.root / (backend + '_tb')), '--no_optional_tie_skips', '--timeout', '5']
            baseline = os.environ.get('IMPUT2T_BASELINE_SCRIPT')
            if backend == 'legacy' and baseline:
                cmd[1] = str(Path(baseline).parent / 'aggregate_paths_info_semi-greedy_greedy-tiebreak_0517.py')
            elif backend == 'legacy':
                cmd.extend(['--tiebreak_backend', 'legacy'])
            with (self.root / (backend + '.log')).open('w') as log:
                result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, (self.root / (backend + '.log')).read_text()[-8000:])
        for suffix in ('.agp', '.edge.log'):
            self.assertEqual((self.root / ('legacy' + suffix)).read_bytes(), (self.root / ('targeted' + suffix)).read_bytes())
        self.assertEqual((self.root / 'legacy_tb/tiebreak_resolutions.tsv').read_bytes(), (self.root / 'targeted_tb/tiebreak_resolutions.tsv').read_bytes())
        self.assertTrue((self.root / 'targeted_tb/targeted/cache.sqlite3').is_file())
        self.assertFalse((self.root / 'targeted_tb/targeted/info').exists())

    def test_resolver_preserves_original_score_distance_and_skips(self):
        pair = ('fragment_0_e', 'fragment_1_e')
        agg = aggregate.EdgeAggregator('', ['donorA', 'donorB'], 2, 2, 2, 2, tiebreak_cfg=self.cfg)
        info = {'scores': [9000, 9100], 'dists': [100, 200], 'relevant_samples': ['donorA', 'donorB'], 'tie_idx': [0, 1], 'unique_tie_samples': ['donorA', 'donorB'], 'tie_mu': 150, 'original_selected_score': 9100}
        agg._tiebreak_lookup = {(pair, 'donorA'): (49000, 999), (pair, 'donorB'): (48000, 888)}
        self.assertEqual(agg._resolve_tie_for_pair(pair, info), (9100, 100, 'donorA', 1))
        info['force_no_subprocess'] = True
        self.assertEqual(agg._resolve_tie_for_pair(pair, info), (9100, 200, 'donorB', 0))


if __name__ == '__main__':
    unittest.main()
