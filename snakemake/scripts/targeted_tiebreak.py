"""Donor-batched terminal scoring and a provenance-checked, per-edge SQLite cache."""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

SCHEMA_VERSION = 1


def canonical_pair(pair):
    pair = tuple(sorted(pair))
    if len(pair) != 2 or any(not node.endswith(("_b", "_e")) for node in pair):
        raise ValueError(f"Invalid tiebreak endpoint pair: {pair!r}")
    return pair


def collect_requests(used_pairs, pending):
    """Preserve edge/donor membership, excluding ties explicitly skipped upstream."""
    requests = {}
    for pair in sorted(used_pairs):
        info = pending.get(pair)
        if info is None or info.get("force_no_subprocess", False):
            continue
        for sample in info["unique_tie_samples"]:
            requests.setdefault(str(sample), set()).add(canonical_pair(pair))
    return requests


def _digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


class EvidenceCache:
    """Only the parent writes. A completed donor batch is committed atomically."""
    def __init__(self, path):
        self.db = sqlite3.connect(path, timeout=60)
        self.db.execute("CREATE TABLE IF NOT EXISTS identities (path TEXT PRIMARY KEY, signature TEXT, digest TEXT)")
        self.db.execute("CREATE TABLE IF NOT EXISTS evidence (fingerprint TEXT, donor TEXT, u TEXT, v TEXT, record TEXT, PRIMARY KEY (fingerprint, donor, u, v))")
        self.memo = {}

    def identity(self, path):
        path = os.path.realpath(path)
        st = os.stat(path)
        signature = json.dumps([st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns])
        key = (path, signature)
        if key in self.memo:
            return self.memo[key]
        old = self.db.execute("SELECT signature, digest FROM identities WHERE path=?", (path,)).fetchone()
        digest = old[1] if old and old[0] == signature else _digest(path)
        after = os.stat(path)
        if (st.st_size, st.st_mtime_ns, st.st_ctime_ns) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise RuntimeError(f"Input changed while hashing: {path}")
        with self.db:
            self.db.execute("INSERT OR REPLACE INTO identities VALUES (?, ?, ?)", (path, signature, digest))
        self.memo[key] = {"path": path, "sha256": digest}
        return self.memo[key]

    def get(self, fingerprint, donor, pair):
        row = self.db.execute("SELECT record FROM evidence WHERE fingerprint=? AND donor=? AND u=? AND v=?", (fingerprint, donor, *pair)).fetchone()
        return json.loads(row[0]) if row else None

    def put_batch(self, fingerprint, donor, records):
        with self.db:
            for pair, record in records.items():
                if record["status"] not in ("scored", "no_connection"):
                    raise ValueError("Execution failures cannot be cached")
                self.db.execute("INSERT OR REPLACE INTO evidence VALUES (?, ?, ?, ?, ?)", (fingerprint, donor, *pair, json.dumps(record, sort_keys=True)))

    def close(self):
        self.db.close()


def _worker(task):
    script, inputs, pairs, map_length, threads, minimap2, scratch_root = task
    watched = list(inputs)
    if os.path.isfile(inputs[1] + ".fai"):
        watched.append(inputs[1] + ".fai")
    def signatures():
        return [(p, os.stat(p).st_size, os.stat(p).st_mtime_ns, os.stat(p).st_ctime_ns) for p in watched]
    before = signatures()
    spec = importlib.util.spec_from_file_location("tiebreak_from_paf", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory(prefix="batch-", dir=scratch_root) as scratch:
        _, rows, windows = module.build_connection_evidence(
            inputs[0], inputs[1], inputs[2], os.path.join(scratch, "terminals"),
            threads=threads, MAP_LENGTH=map_length, requested_pairs=set(pairs),
            minimap2_path=minimap2)
    if signatures() != before:
        raise RuntimeError("Tiebreak inputs changed during scoring; results were not cached")
    return {pair: {
        "status": "scored" if pair in rows else "no_connection",
        "score": int(rows[pair][0]) if pair in rows else None,
        "distance": int(rows[pair][1]) if pair in rows else None,
        "windows": {node: windows.get(node) for node in pair},
    } for pair in pairs}


def _atomic_text(path, writer):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name, dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w") as stream:
            writer(stream)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def export_info(root, records, provenance):
    """Partial four-column tables plus a complete request/status audit table."""
    donors = sorted({sample for _, sample in records})
    for sample in donors:
        if Path(sample).name != sample or sample in (".", ".."):
            raise ValueError(f"Unsafe donor filename: {sample!r}")
        def write_rows(stream, sample=sample):
            for (pair, donor), record in sorted(records.items()):
                if donor == sample and record["status"] == "scored":
                    stream.write(f"{pair[0]}\t{pair[1]}\t{record['score']}\t{record['distance']}\n")
        _atomic_text(Path(root) / "info" / (sample + ".info"), write_rows)
    def write_evidence(stream):
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["node_u", "node_v", "donor", "status", "origin", "score", "distance", "windows", "fingerprint", "error"])
        for (pair, donor), record in sorted(records.items()):
            writer.writerow([*pair, donor, record["status"], record["origin"], record.get("score"), record.get("distance"), json.dumps(record.get("windows", {}), sort_keys=True), record.get("fingerprint", ""), record.get("error", "")])
    _atomic_text(Path(root) / "evidence.tsv", write_evidence)
    _atomic_text(Path(root) / "info.metadata.json", lambda stream: json.dump({"schema": SCHEMA_VERSION, "partial": True, "donors": donors, "requests": len(records), "provenance": provenance}, stream, indent=2, sort_keys=True))


def resolve_requests(cfg, requests):
    """Return current request records; failed computations retain legacy fallback behavior.

    Cache-only mode raises before scoring if any requested record is unavailable.
    The cache stores completed no-connection results, but never execution errors.
    """
    root = Path(cfg.tiebreak_dir) / "targeted"
    root.mkdir(parents=True, exist_ok=True)
    cache = EvidenceCache(root / "cache.sqlite3")
    records, missing, provenance = {}, {}, {}
    try:
        if requests:
            minimap2 = shutil.which("minimap2")
            if minimap2 is None:
                raise FileNotFoundError("minimap2 is required to validate targeted tiebreak provenance")
            version = subprocess.check_output([minimap2, "--version"], text=True).strip()
            script = Path(cfg.from_paf_script).resolve()
            code_files = [Path(__file__).resolve(), script, *sorted((script.parent / "ragtag_utilities").glob("*.py"))]
            common = {
                "schema": SCHEMA_VERSION, "code": [cache.identity(p) for p in code_files],
                "minimap2": {"version": version, "binary": cache.identity(minimap2), "options": ["-cx", "asm5"], "threads": cfg.threads},
                "map_length": cfg.map_length,
                "filter": {"terminal_ratio": 0.05, "min_contig_length": 10000, "max_min_id_ratio": 0.5, "connect_mode": "closest", "blacklist": []},
            }
            for donor, pairs in sorted(requests.items()):
                try:
                    if donor not in cfg.sample_inputs:
                        raise FileNotFoundError(f"Donor {donor} is absent from the tiebreak manifest")
                    inputs = cfg.sample_inputs[donor]
                    data = {**common, "inputs": [cache.identity(p) for p in inputs]}
                    if os.path.isfile(inputs[1] + ".fai"):
                        data["query_index"] = cache.identity(inputs[1] + ".fai")
                    fingerprint = hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
                    provenance[donor] = {"fingerprint": fingerprint, **data}
                    for pair in sorted(pairs):
                        cached = cache.get(fingerprint, donor, pair) if cfg.cache_mode != "refresh" else None
                        if cached is not None:
                            records[(pair, donor)] = {**cached, "origin": "cache", "fingerprint": fingerprint}
                        else:
                            missing.setdefault(donor, []).append(pair)
                except OSError as exc:
                    for pair in pairs:
                        records[(pair, donor)] = {"status": "error", "origin": "unavailable", "error": str(exc)}
            unavailable = [(pair, donor) for donor, pairs in missing.items() for pair in pairs]
            unavailable += [key for key, record in records.items() if record["status"] == "error"]
            if cfg.cache_mode == "only" and unavailable:
                raise RuntimeError("Targeted cache-only mode has missing/stale/unavailable requests: " + "; ".join(f"{donor}:{pair[0]},{pair[1]}" for pair, donor in sorted(unavailable)))
            print(f"Tiebreak targeted: {sum(map(len, requests.values()))} requests, {sum(r['origin'] == 'cache' for r in records.values())} cached, {sum(map(len, missing.values()))} to compute across {len(missing)} donors")
            tasks = {donor: (str(script), cfg.sample_inputs[donor], pairs, cfg.map_length, cfg.threads, minimap2, str(root)) for donor, pairs in missing.items()}
            def accept(donor, result=None, error=None):
                fingerprint = provenance[donor]["fingerprint"]
                if error is not None:
                    print(f"Tiebreak targeted: skipping donor {donor} ({error})")
                    result = {pair: {"status": "error", "error": str(error)} for pair in missing[donor]}
                else:
                    cache.put_batch(fingerprint, donor, result)
                for pair, record in result.items():
                    records[(pair, donor)] = {**record, "origin": "computed" if error is None else "failed", "fingerprint": fingerprint}
            if cfg.jobs <= 1 or len(tasks) <= 1:
                for donor, task in tasks.items():
                    try:
                        result = _worker(task)
                    except (OSError, subprocess.CalledProcessError) as exc:
                        accept(donor, error=exc)
                    else:
                        accept(donor, result=result)
            else:
                with ProcessPoolExecutor(max_workers=cfg.jobs) as executor:
                    futures = {executor.submit(_worker, task): donor for donor, task in tasks.items()}
                    for future in as_completed(futures):
                        donor = futures[future]
                        try:
                            result = future.result()
                        except (OSError, subprocess.CalledProcessError) as exc:
                            accept(donor, error=exc)
                        else:
                            accept(donor, result=result)
        for (pair, donor), record in records.items():
            if record["origin"] == "unavailable":
                print(f"Tiebreak targeted: unavailable {donor} {pair}: {record['error']}")
        if cfg.dump_info:
            export_info(root, records, provenance)
        return records
    finally:
        cache.close()
