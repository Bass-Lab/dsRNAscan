#!/usr/bin/env python3
"""
dsRNAscan - A tool for genome-wide prediction of double-stranded RNA structures
Copyright (C) 2024 Bass Lab
"""

__version__ = '0.7.1'
__author__ = 'Bass Lab'

import argparse
import gzip
import itertools
import logging
import os
import queue
import re
import sys
import threading
import time
import types

import numpy as np
import pandas as pd
from Bio import SeqIO

def available_cpus():
    """CPUs this process may run on: a job's allocation (hyperthreads included), else the machine."""
    if hasattr(os, 'process_cpu_count'):  # Python 3.13+, also honours CPU quotas
        return os.process_cpu_count() or 1
    if hasattr(os, 'sched_getaffinity'):  # Linux: what SLURM/cgroups allow
        return len(os.sched_getaffinity(0)) or 1
    return os.cpu_count() or 1


def default_cpus():
    """A scheduler job gets its whole allocation; a shared machine or laptop gets at most 4."""
    in_job = any(v in os.environ for v in ('SLURM_JOB_ID', 'PBS_JOBID', 'LSB_JOBID', 'JOB_ID'))
    return available_cpus() if in_job else min(4, available_cpus())


BATCH_BP = 5_000_000  # sequences are scanned, scored and written in batches of at least this many bp


# extra= for log lines that go to dsrnascan.log only, not the console
LOG_ONLY = {'log_only': True}


def output_stem(args):
    """Short, filesystem-safe name for output files: the input name without FASTA extensions, plus
    --only_seq if given, or --output_label. The run's settings are in dsrnascan.log."""
    stem = os.path.basename(args.filename)
    if stem.endswith('.gz'):
        stem = stem[:-3]
    root, ext = os.path.splitext(stem)
    if ext.lower() in ('.fa', '.fasta', '.fsa', '.fna', '.fas'):
        stem = root
    # sequence IDs can hold |, : and the like; keep names well under 255 chars
    safe = lambda s, n: re.sub(r'[^A-Za-z0-9._-]+', '_', s).strip('._')[:n]
    if args.output_label != 'header':
        stem = safe(args.output_label, 120)
    else:
        stem = safe(stem, 80)
        if args.only_seq and safe(args.only_seq, 60) != stem:
            stem = f"{stem}.{safe(args.only_seq, 60)}"
    return stem or 'dsrnascan'


def changed_from_defaults(args):
    """Settings that differ from the defaults (the human dsRNAome resource settings), as 'name value (default x)'."""
    as_run = {
        'window': ', '.join(map(str, args.windows)) if args.windows else str(args.w),
        'max_span': 'window' if args.max_span == max(args.windows or [args.w]) else str(args.max_span), 'step': str(args.step), 'score': str(args.score),
        'min_bp': str(args.min_bp), 'paired_cutoff': str(args.paired_cutoff),
        'percent_paired_cutoff': str(args.percent_paired_cutoff), 'match': str(args.match),
        'mismatch': str(args.mismatch), 'gaps': str(args.gaps),
        'fold': 'banded' if args.banded else 'exact', 'temperature': f"{args.t:g} C",
        'minus strand': 'reverse complement' if args.minus_revcomp else 'same windows as plus',
    }
    defaults = {'window': '10000, 3000', 'max_span': 'window', 'step': '150', 'score': '75', 'min_bp': '0',
                'paired_cutoff': '0', 'percent_paired_cutoff': '70', 'match': '3', 'mismatch': '-4',
                'gaps': '12', 'fold': 'banded', 'temperature': '37 C', 'minus strand': 'reverse complement'}
    return [f"{k} {v} (default {defaults[k]})" for k, v in as_run.items() if v != defaults[k]]


def fmt_bp(n):
    return f"{n / 1e6:.1f} Mb" if n >= 1e6 else f"{n / 1e3:.1f} kb" if n >= 1e3 else f"{n} bp"


def fmt_duration(seconds):
    m, sec = divmod(int(round(seconds)), 60)
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m{sec:02d}s" if h else f"{m}m{sec:02d}s" if m else f"{sec}s"


def resource_summary(t_start, threads):
    """One line on CPUs available/used, CPU vs wall time and peak memory, for the log."""
    wall = time.time() - t_start
    try:
        import resource
        own, kids = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = own.ru_utime + own.ru_stime + kids.ru_utime + kids.ru_stime
        peak = max(own.ru_maxrss, kids.ru_maxrss) * (1 if sys.platform == 'darwin' else 1024)  # bytes
        return (f"Resources: {threads} threads used, {available_cpus()} CPUs available; CPU time "
                f"{fmt_duration(cpu)} over {fmt_duration(wall)} wall = {cpu / max(wall, 1e-9):.1f} CPUs busy on "
                f"average; peak memory {peak / 1e9:.1f} GB")
    except ImportError:  # no resource module (Windows)
        return f"Resources: {threads} threads used, {available_cpus()} CPUs available; wall {fmt_duration(wall)}"


def smart_open(filename, mode='rt'):
    """
    Open a file, automatically detecting if it's gzipped based on extension.
    
    Args:
        filename: Path to file
        mode: Mode to open file (default 'rt' for reading text)
        
    Returns:
        File handle (either regular file or gzip file)
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def read_fai(filename):
    """{name: (length, offset, line bases, line width)} from a samtools .fai next to an uncompressed FASTA."""
    fai = filename + ".fai"
    if filename.endswith(".gz") or not os.path.exists(fai) or os.path.getmtime(fai) < os.path.getmtime(filename):
        return None
    with open(fai) as f:
        return {r[0]: tuple(map(int, r[1:5])) for r in (line.split("\t") for line in f) if len(r) >= 5}


def fai_record(filename, name):
    """One sequence read straight from its .fai offset (no parse of the rest of the genome), or None."""
    index = read_fai(filename)
    if not index or name not in index:
        return None
    length, offset, line_bases, line_width = index[name]
    with open(filename, "rb") as f:
        f.seek(offset)
        raw = f.read(length + (length // line_bases + 1) * (line_width - line_bases))
    seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:length].decode()
    return types.SimpleNamespace(id=name, seq=seq)


def _ids(df, first):
    n = pd.Series(range(first, first + len(df)), index=df.index).astype(str).str.zfill(6)
    return "dsRNA_" + df["Chromosome"].astype(str) + "_" + n


def gff3_lines(df, first):
    """GFF3 with mRNA/exon types, so IGV draws the two arms as one linked feature with the loop as
    an intron. `first` numbers the dsRNAs across batches."""
    ids, c, st = _ids(df, first), df["Chromosome"].astype(str), df["Strand"].astype(str)
    i0, i1, j0, j1 = (df[k].astype(int) for k in ("i_start", "i_end", "j_start", "j_end"))
    score = df["Score"].astype(str) if "Score" in df.columns else "."
    attrs = ("ID=" + ids + ";Name=" + ids + ";dG=" + df["dG(kcal/mol)"].astype(str)
             + ";percent_paired=" + df["percent_paired"].astype(str) + ";longest_helix=" + df["longest_helix"].astype(str))
    parent = (c + "\tdsRNAscan\tmRNA\t" + np.minimum(i0, j0).astype(str) + "\t" + np.maximum(i1, j1).astype(str)
              + "\t" + score + "\t" + st + "\t.\t" + attrs)
    arm = lambda a, b, name: (c + "\tdsRNAscan\texon\t" + a.astype(str) + "\t" + b.astype(str) + "\t.\t" + st
                              + "\t.\tParent=" + ids + ";Name=" + name)
    return "".join(l + "\n" for l in np.column_stack([parent, arm(i0, i1, "i_arm"), arm(j0, j1, "j_arm")]).ravel())


def bedpe_lines(df, first):
    """BEDPE: one line per dsRNA with both arms (0-based starts); details are in *.dsrnas.txt."""
    c, st = df["Chromosome"].astype(str), df["Strand"].astype(str)
    score = df["Score"].astype(str) if "Score" in df.columns else "."
    lines = (c + "\t" + (df["i_start"].astype(int) - 1).astype(str) + "\t" + df["i_end"].astype(int).astype(str) + "\t"
             + c + "\t" + (df["j_start"].astype(int) - 1).astype(str) + "\t" + df["j_end"].astype(int).astype(str) + "\t"
             + _ids(df, first) + "\t" + score + "\t" + st + "\t" + st)
    return "".join(l + "\n" for l in lines)


def bp_lines(d, strand):
    """IGV .bp arcs for one strand, coloured by percent paired (header colours 0-2 plus, 3-5 minus)."""
    pp = d["percent_paired"].astype(float)
    color = (pp >= 80).astype(int) + (pp >= 90).astype(int) + (3 if strand == "-" else 0)
    lines = (d["Chromosome"].astype(str) + "\t" + d["i_start"].astype(int).astype(str) + "\t"
             + d["i_end"].astype(int).astype(str) + "\t" + d["j_start"].astype(int).astype(str) + "\t"
             + d["j_end"].astype(int).astype(str) + "\t" + color.astype(str))
    return "".join(l + "\n" for l in lines)


BP_HEADER = ("color:\t100\t149\t237\t70-80% paired (forward strand)\n"
             "color:\t65\t105\t225\t80-90% paired (forward strand)\n"
             "color:\t0\t0\t139\t90-100% paired (forward strand)\n"
             "color:\t205\t92\t92\t70-80% paired (reverse strand)\n"
             "color:\t178\t34\t34\t80-90% paired (reverse strand)\n"
             "color:\t139\t0\t0\t90-100% paired (reverse strand)\n")
BEDPE_HEADER = "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\n"


class OutputFiles:
    """*.dsrnas.txt plus the .bp/.bedpe/.gff3 tracks, appended batch by batch so results are never
    all held in memory. A track file is created when its first row arrives."""

    def __init__(self, output_file, fmt):
        self.txt, self.base, self.fmt = output_file, output_file[:-len('.txt')], fmt
        self.n, self.handles, self.txt_handle = 0, {}, None

    def _handle(self, sfx, header):
        if sfx not in self.handles:
            self.handles[sfx] = open(self.base + sfx, "w")
            self.handles[sfx].write(header)
        return self.handles[sfx]

    def write(self, out):
        if self.txt_handle is None:
            self.txt_handle = open(self.txt, "w")
            out.to_csv(self.txt_handle, sep="\t", index=False)
        else:
            out.to_csv(self.txt_handle, sep="\t", index=False, header=False)
        if out.empty:
            return
        first = self.n + 1
        for sfx, strand in ((".forward.bp", "+"), (".reverse.bp", "-")):
            d = out[out["Strand"] == strand]
            if not d.empty:
                self._handle(sfx, BP_HEADER).write(bp_lines(d, strand))
        if self.fmt in ("bedpe", "both"):
            self._handle(".bedpe", BEDPE_HEADER).write(bedpe_lines(out, first))
        if self.fmt in ("gff3", "both"):
            self._handle(".gff3", f"##gff-version 3\n##source dsRNAscan v{__version__}\n").write(gff3_lines(out, first))
        self.n += len(out)

    def close(self):
        """Close everything; returns the track suffixes written, in a fixed order."""
        for h in list(self.handles.values()) + ([self.txt_handle] if self.txt_handle else []):
            h.close()
        return [s for s in (".forward.bp", ".reverse.bp", ".bedpe", ".gff3") if s in self.handles]


def count_base_pairs(structure):
    """Count actual base pairs in a structure string"""
    if not structure:
        return 0
    # Count '(' characters (each represents a base pair)
    # Don't double count - just count one side
    return structure.count('(')


def scanned_length(args, length):
    """bp of a sequence of `length` that --start/--end leave to scan."""
    if args.start > 0 or args.end > 0:
        a, b = max(0, args.start - 1), (args.end if args.end > 0 else length)
        return max(0, min(b, length) - a)
    return length


def sequence_lengths(args):
    """[(name, length)] of the sequences to scan, in file order, from the .fai or a quick pass over
    the file (headers and line lengths only), so progress can be reported while reading lazily."""
    index = read_fai(args.filename)
    if index is not None:
        lens = [(name, v[0]) for name, v in index.items()]
    else:
        lens, name, n = [], None, 0
        with smart_open(args.filename) as f:
            for line in f:
                if line.startswith(">"):
                    if name is not None:
                        lens.append((name, n))
                    name, n = (line[1:].split(None, 1) or [""])[0], 0
                elif name is not None:
                    n += len(line.strip())
        if name is not None:
            lens.append((name, n))
    if args.only_seq:
        lens = [x for x in lens if x[0] == args.only_seq]
    return lens


def iter_records(args):
    """(name, offset, sequence bytes, full length) one sequence at a time, with --only_seq/--start/--end
    applied; offset is where the scanned part starts in the full sequence."""
    indexed = fai_record(args.filename, args.only_seq) if args.only_seq else None
    with smart_open(args.filename) as handle:
        for record in ([indexed] if indexed else SeqIO.parse(handle, "fasta")):
            if args.only_seq and record.id != args.only_seq:
                continue
            seq = str(record.seq).upper()
            full, offset = len(seq), 0
            if args.start > 0 or args.end > 0:
                offset = max(0, args.start - 1)
                seq = seq[offset:args.end if args.end > 0 else len(seq)]
            yield record.id, offset, seq.encode(), full


ENGINE_COLUMNS = {'Chromosome': 'chromosome', 'Strand': 'strand', 'Score': 'score', 'RawMatch': 'raw_match',
                  'PercMatch': 'match_perc', 'Gaps': 'gap_numb', 'dG(kcal/mol)': 'energy'}


def scan_batches(args, handle):
    """Scan with the compiled engine (dsrnascan._core, built from dsrnascan-core/), reading the FASTA
    lazily, and pass each piece of results to `handle(DataFrame)` in output order: a sequence of at
    least BATCH_BP is streamed a stretch of BATCH_BP at a time; smaller ones are scanned together."""
    from dsrnascan import _core as dsrnascan_core
    logger = logging.getLogger(__name__)
    lens = sequence_lengths(args)
    if not lens:
        sys.exit(f"No sequence named '{args.only_seq}' in {args.filename}" if args.only_seq
                 else f"No sequences in {args.filename}")
    sizes = [scanned_length(args, n) for _, n in lens]
    total = sum(sizes)
    if total == 0:
        sys.exit(f"Nothing to scan: --start {args.start} is past the end of the sequence" if args.start
                 else f"Nothing to scan: the sequences in {args.filename} are empty")
    strands = 'plus' if args.forward_only else 'minus' if args.reverse_only else 'both'
    logger.info(f"Scanning {len(lens):,} sequence{'s' if len(lens) > 1 else ''} ({fmt_bp(total)}), "
                f"{'both strands' if strands == 'both' else strands + ' strand'}, {args.cpus} CPUs")
    # Batches of at least BATCH_BP: one large sequence, or several small ones scanned together.
    # Each sequence is scanned on its own either way, so batching does not change the results.
    plan, cur, cur_bp = [], 0, 0
    for bp in sizes:
        cur, cur_bp = cur + 1, cur_bp + bp
        if cur_bp >= BATCH_BP:
            plan.append(cur)
            cur, cur_bp = 0, 0
    if cur:
        plan.append(cur)

    found = 0

    def piece(cols):
        nonlocal found
        df = pd.DataFrame(cols).rename(columns=ENGINE_COLUMNS)
        found += len(df)
        handle(df)

    def run(batch, stream):
        cols = dsrnascan_core.scan(
            [r[:3] for r in batch], strands=strands, threads=args.cpus, window=args.w,
            max_span=min(args.w, args.max_span), step=args.step, threshold=args.score,
            match_score=args.match, mismatch=args.mismatch, gap=args.gaps,
            min_match_perc=args.paired_cutoff, min_percent_paired=args.percent_paired_cutoff,
            min_bp=args.min_bp, eliminate_nested=not args.no_eliminate_nested, banded=args.banded,
            windows=args.windows, minus_revcomp=args.minus_revcomp, lengths=[r[3] for r in batch],
            temperature=args.t, on_rows=piece if stream else None, chunk_bp=BATCH_BP)
        if not stream:
            piece(cols)

    t0, done = time.time(), 0
    records = iter_records(args)
    for k, size in enumerate(plan):
        # the last batch takes whatever is left, in case the quick length pass and the parser differ
        batch = list(records) if k == len(plan) - 1 else list(itertools.islice(records, size))
        if not batch:
            break
        tb, before = time.time(), found
        run(batch, stream=len(batch) == 1 and len(batch[0][2]) >= BATCH_BP)
        bp = sum(len(r[2]) for r in batch)
        done += bp
        if len(plan) > 1:
            name = batch[0][0] if len(batch) == 1 else f"{batch[0][0]} .. {batch[-1][0]} ({len(batch):,} sequences)"
            elapsed = time.time() - t0
            left = f", ~{fmt_duration(elapsed / done * (total - done))} left" if done < total else ""
            logger.info(f"  {name}: {fmt_bp(bp)}, {found - before:,} dsRNAs in {fmt_duration(time.time() - tb)} "
                        f"({done / total:.0%} done{left})")
        del batch
    logger.info(f"Scan done: {found:,} dsRNAs in {fmt_duration(time.time() - t0)} "
                f"(compiled engine {dsrnascan_core.__version__}{', banded fold' if args.banded else ''})")


def output_table(results, stem):
    """dsRNAscan's output columns for one batch of engine results (sequences and structure last)."""
    n = len(results)

    def col(name, default):
        return results[name] if name in results.columns else pd.Series([default] * n, index=results.index, dtype=object)
    structure = col('structure', '')
    out = pd.DataFrame({
        'Chromosome': col('chromosome', stem), 'Strand': col('strand', '+'),
        'i_start': col('i_start', 0).astype(int), 'i_end': col('i_end', 0).astype(int),
        'j_start': col('j_start', 0).astype(int), 'j_end': col('j_end', 0).astype(int),
        'Score': col('score', 0), 'RawMatch': col('raw_match', '0/0'), 'PercMatch': col('match_perc', 0),
        'Gaps': col('gap_numb', 0), 'dG(kcal/mol)': col('energy', 0.0),
        'base_pairs': results['base_pairs'] if 'base_pairs' in results.columns else structure.map(count_base_pairs),
        'percent_paired': col('percent_paired', 0), 'longest_helix': col('longest_helix', 0),
        'eff_i_start': col('eff_i_start', 0).astype(int), 'eff_i_end': col('eff_i_end', 0).astype(int),
        'eff_j_start': col('eff_j_start', 0).astype(int), 'eff_j_end': col('eff_j_end', 0).astype(int),
    })
    return out, (col('i_seq', ''), col('j_seq', ''), structure)


def run_scan(args):
    """Scan, score and write the outputs, one batch at a time."""
    logger = logging.getLogger(__name__)
    t_start = time.time()
    stem = output_stem(args)

    # Output folder: --output-dir, else <input>_dsrnascan (numbered rather than mixing two runs)
    output_dir = args.output_dir
    if not output_dir:
        output_dir, k = f"{stem}_dsrnascan", 2
        while os.path.isdir(output_dir) and os.listdir(output_dir):
            output_dir, k = f"{stem}_dsrnascan_{k}", k + 1
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, 'dsrnascan.log')

    root_logger = logging.getLogger()
    root_logger.handlers = []
    root_logger.setLevel(logging.INFO)
    # File handler with timestamps; appends, so several runs into one folder each keep their settings
    file_handler = logging.FileHandler(log_file, mode='a')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    root_logger.addHandler(file_handler)
    # Console: messages only, without the LOG_ONLY detail
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(message)s'))
    console_handler.addFilter(lambda record: not getattr(record, 'log_only', False))
    root_logger.addHandler(console_handler)

    logger.info(f"dsRNAscan {__version__}: {args.filename} -> {output_dir}/ (settings in dsrnascan.log)")
    logger.info(f"Command: {' '.join(sys.argv)}", extra=LOG_ONLY)
    logger.info("Settings as run:", extra=LOG_ONLY)
    for arg, value in sorted(vars(args).items()):
        logger.info(f"  {arg}: {value}", extra=LOG_ONLY)
    windows = ', '.join(map(str, args.windows)) if args.windows else args.w
    logger.info(f"Key settings: score {args.score}, min_bp {args.min_bp}, paired_cutoff {args.paired_cutoff}, "
                f"percent_paired_cutoff {args.percent_paired_cutoff}, window {windows}, step {args.step}"
                f"{', banded fold' if args.banded else ''}", extra=LOG_ONLY)
    changed = changed_from_defaults(args)
    if changed:
        logger.info(f"Changed from the default settings: {'; '.join(changed)}")

    status = None
    if not args.no_ml:
        from dsrnascan.ml_scoring import score_results, check_ml_dependencies
        status = check_ml_dependencies()
    use_ml = bool(status and (status['stability'] or status['probing']))
    files = OutputFiles(os.path.join(output_dir, f"{stem}.dsrnas.txt"), args.format)
    counts = {'n': 0}

    def process(results):
        """Score and write one piece of results (runs on the writer thread)."""
        out, (i_seq, j_seq, structure) = output_table(results, stem)
        if out.empty:
            return
        if use_ml:
            out = score_results(out, threads=args.cpus, log=False)
            for c in ('likely_edited', 'likely_forms'):
                if c in out.columns:
                    counts[c] = counts.get(c, 0) + int((out[c] == 'Yes').sum())
        out['i_seq'], out['j_seq'], out['structure'] = i_seq, j_seq, structure
        files.write(out)
        counts['n'] += len(out)

    # The engine hands over results while it keeps scanning; one thread scores and writes them.
    # At most two pieces wait in the queue, so memory stays bounded if writing falls behind.
    pieces, failed = queue.Queue(maxsize=2), []

    def writer():
        while (results := pieces.get()) is not None:
            if not failed:
                try:
                    process(results)
                except BaseException as e:  # reported after the scan stops
                    failed.append(e)

    thread = threading.Thread(target=writer, daemon=True)
    thread.start()

    def handle(results):
        if failed:
            raise failed[0]
        pieces.put(results)

    try:
        scan_batches(args, handle)
    finally:
        pieces.put(None)
        thread.join()
    if failed:
        raise failed[0]
    n = counts['n']
    if files.txt_handle is None:  # no dsRNAs: write the header so the file still opens
        empty, seqs = output_table(pd.DataFrame(), stem)
        empty['i_seq'], empty['j_seq'], empty['structure'] = seqs
        files.write(empty)
    extra = files.close()

    if n and use_ml:
        logger.info("ML scores:")
        if 'likely_edited' in counts:
            logger.info(f"  likely_edited (stability model): {counts['likely_edited']:,} of {n:,} ({counts['likely_edited']/n*100:.1f}%)")
        if 'likely_forms' in counts:
            logger.info(f"  likely_forms (probing model): {counts['likely_forms']:,} of {n:,} ({counts['likely_forms']/n*100:.1f}%)")
    if n and status and status['missing']:
        logger.info(f"  Install {', '.join(status['missing'])} for ML scores")
    logger.info(f"Wrote {stem}.dsrnas.txt" + (f" (+ {', '.join(extra)})" if extra else "") + f" in {output_dir}/")
    logger.info(resource_summary(t_start, args.cpus))
    if n:
        logger.info(f"Browse with: dsrna-browse {output_dir}")
    return n


def main():
    ### Arguments
    parser = argparse.ArgumentParser(
        description='dsRNAscan - A tool for genome-wide prediction of double-stranded RNA structures',
        epilog='Version: {} | Copyright (C) 2024 Bass Lab'.format(__version__)
    )
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))
    parser.add_argument('filename',  type=str,
                        help='input filename')
    parser.add_argument('-t', type=float, default=37,
                        help='Folding temperature in Celsius (default: 37)')
    parser.add_argument('-s', '--step', type=int, default=None,
                        help='Step size (default: 150)')
    parser.add_argument('-w', type=int, default=None,
                        help='Window size (default: windows of 10000 and 3000; -w gives one window of that size)')
    parser.add_argument('--max_span', type=int, default=None,
                        help='Max span of inverted repeat; default = window size')
    parser.add_argument('--min_bp', type=int, default=None,
                        help='Minimum number of base pairs required (overrides --score if set; default: 0)')
    parser.add_argument('--score', type=int, default=None,
                        help='Minimum einverted score (default: 75)')
    parser.add_argument('--gaps', type=int, default=12,
                        help='Gap penalty')
    parser.add_argument('--start', type=int, default=0,
                        help='Starting coordinate for scan; Default = 0')
    parser.add_argument('--end', type=int, default=0,
                        help='Ending coordinate for scan; Default = 0')
    parser.add_argument('-x', '--mismatch', type=int, default=-4,
            help='Mismatch score')
    parser.add_argument('--match', type=int, default=3,
            help='Match score')
    parser.add_argument('--paired_cutoff', type=int, default=None,
                        help='Minimum einverted match percentage before folding (default: 0)')
    parser.add_argument('--percent_paired_cutoff', type=int, default=None,
                        help='Minimum percent paired after folding (default: 70)')
    parser.add_argument('--forward-only', action='store_true', default=False,
                        help='Process forward strand only (default: both strands)')
    parser.add_argument('--reverse-only', action='store_true', default=False,
                        help='Process reverse strand only (default: both strands)')
    parser.add_argument('--output_label', type=str, default='header',
                        help='Label for output files and results (default: use sequence header)')
    parser.add_argument('--only_seq', type=str, default=None,
                        help='Only scan this specific sequence')
    parser.add_argument('-c', '--cpus', type=int, default=None,
                        help='CPUs to use (default: the whole allocation in a SLURM job, otherwise up to 4)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help="Output directory (default: <input>_dsrnascan)")
    parser.add_argument('--no-eliminate-nested', action='store_true', default=False,
                        help='Disable removal of nested dsRNAs (default: removal enabled)')
    parser.add_argument('--v1', action='store_true', default=False,
                        help='On by default: settings used for the human dsRNAome resource')
    parser.add_argument('--exact', action='store_true', default=False,
                        help='Exact fold (slower; the default banded fold gives ~99%% identical results)')
    parser.add_argument('--minus-revcomp', action='store_true', default=False,
                        help=argparse.SUPPRESS)  # part of the default settings; kept for old command lines
    parser.add_argument('--banded', action='store_true', default=False,
                        help=argparse.SUPPRESS)  # the default fold; kept for old command lines
    parser.add_argument('--no-ml', action='store_true', default=False,
                        help='Disable ML scoring of predictions (default: scoring enabled if ydf/sklearn installed)')
    parser.add_argument('--format', type=str, default='bedpe', choices=['bedpe', 'gff3', 'both'],
                        help='Output format: bedpe, gff3, or both (default: bedpe)')

    args = parser.parse_args()
    
    # Get logger for main function
    logger = logging.getLogger(__name__)
    
    # Defaults are the settings used for the human dsRNAome resource (--v1); each flag changes only
    # its own setting (-w gives one window of that size), and the log lists what differs.
    if args.exact and args.banded:
        parser.error("--exact and --banded cannot be combined")
    args.windows = [] if args.w is not None else [10000, 3000]
    args.w = args.w or 10000
    args.step = args.step or 150
    # einverted never searches past the window, so record the span actually used
    args.max_span = min(args.max_span or args.w, max(args.windows or [args.w]))
    if args.paired_cutoff is None:
        args.paired_cutoff = 0
    if args.percent_paired_cutoff is None:
        args.percent_paired_cutoff = 70
    args.minus_revcomp = True
    args.banded = not args.exact

    # min_bp vs score (the resolved values are written to the log)
    if args.min_bp is not None and args.score is None:
        args.score = args.min_bp * args.match
    elif args.min_bp is None and args.score is not None:
        args.min_bp = 0  # a custom score alone means no base-pair floor
    elif args.min_bp is None:
        args.min_bp, args.score = 0, 75
    
    # Validate strand options
    if args.forward_only and args.reverse_only:
        parser.error("Cannot use both --forward-only and --reverse-only at the same time")
    
    # Validate command line arguments
    if args.w <= 0:
        parser.error("Window size must be greater than 0")
    if args.step <= 0:
        parser.error("Step size must be greater than 0")
    if args.step > args.w:
        msg = "Warning: Step size is larger than window size. This may cause gaps in coverage."
        print(msg)
    if args.cpus is None:
        args.cpus = default_cpus()
    if args.cpus <= 0:
        parser.error("Number of CPUs must be greater than 0")
    if args.paired_cutoff < 0 or args.paired_cutoff > 100:
        parser.error("Paired cutoff must be between 0 and 100")
    if args.percent_paired_cutoff < 0 or args.percent_paired_cutoff > 100:
        parser.error("Percent paired cutoff must be between 0 and 100")
    if not 0 <= args.t <= 100:
        parser.error("Folding temperature (-t) must be between 0 and 100 C")
    if args.start < 0:
        parser.error("Start coordinate must be non-negative")
    if args.end < 0:
        parser.error("End coordinate must be non-negative")
    if args.end != 0 and args.end <= args.start:
        parser.error("End coordinate must be greater than start coordinate")
        
    # Check if input file exists and is readable
    if not os.path.exists(args.filename):
        parser.error(f"Input file '{args.filename}' does not exist")
    if not os.access(args.filename, os.R_OK):
        parser.error(f"Input file '{args.filename}' is not readable")
    if os.path.getsize(args.filename) == 0:
        parser.error(f"Input file '{args.filename}' is empty")
    
    # Try to open the file to ensure it's a valid FASTA
    try:
        with smart_open(args.filename) as test_file:
            first_record = next(SeqIO.parse(test_file, "fasta"), None)
            if first_record is None:
                parser.error(f"Input file '{args.filename}' does not appear to be a valid FASTA file")
        index = read_fai(args.filename) if args.only_seq else None
        if args.only_seq and index is not None:
            if args.only_seq not in index:
                parser.error(f"No sequence named '{args.only_seq}' in {args.filename}")
        elif args.only_seq:
            with smart_open(args.filename) as fh:
                if not any(line[1:].split(None, 1)[0] == args.only_seq for line in fh if line.startswith('>') and len(line) > 1):
                    parser.error(f"No sequence named '{args.only_seq}' in {args.filename}")
    except Exception as e:
        parser.error(f"Error reading input file '{args.filename}': {str(e)}")
    
    args.v1 = not changed_from_defaults(args)  # logged: True when the run uses the v1 settings
    run_scan(args)


if __name__ == "__main__":
    main()
