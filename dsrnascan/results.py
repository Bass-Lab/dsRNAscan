"""Query layer over dsRNAscan results: filter, sort, page, summarize and export, without loading the
sequences unless asked. DuckDB runs the queries straight on the *.dsrnas.parquet file (or, for older
runs, the *.dsrnas.txt / *_merged_results.txt file), reading only the columns each query needs.

Used by dsrna-browse (web UI and command-line export)."""
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

RESULT_PATTERNS = ('*.dsrnas.parquet', '*.dsrnas.txt', '*_merged_results.txt')

# derived columns, computed in the view from the stored ones
DERIVED = {
    'loop_length': 'j_start - i_end',
    'i_length': 'i_end - i_start + 1',
    'j_length': 'j_end - j_start + 1',
    'span': 'j_end - i_start + 1',
}
EXPORT_FORMATS = {'tsv': '.tsv', 'txt': '.txt', 'parquet': '.parquet', 'bedpe': '.bedpe', 'bed': '.bed',
                  'gff3': '.gff3', 'bp': '.bp', 'fasta': '.fa', 'dbn': '.dbn'}
BP_HEADER = ("color:\t100\t149\t237\t70-80% paired (forward strand)\n"
             "color:\t65\t105\t225\t80-90% paired (forward strand)\n"
             "color:\t0\t0\t139\t90-100% paired (forward strand)\n"
             "color:\t205\t92\t92\t70-80% paired (reverse strand)\n"
             "color:\t178\t34\t34\t80-90% paired (reverse strand)\n"
             "color:\t139\t0\t0\t90-100% paired (reverse strand)\n")
BEDPE_HEADER = "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\n"


def q(col):
    """Quote a column name for SQL."""
    return '"' + col.replace('"', '""') + '"'


def lit(text):
    """Quote a string literal for SQL."""
    return "'" + str(text).replace("'", "''") + "'"


def find_results(path):
    """Results files under `path` (a folder, or one file). A *.dsrnas.txt next to its parquet is skipped."""
    p = Path(path)
    if p.is_file():
        return [str(p)]
    files = []
    for pat in RESULT_PATTERNS:
        files += sorted(glob.glob(str(p / pat)))
    parquets = {f[:-len('.parquet')] for f in files if f.endswith('.parquet')}
    return [f for f in files if not (f.endswith('.dsrnas.txt') and f[:-len('.txt')] in parquets)]


def natural_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', str(s))]


def parse_editing_file(path):
    """RNA editing sites from a BED (0-based start) or GFF3 file as a DataFrame(chrom, strand, pos, freq);
    pos is 1-based. Frequencies come from the score column or a frequency/freq/editing_level attribute,
    scaled to 0-1."""
    rows = []
    is_gff = str(path).lower().endswith(('.gff', '.gff3', '.gtf'))
    opener = __import__('gzip').open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            t = line.rstrip('\n').split('\t')
            try:
                if is_gff:
                    if len(t) < 9:
                        continue
                    pos, strand, freq = int(t[3]), t[6], 1.0
                    if t[5] != '.':
                        s = float(t[5])
                        freq = s if s <= 1 else s / 100 if s <= 100 else s / 1000
                    for attr in t[8].split(';'):
                        k, _, v = attr.partition('=')
                        if k.strip().lower() in ('frequency', 'freq', 'editing_level', 'confidence'):
                            v = float(v.strip().rstrip('%'))
                            freq = v / 100 if v > 1 else v
                            break
                else:
                    if len(t) < 3:
                        continue
                    pos, strand, freq = int(t[1]) + 1, (t[5] if len(t) > 5 else '+'), 1.0
                    if len(t) > 4 and t[4] not in ('.', ''):
                        s = float(t[4])
                        freq = s if s <= 1 else s / 1000
            except ValueError:
                continue
            rows.append((t[0], strand, pos, freq))
    return pd.DataFrame(rows, columns=['chrom', 'strand', 'pos', 'freq'])


class Results:
    """A dsRNAscan results set with SQL-backed filtering.

    Filters are plain dicts: {'chrom': [...], 'region': (chrom, start, end), 'strand': '+',
    'min': {col: value}, 'max': {col: value}, 'likely_edited': True, 'likely_forms': True,
    'search': text, 'where': 'raw SQL'}. Every key is optional."""

    LIGHT_COLUMNS = ['id', 'Chromosome', 'Strand', 'i_start', 'i_end', 'j_start', 'j_end', 'Score',
                     'dG(kcal/mol)', 'base_pairs', 'percent_paired', 'longest_helix', 'i_length', 'j_length',
                     'loop_length', 'stability_model_score', 'probing_model_score', 'likely_edited',
                     'likely_forms', 'edits_i', 'edits_j', 'edits']

    def __init__(self, path, editing_file=None, threads=None):
        import duckdb
        self.files = find_results(path)
        if not self.files:
            raise FileNotFoundError(f"no dsRNAscan results ({', '.join(RESULT_PATTERNS)}) in {path}")
        self.con = duckdb.connect()
        if threads:
            self.con.execute(f"SET threads = {int(threads)}")
        parquet = [f for f in self.files if f.endswith('.parquet')]
        if len(self.files) == 1 and parquet:
            # one parquet: a lazy view, with the row's position in the file as its id
            self.con.execute(f"CREATE VIEW src AS SELECT file_row_number AS id, * EXCLUDE (file_row_number) "
                             f"FROM read_parquet({lit(self.files[0])}, file_row_number = true)")
        else:
            # several files, or text: materialize once (text has no stable row numbers otherwise)
            parts = []
            for f in self.files:
                parts.append(f"SELECT * FROM read_parquet({lit(f)})" if f.endswith('.parquet')
                             else f"SELECT * FROM read_csv({lit(f)}, delim='\\t', header=true, sample_size=100000)")
            self.con.execute("CREATE TABLE src AS SELECT row_number() OVER () - 1 AS id, * FROM ("
                             + " UNION ALL BY NAME ".join(parts) + ")")
        stored = [r[0] for r in self.con.execute("DESCRIBE src").fetchall()]
        self.stored_columns = [c for c in stored if c != 'id']
        derived = ', '.join(f"{expr} AS {name}" for name, expr in DERIVED.items() if name not in stored)
        self.con.execute(f"CREATE VIEW base AS SELECT *, {derived} FROM src")
        self.editing = None
        if editing_file:
            self._load_editing(editing_file)
        else:
            self.con.execute("CREATE VIEW dsrnas AS SELECT * FROM base")
        self.columns = [(r[0], r[1]) for r in self.con.execute("DESCRIBE dsrnas").fetchall()]
        self.numeric = [c for c, t in self.columns
                        if any(k in t.upper() for k in ('INT', 'DOUBLE', 'FLOAT', 'DECIMAL', 'REAL'))]
        self.names = [c for c, _ in self.columns]

    # ------------------------------------------------------------------ editing sites
    def _load_editing(self, editing_file):
        """Count editing sites per arm (numpy searchsorted per chromosome and strand) and join them in."""
        sites = parse_editing_file(editing_file)
        self.editing = {k: (g['pos'].to_numpy(), g['freq'].to_numpy())
                        for k, g in sites.sort_values('pos').groupby(['chrom', 'strand'])}
        d = self.con.execute('SELECT id, "Chromosome", "Strand", i_start, i_end, j_start, j_end FROM base').fetch_df()
        ei, ej = np.zeros(len(d), dtype=np.int64), np.zeros(len(d), dtype=np.int64)
        for (chrom, strand), g in d.groupby(['Chromosome', 'Strand']):
            pos = self.editing.get((chrom, strand))
            if pos is None:
                continue
            pos = pos[0]
            idx = g.index.to_numpy()
            ei[idx] = np.searchsorted(pos, g['i_end'].to_numpy(), 'right') - np.searchsorted(pos, g['i_start'].to_numpy(), 'left')
            ej[idx] = np.searchsorted(pos, g['j_end'].to_numpy(), 'right') - np.searchsorted(pos, g['j_start'].to_numpy(), 'left')
        edits = pd.DataFrame({'id': d['id'], 'edits_i': ei, 'edits_j': ej})
        self.con.register('edits_df', edits)
        self.con.execute("CREATE TABLE edits AS SELECT * FROM edits_df")
        self.con.execute("CREATE VIEW dsrnas AS SELECT b.*, coalesce(e.edits_i, 0) AS edits_i, "
                         "coalesce(e.edits_j, 0) AS edits_j, coalesce(e.edits_i, 0) + coalesce(e.edits_j, 0) AS edits "
                         "FROM base b LEFT JOIN edits e USING (id)")
        self.n_sites = len(sites)

    def editing_sites_of(self, row):
        """[[position in i_seq + j_seq, frequency], ...] for one dsRNA (dict with eff_* coordinates and
        sequences), for the structure view."""
        if not self.editing:
            return []
        pos, freq = self.editing.get((row['Chromosome'], row['Strand']), (None, None))
        if pos is None:
            return []
        out = []
        offset = 0
        for arm, seq in (('i', row['i_seq']), ('j', row['j_seq'])):
            start, end = int(row.get(f'eff_{arm}_start') or row[f'{arm}_start']), int(row.get(f'eff_{arm}_end') or row[f'{arm}_end'])
            a, b = np.searchsorted(pos, start, 'left'), np.searchsorted(pos, end, 'right')
            for p, f in zip(pos[a:b], freq[a:b]):
                k = int(p - start) if row['Strand'] == '+' else int(end - p)  # minus arms are reverse-complemented
                if 0 <= k < len(seq):
                    out.append([offset + k, float(f)])
            offset += len(seq)
        return out

    # ------------------------------------------------------------------ filters
    def where(self, f):
        """(SQL WHERE clause without the keyword, parameters) for a filter dict."""
        clauses, params = [], []
        chroms = f.get('chrom') or []
        if isinstance(chroms, str):
            chroms = [chroms]
        if chroms:
            clauses.append(f'"Chromosome" IN ({", ".join("?" * len(chroms))})')
            params += list(chroms)
        if f.get('region'):
            chrom, start, end = f['region']
            clauses.append('"Chromosome" = ? AND i_start <= ? AND j_end >= ?')
            params += [chrom, int(end), int(start)]
        if f.get('strand') in ('+', '-'):
            clauses.append('"Strand" = ?')
            params.append(f['strand'])
        for key, op in (('min', '>='), ('max', '<=')):
            for col, val in (f.get(key) or {}).items():
                if col not in self.numeric:
                    raise ValueError(f"unknown numeric column: {col}")
                clauses.append(f"{q(col)} {op} ?")
                params.append(float(val))
        for flag in ('likely_edited', 'likely_forms'):
            if f.get(flag) and flag in self.names:
                clauses.append(f"{q(flag)} = 'Yes'")
        s = (f.get('search') or '').strip()
        if s:
            m = re.match(r'^(.+?):([\d,]+)-([\d,]+)$', s)
            if m:
                clauses.append('"Chromosome" = ? AND i_start <= ? AND j_end >= ?')
                params += [m.group(1), int(m.group(3).replace(',', '')), int(m.group(2).replace(',', ''))]
            elif s.isdigit():
                clauses.append('id = ?')
                params.append(int(s) - 1)
            else:
                clauses.append('"Chromosome" ILIKE ?')
                params.append(f'%{s}%')
        if f.get('where'):
            clauses.append(f"({f['where']})")
        return (' AND '.join(clauses) if clauses else 'TRUE'), params

    @staticmethod
    def filters_from_params(params):
        """A filter dict from URL query parameters ({name: [values]}) or command-line values."""
        g = lambda k: (params.get(k) or [None])[0]
        f = {'chrom': params.get('chrom') or [], 'strand': g('strand'), 'search': g('search'), 'where': g('where'),
             'likely_edited': g('likely_edited') in ('1', 'true', 'yes'), 'likely_forms': g('likely_forms') in ('1', 'true', 'yes'),
             'min': {}, 'max': {}}
        for k, v in params.items():
            if k.startswith('min_') and v and v[0] != '':
                f['min'][k[4:]] = v[0]
            elif k.startswith('max_') and v and v[0] != '':
                f['max'][k[4:]] = v[0]
        if g('region'):
            m = re.match(r'^(.+?):([\d,]+)-([\d,]+)$', g('region').strip())
            if not m:
                raise ValueError("region must look like chr1:1000-2000")
            f['region'] = (m.group(1), int(m.group(2).replace(',', '')), int(m.group(3).replace(',', '')))
        return f

    def order_sql(self, sort, desc):
        if sort not in self.names:
            raise ValueError(f"unknown column: {sort}")
        return f"ORDER BY {q(sort)} {'DESC' if desc else 'ASC'} NULLS LAST, id"

    # ------------------------------------------------------------------ queries
    def count(self, f=None):
        w, p = self.where(f or {})
        return self.con.execute(f"SELECT count(*) FROM dsrnas WHERE {w}", p).fetchone()[0]

    def page(self, f=None, sort='Score', desc=True, offset=0, limit=50, columns=None):
        """A page of rows as a DataFrame, without sequences unless asked."""
        cols = [c for c in (columns or self.LIGHT_COLUMNS) if c in self.names]
        w, p = self.where(f or {})
        sql = f"SELECT {', '.join(map(q, cols))} FROM dsrnas WHERE {w} {self.order_sql(sort, desc)} LIMIT ? OFFSET ?"
        return self.con.execute(sql, p + [int(limit), int(offset)]).fetch_df()

    def get(self, id_):
        d = self.con.execute("SELECT * FROM dsrnas WHERE id = ?", [int(id_)]).fetch_df()
        if d.empty:
            return None
        row = {k: (None if pd.isna(v) else v.item() if hasattr(v, 'item') else v) for k, v in d.iloc[0].items()}
        row['editing_sites'] = self.editing_sites_of(row)
        return row

    def summary(self):
        total = self.con.execute("SELECT count(*) FROM dsrnas").fetchone()[0]
        chroms = self.con.execute('SELECT "Chromosome", count(*) FROM dsrnas GROUP BY 1').fetchall()
        stats = {}
        if total:
            aggs = ', '.join(f"min({q(c)}), max({q(c)})" for c in self.numeric)
            vals = self.con.execute(f"SELECT {aggs} FROM dsrnas").fetchone()
            for k, c in enumerate(self.numeric):
                lo, hi = vals[2 * k], vals[2 * k + 1]
                stats[c] = {'min': None if lo is None else float(lo), 'max': None if hi is None else float(hi)}
        return {
            'files': [os.path.basename(f) for f in self.files], 'total': total,
            'chromosomes': [{'name': c, 'count': n} for c, n in sorted(chroms, key=lambda x: natural_key(x[0]))],
            'columns': [{'name': c, 'type': t, 'numeric': c in self.numeric, **stats.get(c, {})} for c, t in self.columns
                        if c != 'id'],
            'has_editing': self.editing is not None, 'editing_sites': getattr(self, 'n_sites', 0),
        }

    STAT_COLUMNS = ['percent_paired', 'dG(kcal/mol)', 'base_pairs', 'longest_helix', 'Score', 'i_length', 'j_length',
                    'loop_length', 'span', 'stability_model_score', 'probing_model_score', 'edits']

    def histogram(self, col, f=None, bins=30, log=False, lo=None, hi=None):
        """Bin counts over the filtered rows and over all rows with shared edges. The edges span the
        column's full range, or lo..hi if given (values outside are left out); log-spaced over the
        positive values if `log`. `n` and `n_all` count every value of each set, for densities."""
        if col not in self.numeric:
            raise ValueError(f"unknown numeric column: {col}")
        w, p = self.where(f or {})
        expr, cond, cp = q(col), f"{q(col)} IS NOT NULL", []
        if log:
            cond += f" AND {q(col)} > 0"
            expr = f"log10({q(col)})"
        n_all, n = self.con.execute(f"SELECT count(*), count(*) FILTER (WHERE {w}) FROM dsrnas WHERE {cond}", p).fetchone()
        if lo is not None and hi is not None:
            lo, hi = float(lo), float(hi)
            cond += f" AND {q(col)} BETWEEN ? AND ?"
            cp = [lo, hi]
            if log:
                lo, hi = (__import__('math').log10(v) if v > 0 else 0.0 for v in (lo, hi))
        else:
            lo, hi = self.con.execute(f"SELECT min({expr}), max({expr}) FROM dsrnas WHERE {cond}").fetchone()
            if lo is None:
                return {'col': col, 'log': log, 'edges': [], 'counts': [], 'all': [], 'n': n, 'n_all': n_all}
            lo, hi = float(lo), float(hi)
        if hi <= lo:
            hi = lo + 1
        bins = max(1, int(bins))
        rows = self.con.execute(
            f"SELECT CAST(floor(({expr} - ?) / ? * ?) AS INTEGER) AS b, count(*), count(*) FILTER (WHERE {w}) "
            f"FROM dsrnas WHERE {cond} GROUP BY b", [lo, hi - lo, bins] + p + cp).fetchall()  # placeholders in text order
        counts, every = [0] * bins, [0] * bins
        for bi, n_bin_all, n_bin in rows:
            k = min(max(int(bi), 0), bins - 1)
            every[k] += n_bin_all
            counts[k] += n_bin
        edges = [lo + (hi - lo) * k / bins for k in range(bins + 1)]
        if log:
            edges = [10 ** e for e in edges]
        return {'col': col, 'log': log, 'edges': edges, 'counts': counts, 'all': every, 'n': n, 'n_all': n_all}

    def kde(self, col, f=None, log=False, lo=None, hi=None, bandwidth=None, grid=8192, points=2048):
        """Smoothed density of the filtered rows and of all rows: a fine histogram convolved with a
        Gaussian kernel. The bandwidth follows Silverman's rule from each set's spread and size unless
        given (in the column's units, or log10 units on a log axis). Densities are per unit of the
        axis (log10 units if `log`) and normalized by each set's full size, as the histogram's are.
        The curve is computed on `grid` bins and returned averaged down to `points`."""
        h = self.histogram(col, f, grid, log, lo, hi)
        if not h['edges']:
            return {'col': col, 'log': log, 'x': [], 'filtered': [], 'all': [], 'n': h['n'], 'n_all': h['n_all'],
                    'bandwidth': {'filtered': None, 'all': None}}
        edges = np.log10(h['edges']) if log else np.asarray(h['edges'], dtype=float)
        step = float(edges[1] - edges[0])
        centers = (edges[:-1] + edges[1:]) / 2
        down = max(1, len(centers) // points)
        shrink = lambda v: v[:len(v) // down * down].reshape(-1, down).mean(axis=1)
        out = {'col': col, 'log': log, 'x': (10 ** shrink(centers) if log else shrink(centers)).tolist(), 'n': h['n'],
               'n_all': h['n_all'], 'bandwidth': {}}
        for key, counts, total in (('filtered', h['counts'], h['n']), ('all', h['all'], h['n_all'])):
            k = np.asarray(counts, dtype=float)
            m = k.sum()
            if m == 0 or total == 0:
                out[key], out['bandwidth'][key] = [0.0] * (len(centers) // down), None
                continue
            if bandwidth:
                bw = float(bandwidth)
            else:  # Silverman: 0.9 * min(sd, IQR / 1.34) * n^(-1/5), from the binned values
                mean = (k * centers).sum() / m
                sd = float(np.sqrt((k * (centers - mean) ** 2).sum() / m))
                cum = np.cumsum(k) / m
                q1, q3 = (float(np.interp(qv, cum, centers)) for qv in (0.25, 0.75))
                bw = 0.9 * min(sd, (q3 - q1) / 1.34 if q3 > q1 else sd) * m ** (-0.2)
                bw = max(bw, step)
            half = int(np.ceil(4 * bw / step))
            kern = np.exp(-0.5 * (np.arange(-half, half + 1) * step / bw) ** 2)
            kern /= kern.sum() * step  # the kernel integrates to 1 over the axis
            out[key] = shrink(np.convolve(k, kern, mode='same') / total).tolist()
            out['bandwidth'][key] = bw
        return out

    def scatter(self, x, y, f=None, color=None, max_points=20000, bins=80):
        """Points (id, x, y[, color]) when the filtered set is small enough, else a 2-D density grid."""
        for c in (x, y):
            if c not in self.numeric:
                raise ValueError(f"unknown numeric column: {c}")
        w, p = self.where(f or {})
        cond = f"{w} AND {q(x)} IS NOT NULL AND {q(y)} IS NOT NULL"
        n = self.con.execute(f"SELECT count(*) FROM dsrnas WHERE {cond}", p).fetchone()[0]
        if n <= max_points:
            color = color if color in self.names and color not in (x, y) else None
            cols = ['id', x, y] + ([color] if color else [])
            d = self.con.execute(f"SELECT {', '.join(map(q, cols))} FROM dsrnas WHERE {cond}", p).fetch_df()
            tolist = lambda c: d[c].astype(object).where(d[c].notna(), None).tolist()
            return {'mode': 'points', 'n': n, 'x': x, 'y': y, 'color': color, 'id': tolist('id'), 'xv': tolist(x),
                    'yv': tolist(y), 'cv': tolist(color) if color else None}
        xlo, xhi, ylo, yhi = (float(v) for v in self.con.execute(
            f"SELECT min({q(x)}), max({q(x)}), min({q(y)}), max({q(y)}) FROM dsrnas WHERE {cond}", p).fetchone())
        xr, yr = (xhi - xlo) or 1.0, (yhi - ylo) or 1.0
        rows = self.con.execute(
            f"SELECT CAST(floor(({q(x)} - ?) / ? * ?) AS INTEGER), CAST(floor(({q(y)} - ?) / ? * ?) AS INTEGER), count(*) "
            f"FROM dsrnas WHERE {cond} GROUP BY 1, 2", [xlo, xr, bins, ylo, yr, bins] + p).fetchall()
        z = [[0] * bins for _ in range(bins)]
        for bx, by, c in rows:
            z[min(max(int(by), 0), bins - 1)][min(max(int(bx), 0), bins - 1)] += c
        return {'mode': 'density', 'n': n, 'x': x, 'y': y, 'xedges': [xlo + xr * k / bins for k in range(bins + 1)],
                'yedges': [ylo + yr * k / bins for k in range(bins + 1)], 'z': z}

    def genome(self, f=None, bin_size=1_000_000, max_chroms=40):
        """dsRNAs per bin along each chromosome (filtered), chromosomes ordered naturally; the extent of
        each is the last j_end in the full set."""
        w, p = self.where(f or {})
        bin_size = max(1000, int(bin_size))
        ext = dict(self.con.execute('SELECT "Chromosome", max(j_end) FROM dsrnas GROUP BY 1').fetchall())
        rows = self.con.execute(
            f'SELECT "Chromosome", CAST(floor((i_start - 1) / ?) AS INTEGER) AS b, count(*) FROM dsrnas WHERE {w} GROUP BY 1, 2',
            [bin_size] + p).fetchall()
        keep = sorted(ext, key=lambda c: -ext[c])[:max_chroms]
        counts = {c: [0] * (int(ext[c] - 1) // bin_size + 1) for c in keep}
        for c, b, n in rows:
            if c in counts:
                counts[c][min(max(int(b), 0), len(counts[c]) - 1)] += n
        return {'bin': bin_size, 'omitted': len(ext) - len(keep),
                'chromosomes': [{'name': c, 'length': int(ext[c]), 'counts': counts[c]} for c in sorted(keep, key=natural_key)]}

    def groups(self, by, f=None):
        """Summary statistics of the filtered rows per value of a categorical column."""
        if by not in ('Chromosome', 'Strand', 'likely_edited', 'likely_forms') or by not in self.names:
            raise ValueError(f"cannot group by {by}")
        w, p = self.where(f or {})
        metrics = [('n', 'count(*)'), ('median dG', 'median("dG(kcal/mol)")'), ('median % paired', 'median(percent_paired)'),
                   ('median longest helix', 'median(longest_helix)'), ('median i length', 'median(i_length)'),
                   ('median loop length', 'median(loop_length)')]
        if 'stability_model_score' in self.names:
            metrics += [('mean stability score', 'avg(stability_model_score)'),
                        ('% likely edited', "100.0 * avg(CASE WHEN likely_edited = 'Yes' THEN 1 ELSE 0 END)")]
        if 'probing_model_score' in self.names:
            metrics += [('mean probing score', 'avg(probing_model_score)'),
                        ('% likely forms', "100.0 * avg(CASE WHEN likely_forms = 'Yes' THEN 1 ELSE 0 END)")]
        if 'edits' in self.names:
            metrics += [('editing sites', 'sum(edits)'), ('% with editing sites', '100.0 * avg(CASE WHEN edits > 0 THEN 1 ELSE 0 END)')]
        rows = self.con.execute(f"SELECT {q(by)}, {', '.join(e for _, e in metrics)} FROM dsrnas WHERE {w} GROUP BY 1", p).fetchall()
        rows = sorted(rows, key=lambda r: natural_key(r[0]))
        clean = lambda v: None if v is None else float(v) if not isinstance(v, str) else v
        return {'by': by, 'columns': [by] + [m for m, _ in metrics], 'rows': [[clean(v) for v in r] for r in rows]}

    def stats(self, f=None):
        """n, mean, quartiles, min and max of the main columns, for the filtered rows and for all rows."""
        cols = [c for c in self.STAT_COLUMNS if c in self.numeric]
        aggs = ', '.join(f"count({q(c)}), avg({q(c)}), quantile_cont({q(c)}, [0.25, 0.5, 0.75]), min({q(c)}), max({q(c)})" for c in cols)
        out = {'columns': cols}
        for key, (w, p) in (('filtered', self.where(f or {})), ('all', ('TRUE', []))):
            vals = self.con.execute(f"SELECT {aggs} FROM dsrnas WHERE {w}", p).fetchone()
            out[key] = {}
            for k, c in enumerate(cols):
                n, mean, qs, lo, hi = vals[5 * k:5 * k + 5]
                fl = lambda v: None if v is None else float(v)
                out[key][c] = {'n': int(n), 'mean': fl(mean), 'q1': fl(qs[0]) if qs else None, 'median': fl(qs[1]) if qs else None,
                               'q3': fl(qs[2]) if qs else None, 'min': fl(lo), 'max': fl(hi)}
        return out

    # ------------------------------------------------------------------ genome-browser support
    def region_rows(self, f, chrom, start, end, strand=None, limit=20000):
        """Filtered dsRNAs overlapping chrom:start-end (strongest first if more than `limit`)."""
        w, p = self.where(f or {})
        cols = ['id', 'Chromosome', 'Strand', 'i_start', 'i_end', 'j_start', 'j_end', 'Score', 'dG(kcal/mol)', 'percent_paired', 'longest_helix']
        sql = (f'SELECT {", ".join(map(q, cols))} FROM dsrnas WHERE {w} AND "Chromosome" = ? AND i_start <= ? AND j_end >= ?'
               + (' AND "Strand" = ?' if strand in ('+', '-') else '') + ' ORDER BY "Score" DESC LIMIT ?')
        return self.con.execute(sql, p + [chrom, int(end), int(start)] + ([strand] if strand in ('+', '-') else []) + [int(limit)]).fetch_df()

    def editing_sites_in(self, chrom, start, end):
        """[(pos, freq, strand)] of loaded editing sites in chrom:start-end (1-based, inclusive)."""
        out = []
        for strand in ('+', '-'):
            pos, freq = self.editing.get((chrom, strand), (None, None)) if self.editing else (None, None)
            if pos is None:
                continue
            a, b = np.searchsorted(pos, int(start), 'left'), np.searchsorted(pos, int(end), 'right')
            out += [(int(x), float(y), strand) for x, y in zip(pos[a:b], freq[a:b])]
        return sorted(out)


    # ------------------------------------------------------------------ export
    def export(self, path, fmt=None, f=None, sort='id', desc=False, limit=None, batch=100_000):
        """Write the filtered rows to `path` in `fmt` (by extension if not given). Returns the row count."""
        ext = os.path.splitext(path)[1].lower()
        fmt = fmt or {**{v: k for k, v in EXPORT_FORMATS.items()}, '.fasta': 'fasta', '.fna': 'fasta', '.gff': 'gff3'}.get(ext)
        if fmt not in EXPORT_FORMATS:
            raise ValueError(f"unknown format {fmt!r}; use one of {', '.join(EXPORT_FORMATS)}")
        cols = [c for c in self.names if c not in DERIVED]  # the stored columns (+ edits_*), not the derived ones
        w, p = self.where(f or {})
        lim = f" LIMIT {int(limit)}" if limit else ""
        sql = f"SELECT {', '.join(map(q, cols))} FROM dsrnas WHERE {w} {self.order_sql(sort, desc)}{lim}"
        cur = self.con.execute(sql, p)

        def chunked():  # ~batch rows per DataFrame; the last chunk is empty
            while True:
                d = cur.fetch_df_chunk(max(1, batch // 2048))
                if d is None or len(d) == 0:
                    return
                yield d
        chunks = chunked()
        n = 0
        if fmt == 'parquet':
            import pyarrow as pa
            import pyarrow.parquet as pq
            writer = None
            for d in chunks:
                t = pa.Table.from_pandas(d, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(path, t.schema, compression='zstd')
                writer.write_table(t)
                n += len(d)
            if writer is None:  # nothing matched: still write a valid (empty) file
                pq.write_table(pa.Table.from_pandas(cur.fetch_df().iloc[0:0], preserve_index=False), path, compression='zstd')
            else:
                writer.close()
            return n
        with open(path, 'w') as out:
            if fmt in ('tsv', 'txt'):
                out.write('\t'.join(c for c in cols if c != 'id') + '\n')
            elif fmt == 'bedpe':
                out.write(BEDPE_HEADER)
            elif fmt == 'gff3':
                from dsrnascan import __version__
                out.write(f"##gff-version 3\n##source dsRNAscan v{__version__}\n")
            elif fmt == 'bp':
                out.write(BP_HEADER)
            for d in chunks:
                n += len(d)
                out.write(self._lines(d, fmt))
        return n

    @staticmethod
    def _lines(d, fmt):
        c, st = d['Chromosome'].astype(str), d['Strand'].astype(str)
        ids = 'dsRNA_' + c + '_' + (d['id'] + 1).astype(str).str.zfill(6)
        i0, i1, j0, j1 = (d[k].astype(int) for k in ('i_start', 'i_end', 'j_start', 'j_end'))
        S = lambda s: s.astype(str)
        if fmt in ('tsv', 'txt'):
            return d.drop(columns=['id']).to_csv(sep='\t', index=False, header=False)
        if fmt == 'bedpe':
            lines = (c + '\t' + S(i0 - 1) + '\t' + S(i1) + '\t' + c + '\t' + S(j0 - 1) + '\t' + S(j1) + '\t' + ids + '\t'
                     + S(d['Score']) + '\t' + st + '\t' + st)
        elif fmt == 'bed':  # one line per arm
            lines = pd.concat([c + '\t' + S(i0 - 1) + '\t' + S(i1) + '\t' + ids + '_i\t' + S(d['Score']) + '\t' + st,
                               c + '\t' + S(j0 - 1) + '\t' + S(j1) + '\t' + ids + '_j\t' + S(d['Score']) + '\t' + st])
        elif fmt == 'gff3':
            attrs = ('ID=' + ids + ';Name=' + ids + ';dG=' + S(d['dG(kcal/mol)']) + ';percent_paired='
                     + S(d['percent_paired']) + ';longest_helix=' + S(d['longest_helix']))
            parent = (c + '\tdsRNAscan\tmRNA\t' + S(np.minimum(i0, j0)) + '\t' + S(np.maximum(i1, j1)) + '\t' + S(d['Score'])
                      + '\t' + st + '\t.\t' + attrs)
            arm = lambda a, b, name: (c + '\tdsRNAscan\texon\t' + S(a) + '\t' + S(b) + '\t.\t' + st + '\t.\tParent=' + ids
                                      + ';Name=' + name)
            return ''.join(l + '\n' for l in np.column_stack([parent, arm(i0, i1, 'i_arm'), arm(j0, j1, 'j_arm')]).ravel())
        elif fmt == 'bp':
            pp = d['percent_paired'].astype(float)
            color = (pp >= 80).astype(int) + (pp >= 90).astype(int) + (st == '-').astype(int) * 3
            lines = c + '\t' + S(i0) + '\t' + S(i1) + '\t' + S(j0) + '\t' + S(j1) + '\t' + S(color)
        elif fmt == 'fasta':
            loc = c + ':' + S(i0) + '-' + S(j1) + '(' + st + ')'
            return ''.join('>' + ids + '_i ' + loc + '\n' + d['i_seq'].astype(str) + '\n>' + ids + '_j ' + loc + '\n'
                           + d['j_seq'].astype(str) + '\n')
        elif fmt == 'dbn':
            loc = c + ':' + S(i0) + '-' + S(j1) + '(' + st + ')'
            return ''.join('>' + ids + ' ' + loc + ' dG=' + S(d['dG(kcal/mol)']) + '\n' + d['i_seq'].astype(str) + '&'
                           + d['j_seq'].astype(str) + '\n' + d['structure'].astype(str) + '\n')
        return ''.join(l + '\n' for l in lines)


HOSTED_GENOMES_URL = 'https://s3.amazonaws.com/igv.org.genomes/genomes.json'
_hosted = {}


def hosted_genome(genome_id):
    """igv.js's hosted reference for an assembly id (dm6, hg38, mm39, ...), or None. Fetched once."""
    import json
    import urllib.request
    if 'list' not in _hosted:
        try:
            with urllib.request.urlopen(HOSTED_GENOMES_URL, timeout=15) as r:
                _hosted['list'] = {g['id']: g for g in json.load(r) if g.get('fastaURL')}
        except Exception as e:  # offline, or the list moved
            _hosted['list'] = {}
            _hosted['error'] = str(e)
    g = _hosted['list'].get(genome_id)
    if not g:
        return None
    ref = {k: g[k] for k in ('id', 'name', 'fastaURL', 'indexURL', 'cytobandURL', 'aliasURL', 'chromosomeOrder', 'tracks') if k in g}
    return ref


def write_fai(fasta, out):
    """A samtools-style .fai for an uncompressed FASTA (name, length, offset, bases per line, bytes per line)."""
    with open(fasta, 'rb') as f, open(out, 'w') as o:
        name, length, offset, bases_per_line, bytes_per_line, pos, first = None, 0, 0, 0, 0, 0, True
        for line in f:
            if line.startswith(b'>'):
                if name is not None:
                    o.write(f"{name}\t{length}\t{offset}\t{bases_per_line}\t{bytes_per_line}\n")
                name, length, offset, bases_per_line, bytes_per_line, first = line[1:].split(None, 1)[0].decode(), 0, pos + len(line), 0, 0, True
            elif name is not None:
                n = len(line.rstrip(b'\r\n'))
                if first and n:
                    bases_per_line, bytes_per_line, first = n, len(line), False
                length += n
            pos += len(line)
        if name is not None:
            o.write(f"{name}\t{length}\t{offset}\t{bases_per_line}\t{bytes_per_line}\n")
    return out


def input_fasta_from_log(results_path):
    """The FASTA a run was given, from its dsrnascan.log (None if unknown)."""
    folder = Path(results_path) if Path(results_path).is_dir() else Path(results_path).parent
    log = folder / 'dsrnascan.log'
    if not log.exists():
        return None
    with open(log) as f:
        names = [line.split('filename:', 1)[1].strip() for line in f if '  filename: ' in line]
    if not names:
        return None
    name = names[-1]
    for cand in (Path(name), folder / name, folder.parent / name, Path.cwd() / name):
        if cand.is_file():
            return str(cand)
    return name  # not found as a file; still useful for the message
