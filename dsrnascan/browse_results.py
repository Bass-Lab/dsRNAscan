#!/usr/bin/env python3
"""dsrna-browse: explore and export dsRNAscan results.

    dsrna-browse results_folder/                       web UI (filters, sorting, structure view, export)
    dsrna-browse results_folder/ --export hits.gff3    filter on the command line and write a file

Reads *.dsrnas.parquet (or the older *.dsrnas.txt / *_merged_results.txt) through dsrnascan.results,
so filtering a whole-genome result set takes well under a second without loading the sequences.
Optional RNA editing sites (BED or GFF3) are counted per arm and shown on the structure."""
import argparse
import json
import os
import socketserver
import sys
import tempfile
import webbrowser
from http.server import SimpleHTTPRequestHandler
from urllib.parse import parse_qs, urlparse

from dsrnascan.results import BP_HEADER, EXPORT_FORMATS, Results, hosted_genome, input_fasta_from_log, write_fai


class Handler(SimpleHTTPRequestHandler):
    results = None   # set before the server starts
    html = ""
    stem = "dsrnascan"
    genome_opt = None   # --genome: an igv.js assembly id or a FASTA path
    path = '.'
    genome = None       # resolved on first use: {'reference': {...}} or {'reference': None, 'message': ...}

    def log_message(self, *args):
        pass

    def _json(self, data, status=200):
        body = json.dumps(data, ensure_ascii=False).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', len(body))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        url = urlparse(self.path)
        params = parse_qs(url.query, keep_blank_values=True)
        try:
            if url.path in ('/', '/index.html'):
                body = self.html.encode()
                self.send_response(200)
                self.send_header('Content-Type', 'text/html; charset=utf-8')
                self.send_header('Content-Length', len(body))
                self.end_headers()
                self.wfile.write(body)
            elif url.path == '/api/summary':
                self._json(self.results.summary())
            elif url.path == '/api/results':
                self._results(params)
            elif url.path.startswith('/api/result/'):
                row = self.results.get(int(url.path.rsplit('/', 1)[1]))
                self._json(row) if row else self._json({'error': 'not found'}, 404)
            elif url.path == '/api/hist':
                f = Results.filters_from_params(params)
                g = lambda k: (params.get(k) or [''])[0]
                lo, hi = (float(g('lo')), float(g('hi'))) if g('lo') != '' and g('hi') != '' else (None, None)
                self._json(self.results.histogram(params['col'][0], f, int(g('bins') or 30), g('log') == '1', lo, hi))
            elif url.path == '/api/kde':
                f = Results.filters_from_params(params)
                g = lambda k: (params.get(k) or [''])[0]
                lo, hi = (float(g('lo')), float(g('hi'))) if g('lo') != '' and g('hi') != '' else (None, None)
                self._json(self.results.kde(params['col'][0], f, g('log') == '1', lo, hi, float(g('bw')) if g('bw') else None))
            elif url.path == '/api/scatter':
                f = Results.filters_from_params(params)
                self._json(self.results.scatter(params['x'][0], params['y'][0], f, (params.get('color') or [None])[0],
                                                int(params.get('max', ['20000'])[0])))
            elif url.path == '/api/genome':
                self._json(self.results.genome(Results.filters_from_params(params), int(params.get('bin', ['1000000'])[0])))
            elif url.path == '/api/groups':
                self._json(self.results.groups(params['by'][0], Results.filters_from_params(params)))
            elif url.path == '/api/stats':
                self._json(self.results.stats(Results.filters_from_params(params)))
            elif url.path == '/api/export':
                self._export(params)
            elif url.path == '/api/igv/config':
                self._json(self._genome())
            elif url.path == '/api/igv/bp':
                self._bp(params)
            elif url.path == '/api/igv/edits':
                self._edits(params)
            elif url.path in ('/data/genome.fa', '/data/genome.fa.fai'):
                g = self._genome()
                if not g.get('local'):
                    self.send_error(404)
                else:
                    self._file(g['local']['fasta'] if url.path.endswith('.fa') else g['local']['fai'])
            else:
                self.send_error(404)
        except (ValueError, KeyError) as e:
            self._json({'error': str(e)}, 400)
        except Exception as e:  # DuckDB errors from a bad SQL filter, etc.
            self._json({'error': f'{type(e).__name__}: {e}'}, 400)

    def _results(self, params):
        g = lambda k, d: (params.get(k) or [d])[0]
        f = Results.filters_from_params(params)
        page, size = int(g('page', 0)), min(int(g('page_size', 50)), 500)
        sort, desc = g('sort', 'Score'), g('dir', 'desc') != 'asc'
        total = self.results.count(f)
        df = self.results.page(f, sort, desc, page * size, size)
        self._json({'total': total, 'page': page, 'page_size': size, 'results': json.loads(df.to_json(orient='records'))})

    def _export(self, params):
        g = lambda k, d: (params.get(k) or [d])[0]
        fmt = g('format', 'tsv')
        f = Results.filters_from_params(params)
        ext = EXPORT_FORMATS[fmt]
        with tempfile.NamedTemporaryFile(suffix=ext, delete=False) as tmp:
            path = tmp.name
        try:
            n = self.results.export(path, fmt, f, g('sort', 'id'), g('dir', 'asc') == 'desc')
            size = os.path.getsize(path)
            self.send_response(200)
            self.send_header('Content-Type', 'application/octet-stream')
            self.send_header('Content-Disposition', f'attachment; filename="{self.stem}.filtered.{n}{ext}"')
            self.send_header('Content-Length', size)
            self.end_headers()
            with open(path, 'rb') as fh:
                while chunk := fh.read(1 << 20):
                    self.wfile.write(chunk)
        finally:
            os.unlink(path)

    # ---------------------------------------------------------------- genome browser (igv.js)
    def _genome(self):
        """The igv.js reference to use: --genome (hosted id or FASTA path), else the results' name if it
        is a hosted assembly, else the run's input FASTA when it is local and uncompressed."""
        cls = self.__class__
        if cls.genome is not None:
            return cls.genome
        opt, stem, log_fa = cls.genome_opt, cls.stem, input_fasta_from_log(cls.path)
        # assembly ids to try: the results name, its first dotted part, and the input FASTA's name
        guesses = [stem, stem.split('.')[0]]
        if log_fa:
            base = os.path.basename(log_fa)
            for ext in ('.gz', '.fa', '.fasta', '.fna'):
                base = base[:-len(ext)] if base.endswith(ext) else base
            guesses += [log_fa, base]
        seen, cands = set(), []
        for c in ([opt] if opt else []) + guesses:
            if c and c not in seen:
                seen.add(c)
                cands.append(c)
        tried = []
        for cand in cands:
            if os.path.isfile(cand):
                if cand.endswith('.gz'):
                    tried.append(f"{cand} is gzipped (igv.js needs an uncompressed FASTA: gunzip it and pass --genome <path>)")
                    continue
                fai = cand + '.fai'
                if not os.path.exists(fai):
                    alt = os.path.join(cls.path if os.path.isdir(cls.path) else os.path.dirname(cls.path), os.path.basename(cand) + '.fai')
                    try:
                        print(f"indexing {cand} ...")
                        write_fai(cand, fai)
                    except OSError:
                        write_fai(cand, alt)
                        fai = alt
                name = os.path.basename(cand)
                cls.genome = {'reference': {'id': name, 'name': name, 'fastaURL': '/data/genome.fa', 'indexURL': '/data/genome.fa.fai'},
                              'source': f'local FASTA {cand}', 'local': {'fasta': cand, 'fai': fai}}
                return cls.genome
            ref = hosted_genome(cand)
            if ref:
                cls.genome = {'reference': ref, 'source': f"igv.js hosted assembly {cand} ({ref.get('name', '')})"}
                return cls.genome
            tried.append(f"{cand}: not a hosted igv.js assembly id" + ("" if os.path.exists(cand) else " and not a local file"))
        cls.genome = {'reference': None, 'message': "No reference genome for the IGV view. Start dsrna-browse with "
                      "--genome <assembly id> (an igv.js hosted genome such as hg38, mm39, dm6, ce11; needs internet) or "
                      "--genome <path to the uncompressed FASTA the scan used>. Tried: " + "; ".join(tried)}
        return cls.genome

    def _bp(self, params):
        """BP arcs (igv.js 'bp' format) for the filtered dsRNAs in a region; stub=1 gives the colour header only."""
        g = lambda k, d=None: (params.get(k) or [d])[0]
        body = BP_HEADER
        if g('stub') != '1':
            f = Results.filters_from_params(params)
            d = self.results.region_rows(f, g('chr'), int(float(g('start', 0))) + 1, int(float(g('end', 0))), g('strand'))
            if len(d):
                body += Results._lines(d, 'bp')
        self._text(body)

    def _edits(self, params):
        """Editing sites in a region as igv.js features (0-based starts)."""
        g = lambda k, d=None: (params.get(k) or [d])[0]
        sites = self.results.editing_sites_in(g('chr'), int(float(g('start', 0))) + 1, int(float(g('end', 0)))) if g('stub') != '1' else []
        self._json([{'chr': g('chr'), 'start': p - 1, 'end': p, 'name': f'{fq:.0%} on {st}', 'score': round(fq * 1000), 'strand': st,
                     'color': '#0B5345' if fq >= 0.8 else '#148F77' if fq >= 0.5 else '#27AE60' if fq >= 0.3 else '#52BE80'}
                    for p, fq, st in sites])

    def _text(self, body):
        body = body.encode()
        self.send_response(200)
        self.send_header('Content-Type', 'text/plain; charset=utf-8')
        self.send_header('Content-Length', len(body))
        self.end_headers()
        self.wfile.write(body)

    def _file(self, path):
        """Serve a file with HTTP range support (igv.js reads an indexed FASTA by byte ranges)."""
        size = os.path.getsize(path)
        start, end, status = 0, size - 1, 200
        rng = self.headers.get('Range', '')
        if rng.startswith('bytes='):
            a, _, b = rng[6:].partition('-')
            start = int(a) if a else max(0, size - int(b))
            end = min(int(b), size - 1) if b and a else end
            status = 206
        self.send_response(status)
        self.send_header('Content-Type', 'application/octet-stream')
        self.send_header('Accept-Ranges', 'bytes')
        self.send_header('Content-Length', end - start + 1)
        if status == 206:
            self.send_header('Content-Range', f'bytes {start}-{end}/{size}')
        self.end_headers()
        with open(path, 'rb') as fh:
            fh.seek(start)
            left = end - start + 1
            while left > 0:
                chunk = fh.read(min(1 << 20, left))
                if not chunk:
                    break
                self.wfile.write(chunk)
                left -= len(chunk)


HTML = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>dsRNA Browser</title>
<script src="https://d3js.org/d3.v3.min.js"></script>
<script src="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.js"></script>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.css">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f0f2f5; color: #1a1a2e; font-size: 13px; }
.header { background: #1a1a2e; color: white; padding: 10px 18px; display: flex; align-items: baseline; gap: 16px; }
.header h1 { font-size: 18px; font-weight: 600; }
.header .sub { color: #9a9ac0; font-size: 12px; }
.layout { display: grid; grid-template-columns: 260px 1fr 460px; height: calc(100vh - 44px); }
.side { background: white; border-right: 1px solid #e0e0e0; overflow-y: auto; padding: 10px 12px; }
.side h3 { font-size: 11px; text-transform: uppercase; letter-spacing: .04em; color: #777; margin: 12px 0 6px; }
.side input[type=text], .side select, .side textarea { width: 100%; padding: 4px 6px; border: 1px solid #ccc; border-radius: 4px; font-size: 12px; }
.side select[multiple] { height: 118px; }
.range { display: grid; grid-template-columns: 1fr 58px 58px; gap: 4px; align-items: center; margin-bottom: 3px; }
.range label { font-size: 12px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.range input { width: 100%; padding: 3px 4px; border: 1px solid #ccc; border-radius: 3px; font-size: 11px; }
.btn { padding: 5px 10px; border: 1px solid #bbb; border-radius: 4px; background: #f7f7f7; cursor: pointer; font-size: 12px; }
.btn:hover { background: #e9e9e9; }
.btn.primary { background: #2d5be3; color: white; border-color: #2d5be3; }
.btn.primary:hover { background: #2148c0; }
.row { display: flex; gap: 6px; align-items: center; margin-top: 6px; }
.center { display: flex; flex-direction: column; overflow: hidden; background: white; }
.stats { display: flex; gap: 18px; padding: 6px 14px; border-bottom: 1px solid #e0e0e0; align-items: flex-end; background: #fafbfc; }
.stats .count { font-size: 14px; font-weight: 600; min-width: 170px; }
.stats .count small { font-weight: normal; color: #777; }
.hist { display: flex; flex-direction: column; font-size: 10px; color: #666; }
.hist svg { display: block; }
.hist rect { fill: #7a93d8; } .hist rect:hover { fill: #2d5be3; }
.tablewrap { flex: 1; overflow: auto; }
table { width: 100%; border-collapse: collapse; font-size: 12px; }
th { position: sticky; top: 0; background: #f8f9fa; padding: 7px 8px; text-align: left; border-bottom: 2px solid #dee2e6; cursor: pointer; white-space: nowrap; user-select: none; }
th:hover { background: #e9ecef; }
th.asc::after { content: " \25B2"; font-size: 9px; } th.desc::after { content: " \25BC"; font-size: 9px; }
td { padding: 5px 8px; border-bottom: 1px solid #f0f0f0; white-space: nowrap; }
tbody tr { cursor: pointer; } tbody tr:hover { background: #f0f4ff; } tbody tr.sel { background: #d0e0ff; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
.pager { display: flex; align-items: center; justify-content: center; gap: 8px; padding: 6px; border-top: 1px solid #e0e0e0; background: #fafafa; }
.detail { background: white; border-left: 1px solid #e0e0e0; display: flex; flex-direction: column; overflow: hidden; }
.metrics { padding: 10px 12px; border-bottom: 1px solid #e0e0e0; font-size: 12px; line-height: 1.6; overflow-y: auto; max-height: 45%; }
.metrics b { font-family: monospace; }
.metrics .lab { color: #666; }
#forna { flex: 1; min-height: 300px; }
.muted { color: #999; text-align: center; padding: 30px; }
.err { color: #c0392b; font-size: 12px; margin-top: 6px; white-space: pre-wrap; }
.tag { display: inline-block; padding: 1px 6px; border-radius: 3px; font-size: 11px; font-weight: 600; }
.tag.yes { background: #e3f6ea; color: #1e7b45; } .tag.no { background: #f3f3f3; color: #999; }
details summary { cursor: pointer; font-size: 12px; color: #555; margin-top: 8px; }
.tabs { display: flex; gap: 2px; padding: 6px 10px 0; border-bottom: 1px solid #e0e0e0; background: #fafbfc; }
.tabs button { border: 1px solid transparent; border-bottom: none; background: none; padding: 6px 12px; font-size: 12px; cursor: pointer; color: #555; border-radius: 6px 6px 0 0; }
.tabs button.on { background: white; border-color: #e0e0e0; color: #1a1a2e; font-weight: 600; margin-bottom: -1px; }
.panel { display: none; flex: 1; flex-direction: column; overflow: hidden; }
.panel.on { display: flex; }
.ctrl { display: flex; gap: 10px; align-items: center; padding: 8px 12px; border-bottom: 1px solid #eee; font-size: 12px; flex-wrap: wrap; }
.ctrl select, .ctrl input { padding: 3px 6px; border: 1px solid #ccc; border-radius: 4px; font-size: 12px; }
.ctrl .dl { margin-left: auto; display: flex; gap: 4px; }
.chart { flex: 1; min-height: 300px; }
.panel table.stat { font-size: 12px; } .panel table.stat th { position: static; } .panel .tablewrap { padding: 8px; }
.note { font-size: 11px; color: #777; padding: 4px 12px; }
</style>
</head>
<body>
<div class="header"><h1>dsRNA Browser</h1><div class="sub" id="sub">loading...</div></div>
<div class="layout">
  <div class="side">
    <h3>Search</h3>
    <input type="text" id="search" placeholder="chr1:100000-200000, chromosome, or ID">
    <h3>Chromosome <span style="font-weight:normal;color:#999">(ctrl/cmd-click for several)</span></h3>
    <select id="chrom" multiple></select>
    <h3>Strand</h3>
    <select id="strand"><option value="">both</option><option value="+">+ (forward)</option><option value="-">- (reverse)</option></select>
    <div id="flags"></div>
    <h3>Ranges <span style="font-weight:normal;color:#999">(blank = any)</span></h3>
    <div id="ranges"></div>
    <details><summary>Advanced: SQL filter</summary>
      <textarea id="where" rows="3" placeholder='e.g. "dG(kcal/mol)" < -200 AND percent_paired > 85'></textarea>
    </details>
    <div class="row"><button class="btn primary" id="apply">Apply</button><button class="btn" id="reset">Reset</button></div>
    <div class="err" id="err"></div>
    <h3>Export the filtered set</h3>
    <div class="row">
      <select id="fmt">
        <option value="tsv">TSV (all columns)</option><option value="parquet">Parquet</option>
        <option value="bedpe">BEDPE (arm pairs)</option><option value="bed">BED (arms)</option>
        <option value="gff3">GFF3 (IGV)</option><option value="bp">BP arcs (IGV)</option>
        <option value="fasta">FASTA (arm sequences)</option><option value="dbn">Dot-bracket</option>
      </select>
      <button class="btn" id="export">Export</button>
    </div>
    <div style="color:#777;font-size:11px;margin-top:4px" id="exportinfo"></div>
  </div>
  <div class="center">
    <div class="stats"><div class="count" id="count">...</div><div id="hists" style="display:flex;gap:14px"></div></div>
    <div class="tabs" id="tabs">
      <button data-tab="table" class="on">Table</button><button data-tab="dist">Distributions</button><button data-tab="scatter">Scatter</button>
      <button data-tab="genome">Genome</button><button data-tab="groups">Groups</button><button data-tab="stats">Summary</button><button data-tab="igv">IGV</button>
    </div>
    <div class="panel on" id="tab-table">
      <div class="tablewrap"><table><thead><tr id="head"></tr></thead><tbody id="body"></tbody></table></div>
      <div class="pager">
        <button class="btn" id="first">&laquo;</button><button class="btn" id="prev">&lsaquo;</button>
        <span id="pageinfo"></span>
        <button class="btn" id="next">&rsaquo;</button><button class="btn" id="last">&raquo;</button>
        <select id="pagesize"><option>25</option><option selected>50</option><option>100</option><option>200</option></select>
      </div>
    </div>
    <div class="panel" id="tab-dist">
      <div class="ctrl"><label>column <select id="dist-col"></select></label>
        <label>plot <select id="dist-kind"><option value="kde">smoothed (KDE)</option><option value="hist">histogram</option></select></label>
        <label id="dist-y-wrap">y <select id="dist-y"><option value="density">density</option><option value="counts">counts</option></select></label>
        <label id="dist-bins-wrap">bins <input type="number" id="dist-bins" value="40" min="5" max="200" style="width:60px"></label>
        <label id="dist-bw-wrap">bandwidth <input type="number" step="any" id="dist-bw" placeholder="auto" style="width:80px"></label>
        <label>x from <input type="number" step="any" id="dist-lo" placeholder="min" style="width:80px"></label><label>to <input type="number" step="any" id="dist-hi" placeholder="max" style="width:80px"></label>
        <label><input type="checkbox" id="dist-log"> log x</label><label><input type="checkbox" id="dist-ylog"> log y</label><span class="note">filtered set in blue over the full set in grey</span>
        <span class="dl"><button class="btn" onclick="dlImage('dist','png')">PNG</button><button class="btn" onclick="dlImage('dist','svg')">SVG</button><button class="btn" onclick="dlData('dist')">TSV</button></span></div>
      <div class="chart" id="chart-dist"></div>
    </div>
    <div class="panel" id="tab-scatter">
      <div class="ctrl"><label>x <select id="sc-x"></select></label><label>y <select id="sc-y"></select></label><label>color <select id="sc-color"><option value="">none</option></select></label>
        <label>points up to <select id="sc-max"><option>5000</option><option selected>20000</option><option>50000</option><option>100000</option><option>200000</option></select></label>
        <span class="note" id="sc-note">box-select points to turn the box into a filter</span>
        <span class="dl"><button class="btn" onclick="dlImage('scatter','png')">PNG</button><button class="btn" onclick="dlImage('scatter','svg')">SVG</button><button class="btn" onclick="dlData('scatter')">TSV</button></span></div>
      <div class="chart" id="chart-scatter"></div>
    </div>
    <div class="panel" id="tab-genome">
      <div class="ctrl"><label>bin <select id="gen-bin"><option value="100000">100 kb</option><option value="500000">500 kb</option><option value="1000000" selected>1 Mb</option><option value="5000000">5 Mb</option></select></label>
        <label>colors <select id="gen-cs"></select></label><label><input type="checkbox" id="gen-rev"> reverse</label>
        <label>scale <select id="gen-zscale"><option value="linear">linear</option><option value="log">log</option></select></label>
        <label>max <input type="number" step="any" id="gen-zmax" placeholder="auto" style="width:70px"></label>
        <span class="note" id="gen-note">dsRNAs per bin along each chromosome; click a bin to filter to that region</span>
        <span class="dl"><button class="btn" onclick="dlImage('genome','png')">PNG</button><button class="btn" onclick="dlImage('genome','svg')">SVG</button><button class="btn" onclick="dlData('genome')">TSV</button></span></div>
      <div class="chart" id="chart-genome" style="overflow-y:auto"></div>
    </div>
    <div class="panel" id="tab-groups">
      <div class="ctrl"><label>group by <select id="gr-by"></select></label><label>chart <select id="gr-metric"></select></label>
        <span class="dl"><button class="btn" onclick="dlImage('groups','png')">PNG</button><button class="btn" onclick="dlImage('groups','svg')">SVG</button><button class="btn" onclick="dlData('groups')">TSV</button></span></div>
      <div class="chart" id="chart-groups" style="flex:0 0 260px"></div>
      <div class="tablewrap" id="table-groups"></div>
    </div>
    <div class="panel" id="tab-igv">
      <div class="ctrl"><span class="note" id="igv-note">loading...</span><span class="dl"><button class="btn" onclick="igvFit()">go to selection</button></span></div>
      <div class="chart" id="igv-div" style="overflow-y:auto;padding:4px"></div>
    </div>
    <div class="panel" id="tab-stats">
      <div class="ctrl"><span class="note">filtered set next to the full set</span><span class="dl"><button class="btn" onclick="dlData('stats')">TSV</button></span></div>
      <div class="tablewrap" id="table-stats"></div>
    </div>
  </div>
  <div class="detail">
    <div class="metrics" id="metrics"><div class="muted">Select a dsRNA</div></div>
    <div id="forna"><div class="muted"></div></div>
  </div>
</div>
<script>
const S = { summary: null, filters: {}, sort: 'Score', dir: 'desc', page: 0, size: 50, total: 0, selected: null };
const RANGE_ORDER = ['percent_paired', 'dG(kcal/mol)', 'base_pairs', 'longest_helix', 'Score', 'i_length', 'j_length',
                     'loop_length', 'span', 'stability_model_score', 'probing_model_score', 'edits', 'edits_i', 'edits_j'];
const LABELS = { 'dG(kcal/mol)': 'dG (kcal/mol)', percent_paired: '% paired', base_pairs: 'base pairs', longest_helix: 'longest helix',
                 i_length: 'i arm length', j_length: 'j arm length', loop_length: 'loop length', span: 'span',
                 stability_model_score: 'stability score', probing_model_score: 'probing score', edits: 'editing sites',
                 edits_i: 'editing sites, i arm', edits_j: 'editing sites, j arm', Score: 'search score' };
const HISTS = ['dG(kcal/mol)', 'percent_paired', 'i_length'];
const COLS = [ ['loc', 'Location', 'i_start', false], ['Strand', 'Strand', 'Strand', false], ['Score', 'Score', 'Score', true],
  ['dG(kcal/mol)', 'dG', 'dG(kcal/mol)', true], ['base_pairs', 'BP', 'base_pairs', true], ['percent_paired', '% paired', 'percent_paired', true],
  ['longest_helix', 'Helix', 'longest_helix', true], ['arms', 'Arms', 'i_length', false], ['loop_length', 'Loop', 'loop_length', true],
  ['stability_model_score', 'Stability', 'stability_model_score', true], ['probing_model_score', 'Probing', 'probing_model_score', true],
  ['edits', 'Edits', 'edits', true] ];
const fmt = n => n == null ? '' : n.toLocaleString();
const num = (v, d) => v == null ? '-' : Number(v).toFixed(d);
async function api(path, params) {
  const r = await fetch(path + (params ? '?' + params : ''));
  const j = await r.json();
  if (!r.ok) throw new Error(j.error || r.statusText);
  return j;
}
function params(extra) {
  const p = new URLSearchParams();
  const f = S.filters;
  (f.chrom || []).forEach(c => p.append('chrom', c));
  for (const k of ['search', 'strand', 'where']) if (f[k]) p.set(k, f[k]);
  for (const k of ['likely_edited', 'likely_forms']) if (f[k]) p.set(k, '1');
  for (const [c, v] of Object.entries(f.min || {})) p.set('min_' + c, v);
  for (const [c, v] of Object.entries(f.max || {})) p.set('max_' + c, v);
  for (const [k, v] of Object.entries(extra || {})) p.set(k, v);
  return p;
}
function readFilters() {
  const f = { min: {}, max: {} };
  f.search = document.getElementById('search').value.trim();
  f.chrom = [...document.getElementById('chrom').selectedOptions].map(o => o.value);
  f.strand = document.getElementById('strand').value;
  f.where = document.getElementById('where').value.trim();
  for (const el of document.querySelectorAll('#flags input')) f[el.dataset.col] = el.checked;
  for (const el of document.querySelectorAll('#ranges input')) if (el.value !== '') f[el.dataset.kind][el.dataset.col] = el.value;
  S.filters = f; S.page = 0;
}
function has(col) { return S.summary.columns.some(c => c.name === col); }
function colInfo(col) { return S.summary.columns.find(c => c.name === col); }
function buildSidebar() {
  const sum = S.summary;
  const sel = document.getElementById('chrom');
  sum.chromosomes.forEach(c => { const o = document.createElement('option'); o.value = c.name; o.textContent = `${c.name}  (${fmt(c.count)})`; sel.appendChild(o); });
  const flags = document.getElementById('flags');
  let h = '';
  if (has('likely_edited')) h += '<label style="display:block;margin-top:8px"><input type="checkbox" data-col="likely_edited"> likely edited (stability model)</label>';
  if (has('likely_forms')) h += '<label style="display:block"><input type="checkbox" data-col="likely_forms"> likely forms (probing model)</label>';
  flags.innerHTML = h ? '<h3>ML flags</h3>' + h : '';
  const ranges = document.getElementById('ranges');
  const numeric = sum.columns.filter(c => c.numeric && c.name !== 'id');
  const ordered = RANGE_ORDER.filter(n => numeric.some(c => c.name === n)).concat(numeric.map(c => c.name).filter(n => !RANGE_ORDER.includes(n) && !['i_start','i_end','j_start','j_end','eff_i_start','eff_i_end','eff_j_start','eff_j_end','PercMatch','Gaps'].includes(n)));
  ranges.innerHTML = ordered.map(n => { const c = colInfo(n); const ph = v => v == null ? '' : (Number.isInteger(v) ? v : v.toFixed(2));
    return `<div class="range"><label title="${n}">${LABELS[n] || n}</label><input type="number" step="any" data-kind="min" data-col="${n}" placeholder="${ph(c.min)}"><input type="number" step="any" data-kind="max" data-col="${n}" placeholder="${ph(c.max)}"></div>`; }).join('');
  const head = document.getElementById('head');
  head.innerHTML = COLS.filter(c => c[0] === 'loc' || c[0] === 'arms' || has(c[2])).map(c =>
    `<th data-sort="${c[2]}" class="${c[3] ? 'num' : ''}">${c[1]}</th>`).join('');
  head.querySelectorAll('th').forEach(th => th.addEventListener('click', () => {
    const col = th.dataset.sort; if (S.sort === col) S.dir = S.dir === 'desc' ? 'asc' : 'desc'; else { S.sort = col; S.dir = ['Score','base_pairs','percent_paired','longest_helix','stability_model_score','probing_model_score','edits'].includes(col) ? 'desc' : 'asc'; }
    S.page = 0; load(); }));
  document.getElementById('hists').innerHTML = HISTS.filter(has).map(c => `<div class="hist" id="hist-${c.replace(/\W/g,'_')}"><svg width="150" height="42"></svg><span>${LABELS[c] || c}</span></div>`).join('');
}
async function load() {
  document.getElementById('err').textContent = '';
  try {
    const data = await api('/api/results', params({ sort: S.sort, dir: S.dir, page: S.page, page_size: S.size }));
    S.total = data.total;
    render(data.results);
    const pages = Math.max(1, Math.ceil(data.total / S.size));
    document.getElementById('count').innerHTML = `${fmt(data.total)} dsRNAs <small>of ${fmt(S.summary.total)}</small>`;
    document.getElementById('pageinfo').textContent = data.total ? `${fmt(S.page * S.size + 1)}-${fmt(Math.min((S.page + 1) * S.size, data.total))}  (page ${S.page + 1} of ${fmt(pages)})` : 'no matches';
    document.getElementById('first').disabled = document.getElementById('prev').disabled = S.page === 0;
    document.getElementById('next').disabled = document.getElementById('last').disabled = S.page >= pages - 1;
    document.getElementById('exportinfo').textContent = `${fmt(data.total)} dsRNAs will be written`;
    document.querySelectorAll('#head th').forEach(th => th.className = (th.classList.contains('num') ? 'num ' : '') + (th.dataset.sort === S.sort ? S.dir : ''));
    HISTS.filter(has).forEach(c => api('/api/hist', params({ col: c, bins: 30 })).then(drawHist).catch(() => {}));
    loadTab();
  } catch (e) { document.getElementById('err').textContent = e.message; }
}
function render(rows) {
  T.rows = rows;
  const body = document.getElementById('body');
  body.innerHTML = '';
  const score = (v, t) => v == null ? '<td class="num" style="color:#bbb">-</td>' : `<td class="num" style="${v >= t ? 'color:#1e7b45;font-weight:600' : 'color:#999'}">${num(v, 3)}</td>`;
  T.firstRow = rows[0] || T.firstRow;
  for (const r of rows) {
    const tr = document.createElement('tr');
    tr.dataset.id = r.id;
    let h = `<td>${r.Chromosome}:${fmt(r.i_start)}-${fmt(r.j_end)}</td><td style="text-align:center">${r.Strand}</td><td class="num">${r.Score}</td>` +
      `<td class="num">${num(r['dG(kcal/mol)'], 1)}</td><td class="num">${r.base_pairs}</td><td class="num">${num(r.percent_paired, 1)}</td>` +
      `<td class="num">${r.longest_helix}</td><td>${r.i_length}/${r.j_length}</td><td class="num">${fmt(r.loop_length)}</td>`;
    if (has('stability_model_score')) h += score(r.stability_model_score, 0.5526);
    if (has('probing_model_score')) h += score(r.probing_model_score, 0.440);
    if (has('edits')) h += `<td class="num">${r.edits || 0}</td>`;
    tr.innerHTML = h;
    if (r.id === S.selected) tr.classList.add('sel');
    tr.addEventListener('click', () => select(r.id, tr));
    body.appendChild(tr);
  }
}
function drawHist(hd) {
  const el = document.querySelector(`#hist-${hd.col.replace(/\W/g,'_')} svg`);
  if (!el || !hd.counts.length) return;
  const w = 150, h = 42, n = hd.counts.length, max = Math.max(1, ...hd.all), bw = w / n;
  el.innerHTML = hd.counts.map((c, i) => { const bh = c / max * (h - 2), ba = hd.all[i] / max * (h - 2);
    return `<rect x="${i * bw}" y="${h - ba}" width="${Math.max(bw - 1, 1)}" height="${ba}" style="fill:#e4e6ea"></rect>` +
      `<rect x="${i * bw}" y="${h - bh}" width="${Math.max(bw - 1, 1)}" height="${bh}"><title>${hd.edges[i].toFixed(1)} to ${hd.edges[i + 1].toFixed(1)}: ${fmt(c)} of ${fmt(hd.all[i])}</title></rect>`; }).join('');
}
// ---------------------------------------------------------------- analysis tabs
const T = { active: 'table', data: {} };
const PL = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['lasso2d'] };
const LAY = { margin: { l: 60, r: 20, t: 10, b: 50 }, font: { family: 'inherit', size: 12 }, paper_bgcolor: 'white', plot_bgcolor: 'white' };
const label = c => LABELS[c] || c;
function fillSelect(id, cols, chosen) { const el = document.getElementById(id); el.innerHTML = cols.map(c => `<option value="${c}" ${c === chosen ? 'selected' : ''}>${label(c)}</option>`).join(''); }
function buildTabs() {
  const numeric = S.summary.columns.filter(c => c.numeric).map(c => c.name).filter(n => !['i_start','i_end','j_start','j_end','eff_i_start','eff_i_end','eff_j_start','eff_j_end','PercMatch','Gaps'].includes(n));
  const ordered = RANGE_ORDER.filter(n => numeric.includes(n)).concat(numeric.filter(n => !RANGE_ORDER.includes(n)));
  fillSelect('dist-col', ordered, 'dG(kcal/mol)'); fillSelect('sc-x', ordered, 'i_length'); fillSelect('sc-y', ordered, 'dG(kcal/mol)');
  const cats = ['Strand', 'likely_edited', 'likely_forms'].filter(has);
  document.getElementById('sc-color').innerHTML = '<option value="">none</option>' + cats.map(c => `<option value="${c}">${label(c)}</option>`).join('');
  fillSelect('gr-by', ['Chromosome'].concat(cats), 'Chromosome');
  document.getElementById('gen-cs').innerHTML = Object.keys(PALETTES).map(k => `<option>${k}</option>`).join('');
  document.querySelectorAll('#tabs button').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('#tabs button').forEach(x => x.classList.toggle('on', x === b));
    document.querySelectorAll('.panel').forEach(p => p.classList.toggle('on', p.id === 'tab-' + b.dataset.tab));
    T.active = b.dataset.tab; loadTab(); }));
  for (const id of ['dist-col', 'dist-kind', 'dist-bins', 'dist-bw', 'dist-log', 'dist-lo', 'dist-hi', 'sc-x', 'sc-y', 'sc-color', 'sc-max', 'gen-bin', 'gr-by']) document.getElementById(id).addEventListener('change', loadTab);
  document.getElementById('dist-kind').addEventListener('change', () => { document.getElementById('dist-bw').value = ''; distControls(); });
  distControls();
  for (const id of ['dist-y', 'dist-ylog']) document.getElementById(id).addEventListener('change', () => drawDist(T.data.dist));
  for (const id of ['gen-cs', 'gen-rev', 'gen-zscale', 'gen-zmax']) document.getElementById(id).addEventListener('change', () => drawGenome(T.data.genome));
  document.getElementById('dist-col').addEventListener('change', () => { document.getElementById('dist-lo').value = ''; document.getElementById('dist-hi').value = ''; });
  document.getElementById('gr-metric').addEventListener('change', () => drawGroups(T.data.groups));
}
async function loadTab() {
  const t = T.active; if (t === 'table') return;
  try {
    if (t === 'dist') { const kde = document.getElementById('dist-kind').value === 'kde'; const lo = document.getElementById('dist-lo').value, hi = document.getElementById('dist-hi').value, bw = document.getElementById('dist-bw').value;
      const extra = { col: document.getElementById('dist-col').value, bins: document.getElementById('dist-bins').value, log: document.getElementById('dist-log').checked ? 1 : 0 }; if (lo !== '' && hi !== '') { extra.lo = lo; extra.hi = hi; } if (kde && bw !== '') extra.bw = bw;
      const d = await api(kde ? '/api/kde' : '/api/hist', params(extra)); d.kind = kde ? 'kde' : 'hist'; T.data.dist = d; drawDist(d); }
    else if (t === 'scatter') { const d = await api('/api/scatter', params({ x: document.getElementById('sc-x').value, y: document.getElementById('sc-y').value, color: document.getElementById('sc-color').value, max: document.getElementById('sc-max').value })); T.data.scatter = d; drawScatter(d); }
    else if (t === 'genome') { const d = await api('/api/genome', params({ bin: document.getElementById('gen-bin').value })); T.data.genome = d; drawGenome(d); }
    else if (t === 'groups') { const d = await api('/api/groups', params({ by: document.getElementById('gr-by').value })); T.data.groups = d; drawGroups(d); }
    else if (t === 'stats') { const d = await api('/api/stats', params()); T.data.stats = d; drawStats(d); }
    else if (t === 'igv') { await showIgv(); }
  } catch (e) { document.getElementById('err').textContent = e.message; }
}
function centers(e) { return e.slice(0, -1).map((v, i) => (v + e[i + 1]) / 2); }
function densities(d) {  // per set: count / (set size x bin width), widths in log10 units on a log axis
  const w = d.edges.slice(0, -1).map((v, i) => d.log ? Math.log10(d.edges[i + 1]) - Math.log10(v) : d.edges[i + 1] - v);
  return { filtered: d.counts.map((c, i) => d.n ? c / (d.n * w[i]) : 0), all: d.all.map((c, i) => d.n_all ? c / (d.n_all * w[i]) : 0) };
}
function distControls() {
  const kde = document.getElementById('dist-kind').value === 'kde';
  document.getElementById('dist-bins-wrap').style.display = kde ? 'none' : '';
  document.getElementById('dist-bw-wrap').style.display = kde ? '' : 'none';
  document.getElementById('dist-y-wrap').style.display = kde ? 'none' : '';
}
function drawKde(d) {
  if (!d.x.length) { document.getElementById('chart-dist').innerHTML = '<div class="muted">no values</div>'; return; }
  const bw = d.bandwidth, fmtbw = v => v == null ? '' : ` bw ${v.toPrecision(3)}`;
  document.getElementById('dist-bw').placeholder = bw.filtered != null ? `auto (${bw.filtered.toPrecision(3)})` : 'auto';
  const line = (y, name, color, fill) => ({ x: d.x, y, name, type: 'scatter', mode: 'lines', line: { color, width: 2 }, fill: 'tozeroy', fillcolor: fill, hovertemplate: `${label(d.col)} %{x:.4g}<br>${name}: %{y:.3g}<extra></extra>` });
  Plotly.react('chart-dist', [line(d.all, `all (n = ${fmt(d.n_all)}${fmtbw(bw.all)})`, '#b8bcc6', 'rgba(184,188,198,0.35)'), line(d.filtered, `filtered (n = ${fmt(d.n)}${fmtbw(bw.filtered)})`, '#2d5be3', 'rgba(45,91,227,0.35)')],
    { ...LAY, xaxis: { title: label(d.col), type: d.log ? 'log' : 'linear' }, yaxis: { title: d.log ? 'density (per log10 unit)' : 'density', type: document.getElementById('dist-ylog').checked ? 'log' : 'linear' }, legend: { x: 1, xanchor: 'right', y: 1 } }, PL);
}
function drawDist(d) {
  if (!d) return;
  if (d.kind === 'kde') return drawKde(d);
  if (!d.edges.length) { document.getElementById('chart-dist').innerHTML = '<div class="muted">no values</div>'; return; }
  const dens = document.getElementById('dist-y').value === 'density', ds = densities(d);
  const yf = dens ? ds.filtered : d.counts, ya = dens ? ds.all : d.all;
  const x = centers(d.edges), hover = d.edges.slice(0, -1).map((v, i) => `${v.toPrecision(4)} to ${d.edges[i + 1].toPrecision(4)}`);
  const step = (y, name, color, fill) => ({ x: d.edges, y: y.concat([y[y.length - 1]]), name, type: 'scatter', mode: 'lines', line: { shape: 'hv', color, width: 1.5 }, fill: 'tozeroy', fillcolor: fill, hoverinfo: 'skip' });
  const pts = (y, n, name, color) => ({ x, y, name, type: 'scatter', mode: 'markers', marker: { color, size: 1, opacity: 0 }, text: hover.map((h, i) => `${h}<br>${name}: ${fmt(n[i])} dsRNAs` + (dens ? ` (density ${y[i].toPrecision(3)})` : '')), hovertemplate: '%{text}<extra></extra>', showlegend: false });
  Plotly.react('chart-dist', [step(ya, `all (n = ${fmt(d.n_all)})`, '#b8bcc6', 'rgba(184,188,198,0.35)'), step(yf, `filtered (n = ${fmt(d.n)})`, '#2d5be3', 'rgba(45,91,227,0.35)'), pts(ya, d.all, 'all', '#b8bcc6'), pts(yf, d.counts, 'filtered', '#2d5be3')],
    { ...LAY, xaxis: { title: label(d.col), type: d.log ? 'log' : 'linear' }, yaxis: { title: dens ? (d.log ? 'density (per log10 unit)' : 'density') : 'dsRNAs', type: document.getElementById('dist-ylog').checked ? 'log' : 'linear' }, legend: { x: 1, xanchor: 'right', y: 1 } }, PL);
}
function drawScatter(d) {
  const el = document.getElementById('chart-scatter');
  if (d.mode === 'density') {
    document.getElementById('sc-note').textContent = `${fmt(d.n)} dsRNAs: shown as a density grid (more than the point limit; raise it or narrow the filters for points and box selection)`;
    Plotly.react(el, [{ type: 'heatmap', x: centers(d.xedges), y: centers(d.yedges), z: d.z.map(r => r.map(v => v || null)), ...palette('Blues', false), hovertemplate: `${label(d.x)} %{x:.4g}<br>${label(d.y)} %{y:.4g}<br>%{z:,} dsRNAs<extra></extra>`, colorbar: { title: 'dsRNAs', thickness: 12 } }],
      { ...LAY, xaxis: { title: label(d.x) }, yaxis: { title: label(d.y) } }, PL);
    return;
  }
  document.getElementById('sc-note').textContent = `${fmt(d.n)} dsRNAs; box-select points to turn the box into a filter`;
  const traces = [];
  const groups = d.cv ? [...new Set(d.cv)].sort() : [null];
  const colors = ['#2d5be3', '#e3572d', '#2da84a', '#8e44ad', '#f39c12'];
  groups.forEach((g, k) => { const idx = d.cv ? d.cv.map((v, i) => v === g ? i : -1).filter(i => i >= 0) : d.xv.map((_, i) => i);
    traces.push({ type: 'scattergl', mode: 'markers', name: g == null ? 'dsRNAs' : `${label(d.color)} ${g}`, x: idx.map(i => d.xv[i]), y: idx.map(i => d.yv[i]), customdata: idx.map(i => d.id[i]),
      marker: { size: d.n > 5000 ? 3 : 5, opacity: d.n > 5000 ? 0.4 : 0.7, color: colors[k % colors.length] }, hovertemplate: `${label(d.x)} %{x:.4g}<br>${label(d.y)} %{y:.4g}<br>id %{customdata}<extra></extra>` }); });
  Plotly.react(el, traces, { ...LAY, xaxis: { title: label(d.x) }, yaxis: { title: label(d.y) }, dragmode: 'select', showlegend: groups.length > 1, legend: { x: 1, xanchor: 'right', y: 1 } }, PL);
  el.removeAllListeners && el.removeAllListeners('plotly_selected');
  el.on('plotly_selected', ev => { if (!ev || !ev.range) return;
    const set = (kind, col, v) => { const i = document.querySelector(`#ranges input[data-kind=${kind}][data-col="${col}"]`); if (i) i.value = +v.toPrecision(5); };
    set('min', d.x, ev.range.x[0]); set('max', d.x, ev.range.x[1]); set('min', d.y, ev.range.y[0]); set('max', d.y, ev.range.y[1]); readFilters(); load(); });
  el.on('plotly_click', ev => { const id = ev.points[0].customdata; if (id != null) select(id, null); });
}
// light-to-dark palettes (ColorBrewer), so larger values are darker; Viridis/Plasma are Plotly's
const PALETTES = {
  Blues: [[0, '#f7fbff'], [0.25, '#c6dbef'], [0.5, '#6baed6'], [0.75, '#2171b5'], [1, '#08306b']],
  Greens: [[0, '#f7fcf5'], [0.25, '#c7e9c0'], [0.5, '#74c476'], [0.75, '#238b45'], [1, '#00441b']],
  Reds: [[0, '#fff5f0'], [0.25, '#fcbba1'], [0.5, '#fb6a4a'], [0.75, '#cb181d'], [1, '#67000d']],
  Oranges: [[0, '#fff5eb'], [0.25, '#fdd0a2'], [0.5, '#fd8d3c'], [0.75, '#d94801'], [1, '#7f2704']],
  Greys: [[0, '#ffffff'], [0.25, '#d9d9d9'], [0.5, '#969696'], [0.75, '#525252'], [1, '#000000']],
  YlGnBu: [[0, '#ffffd9'], [0.25, '#c7e9b4'], [0.5, '#41b6c4'], [0.75, '#225ea8'], [1, '#081d58']],
  Viridis: 'Viridis', Plasma: 'Plasma',
};
function palette(name, reverse) { const p = PALETTES[name] || PALETTES.Blues; if (typeof p === 'string') return { colorscale: p, reversescale: reverse }; return { colorscale: reverse ? p.map(([v, c], i, a) => [v, a[a.length - 1 - i][1]]) : p }; }
function drawGenome(d) {
  const el = document.getElementById('chart-genome');
  if (!d) return;
  if (!d.chromosomes.length) { el.innerHTML = '<div class="muted">nothing to show</div>'; return; }
  const names = d.chromosomes.map(c => c.name), nb = Math.max(...d.chromosomes.map(c => c.counts.length));
  const counts = d.chromosomes.map(c => c.counts.concat(Array(nb - c.counts.length).fill(null)));
  const mb = d.bin / 1e6, x = Array.from({ length: nb }, (_, i) => (i + 0.5) * mb);
  const logz = document.getElementById('gen-zscale').value === 'log', zmaxIn = parseFloat(document.getElementById('gen-zmax').value);
  const dataMax = Math.max(1, ...counts.flat().filter(v => v != null)), zmax = zmaxIn > 0 ? zmaxIn : dataMax;
  const tf = v => v == null ? null : logz ? Math.log10(1 + v) : v;
  const z = counts.map(r => r.map(tf));
  const ticks = logz ? [0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000].filter(v => v <= zmax) : undefined;
  const colorbar = { title: 'dsRNAs', thickness: 12, ...(logz ? { tickvals: ticks.map(tf), ticktext: ticks.map(v => fmt(v)) } : {}) };
  document.getElementById('gen-note').textContent = `dsRNAs per ${mb >= 1 ? mb + ' Mb' : d.bin / 1e3 + ' kb'} bin` + (d.omitted ? `; ${d.omitted} short sequences not shown` : '') + '; click a bin to filter to it';
  Plotly.react(el, [{ type: 'heatmap', x, y: names, z, customdata: counts, zmin: 0, zmax: tf(zmax), ...palette(document.getElementById('gen-cs').value, document.getElementById('gen-rev').checked),
      hoverongaps: false, hovertemplate: '%{y}: %{x:.2f} Mb<br>%{customdata:,} dsRNAs<extra></extra>', colorbar, xgap: 1, ygap: 1 }],
    { ...LAY, height: Math.max(260, 22 * names.length + 80), xaxis: { title: 'position (Mb)' }, yaxis: { autorange: 'reversed', type: 'category', dtick: 1 } }, PL);
  el.removeAllListeners && el.removeAllListeners('plotly_click');
  el.on('plotly_click', ev => { const p = ev.points[0]; const start = Math.round((p.x - mb / 2) * 1e6) + 1;
    document.getElementById('search').value = `${p.y}:${start}-${start + d.bin - 1}`; readFilters(); load(); });
}
function drawGroups(d) {
  if (!d) return;
  const sel = document.getElementById('gr-metric'); const cur = sel.value;
  sel.innerHTML = d.columns.slice(1).map(c => `<option ${c === cur ? 'selected' : ''}>${c}</option>`).join('');
  const m = d.columns.indexOf(sel.value) > 0 ? d.columns.indexOf(sel.value) : 1;
  Plotly.react('chart-groups', [{ type: 'bar', x: d.rows.map(r => String(r[0])), y: d.rows.map(r => r[m]), marker: { color: '#2d5be3' }, hovertemplate: '%{x}: %{y:,.3~f}<extra></extra>' }],
    { ...LAY, margin: { l: 60, r: 20, t: 10, b: 70 }, xaxis: { title: label(d.by), type: 'category' }, yaxis: { title: d.columns[m] } }, PL);
  const f = v => v == null ? '' : typeof v === 'number' ? (Number.isInteger(v) ? fmt(v) : v.toFixed(2)) : v;
  document.getElementById('table-groups').innerHTML = '<table class="stat"><thead><tr>' + d.columns.map(c => `<th>${c}</th>`).join('') + '</tr></thead><tbody>' +
    d.rows.map(r => '<tr>' + r.map((v, i) => `<td class="${i ? 'num' : ''}">${f(v)}</td>`).join('') + '</tr>').join('') + '</tbody></table>';
}
function drawStats(d) {
  const f = v => v == null ? '' : Number.isInteger(v) ? fmt(v) : v.toFixed(2);
  const keys = ['n', 'mean', 'median', 'q1', 'q3', 'min', 'max'];
  document.getElementById('table-stats').innerHTML = '<table class="stat"><thead><tr><th>column</th><th></th>' + keys.map(k => `<th class="num">${k}</th>`).join('') + '</tr></thead><tbody>' +
    d.columns.map(c => `<tr><td rowspan="2"><b>${label(c)}</b></td><td>filtered</td>` + keys.map(k => `<td class="num">${f(d.filtered[c][k])}</td>`).join('') + '</tr>' +
      `<tr><td style="color:#999">all</td>` + keys.map(k => `<td class="num" style="color:#999">${f(d.all[c][k])}</td>`).join('') + '</tr>').join('') + '</tbody></table>';
}
// ---------------------------------------------------------------- IGV tab (igv.js)
const G = { browser: null, loading: false, tracks: [], locus: null };
function parseBP(text) {  // igv.js 'bp' features, as the site's reader builds them
  const feats = [], colors = [];
  for (const raw of text.split('\n')) {
    const line = raw.trim(); if (!line) continue;
    if (line.startsWith('color:')) { const t = line.split(/\s+/); if (t.length >= 4) colors.push(`rgb(${t[1]},${t[2]},${t[3]})`); continue; }
    const t = line.split('\t'); if (t.length < 6) continue;
    const a0 = +t[1] - 1, a1 = +t[2] - 1, b0 = +t[3], b1 = +t[4], ci = +t[5];
    const f = { chr: t[0], startLeft: Math.min(a0, a1), startRight: Math.max(a0, a1), endLeft: Math.min(b0, b1), endRight: Math.max(b0, b1), color: colors[ci] || 'rgb(0,0,128)', score: ci };
    f.start = f.startLeft; f.end = f.endRight; feats.push(f);
  }
  return feats;
}
const bpSource = strand => ({ getFeatures: async ({ chr, start, end }) => {
  if (!chr || chr.toLowerCase() === 'all') return [];
  const r = await fetch('/api/igv/bp?' + params({ chr, start: Math.max(0, Math.floor(start)), end: Math.ceil(end), strand })); return parseBP(await r.text()); } });
const editsSource = () => ({ getFeatures: async ({ chr, start, end }) => {
  if (!chr || chr.toLowerCase() === 'all') return [];
  return (await fetch('/api/igv/edits?' + params({ chr, start: Math.max(0, Math.floor(start)), end: Math.ceil(end) }))).json(); } });
function locusOf(chr, s, e) { const pad = Math.max(200, Math.round((e - s) * 0.2)); return `${chr}:${Math.max(1, s - pad)}-${e + pad}`; }
async function showIgv() { if (G.browser) await refreshIgv(); else if (!G.loading) await initIgv(); }
async function initIgv() {
  const el = document.getElementById('igv-div'), note = document.getElementById('igv-note');
  G.loading = true;
  try {
    const cfg = await api('/api/igv/config');
    if (!cfg.reference) { note.textContent = 'no genome'; el.innerHTML = `<div class="muted" style="text-align:left;max-width:640px">${cfg.message}</div>`; return; }
    note.textContent = 'loading igv.js ...';
    const igvjs = (await import('https://cdn.jsdelivr.net/npm/igv@3.8.9/dist/igv.esm.min.js')).default;
    const fr = T.firstRow;
    G.locus = G.locus || (fr ? locusOf(fr.Chromosome, fr.i_start, fr.j_end) : undefined);
    G.browser = await igvjs.createBrowser(el, { loadDefaultGenomes: false, reference: cfg.reference, locus: G.locus, queryParametersSupported: false, showAllChromosomes: false, tracks: [] });
    for (const t of cfg.reference.tracks || []) { try { await G.browser.loadTrack(t); } catch (e) { console.warn('gene track', e); } }
    for (const [name, strand, orient, order] of [['dsRNAs + strand', '+', 'UP', -1], ['dsRNAs - strand', '-', 'DOWN', 100]]) {
      const t = await G.browser.loadTrack({ name, format: 'bp', url: '/api/igv/bp?stub=1', arcOrientation: orient, order, height: 140 });
      t.featureSource = bpSource(strand); G.tracks.push(t);
    }
    if (S.summary.has_editing) {
      const t = await G.browser.loadTrack({ name: 'editing sites', type: 'annotation', format: 'bed', url: '/api/igv/edits?stub=1', height: 40, displayMode: 'COLLAPSED', color: '#27ae60' });
      t.featureSource = editsSource(); G.tracks.push(t);
    }
    await refreshIgv();
    note.textContent = `${cfg.source}; arcs follow the current filters (up to 20,000 per view)`;
    // the results' chromosome names must exist in the reference, or every track stays empty
    const known = new Set((G.browser.genome?.chromosomeNames || []).map(c => c.toLowerCase()));
    const mine = S.summary.chromosomes.map(c => c.name);
    if (known.size && !mine.some(c => known.has(c.toLowerCase())))
      note.innerHTML += ` <span class="err" style="display:inline">none of the results' sequence names (${mine.slice(0, 3).join(', ')}${mine.length > 3 ? ', ...' : ''}) match the reference's (${[...known].slice(0, 3).join(', ')}, ...): use --genome with the FASTA the scan used</span>`;
  } catch (e) { note.textContent = ''; el.innerHTML = `<div class="err">IGV: ${e.message}</div>`; }
  finally { G.loading = false; }
}
async function refreshIgv() { for (const t of G.tracks) t.clearCachedFeatures(); await Promise.all(G.tracks.map(t => t.trackView?.updateViews())); }
function igvTab(chr, s, e) { G.locus = locusOf(chr, s, e); document.querySelector('#tabs button[data-tab=igv]').click(); if (G.browser) G.browser.search(G.locus); }
function igvFit() { if (G.browser && G.locus) G.browser.search(G.locus); }
function dlImage(tab, format) { const el = document.getElementById('chart-' + tab); if (el && el.data) Plotly.downloadImage(el, { format, width: 1200, height: tab === 'genome' ? Math.max(500, el.clientHeight) : 700, filename: 'dsrnascan_' + tab }); }
function dlData(tab) {
  const d = T.data[tab]; if (!d) return;
  let rows = [];
  if (tab === 'dist' && d.kind === 'kde') { rows = [[d.col, 'filtered_density', 'all_density']].concat(d.x.map((x, i) => [x, d.filtered[i], d.all[i]])); }
  else if (tab === 'dist') { const ds = densities(d); rows = [['bin_start', 'bin_end', 'filtered_count', 'all_count', 'filtered_density', 'all_density']].concat(d.counts.map((c, i) => [d.edges[i], d.edges[i + 1], c, d.all[i], ds.filtered[i], ds.all[i]])); }
  else if (tab === 'scatter') { rows = d.mode === 'points' ? [['id', d.x, d.y].concat(d.color ? [d.color] : [])].concat(d.id.map((id, i) => [id, d.xv[i], d.yv[i]].concat(d.color ? [d.cv[i]] : []))) : [[d.x + '_bin_start', d.y + '_bin_start', 'count']].concat(d.z.flatMap((r, yi) => r.map((v, xi) => [d.xedges[xi], d.yedges[yi], v]).filter(t => t[2]))); }
  else if (tab === 'genome') { rows = [['chromosome', 'bin_start', 'bin_end', 'count']].concat(d.chromosomes.flatMap(c => c.counts.map((v, i) => [c.name, i * d.bin + 1, Math.min((i + 1) * d.bin, c.length), v]))); }
  else if (tab === 'groups') { rows = [d.columns].concat(d.rows); }
  else if (tab === 'stats') { rows = [['column', 'set', 'n', 'mean', 'median', 'q1', 'q3', 'min', 'max']].concat(d.columns.flatMap(c => ['filtered', 'all'].map(k => [c, k].concat(['n', 'mean', 'median', 'q1', 'q3', 'min', 'max'].map(s => d[k][c][s]))))); }
  const tsv = rows.map(r => r.map(v => v == null ? '' : v).join('\t')).join('\n') + '\n';
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([tsv], { type: 'text/tab-separated-values' })); a.download = 'dsrnascan_' + tab + '.tsv'; a.click(); URL.revokeObjectURL(a.href);
}
async function select(id, tr) {
  document.querySelectorAll('tr.sel').forEach(e => e.classList.remove('sel'));
  if (T.active === 'igv' && G.browser && tr) { const r = T.rows && T.rows.find(x => x.id === id); if (r) { G.locus = locusOf(r.Chromosome, r.i_start, r.j_end); G.browser.search(G.locus); } }
  if (tr) tr.classList.add('sel');
  S.selected = id;
  try { const d = await api('/api/result/' + id); showMetrics(d); showStructure(d); }
  catch (e) { document.getElementById('metrics').innerHTML = `<div class="err">${e.message}</div>`; }
}
function showMetrics(d) {
  const loc = `${d.Chromosome}:${fmt(d.i_start)}-${fmt(d.j_end)}`;
  const tag = (v) => v == null ? '' : `<span class="tag ${v === 'Yes' ? 'yes' : 'no'}">${v === 'Yes' ? 'yes' : 'no'}</span>`;
  const btn = 'style="font-size:11px;padding:2px 7px;margin:0 2px 2px 0"';
  let h = `<div><span class="lab">dsRNA_${d.Chromosome}_${String(d.id + 1).padStart(6, '0')}</span> &nbsp; <b>${loc}</b> (${d.Strand}) &nbsp;` +
    `<button class="btn" ${btn} onclick="igvTab('${d.Chromosome}', ${d.i_start}, ${d.j_end})">IGV tab</button><button class="btn" ${btn} onclick="igv('${loc}')">IGV app</button><button class="btn" ${btn} onclick="copy('${loc}', this)">copy locus</button></div>` +
    `<div><span class="lab">arms</span> <b>${d.i_start}-${d.i_end}</b> and <b>${d.j_start}-${d.j_end}</b> (${d.i_length}/${d.j_length} nt, loop ${fmt(d.loop_length)} nt)</div>` +
    `<div><span class="lab">energy</span> <b>${num(d['dG(kcal/mol)'], 1)} kcal/mol</b> &nbsp; <span class="lab">paired</span> <b>${d.base_pairs} bp (${num(d.percent_paired, 1)}%)</b> &nbsp; <span class="lab">longest helix</span> <b>${d.longest_helix}</b> &nbsp; <span class="lab">search score</span> <b>${d.Score}</b></div>`;
  if (d.stability_model_score != null) h += `<div><span class="lab">stability model</span> <b>${num(d.stability_model_score, 3)}</b> ${tag(d.likely_edited)} &nbsp; <span class="lab">probing model</span> <b>${num(d.probing_model_score, 3)}</b> ${tag(d.likely_forms)}</div>`;
  if (d.edits != null) h += `<div><span class="lab">editing sites</span> <b>${d.edits_i}</b> in i arm, <b>${d.edits_j}</b> in j arm</div>`;
  const parts = (d.structure || '&').split('&');
  const dbn = `>dsRNA_${d.Chromosome}_${String(d.id + 1).padStart(6, '0')} ${loc}(${d.Strand}) dG=${d['dG(kcal/mol)']}\n${d.i_seq}&${d.j_seq}\n${d.structure}`;
  h += `<div style="margin-top:6px"><span class="lab">copy</span> <button class="btn" ${btn} onclick="copyText(this)" data-t="${enc(dbn)}">dot-bracket</button>` +
    `<button class="btn" ${btn} onclick="copyText(this)" data-t="${enc(d.i_seq)}">i seq</button><button class="btn" ${btn} onclick="copyText(this)" data-t="${enc(d.j_seq)}">j seq</button>` +
    `<button class="btn" ${btn} onclick="copyText(this)" data-t="${enc(parts[0])}">i structure</button><button class="btn" ${btn} onclick="copyText(this)" data-t="${enc(parts[1] || '')}">j structure</button></div>`;
  h += '<div style="margin-top:6px;font-size:10px;color:#888;line-height:1.4">Stability model: likelihood of ADAR editing from structure (cutoff 0.553). Probing model: agreement with in vivo structure probing, 0-1 (cutoff 0.440). Both trained on human data.</div>';
  document.getElementById('metrics').innerHTML = h;
}
function enc(s) { return encodeURIComponent(s || ''); }
function copyText(btn) { copy(decodeURIComponent(btn.dataset.t), btn); }
function copy(text, btn) { navigator.clipboard.writeText(text).then(() => { const o = btn.textContent; btn.textContent = 'copied'; setTimeout(() => btn.textContent = o, 1200); }); }
function igv(loc) { fetch('http://localhost:60151/goto?locus=' + loc, { mode: 'no-cors' }).catch(() => {}); }
function showStructure(d) {
  const v = document.getElementById('forna');
  v.innerHTML = '';
  const parts = (d.structure || '').split('&');
  const seq = (d.i_seq || '') + (d.j_seq || ''), str = (parts[0] || '') + (parts[1] || '');
  const maxArm = Math.max((d.i_seq || '').length, (d.j_seq || '').length);
  if (!seq || maxArm > 500) { v.innerHTML = `<div class="muted">${maxArm > 500 ? 'Arms too long for the interactive view (' + maxArm + ' nt); use the copy buttons.' : 'No structure'}</div>`; return; }
  const c = document.createElement('div'); c.id = 'forna-c'; c.style.width = '100%'; c.style.height = '100%'; v.appendChild(c);
  const rna = new fornac.FornaContainer('#forna-c', { applyForce: true, allowPanningAndZooming: true, initialSize: [v.offsetWidth || 450, v.offsetHeight || 400], friction: 0.35, middleCharge: -30, otherCharge: -30 });
  rna.addRNA(str, { structure: str, sequence: seq, labelInterval: 10 });
  if (d.editing_sites && d.editing_sites.length) setTimeout(() => annotate(d.editing_sites), 1500);
  setTimeout(() => { rna.centerView(); rna.zoomToFit(); }, 1000);
}
function annotate(sites) {
  const circles = d3.select('#forna-c svg').selectAll('circle');
  const counts = {}; circles.each(function () { const r = d3.select(this).attr('r'); counts[r] = (counts[r] || 0) + 1; });
  const common = Object.entries(counts).sort((a, b) => b[1] - a[1])[0]?.[0];
  const nucs = []; circles.each(function () { if (d3.select(this).attr('r') === common) nucs.push(this); });
  for (const [pos, f] of sites) if (pos < nucs.length) {
    const c = d3.select(nucs[pos]); c.style('stroke', f >= 0.8 ? '#0B5345' : f >= 0.5 ? '#148F77' : f >= 0.3 ? '#27AE60' : '#52BE80').style('stroke-width', '3px');
    let t = c.select('title'); if (t.empty()) t = c.append('title'); t.text(`editing site, position ${pos + 1} (${(f * 100).toFixed(0)}%)`);
  }
}
document.getElementById('apply').addEventListener('click', () => { readFilters(); load(); });
document.getElementById('reset').addEventListener('click', () => {
  document.getElementById('search').value = ''; document.getElementById('where').value = ''; document.getElementById('strand').value = '';
  document.getElementById('chrom').selectedIndex = -1;
  document.querySelectorAll('#flags input').forEach(e => e.checked = false); document.querySelectorAll('#ranges input').forEach(e => e.value = '');
  readFilters(); load(); });
document.querySelector('.side').addEventListener('keydown', e => { if (e.key === 'Enter' && e.target.tagName !== 'TEXTAREA') { readFilters(); load(); } });
document.getElementById('first').addEventListener('click', () => { S.page = 0; load(); });
document.getElementById('prev').addEventListener('click', () => { S.page = Math.max(0, S.page - 1); load(); });
document.getElementById('next').addEventListener('click', () => { S.page++; load(); });
document.getElementById('last').addEventListener('click', () => { S.page = Math.max(0, Math.ceil(S.total / S.size) - 1); load(); });
document.getElementById('pagesize').addEventListener('change', e => { S.size = +e.target.value; S.page = 0; load(); });
document.getElementById('export').addEventListener('click', () => { window.location.href = '/api/export?' + params({ format: document.getElementById('fmt').value, sort: S.sort, dir: S.dir }); });
(async () => {
  S.summary = await api('/api/summary');
  document.getElementById('sub').textContent = `${S.summary.files.join(', ')}: ${fmt(S.summary.total)} dsRNAs on ${S.summary.chromosomes.length} sequences` + (S.summary.has_editing ? ` \u00b7 ${fmt(S.summary.editing_sites)} editing sites loaded` : '');
  buildSidebar();
  buildTabs();
  load();
})();
</script>
</body>
</html>'''


def main():
    p = argparse.ArgumentParser(description='Browse, filter and export dsRNAscan results',
                                epilog='Export filters: --chrom chr19 --region chr19:1-3000000 --strand + '
                                       '--min percent_paired 80 --max "dG(kcal/mol)" -100 --likely-edited --likely-forms '
                                       '--where "longest_helix > 30"')
    p.add_argument('path', nargs='?', default='.', help='dsRNAscan output folder, or a results file (default: .)')
    p.add_argument('--editing-file', help='BED or GFF3 of RNA editing sites: counted per arm, shown on structures')
    p.add_argument('--large-editing-file', action='store_true', help=argparse.SUPPRESS)  # older option, no longer needed
    p.add_argument('--port', type=int, default=8080, help='web server port (default: 8080)')
    p.add_argument('--genome', help='for the IGV tab: an igv.js assembly id (hg38, mm39, dm6, ce11, ...) or the uncompressed '
                                    'FASTA the scan used (default: the results name if it is a hosted assembly, else the FASTA from dsrnascan.log)')
    p.add_argument('--no-browser', action='store_true', help='do not open the web browser')
    x = p.add_argument_group('export instead of serving')
    x.add_argument('--export', metavar='FILE', help='write the filtered dsRNAs to FILE and exit; the format follows the '
                                                    'extension (' + ', '.join(f'.{v.lstrip(".")}' for v in EXPORT_FORMATS.values()) + ')')
    x.add_argument('--format', choices=sorted(EXPORT_FORMATS), help='export format, if not clear from the extension')
    x.add_argument('--count', action='store_true', help='print the number of matching dsRNAs and exit')
    x.add_argument('--columns', action='store_true', help='list the filterable columns and exit')
    x.add_argument('--chrom', action='append', default=[], help='keep this chromosome (repeatable)')
    x.add_argument('--region', help='keep dsRNAs overlapping chr:start-end')
    x.add_argument('--strand', choices=['+', '-'])
    x.add_argument('--min', nargs=2, action='append', default=[], metavar=('COLUMN', 'VALUE'), help='keep COLUMN >= VALUE (repeatable)')
    x.add_argument('--max', nargs=2, action='append', default=[], metavar=('COLUMN', 'VALUE'), help='keep COLUMN <= VALUE (repeatable)')
    x.add_argument('--likely-edited', action='store_true', help='keep only likely_edited = Yes')
    x.add_argument('--likely-forms', action='store_true', help='keep only likely_forms = Yes')
    x.add_argument('--where', help='extra SQL condition on the columns, e.g. "i_length > 300 AND edits > 0"')
    x.add_argument('--sort', default='id', help='sort column for the export (default: file order)')
    x.add_argument('--desc', action='store_true', help='sort descending')
    x.add_argument('--limit', type=int, help='write at most this many dsRNAs')
    args = p.parse_args()

    try:
        results = Results(args.path, args.editing_file)
    except FileNotFoundError as e:
        sys.exit(f"Error: {e}")
    if args.columns:
        for c, t in results.columns:
            print(f"{c:24s} {t}")
        return
    f = {'chrom': args.chrom, 'strand': args.strand, 'where': args.where, 'likely_edited': args.likely_edited,
         'likely_forms': args.likely_forms, 'min': dict(args.min), 'max': dict(args.max)}
    if args.region:
        f['region'] = Results.filters_from_params({'region': [args.region]})['region']
    if args.count:
        print(results.count(f))
        return
    if args.export:
        n = results.export(args.export, args.format, f, args.sort, args.desc, args.limit)
        print(f"{n:,} dsRNAs written to {args.export}")
        return

    Handler.results = results
    Handler.html = HTML
    Handler.stem = os.path.basename(results.files[0]).split('.dsrnas')[0].split('_merged_results')[0]
    Handler.genome_opt, Handler.path = args.genome, args.path
    print(f"{results.summary()['total']:,} dsRNAs in {', '.join(os.path.basename(x) for x in results.files)}")
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.ThreadingTCPServer(("", args.port), Handler) as httpd:
        print(f"Browser: http://localhost:{args.port}/   (Ctrl+C to stop)")
        if not args.no_browser:
            webbrowser.open(f'http://localhost:{args.port}/')
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nStopped.")


if __name__ == '__main__':
    main()
