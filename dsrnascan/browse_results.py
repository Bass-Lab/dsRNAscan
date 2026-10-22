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

from dsrnascan.results import EXPORT_FORMATS, Results


class Handler(SimpleHTTPRequestHandler):
    results = None   # set before the server starts
    html = ""
    stem = "dsrnascan"

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
                self._json(self.results.histogram(params['col'][0], f, int(params.get('bins', ['30'])[0])))
            elif url.path == '/api/export':
                self._export(params)
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


HTML = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>dsRNA Browser</title>
<script src="https://d3js.org/d3.v3.min.js"></script>
<script src="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.js"></script>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.css">
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
    <div class="tablewrap"><table><thead><tr id="head"></tr></thead><tbody id="body"></tbody></table></div>
    <div class="pager">
      <button class="btn" id="first">&laquo;</button><button class="btn" id="prev">&lsaquo;</button>
      <span id="pageinfo"></span>
      <button class="btn" id="next">&rsaquo;</button><button class="btn" id="last">&raquo;</button>
      <select id="pagesize"><option>25</option><option selected>50</option><option>100</option><option>200</option></select>
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
  } catch (e) { document.getElementById('err').textContent = e.message; }
}
function render(rows) {
  const body = document.getElementById('body');
  body.innerHTML = '';
  const score = (v, t) => v == null ? '<td class="num" style="color:#bbb">-</td>' : `<td class="num" style="${v >= t ? 'color:#1e7b45;font-weight:600' : 'color:#999'}">${num(v, 3)}</td>`;
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
  if (!el) return;
  const w = 150, h = 42, n = hd.counts.length, max = Math.max(1, ...hd.counts), bw = w / n;
  el.innerHTML = hd.counts.map((c, i) => { const bh = c / max * (h - 2); const lo = hd.lo + (hd.hi - hd.lo) * i / n, hi = lo + (hd.hi - hd.lo) / n;
    return `<rect x="${i * bw}" y="${h - bh}" width="${Math.max(bw - 1, 1)}" height="${bh}"><title>${lo.toFixed(1)} to ${hi.toFixed(1)}: ${fmt(c)}</title></rect>`; }).join('');
}
async function select(id, tr) {
  document.querySelectorAll('tr.sel').forEach(e => e.classList.remove('sel'));
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
    `<button class="btn" ${btn} onclick="igv('${loc}')">IGV</button><button class="btn" ${btn} onclick="copy('${loc}', this)">copy locus</button></div>` +
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
