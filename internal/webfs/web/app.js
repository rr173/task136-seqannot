// seqannot frontend — native, no build. All computation is on the Go side;
// this script only issues HTTP calls and renders results.
'use strict';

let currentSeq = '';

const $ = (id) => document.getElementById(id);
const show = (el, text) => { el.textContent = text; };
const json = (v) => JSON.stringify(v, null, 2);

async function api(method, path, body) {
  const opt = { method, headers: {} };
  if (body !== undefined) {
    opt.headers['Content-Type'] = 'application/json';
    opt.body = JSON.stringify(body);
  }
  const resp = await fetch(path, opt);
  const text = await resp.text();
  let data;
  try { data = text ? JSON.parse(text) : null; } catch { data = text; }
  return { status: resp.status, data };
}

function setCur(id) {
  currentSeq = id;
  show($('cur-seq'), id || '（未选择）');
}

async function submitSeq() {
  const name = $('seq-name').value.trim() || 'seq1';
  const residues = $('seq-residues').value.trim();
  const type = $('seq-type').value;
  if (!residues) { show($('seq-status'), '请输入序列'); return; }
  // if it starts with '>', treat as FASTA
  if (residues.startsWith('>')) {
    const r = await api('POST', '/sequences/import-fasta', { fasta: residues });
    if (r.status >= 300) { show($('seq-status'), '导入失败: ' + json(r.data)); return; }
    setCur(r.data.id);
    show($('seq-status'), 'FASTA 导入成功: ' + r.data.id + ' (len=' + r.data.length + ')');
    return;
  }
  const r = await api('POST', '/sequences', { name, residues, type });
  if (r.status >= 300) { show($('seq-status'), '提交失败: ' + json(r.data)); return; }
  setCur(r.data.id);
  show($('seq-status'), '提交成功: ' + r.data.id + ' (len=' + r.data.length + ', gc=' + (r.data.length ? '' : '') + ')');
}

async function importFasta() {
  const fasta = $('seq-residues').value.trim();
  if (!fasta) { show($('seq-status'), '请粘贴 FASTA'); return; }
  const r = await api('POST', '/sequences/import-fasta', { fasta });
  if (r.status >= 300) { show($('seq-status'), '导入失败: ' + json(r.data)); return; }
  setCur(r.data.id);
  show($('seq-status'), 'FASTA 导入: ' + r.data.id);
}

async function runAnalysis(kind) {
  if (!currentSeq) { show($('analyze-out'), '请先提交序列'); return; }
  let path = '/sequences/' + currentSeq + '/' + kind;
  let body = undefined;
  if (kind === 'orfs') path += '?min_aa_len=' + (encodeURIComponent($('min-aa').value || '1'));
  if (kind === 'motif') body = { pattern: $('motif-pattern').value };
  if (kind === 'restriction' || kind === 'digest') {
    body = { enzyme_ids: ($('enzyme-ids').value || '').split(',').map(s=>s.trim()).filter(Boolean) };
  }
  const r = await api('POST', path, body);
  show($('analyze-out'), kind + ' -> ' + r.status + '\n' + json(r.data));
}

async function submitJob() {
  if (!currentSeq) { show($('analyze-out'), '请先提交序列'); return; }
  const body = {
    sequence_id: currentSeq,
    steps: [
      { name: 'composition' },
      { name: 'translate' },
      { name: 'orf', params: { min_aa_len: parseInt($('min-aa').value || '1', 10) } },
    ],
  };
  const r = await api('POST', '/jobs', body);
  show($('analyze-out'), 'job -> ' + r.status + '\n' + json(r.data));
}

async function addMotif() {
  const name = $('motif-name').value.trim();
  const pattern = $('motif-pat').value.trim();
  if (!name || !pattern) { show($('library-out'), '需要 name 和 pattern'); return; }
  const r = await api('POST', '/motifs', { name, pattern });
  show($('library-out'), 'add motif -> ' + r.status + '\n' + json(r.data));
}

async function addEnzyme() {
  const name = $('enz-name').value.trim();
  const site = $('enz-site').value.trim();
  const cut = parseInt($('enz-cut').value || '0', 10);
  if (!name || !site) { show($('library-out'), '需要 name 和 site'); return; }
  const r = await api('POST', '/enzymes', { name, site, cut_offset: cut });
  show($('library-out'), 'add enzyme -> ' + r.status + '\n' + json(r.data));
}

async function listMotifs() {
  const r = await api('GET', '/motifs');
  show($('library-out'), 'motifs -> ' + r.status + '\n' + json(r.data));
}

async function listEnzymes() {
  const r = await api('GET', '/enzymes');
  show($('library-out'), 'enzymes -> ' + r.status + '\n' + json(r.data));
}

async function listJobs() {
  const r = await api('GET', '/jobs');
  show($('job-out'), 'jobs -> ' + r.status + '\n' + json(r.data));
}

async function reconcile() {
  const r = await api('POST', '/admin/reconcile');
  show($('job-out'), 'reconcile -> ' + r.status + '\n' + json(r.data));
}

document.addEventListener('DOMContentLoaded', () => {
  $('seq-submit').addEventListener('click', submitSeq);
  $('seq-import').addEventListener('click', importFasta);
  document.querySelectorAll('button[data-an]').forEach(b => {
    b.addEventListener('click', () => runAnalysis(b.dataset.an));
  });
  document.querySelector('button[data-job]').addEventListener('click', submitJob);
  $('motif-add').addEventListener('click', addMotif);
  $('enz-add').addEventListener('click', addEnzyme);
  $('motif-list').addEventListener('click', listMotifs);
  $('enz-list').addEventListener('click', listEnzymes);
  $('job-list').addEventListener('click', listJobs);
  $('reconcile').addEventListener('click', reconcile);
});
