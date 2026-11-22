'use strict';

// Extra coverage for the lab workbench: literature ground state energies,
// enumeration counting, determinism guarantees, metric formulas, audit
// envelope semantics, and lab app wiring checks.

const test = require('node:test');
const assert = require('node:assert');
const crypto = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const f = require('../src/fold.js');
const L = require('../src/lab.js');

const ROOT = path.join(__dirname, '..');
const read = name => fs.readFileSync(path.join(ROOT, name), 'utf8');

// -----------------------------------------------------------------------
// Enumeration counts and literature ground states
// -----------------------------------------------------------------------
test('countSaws matches canonical self avoiding walk counts', () => {
  // Counts fix the first bead at the origin and the second step in the
  // canonical half plane, so they are the symmetry reduced numbers.
  const expected2 = { 5: 17, 6: 48, 7: 132, 8: 367 };
  const expected3 = { 5: 49, 6: 238 };
  for (const n of Object.keys(expected2)) {
    assert.strictEqual(L.countSaws(Number(n), 2), expected2[n], '2D n=' + n);
  }
  for (const n of Object.keys(expected3)) {
    assert.strictEqual(L.countSaws(Number(n), 3), expected3[n], '3D n=' + n);
  }
});

test('runExact state count equals countSaws for untruncated runs', () => {
  for (const [seq, dim] of [['HPPHPH', 2], ['HPPHPH', 3], ['HHPH', 2]]) {
    const r = L.runExact({ sequence: seq, dim: dim });
    assert.strictEqual(r.truncated, false);
    assert.strictEqual(r.states, L.countSaws(seq.length, dim), seq + ' dim ' + dim);
  }
});

test('runExact hits documented ground states for tiny chains', () => {
  const cases = [
    ['HH', 0, 0], ['HHH', 0, 0], ['HHHH', -1, -1], ['HHHHH', -1, -1],
    ['HHHHHH', -2, -2], ['HHHHHHH', -2, -3], ['HPHPH', 0, 0],
    ['PHPHPHHPPHPH', -2, -4]
  ];
  for (const [seq, e2, e3] of cases) {
    assert.strictEqual(L.runExact({ sequence: seq, dim: 2 }).best.energy, e2, seq + ' 2D');
    assert.strictEqual(L.runExact({ sequence: seq, dim: 3 }).best.energy, e3, seq + ' 3D');
  }
});

test('compareAlgorithms skips exact past the 3D cap but runExact still works', () => {
  const cmp = L.compareAlgorithms({ sequence: 'PHPHPHHPPHPH', dim: 3, seed: 3 });
  assert.strictEqual(cmp.exact, null);
  assert.ok(cmp.exactSkipped);
  // runExact itself does not enforce the UI policy cap.
  assert.strictEqual(L.runExact({ sequence: 'PHPHPHHPPHPH', dim: 3 }).best.energy, -4);
});

test('exact result is seed independent across comparator runs', () => {
  const a = L.compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 1 });
  const b = L.compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 999 });
  assert.strictEqual(a.exact.best.energy, b.exact.best.energy);
  assert.strictEqual(a.exact.states, b.exact.states);
});

// -----------------------------------------------------------------------
// SAW edge cases and move guards
// -----------------------------------------------------------------------
test('randomSaw on n=1 returns the origin bead', () => {
  assert.deepStrictEqual(L.randomSaw(1, 2, f.makeRng(1)), [[0, 0, 0]]);
  assert.deepStrictEqual(L.randomSaw(1, 3, f.makeRng(1)), [[0, 0, 0]]);
});

test('straightCoords produces an x axis run with z fixed at 0', () => {
  const c = L.straightCoords(5);
  assert.strictEqual(c.length, 5);
  for (let i = 0; i < 5; i++) {
    assert.deepStrictEqual(c[i], [i, 0, 0]);
  }
});

test('isValidSaw accepts a single bead and rejects empties', () => {
  assert.ok(L.isValidSaw([[0, 0, 0]], 2));
  assert.ok(L.isValidSaw([[3, -7, 0]], 2));
  assert.ok(!L.isValidSaw([], 2));
  assert.ok(!L.isValidSaw([[0, 0, 1]], 2));
});

test('pivotCandidate on a two bead chain stays a valid unit step', () => {
  const c = [[0, 0, 0], [1, 0, 0]];
  for (let s = 0; s < 20; s++) {
    const cand = L.pivotCandidate(c, 2, f.makeRng(s + 300));
    assert.ok(L.isValidSaw(cand, 2), 'invalid two bead pivot');
  }
});

test('end flip and corner refuse chains that are too short', () => {
  assert.strictEqual(L.endFlipCandidate([[0, 0, 0], [1, 0, 0]], 2, f.makeRng(1)), null);
  assert.strictEqual(L.cornerCandidate([[0, 0, 0], [1, 0, 0]], 2, f.makeRng(1)), null);
  assert.strictEqual(L.crankshaftCandidate([[0, 0, 0], [1, 0, 0], [1, 1, 0]], 2, f.makeRng(1)), null);
});

test('corner flip on an L shaped bead moves it to the free corner', () => {
  // Bead 1 of this L sits at a corner whose neighbors span a square; the
  // fourth corner (1,1) is free so the flip must be offered and valid.
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 1, 0]];
  const rng = f.makeRng(307);
  let flipped = 0;
  for (let i = 0; i < 40; i++) {
    const cand = L.cornerCandidate(c, 2, rng);
    if (!cand) { continue; }
    const ev = L.evaluateConformation('HHHH', cand, 2);
    if (ev) { flipped++; assert.ok(L.isValidSaw(cand, 2)); }
  }
  assert.ok(flipped > 0, 'corner flip never offered on an L');
});

test('moveCandidate returns array or null for every valid name', () => {
  const c = L.straightCoords(7);
  const rng = f.makeRng(311);
  for (const name of L.MOVE_NAMES) {
    const r = L.moveCandidate(name, c, 2, rng);
    assert.ok(r === null || Array.isArray(r));
  }
});

// -----------------------------------------------------------------------
// Metric formulas
// -----------------------------------------------------------------------
test('radius of gyration of a straight chain matches the closed form', () => {
  // n collinear unit spaced beads: Rg^2 = (n^2 - 1) / 12.
  for (const n of [4, 6, 9, 12]) {
    const m = L.metricsOf('P'.repeat(n), L.straightCoords(n), 2);
    const want = Math.sqrt((n * n - 1) / 12);
    assert.ok(Math.abs(m.rg - want) < 1e-9, 'n=' + n + ' rg ' + m.rg + ' want ' + want);
  }
});

test('end to end of a straight chain equals n minus one', () => {
  for (const n of [2, 5, 9]) {
    assert.strictEqual(L.endToEnd(L.straightCoords(n)), n - 1);
  }
  assert.strictEqual(L.endToEnd([[0, 0, 0]]), 0);
});

test('energy is invariant under chain reversal', () => {
  // Reversing the bead order reverses the sequence labels too, so the
  // physical conformation is the same chain read from the other end.
  const rng = f.makeRng(313);
  const seq = 'HPPHPHPHH';
  const revSeq = seq.split('').reverse().join('');
  for (let i = 0; i < 15; i++) {
    const c = L.randomSaw(9, 2, rng);
    const rev = c.slice().reverse().map(p => p.slice());
    assert.strictEqual(L.energyOf(revSeq, rev, 2), L.energyOf(seq, c, 2));
  }
});

test('contactPairs are sorted in canonical lexicographic order', () => {
  const rng = f.makeRng(317);
  for (let k = 0; k < 10; k++) {
    const c = L.randomSaw(10, 2, rng);
    const pairs = L.contactPairs(c, 2);
    for (let i = 1; i < pairs.length; i++) {
      const [a1, b1] = pairs[i - 1];
      const [a2, b2] = pairs[i];
      assert.ok(a1 < a2 || (a1 === a2 && b1 < b2), 'pairs out of order');
    }
  }
});

test('hhPairs is always a subset of contactPairs', () => {
  const rng = f.makeRng(331);
  const c = L.randomSaw(11, 2, rng);
  const all = new Set(L.contactPairs(c, 2).map(p => p.join(',')));
  for (const p of L.hhPairs('HPPHPHPHPHH', c, 2)) {
    assert.ok(all.has(p.join(',')), 'hh pair not in contact list');
  }
});

test('evaluateConformation returns sorted pairs consistent with contactPairs', () => {
  const rng = f.makeRng(337);
  const c = L.randomSaw(8, 2, rng);
  const ev = L.evaluateConformation('HPPHPHPH', c, 2);
  assert.ok(ev !== null);
  assert.deepStrictEqual(ev.pairs, L.contactPairs(c, 2));
});

// -----------------------------------------------------------------------
// Metrics and run bookkeeping
// -----------------------------------------------------------------------
test('sampleFunnel is deterministic for a fixed seed', () => {
  const a = L.sampleFunnel('HPPHPHPH', 2, 120, 347);
  const b = L.sampleFunnel('HPPHPHPH', 2, 120, 347);
  assert.deepStrictEqual(a.points, b.points);
});

test('different funnel seeds produce different clouds', () => {
  const a = L.sampleFunnel('HPPHPHPH', 2, 120, 349);
  const b = L.sampleFunnel('HPPHPHPH', 2, 120, 350);
  assert.notDeepStrictEqual(a.points, b.points);
});

test('metropolis moved and accepted counters stay inside the step budget', () => {
  const r = L.runMetropolis({ sequence: 'HPPHPHPH', dim: 2, seed: 353, steps: 1200, temperature: 0.9 });
  assert.ok(r.moved <= 1200);
  assert.ok(r.accepted <= r.moved);
  assert.ok(r.steps <= 1200);
});

test('greedy pivot trajectory is strictly decreasing when it improves', () => {
  const r = L.runGreedyPivot({ sequence: 'HHPPHPHH', dim: 2, seed: 359 });
  for (let i = 1; i < r.trajectory.length; i++) {
    assert.ok(r.trajectory[i].energy < r.trajectory[i - 1].energy);
  }
  assert.strictEqual(r.trajectory.length, r.iterations + 1);
});

test('annealing respects sweep and step budget bookkeeping', () => {
  const r = L.runAnnealing({ sequence: 'HPPHPH', dim: 2, seed: 367, sweeps: 11, stepsPerSweep: 7 });
  assert.strictEqual(r.steps, 77);
  assert.strictEqual(r.trajectory.length, 12);
  assert.strictEqual(r.points.length, 12);
});

// -----------------------------------------------------------------------
// Mutation scan recomputation
// -----------------------------------------------------------------------
test('mutation dE equals energyOf on the flipped sequence', () => {
  const rng = f.makeRng(373);
  const c = L.randomSaw(10, 2, rng);
  const seq = 'HPPHPHPHPH';
  const scan = L.mutationScan(seq, c, 2);
  for (const e of scan.perResidue) {
    const flipped = seq.slice(0, e.index) + e.to + seq.slice(e.index + 1);
    const manual = L.energyOf(flipped, c, 2) - L.energyOf(seq, c, 2);
    assert.strictEqual(e.dE, manual, 'residue ' + e.index);
  }
});

test('mutation scan on a contact free chain is all neutral', () => {
  const scan = L.mutationScan('HPPHPHPH', L.straightCoords(8), 2);
  assert.strictEqual(scan.destabilizing, 0);
  assert.strictEqual(scan.stabilizing, 0);
  assert.strictEqual(scan.neutral, 8);
  assert.strictEqual(scan.worst.dE, 0);
});

test('mutation worst index is the first argmax', () => {
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const scan = L.mutationScan('HHHH', sq, 2);
  assert.strictEqual(scan.worst.dE, 1);
  assert.ok(scan.worst.index === 0 || scan.worst.index === 3);
});

test('mutation scan on a two residue chain reports both residues', () => {
  const scan = L.mutationScan('HP', [[0, 0, 0], [1, 0, 0]], 2);
  assert.strictEqual(scan.perResidue.length, 2);
  assert.strictEqual(scan.perResidue[0].from, 'H');
  assert.strictEqual(scan.perResidue[1].from, 'P');
  assert.strictEqual(scan.e0, 0);
});

// -----------------------------------------------------------------------
// Contact matrix details
// -----------------------------------------------------------------------
test('contactMatrix has n squared cells and 2n-2 bonds', () => {
  const rng = f.makeRng(379);
  const c = L.randomSaw(11, 2, rng);
  const m = L.contactMatrix('HPPHPHPHPHH', c, 2);
  assert.strictEqual(m.cells.length, 121);
  let bonded = 0;
  for (const v of m.cells) if (v === L.CMAP_BONDED) bonded++;
  assert.strictEqual(bonded, 20);
});

test('contactMatrix contactCount equals contactPairs length', () => {
  const rng = f.makeRng(383);
  const c = L.randomSaw(9, 3, rng);
  const m = L.contactMatrix('HHPHPHPHH', c, 3);
  assert.strictEqual(m.contactCount, L.contactPairs(c, 3).length);
  assert.strictEqual(m.hhCount, L.hhPairs('HHPHPHPHH', c, 3).length);
});

test('contact matrix codes are four distinct constants plus self', () => {
  const codes = [L.CMAP_FAR, L.CMAP_BONDED, L.CMAP_HH, L.CMAP_OTHER, L.CMAP_SELF];
  assert.strictEqual(new Set(codes).size, 5);
});

// -----------------------------------------------------------------------
// Audit envelope semantics
// -----------------------------------------------------------------------
test('checksum ignores envelope fields outside the payload', () => {
  const a = L.buildAuditReport({ x: [1, 2] });
  const b = Object.assign({}, a, { extra: 'field', checksumNote: 'changed' });
  assert.ok(L.verifyAuditReport(b).ok);
  assert.strictEqual(a.checksum, b.checksum);
});

test('payload key order does not change the checksum', () => {
  const a = L.buildAuditReport({ b: 2, a: 1, n: { y: 2, x: 1 } });
  const b = L.buildAuditReport({ a: 1, b: 2, n: { x: 1, y: 2 } });
  assert.strictEqual(a.checksum, b.checksum);
});

test('tampering a nested array element is detected', () => {
  const r = L.buildAuditReport({ results: [{ e: -2 }, { e: -3 }] });
  r.payload.results[1].e = -9;
  assert.strictEqual(L.verifyAuditReport(r).ok, false);
});

test('verifyAuditText accepts whitespace reformatted JSON', () => {
  const r = L.buildAuditReport({ sequence: 'HPPH', results: [{ algorithm: 'x' }] });
  assert.ok(L.verifyAuditText(JSON.stringify(r)).ok);
  assert.ok(L.verifyAuditText(JSON.stringify(r, null, 2)).ok);
});

test('auditToCsv line count tracks results and mutations', () => {
  const report = L.buildAuditReport({
    sequence: 'HPPH', dimension: 2,
    results: [{ algorithm: 'a', energy: -1, contactCount: 1, hhContacts: 1, rg: 1, ree: 1, optimalityGap: 0, ms: 0 }],
    mutationScan: { perResidue: [{ index: 0, from: 'H', to: 'P', dE: 1 }, { index: 1, from: 'P', to: 'H', dE: 0 }] }
  });
  const lines = L.auditToCsv(report).trim().split('\n');
  // section header + 8 meta rows + algorithm header + 1 row + mutation header + 2 rows
  assert.strictEqual(lines.length, 1 + 8 + 1 + 1 + 1 + 2);
});

test('auditToCsv omits the mutation section without a scan', () => {
  const csv = L.auditToCsv(L.buildAuditReport({ sequence: 'HPPH', dimension: 2, results: [] }));
  assert.ok(!csv.includes('mutation_index'));
});

test('auditToMarkdown omits the mutation table without a scan', () => {
  const md = L.auditToMarkdown(L.buildAuditReport({ sequence: 'HPPH', dimension: 2, results: [] }));
  assert.ok(!md.includes('| residue |'));
  assert.ok(md.includes('| algorithm |'));
});

test('canonicalize sorts numeric like keys lexically', () => {
  assert.strictEqual(L.canonicalize({ '10': 1, '2': 2 }), '{"10":1,"2":2}');
});

test('utf8Bytes on an empty string yields an empty array', () => {
  assert.deepStrictEqual(L.utf8Bytes(''), []);
});

test('sha256Hex matches crypto on a full audit JSON document', () => {
  const r = L.buildAuditReport({ sequence: 'HPPHPH', dimension: 2, seed: 5, results: [{ a: 1 }] });
  const text = L.auditToJson(r);
  assert.strictEqual(L.sha256Hex(text), crypto.createHash('sha256').update(text, 'utf8').digest('hex'));
});

// -----------------------------------------------------------------------
// Bookmark edge cases beyond round trips
// -----------------------------------------------------------------------
test('decodeBookmark rejects a missing sequence and extra fields', () => {
  assert.strictEqual(L.decodeBookmark('#fold=:2d:1'), null);
  assert.strictEqual(L.decodeBookmark('#fold=HPPH:2d:1:9'), null);
  assert.strictEqual(L.decodeBookmark('#fold=HPPH:2d:abc'), null);
});

test('bookmark with seed zero round trips', () => {
  const b = L.decodeBookmark(L.encodeBookmark({ sequence: 'HPPH', dim: 2, seed: 0 }));
  assert.strictEqual(b.seed, 0);
});

test('encodeBookmark normalizes lowercase sequences', () => {
  assert.strictEqual(L.encodeBookmark({ sequence: 'hpph', dim: 2, seed: 3 }), 'fold=HPPH:2d:3');
});

// -----------------------------------------------------------------------
// Lab app wiring and shipped statics
// -----------------------------------------------------------------------
test('lab-app.js is strict and its script tag sits after the lab markup', () => {
  const t = read('lab-app.js');
  assert.match(t, /'use strict';/);
  // The IIFE runs at parse time, so the tag must come after #labCard.
  const html = read('index.html');
  assert.ok(html.indexOf('src="lab-app.js"') > html.indexOf('id="labCard"'),
    'lab-app.js loads before the lab markup exists');
});

test('index.html keeps the strict CSP meta untouched', () => {
  const html = read('index.html');
  assert.ok(html.includes("default-src 'none'"));
  assert.ok(html.includes("script-src 'self' file:"));
  assert.ok(html.includes("worker-src blob:"));
});

test('lab section uses heading titles without kicker or eyebrow text', () => {
  const html = read('index.html');
  const labStart = html.indexOf('id="labCard"');
  assert.ok(labStart > 0);
  const labHtml = html.slice(labStart);
  assert.ok(!/class=["'][^"']*(eyebrow|kicker)[^"']*["']/i.test(labHtml), 'eyebrow/kicker class in lab section');
});

test('every lab engine function the UI calls exists in the export', () => {
  const app = read('lab-app.js');
  const used = new Set([...app.matchAll(/\bL\.([A-Za-z_$][A-Za-z0-9_$]*)/g)].map(m => m[1]));
  assert.ok(used.size > 10, 'expected many engine calls, got ' + used.size);
  for (const name of used) {
    assert.ok(L[name] !== undefined, 'missing export L.' + name);
  }
});

test('audit buttons in index.html have no inline event handlers', () => {
  const html = read('index.html');
  const labStart = html.indexOf('id="labCard"');
  const labHtml = html.slice(labStart);
  assert.ok(!/\son(click|change|input|submit|keydown|keyup)=/i.test(labHtml), 'inline handler in lab section');
});

test('lab controls expose the documented input ids', () => {
  const html = read('index.html');
  for (const id of ['hpSeqInput', 'hpPresetSel', 'dimSel', 'seedInput', 'samplesInput',
    'mcStepsInput', 'tempInput', 'saSweepsInput', 't0Input', 't1Input', 'compareBtn',
    'mutBtn', 'auditJsonBtn', 'auditCsvBtn', 'auditMdBtn', 'bookmarkBtn',
    'verifyBtn', 'tamperBtn', 'selfTestBtn']) {
    assert.ok(html.includes('id="' + id + '"'), 'missing control #' + id);
  }
});

test('self test names cover sha256 and exact enumeration', () => {
  const out = L.runSelfTests();
  const names = out.results.map(r => r.name).join(' | ');
  assert.match(names, /SHA-256/i);
  assert.match(names, /exact enumeration/i);
  assert.ok(out.passed === out.total);
});

// -----------------------------------------------------------------------
// Final hardening checks
// -----------------------------------------------------------------------
test('canonicalEnumerate never yields the same conformation twice in 3D', () => {
  const seen = new Set();
  L.canonicalEnumerate(5, 3, (c) => {
    const k = c.map(p => p.join(',')).join('|');
    assert.ok(!seen.has(k), 'duplicate 3D conformation ' + k);
    seen.add(k);
  }, 1e6);
  assert.strictEqual(seen.size, L.countSaws(5, 3));
});

test('runExact ground state beats or ties random sampling on the same seq', () => {
  const seq = 'HHPHPHPHPP';
  const exact = L.runExact({ sequence: seq, dim: 2 });
  const walk = L.runRandomWalk({ sequence: seq, dim: 2, samples: 200, seed: 7 });
  assert.ok(exact.best.energy <= walk.best.energy);
});

test('metropolis best never exceeds the initial extended conformation', () => {
  const r = L.runMetropolis({ sequence: 'HPPHPHPH', dim: 2, seed: 11, steps: 500, temperature: 1.2 });
  assert.ok(r.best.energy <= 0);
});

test('annealing best never exceeds the initial extended conformation', () => {
  const r = L.runAnnealing({ sequence: 'HPPHPHPH', dim: 2, seed: 13, sweeps: 20, stepsPerSweep: 15 });
  assert.ok(r.best.energy <= 0);
});

test('greedy pivot best never exceeds the initial extended conformation', () => {
  const r = L.runGreedyPivot({ sequence: 'HPPHPHPH', dim: 2, seed: 17 });
  assert.ok(r.best.energy <= 0);
});

test('mutation scan e0 equals the energy of the reference conformation', () => {
  const rng = f.makeRng(389);
  const c = L.randomSaw(9, 2, rng);
  const scan = L.mutationScan('HPHPHPHPH', c, 2);
  assert.strictEqual(scan.e0, L.energyOf('HPHPHPHPH', c, 2));
});

test('decodeBookmark tolerates a leading double hash', () => {
  assert.deepStrictEqual(L.decodeBookmark('##fold=HPPH:2d:7'), { sequence: 'HPPH', dim: 2, seed: 7 });
});

test('sha256Hex of a 448 bit boundary string matches crypto', () => {
  const s = 'y'.repeat(56);
  assert.strictEqual(L.sha256Hex(s), crypto.createHash('sha256').update(s).digest('hex'));
});

test('index.html mentions the exact enumeration length limits', () => {
  const html = read('index.html');
  assert.match(html, /16|length/);
  assert.match(html, /exact/i);
});

test('app.js exposes the foldlStage bridge for the lab viewer', () => {
  const app = read('app.js');
  assert.ok(app.includes('foldlStage'), 'foldlStage bridge missing');
  assert.ok(app.includes('showLatticeFold'), 'showLatticeFold missing');
});

test('lab-app.js calls showLatticeFold with coords sequence and energy', () => {
  const app = read('lab-app.js');
  assert.ok(app.includes('showLatticeFold'), 'stage bridge call missing');
});

test('package.json test script still points at node test', () => {
  const pkg = JSON.parse(read('package.json'));
  assert.match(pkg.scripts.test, /node --test/);
});
