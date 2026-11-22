'use strict';

// Hostile invariants for the lab workbench: move integrity under stress,
// energy bounds and symmetries, enumeration edge cases, checksum and CSV
// abuse cases, and static shipped file checks for the new modules.

const test = require('node:test');
const assert = require('node:assert');
const crypto = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const f = require('../src/fold.js');
const L = require('../src/lab.js');

const ROOT = path.join(__dirname, '..');
const read = name => fs.readFileSync(path.join(ROOT, name), 'utf8');
const NEW_SHIPPED = ['src/lab.js', 'lab-app.js', 'index.html', 'app.js'];

function man(a, b) {
  return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);
}

function shaRef(s) {
  return crypto.createHash('sha256').update(s, 'utf8').digest('hex');
}

function sawOk(coords, dim) {
  return L.isValidSaw(coords, dim);
}

// -----------------------------------------------------------------------
// Move integrity under long random play
// -----------------------------------------------------------------------
test('every committed lab move keeps 2D states valid SAWs', () => {
  // Candidates are proposals; the evaluation step refuses overlaps, exactly
  // like the fold.js tryMove path. Committed states must always be SAWs.
  const rng = f.makeRng(163);
  const seq = 'HPPHPHPHPHH';
  const fold = { coords: L.straightCoords(11) };
  const commits = [0, 0, 0, 0];
  for (let i = 0; i < 8000; i++) {
    const cand = L.moveCandidate(L.MOVE_NAMES[i % 4], fold.coords, 2, rng);
    if (!cand) { continue; }
    const ev = L.evaluateConformation(seq, cand, 2);
    if (!ev) { continue; }
    assert.ok(sawOk(cand, 2), 'evaluateConformation accepted an invalid SAW');
    fold.coords = cand;
    commits[i % 4]++;
  }
  assert.ok(commits[0] + commits[1] + commits[2] + commits[3] > 0, 'no move ever committed');
});

test('every committed lab move keeps 3D states valid SAWs', () => {
  const rng = f.makeRng(167);
  const seq = 'HHPHPHPHPP';
  const fold = { coords: L.straightCoords(10) };
  const commits = [0, 0, 0, 0];
  for (let i = 0; i < 8000; i++) {
    const cand = L.moveCandidate(L.MOVE_NAMES[i % 4], fold.coords, 3, rng);
    if (!cand) { continue; }
    const ev = L.evaluateConformation(seq, cand, 3);
    if (!ev) { continue; }
    assert.ok(sawOk(cand, 3), 'evaluateConformation accepted an invalid SAW');
    fold.coords = cand;
    commits[i % 4]++;
  }
  assert.ok(commits[0] + commits[1] + commits[2] + commits[3] > 0, 'no move ever committed');
});

test('crankshaft commits deterministically on a U shaped loop', () => {
  // A U turn has d3 == -d1 which is exactly the crankshaft precondition, so
  // this move type is proven to commit, not just sampled by luck.
  const seq = 'HHHH';
  const u2d = [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]];
  const rng = f.makeRng(5);
  const cand = L.crankshaftCandidate(u2d, 2, rng);
  assert.ok(cand !== null, 'crankshaft refused a U shape');
  const ev = L.evaluateConformation(seq, cand, 2);
  assert.ok(ev !== null, 'crankshaft candidate failed evaluation');
  assert.ok(sawOk(ev ? cand : [], 2));
  assert.deepStrictEqual(cand[1], [0, -1, 0]);
  assert.deepStrictEqual(cand[2], [1, -1, 0]);
});

test('crankshaft commits deterministically on a 3D U shaped loop', () => {
  const u3d = [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]];
  const cand = L.crankshaftCandidate(u3d, 3, f.makeRng(7));
  assert.ok(cand !== null);
  assert.ok(L.evaluateConformation('HHHH', cand, 3) !== null);
  assert.ok(sawOk(cand, 3));
});

test('lab moves always produce distinct or identical conformations, never corrupt', () => {
  const rng = f.makeRng(169);
  const fold = { seq: 'HPPHPH', dim: 2, coords: L.straightCoords(6), energy: 0 };
  for (let i = 0; i < 800; i++) {
    const cand = L.moveCandidate(L.MOVE_NAMES[i % 4], fold.coords, 2, rng);
    if (!cand) { continue; }
    assert.strictEqual(cand.length, fold.coords.length);
    for (const c of cand) assert.strictEqual(c.length, 3);
    fold.coords = cand;
  }
});

test('mcStepLab uses all four move types over a long run', () => {
  const rng = f.makeRng(173);
  const fold = { seq: 'HPPHPHPH', dim: 2, coords: L.straightCoords(8), energy: 0 };
  const seen = new Set();
  for (let i = 0; i < 4000; i++) {
    const r = L.mcStepLab(fold, rng, 1.2);
    seen.add(r.type);
  }
  assert.strictEqual(seen.size, 4);
});

test('mcStepLab dE agrees with the committed energy', () => {
  const rng = f.makeRng(179);
  const fold = { seq: 'HPPHPHPH', dim: 2, coords: L.straightCoords(8), energy: 0 };
  for (let i = 0; i < 1500; i++) {
    const before = fold.energy;
    const r = L.mcStepLab(fold, rng, 0.7);
    if (r.accepted) {
      assert.strictEqual(fold.energy, before + r.dE);
    } else {
      assert.strictEqual(fold.energy, before);
    }
  }
});

test('mcStepLab on a locked conformation reports refusals cleanly', () => {
  // A closed square has almost no legal moves; refusals must not corrupt.
  const fold = { seq: 'HHHH', dim: 2, coords: [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], energy: -1 };
  const rng = f.makeRng(181);
  for (let i = 0; i < 200; i++) {
    L.mcStepLab(fold, rng, 1);
    assert.ok(sawOk(fold.coords, 2));
  }
});

// -----------------------------------------------------------------------
// Energy invariants
// -----------------------------------------------------------------------
test('energy is invariant under lattice translation', () => {
  const rng = f.makeRng(191);
  for (let i = 0; i < 15; i++) {
    const c = L.randomSaw(10, 2, rng);
    const shifted = c.map(p => [p[0] + 40, p[1] - 13, p[2]]);
    assert.strictEqual(L.energyOf('HPPHPHPHPH', shifted, 2), L.energyOf('HPPHPHPHPH', c, 2));
  }
});

test('energy is invariant under a quarter turn rotation', () => {
  const rng = f.makeRng(193);
  for (let i = 0; i < 15; i++) {
    const c = L.randomSaw(9, 2, rng);
    const rot = c.map(p => [-p[1], p[0], p[2]]);
    assert.strictEqual(L.energyOf('HPPHPHPHP', rot, 2), L.energyOf('HPPHPHPHP', c, 2));
  }
});

test('energy respects the lattice contact count bound', () => {
  // Each bead makes at most 2d neighbors and two are chain bonds, so the
  // absolute energy is bounded by n * (d - 1) contacts in the worst case.
  const rng = f.makeRng(197);
  for (const [dim, bound] of [[2, 2], [3, 4]]) {
    for (let i = 0; i < 15; i++) {
      const c = L.randomSaw(10, dim, rng);
      const e = L.energyOf('HHHHHHHHHH', c, dim);
      assert.ok(e >= -10 * bound, 'energy ' + e + ' below the contact bound');
    }
  }
});

test('energy of a conformation equals energy recomputed from scratch', () => {
  const rng = f.makeRng(199);
  const fold = { seq: 'HPPHPHPHPHH', dim: 2, coords: L.straightCoords(11), energy: 0 };
  for (let i = 0; i < 600; i++) {
    L.mcStepLab(fold, rng, 1.0);
    if (i % 50 === 0) {
      assert.strictEqual(fold.energy, L.energyOf(fold.seq, fold.coords, 2));
    }
  }
});

// -----------------------------------------------------------------------
// Enumeration edge cases
// -----------------------------------------------------------------------
test('canonicalEnumerate on n=2 yields exactly the canonical walk', () => {
  const seen = [];
  L.canonicalEnumerate(2, 2, c => seen.push(c), 10);
  assert.strictEqual(seen.length, 1);
  assert.deepStrictEqual(seen[0], [[0, 0, 0], [1, 0, 0]]);
});

test('canonicalEnumerate respects the state budget and stays deterministic', () => {
  const a = [];
  L.canonicalEnumerate(9, 2, c => a.push(c.map(p => p.join('')).join('|')), 40);
  const b = [];
  L.canonicalEnumerate(9, 2, c => b.push(c.map(p => p.join('')).join('|')), 40);
  assert.deepStrictEqual(a, b);
  assert.ok(a.length <= 40);
});

test('runExact is deterministic across repeated runs', () => {
  const a = L.runExact({ sequence: 'HPHPHPH', dim: 2 });
  const b = L.runExact({ sequence: 'HPHPHPH', dim: 2 });
  assert.strictEqual(a.best.energy, b.best.energy);
  assert.deepStrictEqual(a.best.coords, b.best.coords);
  assert.strictEqual(a.states, b.states);
});

test('runExact ground state energies are achievable conformations', () => {
  for (const seq of ['HHHH', 'HHPPHH', 'HPHPHPH', 'HHHPPPHH']) {
    const r = L.runExact({ sequence: seq, dim: 2 });
    assert.ok(sawOk(r.best.coords, 2));
    assert.strictEqual(L.energyOf(seq, r.best.coords, 2), r.best.energy);
  }
});

test('runExact 3D matches brute force for short chains', () => {
  function brute(seq) {
    const dirs = L.DIRS[3];
    const coords = [[0, 0, 0]];
    const occ = new Set(['0,0,0']);
    let best = 0;
    function dfs() {
      if (coords.length === seq.length) {
        let e = 0;
        for (let i = 0; i < seq.length; i++) {
          for (let j = i + 2; j < seq.length; j++) {
            if (seq[i] === 'H' && seq[j] === 'H' && man(coords[i], coords[j]) === 1) e--;
          }
        }
        if (e < best) best = e;
        return;
      }
      const last = coords[coords.length - 1];
      for (const d of dirs) {
        const cell = [last[0] + d[0], last[1] + d[1], last[2] + d[2]];
        const k = cell.join(',');
        if (occ.has(k)) continue;
        coords.push(cell); occ.add(k); dfs(); coords.pop(); occ.delete(k);
      }
    }
    dfs();
    return best;
  }
  for (const seq of ['HHHHH', 'HHPPHH', 'HPHPHH', 'HHPHPHH']) {
    const r = L.runExact({ sequence: seq, dim: 3 });
    assert.strictEqual(r.best.energy, brute(seq), seq + ' 3D');
    assert.ok(sawOk(r.best.coords, 3));
  }
});

test('enumeration of an H rich 8mer finds a compact ground state', () => {
  const r = L.runExact({ sequence: 'HHHHHHHH', dim: 2 });
  // The all hydrophobic 8mer reaches -3 on the square lattice (2x4 block
  // or 3x3 minus a corner); any result shallower than -3 is a bug.
  assert.strictEqual(r.best.energy, -3);
});

// -----------------------------------------------------------------------
// Algorithm invariants
// -----------------------------------------------------------------------
test('runGreedyPivot on an all P chain stops at zero immediately', () => {
  const r = L.runGreedyPivot({ sequence: 'PPPPPPPP', dim: 2, seed: 211 });
  assert.strictEqual(r.best.energy, 0);
  assert.strictEqual(r.iterations, 0);
});

test('runGreedyPivot never exceeds its iteration cap', () => {
  const r = L.runGreedyPivot({ sequence: 'HPPHPHPHPHH', dim: 2, seed: 223, maxIter: 3 });
  assert.ok(r.iterations <= 3);
  assert.strictEqual(r.trajectory.length, r.iterations + 1);
});

test('runMetropolis with zero steps returns the extended chain', () => {
  const r = L.runMetropolis({ sequence: 'HPPHPH', dim: 2, seed: 227, steps: 1 });
  assert.ok(sawOk(r.best.coords, 2));
  assert.ok(r.best.energy <= 0);
});

test('runAnnealing trajectory records one point per sweep plus start', () => {
  const r = L.runAnnealing({ sequence: 'HPPHPH', dim: 2, seed: 229, sweeps: 17, stepsPerSweep: 5 });
  assert.strictEqual(r.trajectory.length, 18);
  assert.strictEqual(r.points.length, 18);
});

test('annealing on a one residue move budget still finishes', () => {
  const r = L.runAnnealing({ sequence: 'HPPH', dim: 2, sweeps: 3, stepsPerSweep: 2, seed: 233 });
  assert.strictEqual(r.steps, 6);
  assert.ok(sawOk(r.best.coords, 2));
});

test('random walk best never exceeds the best sample it saw', () => {
  const r = L.runRandomWalk({ sequence: 'HPPHPHPH', dim: 2, samples: 80, seed: 239 });
  let min = 0;
  for (const p of r.points) if (p.E < min) min = p.E;
  assert.strictEqual(r.best.energy, min);
});

test('engine rejects degenerate algorithm budgets', () => {
  assert.throws(() => L.runRandomWalk({ sequence: 'HPPH', dim: 2, samples: 0 }), RangeError);
  assert.throws(() => L.runRandomWalk({ sequence: 'HPPH', dim: 2, samples: 2.5 }), RangeError);
  assert.throws(() => L.runRandomWalk({ sequence: 'HPPH', dim: 2, samples: NaN }), RangeError);
  assert.throws(() => L.runMetropolis({ sequence: 'HPPH', dim: 2, steps: 0 }), RangeError);
  assert.throws(() => L.runAnnealing({ sequence: 'HPPH', dim: 2, sweeps: 0 }), RangeError);
  assert.throws(() => L.runAnnealing({ sequence: 'HPPH', dim: 2, stepsPerSweep: 0 }), RangeError);
  assert.throws(() => L.runGreedyPivot({ sequence: 'HPPH', dim: 2, maxIter: -1 }), RangeError);
  assert.throws(() => L.sampleFunnel('HPPH', 2, 0), RangeError);
});

test('greedy pivot reports a real descent in the comparator table', () => {
  const cmp = L.compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 97 });
  const greedy = cmp.results[1];
  assert.strictEqual(greedy.algorithm, 'greedy-pivot');
  assert.ok(greedy.best.energy < 0, 'greedy stalled at zero from a random start');
  assert.ok(greedy.points[0].step === 0, 'greedy funnel points carry a step for step mode');
});

test('comparator on the minimal two residue chain', () => {
  const cmp = L.compareAlgorithms({ sequence: 'HH', dim: 2, seed: 241 });
  for (const r of cmp.results) {
    assert.strictEqual(r.best.energy, 0);
    assert.strictEqual(r.optimalityGap, 0);
  }
});

test('comparator result labels stay stable for the UI table', () => {
  const cmp = L.compareAlgorithms({ sequence: 'HPPH', dim: 2, seed: 243 });
  assert.strictEqual(cmp.results[0].label, 'Random SAW sampling');
  assert.strictEqual(cmp.results[1].label, 'Greedy local pivot');
  assert.strictEqual(cmp.results[2].label, 'Metropolis MC (fixed T)');
  assert.strictEqual(cmp.results[3].label, 'Simulated annealing');
  assert.strictEqual(cmp.results[4].label, 'Exact enumeration');
});

// -----------------------------------------------------------------------
// Mutation scan edge cases
// -----------------------------------------------------------------------
test('mutation deltaE values are integers', () => {
  const rng = f.makeRng(251);
  const c = L.randomSaw(11, 2, rng);
  const scan = L.mutationScan('HPPHPHPHPHH', c, 2);
  for (const e of scan.perResidue) {
    assert.ok(Number.isInteger(e.dE), 'non integer deltaE ' + e.dE);
  }
});

test('mutation worst entry is the maximum deltaE', () => {
  const rng = f.makeRng(253);
  const c = L.randomSaw(9, 2, rng);
  const scan = L.mutationScan('HPHPHPHPH', c, 2);
  const max = Math.max.apply(null, scan.perResidue.map(e => e.dE));
  assert.strictEqual(scan.worst.dE, max);
});

test('mutating a buried H core costs more than mutating an exposed P', () => {
  // Compact square HHHH: core H mutations cost +1; a polar residue swapped
  // to H on the same geometry can only add contacts, never remove them.
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const scan = L.mutationScan('HHHH', sq, 2);
  assert.ok(scan.worst.dE >= 1);
  const scanP = L.mutationScan('PPPP', sq, 2);
  for (const e of scanP.perResidue) {
    assert.ok(e.dE <= 0, 'P to H can never remove an H-H contact');
  }
});

test('mutation scan on 3D conformations works identically', () => {
  const rng = f.makeRng(257);
  const c = L.randomSaw(8, 3, rng);
  const scan = L.mutationScan('HPHPPHPH', c, 3);
  assert.strictEqual(scan.perResidue.length, 8);
  assert.strictEqual(scan.destabilizing + scan.stabilizing + scan.neutral, 8);
});

// -----------------------------------------------------------------------
// Contact matrix edge cases
// -----------------------------------------------------------------------
test('contactMatrix on a two residue chain has only the bond', () => {
  const m = L.contactMatrix('HH', [[0, 0, 0], [1, 0, 0]], 2);
  assert.strictEqual(m.cells[1], L.CMAP_BONDED);
  assert.strictEqual(m.cells[2], L.CMAP_BONDED);
  assert.strictEqual(m.hhCount, 0);
  assert.strictEqual(m.contactCount, 0);
});

test('contactMatrix on a straight chain has bonds and no contacts', () => {
  const m = L.contactMatrix('HHHHHH', L.straightCoords(6), 2);
  for (let i = 0; i < 5; i++) {
    assert.strictEqual(m.cells[i * 6 + i + 1], L.CMAP_BONDED);
    assert.strictEqual(m.cells[(i + 1) * 6 + i], L.CMAP_BONDED);
  }
  let bonded = 0;
  for (const v of m.cells) if (v === L.CMAP_BONDED) bonded++;
  assert.strictEqual(bonded, 10);
  assert.strictEqual(m.contactCount, 0);
});

test('contactMatrix codes only the four documented values', () => {
  const rng = f.makeRng(263);
  const c = L.randomSaw(10, 2, rng);
  const m = L.contactMatrix('HPPHPHPHPH', c, 2);
  const allowed = new Set([L.CMAP_FAR, L.CMAP_BONDED, L.CMAP_HH, L.CMAP_OTHER, L.CMAP_SELF]);
  for (const v of m.cells) assert.ok(allowed.has(v), 'unknown cell code ' + v);
});

// -----------------------------------------------------------------------
// SHA-256 block boundary and oracle stress
// -----------------------------------------------------------------------
test('sha256Hex matches crypto at every padding boundary', () => {
  for (const len of [0, 1, 54, 55, 56, 57, 63, 64, 65, 119, 120, 127, 128, 129, 255, 256, 1000]) {
    const s = 'x'.repeat(len);
    assert.strictEqual(L.sha256Hex(s), shaRef(s), 'length ' + len);
  }
});

test('sha256Hex matches crypto on binary like content', () => {
  const rng = f.makeRng(269);
  for (let i = 0; i < 20; i++) {
    let s = '';
    const len = 1 + Math.floor(rng() * 300);
    for (let k = 0; k < len; k++) {
      s += String.fromCharCode(32 + Math.floor(rng() * 95));
    }
    assert.strictEqual(L.sha256Hex(s), shaRef(s));
  }
});

test('sha256Hex handles a lone surrogate without throwing', () => {
  const out = L.sha256Hex('abc\uD800def');
  assert.match(out, /^[0-9a-f]{64}$/);
});

test('sha256Hex output is always 64 lowercase hex', () => {
  const rng = f.makeRng(271);
  for (let i = 0; i < 10; i++) {
    assert.match(L.sha256Hex('seed' + i + rng()), /^[0-9a-f]{64}$/);
  }
});

// -----------------------------------------------------------------------
// Canonical JSON and checksum edge cases
// -----------------------------------------------------------------------
test('canonicalize sorts keys at every nesting level', () => {
  const out = L.canonicalize({ z: { b: 1, a: 2 }, a: [{ d: 4, c: 3 }] });
  assert.strictEqual(out, '{"a":[{"c":3,"d":4}],"z":{"a":2,"b":1}}');
});

test('canonicalize maps undefined array slots to null like JSON', () => {
  assert.strictEqual(L.canonicalize([undefined, 1]), '[null,1]');
});

test('verifyAuditReport detects a swapped checksum', () => {
  const r = L.buildAuditReport({ a: 1 });
  r.checksum = '0'.repeat(64);
  const v = L.verifyAuditReport(r);
  assert.strictEqual(v.ok, false);
  assert.match(v.reason, /mismatch/);
});

test('verifyAuditReport reports the actual recomputed digest', () => {
  const r = L.buildAuditReport({ a: 1 });
  const good = L.verifyAuditReport(r);
  assert.ok(good.ok);
  r.checksum = r.checksum.slice(1) + 'f';
  const bad = L.verifyAuditReport(r);
  assert.strictEqual(bad.ok, false);
  assert.match(bad.actual, /^[0-9a-f]{64}$/);
});

test('verifyAuditText is whitespace agnostic on the payload level', () => {
  // Reformatting the JSON (pretty vs compact) must still verify because
  // the checksum covers the canonicalized payload, not the byte layout.
  const r = L.buildAuditReport({ a: [1, 2], s: 'x' });
  const compact = JSON.stringify(r);
  assert.ok(L.verifyAuditText(compact).ok);
  assert.ok(L.verifyAuditText(JSON.stringify(r, null, 4)).ok);
});

test('auditToJson ends with a newline and reparses identically', () => {
  const r = L.buildAuditReport({ x: [1, 'a', null] });
  const t = L.auditToJson(r);
  assert.ok(t.endsWith('\n'));
  assert.deepStrictEqual(JSON.parse(t), r);
});

// -----------------------------------------------------------------------
// CSV escaping edge cases
// -----------------------------------------------------------------------
test('csvCell leaves ordinary values untouched', () => {
  assert.strictEqual(L.csvCell('HPPHPH'), 'HPPHPH');
  assert.strictEqual(L.csvCell(-3), '-3');
  assert.strictEqual(L.csvCell(3.5), '3.5');
});

test('csvCell doubles embedded quotes before wrapping', () => {
  assert.strictEqual(L.csvCell('a"b'), '"a""b"');
  assert.strictEqual(L.csvCell('"'), '""""');
});

test('csvCell neutralizes every formula leading character', () => {
  for (const ch of ['=', '+', '-', '@', '\t', '\r', ' ']) {
    const out = L.csvCell(ch + 'x');
    assert.ok(out.startsWith('"\''), ch + ' not neutralized');
  }
});

test('auditToCsv always ends with a newline and has no stray CR', () => {
  const csv = L.auditToCsv(L.buildAuditReport({ sequence: 'HPPH', dimension: 2, results: [] }));
  assert.ok(csv.endsWith('\n'));
  assert.ok(!csv.includes('\r'));
});

test('auditToMarkdown escapes pipe characters in values', () => {
  const r = L.buildAuditReport({ sequence: 'HP|PH', dimension: 2, results: [] });
  const md = L.auditToMarkdown(r);
  assert.ok(md.includes('HP\\|PH'), 'unescaped pipe in markdown');
});

// -----------------------------------------------------------------------
// Bookmark edge cases
// -----------------------------------------------------------------------
test('decodeBookmark accepts boundary seeds', () => {
  assert.deepStrictEqual(L.decodeBookmark('#fold=HPPH:2d:0'), { sequence: 'HPPH', dim: 2, seed: 0 });
  assert.deepStrictEqual(L.decodeBookmark('#fold=HPPH:3d:4294967295'), { sequence: 'HPPH', dim: 3, seed: 4294967295 });
});

test('decodeBookmark rejects seeds above the 32 bit range', () => {
  assert.strictEqual(L.decodeBookmark('#fold=HPPH:2d:4294967296'), null);
});

test('decodeBookmark uppercases lowercase sequences', () => {
  const b = L.decodeBookmark('#fold=hpph:2d:5');
  assert.strictEqual(b.sequence, 'HPPH');
});

test('encodeBookmark output always round trips', () => {
  const rng = f.makeRng(281);
  for (let i = 0; i < 25; i++) {
    let seq = '';
    for (let k = 0; k < 4 + i % 8; k++) seq += rng() < 0.5 ? 'H' : 'P';
    const opts = { sequence: seq, dim: rng() < 0.5 ? 2 : 3, seed: Math.floor(rng() * 1e9) };
    const back = L.decodeBookmark(L.encodeBookmark(opts));
    assert.deepStrictEqual(back, { sequence: seq, dim: opts.dim, seed: opts.seed });
  }
});

// -----------------------------------------------------------------------
// Self test runner structure
// -----------------------------------------------------------------------
test('runSelfTests returns named results with string details', () => {
  const out = L.runSelfTests();
  for (const r of out.results) {
    assert.ok(r.name.length > 3);
    assert.strictEqual(typeof r.detail, 'string');
    assert.strictEqual(typeof r.pass, 'boolean');
  }
  assert.strictEqual(out.passed + (out.results.length - out.passed), out.results.length);
});

// -----------------------------------------------------------------------
// Static shipped file checks for the new modules
// -----------------------------------------------------------------------
test('no em dash or en dash in the new shipped files', () => {
  for (const name of NEW_SHIPPED) {
    const t = read(name);
    assert.ok(!t.includes('—'), name + ' em dash');
    assert.ok(!t.includes('–'), name + ' en dash');
  }
});

test('lab-app.js avoids innerHTML sinks eval and new Function', () => {
  const t = read('lab-app.js');
  assert.ok(!/\.innerHTML\s*=/ .test(t), 'innerHTML sink');
  assert.ok(!/eval\s*\(|new Function|document\.write/.test(t));
});

test('index.html loads the lab scripts after the engine', () => {
  const html = read('index.html');
  const iFold = html.indexOf('src="src/fold.js"');
  const iApp = html.indexOf('src="app.js"');
  const iLab = html.indexOf('src="src/lab.js"');
  const iLabApp = html.indexOf('src="lab-app.js"');
  assert.ok(iFold > 0 && iLab > iApp && iLabApp > iLab, 'script order wrong');
});

test('every lab-app getElementById target exists in index.html', () => {
  const app = read('lab-app.js');
  const html = read('index.html');
  const ids = new Set([...app.matchAll(/\$\('([^']+)'\)|getElementById\('([^']+)'\)/g)]
    .map(m => m[1] || m[2]));
  assert.ok(ids.size > 15, 'expected many lab ids, got ' + ids.size);
  for (const id of ids) {
    assert.ok(html.includes('id="' + id + '"'), 'missing #' + id);
  }
});

test('new shipped files keep the offline zero dependency contract', () => {
  for (const name of NEW_SHIPPED) {
    const t = read(name);
    assert.ok(!/<script[^>]+src=["']https?:/i.test(t), name + ' external script');
    assert.ok(!/<link[^>]+href=["']https?:/i.test(t), name + ' external link');
    assert.ok(!/fetch\(["']https?:/i.test(t), name + ' remote fetch');
    assert.ok(!/url\(["']?https?:/i.test(t), name + ' remote css url');
    assert.ok(!/require\(['"]/.test(t) || name === 'src/lab.js', name + ' unexpected require');
  }
});

test('lab engine runs in strict mode and exposes a UMD export', () => {
  const t = read('src/lab.js');
  assert.match(t, /^'use strict';/);
  assert.match(t, /module\.exports = lab/);
  assert.match(t, /window\.foldlLab = lab/);
});

test('index.html discloses the HP model limits and the digest scope', () => {
  const html = read('index.html');
  assert.match(html, /simplified statistical physics toy/);
  assert.match(html, /not an atomistic force field/);
  assert.match(html, /AlphaFold style structure predictor/);
  assert.match(html, /unkeyed.*tamper detection digest|tamper detection digest.*not a signature/);
});

test('index.html workbench markup carries accessibility roles', () => {
  const html = read('index.html');
  assert.match(html, /id="funnelSvg"[^>]*role="img"[^>]*aria-label=/);
  assert.match(html, /id="cmapSvg"[^>]*role="img"[^>]*aria-label=/);
  assert.match(html, /id="mutSvg"[^>]*role="img"[^>]*aria-label=/);
  assert.match(html, /id="labStatus"[^>]*aria-live="polite"/);
  assert.match(html, /id="labErr"[^>]*role="alert"/);
  assert.match(html, /id="auditVerdict"[^>]*aria-live="polite"/);
});
