'use strict';

// Lab workbench suite: HP model on 2D and 3D lattices, the five algorithm
// comparator, exact enumeration, mutation scans, contact matrices,
// SHA-256, checksummed audits, bookmarks and self tests.

const test = require('node:test');
const assert = require('node:assert');
const crypto = require('node:crypto');
const f = require('../src/fold.js');
const L = require('../src/lab.js');

const SEQ2 = 'HPPHPHPHPHH';
const SEQ3 = 'HHPHPHPHPP';

function man(a, b) {
  return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);
}

function shaRef(s) {
  return crypto.createHash('sha256').update(s, 'utf8').digest('hex');
}

// -----------------------------------------------------------------------
// HP sequence and dimension validation
// -----------------------------------------------------------------------
test('validateHpSequence normalizes case and accepts H and P', () => {
  assert.strictEqual(L.validateHpSequence('hpph'), 'HPPH');
  assert.strictEqual(L.validateHpSequence('HPPHPH'), 'HPPHPH');
  assert.strictEqual(L.validateHpSequence('HP HP\n'), 'HPHP');
});

test('validateHpSequence rejects non HP letters', () => {
  assert.throws(() => L.validateHpSequence('HPXH'), /only H and P/);
  assert.throws(() => L.validateHpSequence('HP1H'), /only H and P/);
  assert.throws(() => L.validateHpSequence(''), /at least 2/);
  assert.throws(() => L.validateHpSequence('H'), /at least 2/);
});

test('validateHpSequence rejects non strings', () => {
  assert.throws(() => L.validateHpSequence(12), /must be a string/);
  assert.throws(() => L.validateHpSequence(['H', 'P']), /must be a string/);
  assert.throws(() => L.validateHpSequence(null), /must be a string/);
});

test('validateDim accepts exactly 2 and 3', () => {
  assert.strictEqual(L.validateDim(2), 2);
  assert.strictEqual(L.validateDim(3), 3);
  assert.throws(() => L.validateDim(1), /2.*3/);
  assert.throws(() => L.validateDim(4), /2.*3/);
  assert.throws(() => L.validateDim('2'), /2.*3/);
});

// -----------------------------------------------------------------------
// SAW validity
// -----------------------------------------------------------------------
test('isValidSaw accepts a straight walk in both dims', () => {
  const c = L.straightCoords(7);
  assert.ok(L.isValidSaw(c, 2));
  assert.ok(L.isValidSaw(c, 3));
});

test('isValidSaw rejects overlapping cells', () => {
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0]];
  assert.strictEqual(L.isValidSaw(c, 2), false);
});

test('isValidSaw rejects non unit steps', () => {
  assert.strictEqual(L.isValidSaw([[0, 0, 0], [2, 0, 0], [3, 0, 0]], 2), false);
  assert.strictEqual(L.isValidSaw([[0, 0, 0], [1, 1, 0], [2, 1, 0]], 2), false);
});

test('isValidSaw rejects non integer and out of range cells', () => {
  assert.strictEqual(L.isValidSaw([[0, 0, 0], [0.5, 0, 0], [1.5, 0, 0]], 2), false);
  assert.strictEqual(L.isValidSaw([[0, 0, 0], [1, 0, 0], [100001, 0, 0]], 2), false);
});

test('isValidSaw enforces z = 0 on the 2D lattice', () => {
  assert.strictEqual(L.isValidSaw([[0, 0, 0], [1, 0, 0], [1, 0, 1]], 2), false);
  assert.strictEqual(L.isValidSaw([[0, 0, 0], [1, 0, 0], [1, 0, 1]], 3), true);
});

test('isValidSaw rejects malformed input', () => {
  assert.strictEqual(L.isValidSaw(null, 2), false);
  assert.strictEqual(L.isValidSaw([], 2), false);
  assert.strictEqual(L.isValidSaw('abc', 2), false);
  assert.throws(() => L.isValidSaw([[0, 0, 0]], 4), /dimension/);
});

// -----------------------------------------------------------------------
// Contacts and energy
// -----------------------------------------------------------------------
test('straight chain has no non bonded contacts', () => {
  assert.strictEqual(L.contactPairs(L.straightCoords(8), 2).length, 0);
  assert.strictEqual(L.energyOf('HHHHHHHH', L.straightCoords(8), 2), 0);
});

test('a closed 2D square makes exactly one contact', () => {
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  assert.deepStrictEqual(L.contactPairs(sq, 2), [[0, 3]]);
});

test('energyOf counts only H-H non bonded pairs', () => {
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  assert.strictEqual(L.energyOf('HHHH', sq, 2), -1);
  assert.strictEqual(L.energyOf('HPPH', sq, 2), -1);  // pair (0,3) is H-H
  assert.strictEqual(L.energyOf('HPHP', sq, 2), 0);   // pair (0,3) is H-P
  assert.strictEqual(L.energyOf('PPPP', sq, 2), 0);
});

test('energy is always non positive and matches hhPairs', () => {
  const rng = f.makeRng(31);
  for (let i = 0; i < 30; i++) {
    const c = L.randomSaw(11, 2, rng);
    const e = L.energyOf(SEQ2, c, 2);
    assert.ok(e <= 0);
    assert.strictEqual(e, -L.hhPairs(SEQ2, c, 2).length || 0);
  }
});

test('contactPairs never lists sequential neighbors or duplicates', () => {
  const rng = f.makeRng(37);
  for (let k = 0; k < 20; k++) {
    const c = L.randomSaw(9, 3, rng);
    const seen = new Set();
    for (const [i, j] of L.contactPairs(c, 3)) {
      assert.ok(j > i + 1);
      assert.strictEqual(man(c[i], c[j]), 1);
      assert.ok(!seen.has(i + ',' + j));
      seen.add(i + ',' + j);
    }
  }
});

test('metricsOf bundles energy contacts rg and ree', () => {
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const m = L.metricsOf('HHHH', sq, 2);
  assert.strictEqual(m.energy, -1);
  assert.strictEqual(m.hhContacts, 1);
  assert.strictEqual(m.contactCount, 1);
  assert.ok(Math.abs(m.ree - 1) < 1e-9);
  assert.ok(m.rg > 0);
});

test('endToEnd is Euclidean distance between termini', () => {
  const c = [[0, 0, 0], [1, 0, 0], [2, 0, 0]];
  assert.strictEqual(L.endToEnd(c), 2);
  const bent = [[0, 0, 0], [1, 0, 0], [1, 1, 0]];
  assert.ok(Math.abs(L.endToEnd(bent) - Math.SQRT2) < 1e-12);
});

test('evaluateConformation rejects invalid candidates and scores valid ones', () => {
  assert.strictEqual(L.evaluateConformation('HHHH', [[0, 0, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]], 2), null);
  const built = L.evaluateConformation('HHHH', [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], 2);
  assert.strictEqual(built.energy, -1);
  assert.deepStrictEqual(built.pairs, [[0, 3]]);
});

// -----------------------------------------------------------------------
// Moves
// -----------------------------------------------------------------------
test('pivotAt 2D performs a quarter turn about the pivot bead', () => {
  const c = [[0, 0, 0], [1, 0, 0], [2, 0, 0]];
  const r = L.pivotAt(c, 1, 1, 2);
  assert.deepStrictEqual(r[0], [0, 0, 0]);
  assert.deepStrictEqual(r[1], [1, 0, 0]);
  assert.deepStrictEqual(r[2], [1, 1, 0]); // (dx,dy) = (1,0) -> (0,1)
});

test('pivotAt 2D double turn mirrors the tail', () => {
  const c = [[0, 0, 0], [1, 0, 0], [2, 0, 0]];
  const r = L.pivotAt(c, 0, 2, 2);
  assert.deepStrictEqual(r[1], [-1, 0, 0]);
  assert.deepStrictEqual(r[2], [-2, 0, 0]);
});

test('pivotAt 3D keeps the tail on the integer lattice', () => {
  const c = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]];
  for (let t = 1; t <= 3; t++) {
    const r = L.pivotAt(c, 0, t, 3);
    for (const p of r) {
      assert.ok(p.every(Number.isInteger));
    }
  }
  // The tail lies on the rotation axis, so a straight chain is invariant.
  assert.deepStrictEqual(L.pivotAt(c, 0, 2, 3)[3], [3, 0, 0]);
});

test('pivotAt 3D rotates a bent tail a half turn about the bond axis', () => {
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 2, 0]];
  const r = L.pivotAt(c, 0, 2, 3);
  // Bead 1 sits on the rotation axis and stays put; the tail mirrors.
  assert.deepStrictEqual(r, [[0, 0, 0], [1, 0, 0], [1, -1, 0], [1, -2, 0]]);
});

test('pivotAt at the last bead is a safe no op', () => {
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0]];
  assert.deepStrictEqual(L.pivotAt(c, 2, 1, 3), c);
  assert.deepStrictEqual(L.pivotAt(c, -1, 1, 2), c);
});

test('pivot candidates produce valid SAWs on a compact fold', () => {
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const rng = f.makeRng(41);
  let moved = 0;
  for (let i = 0; i < 60; i++) {
    const cand = L.pivotCandidate(c, 2, rng);
    if (L.isValidSaw(cand, 2)) moved++;
  }
  assert.ok(moved > 0, 'no valid pivot candidate in 60 tries');
});

test('corner candidate is null on straight chains', () => {
  for (let i = 0; i < 20; i++) {
    assert.strictEqual(L.cornerCandidate(L.straightCoords(8), 2, f.makeRng(i)), null);
  }
});

test('corner candidate proposes or refuses on a closed square', () => {
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const rng = f.makeRng(43);
  // bead 2 at (1,1): its neighbors (1,0) and (0,1) span a square, fourth
  // corner is (0,0) which is occupied; bead 1 at (1,0) also squares with
  // (0,0) and (1,1) -> fourth corner (0,1) occupied. Proposals that would
  // overlap are refused by the evaluation step, like the fold.js moves.
  const cand = L.cornerCandidate(c, 2, rng);
  if (cand !== null) {
    assert.strictEqual(L.evaluateConformation('HHHH', cand, 2), null);
  }
});

test('end flip only moves a terminal bead', () => {
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const rng = f.makeRng(47);
  for (let i = 0; i < 40; i++) {
    const cand = L.endFlipCandidate(c, 2, rng);
    if (!cand) { continue; }
    let diffs = 0;
    for (let k = 0; k < c.length; k++) {
      if (!(cand[k][0] === c[k][0] && cand[k][1] === c[k][1])) diffs++;
    }
    assert.strictEqual(diffs, 1);
    assert.ok(cand[0][0] !== c[0][0] || cand[0][1] !== c[0][1]
      ? diffs === 1 && man(cand[0], c[1]) === 1
      : man(cand[cand.length - 1], c[c.length - 2]) === 1);
  }
});

test('crankshaft candidate refuses non square loops', () => {
  for (let i = 0; i < 25; i++) {
    assert.strictEqual(L.crankshaftCandidate(L.straightCoords(9), 2, f.makeRng(i + 50)), null);
  }
});

test('moveCandidate dispatches and rejects unknown names', () => {
  const c = L.straightCoords(6);
  const rng = f.makeRng(53);
  assert.ok(L.moveCandidate('pivot', c, 2, rng) instanceof Array);
  assert.throws(() => L.moveCandidate('teleport', c, 2, rng), /unknown lab move/);
});

test('mcStepLab at T=0 never accepts an uphill move', () => {
  const fold = { seq: 'HPPHPHPH', dim: 2, coords: L.straightCoords(8), energy: 0 };
  const rng = f.makeRng(59);
  let prev = fold.energy;
  for (let i = 0; i < 1500; i++) {
    const r = L.mcStepLab(fold, rng, 0);
    assert.ok(!r.accepted || r.dE <= 0);
    assert.ok(fold.energy <= prev);
    prev = fold.energy;
  }
});

test('mcStepLab keeps every accepted state a valid SAW', () => {
  const rng = f.makeRng(61);
  const fold = { seq: SEQ2, dim: 2, coords: L.straightCoords(SEQ2.length), energy: 0 };
  for (let i = 0; i < 2000; i++) {
    L.mcStepLab(fold, rng, 0.9);
    assert.ok(L.isValidSaw(fold.coords, 2), 'invalid SAW at step ' + i);
  }
});

test('mcStepLab rejects negative and non finite temperatures', () => {
  const fold = { seq: 'HPPH', dim: 2, coords: L.straightCoords(4), energy: 0 };
  const rng = f.makeRng(67);
  assert.throws(() => L.mcStepLab(fold, rng, -0.5), RangeError);
  assert.throws(() => L.mcStepLab(fold, rng, NaN), RangeError);
  assert.throws(() => L.mcStepLab(fold, rng, Infinity), RangeError);
});

test('mcStepLab at high T accepts some uphill moves', () => {
  const rng = f.makeRng(71);
  const fold = { seq: SEQ2, dim: 2, coords: L.straightCoords(SEQ2.length), energy: 0 };
  let uphill = 0;
  for (let i = 0; i < 4000; i++) {
    const r = L.mcStepLab(fold, rng, 40);
    if (r.accepted && r.dE > 0) uphill++;
  }
  assert.ok(uphill > 0, 'no uphill acceptance at T=40');
});

// -----------------------------------------------------------------------
// Canonical enumeration
// -----------------------------------------------------------------------
test('countSaws counts canonical walks for tiny chains', () => {
  assert.strictEqual(L.countSaws(1, 2), 1);
  assert.strictEqual(L.countSaws(2, 2), 1);
  assert.strictEqual(L.countSaws(3, 2), 2);
  assert.strictEqual(L.countSaws(4, 2), 6);
  assert.strictEqual(L.countSaws(3, 3), 2);
  assert.strictEqual(L.countSaws(4, 3), 10);
});

test('canonical enumeration visits every conformation exactly once', () => {
  const seen = new Set();
  L.canonicalEnumerate(6, 2, (coords) => {
    const key = coords.map(c => c.join(',')).join('|');
    assert.ok(!seen.has(key), 'duplicate conformation ' + key);
    seen.add(key);
  }, 1e6);
  assert.strictEqual(seen.size, L.countSaws(6, 2));
});

test('canonical walks all start with the fixed first two beads', () => {
  L.canonicalEnumerate(5, 3, (coords) => {
    assert.deepStrictEqual(coords[0], [0, 0, 0]);
    assert.deepStrictEqual(coords[1], [1, 0, 0]);
    assert.ok(
      (coords[2][0] === 2 && coords[2][1] === 0 && coords[2][2] === 0) ||
      (coords[2][0] === 1 && coords[2][1] === 1 && coords[2][2] === 0),
      'second step outside the canonical set'
    );
  }, 1e6);
});

test('canonical enumeration only yields valid SAWs', () => {
  L.canonicalEnumerate(7, 2, (coords) => {
    assert.ok(L.isValidSaw(coords, 2));
  }, 1e6);
});

// -----------------------------------------------------------------------
// Exact ground state enumeration
// -----------------------------------------------------------------------
test('runExact finds the hand verified ground state of a square', () => {
  // HHHH in 2D folds into the unit square: exactly one H-H contact.
  const r = L.runExact({ sequence: 'HHHH', dim: 2 });
  assert.strictEqual(r.best.energy, -1);
  assert.strictEqual(r.truncated, false);
  assert.ok(L.isValidSaw(r.best.coords, 2));
});

test('runExact on an all P chain bottoms at zero', () => {
  const r = L.runExact({ sequence: 'PPPPPP', dim: 2 });
  assert.strictEqual(r.best.energy, 0);
});

test('runExact matches brute force enumeration without symmetry', () => {
  // Independent DFS with no symmetry reduction for n <= 6.
  function brute(seq, dim) {
    const dirs = L.DIRS[dim];
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
  for (const dim of [2, 3]) {
    for (const seq of ['HHHH', 'HPPH', 'HPHPH', 'HHHPPH']) {
      const exact = L.runExact({ sequence: seq, dim: dim });
      assert.strictEqual(exact.best.energy, brute(seq, dim), seq + ' dim ' + dim);
    }
  }
});

test('runExact marks truncation when the state budget runs out', () => {
  const r = L.runExact({ sequence: 'HHHHHHHHHH', dim: 2, maxStates: 5 });
  assert.strictEqual(r.truncated, true);
  assert.ok(r.states >= 5);
});

test('exact ground state is never beaten by the seeded heuristics', () => {
  const cmp = L.compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 19, samples: 60, mcSteps: 500, saSweeps: 25 });
  for (const r of cmp.results) {
    assert.ok(r.best.energy >= cmp.exact.best.energy, r.algorithm + ' beat exact');
    assert.strictEqual(r.optimalityGap, r.best.energy - cmp.exact.best.energy);
  }
});

test('exactAllowed respects the per dimension caps', () => {
  assert.ok(L.exactAllowed(12, 2));
  assert.ok(L.exactAllowed(16, 2));
  assert.ok(!L.exactAllowed(17, 2));
  assert.ok(L.exactAllowed(10, 3));
  assert.ok(!L.exactAllowed(11, 3));
});

// -----------------------------------------------------------------------
// Search algorithms
// -----------------------------------------------------------------------
test('randomSaw returns valid walks of the requested length', () => {
  const rng = f.makeRng(79);
  for (let i = 0; i < 30; i++) {
    const c = L.randomSaw(12, 2, rng);
    assert.strictEqual(c.length, 12);
    assert.ok(L.isValidSaw(c, 2));
  }
  for (let i = 0; i < 30; i++) {
    const c = L.randomSaw(10, 3, rng);
    assert.ok(L.isValidSaw(c, 3));
  }
});

test('runRandomWalk samples and keeps the best conformation', () => {
  const r = L.runRandomWalk({ sequence: SEQ2, dim: 2, samples: 120, seed: 83 });
  assert.strictEqual(r.points.length, 120);
  assert.strictEqual(r.trajectory.length, 120);
  assert.ok(L.isValidSaw(r.best.coords, 2));
  let prev = 0;
  for (const t of r.trajectory) {
    assert.ok(t.energy <= prev + 0 || t.energy <= prev, 'best energy must be monotone');
    prev = t.energy;
  }
});

test('runRandomWalk is deterministic for a fixed seed', () => {
  const a = L.runRandomWalk({ sequence: SEQ2, dim: 2, samples: 60, seed: 89 });
  const b = L.runRandomWalk({ sequence: SEQ2, dim: 2, samples: 60, seed: 89 });
  assert.strictEqual(a.best.energy, b.best.energy);
  assert.deepStrictEqual(a.best.coords, b.best.coords);
});

test('runGreedyPivot descends monotonically and stops at a local optimum', () => {
  const r = L.runGreedyPivot({ sequence: 'HPPHPHPH', dim: 2, seed: 97 });
  assert.ok(L.isValidSaw(r.best.coords, 2));
  // From a seeded random start the descent does real work: two improving
  // pivots take this sequence from 0 to -3, still short of the -5 ground
  // state, which is exactly the teaching point of the comparator.
  assert.ok(r.iterations > 0, 'greedy pivot must do real work, not stall at start');
  assert.ok(r.best.energy < 0);
  for (let i = 1; i < r.trajectory.length; i++) {
    assert.ok(r.trajectory[i].energy < r.trajectory[i - 1].energy,
      'greedy steps must strictly improve');
  }
});

test('runGreedyPivot never worsens its start energy', () => {
  const rng = f.makeRng(401);
  for (let s = 0; s < 12; s++) {
    const r = L.runGreedyPivot({ sequence: 'HPPHPHPH', dim: 2, seed: Math.floor(rng() * 1e6) });
    assert.ok(r.best.energy <= r.trajectory[0].energy,
      'final energy above the start on seed ' + s);
    assert.ok(L.isValidSaw(r.best.coords, 2));
  }
});

test('runGreedyPivot is deterministic and validates its output', () => {
  const a = L.runGreedyPivot({ sequence: 'HPHPHPH', dim: 2, seed: 101 });
  const b = L.runGreedyPivot({ sequence: 'HPHPHPH', dim: 2, seed: 101 });
  assert.strictEqual(a.best.energy, b.best.energy);
  assert.deepStrictEqual(a.best.coords, b.best.coords);
});

test('runMetropolis reports stats and a monotone best trace', () => {
  const r = L.runMetropolis({ sequence: SEQ2, dim: 2, seed: 103, steps: 2500, temperature: 0.8 });
  assert.ok(L.isValidSaw(r.best.coords, 2));
  assert.strictEqual(r.temperature, 0.8);
  assert.ok(r.moved >= r.accepted && r.accepted > 0);
  let prev = r.trajectory[0].energy;
  for (const t of r.trajectory) {
    assert.ok(t.energy <= prev, 'best trace regressed');
    prev = t.energy;
  }
});

test('runMetropolis is deterministic for a fixed seed', () => {
  const a = L.runMetropolis({ sequence: 'HPPHPHPH', dim: 2, seed: 107, steps: 800, temperature: 1 });
  const b = L.runMetropolis({ sequence: 'HPPHPHPH', dim: 2, seed: 107, steps: 800, temperature: 1 });
  assert.strictEqual(a.best.energy, b.best.energy);
  assert.deepStrictEqual(a.best.coords, b.best.coords);
});

test('runMetropolis validates its temperature', () => {
  assert.throws(() => L.runMetropolis({ sequence: 'HPPH', dim: 2, temperature: -1 }), RangeError);
  assert.throws(() => L.runMetropolis({ sequence: 'HPPH', dim: 2, temperature: NaN }), RangeError);
});

test('runAnnealing cools geometrically and keeps a valid best', () => {
  const r = L.runAnnealing({ sequence: SEQ2, dim: 2, seed: 109, sweeps: 60, stepsPerSweep: 20, T0: 1.5, T1: 0.05 });
  assert.ok(L.isValidSaw(r.best.coords, 2));
  assert.strictEqual(r.steps, 60 * 20);
  assert.strictEqual(r.trajectory.length, 61);
  let prev = r.trajectory[0].energy;
  for (const t of r.trajectory) {
    assert.ok(t.energy <= prev);
    prev = t.energy;
  }
});

test('runAnnealing is deterministic and seed sensitive', () => {
  const a = L.runAnnealing({ sequence: 'HPPHPHPH', dim: 2, seed: 113, sweeps: 30 });
  const b = L.runAnnealing({ sequence: 'HPPHPHPH', dim: 2, seed: 113, sweeps: 30 });
  assert.strictEqual(a.best.energy, b.best.energy);
  const c = L.runAnnealing({ sequence: 'HPPHPHPH', dim: 2, seed: 114, sweeps: 30 });
  assert.ok(JSON.stringify(a.best.coords) !== JSON.stringify(c.best.coords) || a.best.energy !== c.best.energy);
});

test('annealing validates its schedule bounds', () => {
  assert.throws(() => L.runAnnealing({ sequence: 'HPPH', dim: 2, T0: 0 }), RangeError);
  assert.throws(() => L.runAnnealing({ sequence: 'HPPH', dim: 2, T0: 1, T1: -0.1 }), RangeError);
});

// -----------------------------------------------------------------------
// Comparator driver
// -----------------------------------------------------------------------
test('compareAlgorithms returns five seeded results plus exact', () => {
  const cmp = L.compareAlgorithms({ sequence: SEQ2, dim: 2, seed: 5 });
  assert.strictEqual(cmp.sequence, SEQ2);
  assert.strictEqual(cmp.dim, 2);
  const names = cmp.results.map(r => r.algorithm);
  assert.deepStrictEqual(names, ['random-walk', 'greedy-pivot', 'metropolis', 'annealing', 'exact']);
  assert.ok(cmp.exact && !cmp.exact.truncated);
});

test('comparator results all carry metrics and nonnegative gaps', () => {
  const cmp = L.compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 13 });
  for (const r of cmp.results) {
    assert.ok(r.metrics.rg >= 0);
    assert.ok(r.metrics.ree >= 1);
    assert.strictEqual(r.metrics.energy, r.best.energy);
    assert.ok(r.optimalityGap >= 0);
    assert.ok(r.ms >= 0);
  }
});

test('comparator marks the exact result and skips it beyond the cap', () => {
  const cmp = L.compareAlgorithms({ sequence: 'H'.repeat(20), dim: 3, seed: 3 });
  assert.strictEqual(cmp.exact, null);
  assert.ok(cmp.exactSkipped);
  assert.strictEqual(cmp.results.length, 4);
  for (const r of cmp.results) assert.strictEqual(r.optimalityGap, null);
});

test('comparator runs 3D chains with valid conformations', () => {
  const cmp = L.compareAlgorithms({ sequence: SEQ3, dim: 3, seed: 23 });
  for (const r of cmp.results) {
    assert.ok(L.isValidSaw(r.best.coords, 3), r.algorithm + ' produced an invalid 3D walk');
  }
  assert.ok(cmp.exact.best.energy <= 0);
});

test('comparator is fully deterministic for a fixed seed', () => {
  const a = L.compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 29 });
  const b = L.compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 29 });
  assert.strictEqual(a.results.length, b.results.length);
  for (let i = 0; i < a.results.length; i++) {
    assert.strictEqual(a.results[i].best.energy, b.results[i].best.energy);
    assert.deepStrictEqual(a.results[i].best.coords, b.results[i].best.coords);
  }
});

test('includeExact false suppresses the exact run', () => {
  const cmp = L.compareAlgorithms({ sequence: 'HPPHPH', dim: 2, seed: 2, includeExact: false });
  assert.strictEqual(cmp.results.length, 4);
  assert.strictEqual(cmp.exact, null);
});

// -----------------------------------------------------------------------
// Funnel sampling
// -----------------------------------------------------------------------
test('sampleFunnel returns the requested cloud with sane extrema', () => {
  const r = L.sampleFunnel(SEQ2, 2, 200, 127);
  assert.strictEqual(r.points.length, 200);
  for (const p of r.points) {
    assert.ok(p.E <= 0);
    assert.ok(p.Rg >= 0);
    assert.ok(p.Ree >= 1);
  }
  assert.ok(r.minEnergy <= r.maxEnergy);
  let recomputeMin = 0;
  for (const p of r.points) if (p.E < recomputeMin) recomputeMin = p.E;
  assert.strictEqual(r.minEnergy, recomputeMin);
});

// -----------------------------------------------------------------------
// Mutation stress scan
// -----------------------------------------------------------------------
test('flipHp swaps H and P', () => {
  assert.strictEqual(L.flipHp('H'), 'P');
  assert.strictEqual(L.flipHp('P'), 'H');
});

test('mutationScan reports one entry per residue with correct flips', () => {
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const scan = L.mutationScan('HHHH', sq, 2);
  assert.strictEqual(scan.perResidue.length, 4);
  for (const e of scan.perResidue) {
    assert.strictEqual(e.from, 'H');
    assert.strictEqual(e.to, 'P');
  }
});

test('mutation deltaE counts lost and gained H-H contacts', () => {
  // HHHH square: energy -1 via the (0,3) contact. Mutating residue 0 or 3
  // to P destroys it (+1 destabilizing); residues 1 and 2 change nothing.
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const scan = L.mutationScan('HHHH', sq, 2);
  assert.strictEqual(scan.e0, -1);
  assert.strictEqual(scan.perResidue[0].dE, 1);
  assert.strictEqual(scan.perResidue[3].dE, 1);
  assert.strictEqual(scan.perResidue[1].dE, 0);
  assert.strictEqual(scan.perResidue[2].dE, 0);
  assert.strictEqual(scan.destabilizing, 2);
  assert.strictEqual(scan.neutral, 2);
  assert.strictEqual(scan.worst.dE, 1);
});

test('mutation categories sum to the chain length', () => {
  const rng = f.makeRng(131);
  const c = L.randomSaw(10, 2, rng);
  const scan = L.mutationScan('HPPHPHPHPH', c, 2);
  assert.strictEqual(scan.destabilizing + scan.stabilizing + scan.neutral, 10);
});

test('P to H mutations on a contact free chain never destabilize', () => {
  const c = L.straightCoords(8);
  const scan = L.mutationScan('PPPPPPPP', c, 2);
  assert.strictEqual(scan.e0, 0);
  assert.strictEqual(scan.destabilizing, 0);
});

test('mutationScan validates the HP sequence', () => {
  assert.throws(() => L.mutationScan('HXH', L.straightCoords(3), 2), /only H and P/);
});

// -----------------------------------------------------------------------
// Contact matrix
// -----------------------------------------------------------------------
test('contactMatrix codes bonds, contacts and the diagonal', () => {
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const m = L.contactMatrix('HHHH', sq, 2);
  assert.strictEqual(m.n, 4);
  assert.strictEqual(m.cells[0 * 4 + 0], L.CMAP_SELF);
  assert.strictEqual(m.cells[0 * 4 + 1], L.CMAP_BONDED);
  assert.strictEqual(m.cells[0 * 4 + 3], L.CMAP_HH);
  assert.strictEqual(m.cells[1 * 4 + 3], L.CMAP_FAR);
  assert.strictEqual(m.hhCount, 1);
  assert.strictEqual(m.contactCount, 1);
});

test('contactMatrix is symmetric with a full diagonal', () => {
  const rng = f.makeRng(137);
  const c = L.randomSaw(11, 2, rng);
  const m = L.contactMatrix(SEQ2, c, 2);
  for (let a = 0; a < m.n; a++) {
    for (let b = 0; b < m.n; b++) {
      assert.strictEqual(m.cells[a * m.n + b], m.cells[b * m.n + a]);
    }
    assert.strictEqual(m.cells[a * m.n + a], L.CMAP_SELF);
  }
});

test('contactMatrix H-H count equals minus the energy', () => {
  const rng = f.makeRng(139);
  for (let i = 0; i < 15; i++) {
    const c = L.randomSaw(9, 3, rng);
    const m = L.contactMatrix('HHPHPHPHH', c, 3);
    assert.strictEqual(m.hhCount, Math.abs(L.energyOf('HHPHPHPHH', c, 3)));
    let hh = 0;
    for (const v of m.cells) if (v === L.CMAP_HH) hh++;
    assert.strictEqual(hh / 2, m.hhCount);
  }
});

test('contactMatrix distinguishes H-H from other contacts', () => {
  const sq = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  const m = L.contactMatrix('HPHP', sq, 2); // (0,3) is H-P
  assert.strictEqual(m.cells[0 * 4 + 3], L.CMAP_OTHER);
  assert.strictEqual(m.hhCount, 0);
  assert.strictEqual(m.contactCount, 1);
});

// -----------------------------------------------------------------------
// SHA-256 (cross validated against the platform crypto)
// -----------------------------------------------------------------------
test('sha256Hex matches the FIPS 180-4 empty string vector', () => {
  assert.strictEqual(L.sha256Hex(''), 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855');
});

test('sha256Hex matches the FIPS 180-4 abc vector', () => {
  assert.strictEqual(L.sha256Hex('abc'), 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
});

test('sha256Hex matches the 448 bit two block vector', () => {
  assert.strictEqual(
    L.sha256Hex('abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'),
    '248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1');
});

test('sha256Hex matches the quick brown fox vector', () => {
  assert.strictEqual(L.sha256Hex('The quick brown fox jumps over the lazy dog'),
    'd7a8fbb307d7809469ca9abcb0082e4f8d5651e46d3cdb762d02d0bf37c9e592');
});

test('sha256Hex matches the one million a vector', () => {
  assert.strictEqual(L.sha256Hex('a'.repeat(1000000)),
    'cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0');
});

test('sha256Hex agrees with platform crypto on arbitrary strings', () => {
  const rng = f.makeRng(149);
  const alphabet = 'HPXYZ ,:{}[]=#@!$%^&*()0123456789';
  for (let i = 0; i < 40; i++) {
    let s = '';
    const len = 1 + Math.floor(rng() * 200);
    for (let k = 0; k < len; k++) s += alphabet[Math.floor(rng() * alphabet.length)];
    assert.strictEqual(L.sha256Hex(s), shaRef(s), 'mismatch on ' + JSON.stringify(s.slice(0, 30)));
  }
});

test('sha256Hex agrees with platform crypto on multibyte input', () => {
  for (const s of ['é', 'español', '日本語', 'emoji 🙂 fold', 'Ω≈ç√∫', 'a̐éö̲۬']) {
    assert.strictEqual(L.sha256Hex(s), shaRef(s), 'mismatch on ' + s);
  }
});

test('utf8Bytes encodes ASCII, two, three and four byte sequences', () => {
  assert.deepStrictEqual(L.utf8Bytes('ABC'), [65, 66, 67]);
  assert.deepStrictEqual(L.utf8Bytes('é'), [0xC3, 0xA9]);
  assert.deepStrictEqual(L.utf8Bytes('€'), [0xE2, 0x82, 0xAC]);
  assert.deepStrictEqual(L.utf8Bytes('🙂'), [0xF0, 0x9F, 0x99, 0x82]);
  assert.deepStrictEqual(L.utf8Bytes('Aé€🙂'), Buffer.from('Aé€🙂', 'utf8').toJSON().data);
});

// -----------------------------------------------------------------------
// Canonical JSON and checksummed reports
// -----------------------------------------------------------------------
test('canonicalize is key order independent', () => {
  const a = L.canonicalize({ b: 1, a: [3, { y: 2, x: 1 }] });
  const b = L.canonicalize({ a: [3, { x: 1, y: 2 }], b: 1 });
  assert.strictEqual(a, b);
  assert.strictEqual(a, '{"a":[3,{"x":1,"y":2}],"b":1}');
});

test('canonicalize handles scalars arrays and nulls', () => {
  assert.strictEqual(L.canonicalize(5), '5');
  assert.strictEqual(L.canonicalize('x'), '"x"');
  assert.strictEqual(L.canonicalize(null), 'null');
  assert.strictEqual(L.canonicalize([1, 'a', null]), '[1,"a",null]');
  assert.strictEqual(L.canonicalize({}), '{}');
});

test('buildAuditReport stamps a 64 hex checksum and metadata', () => {
  const r = L.buildAuditReport({ sequence: 'HPPH', dimension: 2 });
  assert.match(r.checksum, /^[0-9a-f]{64}$/);
  assert.strictEqual(r.tool, 'fold-lattice');
  assert.strictEqual(r.toolVersion, f.VERSION);
  assert.strictEqual(r.kind, 'hp-folding-audit');
  assert.match(r.checksumAlgorithm, /SHA-256/);
});

test('audit checksum is deterministic and payload sensitive', () => {
  const a = L.buildAuditReport({ sequence: 'HPPH', x: 1 });
  const b = L.buildAuditReport({ sequence: 'HPPH', x: 1 });
  const c = L.buildAuditReport({ sequence: 'HPPH', x: 2 });
  assert.strictEqual(a.checksum, b.checksum);
  assert.notStrictEqual(a.checksum, c.checksum);
});

test('verifyAuditReport accepts a fresh report and rejects edits', () => {
  const r = L.buildAuditReport({ sequence: 'HPPH', results: [{ algorithm: 'x', energy: -2 }] });
  assert.ok(L.verifyAuditReport(r).ok);
  r.payload.sequence = 'HPHP';
  assert.strictEqual(L.verifyAuditReport(r).ok, false);
});

test('verifyAuditReport rejects malformed reports', () => {
  assert.strictEqual(L.verifyAuditReport(null).ok, false);
  assert.strictEqual(L.verifyAuditReport({}).ok, false);
  assert.strictEqual(L.verifyAuditReport({ checksum: 'zzz', payload: {} }).ok, false);
});

test('verifyAuditText parses, verifies and reports bad JSON', () => {
  const r = L.buildAuditReport({ a: 1 });
  const text = L.auditToJson(r);
  assert.ok(L.verifyAuditText(text).ok);
  const bad = L.verifyAuditText('{not json');
  assert.strictEqual(bad.ok, false);
  assert.match(bad.reason, /not valid JSON/);
});

test('tamperAuditText flips one checksum digit', () => {
  const r = L.buildAuditReport({ a: 1 });
  const t = L.tamperAuditText(L.auditToJson(r));
  assert.ok(!L.verifyAuditText(t).ok);
  assert.strictEqual(L.tamperAuditText('no hex here'), 'no hex here');
});

test('audit JSON round trips through verifyAuditText', () => {
  const payload = { sequence: SEQ2, dimension: 2, results: [{ algorithm: 'exact', energy: -5 }] };
  const text = L.auditToJson(L.buildAuditReport(payload));
  const v = L.verifyAuditText(text);
  assert.ok(v.ok);
});

// -----------------------------------------------------------------------
// CSV export and formula injection escaping
// -----------------------------------------------------------------------
test('csvCell escapes commas quotes and newlines', () => {
  assert.strictEqual(L.csvCell('a,b'), '"a,b"');
  assert.strictEqual(L.csvCell('say "hi"'), '"say ""hi"""');
  assert.strictEqual(L.csvCell('two\nlines'), '"two\nlines"');
  assert.strictEqual(L.csvCell('plain'), 'plain');
  assert.strictEqual(L.csvCell(42), '42');
  assert.strictEqual(L.csvCell(null), '');
});

test('csvCell neutralizes spreadsheet formula prefixes', () => {
  for (const bad of ['=1+1', '+cmd', '-1+1', '-x', '@SUM(A1)', '\t=x', ' =x']) {
    const cell = L.csvCell(bad);
    assert.ok(cell.startsWith('"\''), bad + ' not neutralized: ' + cell);
  }
  // Plain numbers are values, not formulas, and stay untouched.
  assert.strictEqual(L.csvCell('-10'), '-10');
  assert.strictEqual(L.csvCell('+3.5'), '+3.5');
  assert.strictEqual(L.csvCell("=HYPERLINK(\"http://x\")"), '"\'=HYPERLINK(""http://x"")"');
});

test('auditToCsv emits meta algorithm and mutation sections', () => {
  const report = L.buildAuditReport({
    sequence: SEQ2, dimension: 2, seed: 7, groundEnergy: -5,
    results: [{ algorithm: 'exact', energy: -5, contactCount: 7, hhContacts: 5, rg: 1.4, ree: 2, optimalityGap: 0, ms: 3 }],
    mutationScan: { perResidue: [{ index: 0, from: 'H', to: 'P', dE: 1 }] }
  });
  const csv = L.auditToCsv(report);
  assert.match(csv, /section,field,value/);
  assert.match(csv, /algorithm,best_energy,contact_count/);
  assert.match(csv, /exact,-5,7,5,1.4,2,0,3/);
  assert.match(csv, /mutation_index,from,to,delta_e/);
  assert.match(csv, /0,H,P,1/);
});

test('auditToCsv escapes hostile sequence names', () => {
  const report = L.buildAuditReport({ sequence: '=HPPH+1', dimension: 2, results: [] });
  const csv = L.auditToCsv(report);
  assert.ok(csv.includes("'=HPPH+1"), 'formula sequence not escaped');
  assert.ok(!/\n=HPPH/.test(csv), 'raw formula cell on a line');
});

// -----------------------------------------------------------------------
// Markdown export
// -----------------------------------------------------------------------
test('auditToMarkdown produces a table and the disclosure', () => {
  const report = L.buildAuditReport({
    sequence: SEQ2, dimension: 2, seed: 7, groundEnergy: -5,
    results: [{ algorithm: 'annealing', energy: -5, contactCount: 7, hhContacts: 5, rg: 1.4, ree: 2, optimalityGap: 0, ms: 3 }],
    mutationScan: { destabilizing: 2, stabilizing: 0, neutral: 9, worst: { index: 3, from: 'H', to: 'P', dE: 2 }, perResidue: [] }
  });
  const md = L.auditToMarkdown(report);
  assert.match(md, /# fold-lattice folding audit/);
  assert.match(md, /\| algorithm \| best E \|/);
  assert.match(md, /\| annealing \| -5 \|/);
  assert.match(md, /Dill HP lattice model/);
  assert.match(md, new RegExp(report.checksum));
});

// -----------------------------------------------------------------------
// URL hash bookmarks
// -----------------------------------------------------------------------
test('encodeBookmark formats fold=SEQ:DIMd:SEED', () => {
  assert.strictEqual(L.encodeBookmark({ sequence: 'hpph', dim: 2, seed: 42 }), 'fold=HPPH:2d:42');
  assert.strictEqual(L.encodeBookmark({ sequence: 'HPPH', dim: 3, seed: 7 }), 'fold=HPPH:3d:7');
});

test('decodeBookmark parses with or without the hash mark', () => {
  assert.deepStrictEqual(L.decodeBookmark('#fold=HPPH:2d:42'), { sequence: 'HPPH', dim: 2, seed: 42 });
  assert.deepStrictEqual(L.decodeBookmark('fold=hpph:3d:9'), { sequence: 'HPPH', dim: 3, seed: 9 });
});

test('decodeBookmark rejects malformed hashes', () => {
  assert.strictEqual(L.decodeBookmark('#fold=HPXH:2d:1'), null);
  assert.strictEqual(L.decodeBookmark('#fold=HPPH:4d:1'), null);
  assert.strictEqual(L.decodeBookmark('#fold=HPPH:2d'), null);
  assert.strictEqual(L.decodeBookmark('#other=HPPH:2d:1'), null);
  assert.strictEqual(L.decodeBookmark(''), null);
  assert.strictEqual(L.decodeBookmark(null), null);
  assert.strictEqual(L.decodeBookmark('#fold=HPPH:2d:99999999999'), null);
});

test('bookmark round trips through encode then decode', () => {
  const mark = L.encodeBookmark({ sequence: 'HHPPHPH', dim: 3, seed: 1234 });
  const back = L.decodeBookmark('#' + mark);
  assert.deepStrictEqual(back, { sequence: 'HHPPHPH', dim: 3, seed: 1234 });
});

test('encodeBookmark validates its inputs', () => {
  assert.throws(() => L.encodeBookmark({ sequence: 'HX', dim: 2, seed: 1 }), /only H and P|at least/);
  assert.throws(() => L.encodeBookmark({ sequence: 'HPPH', dim: 4, seed: 1 }), /dimension/);
  assert.throws(() => L.encodeBookmark({ sequence: 'HPPH', dim: 2, seed: -3 }), /seed/);
  assert.throws(() => L.encodeBookmark({ sequence: 'HPPH', dim: 2, seed: 4294967296 }), /seed/);
  assert.throws(() => L.encodeBookmark({ sequence: 'H'.repeat(121), dim: 2, seed: 1 }), /too long/);
});

test('decodeBookmark rejects overlong sequences', () => {
  assert.strictEqual(L.decodeBookmark('#fold=' + 'H'.repeat(121) + ':2d:1'), null);
  assert.ok(L.decodeBookmark('#fold=' + 'H'.repeat(120) + ':2d:1') !== null);
});

// -----------------------------------------------------------------------
// HP to AA bridge and self tests
// -----------------------------------------------------------------------
test('hpToAa maps H to hydrophobic A and P to polar S', () => {
  assert.strictEqual(L.hpToAa('HPPHPH'), 'ASSASA');
  for (const ch of L.hpToAa('HHPPHPHH')) {
    assert.ok(ch === 'A' || ch === 'S');
  }
});

test('hpToAa bridge lands inside the engine hydrophobic set', () => {
  const aa = L.hpToAa('HPPHPH');
  assert.ok(f.isHydrophobic(aa[0]));
  assert.ok(!f.isHydrophobic(aa[1]));
});

test('lab energy agrees with the LatticeFolder HP energy on the same fold', () => {
  // Cross engine check: a lab conformation scored identically when the
  // sequence is mapped into A/S letters through the hpToAa bridge.
  const rng = f.makeRng(151);
  for (let i = 0; i < 12; i++) {
    const coords = L.randomSaw(9, 3, rng);
    const labE = L.energyOf('HHPHPHPHH', coords, 3);
    const lf = new f.LatticeFolder(L.hpToAa('HHPHPHPHH'), { mode: 'hp' });
    const occ = new Map();
    for (let k = 0; k < coords.length; k++) {
      occ.set(f.LatticeFolder.key(coords[k][0], coords[k][1], coords[k][2]), k);
    }
    lf._commit(coords, occ, null, null);
    assert.strictEqual(lf.energy(), labE);
  }
});

test('in page self tests all pass', () => {
  const out = L.runSelfTests();
  assert.strictEqual(out.passed, out.total, out.results.filter(r => !r.pass).map(r => r.name + ': ' + r.detail).join('; '));
  assert.ok(out.total >= 12);
});

test('LAB_VERSION is exposed', () => {
  assert.strictEqual(L.LAB_VERSION, '2.0.0');
});
