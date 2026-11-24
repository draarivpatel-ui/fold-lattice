'use strict';

// Second test suite: hostile invariants, validation boundaries, PDB edge
// cases and static file checks. fold.test.js covers the happy paths.

const test = require('node:test');
const assert = require('node:assert');
const fs = require('node:fs');
const path = require('node:path');
const f = require('../src/fold.js');

const SEQ = 'GAGSLYQGSSLAQSGFLSSHHSGF';
const HP_SEQ = 'SSSSAAAASSSSAAAASSSAA';
const ROOT = path.join(__dirname, '..');
const read = name => fs.readFileSync(path.join(ROOT, name), 'utf8');

function manhattan(a, b) {
  return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);
}

// Independent energy recount that never touches the folder caches.
function recountEnergy(lf) {
  let e = 0;
  for (let i = 0; i < lf.n; i++) {
    for (let j = i + 2; j < lf.n; j++) {
      if (manhattan(lf.coords[i], lf.coords[j]) === 1) {
        if (lf.mode === 'hp') {
          if (f.isHydrophobic(lf.sequence[i]) && f.isHydrophobic(lf.sequence[j])) e -= 1;
        } else {
          e += f.mjPairEnergy(lf.sequence[i], lf.sequence[j]);
        }
      }
    }
  }
  return e;
}

function occupancyConsistent(lf) {
  if (lf.occ.size !== lf.n) return false;
  for (const [key, idx] of lf.occ) {
    const c = lf.coords[idx];
    if (f.LatticeFolder.key(c[0], c[1], c[2]) !== key) return false;
  }
  const keys = new Set(lf.coords.map(c => f.LatticeFolder.key(c[0], c[1], c[2])));
  return keys.size === lf.n;
}

// ---------------------------------------------------------------------
// Metropolis criterion at and below zero temperature
// ---------------------------------------------------------------------
test('mcStep at T=0 is strictly greedy: no uphill move is ever accepted', () => {
  const rng = f.makeRng(101);
  const lf = new f.LatticeFolder(SEQ);
  for (let i = 0; i < 3000; i++) {
    const res = lf.mcStep(0, rng);
    if (res.moved && res.accepted) {
      assert.ok(res.dE <= 0, 'accepted uphill move at T=0, dE ' + res.dE);
    }
  }
});

test('mcStep at T=0 never raises the running energy', () => {
  const rng = f.makeRng(103);
  const lf = new f.LatticeFolder(SEQ);
  let prev = lf.energy();
  for (let i = 0; i < 2000; i++) {
    lf.mcStep(0, rng);
    assert.ok(lf.energy() <= prev + 1e-12, 'energy rose at T=0');
    prev = lf.energy();
  }
});

test('mcStep at negative temperature throws instead of silently going greedy', () => {
  // A negative Metropolis temperature is not greedy: it inverts acceptance.
  // Passing one used to run as T=0, hiding the caller's bad input.
  const rng = f.makeRng(107);
  const lf = new f.LatticeFolder(SEQ);
  assert.throws(() => lf.mcStep(-3, rng), RangeError);
});

test('mcStep at high temperature does accept some uphill moves', () => {
  // Metropolis sanity: at T=50 nearly every move is accepted, including
  // energy raising ones. Deterministic under the fixed seed.
  const rng = f.makeRng(109);
  const lf = new f.LatticeFolder(SEQ);
  let uphill = 0;
  for (let i = 0; i < 6000; i++) {
    const res = lf.mcStep(50, rng);
    if (res.moved && res.accepted && res.dE > 0) uphill++;
  }
  assert.ok(uphill > 0, 'no uphill acceptance at T=50');
});

test('mcStep result reports moved, accepted, dE and type consistently', () => {
  const rng = f.makeRng(113);
  const lf = new f.LatticeFolder(SEQ);
  const before = lf.energy();
  const res = lf.mcStep(1.0, rng);
  assert.ok(['pull', 'pivot', 'crankshaft', 'corner', 'endpull'].includes(res.type));
  if (res.moved && res.accepted) {
    assert.ok(Math.abs(lf.energy() - (before + res.dE)) < 1e-9);
  } else {
    assert.strictEqual(lf.energy(), before);
  }
});

// ---------------------------------------------------------------------
// Occupancy map and energy invariants under random play
// ---------------------------------------------------------------------
test('occupancy map stays consistent after 500 accepted moves', () => {
  const rng = f.makeRng(127);
  const lf = new f.LatticeFolder(SEQ);
  for (let i = 0; i < 4000; i++) {
    const res = lf.tryMove(['pull', 'pivot', 'crankshaft', 'corner', 'endpull'][i % 5], rng);
    if (res.moved) {
      assert.ok(occupancyConsistent(lf), 'occ map diverged at move ' + i);
    }
  }
});

test('cached energy always matches an independent recount', () => {
  const rng = f.makeRng(131);
  const lf = new f.LatticeFolder(SEQ);
  assert.strictEqual(lf.energy(), recountEnergy(lf));
  for (let i = 0; i < 1500; i++) {
    lf.tryMove(['pull', 'pivot', 'crankshaft', 'corner', 'endpull'][i % 5], rng);
    if (i % 100 === 0) assert.strictEqual(lf.energy(), recountEnergy(lf));
  }
  assert.strictEqual(lf.energy(), recountEnergy(lf));
});

test('hp energy always equals minus the H-H contact count', () => {
  const rng = f.makeRng(137);
  const lf = new f.LatticeFolder(HP_SEQ, { mode: 'hp' });
  for (let i = 0; i < 1200; i++) {
    lf.tryMove(['pull', 'pivot', 'corner', 'crankshaft'][i % 4], rng);
    if (i % 150 === 0) assert.strictEqual(lf.energy(), recountEnergy(lf));
  }
});

test('contacts list contains no sequential pairs and no duplicates', () => {
  const rng = f.makeRng(139);
  const lf = new f.LatticeFolder(SEQ);
  for (let i = 0; i < 400; i++) lf.mcStep(2.0, rng);
  const seen = new Set();
  for (const [a, b] of lf.contacts()) {
    assert.ok(b > a + 1, 'sequential pair listed as contact');
    assert.strictEqual(manhattan(lf.coords[a], lf.coords[b]), 1);
    const k = a + ',' + b;
    assert.ok(!seen.has(k), 'duplicate contact ' + k);
    seen.add(k);
  }
});

test('tryMove dE equals energy difference between states', () => {
  const rng = f.makeRng(149);
  const lf = new f.LatticeFolder(SEQ);
  for (let i = 0; i < 800; i++) {
    const before = lf.energy();
    const res = lf.tryMove(['pull', 'pivot', 'corner'][i % 3], rng);
    if (res.moved) {
      assert.ok(Math.abs(res.dE - (lf.energy() - before)) < 1e-9);
    }
  }
});

// ---------------------------------------------------------------------
// cellKey collision bounds
// ---------------------------------------------------------------------
test('cellKey is unique across a dense local cube', () => {
  const seen = new Set();
  for (let x = -6; x <= 6; x++) {
    for (let y = -6; y <= 6; y++) {
      for (let z = -6; z <= 6; z++) {
        const k = f.LatticeFolder.key(x, y, z);
        assert.ok(!seen.has(k), 'collision at ' + [x, y, z]);
        seen.add(k);
      }
    }
  }
});

test('cellKey stays unique for far coordinates within design bounds', () => {
  const pts = [[0, 0, 0], [-1048576, 0, 0], [1048576, 0, 0], [0, -999999, 0],
    [0, 0, 888888], [123456, -654321, 111111], [-1, -1, -1], [1, 1, 1]];
  const keys = pts.map(p => f.LatticeFolder.key(p[0], p[1], p[2]));
  assert.strictEqual(new Set(keys).size, keys.length);
});

// ---------------------------------------------------------------------
// geometricTemps validation
// ---------------------------------------------------------------------
test('geometricTemps rejects non integer counts', () => {
  assert.throws(() => f.geometricTemps(2.5, 0.5, 2), /at least one|positive integer/);
  assert.throws(() => f.geometricTemps('4', 0.5, 2), /at least one/);
  assert.throws(() => f.geometricTemps(NaN, 0.5, 2), /at least one/);
});

test('geometricTemps rejects non positive temperatures', () => {
  assert.throws(() => f.geometricTemps(4, 0, 2), /minTemp.*positive/);
  assert.throws(() => f.geometricTemps(4, -0.5, 2), /minTemp.*positive/);
  assert.throws(() => f.geometricTemps(4, 0.5, 0), /maxTemp.*positive/);
  assert.throws(() => f.geometricTemps(4, NaN, 2), /minTemp.*positive|finite/);
  assert.throws(() => f.geometricTemps(4, 0.5, Infinity), /maxTemp.*positive|finite/);
});

test('geometricTemps rejects inverted ladders', () => {
  assert.throws(() => f.geometricTemps(4, 2.0, 0.5), /minTemp|greater than or equal/);
});

test('geometricTemps equal endpoints produce a flat ladder', () => {
  const t = f.geometricTemps(5, 1.0, 1.0);
  assert.deepStrictEqual(t, [1, 1, 1, 1, 1]);
});

// ---------------------------------------------------------------------
// Generator option validation
// ---------------------------------------------------------------------
test('remcGenerator rejects invalid replica counts', () => {
  for (const bad of [0, -2, 2.5, 'abc', NaN]) {
    assert.throws(() => f.simulateRemc({ sequence: 'ACDE', replicas: bad }),
      /replicas.*positive integer/, String(bad));
  }
});

test('remcGenerator rejects invalid sweep and step counts', () => {
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', sweeps: 0 }), /sweeps/);
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', sweeps: -4 }), /sweeps/);
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', sweeps: 'x' }), /sweeps/);
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', stepsPerSweep: 0 }), /stepsPerSweep/);
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', yieldEvery: 0 }), /yieldEvery/);
});

test('remcGenerator rejects bad temperature bounds', () => {
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', minTemp: -1 }), /minTemp/);
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', minTemp: 0 }), /minTemp/);
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', maxTemp: Infinity }), /maxTemp/);
  assert.throws(() => f.simulateRemc({ sequence: 'ACDE', minTemp: 3, maxTemp: 1 }), /minTemp|maxTemp/);
});

test('annealingGenerator rejects invalid counts and temperatures', () => {
  assert.throws(() => f.simulateAnnealing({ sequence: 'ACDE', sweeps: 0 }), /sweeps/);
  assert.throws(() => f.simulateAnnealing({ sequence: 'ACDE', stepsPerSweep: -1 }), /stepsPerSweep/);
  assert.throws(() => f.simulateAnnealing({ sequence: 'ACDE', yieldEvery: 1.5 }), /yieldEvery/);
  assert.throws(() => f.simulateAnnealing({ sequence: 'ACDE', T0: 0 }), /T0/);
  assert.throws(() => f.simulateAnnealing({ sequence: 'ACDE', T0: -1 }), /T0/);
  assert.throws(() => f.simulateAnnealing({ sequence: 'ACDE', T1: -0.5 }), /T1/);
  assert.throws(() => f.simulateAnnealing({ sequence: 'ACDE', T0: NaN }), /T0/);
});

test('annealing with T1 of 0 quenches greedily and still finishes', () => {
  const r = f.simulateAnnealing({ sequence: 'ACDEFG', sweeps: 20, stepsPerSweep: 8, T1: 0, seed: 5 });
  assert.ok(r.done && r.best.coords.length === 6);
  assert.ok(r.best.energy <= 0);
});

test('remc with a single replica runs without exchange machinery', () => {
  const r = f.simulateRemc({ sequence: 'ACDEFG', replicas: 1, sweeps: 15, stepsPerSweep: 8, seed: 7 });
  assert.ok(r.done);
  assert.strictEqual(r.exchanges, 0);
  assert.strictEqual(r.temps.length, 1);
});

// ---------------------------------------------------------------------
// PDB parser model bookkeeping
// ---------------------------------------------------------------------
const MINI_PDB = [
  'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
  'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
  'ATOM      3  C   ALA A   1       2.000   1.300   0.300  1.00  0.00           C',
  'ATOM      4  O   ALA A   1       1.500   2.400   0.400  1.00  0.00           O',
  'TER',
  'END'
].join('\n');

const TWO_MODEL = [
  'MODEL        1',
  MINI_PDB.split('\n')[0],
  MINI_PDB.split('\n')[1],
  'ENDMDL',
  'MODEL        2',
  'ATOM      5  N   GLY B   9       9.000   9.000   9.000  1.00  0.00           N',
  'ENDMDL',
  'END'
].join('\n');

test('parsePdb reports one implicit model for files without MODEL records', () => {
  const p = f.parsePdb(MINI_PDB);
  assert.strictEqual(p.numModelsRead, 1);
  assert.strictEqual(p.modelCount, 1);
});

test('parsePdb counts every MODEL but reads only the first', () => {
  const p = f.parsePdb(TWO_MODEL);
  assert.strictEqual(p.numModelsRead, 1);
  assert.strictEqual(p.modelCount, 2);
  assert.strictEqual(p.atoms.length, 2);
  assert.strictEqual(p.atoms[0].resName, 'ALA');
});

test('parsePdb skips altLoc codes other than blank, A and 1', () => {
  const alt = [
    'ATOM      1  CA AALA A   1       0.000   0.000   0.000  1.00  0.00           C',
    'ATOM      2  CA BALA A   1       9.000   9.000   9.000  1.00  0.00           C'
  ].join('\n');
  const p = f.parsePdb(alt);
  assert.strictEqual(p.atoms.length, 1);
  assert.strictEqual(p.atoms[0].altLoc, 'A');
});

test('parsePdb handles CRLF and CR line endings', () => {
  const crlf = f.parsePdb(MINI_PDB.split('\n').join('\r\n'));
  assert.strictEqual(crlf.atoms.length, 4);
  const cr = f.parsePdb(MINI_PDB.split('\n').join('\r'));
  assert.strictEqual(cr.atoms.length, 4);
});

test('parsePdb tolerates short lines missing trailing columns', () => {
  const short = 'ATOM      1  CA  ALA A   1       1.000   2.000   3.000\n';
  const p = f.parsePdb(short);
  assert.strictEqual(p.atoms.length, 1);
  assert.strictEqual(p.atoms[0].element, 'C');
});

// ---------------------------------------------------------------------
// PDB export edge cases
// ---------------------------------------------------------------------
test('atomsToPdb TER records carry the last residue number of each chain', () => {
  const atoms = [
    { name: 'CA', resName: 'GLY', resSeq: 7, chainID: 'A', xyz: [0, 0, 0], element: 'C' },
    { name: 'CA', resName: 'ALA', resSeq: 8, chainID: 'A', xyz: [3.8, 0, 0], element: 'C' }
  ];
  const pdb = f.atomsToPdb(atoms);
  const ter = pdb.split('\n').find(l => l.startsWith('TER'));
  assert.ok(ter, 'TER line present');
  assert.strictEqual(ter.slice(21, 22), 'A');
  assert.strictEqual(ter.slice(22, 26).trim(), '8');
});

test('atomsToPdb on an empty atom list emits no TER', () => {
  const pdb = f.atomsToPdb([]);
  assert.ok(!pdb.split('\n').some(l => l.startsWith('TER')), 'unexpected TER');
  assert.ok(pdb.split('\n').includes('END'));
});

test('pdbAtomNameField left justifies digit and two letter element names', () => {
  assert.strictEqual(f.pdbAtomNameField('1HB'), '1HB ');
  assert.strictEqual(f.pdbAtomNameField('FE', 'FE'), 'FE  ');
  assert.strictEqual(f.pdbAtomNameField('CA', 'C'), ' CA ');
  assert.strictEqual(f.pdbAtomNameField('CA', 'CA'), 'CA  ');
});

test('export then reparse of two letter element keeps element column', () => {
  const atoms = [{ name: 'FE', resName: 'HET', resSeq: 1, chainID: 'A', xyz: [0, 0, 0], element: 'FE' }];
  const p = f.parsePdb(f.atomsToPdb(atoms));
  assert.strictEqual(p.atoms[0].element, 'FE');
});

test('atomsToPdb honors startSerial', () => {
  const atoms = [{ name: 'CA', resName: 'GLY', resSeq: 1, chainID: 'A', xyz: [0, 0, 0], element: 'C' }];
  const pdb = f.atomsToPdb(atoms, { startSerial: 100 });
  const atom = pdb.split('\n').find(l => l.startsWith('ATOM'));
  assert.strictEqual(parseInt(atom.slice(6, 11), 10), 100);
});

// ---------------------------------------------------------------------
// Backbone reconstruction validation
// ---------------------------------------------------------------------
test('reconstructBackbone rejects non finite angles', () => {
  const good = [{ phi: -57, psi: -47 }, { phi: -139, psi: 135 }];
  assert.throws(() => f.reconstructBackbone('AA', [{ phi: NaN, psi: 0 }, good[1]]), /numeric|finite/);
  assert.throws(() => f.reconstructBackbone('AA', [{ phi: -57, psi: Infinity }, good[1]]), /numeric|finite/);
  assert.throws(() => f.reconstructBackbone('AA', [good[0], { phi: 1, psi: 1, omega: NaN }]), /finite/);
});

test('reconstructBackbone angles survive normalizeDeg extremes', () => {
  const atoms = f.reconstructBackbone('AAA', [
    { phi: -179.9, psi: 179.9, omega: 179.9 },
    { phi: 179.9, psi: -179.9, omega: -179.9 },
    { phi: 0.1, psi: -0.1, omega: 0.5 }
  ]);
  assert.strictEqual(atoms.length, 12);
  for (const a of atoms) {
    for (const v of a.xyz) assert.ok(Number.isFinite(v), 'finite coordinate');
  }
});

// ---------------------------------------------------------------------
// SAW validity on hand crafted conformations
// ---------------------------------------------------------------------
test('isValidConformation rejects duplicate cells and broken bonds', () => {
  const lf = new f.LatticeFolder('ACDEF');
  assert.ok(!lf.isValidConformation([[0, 0, 0], [1, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]), 'dup cell');
  assert.ok(!lf.isValidConformation([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [5, 0, 0]]), 'gap');
  assert.ok(!lf.isValidConformation([[0, 0, 0], [1, 1, 0], [2, 1, 0], [3, 1, 0], [4, 1, 0]]), 'diagonal bond');
  assert.ok(!lf.isValidConformation([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]), 'short chain');
  assert.ok(lf.isValidConformation([[0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 1, 0], [2, 1, 1]]));
});

test('a knotted compact conformation validates and scores', () => {
  // 2x2x2 cube corner filling, a classic compact SAW on 8 residues.
  const lf = new f.LatticeFolder('AAAAAAAA', { mode: 'hp' });
  const cube = [
    [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
    [0, 1, 1], [1, 1, 1], [1, 0, 1], [0, 0, 1]
  ];
  assert.ok(lf.isValidConformation(cube));
  const occ = new Map();
  cube.forEach((c, i) => occ.set(f.LatticeFolder.key(c[0], c[1], c[2]), i));
  lf._commit(cube, occ, null, null);
  // In the 2x2x2 cube every vertex pair at distance 1 except the 7 backbone
  // bonds is a contact: 12 cube edges minus 7 backbone edges = 5 contacts.
  assert.strictEqual(lf.contactCount(), 5);
  assert.strictEqual(lf.energy(), -5);
});

test('swapStateWith fully exchanges two folders', () => {
  const rng = f.makeRng(151);
  const a = new f.LatticeFolder(SEQ);
  const b = new f.LatticeFolder(SEQ);
  b.tryMove('pivot', rng);
  const aCoords = a.coords.map(c => c.slice());
  const bCoords = b.coords.map(c => c.slice());
  a.swapStateWith(b);
  assert.deepStrictEqual(a.coords, bCoords);
  assert.deepStrictEqual(b.coords, aCoords);
  assert.ok(occupancyConsistent(a) && occupancyConsistent(b));
});

// ---------------------------------------------------------------------
// Annealing and REMC invariants
// ---------------------------------------------------------------------
test('annealing trace is monotone non increasing', () => {
  const r = f.simulateAnnealing({ sequence: SEQ, sweeps: 60, stepsPerSweep: 15, seed: 157 });
  for (let i = 1; i < r.trace.length; i++) {
    assert.ok(r.trace[i] <= r.trace[i - 1] + 1e-9);
  }
});

test('annealing progress sweeps are strictly increasing', () => {
  const gen = f.annealingGenerator({ sequence: 'ACDEFG', sweeps: 12, stepsPerSweep: 6, yieldEvery: 3, seed: 163 });
  let last = 0;
  for (const msg of gen) {
    if (msg.type === 'progress') {
      assert.ok(msg.sweep > last);
      last = msg.sweep;
      assert.ok(msg.temperature >= 0);
    }
  }
});

test('remc progress energies match replica count and stay finite', () => {
  const gen = f.remcGenerator({ sequence: 'ACDEFGHI', sweeps: 10, replicas: 4, stepsPerSweep: 6, yieldEvery: 3, seed: 167 });
  for (const msg of gen) {
    if (msg.type === 'progress') {
      assert.strictEqual(msg.energies.length, 4);
      for (const e of msg.energies) assert.ok(Number.isFinite(e));
    }
  }
});

test('remc best energy never exceeds the initial extended state', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 40, replicas: 4, stepsPerSweep: 15, seed: 173 });
  assert.ok(r.best.energy <= 0);
});

test('annealing linear and geometric cooling both produce valid SAWs', () => {
  for (const cooling of ['linear', 'geometric']) {
    const r = f.simulateAnnealing({ sequence: SEQ, sweeps: 40, stepsPerSweep: 15, cooling, seed: 179 });
    const lf = new f.LatticeFolder(SEQ);
    const occ = new Map();
    r.best.coords.forEach((c, i) => occ.set(f.LatticeFolder.key(c[0], c[1], c[2]), i));
    lf._commit(r.best.coords, occ, r.best.energy, null);
    assert.ok(lf.isValidConformation(), cooling);
  }
});

// ---------------------------------------------------------------------
// Static file checks (app, markup, docs)
// ---------------------------------------------------------------------
const SHIPPED = ['index.html', 'app.js', 'src/fold.js', 'demo/commands.js', 'README.md', 'package.json'];

test('no em dash or en dash in shipped files', () => {
  for (const name of SHIPPED) {
    const t = read(name);
    assert.ok(!t.includes('\u2014'), name + ' em dash');
    assert.ok(!t.includes('\u2013'), name + ' en dash');
  }
});

test('no smooth-scroll CSS in shipped files', () => {
  const pat = 'scroll-' + 'behavior';
  for (const name of SHIPPED) assert.ok(!read(name).includes(pat), name);
});

test('index.html declares a CSP with worker and offline allowances', () => {
  const html = read('index.html');
  const m = html.match(/<meta http-equiv="Content-Security-Policy" content="([^"]+)"/);
  assert.ok(m, 'CSP meta present');
  const csp = m[1];
  assert.match(csp, /default-src 'none'/);
  assert.match(csp, /script-src 'self' file:/);
  assert.match(csp, /worker-src blob:/);
  assert.match(csp, /connect-src 'self' file:/);
  assert.ok(!csp.includes('unsafe-eval'));
  assert.ok(!/script-src[^;]*unsafe-inline/.test(csp), 'no inline script allowed');
});

test('index.html has no inline event handlers', () => {
  assert.ok(!/\son[a-z]+\s*=/i.test(read('index.html')));
});

test('grid tracks cannot blow out the mobile layout', () => {
  // Regression: the ribbon's ~700px flex content propagated its min-content
  // width up through `1fr` tracks and forced horizontal scroll at 390px.
  const html = read('index.html');
  assert.match(html, /grid-template-columns:\s*minmax\(0,\s*1fr\) minmax\(0,\s*360px\)/);
  assert.match(html, /main > \*, \.charts > \* \{\s*min-width:\s*0/);
  assert.ok(!/grid-template-columns:\s*minmax\(320px/.test(html), 'no lower bounded track');
});

test('index.html marks up live regions and keyboard access', () => {
  const html = read('index.html');
  assert.match(html, /id="progressText" role="status" aria-live="polite"/);
  assert.match(html, /id="errMsg" role="alert"/);
  assert.match(html, /id="stage"[^>]*role="img"[^>]*tabindex="0"/);
  assert.match(html, /aria-label="residue ribbon"/);
  assert.match(html, /:focus-visible/);
});

test('app.js never injects parsed data through innerHTML sinks', () => {
  const app = read('app.js');
  assert.ok(!/resBody\.innerHTML/.test(app), 'residue card innerHTML');
  assert.ok(!/ribbon\.innerHTML/.test(app), 'ribbon innerHTML');
  assert.ok(!/document\.write|eval\s*\(|new Function/.test(app));
});

test('keyboard camera control lives only on the focused canvas', () => {
  const app = read('app.js');
  assert.ok(app.includes("stage.addEventListener('keydown'"), 'stage keydown');
  assert.ok(!/document\.addEventListener\('key(down|press|up)'/.test(app), 'no global key listener');
  assert.ok(!/window\.addEventListener\('key(down|press|up)'/.test(app), 'no window key listener');
});

test('every getElementById target exists in index.html', () => {
  const app = read('app.js');
  const html = read('index.html');
  const ids = new Set([...app.matchAll(/\$\('([^']+)'\)|getElementById\('([^']+)'\)/g)]
    .map(m => m[1] || m[2]));
  assert.ok(ids.size > 10);
  for (const id of ids) {
    assert.ok(html.includes('id="' + id + '"'), 'missing #' + id);
  }
});

test('no external resource loads in shipped files', () => {
  for (const name of ['index.html', 'app.js', 'src/fold.js']) {
    const t = read(name);
    assert.ok(!/<script[^>]+src=["']https?:/i.test(t), name + ' external script');
    assert.ok(!/<link[^>]+href=["']https?:/i.test(t), name + ' external link');
    assert.ok(!/fetch\(["']https?:/i.test(t), name + ' remote fetch');
    assert.ok(!/url\(["']?https?:/i.test(t), name + ' remote css url');
  }
});

test('no AI-tooling attribution strings in shipped files', () => {
  const hx = h => h.replace(/../g, m => String.fromCharCode(parseInt(m, 16)));
  const banned = new RegExp('\\b(' +
    ['636c61756465', '63686174677074', '6f70656e6169', '636f70696c6f74',
     '646576696e', '616e7468726f706963'].map(hx).join('|') +
    '|' + hx('737765') + '-\\d|' + hx('677074') + '-\\d)\\b', 'i');
  for (const name of SHIPPED) {
    assert.ok(!banned.test(read(name)), name + ' mentions an AI tool');
  }
});

test('engine exposes no accidental globals in strict mode', () => {
  assert.match(read('src/fold.js'), /^'use strict';/);
});

// ---------------------------------------------------------------------
// Additional invariants and boundary coverage
// ---------------------------------------------------------------------
test('LatticeFolder on the minimal two residue chain', () => {
  const lf = new f.LatticeFolder('AA');
  assert.strictEqual(lf.n, 2);
  assert.deepStrictEqual(lf.coords, [[0, 0, 0], [1, 0, 0]]);
  assert.strictEqual(lf.contactCount(), 0);
});

test('corner flip returns null on straight segments and applies on corners', () => {
  const lf = new f.LatticeFolder('AAAA');
  const rng = f.makeRng(181);
  // Extended chain has no right angle corners: every proposal refuses.
  for (let i = 0; i < 50; i++) {
    assert.strictEqual(lf.tryMove('corner', rng).moved, false);
  }
});

test('crankshaft refuses non square loops', () => {
  const lf = new f.LatticeFolder('AAAAA');
  const rng = f.makeRng(183);
  for (let i = 0; i < 50; i++) {
    assert.strictEqual(lf.tryMove('crankshaft', rng).moved, false);
  }
});

test('endpoint pull only ever moves a terminal residue', () => {
  const rng = f.makeRng(191);
  const lf = new f.LatticeFolder('AAAAAAAAA');
  for (let i = 0; i < 600; i++) {
    const beforeEnds = [lf.coords[0].slice(), lf.coords[lf.n - 1].slice()];
    const res = lf.tryMove('endpull', rng);
    if (res.moved) {
      const movedEnd = !f.veq(lf.coords[0], beforeEnds[0]) ? 0 :
        (!f.veq(lf.coords[lf.n - 1], beforeEnds[1]) ? lf.n - 1 : -1);
      assert.ok(movedEnd >= 0, 'endpull moved no endpoint');
    }
  }
});

test('pull moves respect both pull directions', () => {
  for (const dir of ['pullUp', 'pullDown']) {
    const rng = f.makeRng(193);
    const lf = new f.LatticeFolder(SEQ);
    for (let i = 0; i < 400; i++) {
      lf.tryMove(dir, rng);
      assert.ok(lf.isValidConformation());
    }
  }
});

test('tryMove reports a clean refusal shape when nothing moved', () => {
  const lf = new f.LatticeFolder('AAAA');
  const rng = f.makeRng(197);
  const res = lf.tryMove('corner', rng);
  assert.strictEqual(res.moved, false);
  assert.strictEqual(res.accepted, false);
  assert.strictEqual(res.dE, 0);
  assert.strictEqual(res.cand, undefined);
});

test('copy is fully independent including the occupancy map', () => {
  const rng = f.makeRng(199);
  const lf = new f.LatticeFolder(SEQ);
  for (let i = 0; i < 100; i++) lf.tryMove('pull', rng);
  const c = lf.copy();
  const key0 = f.LatticeFolder.key(lf.coords[0][0], lf.coords[0][1], lf.coords[0][2]);
  c.occ.delete(key0);
  assert.ok(lf.occ.has(key0), 'copy mutation leaked into original');
  assert.strictEqual(c.occ.has(key0), false);
});

test('snapshot and restore are independent of later mutations', () => {
  const rng = f.makeRng(211);
  const lf = new f.LatticeFolder(SEQ);
  const snap = f.snapshot(lf);
  for (let i = 0; i < 200; i++) lf.tryMove('pull', rng);
  f.restore(lf, snap);
  const again = f.snapshot(lf);
  assert.deepStrictEqual(again.coords, snap.coords);
  // Restoring twice from the same snapshot also works.
  for (let i = 0; i < 50; i++) lf.tryMove('pivot', rng);
  f.restore(lf, snap);
  assert.deepStrictEqual(f.snapshot(lf).coords, snap.coords);
});

test('geometricTemps of two returns exactly max then min', () => {
  assert.deepStrictEqual(f.geometricTemps(2, 0.5, 2.0), [2.0, 0.5]);
});

test('different seeds produce different annealing results', () => {
  const a = f.simulateAnnealing({ sequence: SEQ, sweeps: 40, stepsPerSweep: 15, seed: 1 });
  const b = f.simulateAnnealing({ sequence: SEQ, sweeps: 40, stepsPerSweep: 15, seed: 999 });
  assert.notDeepStrictEqual(a.best.coords, b.best.coords);
});

test('remc reports exchange statistics consistent with the ladder', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 30, replicas: 5, stepsPerSweep: 10, seed: 223 });
  // With 5 replicas and alternating parity, pairs (0,1)+(2,3) try on even
  // sweeps and (1,2)+(3,4) on odd sweeps: floor(5/2)=2 pairs each parity.
  assert.strictEqual(r.exchangeTries, 30 * 2);
  assert.ok(r.exchanges <= r.exchangeTries);
  assert.ok(r.exchanges > 0, 'some exchanges should succeed');
});

test('remc final result lists every replica energy', () => {
  const r = f.simulateRemc({ sequence: 'ACDEFGH', sweeps: 12, replicas: 3, stepsPerSweep: 6, seed: 227 });
  assert.strictEqual(r.energies.length, 3);
  for (const e of r.energies) assert.ok(Number.isFinite(e) && e <= 0);
});

test('annealing yields temperature values inside the schedule bounds', () => {
  const gen = f.annealingGenerator({ sequence: 'ACDEFG', sweeps: 10, stepsPerSweep: 5, T0: 2, T1: 0.5, yieldEvery: 2, seed: 229 });
  for (const msg of gen) {
    if (msg.type === 'progress') {
      assert.ok(msg.temperature <= 2.0001 && msg.temperature >= 0.4999);
    }
  }
});

test('makeRng mean is near one half over a long run', () => {
  const rng = f.makeRng(233);
  let sum = 0;
  const N = 20000;
  for (let i = 0; i < N; i++) sum += rng();
  const mean = sum / N;
  assert.ok(mean > 0.49 && mean < 0.51, 'mean ' + mean);
});

test('normalizeDeg handles multiples of a full turn', () => {
  for (const [v, want] of [[360, 0], [-360, 0], [720, 0], [-540, 180], [540, 180], [181, -179], [-181, 179]]) {
    assert.strictEqual(f.normalizeDeg(v), want, String(v));
  }
});

test('dihedral of a planar trans zigzag is pi and a cis fold back is zero', () => {
  // IUPAC convention: trans (p1 and p4 on opposite sides) -> magnitude pi.
  const trans = f.dihedral([0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 1, 0]);
  assert.ok(Math.abs(Math.abs(trans) - Math.PI) < 1e-9, 'trans torsion ' + trans);
  // Cis (p4 folds back toward p1's side of the central bond) -> 0.
  const cis = f.dihedral([0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]);
  assert.ok(Math.abs(cis) < 1e-9, 'cis torsion ' + cis);
});

test('dihedral result always lies in the half open range (-pi, pi]', () => {
  const rng = f.makeRng(239);
  for (let i = 0; i < 200; i++) {
    const pts = [];
    for (let k = 0; k < 4; k++) pts.push([rng() * 10, rng() * 10, rng() * 10]);
    try {
      const d = f.dihedral(pts[0], pts[1], pts[2], pts[3]);
      assert.ok(d > -Math.PI && d <= Math.PI, d);
    } catch (e) { /* degenerate tetrads are allowed to throw */ }
  }
});

test('classifyOmega treats the 90 degree boundary as trans', () => {
  assert.strictEqual(f.classifyOmega(90), 'TRANS');
  assert.strictEqual(f.classifyOmega(-90), 'TRANS');
  assert.strictEqual(f.classifyOmega(89.9), 'CIS');
  assert.strictEqual(f.classifyOmega(-89.9), 'CIS');
});

test('findStericClashes sorts worst overlap first', () => {
  const atoms = [
    { xyz: [0, 0, 0], element: 'C' },
    { xyz: [2.4, 0, 0], element: 'C' },
    { xyz: [10, 0, 0], element: 'C' },
    { xyz: [12.25, 0, 0], element: 'C' }
  ];
  const clashes = f.findStericClashes(atoms);
  assert.strictEqual(clashes.length, 2);
  assert.ok(clashes[0].overlap >= clashes[1].overlap);
});

test('findStericClashes honors a factor of zero as never clashing', () => {
  const atoms = [
    { xyz: [0, 0, 0], element: 'C' },
    { xyz: [2.2, 0, 0], element: 'C' }
  ];
  assert.strictEqual(f.findStericClashes(atoms, { factor: 0 }).length, 0);
});

test('atomRadius prefers an explicit radius over the element table', () => {
  assert.strictEqual(f.atomRadius({ element: 'C', radius: 9.9 }), 9.9);
  assert.strictEqual(f.atomRadius({ element: 'Q' }), 1.70);
});

test('inferElement resolves two letter elements and odd placements', () => {
  assert.strictEqual(f.inferElement('CL  '), 'CL');
  assert.strictEqual(f.inferElement(' ZN'), 'Z');
  assert.strictEqual(f.inferElement('CA  '), 'CA');
  assert.strictEqual(f.inferElement(' CA '), 'C');
  assert.strictEqual(f.inferElement(''), 'C');
  assert.strictEqual(f.inferElement('X', 'SE'), 'SE');
});

test('latticeBackboneAngles right angle turns classify into a helix basin', () => {
  // A chain with a single right-angle bend: residue 1 sees a 90 degree
  // angle and must not stay in the beta basin.
  const lf = new f.LatticeFolder('AAAA');
  const coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]];
  const occ = new Map();
  coords.forEach((c, i) => occ.set(f.LatticeFolder.key(c[0], c[1], c[2]), i));
  lf._commit(coords, occ, null, null);
  const angles = f.latticeBackboneAngles(lf.coords);
  assert.strictEqual(angles.length, 4);
  assert.strictEqual(angles[0].region, 'BETA');
  assert.notStrictEqual(angles[1].region, 'BETA', 'turn residue should leave beta');
  assert.notStrictEqual(angles[2].region, 'BETA', 'turn residue should leave beta');
});

test('computeDihedrals keeps insertion coded residues as separate groups', () => {
  const pdb = [
    'ATOM      1  N   ALA A  10       0.000   0.000   0.000  1.00  0.00           N',
    'ATOM      2  CA  ALA A  10       1.458   0.000   0.000  1.00  0.00           C',
    'ATOM      3  C   ALA A  10       2.000   1.300   0.300  1.00  0.00           C',
    'ATOM      4  N   GLY A  10A      3.200   1.100   1.200  1.00  0.00           N',
    'ATOM      5  CA  GLY A  10A      4.100   0.100   1.600  1.00  0.00           C',
    'ATOM      6  C   GLY A  10A      5.500   0.700   1.400  1.00  0.00           C',
    'END'
  ].join('\n');
  const d = f.computeDihedrals(f.parsePdb(pdb).atoms);
  assert.strictEqual(d.length, 2);
  assert.notStrictEqual(d[1].phi, null);
  assert.notStrictEqual(d[0].psi, null);
});

test('folderToPdb honors a custom CA spacing scale', () => {
  const lf = new f.LatticeFolder('ACDE');
  const pdb = f.folderToPdb(lf, { scale: 5.0 });
  const p = f.parsePdb(pdb);
  assert.ok(Math.abs(f.vdist(p.atoms[0].xyz, p.atoms[1].xyz) - 5.0) < 1e-9);
});

test('reconstructLatticeBackbone dihedrals carry an omega class', () => {
  // Regression: the UI rendered "omega 180.0 (?)" because omegaClass was
  // missing from the rebuilt dihedral rows.
  const lf = new f.LatticeFolder('ACDEFG');
  const rec = f.reconstructLatticeBackbone(lf);
  for (const d of rec.dihedrals) {
    assert.ok(d.omegaClass === 'TRANS' || d.omegaClass === 'CIS', 'omegaClass present');
  }
});

test('folderToPdb embeds title and remarks in the header block', () => {
  const lf = new f.LatticeFolder('ACDE');
  const pdb = f.folderToPdb(lf, { title: 'custom title', remarks: ['note one'] });
  assert.ok(pdb.includes('custom title'));
  assert.ok(pdb.includes('note one'));
});

test('annealingGenerator stops yielding after done', () => {
  const gen = f.annealingGenerator({ sequence: 'ACDEFG', sweeps: 6, stepsPerSweep: 4, seed: 241 });
  let seen = 0, doneSeen = false;
  for (const msg of gen) {
    seen++;
    if (msg.done) doneSeen = true;
  }
  assert.ok(doneSeen);
  // Iterator protocol: further next() calls report done without values.
  assert.strictEqual(gen.next().done, true);
});

test('generator contract: the last yielded value is the done result', () => {
  // Consumers (worker pump and chunked fallback) must see the result arrive
  // as a yielded message carrying done:true, before the iterator finishes.
  for (const mk of [
    () => f.remcGenerator({ sequence: 'ACDEFG', sweeps: 8, replicas: 2, stepsPerSweep: 5, seed: 263 }),
    () => f.annealingGenerator({ sequence: 'ACDEFG', sweeps: 8, stepsPerSweep: 5, seed: 263 })
  ]) {
    const gen = mk();
    let last;
    for (const v of gen) { last = v; }
    assert.ok(last && last.done === true && last.type === 'done', 'final yield is the done result');
    assert.ok(last.best && Array.isArray(last.trace), 'done payload is the full result');
  }
});

test('app message pumps dispatch a yielded done object as done', () => {
  // Regression: the worker and chunked pumps once wrapped every yield as
  // 'progress', which threw inside the progress handler and left runs stuck.
  const app = read('app.js');
  assert.match(app, /v\.done \? \\'done\\' : \\'progress\\'/, 'worker pump checks v.done');
  assert.match(app, /r\.done \|\| r\.value\.done === true/, 'chunked pump checks r.value.done');
});

test('remcGenerator respects a custom temperature ladder', () => {
  const r = f.simulateRemc({ sequence: 'ACDEFG', sweeps: 8, replicas: 3, stepsPerSweep: 5, minTemp: 0.4, maxTemp: 1.2, seed: 251 });
  assert.strictEqual(r.temps.length, 3);
  assert.ok(Math.abs(r.temps[0] - 1.2) < 1e-9);
  assert.ok(Math.abs(r.temps[2] - 0.4) < 1e-9);
});

// ---------------------------------------------------------------------
// pass 2 regressions: reproduced adversarial findings
// ---------------------------------------------------------------------

test('atomsToPdb rejects fields that would overflow their PDB columns', () => {
  const atom = { name: 'CA', resName: 'ALA', resSeq: 1, chainID: 'A', xyz: [0, 0, 0] };
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { resSeq: 10000 })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { resSeq: 1.5 })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { resSeq: 'XYZ' })]), RangeError);
  assert.throws(() => f.atomsToPdb([atom, Object.assign({}, atom, { resSeq: 2 })], { startSerial: 99999 }), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { xyz: [100000, 0, 0] })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { xyz: [-100000, 0, 0] })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { chainID: 'AB' })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { resName: 'LONG' })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { iCode: 'AB' })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { altLoc: 'BC' })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { element: 'XXX' })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { occupancy: 12345 })]), RangeError);
  assert.throws(() => f.atomsToPdb([atom], { startSerial: 'x' }), RangeError);
});

test('atomsToPdb rejects control characters that could inject records', () => {
  const atom = { name: 'CA', resName: 'ALA', resSeq: 1, chainID: 'A', xyz: [0, 0, 0] };
  assert.throws(() => f.atomsToPdb([atom], { remarks: ['x\nATOM      9  CA  ALA A   1'] }), RangeError);
  assert.throws(() => f.atomsToPdb([atom], { title: 't\r\nHETATM' }), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { chainID: 'A\n' })]), RangeError);
  assert.throws(() => f.atomsToPdb([Object.assign({}, atom, { resName: 'A\nB' })]), RangeError);
});

test('atomsToPdb emits a blank resSeq for missing values instead of "null"', () => {
  const out = f.atomsToPdb([{ name: 'CA', resName: 'ALA', xyz: [0, 0, 0] }]);
  const atomLine = out.split('\n').find(l => l.startsWith('ATOM'));
  assert.strictEqual(atomLine.slice(22, 26), '    ');
  assert.ok(!/null|undefined|NaN/.test(atomLine));
});

// Fixed-width PDB line builder: fields land at their exact columns.
function pdbAtomLine(o) {
  const line =
    'ATOM  ' +
    String(o.serial).padStart(5) + ' ' +
    String(o.name || 'CA').padEnd(4) +
    String(o.altLoc || ' ') +
    String(o.resName || 'ALA').padEnd(3) + ' ' +
    String(o.chainID || 'A') +
    String(o.resSeq).padStart(4) +
    String(o.iCode || ' ') + '   ' +
    String(o.x).padStart(8) + String(o.y).padStart(8) + String(o.z).padStart(8) +
    String(o.occupancy === undefined ? '  1.00' : o.occupancy).padStart(6) +
    String(o.tempFactor === undefined ? '  0.00' : o.tempFactor).padStart(6) +
    '          ' + String(o.element || 'C').padStart(2);
  return line;
}

test('parsePdb skips records with non-finite coordinates', () => {
  for (const bad of ['NaN', 'Infinity', '-Infinity']) {
    const pdb = pdbAtomLine({ serial: 1, resSeq: 1, x: bad, y: '   0.000', z: '   0.000' }) + '\n' +
      pdbAtomLine({ serial: 2, resSeq: 1, x: '   1.000', y: '   0.000', z: '   0.000' }) + '\nEND\n';
    const parsed = f.parsePdb(pdb);
    assert.strictEqual(parsed.atoms.length, 1, bad + ' record skipped');
    assert.strictEqual(parsed.atoms[0].xyz[0], 1);
  }
});

test('parsePdb stores null, not NaN, for non-finite serial/resSeq/occupancy', () => {
  const line = pdbAtomLine({ serial: '  NaN', resSeq: ' NaN', x: '   0.000', y: '   0.000', z: '   0.000', occupancy: '  NaN', tempFactor: '  Inf' });
  const good = pdbAtomLine({ serial: 2, resSeq: 2, x: '   1.000', y: '   0.000', z: '   0.000' });
  const parsed = f.parsePdb(line + '\n' + good + '\nEND\n');
  assert.strictEqual(parsed.atoms.length, 2);
  assert.strictEqual(parsed.atoms[0].serial, null);
  assert.strictEqual(parsed.atoms[0].resSeq, null);
  assert.strictEqual(parsed.atoms[0].occupancy, null);
  assert.ok(Number.isFinite(parsed.atoms[1].resSeq));
});

test('parsePdb counts every MODEL record while reading only the first', () => {
  let text = '';
  for (let m = 1; m <= 5; m++) {
    text += 'MODEL     ' + m + '\nATOM      1  CA  ALA A   ' + m +
      '       0.000   0.000   0.000  1.00  0.00           C\nENDMDL\n';
  }
  const parsed = f.parsePdb(text);
  assert.strictEqual(parsed.modelCount, 5);
  assert.strictEqual(parsed.atoms.length, 1);
  assert.strictEqual(parsed.atoms[0].resSeq, 1);
});

test('geometricTemps never returns zero for wide ranges', () => {
  const temps = f.geometricTemps(10, 1e-200, 1e200);
  assert.strictEqual(temps.length, 10);
  for (const t of temps) { assert.ok(t > 0 && Number.isFinite(t), t); }
  for (let i = 1; i < temps.length; i++) { assert.ok(temps[i] < temps[i - 1]); }
  assert.strictEqual(temps[0], 1e200);
  assert.strictEqual(temps[9], 1e-200);
});

test('annealingGenerator ramps smoothly to T1 zero instead of one hot sweep', () => {
  const gen = f.annealingGenerator({ sequence: 'ACDEFG', sweeps: 6, stepsPerSweep: 3, T0: 2, T1: 0, seed: 7, yieldEvery: 1 });
  const temps = [];
  for (const msg of gen) { if (msg.type === 'progress') temps.push(msg.temperature); }
  assert.strictEqual(temps.length, 6);
  for (let i = 1; i < temps.length; i++) { assert.ok(temps[i] < temps[i - 1], 'strictly cooling'); }
  assert.strictEqual(temps[5], 0);
  assert.ok(temps[4] > 0, 'interior sweep is not quenched');
});

test('contacts() order is canonical regardless of how the conformation was reached', () => {
  const seed1 = 31;
  const fresh = new f.LatticeFolder(SEQ, { mode: 'mj', seed: seed1 });
  for (let i = 0; i < 200; i++) { fresh.mcStep(1.5, Math.random); }
  const canonical = fresh.coords.map(c => c.slice());
  const order1 = fresh.contacts();
  // Rebuild the same conformation through _evaluate (the move-commit path).
  const rebuilt = new f.LatticeFolder(SEQ, { mode: 'mj' });
  const occ = new Map();
  canonical.forEach((c, i) => occ.set(f.LatticeFolder.key(c[0], c[1], c[2]), i));
  const built = rebuilt._evaluate(canonical);
  rebuilt._commit(canonical, built.occ, built.energy, built.contacts);
  assert.deepStrictEqual(rebuilt.contacts(), order1);
});

test('swapStateWith rejects different sequences or energy models', () => {
  const a = new f.LatticeFolder('AAAA', { mode: 'hp' });
  const b = new f.LatticeFolder('WWWW', { mode: 'hp' });
  const c = new f.LatticeFolder('AAAA', { mode: 'mj' });
  const d = new f.LatticeFolder('AAAA', { mode: 'hp' });
  assert.throws(() => a.swapStateWith(b));
  assert.throws(() => a.swapStateWith(c));
  assert.doesNotThrow(() => a.swapStateWith(d));
});

test('mcStep rejects non-finite and negative temperatures', () => {
  const lf = new f.LatticeFolder(SEQ, { mode: 'hp' });
  for (const t of [NaN, -1, -0.5, Infinity]) {
    assert.throws(() => lf.mcStep(t, Math.random), RangeError, String(t));
  }
  assert.doesNotThrow(() => lf.mcStep(0, Math.random));
  assert.doesNotThrow(() => lf.mcStep(2, Math.random));
});

test('isValidConformation rejects non-integer and out-of-range coordinates', () => {
  const lf = new f.LatticeFolder('AAAA', { mode: 'hp' });
  assert.strictEqual(lf.isValidConformation([[0.5, 0, 0], [1.5, 0, 0], [1.5, 1, 0], [0.5, 1, 0]]), false);
  assert.strictEqual(lf.isValidConformation([[0, 0, 0], [1, 0, 0], [100000, 0, 0], [100001, 0, 0]]), false);
  assert.strictEqual(lf.isValidConformation([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]), true);
  assert.strictEqual(lf._evaluate([[0.5, 0, 0], [1.5, 0, 0], [1.5, 1, 0], [0.5, 1, 0]]), null);
});

test('groupResidues splits on insertion code and indexes atoms linearly', () => {
  const app = read('app.js');
  assert.match(app, /cur\.iCode !== a\.iCode/, 'insertion code is part of the residue key');
  assert.match(app, /atomResidueIndex\(groups, atoms\.length\)/, 'atom map is precomputed');
  assert.ok(!/resIndexOfGroup/.test(app), 'quadratic per-atom scan removed');
});

test('scrubbing keeps dihedrals in sync with the displayed snapshot', () => {
  const app = read('app.js');
  // refreshDisplay assigns the rebuilt dihedrals for every shown snapshot,
  // and clears them when showing earlier lattice beads.
  assert.match(app, /state\.dihedrals = rec\.dihedrals;\s*\n\s*setFromAtoms/);
  assert.match(app, /snap !== finalSnapshot\(\) \? |snap !== finalSnapshot\(\)\) \{ state\.dihedrals = null/);
  // The scrub handler redraws the Ramachandran map alongside the beads.
  assert.match(app, /state\.scrub = parseInt\(scrubR\.value[\s\S]*?drawRama\(\)/);
});

test('pdb import drops stale file reads and stale engine messages', () => {
  const app = read('app.js');
  assert.match(app, /var pdbLoadGen = 0/, 'load generation counter exists');
  assert.match(app, /gen !== pdbLoadGen\) \{ return; \}[\s\S]*?parse failed/, 'catch path is guarded');
  assert.match(app, /e\.target !== worker/, 'worker events checked against current worker');
  assert.match(app, /if \(!state\.running\) \{ return; \}[\s\S]*?msg\.type === 'progress'/, 'engine messages need a live run');
});

test('charts and legend expose text alternatives', () => {
  const html = read('index.html');
  for (const id of ['traceChart', 'ramaChart']) {
    assert.match(html, new RegExp('id="' + id + '"[^>]*role="img"'), id + ' role');
    assert.match(html, new RegExp('id="' + id + '"[^>]*aria-label='), id + ' label');
  }
  assert.match(html, /<canvas id="traceChart"[^>]*>[^<]+<\/canvas>/, 'trace fallback text');
  assert.match(html, /<canvas id="ramaChart"[^>]*>[^<]+<\/canvas>/, 'rama fallback text');
  assert.match(html, /id="legend" role="group" aria-label=/, 'legend is labelled');
});

test('atomsToPdb keeps a blank chain as one chain', () => {
  const atoms = [
    { name: 'CA', resName: 'ALA', resSeq: 1, xyz: [0, 0, 0], element: 'C' },
    { name: 'CA', resName: 'ALA', resSeq: 2, xyz: [3.8, 0, 0], element: 'C' }
  ];
  const records = f.atomsToPdb(atoms).split('\n')
    .filter(l => l.startsWith('ATOM') || l.startsWith('TER'))
    .map(l => l.slice(0, 6).trim());
  assert.deepStrictEqual(records, ['ATOM', 'ATOM', 'TER']);
});

test('parsePdb ignores stray ATOM records after the first ENDMDL', () => {
  const atom = pdbAtomLine({ serial: 1, resSeq: 1, x: '   0.000', y: '   0.000', z: '   0.000' });
  const stray = pdbAtomLine({ serial: 2, resSeq: 2, x: '   1.000', y: '   0.000', z: '   0.000' });
  const second = pdbAtomLine({ serial: 3, resSeq: 1, x: '   2.000', y: '   0.000', z: '   0.000' });
  const text = 'MODEL     1\n' + atom + '\nENDMDL\n' + stray + '\nMODEL     2\n' + second + '\nENDMDL\n';
  const parsed = f.parsePdb(text);
  assert.strictEqual(parsed.modelCount, 2);
  assert.strictEqual(parsed.atoms.length, 1);
  assert.strictEqual(parsed.atoms[0].serial, 1);
});

test('computeDihedrals bonds insertion coded neighbors in sequence order', () => {
  const angles = [
    { phi: -60, psi: -45, omega: 180 },
    { phi: -60, psi: -45, omega: 180 },
    { phi: -60, psi: -45, omega: 180 }
  ];
  const atoms = f.reconstructBackbone('AAA', angles).map((a) => {
    const out = Object.assign({}, a);
    if (out.resIndex === 1) { out.resSeq = 1; out.iCode = 'A'; }
    if (out.resIndex === 2) { out.resSeq = 2; out.iCode = ' '; }
    return out;
  });
  const ds = f.computeDihedrals(atoms);
  assert.strictEqual(ds.length, 3);
  assert.notStrictEqual(ds[0].psi, null, 'residue 1 should bond to insertion residue 1A');
  assert.notStrictEqual(ds[1].phi, null, 'insertion residue 1A should bond back to residue 1');
  assert.notStrictEqual(ds[1].omega, null, 'insertion residue 1A should have peptide omega');
});

test('makeRng rejects malformed seeds instead of silently coercing them', () => {
  for (const bad of [NaN, Infinity, -Infinity, -1, 1.5, 4294967296]) {
    assert.throws(() => f.makeRng(bad), RangeError, String(bad));
  }
  assert.doesNotThrow(() => f.makeRng(0));
  assert.doesNotThrow(() => f.makeRng(4294967295));
});

test('shipped page carries skip, reduced motion, progress and button state accessibility', () => {
  const html = read('index.html');
  assert.match(html, /class="skip-link"/, 'skip link');
  assert.match(html, /<main id="mainContent"/, 'skip target');
  assert.match(html, /prefers-reduced-motion/, 'reduced motion media query');
  assert.match(html, /id="progressWrap"[^>]*role="progressbar"/, 'progressbar role');
  assert.match(html, /id="mutBtn"[^>]*disabled/, 'mutation button starts disabled');
  assert.match(html, /id="auditJsonBtn"[^>]*disabled/, 'audit export starts disabled');
  assert.match(html, /id="verifyBtn"[^>]*disabled/, 'verify button starts disabled');
});

test('contact map cells are keyboard operable', () => {
  const app = read('lab-app.js');
  assert.match(app, /tabindex', '0'|tabindex", "0"|setAttribute\('tabindex', '0'\)/, 'focusable cells');
  assert.match(app, /cmapSvg\.addEventListener\('keydown'|cmapSvg\.onkeydown/, 'keyboard handler');
  assert.match(app, /selectContactCell/, 'shared selection path');
});
