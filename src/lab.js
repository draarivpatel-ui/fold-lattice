'use strict';

/**
 * fold-lattice lab: a 2D/3D Dill HP model workbench that runs beside the 3D
 * lattice engine in src/fold.js.
 *
 * Contents:
 *   - HP folding on square (2D) and cubic (3D) lattices as self avoiding
 *     walks of hydrophobic (H) and polar (P) beads, scored at -1 per
 *     non bonded H-H lattice contact.
 *   - Deterministic seeded algorithms side by side: random SAW sampling,
 *     greedy local pivot descent, fixed temperature Metropolis Monte Carlo,
 *     simulated annealing, and exact exhaustive enumeration of the ground
 *     state over symmetry reduced self avoiding walks.
 *   - Metrics: contact energy E, H-H contact count, radius of gyration R_g,
 *     end to end distance R_ee, optimality gap against the exact ground
 *     state, and a point mutation stress scan (H <-> P, delta E).
 *   - L x L contact matrix, folding funnel sampling, FIPS 180-4 SHA-256 in
 *     pure JavaScript, canonical JSON checksums, checksummed audit reports
 *     (JSON, CSV, Markdown) with CSV formula injection escaping, URL hash
 *     bookmarks, and an in page self test runner.
 *
 * The Dill HP lattice model is a simplified statistical physics toy model:
 * it demonstrates folding thermodynamics and search algorithms, it is not
 * an atomistic force field and not an AlphaFold style structure predictor.
 */

var F = (typeof module === 'object' && module.exports)
  ? require('./fold.js')
  : window.foldl;

var LAB_VERSION = '2.0.0';

var makeRng = F.makeRng;

// ---------------------------------------------------------------------------
// Lattice directions. Lab coordinates are always length 3 arrays; 2D
// conformations keep z = 0 so the vector helpers and radius of gyration from
// the main engine apply unchanged.
// ---------------------------------------------------------------------------
var DIRS = {
  2: [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
  3: [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
};

function validateDim(dim) {
  if (dim !== 2 && dim !== 3) {
    throw new Error('dimension must be 2 (square lattice) or 3 (cubic lattice)');
  }
  return dim;
}

function validateHpSequence(sequence) {
  if (typeof sequence !== 'string') {
    throw new Error('HP sequence must be a string');
  }
  var seq = sequence.toUpperCase().replace(/\s/g, '');
  if (seq.length < 2) {
    throw new Error('HP sequence needs at least 2 residues');
  }
  for (var i = 0; i < seq.length; i++) {
    var ch = seq[i];
    if (ch !== 'H' && ch !== 'P') {
      throw new Error('HP sequence holds only H and P letters, found "' + ch + '" at position ' + i);
    }
  }
  return seq;
}

function seedOption(value, fallback) {
  return F.validateSeed(value === undefined ? fallback : value);
}

// Occupancy keys; identical injective scheme to the main engine.
var LKEY = 200001;
var LOFF = 100000;
function labKey(x, y, z) {
  return (x + LOFF) + LKEY * (y + LOFF) + LKEY * LKEY * (z + LOFF);
}

function cloneCoords(coords) {
  return coords.map(function (c) { return c.slice(); });
}

function manhattan(a, b) {
  return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);
}

/**
 * Full conformation validation: integer cells inside the key range, self
 * avoidance, unit step adjacency, and z = 0 everywhere for the 2D lattice.
 */
function isValidSaw(coords, dim) {
  validateDim(dim);
  if (!Array.isArray(coords) || coords.length < 1) { return false; }
  var seen = {};
  for (var i = 0; i < coords.length; i++) {
    var c = coords[i];
    if (!c || c.length !== 3) { return false; }
    if (!Number.isInteger(c[0]) || !Number.isInteger(c[1]) || !Number.isInteger(c[2])) {
      return false;
    }
    if (dim === 2 && c[2] !== 0) { return false; }
    if (Math.abs(c[0]) >= LOFF || Math.abs(c[1]) >= LOFF || Math.abs(c[2]) >= LOFF) {
      return false;
    }
    var k = labKey(c[0], c[1], c[2]);
    if (seen[k]) { return false; }
    seen[k] = true;
  }
  for (var j = 1; j < coords.length; j++) {
    if (manhattan(coords[j], coords[j - 1]) !== 1) { return false; }
  }
  return true;
}

function occupancyOf(coords) {
  var occ = new Map();
  for (var i = 0; i < coords.length; i++) {
    occ.set(labKey(coords[i][0], coords[i][1], coords[i][2]), i);
  }
  return occ;
}

/**
 * All non bonded lattice adjacent residue pairs [i, j] with j > i + 1,
 * returned in canonical sorted order.
 */
function contactPairs(coords, dim) {
  var occ = occupancyOf(coords);
  var dirs = DIRS[dim];
  var pairs = [];
  for (var i = 0; i < coords.length; i++) {
    var p = coords[i];
    for (var d = 0; d < dirs.length; d++) {
      var ax = dirs[d];
      var j = occ.get(labKey(p[0] + ax[0], p[1] + ax[1], p[2] + ax[2]));
      if (j !== undefined && j > i + 1) { pairs.push([i, j]); }
    }
  }
  pairs.sort(function (a, b) { return a[0] - b[0] || a[1] - b[1]; });
  return pairs;
}

/** Non bonded adjacent pairs where both residues are hydrophobic. */
function hhPairs(seq, coords, dim) {
  return contactPairs(coords, dim).filter(function (pr) {
    return seq[pr[0]] === 'H' && seq[pr[1]] === 'H';
  });
}

/** Dill HP contact energy: -1 per non bonded H-H lattice contact. */
function energyOf(seq, coords, dim) {
  // The || 0 folds -0 into +0 so strict equality checks see a plain zero.
  return -hhPairs(seq, coords, dim).length || 0;
}

function endToEnd(coords) {
  var a = coords[0], b = coords[coords.length - 1];
  var dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Metric bundle for one conformation:
 *   energy       H-H non bonded contact energy (<= 0)
 *   hhContacts   H-H contacts that drive the energy
 *   contactCount every non bonded lattice adjacent pair
 *   rg           radius of gyration in lattice units
 *   ree          end to end Euclidean distance in lattice units
 */
function metricsOf(seq, coords, dim) {
  var pairs = contactPairs(coords, dim);
  var hh = 0;
  for (var i = 0; i < pairs.length; i++) {
    if (seq[pairs[i][0]] === 'H' && seq[pairs[i][1]] === 'H') { hh++; }
  }
  return {
    energy: -hh || 0,
    hhContacts: hh,
    contactCount: pairs.length,
    rg: F.radiusOfGyration(coords),
    ree: endToEnd(coords)
  };
}

/** Evaluate a candidate: null when invalid, else { energy, pairs }. */
function evaluateConformation(seq, cand, dim) {
  if (!isValidSaw(cand, dim)) { return null; }
  var pairs = contactPairs(cand, dim);
  var e = 0;
  for (var i = 0; i < pairs.length; i++) {
    if (seq[pairs[i][0]] === 'H' && seq[pairs[i][1]] === 'H') { e -= 1; }
  }
  return { energy: e, pairs: pairs };
}

// ---------------------------------------------------------------------------
// Move set: pivot, corner flip, crankshaft, end flip. Candidates are plain
// coordinate arrays; callers evaluate and commit explicitly.
// ---------------------------------------------------------------------------
function pivotAt(coords, i, turns, dim) {
  var n = coords.length;
  if (i < 0 || i >= n - 1) { return cloneCoords(coords); }
  var p = coords[i];
  var cand = new Array(n);
  for (var k = 0; k <= i; k++) { cand[k] = coords[k].slice(); }
  if (dim === 3) {
    var u = [coords[i + 1][0] - p[0], coords[i + 1][1] - p[1], coords[i + 1][2] - p[2]];
    var R = F.rotationAboutAxis(u, turns);
    for (var m = i + 1; m < n; m++) {
      var v = [coords[m][0] - p[0], coords[m][1] - p[1], coords[m][2] - p[2]];
      cand[m] = [
        p[0] + R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        p[1] + R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        p[2] + R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]
      ];
    }
  } else {
    for (var q = i + 1; q < n; q++) {
      var dx = coords[q][0] - p[0];
      var dy = coords[q][1] - p[1];
      var t = ((turns % 4) + 4) % 4;
      if (t === 1) { cand[q] = [p[0] - dy, p[1] + dx, p[2]]; }
      else if (t === 2) { cand[q] = [p[0] - dx, p[1] - dy, p[2]]; }
      else if (t === 3) { cand[q] = [p[0] + dy, p[1] - dx, p[2]]; }
      else { cand[q] = coords[q].slice(); }
    }
  }
  return cand;
}

function pivotCandidate(coords, dim, rng) {
  var n = coords.length;
  var i = Math.floor(rng() * (n - 1));
  var turns = 1 + Math.floor(rng() * 3);
  return pivotAt(coords, i, turns, dim);
}

function cornerCandidate(coords, dim, rng) {
  var n = coords.length;
  if (n < 3) { return null; }
  var i = 1 + Math.floor(rng() * (n - 2));
  var p = coords[i];
  var c = [
    coords[i - 1][0] + coords[i + 1][0] - p[0],
    coords[i - 1][1] + coords[i + 1][1] - p[1],
    coords[i - 1][2] + coords[i + 1][2] - p[2]
  ];
  if (manhattan(c, p) !== 2) { return null; }
  var cand = cloneCoords(coords);
  cand[i] = c;
  return cand;
}

function crankshaftCandidate(coords, dim, rng) {
  var n = coords.length;
  if (n < 4) { return null; }
  var j = Math.floor(rng() * (n - 3));
  var p0 = coords[j];
  var d1 = [coords[j + 1][0] - p0[0], coords[j + 1][1] - p0[1], coords[j + 1][2] - p0[2]];
  var d2 = [coords[j + 2][0] - coords[j + 1][0], coords[j + 2][1] - coords[j + 1][1], coords[j + 2][2] - coords[j + 1][2]];
  var d3 = [coords[j + 3][0] - coords[j + 2][0], coords[j + 3][1] - coords[j + 2][1], coords[j + 3][2] - coords[j + 2][2]];
  if (manhattan(d1, [0, 0, 0]) !== 1 || manhattan(d2, [0, 0, 0]) !== 1 || manhattan(d3, [0, 0, 0]) !== 1) {
    return null;
  }
  if (!(d3[0] === -d1[0] && d3[1] === -d1[1] && d3[2] === -d1[2])) { return null; }
  var a = [p0[0] - d1[0], p0[1] - d1[1], p0[2] - d1[2]];
  var b = [a[0] + d2[0], a[1] + d2[1], a[2] + d2[2]];
  var cand = cloneCoords(coords);
  cand[j + 1] = a;
  cand[j + 2] = b;
  return cand;
}

function endFlipCandidate(coords, dim, rng) {
  var n = coords.length;
  if (n < 3) { return null; }
  var atEnd = rng() < 0.5;
  var e = atEnd ? n - 1 : 0;
  var nb = atEnd ? n - 2 : 1;
  var occ = occupancyOf(coords);
  var free = [];
  var dirs = DIRS[dim];
  for (var d = 0; d < dirs.length; d++) {
    var ax = dirs[d];
    var cell = [coords[nb][0] + ax[0], coords[nb][1] + ax[1], coords[nb][2] + ax[2]];
    if (!occ.has(labKey(cell[0], cell[1], cell[2]))) { free.push(cell); }
  }
  if (!free.length) { return null; }
  var cand = cloneCoords(coords);
  cand[e] = free[Math.floor(rng() * free.length)];
  return cand;
}

var MOVE_NAMES = ['pivot', 'corner', 'crankshaft', 'endflip'];

function moveCandidate(name, coords, dim, rng) {
  switch (name) {
    case 'pivot': return pivotCandidate(coords, dim, rng);
    case 'corner': return cornerCandidate(coords, dim, rng);
    case 'crankshaft': return crankshaftCandidate(coords, dim, rng);
    case 'endflip': return endFlipCandidate(coords, dim, rng);
    default: throw new Error('unknown lab move: ' + name);
  }
}

/**
 * One Metropolis step against a { seq, dim, coords, energy } fold state.
 * Move mix: pivot 45%, corner 20%, crankshaft 15%, end flip 20%.
 */
function mcStepLab(fold, rng, temperature) {
  if (!Number.isFinite(temperature) || temperature < 0) {
    throw new RangeError('temperature must be a nonnegative finite number');
  }
  var r = rng();
  var type;
  if (r < 0.45) { type = 'pivot'; }
  else if (r < 0.65) { type = 'corner'; }
  else if (r < 0.80) { type = 'crankshaft'; }
  else { type = 'endflip'; }
  var cand = moveCandidate(type, fold.coords, fold.dim, rng);
  if (!cand) { return { moved: false, accepted: false, dE: 0, type: type }; }
  var built = evaluateConformation(fold.seq, cand, fold.dim);
  if (!built) { return { moved: false, accepted: false, dE: 0, type: type }; }
  var dE = built.energy - fold.energy;
  var accept = dE <= 0 || (temperature > 0 && rng() < Math.exp(-dE / temperature));
  if (accept) {
    fold.coords = cand;
    fold.energy = built.energy;
  }
  return { moved: true, accepted: accept, dE: dE, type: type, energy: built.energy };
}

// ---------------------------------------------------------------------------
// Exact enumeration: depth first generation of self avoiding walks with the
// first step fixed to +x and the second step restricted to {+x, +y}. This
// quotients rotations and translations but keeps mirror twins, so every
// symmetry class is visited at least once; extra copies do not change the
// minimum, so the reported ground state is still exact.
// ---------------------------------------------------------------------------
function canonicalEnumerate(n, dim, visit, maxStates) {
  validateDim(dim);
  var dirs = DIRS[dim];
  var coords = new Array(n);
  var occ = new Set();
  var state = { states: 0, truncated: false };
  coords[0] = [0, 0, 0];
  occ.add(labKey(0, 0, 0));
  if (n === 1) {
    state.states++;
    visit(coords.slice(0, 1));
    return state;
  }
  coords[1] = [1, 0, 0];
  occ.add(labKey(1, 0, 0));
  if (n === 2) {
    state.states++;
    visit(coords.slice(0, 2));
    return state;
  }
  // Second step: straight (+x) or the canonical perpendicular (+y).
  var SECOND = [[1, 0, 0], [0, 1, 0]];

  function dfs(k) {
    if (state.truncated) { return; }
    if (k === n) {
      state.states++;
      visit(coords.slice());
      if (state.states >= maxStates) { state.truncated = true; }
      return;
    }
    var choices = k === 2 ? SECOND : dirs;
    var base = coords[k - 1];
    for (var d = 0; d < choices.length; d++) {
      var ax = choices[d];
      var cell = [base[0] + ax[0], base[1] + ax[1], base[2] + ax[2]];
      var key = labKey(cell[0], cell[1], cell[2]);
      if (occ.has(key)) { continue; }
      coords[k] = cell;
      occ.add(key);
      dfs(k + 1);
      occ.delete(key);
      if (state.truncated) { return; }
    }
  }
  dfs(2);
  return state;
}

/** Number of canonical self avoiding walks of length n (for tests). */
function countSaws(n, dim, maxStates) {
  var count = { c: 0 };
  canonicalEnumerate(n, dim, function () { count.c++; },
    maxStates === undefined ? 1e9 : maxStates);
  return count.c;
}

/**
 * Exact exhaustive enumeration of the HP ground state. Returns
 * { energy, coords, states, truncated }: truncated means the state budget
 * was hit so energy is only the best seen so far, not proven optimal.
 */
function runExact(options) {
  var opts = options || {};
  var seq = validateHpSequence(opts.sequence);
  var dim = validateDim(opts.dim === undefined ? 2 : opts.dim);
  var maxStates = opts.maxStates === undefined ? 3000000 : opts.maxStates;
  if (!Number.isInteger(maxStates) || maxStates < 1) {
    throw new RangeError('maxStates must be a positive integer');
  }
  var t0 = Date.now();
  var bestE = Infinity;
  var bestCoords = null;
  var st = canonicalEnumerate(seq.length, dim, function (coords) {
    var e = energyOf(seq, coords, dim);
    if (e < bestE) { bestE = e; bestCoords = coords; }
  }, maxStates);
  if (!Number.isFinite(bestE)) { bestE = 0; bestCoords = straightCoords(seq.length); }
  return {
    algorithm: 'exact',
    label: 'Exact enumeration',
    best: { coords: bestCoords, energy: bestE },
    states: st.states,
    truncated: st.truncated,
    ms: Date.now() - t0
  };
}

// Largest chain length the exact enumerator is allowed to attempt by
// default. 2D square lattice: canonical state count stays small well past
// 12. 3D cubic: about half a million canonical states at length 10.
var EXACT_MAX_LEN = { 2: 16, 3: 10 };
function exactAllowed(n, dim) {
  return n <= EXACT_MAX_LEN[dim];
}

// ---------------------------------------------------------------------------
// Sampling and search algorithms
// ---------------------------------------------------------------------------
function straightCoords(n) {
  var coords = [];
  for (var i = 0; i < n; i++) { coords.push([i, 0, 0]); }
  return coords;
}

/** Grow one random self avoiding walk; restarts on dead ends. */
function randomSaw(n, dim, rng, maxAttempts) {
  var attempts = maxAttempts === undefined ? 10000 : maxAttempts;
  var dirs = DIRS[dim];
  for (var a = 0; a < attempts; a++) {
    var coords = [[0, 0, 0]];
    var occ = new Set([labKey(0, 0, 0)]);
    var ok = true;
    for (var k = 1; k < n; k++) {
      var last = coords[k - 1];
      var free = [];
      for (var d = 0; d < dirs.length; d++) {
        var ax = dirs[d];
        var cell = [last[0] + ax[0], last[1] + ax[1], last[2] + ax[2]];
        if (!occ.has(labKey(cell[0], cell[1], cell[2]))) { free.push(cell); }
      }
      if (!free.length) { ok = false; break; }
      var pick = free[Math.floor(rng() * free.length)];
      coords.push(pick);
      occ.add(labKey(pick[0], pick[1], pick[2]));
    }
    if (ok) { return coords; }
  }
  throw new Error('random SAW growth failed after ' + attempts + ' attempts');
}

function foldState(seq, dim, coords) {
  return { seq: seq, dim: dim, coords: coords, energy: energyOf(seq, coords, dim) };
}

function snapshotMetrics(seq, coords, dim, step) {
  var m = metricsOf(seq, coords, dim);
  m.step = step;
  return m;
}

/**
 * Random SAW sampling baseline: independent growth samples, keep the best.
 * Not a uniform sampler; documented honestly as a random baseline.
 */
function runRandomWalk(options) {
  var opts = options || {};
  var seq = validateHpSequence(opts.sequence);
  var dim = validateDim(opts.dim === undefined ? 2 : opts.dim);
  var samples = opts.samples === undefined ? 250 : opts.samples;
  if (!Number.isInteger(samples) || samples < 1) {
    throw new RangeError('samples must be a positive integer');
  }
  var rng = makeRng(seedOption(opts.seed, 1));
  var t0 = Date.now();
  var best = null;
  var points = [];
  var trajectory = [];
  for (var s = 0; s < samples; s++) {
    var coords = randomSaw(seq.length, dim, rng);
    var m = metricsOf(seq, coords, dim);
    points.push({ E: m.energy, Rg: m.rg, Ree: m.ree, step: s });
    if (!best || m.energy < best.energy) {
      best = { coords: coords, energy: m.energy };
    }
    trajectory.push({ step: s, energy: best.energy });
  }
  return {
    algorithm: 'random-walk',
    label: 'Random SAW sampling',
    best: best,
    points: points,
    trajectory: trajectory,
    samples: samples,
    ms: Date.now() - t0
  };
}

/**
 * Greedy local pivot descent: every iteration scores all single pivot
 * moves and commits the lowest energy one (ties broken by the seeded rng).
 * The walk starts from a seeded random conformation, not the straight
 * chain: on a straight chain no single pivot can ever improve on zero
 * contacts, so the descent would stop before it began. It still stalls at
 * the first local optimum; the teaching point is exactly that this gets
 * stuck where annealing does not.
 */
function runGreedyPivot(options) {
  var opts = options || {};
  var seq = validateHpSequence(opts.sequence);
  var dim = validateDim(opts.dim === undefined ? 2 : opts.dim);
  var maxIter = opts.maxIter === undefined ? 12 * seq.length : opts.maxIter;
  if (!Number.isInteger(maxIter) || maxIter < 0) {
    throw new RangeError('maxIter must be a nonnegative integer');
  }
  var rng = makeRng(seedOption(opts.seed, 1));
  var t0 = Date.now();
  var fold = foldState(seq, dim, randomSaw(seq.length, dim, rng));
  var trajectory = [snapshotMetrics(seq, fold.coords, dim, 0)];
  var iter = 0;
  while (iter < maxIter) {
    var bestCand = [];
    var bestE = fold.energy;
    for (var i = 0; i < seq.length - 1; i++) {
      for (var t = 1; t <= 3; t++) {
        var cand = pivotAt(fold.coords, i, t, dim);
        var built = evaluateConformation(seq, cand, dim);
        if (!built) { continue; }
        if (built.energy < bestE) {
          bestE = built.energy;
          bestCand = [cand];
        } else if (built.energy === bestE && built.energy < fold.energy) {
          bestCand.push(cand);
        }
      }
    }
    if (!bestCand.length) { break; }
    fold.coords = bestCand[Math.floor(rng() * bestCand.length)];
    fold.energy = bestE;
    iter++;
    trajectory.push(snapshotMetrics(seq, fold.coords, dim, iter));
  }
  return {
    algorithm: 'greedy-pivot',
    label: 'Greedy local pivot',
    best: { coords: fold.coords, energy: fold.energy },
    trajectory: trajectory,
    points: trajectory.map(function (m) {
      return { E: m.energy, Rg: m.rg, Ree: m.ree, step: m.step };
    }),
    iterations: iter,
    ms: Date.now() - t0
  };
}

function recordEvery(steps, cap) {
  var every = Math.ceil(steps / (cap === undefined ? 160 : cap));
  return Math.max(1, every);
}

/**
 * Fixed temperature Metropolis Monte Carlo on the lab move set.
 * Records periodic (step, E, R_g, R_ee) samples for the funnel chart.
 */
function runMetropolis(options) {
  var opts = options || {};
  var seq = validateHpSequence(opts.sequence);
  var dim = validateDim(opts.dim === undefined ? 2 : opts.dim);
  var steps = opts.steps === undefined ? 3000 : opts.steps;
  if (!Number.isInteger(steps) || steps < 1) {
    throw new RangeError('steps must be a positive integer');
  }
  var temp = opts.temperature === undefined ? 0.8 : opts.temperature;
  if (!Number.isFinite(temp) || temp < 0) {
    throw new RangeError('temperature must be a nonnegative finite number');
  }
  var rng = makeRng(seedOption(opts.seed, 1));
  var t0 = Date.now();
  var fold = foldState(seq, dim, straightCoords(seq.length));
  var best = { coords: cloneCoords(fold.coords), energy: fold.energy };
  var trajectory = [{ step: 0, energy: fold.energy }];
  var points = [snapshotMetrics(seq, fold.coords, dim, 0)];
  var accepted = 0, moved = 0;
  var every = recordEvery(steps);
  for (var s = 0; s < steps; s++) {
    var res = mcStepLab(fold, rng, temp);
    if (res.moved) { moved++; }
    if (res.accepted) { accepted++; }
    if (fold.energy < best.energy) {
      best = { coords: cloneCoords(fold.coords), energy: fold.energy };
    }
    if ((s + 1) % every === 0 || s === steps - 1) {
      trajectory.push({ step: s + 1, energy: best.energy });
      points.push(snapshotMetrics(seq, fold.coords, dim, s + 1));
    }
  }
  return {
    algorithm: 'metropolis',
    label: 'Metropolis MC (fixed T)',
    temperature: temp,
    best: best,
    trajectory: trajectory,
    points: points,
    accepted: accepted,
    moved: moved,
    steps: steps,
    ms: Date.now() - t0
  };
}

/**
 * Simulated annealing: geometric cooling T0 -> T1 over sweeps of
 * stepsPerSweep Metropolis steps. Same move set and recording as
 * runMetropolis so the two stay directly comparable.
 */
function runAnnealing(options) {
  var opts = options || {};
  var seq = validateHpSequence(opts.sequence);
  var dim = validateDim(opts.dim === undefined ? 2 : opts.dim);
  var sweeps = opts.sweeps === undefined ? 80 : opts.sweeps;
  var stepsPerSweep = opts.stepsPerSweep === undefined ? Math.max(seq.length, 20) : opts.stepsPerSweep;
  if (!Number.isInteger(sweeps) || sweeps < 1) {
    throw new RangeError('sweeps must be a positive integer');
  }
  if (!Number.isInteger(stepsPerSweep) || stepsPerSweep < 1) {
    throw new RangeError('stepsPerSweep must be a positive integer');
  }
  var T0 = opts.T0 === undefined ? 1.5 : opts.T0;
  var T1 = opts.T1 === undefined ? 0.05 : opts.T1;
  if (!Number.isFinite(T0) || T0 <= 0) { throw new RangeError('T0 must be positive'); }
  if (!Number.isFinite(T1) || T1 < 0) { throw new RangeError('T1 must be nonnegative'); }
  var rng = makeRng(seedOption(opts.seed, 1));
  var t0 = Date.now();
  var fold = foldState(seq, dim, straightCoords(seq.length));
  var best = { coords: cloneCoords(fold.coords), energy: fold.energy };
  var trajectory = [{ step: 0, energy: fold.energy }];
  var points = [snapshotMetrics(seq, fold.coords, dim, 0)];
  var accepted = 0, moved = 0;
  var step = 0;
  for (var s = 0; s < sweeps; s++) {
    var frac = sweeps === 1 ? 1 : s / (sweeps - 1);
    var T = (T1 === 0) ? T0 + (T1 - T0) * frac : T0 * Math.pow(T1 / T0, frac);
    for (var st = 0; st < stepsPerSweep; st++) {
      var res = mcStepLab(fold, rng, T);
      if (res.moved) { moved++; }
      if (res.accepted) { accepted++; }
      if (fold.energy < best.energy) {
        best = { coords: cloneCoords(fold.coords), energy: fold.energy };
      }
      step++;
    }
    trajectory.push({ step: step, energy: best.energy });
    points.push(snapshotMetrics(seq, fold.coords, dim, step));
  }
  return {
    algorithm: 'annealing',
    label: 'Simulated annealing',
    T0: T0, T1: T1,
    best: best,
    trajectory: trajectory,
    points: points,
    accepted: accepted,
    moved: moved,
    steps: step,
    ms: Date.now() - t0
  };
}

// ---------------------------------------------------------------------------
// Comparator driver
// ---------------------------------------------------------------------------
/**
 * Run the full algorithm bench on one HP sequence. Options:
 *   sequence, dim, seed, samples, mcSteps, temperature, saSweeps,
 *   saStepsPerSweep, T0, T1, greedyMaxIter, exactMaxStates, includeExact
 * Returns { sequence, dim, seed, results, exact }: results in the order
 * random-walk, greedy-pivot, metropolis, annealing, then exact when run.
 * Each result gains metrics, optimalityGap (when the exact ground state is
 * proven) and a boolean exact flag.
 */
function compareAlgorithms(options) {
  var opts = options || {};
  var seq = validateHpSequence(opts.sequence);
  var dim = validateDim(opts.dim === undefined ? 2 : opts.dim);
  var seed = seedOption(opts.seed, 1);
  var results = [];
  results.push(runRandomWalk({ sequence: seq, dim: dim, seed: seed, samples: opts.samples }));
  results.push(runGreedyPivot({ sequence: seq, dim: dim, seed: (seed + 1) >>> 0, maxIter: opts.greedyMaxIter }));
  results.push(runMetropolis({ sequence: seq, dim: dim, seed: (seed + 2) >>> 0, steps: opts.mcSteps, temperature: opts.temperature }));
  results.push(runAnnealing({
    sequence: seq, dim: dim, seed: (seed + 3) >>> 0,
    sweeps: opts.saSweeps, stepsPerSweep: opts.saStepsPerSweep,
    T0: opts.T0, T1: opts.T1
  }));
  var exact = null;
  var wantExact = opts.includeExact !== false && exactAllowed(seq.length, dim);
  if (wantExact) {
    exact = runExact({ sequence: seq, dim: dim, maxStates: opts.exactMaxStates });
    results.push(exact);
  }
  for (var i = 0; i < results.length; i++) {
    var r = results[i];
    r.metrics = metricsOf(seq, r.best.coords, dim);
    r.optimalityGap = (exact && !exact.truncated) ? (r.best.energy - exact.best.energy) : null;
    r.exact = (r === exact);
  }
  return {
    sequence: seq,
    dim: dim,
    seed: seed,
    results: results,
    exact: exact,
    exactSkipped: !wantExact
  };
}

/**
 * Cloud of random SAW metric points used as the folding funnel backdrop.
 * Returns { points: [{ E, Rg, Ree }], minEnergy, maxEnergy }.
 */
function sampleFunnel(seq, dim, count, seed) {
  var s = validateHpSequence(seq);
  var d = validateDim(dim);
  if (!Number.isInteger(count) || count < 1) {
    throw new RangeError('count must be a positive integer');
  }
  var rng = makeRng(seedOption(seed, 1));
  var points = [];
  var lo = Infinity, hi = -Infinity;
  for (var i = 0; i < count; i++) {
    var coords = randomSaw(s.length, d, rng);
    var m = metricsOf(s, coords, d);
    points.push({ E: m.energy, Rg: m.rg, Ree: m.ree });
    if (m.energy < lo) { lo = m.energy; }
    if (m.energy > hi) { hi = m.energy; }
  }
  return { points: points, minEnergy: lo, maxEnergy: hi };
}

// ---------------------------------------------------------------------------
// Mutation stress scan
// ---------------------------------------------------------------------------
function flipHp(ch) {
  return ch === 'H' ? 'P' : 'H';
}

/**
 * Point mutation stress scan: at each residue, flip H <-> P on the SAME
 * conformation and rescore. deltaE > 0 means the mutant is less stable
 * (a hydrophobic core contact lost); deltaE < 0 means a stabilizing swap.
 * Returns { e0, perResidue, destabilizing, stabilizing, neutral, worst }.
 */
function mutationScan(seq, coords, dim) {
  var s = validateHpSequence(seq);
  if (!isValidSaw(coords, dim) || coords.length !== s.length) {
    throw new Error('coords must be a valid self avoiding walk matching the sequence length');
  }
  var e0 = energyOf(s, coords, dim);
  var per = [];
  var destabilizing = 0, stabilizing = 0, neutral = 0;
  var worst = null;
  for (var i = 0; i < s.length; i++) {
    var mut = s.slice(0, i) + flipHp(s[i]) + s.slice(i + 1);
    var em = energyOf(mut, coords, dim);
    var dE = em - e0;
    if (dE > 0) { destabilizing++; }
    else if (dE < 0) { stabilizing++; }
    else { neutral++; }
    var entry = { index: i, from: s[i], to: flipHp(s[i]), dE: dE };
    per.push(entry);
    if (!worst || dE > worst.dE) { worst = entry; }
  }
  return {
    e0: e0,
    perResidue: per,
    destabilizing: destabilizing,
    stabilizing: stabilizing,
    neutral: neutral,
    worst: worst
  };
}

// ---------------------------------------------------------------------------
// Contact matrix
// ---------------------------------------------------------------------------
var CMAP_FAR = 0;        // not lattice adjacent
var CMAP_BONDED = 1;     // chain neighbors (j = i +/- 1)
var CMAP_HH = 2;         // non bonded adjacent H-H contact (native contact)
var CMAP_OTHER = 3;      // non bonded adjacent pair, at least one P
var CMAP_SELF = 4;       // diagonal

/**
 * L x L contact matrix as a flat array of CMAP_* codes.
 * { n, cells, hhCount, contactCount }
 */
function contactMatrix(seq, coords, dim) {
  var s = validateHpSequence(seq);
  if (!isValidSaw(coords, dim) || coords.length !== s.length) {
    throw new Error('coords must be a valid self avoiding walk matching the sequence length');
  }
  var n = coords.length;
  var cells = new Array(n * n);
  for (var i = 0; i < n * n; i++) { cells[i] = CMAP_FAR; }
  for (var k = 0; k < n; k++) { cells[k * n + k] = CMAP_SELF; }
  var hh = 0;
  var cc = 0;
  for (var a = 0; a < n; a++) {
    for (var b = a + 1; b < n; b++) {
      if (manhattan(coords[a], coords[b]) !== 1) { continue; }
      var bonded = (b === a + 1);
      var code;
      if (bonded) { code = CMAP_BONDED; }
      else {
        cc++;
        var isHh = s[a] === 'H' && s[b] === 'H';
        if (isHh) { hh++; }
        code = isHh ? CMAP_HH : CMAP_OTHER;
      }
      cells[a * n + b] = code;
      cells[b * n + a] = code;
    }
  }
  return { n: n, cells: cells, hhCount: hh, contactCount: cc };
}

// ---------------------------------------------------------------------------
// FIPS 180-4 SHA-256, pure JavaScript, unkeyed digest.
// ---------------------------------------------------------------------------
function utf8Bytes(str) {
  var out = [];
  for (var i = 0; i < str.length; i++) {
    var cp = str.charCodeAt(i);
    if (cp >= 0xD800 && cp <= 0xDBFF && i + 1 < str.length) {
      var lo = str.charCodeAt(i + 1);
      if (lo >= 0xDC00 && lo <= 0xDFFF) {
        cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
        i++;
      }
    }
    if (cp < 0x80) { out.push(cp); }
    else if (cp < 0x800) {
      out.push(0xC0 | (cp >>> 6), 0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
      out.push(0xE0 | (cp >>> 12), 0x80 | ((cp >>> 6) & 0x3F), 0x80 | (cp & 0x3F));
    } else {
      out.push(0xF0 | (cp >>> 18), 0x80 | ((cp >>> 12) & 0x3F),
        0x80 | ((cp >>> 6) & 0x3F), 0x80 | (cp & 0x3F));
    }
  }
  return out;
}

var SHA256_H = [
  0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
  0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
];

var SHA256_K = [
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
  0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
  0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
  0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
  0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
  0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
  0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
  0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
  0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
];

function sha256Hex(str) {
  var bytes = utf8Bytes(str);
  var bitLen = bytes.length * 8;
  // Padding: append 0x80, then zeros to 56 mod 64, then 64-bit length.
  bytes.push(0x80);
  while (bytes.length % 64 !== 56) { bytes.push(0); }
  var hi = Math.floor(bitLen / 0x100000000);
  var lo = bitLen >>> 0;
  bytes.push((hi >>> 24) & 255, (hi >>> 16) & 255, (hi >>> 8) & 255, hi & 255);
  bytes.push((lo >>> 24) & 255, (lo >>> 16) & 255, (lo >>> 8) & 255, lo & 255);

  var h = SHA256_H.slice();
  var w = new Array(64);
  function rotr(x, n) { return (x >>> n) | (x << (32 - n)); }
  for (var off = 0; off < bytes.length; off += 64) {
    for (var t = 0; t < 16; t++) {
      w[t] = ((bytes[off + t * 4] << 24) | (bytes[off + t * 4 + 1] << 16) |
        (bytes[off + t * 4 + 2] << 8) | bytes[off + t * 4 + 3]) >>> 0;
    }
    for (var u = 16; u < 64; u++) {
      var s0 = rotr(w[u - 15], 7) ^ rotr(w[u - 15], 18) ^ (w[u - 15] >>> 3);
      var s1 = rotr(w[u - 2], 17) ^ rotr(w[u - 2], 19) ^ (w[u - 2] >>> 10);
      w[u] = (w[u - 16] + s0 + w[u - 7] + s1) >>> 0;
    }
    var a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], g = h[5], hh = h[6], i2 = h[7];
    for (var r = 0; r < 64; r++) {
      var S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      var ch = (e & g) ^ (~e & hh);
      var t1 = (i2 + S1 + ch + SHA256_K[r] + w[r]) >>> 0;
      var S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      var maj = (a & b) ^ (a & c) ^ (b & c);
      var t2 = (S0 + maj) >>> 0;
      i2 = hh; hh = g; g = e; e = (d + t1) >>> 0;
      d = c; c = b; b = a; a = (t1 + t2) >>> 0;
    }
    h[0] = (h[0] + a) >>> 0; h[1] = (h[1] + b) >>> 0;
    h[2] = (h[2] + c) >>> 0; h[3] = (h[3] + d) >>> 0;
    h[4] = (h[4] + e) >>> 0; h[5] = (h[5] + g) >>> 0;
    h[6] = (h[6] + hh) >>> 0; h[7] = (h[7] + i2) >>> 0;
  }
  var out = '';
  for (var q = 0; q < 8; q++) {
    out += ('0000000' + h[q].toString(16)).slice(-8);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Canonical JSON + checksummed audit reports
// ---------------------------------------------------------------------------
function canonicalize(value) {
  if (value === null || typeof value !== 'object') {
    return JSON.stringify(value === undefined ? null : value);
  }
  if (Array.isArray(value)) {
    return '[' + value.map(canonicalize).join(',') + ']';
  }
  var keys = Object.keys(value).sort();
  var parts = [];
  for (var i = 0; i < keys.length; i++) {
    parts.push(JSON.stringify(keys[i]) + ':' + canonicalize(value[keys[i]]));
  }
  return '{' + parts.join(',') + '}';
}

function checksumPayload(payload) {
  return sha256Hex(canonicalize(payload));
}

function comparisonAuditPayload(cmp, scan) {
  if (!cmp || typeof cmp !== 'object' || !Array.isArray(cmp.results)) {
    throw new Error('comparison audit needs a comparator result object');
  }
  var results = cmp.results.map(function (r) {
    if (!r || !r.best || !r.metrics) {
      throw new Error('comparison audit result is missing best metrics');
    }
    return {
      algorithm: r.algorithm,
      label: r.label,
      energy: r.best.energy,
      contactCount: r.metrics.contactCount,
      hhContacts: r.metrics.hhContacts,
      rg: Math.round(r.metrics.rg * 10000) / 10000,
      ree: Math.round(r.metrics.ree * 10000) / 10000,
      optimalityGap: r.optimalityGap === undefined ? null : r.optimalityGap,
      states: r.states === undefined ? null : r.states
    };
  });
  return {
    sequence: cmp.sequence,
    dimension: cmp.dim,
    seed: cmp.seed,
    groundEnergy: cmp.exact && !cmp.exact.truncated ? cmp.exact.best.energy : null,
    results: results,
    mutationScan: scan ? {
      algorithm: scan.algorithm,
      e0: scan.e0,
      destabilizing: scan.destabilizing,
      stabilizing: scan.stabilizing,
      neutral: scan.neutral,
      worst: scan.worst,
      perResidue: scan.perResidue
    } : null,
    modelDisclosure: 'Dill HP lattice model: simplified statistical physics toy, not atomistic prediction',
    timingDisclosure: 'wall clock timings are excluded so this audit is deterministic for fixed inputs'
  };
}

/**
 * Build a checksummed audit report. The checksum covers the canonicalized
 * payload only; it is an unkeyed tamper detection digest, not a signature.
 */
function buildAuditReport(payload) {
  return {
    tool: 'fold-lattice',
    toolVersion: F.VERSION,
    labVersion: LAB_VERSION,
    kind: 'hp-folding-audit',
    checksumAlgorithm: 'SHA-256 (FIPS 180-4, unkeyed digest)',
    checksum: checksumPayload(payload),
    payload: payload
  };
}

function verifyAuditReport(report) {
  if (!report || typeof report !== 'object') {
    return { ok: false, reason: 'report is not an object' };
  }
  if (typeof report.checksum !== 'string' || !/^[0-9a-f]{64}$/.test(report.checksum)) {
    return { ok: false, reason: 'missing or malformed checksum' };
  }
  if (!Object.prototype.hasOwnProperty.call(report, 'payload') ||
      !report.payload || typeof report.payload !== 'object') {
    return { ok: false, reason: 'missing or malformed payload' };
  }
  var actual = checksumPayload(report.payload);
  if (actual !== report.checksum) {
    return { ok: false, reason: 'checksum mismatch', expected: report.checksum, actual: actual };
  }
  return { ok: true, reason: 'checksum verified' };
}

function auditToJson(report) {
  return JSON.stringify(report, null, 2) + '\n';
}

function verifyAuditText(text) {
  var parsed;
  try {
    parsed = JSON.parse(text);
  } catch (err) {
    return { ok: false, reason: 'not valid JSON: ' + err.message };
  }
  return verifyAuditReport(parsed);
}

/**
 * CSV cell escaping with spreadsheet formula injection protection: values
 * that would be read as a formula (= + - @ or leading whitespace control
 * characters) get a leading single quote, and anything holding a comma,
 * quote, newline or leading/trailing space is double quoted with internal
 * quotes doubled.
 */
function csvCell(value) {
  var s = value === null || value === undefined ? '' : String(value);
  // A leading formula character makes a dangerous cell, except on a plain
  // number like -5 or +3.5 which a spreadsheet reads as a value anyway.
  var plainNumber = /^[+-]?\d+(\.\d+)?([eE][+-]?\d+)?$/.test(s);
  var danger = !plainNumber && /^[=+\-@\t\v\f\r ]/.test(s);
  if (danger) { s = "'" + s; }
  if (/[",\r\n]/.test(s) || danger || s !== s.trim()) {
    return '"' + s.replace(/"/g, '""') + '"';
  }
  return s;
}

function auditToCsv(report) {
  var p = report.payload || {};
  var lines = [];
  lines.push('section,' + ['field', 'value'].join(','));
  lines.push(['meta', 'tool', csvCell(report.tool)].join(','));
  lines.push(['meta', 'toolVersion', csvCell(report.toolVersion)].join(','));
  lines.push(['meta', 'labVersion', csvCell(report.labVersion)].join(','));
  lines.push(['meta', 'sequence', csvCell(p.sequence)].join(','));
  lines.push(['meta', 'dimension', csvCell(p.dimension)].join(','));
  lines.push(['meta', 'seed', csvCell(p.seed)].join(','));
  if (p.groundEnergy !== undefined && p.groundEnergy !== null) {
    lines.push(['meta', 'groundEnergy', csvCell(p.groundEnergy)].join(','));
  }
  lines.push(['meta', 'checksumAlgorithm', csvCell(report.checksumAlgorithm)].join(','));
  lines.push(['meta', 'checksum', csvCell(report.checksum)].join(','));
  if (p.timingDisclosure) { lines.push(['meta', 'timingDisclosure', csvCell(p.timingDisclosure)].join(',')); }
  lines.push('algorithm,best_energy,contact_count,hh_contacts,radius_of_gyration,end_to_end,optimality_gap');
  var results = p.results || [];
  for (var i = 0; i < results.length; i++) {
    var r = results[i];
    lines.push([
      csvCell(r.algorithm), csvCell(r.energy), csvCell(r.contactCount),
      csvCell(r.hhContacts), csvCell(r.rg), csvCell(r.ree),
      csvCell(r.optimalityGap === null || r.optimalityGap === undefined ? 'n/a' : r.optimalityGap)
    ].join(','));
  }
  var mut = p.mutationScan;
  if (mut && mut.perResidue) {
    lines.push('mutation_index,from,to,delta_e');
    for (var m = 0; m < mut.perResidue.length; m++) {
      var e = mut.perResidue[m];
      lines.push([csvCell(e.index), csvCell(e.from), csvCell(e.to), csvCell(e.dE)].join(','));
    }
  }
  return lines.join('\n') + '\n';
}

function mdEscape(s) {
  return String(s).replace(/\|/g, '\\|');
}

function auditToMarkdown(report) {
  var p = report.payload || {};
  var out = [];
  out.push('# fold-lattice folding audit');
  out.push('');
  out.push('- sequence: `' + mdEscape(p.sequence || '') + '` (' +
    (p.dimension === 3 ? '3D cubic' : '2D square') + ' lattice, Dill HP model)');
  out.push('- seed: ' + (p.seed === undefined ? 'n/a' : p.seed));
  if (p.groundEnergy !== undefined && p.groundEnergy !== null) {
    out.push('- exact ground state energy: ' + p.groundEnergy);
  }
  out.push('- checksum (' + mdEscape(report.checksumAlgorithm) + '):');
  out.push('  `' + report.checksum + '`');
  out.push('');
  out.push('| algorithm | best E | contacts | H-H | R_g | R_ee | optimality gap |');
  out.push('|---|---|---|---|---|---|---|');
  var results = p.results || [];
  for (var i = 0; i < results.length; i++) {
    var r = results[i];
    var gap = r.optimalityGap === null || r.optimalityGap === undefined ? 'n/a' : String(r.optimalityGap);
    out.push('| ' + mdEscape(r.algorithm) + ' | ' + r.energy + ' | ' + r.contactCount +
      ' | ' + r.hhContacts + ' | ' + r.rg + ' | ' + r.ree + ' | ' + gap + ' |');
  }
  var mut = p.mutationScan;
  if (mut && mut.perResidue) {
    out.push('');
    out.push('## Point mutation scan (H <-> P on the best conformation)');
    out.push('');
    out.push('destabilizing: ' + mut.destabilizing + ', stabilizing: ' + mut.stabilizing +
      ', neutral: ' + mut.neutral + ', worst: residue ' +
      (mut.worst ? (mut.worst.index + 1) + ' (' + mut.worst.from + ' to ' + mut.worst.to +
        ', deltaE ' + mut.worst.dE + ')' : 'n/a'));
  }
  out.push('');
  out.push('The Dill HP lattice model is a simplified statistical physics toy model.');
  if (p.timingDisclosure) { out.push(p.timingDisclosure + '.'); }
  return out.join('\n') + '\n';
}

/** Flip one hexadecimal digit of the checksum in an audit JSON text. */
function tamperAuditText(text) {
  var i = text.search(/[0-9a-f]{64}/);
  if (i < 0) { return text; }
  var ch = text[i];
  var rep = ch === '0' ? '1' : '0';
  return text.slice(0, i) + rep + text.slice(i + 1);
}

// ---------------------------------------------------------------------------
// URL hash bookmarking: #fold=SEQ:DIMd:SEED
// ---------------------------------------------------------------------------
var BOOKMARK_MAX_SEQ = 120;

function encodeBookmark(opts) {
  var seq = validateHpSequence(opts.sequence);
  if (seq.length > BOOKMARK_MAX_SEQ) {
    throw new Error('sequence too long for a bookmark (max ' + BOOKMARK_MAX_SEQ + ')');
  }
  var dim = validateDim(opts.dim === undefined ? 2 : opts.dim);
  var seed = seedOption(opts.seed, 1);
  return 'fold=' + seq + ':' + dim + 'd:' + seed;
}

function decodeBookmark(str) {
  if (typeof str !== 'string') { return null; }
  var s = str.trim().replace(/^#+/, '');
  var m = /^fold=([HPhp]+):([23])d:(\d{1,10})$/.exec(s);
  if (!m || m[1].length > BOOKMARK_MAX_SEQ) { return null; }
  var seed = parseInt(m[3], 10);
  if (!Number.isFinite(seed) || seed > 4294967295) { return null; }
  return {
    sequence: m[1].toUpperCase(),
    dim: parseInt(m[2], 10),
    seed: seed
  };
}

// ---------------------------------------------------------------------------
// HP <-> AA bridge into the 3D viewer (H -> A hydrophobic, P -> S polar).
// ---------------------------------------------------------------------------
function hpToAa(seq) {
  return validateHpSequence(seq).split('').map(function (ch) {
    return ch === 'H' ? 'A' : 'S';
  }).join('');
}

// ---------------------------------------------------------------------------
// In-page self tests: fast invariant checks a visitor can run in the page.
// ---------------------------------------------------------------------------
function runSelfTests() {
  var results = [];
  function check(name, fn) {
    try {
      var detail = fn();
      results.push({ name: name, pass: true, detail: detail || 'ok' });
    } catch (err) {
      results.push({ name: name, pass: false, detail: String(err && err.message ? err.message : err) });
    }
  }
  function eq(a, b, msg) {
    if (a !== b) { throw new Error((msg || 'mismatch') + ': ' + a + ' !== ' + b); }
  }
  function ok(v, msg) { if (!v) { throw new Error(msg || 'assertion failed'); } }

  check('SHA-256 known answer: empty string', function () {
    eq(sha256Hex(''), 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855');
    return 'FIPS 180-4 vector matches';
  });
  check('SHA-256 known answer: abc', function () {
    eq(sha256Hex('abc'), 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
    return 'FIPS 180-4 vector matches';
  });
  check('SHA-256 known answer: quick brown fox', function () {
    eq(sha256Hex('The quick brown fox jumps over the lazy dog'),
      'd7a8fbb307d7809469ca9abcb0082e4f8d5651e46d3cdb762d02d0bf37c9e592');
    return 'FIPS 180-4 vector matches';
  });

  check('self avoiding walk invariants on all algorithms', function () {
    var cmp = compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 7, samples: 40, mcSteps: 400, saSweeps: 20 });
    for (var i = 0; i < cmp.results.length; i++) {
      ok(isValidSaw(cmp.results[i].best.coords, 2), cmp.results[i].algorithm + ' best is not a SAW');
    }
    return cmp.results.length + ' algorithms returned valid walks';
  });

  check('HP energy is non positive and equals minus H-H contacts', function () {
    var rng = makeRng(3);
    for (var i = 0; i < 25; i++) {
      var c = randomSaw(10, 2, rng);
      var e = energyOf('HPPHPHPHPH', c, 2);
      ok(e <= 0, 'positive energy');
      eq(e, -hhPairs('HPPHPHPHPH', c, 2).length, 'energy mismatch');
    }
    return '25 random walks checked';
  });

  check('energy is invariant under chain reversal', function () {
    var rng = makeRng(5);
    for (var i = 0; i < 15; i++) {
      var c = randomSaw(9, 2, rng);
      var seq = 'HPPHPHPHP';
      var rev = seq.split('').reverse().join('');
      eq(energyOf(seq, c, 2), energyOf(rev, cloneCoords(c).reverse(), 2), 'reversal changed energy');
    }
    return 'reversal symmetry holds';
  });

  check('exact enumeration beats or matches every heuristic', function () {
    var cmp = compareAlgorithms({ sequence: 'HPPHPHPH', dim: 2, seed: 11, samples: 40, mcSteps: 400, saSweeps: 20 });
    ok(cmp.exact && !cmp.exact.truncated, 'exact run missing');
    for (var i = 0; i < cmp.results.length; i++) {
      ok(cmp.results[i].best.energy >= cmp.exact.best.energy,
        cmp.results[i].algorithm + ' beat the exact ground state');
    }
    return 'ground state ' + cmp.exact.best.energy + ' verified over ' + cmp.exact.states + ' states';
  });

  check('radius of gyration and end to end bounds', function () {
    var rng = makeRng(9);
    var rgExtended = F.radiusOfGyration(straightCoords(12));
    for (var i = 0; i < 20; i++) {
      var c = randomSaw(12, 2, rng);
      var rg = F.radiusOfGyration(c);
      ok(rg >= 0, 'negative Rg');
      ok(rg <= rgExtended + 1e-9, 'folded Rg exceeds extended');
      var ree = endToEnd(c);
      ok(ree >= 1 - 1e-9 && ree <= 11 + 1e-9, 'R_ee out of bounds');
    }
    return 'Rg and R_ee inside bounds for 20 walks';
  });

  check('contact matrix symmetry and H-H count', function () {
    var rng = makeRng(13);
    var c = randomSaw(10, 2, rng);
    var seq = 'HPPHPHPHPH';
    var m = contactMatrix(seq, c, 2);
    var hh = 0;
    for (var a = 0; a < m.n; a++) {
      for (var b = 0; b < m.n; b++) {
        eq(m.cells[a * m.n + b], m.cells[b * m.n + a], 'matrix not symmetric');
        if (m.cells[a * m.n + b] === CMAP_HH) { hh++; }
      }
    }
    eq(hh / 2, -energyOf(seq, c, 2), 'H-H count disagrees with energy');
    return 'matrix symmetric, ' + (hh / 2) + ' H-H contacts';
  });

  check('CSV formula injection cells are escaped', function () {
    var cell = csvCell('=1+1');
    ok(cell.charAt(0) === '"' && cell.indexOf("'") === 1, 'dangerous cell not quoted');
    ok(csvCell('+cmd').indexOf("'") >= 0, 'plus cell not escaped');
    return 'formula cells neutralized';
  });

  check('checksummed report verifies and tampering is detected', function () {
    var report = buildAuditReport({ sequence: 'HPPHPH', dimension: 2, results: [{ algorithm: 'x', energy: -1 }] });
    var text = auditToJson(report);
    ok(verifyAuditText(text).ok, 'fresh report failed verification');
    ok(!verifyAuditText(tamperAuditText(text)).ok, 'tampered report still verified');
    return 'verify + tamper detection both pass';
  });

  check('URL hash bookmark round trips', function () {
    var mark = encodeBookmark({ sequence: 'hpph', dim: 2, seed: 42 });
    var back = decodeBookmark('#' + mark);
    eq(back.sequence, 'HPPH'); eq(back.dim, 2); eq(back.seed, 42);
    ok(decodeBookmark('#fold=BAD!') === null, 'junk bookmark accepted');
    return mark;
  });

  check('mutation scan conservation and all P sanity', function () {
    var seq = 'HPHPHPHP';
    var c = randomSaw(8, 2, makeRng(17));
    var scan = mutationScan(seq, c, 2);
    eq(scan.perResidue.length, seq.length);
    var allP = 'PPPPPPPP';
    eq(energyOf(allP, c, 2), 0, 'all P chain has nonzero energy');
    eq(scan.destabilizing + scan.stabilizing + scan.neutral, seq.length, 'categories do not sum');
    return 'worst deltaE ' + (scan.worst ? scan.worst.dE : 'n/a');
  });

  check('seeded runs are deterministic', function () {
    var a = runAnnealing({ sequence: 'HPPHPHPH', dim: 2, seed: 21, sweeps: 15, stepsPerSweep: 10 });
    var b = runAnnealing({ sequence: 'HPPHPHPH', dim: 2, seed: 21, sweeps: 15, stepsPerSweep: 10 });
    eq(a.best.energy, b.best.energy, 'same seed, different energy');
    return 'seed 21 reproducible, E ' + a.best.energy;
  });

  var passed = results.filter(function (r) { return r.pass; }).length;
  return { results: results, passed: passed, total: results.length };
}

// ---------------------------------------------------------------------------
var lab = {
  LAB_VERSION: LAB_VERSION,
  DIRS: DIRS,
  validateHpSequence: validateHpSequence,
  validateDim: validateDim,
  isValidSaw: isValidSaw,
  contactPairs: contactPairs,
  hhPairs: hhPairs,
  energyOf: energyOf,
  endToEnd: endToEnd,
  metricsOf: metricsOf,
  evaluateConformation: evaluateConformation,
  pivotAt: pivotAt,
  pivotCandidate: pivotCandidate,
  cornerCandidate: cornerCandidate,
  crankshaftCandidate: crankshaftCandidate,
  endFlipCandidate: endFlipCandidate,
  MOVE_NAMES: MOVE_NAMES,
  moveCandidate: moveCandidate,
  mcStepLab: mcStepLab,
  canonicalEnumerate: canonicalEnumerate,
  countSaws: countSaws,
  runExact: runExact,
  EXACT_MAX_LEN: EXACT_MAX_LEN,
  exactAllowed: exactAllowed,
  BOOKMARK_MAX_SEQ: BOOKMARK_MAX_SEQ,
  straightCoords: straightCoords,
  randomSaw: randomSaw,
  runRandomWalk: runRandomWalk,
  runGreedyPivot: runGreedyPivot,
  runMetropolis: runMetropolis,
  runAnnealing: runAnnealing,
  compareAlgorithms: compareAlgorithms,
  sampleFunnel: sampleFunnel,
  flipHp: flipHp,
  mutationScan: mutationScan,
  CMAP_FAR: CMAP_FAR,
  CMAP_BONDED: CMAP_BONDED,
  CMAP_HH: CMAP_HH,
  CMAP_OTHER: CMAP_OTHER,
  CMAP_SELF: CMAP_SELF,
  contactMatrix: contactMatrix,
  utf8Bytes: utf8Bytes,
  sha256Hex: sha256Hex,
  canonicalize: canonicalize,
  checksumPayload: checksumPayload,
  comparisonAuditPayload: comparisonAuditPayload,
  buildAuditReport: buildAuditReport,
  verifyAuditReport: verifyAuditReport,
  auditToJson: auditToJson,
  verifyAuditText: verifyAuditText,
  csvCell: csvCell,
  auditToCsv: auditToCsv,
  auditToMarkdown: auditToMarkdown,
  tamperAuditText: tamperAuditText,
  encodeBookmark: encodeBookmark,
  decodeBookmark: decodeBookmark,
  hpToAa: hpToAa,
  runSelfTests: runSelfTests
};

/* global module */
if (typeof module === 'object' && module.exports) {
  module.exports = lab;
} else {
  /* global window */
  window.foldlLab = lab;
}
