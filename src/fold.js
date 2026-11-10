'use strict';

/**
 * fold-lattice: a zero-dependency structural bioinformatics and statistical
 * mechanics toolkit for the browser and Node.js.
 *
 * Contents:
 *   - 3D cubic lattice protein folder: self-avoiding walk chains scored with
 *     either the Dill HP model or the Miyazawa-Jernigan (1996) contact
 *     potential, moved by pull, pivot, crankshaft, corner flip and endpoint
 *     pull moves, and searched with Metropolis simulated annealing and
 *     parallel-tempering replica-exchange Monte Carlo (REMC).
 *   - All-atom backbone reconstruction (N, CA, C, O) via NeRF placement,
 *     Ramachandran (phi, psi) region classification, omega cis/trans checks,
 *     a steric clash detector, radius of gyration, and RCSB PDB parsing and
 *     export.
 *
 * Contact energies (RT units) come from Table 3 of:
 *   Miyazawa, S. and Jernigan, R.L. (1996)
 *   "Residue-residue potentials with a favorable contact pair term and an
 *   unfavorable high packing density term, for simulation and threading."
 *   Journal of Molecular Biology 256, 623-644.
 */

var VERSION = '1.0.0';

// ---------------------------------------------------------------------------
// Seeded random number generator (mulberry32). Deterministic across runs.
// ---------------------------------------------------------------------------
function makeRng(seed) {
  var t = seed >>> 0;
  return function rng() {
    t = (t + 0x6D2B79F5) | 0;
    var z = t;
    z = Math.imul(z ^ (z >>> 15), z | 1);
    z ^= z + Math.imul(z ^ (z >>> 7), z | 61);
    return ((z ^ (z >>> 14)) >>> 0) / 4294967296;
  };
}

// ---------------------------------------------------------------------------
// Small vector helpers. Points are length-3 arrays [x, y, z].
// ---------------------------------------------------------------------------
function vadd(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
function vsub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function vscale(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }
function vdot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function vcross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}
function vlen(a) { return Math.sqrt(vdot(a, a)); }
function vdist(a, b) { return vlen(vsub(a, b)); }
function vnorm(a) {
  var l = vlen(a);
  if (l === 0) { throw new Error('cannot normalize zero vector'); }
  return vscale(a, 1 / l);
}
function veq(a, b) { return a[0] === b[0] && a[1] === b[1] && a[2] === b[2]; }

var AXIS_VECTORS = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];

function isUnitAxis(v) {
  var s = Math.abs(v[0]) + Math.abs(v[1]) + Math.abs(v[2]);
  return s === 1;
}

// ---------------------------------------------------------------------------
// Amino acid data
// ---------------------------------------------------------------------------
var AA20 = 'ARNDCQEGHILKMFPSTWYV';

var AA1_TO_3 = {
  A: 'ALA', R: 'ARG', N: 'ASN', D: 'ASP', C: 'CYS', Q: 'GLN', E: 'GLU',
  G: 'GLY', H: 'HIS', I: 'ILE', L: 'LEU', K: 'LYS', M: 'MET', F: 'PHE',
  P: 'PRO', S: 'SER', T: 'THR', W: 'TRP', Y: 'TYR', V: 'VAL'
};

var AA3_TO_1 = (function () {
  var m = {};
  for (var k in AA1_TO_3) { m[AA1_TO_3[k]] = k; }
  return m;
})();

// Miyazawa-Jernigan matrix row/column order (single letter codes).
var MJ_ORDER = 'CMFILVWYAGTSNQDEHRKP';

var MJ_INDEX = (function () {
  var m = {};
  for (var i = 0; i < MJ_ORDER.length; i++) { m[MJ_ORDER[i]] = i; }
  return m;
})();

/**
 * e_ij contact energies in RT units, upper half (including diagonal) of
 * Table 3 of Miyazawa and Jernigan (1996), mirrored to a full symmetric
 * 20x20 matrix. Row/column order is MJ_ORDER (C, M, F, I, L, V, W, Y, A,
 * G, T, S, N, Q, D, E, H, R, K, P).
 */
var MJ1996 = [
  [-5.44, -4.99, -5.80, -5.50, -5.83, -4.96, -4.95, -4.16, -3.57, -3.16, -3.11, -2.86, -2.59, -2.85, -2.41, -2.27, -3.60, -2.57, -1.95, -3.07],
  [-4.99, -5.46, -6.56, -6.02, -6.41, -5.32, -5.55, -4.91, -3.94, -3.39, -3.51, -3.03, -2.95, -3.30, -2.57, -2.89, -3.98, -3.12, -2.48, -3.45],
  [-5.80, -6.56, -7.26, -6.84, -7.28, -6.29, -6.16, -5.66, -4.81, -4.13, -4.28, -4.02, -3.75, -4.10, -3.48, -3.56, -4.77, -3.98, -3.36, -4.25],
  [-5.50, -6.02, -6.84, -6.54, -7.04, -6.05, -5.78, -5.25, -4.58, -3.78, -4.03, -3.52, -3.24, -3.67, -3.17, -3.27, -4.14, -3.63, -3.01, -3.76],
  [-5.83, -6.41, -7.28, -7.04, -7.37, -6.48, -6.14, -5.67, -4.91, -4.16, -4.34, -3.92, -3.74, -4.04, -3.40, -3.59, -4.54, -4.03, -3.37, -4.20],
  [-4.96, -5.32, -6.29, -6.05, -6.48, -5.52, -5.18, -4.62, -4.04, -3.38, -3.46, -3.05, -2.83, -3.07, -2.48, -2.67, -3.58, -3.07, -2.49, -3.32],
  [-4.95, -5.55, -6.16, -5.78, -6.14, -5.18, -5.06, -4.66, -3.82, -3.42, -3.22, -2.99, -3.07, -3.11, -2.84, -2.99, -3.98, -3.41, -2.69, -3.73],
  [-4.16, -4.91, -5.66, -5.25, -5.67, -4.62, -4.66, -4.17, -3.36, -3.01, -3.01, -2.78, -2.76, -2.97, -2.76, -2.79, -3.52, -3.16, -2.60, -3.19],
  [-3.57, -3.94, -4.81, -4.58, -4.91, -4.04, -3.82, -3.36, -2.72, -2.31, -2.32, -2.01, -1.84, -1.89, -1.70, -1.51, -2.41, -1.83, -1.31, -2.03],
  [-3.16, -3.39, -4.13, -3.78, -4.16, -3.38, -3.42, -3.01, -2.31, -2.24, -2.08, -1.82, -1.74, -1.66, -1.59, -1.22, -2.15, -1.72, -1.15, -1.87],
  [-3.11, -3.51, -4.28, -4.03, -4.34, -3.46, -3.22, -3.01, -2.32, -2.08, -2.12, -1.96, -1.88, -1.90, -1.80, -1.74, -2.42, -1.90, -1.31, -1.90],
  [-2.86, -3.03, -4.02, -3.52, -3.92, -3.05, -2.99, -2.78, -2.01, -1.82, -1.96, -1.67, -1.58, -1.49, -1.63, -1.48, -2.11, -1.62, -1.05, -1.57],
  [-2.59, -2.95, -3.75, -3.24, -3.74, -2.83, -3.07, -2.76, -1.84, -1.74, -1.88, -1.58, -1.68, -1.71, -1.68, -1.51, -2.08, -1.64, -1.21, -1.53],
  [-2.85, -3.30, -4.10, -3.67, -4.04, -3.07, -3.11, -2.97, -1.89, -1.66, -1.90, -1.49, -1.71, -1.54, -1.46, -1.42, -1.98, -1.80, -1.29, -1.73],
  [-2.41, -2.57, -3.48, -3.17, -3.40, -2.48, -2.84, -2.76, -1.70, -1.59, -1.80, -1.63, -1.68, -1.46, -1.21, -1.02, -2.32, -2.29, -1.68, -1.33],
  [-2.27, -2.89, -3.56, -3.27, -3.59, -2.67, -2.99, -2.79, -1.51, -1.22, -1.74, -1.48, -1.51, -1.42, -1.02, -0.91, -2.15, -2.27, -1.80, -1.26],
  [-3.60, -3.98, -4.77, -4.14, -4.54, -3.58, -3.98, -3.52, -2.41, -2.15, -2.42, -2.11, -2.08, -1.98, -2.32, -2.15, -3.05, -2.16, -1.35, -2.25],
  [-2.57, -3.12, -3.98, -3.63, -4.03, -3.07, -3.41, -3.16, -1.83, -1.72, -1.90, -1.62, -1.64, -1.80, -2.29, -2.27, -2.16, -1.55, -0.59, -1.70],
  [-1.95, -2.48, -3.36, -3.01, -3.37, -2.49, -2.69, -2.60, -1.31, -1.15, -1.31, -1.05, -1.21, -1.29, -1.68, -1.80, -1.35, -0.59, -0.12, -0.97],
  [-3.07, -3.45, -4.25, -3.76, -4.20, -3.32, -3.73, -3.19, -2.03, -1.87, -1.90, -1.57, -1.53, -1.73, -1.33, -1.26, -2.25, -1.70, -0.97, -1.75]
];

// Hydrophobic set used by the HP mode (Dill style). Everything else is P.
var HYDROPHOBIC = 'ACFILMVWY';

function isHydrophobic(aa1) { return HYDROPHOBIC.indexOf(aa1) >= 0; }

var HP_MATRIX = (function () {
  var m = [];
  for (var i = 0; i < 20; i++) {
    var row = [];
    for (var j = 0; j < 20; j++) {
      row.push(isHydrophobic(AA20[i]) && isHydrophobic(AA20[j]) ? -1.0 : 0.0);
    }
    m.push(row);
  }
  return m;
})();

function mjPairEnergy(aa1, aa2) {
  var i = MJ_INDEX[aa1];
  var j = MJ_INDEX[aa2];
  if (i === undefined || j === undefined) {
    throw new Error('unknown residue letter: ' + aa1 + ' or ' + aa2);
  }
  return MJ1996[i][j];
}

// ---------------------------------------------------------------------------
// Bond geometry helpers (angles and dihedrals in radians internally)
// ---------------------------------------------------------------------------
function bondAngle(p1, p2, p3) {
  // Angle at p2 between p1-p2 and p3-p2.
  var a = vsub(p1, p2);
  var b = vsub(p3, p2);
  var denom = vlen(a) * vlen(b);
  if (denom === 0) { throw new Error('bondAngle with coincident points'); }
  var c = vdot(a, b) / denom;
  if (c > 1) { c = 1; } else if (c < -1) { c = -1; }
  return Math.acos(c);
}

function dihedral(p1, p2, p3, p4) {
  // IUPAC torsion angle of the p1-p2-p3-p4 series, in radians in (-pi, pi].
  var b1 = vsub(p2, p1);
  var b2 = vsub(p3, p2);
  var b3 = vsub(p4, p3);
  var n1 = vcross(b1, b2);
  var n2 = vcross(b2, b3);
  var b2len = vlen(b2);
  if (b2len === 0 || vlen(n1) === 0 || vlen(n2) === 0) {
    throw new Error('dihedral undefined for colinear or coincident points');
  }
  var m = vcross(n1, vscale(b2, 1 / b2len));
  var x = vdot(n1, n2);
  var y = vdot(m, n2);
  return Math.atan2(y, x);
}

/**
 * NeRF (Natural Extension Reference Frame) placement.
 * Given placed atoms a, b, c, returns the position d such that:
 *   |d - c|                = bondLength
 *   angle(b, c, d)         = angle      (radians)
 *   dihedral(a, b, c, d)   = torsion    (radians)
 */
function placeAtom(a, b, c, bondLength, angle, torsion) {
  var b2 = vsub(c, b);
  var b2len = vlen(b2);
  if (b2len === 0) { throw new Error('placeAtom: b and c coincide'); }
  var bhat = vscale(b2, 1 / b2len);
  var b1 = vsub(b, a);
  var n1 = vcross(b1, b2);
  var n1len = vlen(n1);
  if (n1len === 0) { throw new Error('placeAtom: a, b, c are colinear'); }
  var n1hat = vscale(n1, 1 / n1len);
  var mhat = vcross(n1hat, bhat); // in-plane perpendicular, toward the a side
  var sinA = Math.sin(angle);
  var local = [
    -bondLength * Math.cos(angle),
    bondLength * sinA * Math.cos(torsion),
    -bondLength * sinA * Math.sin(torsion)
  ];
  // d = c + local[0]*bhat + local[1]*mhat + local[2]*n1hat
  return [
    c[0] + local[0] * bhat[0] + local[1] * mhat[0] + local[2] * n1hat[0],
    c[1] + local[0] * bhat[1] + local[1] * mhat[1] + local[2] * n1hat[1],
    c[2] + local[0] * bhat[2] + local[1] * mhat[2] + local[2] * n1hat[2]
  ];
}

// ---------------------------------------------------------------------------
// Ramachandran classification
// ---------------------------------------------------------------------------
var DEG = 180 / Math.PI;
var RAD = Math.PI / 180;

function normalizeDeg(d) {
  var x = ((d + 180) % 360 + 360) % 360 - 180;
  if (x === -180) { x = 180; }
  return x;
}

/**
 * Ramachandran region definitions. Order matters: the classifier returns the
 * first region whose phi and psi windows contain the point. Bounds are in
 * degrees, inclusive. Colors are hex strings used by the web viewer.
 */
var RAMACHANDRAN_REGIONS = [
  { name: 'ALPHA_R', label: 'alpha right',   phi: [-165, -20], psi: [-125, 55],  color: '#e0564e' },
  { name: 'ALPHA_L', label: 'alpha left',    phi: [35, 100],   psi: [-45, 95],   color: '#e0913e' },
  { name: 'PP_II',   label: 'polyproline II', phi: [-105, -35], psi: [105, 165], color: '#3ec6a8' },
  { name: 'EPSILON', label: 'beta bridge',   phi: [-180, -60], psi: [-155, -55], color: '#9a6ee0' },
  { name: 'BETA',    label: 'beta sheet',    phi: [-180, -35], psi: [-180, -140], color: '#4e7de0', psiAlt: [55, 180] },
  { name: 'GAMMA',   label: 'gamma turn',    phi: [40, 100],   psi: [-95, -15], color: '#c6b93a' }
];

function inRange(v, range) { return v >= range[0] && v <= range[1]; }

/**
 * Classify a (phi, psi) pair in degrees into a Ramachandran region name.
 * Returns one of ALPHA_R, ALPHA_L, PP_II, EPSILON, BETA, GAMMA, OUTLIER.
 */
function classifyRamachandran(phi, psi) {
  var p = normalizeDeg(phi);
  var s = normalizeDeg(psi);
  for (var i = 0; i < RAMACHANDRAN_REGIONS.length; i++) {
    var r = RAMACHANDRAN_REGIONS[i];
    if (!inRange(p, r.phi)) { continue; }
    if (inRange(s, r.psi)) { return r.name; }
    if (r.psiAlt && inRange(s, r.psiAlt)) { return r.name; }
  }
  return 'OUTLIER';
}

/** Peptide bond omega: cis below 90 degrees of |omega|, trans otherwise. */
function classifyOmega(omegaDeg) {
  var w = Math.abs(normalizeDeg(omegaDeg));
  return w < 90 ? 'CIS' : 'TRANS';
}

// ---------------------------------------------------------------------------
// Steric clash detection
// ---------------------------------------------------------------------------
var VDW_RADII = { H: 1.10, C: 1.70, N: 1.55, O: 1.52, S: 1.80, P: 1.80 };
var DEFAULT_VDW_RADIUS = 1.70;

var TWO_LETTER_ELEMENTS = ['CL', 'NA', 'MG', 'FE', 'ZN', 'SE', 'MN', 'NI', 'CU', 'BR', 'CA'];

function inferElement(atomName, declared) {
  if (declared && declared.length >= 1 && /^[A-Za-z]{1,2}$/.test(declared.trim())) {
    var up = declared.trim().toUpperCase();
    if (VDW_RADII[up] || TWO_LETTER_ELEMENTS.indexOf(up) >= 0 || up.length === 1) { return up; }
  }
  // Work on the raw PDB atom name field, which is right justified for
  // single letter elements: the alpha carbon field is " CA " while calcium
  // is written left justified as "CA  ". The element symbol therefore
  // starts at the first alphabetic character, and a two letter element is
  // only possible when that character sits in the first column.
  var field = String(atomName || '');
  var k = -1;
  for (var i = 0; i < field.length; i++) {
    if (/[A-Za-z]/.test(field[i])) { k = i; break; }
  }
  if (k < 0) { return 'C'; }
  if (k === 0 && field.length >= 2) {
    var two = field.slice(0, 2).toUpperCase();
    if (TWO_LETTER_ELEMENTS.indexOf(two) >= 0) { return two; }
  }
  return field[k].toUpperCase();
}

function atomRadius(atom) {
  if (atom && typeof atom.radius === 'number') { return atom.radius; }
  var el = atom ? (atom.element || inferElement(atom.name)) : 'C';
  return VDW_RADII[el] || DEFAULT_VDW_RADIUS;
}

/**
 * Find steric clashes between atoms. Each atom needs xyz (and optionally
 * element, name or radius). A clash is a pair closer than
 * factor * (radius_i + radius_j). Pairs closer than minBond are treated as
 * covalently bonded and skipped, and pairs for which excludePair(i, j)
 * returns true are skipped as well. Returns an array of
 * { i, j, distance, overlap } sorted by worst overlap first.
 */
function findStericClashes(atoms, options) {
  var opts = options || {};
  var factor = typeof opts.factor === 'number' ? opts.factor : 0.75;
  var minBond = typeof opts.minBond === 'number' ? opts.minBond : 2.1;
  var exclude = typeof opts.excludePair === 'function' ? opts.excludePair : null;
  var radii = atoms.map(atomRadius);
  var clashes = [];
  for (var i = 0; i < atoms.length; i++) {
    for (var j = i + 1; j < atoms.length; j++) {
      var d = vdist(atoms[i].xyz, atoms[j].xyz);
      if (d < minBond) { continue; }
      if (exclude && exclude(i, j)) { continue; }
      var limit = factor * (radii[i] + radii[j]);
      if (d < limit) {
        clashes.push({ i: i, j: j, distance: d, overlap: limit - d });
      }
    }
  }
  clashes.sort(function (a, b) { return b.overlap - a.overlap; });
  return clashes;
}

// ---------------------------------------------------------------------------
// Radius of gyration
// ---------------------------------------------------------------------------
function radiusOfGyration(points) {
  var n = points.length;
  if (n === 0) { throw new Error('radiusOfGyration of empty point set'); }
  var cm = [0, 0, 0];
  for (var i = 0; i < n; i++) { cm = vadd(cm, points[i]); }
  cm = vscale(cm, 1 / n);
  var sum = 0;
  for (var j = 0; j < n; j++) {
    var d = vsub(points[j], cm);
    sum += vdot(d, d);
  }
  return Math.sqrt(sum / n);
}

// ---------------------------------------------------------------------------
// PDB parsing and export
// ---------------------------------------------------------------------------
function parsePdbLine(line) {
  var record = line.slice(0, 6).trim();
  if (record !== 'ATOM' && record !== 'HETATM') { return null; }
  function num(start, end) {
    var s = line.slice(start, end).trim();
    return s === '' ? null : parseFloat(s);
  }
  var serial = num(6, 11);
  var name = line.slice(12, 16);
  var altLoc = line.length > 16 ? line[16] : ' ';
  var resName = line.slice(17, 20).trim();
  var chainID = line.length > 21 ? line[21] : ' ';
  var resSeq = num(22, 26);
  var iCode = line.length > 26 ? line[26] : ' ';
  var x = num(30, 38);
  var y = num(38, 46);
  var z = num(46, 54);
  if (x === null || y === null || z === null) { return null; }
  var occupancy = num(54, 60);
  var tempFactor = num(60, 66);
  var element = line.length >= 78 ? line.slice(76, 78).trim() : '';
  return {
    record: record,
    serial: serial,
    name: name,
    altLoc: altLoc,
    resName: resName,
    chainID: chainID,
    resSeq: resSeq,
    iCode: iCode,
    xyz: [x, y, z],
    occupancy: occupancy,
    tempFactor: tempFactor,
    element: element ? element.toUpperCase() : ''
  };
}

/**
 * Parse the ATOM records of a PDB file. By default only the first MODEL is
 * read and alternate location indicators other than blank and A are skipped.
 * HETATM records are included only with includeHetatm: true.
 */
function parsePdb(text, options) {
  var opts = options || {};
  var lines = String(text).split(/\r\n|\r|\n/);
  var atoms = [];
  var modelIndex = 0;
  var firstModel = -1;
  var remarks = [];
  for (var i = 0; i < lines.length; i++) {
    var line = lines[i];
    if (line.slice(0, 6) === 'ENDMDL') { break; }
    if (line.slice(0, 6) === 'MODEL ') {
      modelIndex++;
      if (firstModel < 0) { firstModel = modelIndex; }
      if (modelIndex > firstModel) { break; }
      continue;
    }
    if (line.slice(0, 6) === 'REMARK') { remarks.push(line.slice(6).trim()); continue; }
    var parsed = parsePdbLine(line);
    if (!parsed) { continue; }
    if (parsed.record === 'HETATM' && !opts.includeHetatm) { continue; }
    if (parsed.altLoc !== ' ' && parsed.altLoc !== 'A' && parsed.altLoc !== '1') { continue; }
    if (!parsed.element) { parsed.element = inferElement(parsed.name); }
    atoms.push(parsed);
  }
  if (atoms.length === 0) {
    throw new Error('no ATOM records found in PDB input');
  }
  return { atoms: atoms, remarks: remarks, numModelsRead: modelIndex };
}

function padLeft(s, n) { s = String(s); while (s.length < n) { s = ' ' + s; } return s; }
function padRight(s, n) { s = String(s); while (s.length < n) { s = s + ' '; } return s; }
function fmtCoord(v) {
  var s = v.toFixed(3);
  if (s.length > 8) { s = v.toFixed(2); }
  return padLeft(s, 8);
}

function pdbAtomNameField(name) {
  // Standard PDB column rules: single and two character atom names start in
  // column 14 (index 13) so the element symbol lines up on the right.
  var n = String(name).trim();
  if (n.length >= 4) { return n.slice(0, 4); }
  if (n.length === 3) { return ' ' + n; }
  if (n.length === 2) { return ' ' + n + ' '; }
  return ' ' + n + '  ';
}

/**
 * Export atoms to PDB ATOM record text. Atoms need name, resName, resSeq,
 * chainID (optional) and xyz. Element is inferred when absent.
 */
function atomsToPdb(atoms, options) {
  var opts = options || {};
  var out = [];
  if (opts.title) {
    out.push('HEADER    ' + String(opts.title).slice(0, 50));
  }
  out.push('REMARK   2 fold-lattice ' + VERSION + ' pdb export');
  if (opts.remarks) {
    for (var r = 0; r < opts.remarks.length; r++) {
      out.push('REMARK   3 ' + String(opts.remarks[r]).slice(0, 66));
    }
  }
  var serial = opts.startSerial || 1;
  var prevChain = null;
  var prevResName = null;
  for (var i = 0; i < atoms.length; i++) {
    var a = atoms[i];
    if (prevChain !== null && a.chainID !== prevChain) {
      out.push('TER   ' + padLeft(serial, 5) + '      ' + padRight(prevResName || '', 3) + ' ' + (prevChain || ' '));
      serial++;
    }
    prevChain = a.chainID || ' ';
    prevResName = a.resName || '';
    var element = (a.element || inferElement(a.name)).toUpperCase();
    var line =
      'ATOM  ' + padLeft(serial, 5) + ' ' + pdbAtomNameField(a.name) +
      (a.altLoc || ' ') + padRight(a.resName || 'UNK', 3) + ' ' + (a.chainID || ' ') +
      padLeft(a.resSeq, 4) + (a.iCode || ' ') + '   ' +
      fmtCoord(a.xyz[0]) + fmtCoord(a.xyz[1]) + fmtCoord(a.xyz[2]) +
      padLeft(typeof a.occupancy === 'number' ? a.occupancy.toFixed(2) : '1.00', 6) +
      padLeft(typeof a.tempFactor === 'number' ? a.tempFactor.toFixed(2) : '0.00', 6) +
      '          ' + padLeft(element, 2);
    out.push(line);
    serial++;
  }
  out.push('TER   ' + padLeft(serial, 5) + '      ' + padRight(prevResName || '', 3) + ' ' + (prevChain || ' '));
  out.push('END');
  return out.join('\n') + '\n';
}

// ---------------------------------------------------------------------------
// All-atom backbone reconstruction from (phi, psi, omega) angles
// ---------------------------------------------------------------------------
var BACKBONE_GEOM = {
  L_N_CA: 1.458,   // N - CA bond length in angstrom
  L_CA_C: 1.525,   // CA - C bond length in angstrom
  L_C_N: 1.329,    // C - N peptide bond length in angstrom
  L_C_O: 1.231,    // C = O bond length in angstrom
  A_N_CA_C: 111.2, // degrees, angle at CA
  A_CA_C_N: 116.2, // degrees, angle at C
  A_C_N_CA: 121.7, // degrees, angle at N
  A_CA_C_O: 120.8  // degrees, angle at C for the carbonyl oxygen
};

/**
 * Build an all-atom protein backbone (N, CA, C, O per residue) from a
 * sequence and per-residue dihedral angles given in degrees.
 * angles[i] = { phi, psi, omega } where omega defaults to 180 (trans).
 * The first residue uses its phi for nothing (no predecessor) but its psi
 * places the carbonyl oxygen. Returns an array of atom objects:
 * { name, resIndex, resName, element, xyz }.
 */
function reconstructBackbone(sequence, angles) {
  var seq = validateSequence(sequence);
  var n = seq.length;
  if (!Array.isArray(angles) || angles.length !== n) {
    throw new Error('angles array length must equal sequence length');
  }
  for (var v = 0; v < n; v++) {
    var entry = angles[v];
    if (!entry || typeof entry.phi !== 'number' || typeof entry.psi !== 'number') {
      throw new Error('angles[' + v + '] must provide numeric phi and psi');
    }
  }
  var g = BACKBONE_GEOM;
  var atoms = [];
  function pushAtom(name, resIndex, element, xyz) {
    atoms.push({
      name: name,
      resIndex: resIndex,
      resName: AA1_TO_3[seq[resIndex]],
      element: element,
      chainID: 'A',
      resSeq: resIndex + 1,
      xyz: xyz
    });
  }
  // Residue 0 seed geometry: N at origin, CA along +x, C in the xy plane.
  var nPrev = [0, 0, 0];
  var caPrev = [g.L_N_CA, 0, 0];
  var caCdir = [
    -Math.cos(g.A_N_CA_C * RAD),
    Math.sin(g.A_N_CA_C * RAD),
    0
  ];
  var cPrev = vadd(caPrev, vscale(caCdir, g.L_CA_C));
  pushAtom('N', 0, 'N', nPrev);
  pushAtom('CA', 0, 'C', caPrev);
  pushAtom('C', 0, 'C', cPrev);
  var psi0 = normalizeDeg(angles[0].psi);
  var o0 = placeAtom(nPrev, caPrev, cPrev, g.L_C_O, g.A_CA_C_O * RAD, (psi0 + 180) * RAD);
  pushAtom('O', 0, 'O', o0);

  for (var i = 1; i < n; i++) {
    var a = angles[i];
    var psiPrev = normalizeDeg(angles[i - 1].psi);
    var omega = a.omega === undefined ? 180 : normalizeDeg(a.omega);
    var phi = normalizeDeg(a.phi);
    var psi = normalizeDeg(a.psi);
    // N_i: dihedral psi of residue i-1 around CA(i-1)-C(i-1)
    var nI = placeAtom(nPrev, caPrev, cPrev, g.L_C_N, g.A_CA_C_N * RAD, psiPrev * RAD);
    // CA_i: dihedral omega around C(i-1)-N_i
    var caI = placeAtom(caPrev, cPrev, nI, g.L_N_CA, g.A_C_N_CA * RAD, omega * RAD);
    // C_i: dihedral phi of residue i around N_i-CA_i
    var cI = placeAtom(cPrev, nI, caI, g.L_CA_C, g.A_N_CA_C * RAD, phi * RAD);
    // O_i: dihedral psi + 180 around CA_i-C_i
    var oI = placeAtom(nI, caI, cI, g.L_C_O, g.A_CA_C_O * RAD, (psi + 180) * RAD);
    pushAtom('N', i, 'N', nI);
    pushAtom('CA', i, 'C', caI);
    pushAtom('C', i, 'C', cI);
    pushAtom('O', i, 'O', oI);
    nPrev = nI; caPrev = caI; cPrev = cI;
  }
  return atoms;
}

/**
 * Compute per-residue backbone dihedrals from parsed PDB atoms.
 * Returns one entry per residue: { resIndex, resName, resSeq, chainID,
 * phi, psi, omega, region } with angles in degrees, or null where the
 * neighboring residue needed for a torsion is missing.
 */
function computeDihedrals(atoms) {
  var groups = [];
  var byKey = {};
  for (var i = 0; i < atoms.length; i++) {
    var a = atoms[i];
    var iCode = a.iCode && a.iCode !== ' ' ? a.iCode : '';
    var key = a.chainID + ':' + a.resSeq + ':' + iCode;
    if (!byKey[key]) {
      byKey[key] = {
        chainID: a.chainID, resSeq: a.resSeq, resName: a.resName,
        N: null, CA: null, C: null
      };
      groups.push(byKey[key]);
    }
    var g = byKey[key];
    var nm = a.name.trim();
    if (nm === 'N' && !g.N) { g.N = a.xyz; }
    if (nm === 'CA' && !g.CA) { g.CA = a.xyz; }
    if (nm === 'C' && !g.C) { g.C = a.xyz; }
  }
  // Consecutive groups count as bonded only when they share the chain and
  // residue numbering is continuous; otherwise there is a chain break.
  function bonded(g1, g2) {
    return g1 && g2 && g1.chainID === g2.chainID &&
      (g2.resSeq === g1.resSeq + 1);
  }
  var out = [];
  for (var r = 0; r < groups.length; r++) {
    var cur = groups[r];
    var prev = r > 0 ? groups[r - 1] : null;
    var next = r < groups.length - 1 ? groups[r + 1] : null;
    var entry = {
      resIndex: r,
      resSeq: cur.resSeq,
      resName: cur.resName,
      chainID: cur.chainID,
      phi: null, psi: null, omega: null, region: null, omegaClass: null
    };
    var ok = cur.N && cur.CA && cur.C;
    if (ok && bonded(prev, cur) && prev.C) {
      try { entry.phi = dihedral(prev.C, cur.N, cur.CA, cur.C) * DEG; }
      catch (err) { entry.phi = null; }
    }
    if (ok && bonded(cur, next) && next.N) {
      try { entry.psi = dihedral(cur.N, cur.CA, cur.C, next.N) * DEG; }
      catch (err) { entry.psi = null; }
    }
    if (ok && bonded(prev, cur) && prev.CA && prev.C) {
      try { entry.omega = dihedral(prev.CA, prev.C, cur.N, cur.CA) * DEG; }
      catch (err) { entry.omega = null; }
    }
    if (entry.phi !== null && entry.psi !== null) {
      entry.region = classifyRamachandran(entry.phi, entry.psi);
    }
    if (entry.omega !== null) { entry.omegaClass = classifyOmega(entry.omega); }
    out.push(entry);
  }
  return out;
}

/**
 * Map a lattice CA trace (one unit per residue) to per-residue backbone
 * dihedral angles. Heuristic used to feed lattice conformations into
 * all-atom rebuilders: straight CA segments map to the beta basin, and
 * right-angle turns map to the alpha right or alpha left basins according
 * to the sign of the CA trace pseudodihedral over (i-1, i, i+1, i+2), whose
 * positive sign corresponds to the right-handed basin in the standard CA
 * pseudodihedral convention. Near-zero pseudodihedrals (local reversals)
 * map to the polyproline II basin. Chain ends without a full angle default
 * to beta.
 */
var LATTICE_ANGLE_MAP = {
  BETA: { phi: -139, psi: 135 },
  ALPHA_R: { phi: -57, psi: -47 },
  ALPHA_L: { phi: 57, psi: 47 },
  PP_II: { phi: -75, psi: 145 }
};

function latticeBackboneAngles(coords) {
  var n = coords.length;
  var out = [];
  for (var i = 0; i < n; i++) {
    var cls;
    var hasPrev = i > 0;
    var hasNext = i < n - 1;
    if (!hasPrev || !hasNext) {
      cls = 'BETA';
    } else {
      var ang = bondAngle(coords[i - 1], coords[i], coords[i + 1]) * DEG;
      if (ang > 155) {
        cls = 'BETA';
      } else if (i <= n - 3) {
        var pd = null;
        try {
          pd = dihedral(coords[i - 1], coords[i], coords[i + 1], coords[i + 2]) * DEG;
        } catch (err) {
          // Pseudodihedral undefined when the next bond runs colinear with
          // the turn axis; fall through to the turn-only assignment below.
          pd = null;
        }
        if (pd === null) { cls = ang > 90 ? 'ALPHA_R' : 'PP_II'; }
        else if (pd > 25) { cls = 'ALPHA_R'; }
        else if (pd < -25) { cls = 'ALPHA_L'; }
        else { cls = 'PP_II'; }
      } else {
        cls = 'PP_II';
      }
    }
    var a = LATTICE_ANGLE_MAP[cls];
    out.push({ phi: a.phi, psi: a.psi, omega: 180, region: cls });
  }
  return out;
}

function validateSequence(sequence) {
  if (typeof sequence !== 'string') { throw new Error('sequence must be a string'); }
  var seq = sequence.toUpperCase().replace(/[^A-Z]/g, '');
  if (seq.length !== sequence.replace(/\s/g, '').length) {
    throw new Error('sequence contains invalid characters: only the 20 standard amino acid letters are allowed');
  }
  for (var i = 0; i < seq.length; i++) {
    if (AA20.indexOf(seq[i]) < 0) {
      throw new Error('sequence contains invalid residue "' + seq[i] + '" at position ' + i);
    }
  }
  if (seq.length < 2) { throw new Error('sequence must contain at least 2 residues'); }
  return seq;
}

// ---------------------------------------------------------------------------
// 3D cubic lattice folder
// ---------------------------------------------------------------------------
var KEY_OFFSET = 2048;
var KEY_K = 4096;

function cellKey(x, y, z) {
  return (x + KEY_OFFSET) + KEY_K * (y + KEY_OFFSET) + KEY_K * KEY_K * (z + KEY_OFFSET);
}

/**
 * A self-avoiding walk on the 3D cubic lattice scored by a contact matrix.
 * Coordinates are integer lattice cells; one lattice step is one residue.
 */
class LatticeFolder {
  constructor(sequence, options) {
    var opts = options || {};
    var seq = validateSequence(sequence);
    this.sequence = seq;
    this.n = seq.length;
    this.mode = opts.mode === 'hp' ? 'hp' : 'mj';
    this.matrix = this.mode === 'hp' ? HP_MATRIX : MJ1996;
    var letters = [];
    var mjIdx = [];
    for (var i = 0; i < this.n; i++) {
      letters.push(seq[i]);
      mjIdx.push(this.mode === 'hp'
        ? (isHydrophobic(seq[i]) ? 0 : 1)
        : MJ_INDEX[seq[i]]);
    }
    this._letters = letters;
    this._mjIdx = mjIdx;
    if (this.mode === 'hp') {
      this._pair = function (i, j) {
        return isHydrophobic(letters[i]) && isHydrophobic(letters[j]) ? -1 : 0;
      };
    } else {
      this._pair = function (i, j) {
        return MJ1996[mjIdx[i]][mjIdx[j]];
      };
    }
    // Straight line start conformation along +x.
    this.coords = [];
    this.occ = new Map();
    for (var k = 0; k < this.n; k++) {
      this.coords.push([k, 0, 0]);
      this.occ.set(cellKey(k, 0, 0), k);
    }
    this._energy = null;
    this._contacts = null;
  }

  static key(x, y, z) { return cellKey(x, y, z); }

  indexAt(x, y, z) { return this.occ.get(cellKey(x, y, z)); }

  /** Full validation: self avoidance plus chain connectivity. */
  isValidConformation(coords) {
    var c = coords || this.coords;
    var n = this.n;
    if (c.length !== n) { return false; }
    var seen = new Set();
    for (var i = 0; i < n; i++) {
      var k = cellKey(c[i][0], c[i][1], c[i][2]);
      if (seen.has(k)) { return false; }
      seen.add(k);
    }
    for (var j = 1; j < n; j++) {
      var dx = Math.abs(c[j][0] - c[j - 1][0]);
      var dy = Math.abs(c[j][1] - c[j - 1][1]);
      var dz = Math.abs(c[j][2] - c[j - 1][2]);
      if (dx + dy + dz !== 1) { return false; }
    }
    return true;
  }

  contacts() {
    if (this._contacts) { return this._contacts; }
    var list = [];
    var coords = this.coords;
    for (var i = 0; i < this.n; i++) {
      var p = coords[i];
      for (var d = 0; d < 6; d++) {
        var ax = AXIS_VECTORS[d];
        var j = this.occ.get(cellKey(p[0] + ax[0], p[1] + ax[1], p[2] + ax[2]));
        if (j !== undefined && j > i + 1) { list.push([i, j]); }
      }
    }
    list.sort(function (a, b) { return a[0] - b[0] || a[1] - b[1]; });
    this._contacts = list;
    return list;
  }

  contactCount() { return this.contacts().length; }

  energy() {
    if (this._energy !== null) { return this._energy; }
    var list = this.contacts();
    var e = 0;
    var pair = this._pair;
    for (var i = 0; i < list.length; i++) { e += pair(list[i][0], list[i][1]); }
    this._energy = e;
    return e;
  }

  radiusOfGyration() { return radiusOfGyration(this.coords); }

  copy() {
    var f = Object.create(LatticeFolder.prototype);
    f.sequence = this.sequence;
    f.n = this.n;
    f.mode = this.mode;
    f.matrix = this.matrix;
    f._letters = this._letters.slice();
    f._mjIdx = this._mjIdx.slice();
    f._pair = this._pair;
    f.coords = this.coords.map(function (c) { return c.slice(); });
    f.occ = new Map(this.occ);
    f._energy = this._energy;
    f._contacts = this._contacts ? this._contacts.map(function (c) { return c.slice(); }) : null;
    return f;
  }

  swapStateWith(other) {
    var t;
    t = this.coords; this.coords = other.coords; other.coords = t;
    t = this.occ; this.occ = other.occ; other.occ = t;
    t = this._energy; this._energy = other._energy; other._energy = t;
    t = this._contacts; this._contacts = other._contacts; other._contacts = t;
  }

  _commit(coords, occ, energy, contacts) {
    this.coords = coords;
    this.occ = occ;
    this._energy = energy;
    this._contacts = contacts;
  }

  /**
   * Build occupancy for a candidate conformation and score it without
   * committing. Returns { occ, energy, contacts } or null when the candidate
   * is invalid (overlap or broken connectivity).
   */
  _evaluate(cand) {
    var occ = new Map();
    for (var i = 0; i < this.n; i++) {
      var k = cellKey(cand[i][0], cand[i][1], cand[i][2]);
      if (occ.has(k)) { return null; }
      occ.set(k, i);
    }
    for (var j = 1; j < this.n; j++) {
      var d = Math.abs(cand[j][0] - cand[j - 1][0]) +
        Math.abs(cand[j][1] - cand[j - 1][1]) +
        Math.abs(cand[j][2] - cand[j - 1][2]);
      if (d !== 1) { return null; }
    }
    var e = 0;
    var contacts = [];
    var pair = this._pair;
    for (var a = 0; a < this.n; a++) {
      var p = cand[a];
      for (var d2 = 0; d2 < 6; d2++) {
        var ax = AXIS_VECTORS[d2];
        var b = occ.get(cellKey(p[0] + ax[0], p[1] + ax[1], p[2] + ax[2]));
        if (b !== undefined && b > a + 1) {
          contacts.push([a, b]);
          e += pair(a, b);
        }
      }
    }
    return { occ: occ, energy: e, contacts: contacts };
  }

  _tryCommit(cand) {
    var built = this._evaluate(cand);
    if (!built) { return null; }
    var dE = built.energy - this.energy();
    this._commit(cand, built.occ, built.energy, built.contacts);
    return dE;
  }

  /**
   * Pivot move: rotate the tail beyond bond i about the axis through
   * residue i directed toward residue i+1, by 1 to 3 quarter turns.
   */
  _pivotCandidate(rng) {
    var n = this.n;
    var i = Math.floor(rng() * (n - 1)); // bond between i and i+1
    var u = vsub(this.coords[i + 1], this.coords[i]);
    if (!isUnitAxis(u)) { return null; }
    var turns = 1 + Math.floor(rng() * 3);
    var R = rotationAboutAxis(u, turns);
    var pivot = this.coords[i];
    var cand = new Array(n);
    for (var k = 0; k < n; k++) {
      if (k <= i) { cand[k] = this.coords[k].slice(); continue; }
      var v = vsub(this.coords[k], pivot);
      var w = [
        R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]
      ];
      cand[k] = [pivot[0] + w[0], pivot[1] + w[1], pivot[2] + w[2]];
    }
    return cand;
  }

  /**
   * Corner flip: move an interior residue to the fourth corner of the unit
   * square spanned by its two chain neighbors (a one-bead crankshaft).
   */
  _cornerFlipCandidate(rng) {
    var n = this.n;
    if (n < 3) { return null; }
    var i = 1 + Math.floor(rng() * (n - 2));
    var p = this.coords[i];
    var c = [
      this.coords[i - 1][0] + this.coords[i + 1][0] - p[0],
      this.coords[i - 1][1] + this.coords[i + 1][1] - p[1],
      this.coords[i - 1][2] + this.coords[i + 1][2] - p[2]
    ];
    var d = Math.abs(c[0] - p[0]) + Math.abs(c[1] - p[1]) + Math.abs(c[2] - p[2]);
    if (d !== 2) { return null; } // not a right-angle corner
    var cand = this.coords.map(function (q) { return q.slice(); });
    cand[i] = c;
    return cand;
  }

  /**
   * Two-bead crankshaft: when residues j..j+3 form a square loop (residues
   * j and j+3 are lattice adjacent), rotate the middle two residues by 180
   * degrees about the axis through the two ends.
   */
  _crankshaftCandidate(rng) {
    var n = this.n;
    if (n < 4) { return null; }
    var j = Math.floor(rng() * (n - 3));
    var p0 = this.coords[j];
    var d1 = vsub(this.coords[j + 1], p0);
    var d2 = vsub(this.coords[j + 2], this.coords[j + 1]);
    var d3 = vsub(this.coords[j + 3], this.coords[j + 2]);
    if (!isUnitAxis(d1) || !isUnitAxis(d2) || !isUnitAxis(d3)) { return null; }
    if (!veq(d3, vscale(d1, -1))) { return null; } // not a square loop
    var a = [p0[0] - d1[0], p0[1] - d1[1], p0[2] - d1[2]];
    var b = [a[0] + d2[0], a[1] + d2[1], a[2] + d2[2]];
    var cand = this.coords.map(function (q) { return q.slice(); });
    cand[j + 1] = a;
    cand[j + 2] = b;
    return cand;
  }

  /**
   * Pull move (Lesh, Mitzenmacher and Whitesides, adapted to the cubic
   * lattice). Residue i is pulled into a free cell L that is adjacent to its
   * chain neighbor and diagonally adjacent to the residue itself; the corner
   * cell C of that unit square must be free or already occupied by the next
   * residue on the far side, and the chain is then pulled two cells at a
   * time until connectivity is restored.
   */
  _pullCandidate(rng, dir) {
    var n = this.n;
    var i = Math.floor(rng() * n);
    // dir = -1 pulls the lower side (neighbor j = i+1), dir = +1 pulls the
    // upper side (neighbor j = i-1).
    if (dir < 0 && i >= n - 1) { return null; }
    if (dir > 0 && i < 1) { return null; }
    var j = dir < 0 ? i + 1 : i - 1;
    var pi = this.coords[i];
    var pj = this.coords[j];
    var A = vsub(pj, pi); // unit axis from i toward its neighbor
    if (!isUnitAxis(A)) { return null; }
    // The two axes perpendicular to A.
    var perp = [];
    for (var ax = 0; ax < 3; ax++) {
      if (A[ax] === 0) { perp.push(ax); }
    }
    var perpIndex = Math.floor(rng() * 2);
    var sign = rng() < 0.5 ? 1 : -1;
    var B = [0, 0, 0];
    B[perp[perpIndex]] = sign;
    var L = [pj[0] + B[0], pj[1] + B[1], pj[2] + B[2]];
    var C = [pi[0] + B[0], pi[1] + B[1], pi[2] + B[2]];
    var keyL = cellKey(L[0], L[1], L[2]);
    if (this.occ.has(keyL)) { return null; }
    var cand = this.coords.map(function (q) { return q.slice(); });
    if (dir < 0) {
      var far = i - 1;
      if (i >= 1) {
        var keyC = cellKey(C[0], C[1], C[2]);
        var occupant = this.occ.get(keyC);
        if (occupant !== far) {
          if (occupant !== undefined) { return null; }
          cand[i - 1] = C;
        }
      }
      cand[i] = L;
      // Pull two-by-two down the chain until connectivity is restored.
      for (var k = i - 2; k >= 0; k--) {
        var d = Math.abs(cand[k][0] - cand[k + 1][0]) +
          Math.abs(cand[k][1] - cand[k + 1][1]) +
          Math.abs(cand[k][2] - cand[k + 1][2]);
        if (d === 1) { break; }
        cand[k] = this.coords[k + 2].slice();
      }
    } else {
      var far2 = i + 1;
      if (i <= n - 2) {
        var keyC2 = cellKey(C[0], C[1], C[2]);
        var occupant2 = this.occ.get(keyC2);
        if (occupant2 !== far2) {
          if (occupant2 !== undefined) { return null; }
          cand[i + 1] = C;
        }
      }
      cand[i] = L;
      for (var k2 = i + 2; k2 < n; k2++) {
        var d2 = Math.abs(cand[k2][0] - cand[k2 - 1][0]) +
          Math.abs(cand[k2][1] - cand[k2 - 1][1]) +
          Math.abs(cand[k2][2] - cand[k2 - 1][2]);
        if (d2 === 1) { break; }
        cand[k2] = this.coords[k2 - 2].slice();
      }
    }
    return cand;
  }

  /**
   * Endpoint pull: move the terminal residue and its neighbor into a path
   * of two free cells attached to the tip, then pull the rest of the chain.
   */
  _endpointPullCandidate(rng) {
    var n = this.n;
    if (n < 3) { return null; }
    var atEnd = rng() < 0.5;
    var e = atEnd ? n - 1 : 0;
    var nb = atEnd ? n - 2 : 1;
    var pe = this.coords[e];
    // Free cells adjacent to the tip.
    var freeAdj = [];
    for (var d = 0; d < 6; d++) {
      var ax = AXIS_VECTORS[d];
      var f1 = [pe[0] + ax[0], pe[1] + ax[1], pe[2] + ax[2]];
      if (!this.occ.has(cellKey(f1[0], f1[1], f1[2]))) { freeAdj.push(f1); }
    }
    if (freeAdj.length === 0) { return null; }
    var F1 = freeAdj[Math.floor(rng() * freeAdj.length)];
    var freeAdj2 = [];
    for (var d2 = 0; d2 < 6; d2++) {
      var ax2 = AXIS_VECTORS[d2];
      var f2 = [F1[0] + ax2[0], F1[1] + ax2[1], F1[2] + ax2[2]];
      var k2 = cellKey(f2[0], f2[1], f2[2]);
      if (this.occ.has(k2)) { continue; }
      if (f2[0] === pe[0] && f2[1] === pe[1] && f2[2] === pe[2]) { continue; }
      freeAdj2.push(f2);
    }
    if (freeAdj2.length === 0) { return null; }
    var F2 = freeAdj2[Math.floor(rng() * freeAdj2.length)];
    var cand = this.coords.map(function (q) { return q.slice(); });
    cand[e] = F2;
    cand[nb] = F1;
    if (atEnd) {
      for (var k = e - 2; k >= 0; k--) {
        var dd = Math.abs(cand[k][0] - cand[k + 1][0]) +
          Math.abs(cand[k][1] - cand[k + 1][1]) +
          Math.abs(cand[k][2] - cand[k + 1][2]);
        if (dd === 1) { break; }
        cand[k] = this.coords[k + 2].slice();
      }
    } else {
      for (var k3 = 2; k3 < n; k3++) {
        var dd3 = Math.abs(cand[k3][0] - cand[k3 - 1][0]) +
          Math.abs(cand[k3][1] - cand[k3 - 1][1]) +
          Math.abs(cand[k3][2] - cand[k3 - 1][2]);
        if (dd3 === 1) { break; }
        cand[k3] = this.coords[k3 - 2].slice();
      }
    }
    return cand;
  }

  /**
   * Attempt one move of the requested type. With { trial: true } the
   * candidate is scored but not committed, and the result carries the
   * candidate data so the caller can decide whether to accept it.
   */
  tryMove(type, rng, options) {
    var trial = !!(options && options.trial);
    var cand = null;
    switch (type) {
      case 'pivot': cand = this._pivotCandidate(rng); break;
      case 'pull': cand = this._pullCandidate(rng, rng() < 0.5 ? -1 : 1); break;
      case 'pullDown': cand = this._pullCandidate(rng, -1); break;
      case 'pullUp': cand = this._pullCandidate(rng, 1); break;
      case 'crankshaft': cand = this._crankshaftCandidate(rng); break;
      case 'corner': cand = this._cornerFlipCandidate(rng); break;
      case 'endpull': cand = this._endpointPullCandidate(rng); break;
      default: throw new Error('unknown move type: ' + type);
    }
    if (!cand) { return { moved: false, accepted: false, dE: 0, type: type }; }
    var built = this._evaluate(cand);
    if (!built) { return { moved: false, accepted: false, dE: 0, type: type }; }
    var dE = built.energy - this.energy();
    if (!trial) { this._commit(cand, built.occ, built.energy, built.contacts); }
    return {
      moved: true, accepted: true, dE: dE, type: type,
      cand: cand, occ: built.occ, energy: built.energy, contacts: built.contacts
    };
  }

  /**
   * One Metropolis Monte Carlo step at the given temperature. Move types
   * are sampled with fixed weights; the candidate is committed only when
   * the Metropolis criterion accepts it, so rejected proposals leave the
   * conformation unchanged.
   */
  mcStep(temperature, rng) {
    var r = rng();
    var type;
    if (r < 0.40) { type = 'pull'; }
    else if (r < 0.60) { type = 'pivot'; }
    else if (r < 0.75) { type = 'crankshaft'; }
    else if (r < 0.90) { type = 'corner'; }
    else { type = 'endpull'; }
    var res = this.tryMove(type, rng, { trial: true });
    if (!res.moved) {
      return { moved: false, accepted: false, dE: 0, type: type, proposed: true };
    }
    var accept = res.dE <= 0 || temperature <= 0 || rng() < Math.exp(-res.dE / temperature);
    if (accept) { this._commit(res.cand, res.occ, res.energy, res.contacts); }
    return { moved: true, accepted: accept, dE: res.dE, type: type, proposed: true };
  }
}

// ---------------------------------------------------------------------------
// Lattice rotation matrices for pivot moves
// ---------------------------------------------------------------------------
var QUARTER_ABOUT_X = [
  [1, 0, 0],
  [0, 0, -1],
  [0, 1, 0]
];
var QUARTER_ABOUT_Y = [
  [0, 0, 1],
  [0, 1, 0],
  [-1, 0, 0]
];
var QUARTER_ABOUT_Z = [
  [0, -1, 0],
  [1, 0, 0],
  [0, 0, 1]
];

function matMul(A, B) {
  var C = [];
  for (var i = 0; i < 3; i++) {
    var row = [];
    for (var j = 0; j < 3; j++) {
      var s = 0;
      for (var k = 0; k < 3; k++) { s += A[i][k] * B[k][j]; }
      row.push(s);
    }
    C.push(row);
  }
  return C;
}

function matVec(A, v) {
  return [
    A[0][0] * v[0] + A[0][1] * v[1] + A[0][2] * v[2],
    A[1][0] * v[0] + A[1][1] * v[1] + A[1][2] * v[2],
    A[2][0] * v[0] + A[2][1] * v[1] + A[2][2] * v[2]
  ];
}

/**
 * Integer rotation matrix for `turns` quarter turns about the axis u
 * (u must be one of the six unit axis vectors). Rotation about a negative
 * axis is folded into the turn count so the matrix stays integral.
 */
function rotationAboutAxis(u, turns) {
  var base;
  var axisIndex = u[0] !== 0 ? 0 : (u[1] !== 0 ? 1 : 2);
  var sign = u[axisIndex] > 0 ? 1 : -1;
  var t = ((sign > 0 ? turns : -turns) % 4 + 4) % 4;
  if (axisIndex === 0) { base = QUARTER_ABOUT_X; }
  else if (axisIndex === 1) { base = QUARTER_ABOUT_Y; }
  else { base = QUARTER_ABOUT_Z; }
  if (t === 0) { return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]; }
  if (t === 1) { return base; }
  if (t === 2) { return matMul(base, base); }
  return matMul(matMul(base, base), base);
}

// ---------------------------------------------------------------------------
// Replica-exchange Monte Carlo (parallel tempering)
// ---------------------------------------------------------------------------

function geometricTemps(count, minTemp, maxTemp) {
  if (count < 1) { throw new Error('need at least one replica'); }
  if (count === 1) { return [minTemp]; }
  var out = [];
  var ratio = Math.pow(minTemp / maxTemp, 1 / (count - 1));
  for (var i = 0; i < count; i++) {
    out.push(maxTemp * Math.pow(ratio, i));
  }
  return out;
}

/**
 * Generator driving a replica-exchange simulation. Yields progress
 * snapshots; the final yield has done: true and carries the result. Pass
 * `yieldEvery` to control snapshot frequency.
 */
function* remcGenerator(options) {
  var opts = options || {};
  var seq = validateSequence(opts.sequence);
  var n = seq.length;
  var mode = opts.mode === 'hp' ? 'hp' : 'mj';
  var replicas = opts.replicas || 8;
  var sweeps = opts.sweeps || 1500;
  var stepsPerSweep = opts.stepsPerSweep || n;
  var minTemp = typeof opts.minTemp === 'number' ? opts.minTemp : 0.2;
  var maxTemp = typeof opts.maxTemp === 'number' ? opts.maxTemp : 2.5;
  var yieldEvery = opts.yieldEvery || 25;
  var seed = opts.seed === undefined ? 12345 : opts.seed;
  var rng = makeRng(seed);
  var temps = geometricTemps(replicas, minTemp, maxTemp);
  var folders = [];
  for (var r = 0; r < replicas; r++) {
    folders.push(new LatticeFolder(seq, { mode: mode }));
  }
  var exchanges = 0, exchangeTries = 0, accepted = 0, proposed = 0;
  var best = null;
  var trace = [];
  var started = Date.now();

  function considerBest(f) {
    var e = f.energy();
    if (!best || e < best.energy) {
      best = { energy: e, coords: f.coords.map(function (c) { return c.slice(); }), replica: 0 };
    }
  }
  for (var r0 = 0; r0 < replicas; r0++) { considerBest(folders[r0]); }

  for (var s = 0; s < sweeps; s++) {
    for (var r1 = 0; r1 < replicas; r1++) {
      var f = folders[r1];
      for (var st = 0; st < stepsPerSweep; st++) {
        var res = f.mcStep(temps[r1], rng);
        proposed++;
        if (res.accepted) { accepted++; }
      }
      considerBest(f);
    }
    // Replica exchanges between adjacent temperature slots, alternating
    // even and odd pairs each sweep.
    var offset = s % 2;
    for (var p = offset; p + 1 < replicas; p += 2) {
      exchangeTries++;
      var fa = folders[p];
      var fb = folders[p + 1];
      var ea = fa.energy();
      var eb = fb.energy();
      var betaA = 1 / temps[p];
      var betaB = 1 / temps[p + 1];
      var lnA = (betaA - betaB) * (ea - eb);
      var doSwap = lnA >= 0 || rng() < Math.exp(lnA);
      if (doSwap) {
        fa.swapStateWith(fb);
        exchanges++;
      }
    }
    var bestE = best.energy;
    trace.push(bestE);
    if (s % yieldEvery === 0 || s === sweeps - 1) {
      yield {
        type: 'progress',
        sweep: s + 1,
        sweeps: sweeps,
        bestEnergy: bestE,
        bestCoords: best.coords,
        energies: folders.map(function (f2) { return f2.energy(); }),
        exchanges: exchanges,
        accepted: accepted,
        proposed: proposed
      };
    }
  }

  var bestFolder = new LatticeFolder(seq, { mode: mode });
  bestFolder._commit(
    best.coords,
    (function () {
      var m = new Map();
      for (var i = 0; i < best.coords.length; i++) {
        m.set(cellKey(best.coords[i][0], best.coords[i][1], best.coords[i][2]), i);
      }
      return m;
    })(),
    best.energy,
    null
  );
  yield {
    type: 'done',
    done: true,
    best: {
      coords: best.coords,
      energy: best.energy,
      contacts: bestFolder.contacts(),
      contactCount: bestFolder.contactCount(),
      radiusOfGyration: bestFolder.radiusOfGyration(),
      sequence: seq,
      mode: mode
    },
    energies: folders.map(function (f3) { return f3.energy(); }),
    trace: trace,
    temps: temps,
    exchanges: exchanges,
    exchangeTries: exchangeTries,
    accepted: accepted,
    proposed: proposed,
    sweeps: sweeps,
    replicas: replicas,
    elapsedMs: Date.now() - started
  };
}

/** Convenience wrapper that drains remcGenerator and returns the result. */
function simulateRemc(options) {
  var it = remcGenerator(options);
  var last;
  while (true) {
    var r = it.next();
    if (r.done) { break; }
    last = r.value;
  }
  return last;
}

/**
 * Simulated annealing: a single walker whose temperature decays geometric
 * (or linearly) from T0 to T1 over the run. Same yield/done protocol as
 * remcGenerator.
 */
function* annealingGenerator(options) {
  var opts = options || {};
  var seq = validateSequence(opts.sequence);
  var n = seq.length;
  var mode = opts.mode === 'hp' ? 'hp' : 'mj';
  var sweeps = opts.sweeps || 800;
  var stepsPerSweep = opts.stepsPerSweep || n;
  var T0 = typeof opts.T0 === 'number' ? opts.T0 : 2.0;
  var T1 = typeof opts.T1 === 'number' ? opts.T1 : 0.05;
  var cooling = opts.cooling === 'linear' ? 'linear' : 'geometric';
  var yieldEvery = opts.yieldEvery || 20;
  var seed = opts.seed === undefined ? 12345 : opts.seed;
  var rng = makeRng(seed);
  var folder = new LatticeFolder(seq, { mode: mode });
  var accepted = 0, proposed = 0;
  var best = { energy: folder.energy(), coords: folder.coords.map(function (c) { return c.slice(); }) };
  var trace = [];
  var started = Date.now();
  for (var s = 0; s < sweeps; s++) {
    var frac = sweeps === 1 ? 1 : s / (sweeps - 1);
    var T = cooling === 'linear' ? T0 + (T1 - T0) * frac : T0 * Math.pow(T1 / T0, frac);
    for (var st = 0; st < stepsPerSweep; st++) {
      var res = folder.mcStep(T, rng);
      proposed++;
      if (res.accepted) { accepted++; }
      if (folder.energy() < best.energy) {
        best = { energy: folder.energy(), coords: folder.coords.map(function (c) { return c.slice(); }) };
      }
    }
    trace.push(best.energy);
    if (s % yieldEvery === 0 || s === sweeps - 1) {
      yield {
        type: 'progress', sweep: s + 1, sweeps: sweeps,
        bestEnergy: best.energy, energy: folder.energy(),
        temperature: T, accepted: accepted, proposed: proposed
      };
    }
  }
  var bf = new LatticeFolder(seq, { mode: mode });
  var occ = new Map();
  for (var i = 0; i < best.coords.length; i++) {
    occ.set(cellKey(best.coords[i][0], best.coords[i][1], best.coords[i][2]), i);
  }
  bf._commit(best.coords, occ, best.energy, null);
  yield {
    type: 'done', done: true,
    best: {
      coords: best.coords, energy: best.energy,
      contacts: bf.contacts(), contactCount: bf.contactCount(),
      radiusOfGyration: bf.radiusOfGyration(), sequence: seq, mode: mode
    },
    energy: folder.energy(),
    trace: trace, accepted: accepted, proposed: proposed,
    sweeps: sweeps, elapsedMs: Date.now() - started
  };
}

function simulateAnnealing(options) {
  var it = annealingGenerator(options);
  var last;
  while (true) {
    var r = it.next();
    if (r.done) { break; }
    last = r.value;
  }
  return last;
}

/** Export a folder conformation as a CA-trace PDB (one ATOM per residue). */
function folderToPdb(folder, options) {
  var opts = options || {};
  var scale = typeof opts.scale === 'number' ? opts.scale : 3.8;
  var atoms = [];
  for (var i = 0; i < folder.n; i++) {
    var c = folder.coords[i];
    atoms.push({
      name: 'CA',
      resName: AA1_TO_3[folder.sequence[i]],
      resSeq: i + 1,
      chainID: opts.chainID || 'A',
      xyz: [c[0] * scale, c[1] * scale, c[2] * scale],
      element: 'C'
    });
  }
  return atomsToPdb(atoms, { title: opts.title || 'fold-lattice CA trace', remarks: opts.remarks });
}

/** Rebuild an all-atom backbone from a lattice conformation. */
function reconstructLatticeBackbone(folder) {
  var angles = latticeBackboneAngles(folder.coords);
  return {
    angles: angles,
    atoms: reconstructBackbone(folder.sequence, angles),
    dihedrals: angles.map(function (a, idx) {
      return {
        resIndex: idx,
        phi: a.phi, psi: a.psi, omega: a.omega,
        region: a.region
      };
    })
  };
}

// Revert support for mcStep: snapshot and restore helpers used by tests.
function snapshot(folder) {
  return {
    coords: folder.coords.map(function (c) { return c.slice(); }),
    occ: new Map(folder.occ),
    energy: folder.energy(),
    contacts: folder.contacts().map(function (c) { return c.slice(); })
  };
}
function restore(folder, snap) {
  folder._commit(snap.coords, snap.occ, snap.energy, snap.contacts);
}

var foldl = {
  VERSION: VERSION,
  AA20: AA20,
  AA1_TO_3: AA1_TO_3,
  AA3_TO_1: AA3_TO_1,
  MJ_ORDER: MJ_ORDER,
  MJ1996: MJ1996,
  MJ_INDEX: MJ_INDEX,
  HYDROPHOBIC: HYDROPHOBIC,
  isHydrophobic: isHydrophobic,
  HP_MATRIX: HP_MATRIX,
  mjPairEnergy: mjPairEnergy,
  makeRng: makeRng,
  vadd: vadd, vsub: vsub, vscale: vscale, vdot: vdot, vcross: vcross,
  vlen: vlen, vdist: vdist, vnorm: vnorm, veq: veq,
  AXIS_VECTORS: AXIS_VECTORS,
  bondAngle: bondAngle,
  dihedral: dihedral,
  placeAtom: placeAtom,
  DEG: DEG, RAD: RAD,
  normalizeDeg: normalizeDeg,
  RAMACHANDRAN_REGIONS: RAMACHANDRAN_REGIONS,
  classifyRamachandran: classifyRamachandran,
  classifyOmega: classifyOmega,
  VDW_RADII: VDW_RADII,
  inferElement: inferElement,
  atomRadius: atomRadius,
  findStericClashes: findStericClashes,
  radiusOfGyration: radiusOfGyration,
  parsePdb: parsePdb,
  parsePdbLine: parsePdbLine,
  atomsToPdb: atomsToPdb,
  pdbAtomNameField: pdbAtomNameField,
  BACKBONE_GEOM: BACKBONE_GEOM,
  reconstructBackbone: reconstructBackbone,
  computeDihedrals: computeDihedrals,
  LATTICE_ANGLE_MAP: LATTICE_ANGLE_MAP,
  latticeBackboneAngles: latticeBackboneAngles,
  validateSequence: validateSequence,
  LatticeFolder: LatticeFolder,
  rotationAboutAxis: rotationAboutAxis,
  geometricTemps: geometricTemps,
  remcGenerator: remcGenerator,
  simulateRemc: simulateRemc,
  annealingGenerator: annealingGenerator,
  simulateAnnealing: simulateAnnealing,
  folderToPdb: folderToPdb,
  reconstructLatticeBackbone: reconstructLatticeBackbone,
  snapshot: snapshot,
  restore: restore
};

/* global module */
if (typeof module === 'object' && module.exports) {
  module.exports = foldl;
} else {
  /* global window */
  window.foldl = foldl;
}
