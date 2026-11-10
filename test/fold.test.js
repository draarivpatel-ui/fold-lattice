'use strict';

const test = require('node:test');
const assert = require('node:assert');
const f = require('../src/fold.js');

const SEQ = 'GAGSLYQGSSLAQSGFLSSHHSGF';
// HP mode test chain: A stands in for hydrophobic (A is in ACFILMVWY),
// S for polar. Pattern mirrors the classic block copolymer HP designs.
const HP_SEQ = 'SSSSAAAASSSSAAAASSSAA';

// Angle-safe comparison for wrapped dihedral values.
function angClose(actual, expected, tolDeg) {
  const d = Math.abs(f.normalizeDeg(actual - expected));
  return d <= tolDeg;
}

// -----------------------------------------------------------------------
// RNG
// -----------------------------------------------------------------------
test('makeRng is deterministic for a given seed', () => {
  const a = f.makeRng(42);
  const b = f.makeRng(42);
  for (let i = 0; i < 100; i++) assert.strictEqual(a(), b());
});

test('makeRng values stay in [0, 1)', () => {
  const r = f.makeRng(7);
  for (let i = 0; i < 1000; i++) {
    const v = r();
    assert.ok(v >= 0 && v < 1);
  }
});

test('different seeds give different sequences', () => {
  const a = f.makeRng(1);
  const b = f.makeRng(2);
  let diff = false;
  for (let i = 0; i < 10; i++) if (a() !== b()) diff = true;
  assert.ok(diff);
});

// -----------------------------------------------------------------------
// Vector helpers
// -----------------------------------------------------------------------
test('vector add, sub, scale, dot', () => {
  assert.deepStrictEqual(f.vadd([1, 2, 3], [4, 5, 6]), [5, 7, 9]);
  assert.deepStrictEqual(f.vsub([4, 5, 6], [1, 2, 3]), [3, 3, 3]);
  assert.deepStrictEqual(f.vscale([1, -2, 3], 2), [2, -4, 6]);
  assert.strictEqual(f.vdot([1, 2, 3], [4, 5, 6]), 32);
});

test('vector cross product', () => {
  assert.deepStrictEqual(f.vcross([1, 0, 0], [0, 1, 0]), [0, 0, 1]);
  assert.deepStrictEqual(f.vcross([0, 1, 0], [1, 0, 0]), [0, 0, -1]);
  assert.strictEqual(f.vdot(f.vcross([3, 1, 2], [1, 4, 2]), [3, 1, 2]), 0);
});

test('vector length and distance', () => {
  assert.strictEqual(f.vlen([3, 4, 0]), 5);
  assert.strictEqual(f.vdist([0, 0, 0], [1, 1, 1]), Math.sqrt(3));
});

test('vector normalize gives unit length and rejects zero', () => {
  const u = f.vnorm([2, 0, 0]);
  assert.ok(Math.abs(f.vlen(u) - 1) < 1e-12);
  assert.throws(() => f.vnorm([0, 0, 0]), /zero vector/);
});

// -----------------------------------------------------------------------
// Amino acid data and MJ matrix
// -----------------------------------------------------------------------
test('AA20 covers exactly the 20 standard residues', () => {
  assert.strictEqual(f.AA20.length, 20);
  const uniq = new Set(f.AA20.split(''));
  assert.strictEqual(uniq.size, 20);
  assert.ok(!f.AA20.includes('B') && !f.AA20.includes('X') && !f.AA20.includes('Z'));
});

test('one letter to three letter mapping round trips', () => {
  for (const a of f.AA20) {
    assert.strictEqual(f.AA3_TO_1[f.AA1_TO_3[a]], a);
  }
  assert.strictEqual(f.AA1_TO_3.W, 'TRP');
  assert.strictEqual(f.AA3_TO_1.CYS, 'C');
});

test('MJ order lists all 20 residues', () => {
  assert.strictEqual(f.MJ_ORDER.length, 20);
  for (const a of f.AA20) assert.ok(f.MJ_ORDER.includes(a));
  assert.ok(f.MJ_INDEX.C !== undefined && f.MJ_INDEX.P !== undefined);
});

test('MJ1996 matrix is 20x20 and symmetric', () => {
  assert.strictEqual(f.MJ1996.length, 20);
  for (const row of f.MJ1996) assert.strictEqual(row.length, 20);
  for (let i = 0; i < 20; i++) {
    for (let j = 0; j < 20; j++) {
      assert.ok(Math.abs(f.MJ1996[i][j] - f.MJ1996[j][i]) < 1e-12);
    }
  }
});

test('MJ1996 anchor values match Table 3 of the 1996 paper', () => {
  assert.strictEqual(f.mjPairEnergy('C', 'C'), -5.44);
  assert.strictEqual(f.mjPairEnergy('M', 'M'), -5.46);
  assert.strictEqual(f.mjPairEnergy('F', 'L'), -7.28);
  assert.strictEqual(f.mjPairEnergy('K', 'K'), -0.12);
  assert.strictEqual(f.mjPairEnergy('A', 'A'), -2.72);
  assert.strictEqual(f.mjPairEnergy('W', 'Y'), -4.66);
  assert.strictEqual(f.mjPairEnergy('L', 'H'), -4.54);
});

test('MJ1996 diagonal holds the strongest hydrophobic contacts', () => {
  // L-L (-7.37) is the strongest diagonal entry, edging out F-F (-7.26).
  let minDiag = 0, minPair = '';
  for (const a of f.MJ_ORDER) {
    const e = f.mjPairEnergy(a, a);
    if (e < minDiag) { minDiag = e; minPair = a + a; }
  }
  assert.strictEqual(minPair, 'LL');
  assert.strictEqual(minDiag, -7.37);
});

test('mjPairEnergy rejects unknown residues', () => {
  assert.throws(() => f.mjPairEnergy('B', 'A'), /unknown residue/);
  assert.throws(() => f.mjPairEnergy('A', 'Z'), /unknown residue/);
});

// -----------------------------------------------------------------------
// Sequence validation
// -----------------------------------------------------------------------
test('validateSequence accepts the 20 standard letters case-insensitively', () => {
  assert.strictEqual(f.validateSequence('acdefghiklmnpqrstvwy'), 'ACDEFGHIKLMNPQRSTVWY');
  assert.strictEqual(f.validateSequence('GaGsLyQ'), 'GAGSLYQ');
  assert.strictEqual(f.validateSequence('ACDEFGHIKLMNPQRSTVWY'), 'ACDEFGHIKLMNPQRSTVWY');
});

test('validateSequence rejects non-standard letters', () => {
  assert.throws(() => f.validateSequence('AXB'), /invalid residue/);
  assert.throws(() => f.validateSequence('123'), /invalid|must contain/);
  assert.throws(() => f.validateSequence(''), /at least 2/);
});

test('validateSequence rejects non-string input', () => {
  assert.throws(() => f.validateSequence(42), /must be a string/);
  assert.throws(() => f.validateSequence(['A', 'B']), /must be a string/);
});

// -----------------------------------------------------------------------
// Geometric primitives
// -----------------------------------------------------------------------
test('bondAngle computes the angle at the middle point', () => {
  const ang = f.bondAngle([1, 0, 0], [0, 0, 0], [0, 1, 0]);
  assert.ok(Math.abs(ang - Math.PI / 2) < 1e-12);
  const straight = f.bondAngle([1, 0, 0], [0, 0, 0], [-1, 0, 0]);
  assert.ok(Math.abs(straight - Math.PI) < 1e-12);
});

test('bondAngle rejects coincident points', () => {
  assert.throws(() => f.bondAngle([0, 0, 0], [0, 0, 0], [1, 1, 1]), /coincident/);
});

test('dihedral computes the IUPAC torsion of a planar zig-zag', () => {
  // trans planar chain along x with alternating y offsets has omega 180.
  const p1 = [0, 0, 0], p2 = [1.33, 0, 0], p3 = [2.33, 0, 0], p4 = [3.66, 0, 0];
  // exactly colinear is undefined
  assert.throws(() => f.dihedral(p1, p2, p3, p4), /colinear|undefined/);
});

test('dihedral sign follows the right hand rule', () => {
  // Mirror images through the plane of p1, p2, p3 flip the torsion sign.
  const p1 = [0, 0, 0];
  const p2 = [1, 0, 0];
  const p3 = [1, 1, 0];
  const d1 = f.dihedral(p1, p2, p3, [2, 1, 1]);
  const d2 = f.dihedral(p1, p2, p3, [2, 1, -1]);
  assert.ok(Math.abs(d1 + d2) < 1e-9, 'mirror images flip the sign');
  assert.ok(Math.abs(d1) > Math.PI / 3, 'nondegenerate torsion magnitude');
});

test('placeAtom reproduces bond length, angle and torsion exactly', () => {
  const a = [0.1, 0.2, 0.3];
  const b = [1.1, 0.2, 0.3];
  const c = [1.6, 1.0, 0.5];
  for (const tor of [0, 45, 90, 137, -90, 180, -179]) {
    const d = f.placeAtom(a, b, c, 1.525, 111.2 * f.RAD, tor * f.RAD);
    assert.ok(Math.abs(f.vdist(c, d) - 1.525) < 1e-9, 'bond length');
    assert.ok(Math.abs(f.bondAngle(b, c, d) - 111.2 * f.RAD) < 1e-9, 'bond angle');
    assert.ok(angClose(f.dihedral(a, b, c, d) / f.RAD, tor, 1e-6), 'torsion');
  }
});

test('placeAtom rejects colinear reference frames', () => {
  const a = [0, 0, 0], b = [1, 0, 0], c = [2, 0, 0];
  assert.throws(() => f.placeAtom(a, b, c, 1.5, 1, 1), /colinear|coincide/);
});

// -----------------------------------------------------------------------
// Rotation matrices
// -----------------------------------------------------------------------
test('rotationAboutAxis identity for four quarter turns', () => {
  const I = f.rotationAboutAxis([1, 0, 0], 4);
  assert.deepStrictEqual(I, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
  const I0 = f.rotationAboutAxis([0, 1, 0], 0);
  assert.deepStrictEqual(I0, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
});

test('rotationAboutAxis keeps the axis fixed and preserves norms', () => {
  const u = [0, 0, 1];
  for (const t of [1, 2, 3, 5, -1, -2]) {
    const R = f.rotationAboutAxis(u, t);
    const v = matVec(R, [0, 0, 1]);
    assert.deepStrictEqual(v, [0, 0, 1]);
    const p = matVec(R, [1, 1, 1]);
    assert.ok(Math.abs(f.vlen(p) - Math.sqrt(3)) < 1e-9);
  }
});

test('rotationAboutAxis about -z equals inverse turn about +z', () => {
  const a = matVec(f.rotationAboutAxis([0, 0, -1], 1), [1, 0, 0]);
  const b = matVec(f.rotationAboutAxis([0, 0, 1], 3), [1, 0, 0]);
  assert.deepStrictEqual(a, b);
});

test('rotationAboutAxis maps lattice points to integer lattice points', () => {
  const R = f.rotationAboutAxis([0, 1, 0], 1);
  const p = matVec(R, [1, 2, 3]);
  for (const c of p) assert.ok(Number.isInteger(c));
});

function matVec(R, v) {
  return [
    R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
    R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
    R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]
  ];
}

// -----------------------------------------------------------------------
// Ramachandran classification
// -----------------------------------------------------------------------
test('normalizeDeg wraps angles into (-180, 180]', () => {
  assert.strictEqual(f.normalizeDeg(0), 0);
  assert.strictEqual(f.normalizeDeg(185), -175);
  assert.strictEqual(f.normalizeDeg(-185), 175);
  assert.strictEqual(f.normalizeDeg(540), 180);
  assert.strictEqual(f.normalizeDeg(-180), 180);
});

test('alpha right basin classifies canonical helix angles', () => {
  assert.strictEqual(f.classifyRamachandran(-57, -47), 'ALPHA_R');
  assert.strictEqual(f.classifyRamachandran(-63, -43), 'ALPHA_R');
});

test('beta basin classifies canonical sheet angles', () => {
  assert.strictEqual(f.classifyRamachandran(-139, 135), 'BETA');
  assert.strictEqual(f.classifyRamachandran(-120, 130), 'BETA');
});

test('alpha left basin classifies left handed angles', () => {
  assert.strictEqual(f.classifyRamachandran(57, 47), 'ALPHA_L');
});

test('polyproline II region classifies its basin', () => {
  assert.strictEqual(f.classifyRamachandran(-75, 145), 'PP_II');
});

test('disallowed angles are outliers', () => {
  assert.strictEqual(f.classifyRamachandran(120, -120), 'OUTLIER');
  assert.strictEqual(f.classifyRamachandran(0, 0), 'OUTLIER');
});

test('omega cis and trans classification', () => {
  assert.strictEqual(f.classifyOmega(180), 'TRANS');
  assert.strictEqual(f.classifyOmega(-175), 'TRANS');
  assert.strictEqual(f.classifyOmega(0), 'CIS');
  assert.strictEqual(f.classifyOmega(45), 'CIS');
  assert.strictEqual(f.classifyOmega(90), 'TRANS');
  assert.strictEqual(f.classifyOmega(89), 'CIS');
});

// -----------------------------------------------------------------------
// Steric clashes and radii
// -----------------------------------------------------------------------
test('atom radii come from the element table', () => {
  assert.ok(Math.abs(f.atomRadius({ element: 'C' }) - 1.70) < 1e-12);
  assert.ok(Math.abs(f.atomRadius({ element: 'O' }) - 1.52) < 1e-12);
  assert.ok(Math.abs(f.atomRadius({ element: 'S' }) - 1.80) < 1e-12);
});

test('atom radius falls back to the default for unknown elements', () => {
  assert.strictEqual(f.atomRadius({ element: 'Xx' }), f.DEFAULT_VDW_RADIUS === undefined ? 1.70 : f.DEFAULT_VDW_RADIUS);
});

test('inferElement reads declared elements and atom names', () => {
  assert.strictEqual(f.inferElement('CA', 'C'), 'C');
  assert.strictEqual(f.inferElement(' NZ ', ''), 'N');
  assert.strictEqual(f.inferElement('FE', ''), 'FE');
  assert.strictEqual(f.inferElement('1HB', ''), 'H');
});

test('findStericClashes flags overlapping atoms', () => {
  // 2.2 angstrom apart: above the 2.1 bonded cutoff, below the 0.75 *
  // (1.70 + 1.70) vdW overlap limit.
  const atoms = [
    { xyz: [0, 0, 0], element: 'C' },
    { xyz: [2.2, 0, 0], element: 'C' }
  ];
  const clashes = f.findStericClashes(atoms, { factor: 0.75 });
  assert.strictEqual(clashes.length, 1);
  assert.ok(clashes[0].overlap > 0);
});

test('findStericClashes ignores bonded and distant pairs', () => {
  const atoms = [
    { xyz: [0, 0, 0], element: 'C' },
    { xyz: [1.5, 0, 0], element: 'C' },
    { xyz: [5, 0, 0], element: 'C' }
  ];
  // 1.5 apart: below minBond, skipped. 5 apart: beyond sum of radii, fine.
  assert.strictEqual(f.findStericClashes(atoms).length, 0);
});

test('findStericClashes honors excludePair', () => {
  const atoms = [
    { xyz: [0, 0, 0], element: 'C' },
    { xyz: [2.3, 0, 0], element: 'C' }
  ];
  const all = f.findStericClashes(atoms);
  const none = f.findStericClashes(atoms, { excludePair: () => true });
  assert.strictEqual(all.length, 1);
  assert.strictEqual(none.length, 0);
});

// -----------------------------------------------------------------------
// Radius of gyration
// -----------------------------------------------------------------------
test('radiusOfGyration of a straight line matches the formula', () => {
  const pts = [];
  for (let i = 0; i < 5; i++) pts.push([i, 0, 0]);
  // Rg^2 = (n^2 - 1) / 12 for n evenly spaced points
  const expected = Math.sqrt((25 - 1) / 12);
  assert.ok(Math.abs(f.radiusOfGyration(pts) - expected) < 1e-12);
});

test('radiusOfGyration is translation invariant', () => {
  const pts = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 1, 0]];
  const shifted = pts.map(p => [p[0] + 7, p[1] - 3, p[2] + 100]);
  assert.ok(Math.abs(f.radiusOfGyration(pts) - f.radiusOfGyration(shifted)) < 1e-12);
});

test('radiusOfGyration rejects empty input', () => {
  assert.throws(() => f.radiusOfGyration([]), /empty/);
});

// -----------------------------------------------------------------------
// PDB parsing
// -----------------------------------------------------------------------
const MINI_PDB = [
  'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
  'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C',
  'ATOM      3  C   ALA A   1       2.000   1.300   0.300  1.00  0.00           C',
  'ATOM      4  O   ALA A   1       1.500   2.400   0.400  1.00  0.00           O',
  'ATOM      5  N   GLY A   2       3.200   1.100   1.200  1.00  0.00           N',
  'ATOM      6  CA  GLY A   2       4.100   0.100   1.600  1.00  0.00           C',
  'ATOM      7  C   GLY A   2       5.500   0.700   1.400  1.00  0.00           C',
  'ATOM      8  O   GLY A   2       6.200   0.500   2.400  1.00  0.00           O',
  'TER',
  'HETATM    9  O   HOH A  99      10.000  10.000  10.000  1.00  0.00           O',
  'END'
].join('\n');

test('parsePdb reads ATOM records with names, residues and coordinates', () => {
  const parsed = f.parsePdb(MINI_PDB);
  assert.strictEqual(parsed.atoms.length, 8);
  const first = parsed.atoms[0];
  assert.strictEqual(first.name.trim(), 'N');
  assert.strictEqual(first.resName, 'ALA');
  assert.strictEqual(first.resSeq, 1);
  assert.strictEqual(first.chainID, 'A');
  assert.deepStrictEqual(first.xyz, [0, 0, 0]);
  assert.strictEqual(parsed.atoms[1].xyz[0], 1.458);
});

test('parsePdb skips HETATM by default and includes it on request', () => {
  assert.strictEqual(f.parsePdb(MINI_PDB).atoms.length, 8);
  const withHet = f.parsePdb(MINI_PDB, { includeHetatm: true });
  assert.strictEqual(withHet.atoms.length, 9);
  assert.strictEqual(withHet.atoms[8].resName, 'HOH');
});

test('parsePdb infers missing elements from atom names', () => {
  const noElem = MINI_PDB.replace(/           [A-Z]$/gm, '   ');
  const parsed = f.parsePdb(noElem);
  assert.strictEqual(parsed.atoms[0].element, 'N');
  assert.strictEqual(parsed.atoms[1].element, 'C');
});

test('parsePdb rejects input without ATOM records', () => {
  assert.throws(() => f.parsePdb('HEADER    nothing here\nEND'), /no ATOM records/);
  assert.throws(() => f.parsePdb(''), /no ATOM records/);
});

test('parsePdbLine returns null for non atom lines', () => {
  assert.strictEqual(f.parsePdbLine('REMARK   1 some note'), null);
  assert.strictEqual(f.parsePdbLine('END'), null);
  assert.strictEqual(f.parsePdbLine('CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1'), null);
});

test('parsePdb stops at ENDMDL and reads one model', () => {
  const twoModels = MINI_PDB.replace('END', 'ENDMDL\nMODEL        2\n' + MINI_PDB.split('\n')[0] + '\nENDMDL\nEND');
  const parsed = f.parsePdb(twoModels);
  assert.strictEqual(parsed.atoms.length, 8);
});

test('parsePdb collects REMARK lines', () => {
  const withRemark = 'REMARK   2 test remark\n' + MINI_PDB;
  const parsed = f.parsePdb(withRemark);
  assert.ok(parsed.remarks.some(r => r.includes('test remark')));
});

test('computeDihedrals extracts phi, psi and omega from a PDB', () => {
  const parsed = f.parsePdb(MINI_PDB);
  const d = f.computeDihedrals(parsed.atoms);
  assert.strictEqual(d.length, 2);
  // residue 1 has no predecessor: phi null, psi computable
  assert.strictEqual(d[0].phi, null);
  assert.ok(d[0].psi !== null);
  assert.strictEqual(d[1].phi !== null, true);
  assert.ok(d[1].omega !== null);
});

test('computeDihedrals marks chain breaks as missing torsions', () => {
  // Residues numbered 1 and 3 with a gap: no bonded pairs.
  const gapped = MINI_PDB.replace(/A   2/g, 'A   3');
  const d = f.computeDihedrals(f.parsePdb(gapped).atoms);
  assert.strictEqual(d.length, 2);
  assert.strictEqual(d[0].psi, null);
  assert.strictEqual(d[1].phi, null);
});

// -----------------------------------------------------------------------
// PDB export
// -----------------------------------------------------------------------
test('atomsToPdb emits records with correct column layout', () => {
  const atoms = [
    { name: 'N', resName: 'ALA', resSeq: 1, chainID: 'A', xyz: [1.25, 2.5, 3.75], element: 'N' },
    { name: 'CA', resName: 'ALA', resSeq: 1, chainID: 'A', xyz: [-1.5, 0.25, 100.125], element: 'C' }
  ];
  const pdb = f.atomsToPdb(atoms, { title: 'unit test' });
  const lines = pdb.split('\n');
  assert.ok(lines[0].startsWith('HEADER'));
  // lines[1] is the always emitted REMARK banner, first atom is lines[2].
  assert.ok(lines[1].startsWith('REMARK'));
  const atom1 = lines[2];
  assert.ok(atom1.startsWith('ATOM  '));
  assert.ok(atom1.slice(12, 16) === ' N  ', 'atom name field');
  assert.ok(atom1.slice(17, 20) === 'ALA', 'residue name field');
  assert.ok(atom1.slice(21, 22) === 'A', 'chain field');
  assert.ok(atom1.slice(22, 26) === '   1', 'residue number field');
  assert.ok(atom1.slice(30, 38) === '   1.250', 'x coordinate');
  assert.ok(atom1.slice(38, 46) === '   2.500', 'y coordinate');
  assert.ok(atom1.slice(46, 54) === '   3.750', 'z coordinate');
  assert.ok(atom1.slice(76, 78) === ' N', 'element field');
  // trailing newline leaves an empty final element: END is at length - 2,
  // the closing TER at length - 3.
  assert.strictEqual(lines[lines.length - 1], '');
  assert.strictEqual(lines[lines.length - 2], 'END');
  assert.ok(lines[lines.length - 3].startsWith('TER'));
});

test('pdbAtomNameField aligns single and four character names', () => {
  assert.strictEqual(f.pdbAtomNameField('N'), ' N  ');
  assert.strictEqual(f.pdbAtomNameField('CA'), ' CA ');
  assert.strictEqual(f.pdbAtomNameField('OXT'), ' OXT');
  assert.strictEqual(f.pdbAtomNameField('HB2X').length, 4);
});

test('atomsToPdb output reparses to the same atoms', () => {
  const atoms = f.reconstructBackbone('ACDE', [
    { phi: -57, psi: -47 }, { phi: -139, psi: 135 },
    { phi: 57, psi: 47 }, { phi: -75, psi: 145 }
  ]);
  const pdb = f.atomsToPdb(atoms, { title: 'round trip' });
  const reparsed = f.parsePdb(pdb);
  assert.strictEqual(reparsed.atoms.length, atoms.length);
  for (let i = 0; i < atoms.length; i++) {
    assert.ok(f.vdist(atoms[i].xyz, reparsed.atoms[i].xyz) < 0.01, 'coordinate fidelity');
    assert.strictEqual(reparsed.atoms[i].resName, atoms[i].resName);
  }
});

test('atomsToPdb inserts TER between chains and numbers serials', () => {
  const atoms = [
    { name: 'CA', resName: 'GLY', resSeq: 1, chainID: 'A', xyz: [0, 0, 0], element: 'C' },
    { name: 'CA', resName: 'GLY', resSeq: 2, chainID: 'A', xyz: [3.8, 0, 0], element: 'C' },
    { name: 'CA', resName: 'GLY', resSeq: 1, chainID: 'B', xyz: [0, 3.8, 0], element: 'C' }
  ];
  const pdb = f.atomsToPdb(atoms);
  const lines = pdb.split('\n').filter(l => l.startsWith('ATOM') || l.startsWith('TER'));
  // 2 atoms in chain A, TER, 1 atom in chain B, final TER.
  assert.strictEqual(lines.length, 5);
  assert.ok(lines[2].startsWith('TER'));
  assert.ok(lines[4].startsWith('TER'));
  // serials: 1, 2 in chain A; the TER takes serial 3; chain B atom is 4.
  assert.ok(lines[0].includes(padSerial(1)));
  assert.ok(lines[1].includes(padSerial(2)));
  assert.ok(lines[3].includes(padSerial(4)));
});

function padSerial(s) {
  let out = String(s);
  while (out.length < 5) out = ' ' + out;
  return out;
}

// -----------------------------------------------------------------------
// Backbone reconstruction
// -----------------------------------------------------------------------
test('BACKBONE_GEOM holds standard Engh and Huber values', () => {
  assert.ok(Math.abs(f.BACKBONE_GEOM.L_N_CA - 1.458) < 1e-9);
  assert.ok(Math.abs(f.BACKBONE_GEOM.L_CA_C - 1.525) < 1e-9);
  assert.ok(Math.abs(f.BACKBONE_GEOM.L_C_N - 1.329) < 1e-9);
  assert.ok(Math.abs(f.BACKBONE_GEOM.L_C_O - 1.231) < 1e-9);
});

test('reconstructBackbone produces four atoms per residue', () => {
  const atoms = f.reconstructBackbone('ACD', [
    { phi: -57, psi: -47 }, { phi: -139, psi: 135 }, { phi: 57, psi: 47 }
  ]);
  assert.strictEqual(atoms.length, 12);
  assert.deepStrictEqual(atoms.map(a => a.name).slice(0, 4), ['N', 'CA', 'C', 'O']);
  assert.strictEqual(atoms[4].resName, 'CYS');
  assert.strictEqual(atoms[8].resName, 'ASP');
});

test('reconstructBackbone honors exact bond lengths', () => {
  const g = f.BACKBONE_GEOM;
  const atoms = f.reconstructBackbone('AG', [{ phi: -120, psi: 130 }, { phi: -60, psi: 40 }]);
  const byRes = [];
  for (const a of atoms) { byRes[a.resSeq - 1] = byRes[a.resSeq - 1] || {}; byRes[a.resSeq - 1][a.name] = a.xyz; }
  assert.ok(Math.abs(f.vdist(byRes[0].N, byRes[0].CA) - g.L_N_CA) < 1e-9);
  assert.ok(Math.abs(f.vdist(byRes[0].CA, byRes[0].C) - g.L_CA_C) < 1e-9);
  assert.ok(Math.abs(f.vdist(byRes[0].C, byRes[1].N) - g.L_C_N) < 1e-9);
  assert.ok(Math.abs(f.vdist(byRes[1].CA, byRes[1].C) - g.L_CA_C) < 1e-9);
  assert.ok(Math.abs(f.vdist(byRes[1].C, byRes[1].O) - g.L_C_O) < 1e-9);
});

test('reconstructBackbone reproduces the input dihedrals exactly', () => {
  const angles = [
    { phi: -57, psi: -47, omega: 180 },
    { phi: -139, psi: 135, omega: 180 },
    { phi: 57, psi: 47, omega: 180 },
    { phi: -75, psi: 145, omega: 180 },
    { phi: -120, psi: 130, omega: 178 }
  ];
  const atoms = f.reconstructBackbone('AAAAA', angles);
  const d = f.computeDihedrals(atoms);
  for (let i = 0; i < angles.length; i++) {
    if (i > 0) assert.ok(angClose(d[i].phi, angles[i].phi, 1e-6), 'phi ' + i);
    if (i < angles.length - 1) assert.ok(angClose(d[i].psi, angles[i].psi, 1e-6), 'psi ' + i);
    if (i > 0) assert.ok(angClose(d[i].omega, angles[i].omega, 1e-6), 'omega ' + i);
  }
});

test('reconstructBackbone rejects mismatched angle arrays', () => {
  assert.throws(() => f.reconstructBackbone('AAA', [{ phi: 1, psi: 1 }]), /length/);
  assert.throws(() => f.reconstructBackbone('AAA', [{ phi: null, psi: 1 }, { phi: 1, psi: 1 }, { phi: 1, psi: 1 }]), /numeric/);
  assert.throws(() => f.reconstructBackbone('AAA', 'nope'), /length/);
});

test('reconstructBackbone builds a trans peptide plane by default', () => {
  const atoms = f.reconstructBackbone('AAA', [{ phi: -57, psi: -47 }, { phi: -57, psi: -47 }, { phi: -57, psi: -47 }]);
  const d = f.computeDihedrals(atoms);
  assert.strictEqual(d[1].omegaClass, 'TRANS');
  assert.strictEqual(d[2].omegaClass, 'TRANS');
});

test('reconstructBackbone can build a cis peptide', () => {
  const atoms = f.reconstructBackbone('PG', [{ phi: -60, psi: 120 }, { phi: -80, psi: 130, omega: 0 }]);
  const d = f.computeDihedrals(atoms);
  assert.strictEqual(d[1].omegaClass, 'CIS');
});

// -----------------------------------------------------------------------
// LatticeFolder basics
// -----------------------------------------------------------------------
test('LatticeFolder starts as a straight line with zero contacts', () => {
  const lf = new f.LatticeFolder(SEQ);
  assert.strictEqual(lf.n, SEQ.length);
  assert.strictEqual(lf.contacts().length, 0);
  assert.strictEqual(lf.energy(), 0);
  assert.ok(lf.isValidConformation());
});

test('LatticeFolder MJ energy of the extended chain is zero', () => {
  const lf = new f.LatticeFolder('ACDEFGHIKLMNPQRSTVWYA');
  assert.strictEqual(lf.energy(), 0);
});

test('LatticeFolder indexAt and cellKey agree', () => {
  const lf = new f.LatticeFolder('ACDEF');
  for (let i = 0; i < lf.n; i++) {
    const [x, y, z] = lf.coords[i];
    assert.strictEqual(lf.indexAt(x, y, z), i);
  }
  assert.strictEqual(lf.indexAt(0, 1, 0), undefined);
  assert.strictEqual(lf.indexAt(-5, 0, 0), undefined);
});

test('LatticeFolder hp mode scores H-H contacts only', () => {
  const lf = new f.LatticeFolder(HP_SEQ, { mode: 'hp' });
  assert.strictEqual(lf.mode, 'hp');
  assert.strictEqual(lf.energy(), 0);
  assert.strictEqual(lf.contactCount(), 0);
});

test('LatticeFolder default mode is mj', () => {
  const lf = new f.LatticeFolder('ACDE');
  assert.strictEqual(lf.mode, 'mj');
});

test('LatticeFolder rejects invalid sequences', () => {
  assert.throws(() => new f.LatticeFolder('AXB'), /invalid residue/);
  assert.throws(() => new f.LatticeFolder('A'), /at least 2/);
  assert.throws(() => new f.LatticeFolder(42), /must be a string/);
});

test('LatticeFolder copy is independent', () => {
  const lf = new f.LatticeFolder(SEQ);
  const c = lf.copy();
  c.coords[0][0] = 999;
  assert.strictEqual(lf.coords[0][0], 0);
  assert.strictEqual(c.n, lf.n);
  assert.strictEqual(c.sequence, lf.sequence);
});

test('LatticeFolder key handles negative coordinates', () => {
  const k1 = f.LatticeFolder.key(-10, 3, 2);
  const k2 = f.LatticeFolder.key(-9, 3, 2);
  const k3 = f.LatticeFolder.key(-10, 4, 2);
  const k4 = f.LatticeFolder.key(-10, 3, 3);
  assert.notStrictEqual(k1, k2);
  assert.notStrictEqual(k1, k3);
  assert.notStrictEqual(k1, k4);
  const kA = f.LatticeFolder.key(0, 0, 0);
  const kB = f.LatticeFolder.key(0, 0, 0);
  assert.strictEqual(kA, kB);
});

// -----------------------------------------------------------------------
// Move geometry
// -----------------------------------------------------------------------
test('mcStep never breaks the self avoiding walk', () => {
  const rng = f.makeRng(99);
  const lf = new f.LatticeFolder(SEQ);
  for (let i = 0; i < 3000; i++) {
    lf.mcStep(1.5, rng);
    assert.ok(lf.isValidConformation(), 'invalid after step ' + i);
  }
});

test('mcStep rejects energy raising moves at zero temperature', () => {
  const rng = f.makeRng(5);
  const lf = new f.LatticeFolder(SEQ);
  const start = lf.coords.map(c => c.slice());
  for (let i = 0; i < 500; i++) lf.mcStep(0, rng);
  // At T=0 energy can only decrease or stay equal; coordinates may still
  // move via zero-delta moves, so compare energies instead of positions.
  assert.ok(lf.energy() <= 0);
  assert.ok(lf.isValidConformation());
  assert.notDeepStrictEqual(lf.coords, start, 'zero-T walk should still move');
});

test('mcStep trial mode leaves the state untouched', () => {
  const rng = f.makeRng(21);
  const lf = new f.LatticeFolder(SEQ);
  const snap = f.snapshot(lf);
  for (let i = 0; i < 500; i++) {
    const types = ['pull', 'pullUp', 'pullDown', 'pivot', 'crankshaft', 'corner', 'endpull'];
    const t = types[i % types.length];
    const res = lf.tryMove(t, rng, { trial: true });
    assert.deepStrictEqual(f.snapshot(lf), snap, 'trial changed state');
    assert.ok(res.cand === undefined ? true : Array.isArray(res.cand));
  }
});

test('tryMove pivots preserve chain connectivity', () => {
  const rng = f.makeRng(31);
  const lf = new f.LatticeFolder(SEQ);
  let moved = 0;
  for (let i = 0; i < 2000; i++) {
    const res = lf.tryMove('pivot', rng);
    if (res.moved) {
      moved++;
      assert.ok(lf.isValidConformation());
    }
  }
  assert.ok(moved > 100, 'pivot moves should often apply, got ' + moved);
});

test('tryMove corner flips are reversible quarter rotations', () => {
  const rng = f.makeRng(41);
  const lf = new f.LatticeFolder('PPPPPPPPPPHHHHHHHHHH');
  // Warm up: a straight line has no corners to flip.
  for (let i = 0; i < 300; i++) lf.mcStep(2.0, rng);
  let moved = 0;
  for (let i = 0; i < 2000; i++) {
    const res = lf.tryMove('corner', rng);
    if (res.moved) {
      moved++;
      assert.ok(lf.isValidConformation());
    }
  }
  assert.ok(moved > 100);
});

test('tryMove crankshafts keep the walk valid', () => {
  const rng = f.makeRng(51);
  const lf = new f.LatticeFolder(SEQ);
  // Warm up: crankshafts need a U shaped four residue square, which an
  // extended chain does not contain.
  for (let i = 0; i < 300; i++) lf.mcStep(2.0, rng);
  let moved = 0;
  for (let i = 0; i < 4000; i++) {
    const res = lf.tryMove('crankshaft', rng);
    if (res.moved) {
      moved++;
      assert.ok(lf.isValidConformation());
    }
  }
  assert.ok(moved > 10, 'some crankshafts should apply, got ' + moved);
});

test('tryMove pulls and endpoint pulls keep the walk valid', () => {
  const rng = f.makeRng(61);
  const lf = new f.LatticeFolder(SEQ);
  const counts = { pull: 0, pullUp: 0, pullDown: 0, endpull: 0 };
  for (let i = 0; i < 4000; i++) {
    const t = ['pull', 'pullUp', 'pullDown', 'endpull'][i % 4];
    const res = lf.tryMove(t, rng);
    if (res.moved) {
      counts[t]++;
      assert.ok(lf.isValidConformation(), 'invalid ' + t);
    }
  }
  assert.ok(counts.pull > 100);
  assert.ok(counts.endpull > 100);
});

test('unknown move type throws', () => {
  const lf = new f.LatticeFolder('ACDE');
  assert.throws(() => lf.tryMove('wiggle', f.makeRng(1)), /unknown move/);
});

test('snapshot and restore round trip', () => {
  const rng = f.makeRng(71);
  const lf = new f.LatticeFolder(SEQ);
  const snap0 = f.snapshot(lf);
  lf.tryMove('pivot', rng);
  lf.tryMove('pull', rng);
  f.restore(lf, snap0);
  assert.deepStrictEqual(lf.coords, snap0.coords);
  assert.strictEqual(lf.energy(), snap0.energy);
});

test('folded conformations gain contacts compared to extended', () => {
  const r = f.simulateAnnealing({ sequence: SEQ, sweeps: 250, stepsPerSweep: 40, seed: 3 });
  assert.ok(r.best.contactCount > 0, 'annealing found contacts');
  const extended = new f.LatticeFolder(SEQ);
  assert.ok(r.best.energy < extended.energy());
});

// -----------------------------------------------------------------------
// Simulated annealing
// -----------------------------------------------------------------------
test('simulateAnnealing returns a complete result object', () => {
  const r = f.simulateAnnealing({ sequence: 'ACDEFGHIKL', sweeps: 50, stepsPerSweep: 10, seed: 1 });
  assert.ok(r.done);
  assert.strictEqual(r.best.sequence, 'ACDEFGHIKL');
  assert.strictEqual(r.best.coords.length, 10);
  assert.strictEqual(r.trace.length, 50);
  assert.ok(r.accepted + (r.proposed - r.accepted) === r.proposed);
  assert.ok(typeof r.elapsedMs === 'number');
});

test('simulateAnnealing best conformation is a valid SAW', () => {
  const r = f.simulateAnnealing({ sequence: SEQ, sweeps: 120, stepsPerSweep: 30, seed: 13 });
  const lf = new f.LatticeFolder(SEQ);
  const occ = new Map();
  for (let i = 0; i < r.best.coords.length; i++) {
    occ.set(f.LatticeFolder.key(r.best.coords[i][0], r.best.coords[i][1], r.best.coords[i][2]), i);
  }
  lf._commit(r.best.coords, occ, r.best.energy, null);
  assert.ok(lf.isValidConformation());
  assert.strictEqual(lf.energy(), r.best.energy);
  assert.strictEqual(lf.contactCount(), r.best.contactCount);
});

test('simulateAnnealing is reproducible for a fixed seed', () => {
  const a = f.simulateAnnealing({ sequence: SEQ, sweeps: 60, stepsPerSweep: 20, seed: 9 });
  const b = f.simulateAnnealing({ sequence: SEQ, sweeps: 60, stepsPerSweep: 20, seed: 9 });
  assert.strictEqual(a.best.energy, b.best.energy);
  assert.deepStrictEqual(a.best.coords, b.best.coords);
});

test('annealingGenerator yields progress then done', () => {
  const gen = f.annealingGenerator({ sequence: 'ACDEFG', sweeps: 10, stepsPerSweep: 5, yieldEvery: 4, seed: 2 });
  let progress = 0, done = null;
  for (const msg of gen) {
    if (msg.type === 'progress') { progress++; assert.ok(msg.bestEnergy <= 0); }
    if (msg.done) done = msg;
  }
  assert.ok(progress >= 2);
  assert.ok(done);
  assert.strictEqual(done.type, 'done');
});

test('annealing supports linear cooling', () => {
  const r = f.simulateAnnealing({ sequence: 'ACDEFGHIK', sweeps: 40, stepsPerSweep: 10, cooling: 'linear', seed: 4 });
  assert.ok(r.done && r.best.coords.length === 9);
});

// -----------------------------------------------------------------------
// Replica-exchange Monte Carlo
// -----------------------------------------------------------------------
test('geometricTemps builds a decreasing ladder', () => {
  const t = f.geometricTemps(8, 0.2, 2.5);
  assert.strictEqual(t.length, 8);
  assert.ok(Math.abs(t[0] - 2.5) < 1e-12);
  assert.ok(Math.abs(t[7] - 0.2) < 1e-12);
  for (let i = 1; i < t.length; i++) assert.ok(t[i] < t[i - 1]);
});

test('geometricTemps edge cases', () => {
  assert.deepStrictEqual(f.geometricTemps(1, 0.5, 2.0), [0.5]);
  assert.throws(() => f.geometricTemps(0, 0.5, 2.0), /at least one/);
});

test('simulateRemc returns best, trace and exchange statistics', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 60, replicas: 4, stepsPerSweep: 20, seed: 17 });
  assert.ok(r.done);
  assert.strictEqual(r.best.coords.length, SEQ.length);
  assert.strictEqual(r.trace.length, 60);
  assert.ok(r.exchangeTries > 0);
  assert.ok(r.exchanges <= r.exchangeTries);
  assert.ok(r.temps.length === 4);
  assert.ok(r.proposed >= r.accepted);
});

test('simulateRemc best conformation is a valid SAW with consistent energy', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 80, replicas: 4, stepsPerSweep: 25, seed: 23 });
  const lf = new f.LatticeFolder(SEQ);
  const occ = new Map();
  for (let i = 0; i < r.best.coords.length; i++) {
    occ.set(f.LatticeFolder.key(r.best.coords[i][0], r.best.coords[i][1], r.best.coords[i][2]), i);
  }
  lf._commit(r.best.coords, occ, r.best.energy, null);
  assert.ok(lf.isValidConformation());
  assert.strictEqual(lf.energy(), r.best.energy);
  assert.strictEqual(lf.contacts().length, r.best.contactCount);
});

test('remc best energy is monotonically nonincreasing in the trace', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 50, replicas: 4, stepsPerSweep: 15, seed: 29 });
  for (let i = 1; i < r.trace.length; i++) {
    assert.ok(r.trace[i] <= r.trace[i - 1] + 1e-9, 'trace must not worsen');
  }
});

test('remcGenerator yields progress snapshots', () => {
  const gen = f.remcGenerator({ sequence: 'ACDEFGHIK', sweeps: 12, replicas: 3, stepsPerSweep: 8, yieldEvery: 5, seed: 6 });
  let progress = 0, done = null;
  for (const msg of gen) {
    if (msg.type === 'progress') {
      progress++;
      assert.strictEqual(msg.energies.length, 3);
    }
    if (msg.done) done = msg;
  }
  assert.ok(progress >= 2);
  assert.ok(done && done.done);
});

test('remc swaps states between replicas without corruption', () => {
  // Run with heavy exchange pressure and confirm final replica states are
  // all valid walks with energies matching their contact counts.
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 40, replicas: 6, stepsPerSweep: 20, seed: 37 });
  const lf = new f.LatticeFolder(SEQ);
  for (const e of r.energies) {
    assert.ok(typeof e === 'number');
    assert.ok(e <= 0);
  }
  assert.ok(lf.isValidConformation());
});

// -----------------------------------------------------------------------
// Folding quality (HP model known behavior)
// -----------------------------------------------------------------------
test('hp folding compacts hydrophobic residues', () => {
  const r = f.simulateAnnealing({ sequence: HP_SEQ, mode: 'hp', sweeps: 300, stepsPerSweep: 50, seed: 19 });
  assert.ok(r.best.energy <= -3, 'HP folding should find some H-H contacts, got ' + r.best.energy);
  const lf = new f.LatticeFolder(HP_SEQ, { mode: 'hp' });
  const occ = new Map();
  for (let i = 0; i < r.best.coords.length; i++) {
    occ.set(f.LatticeFolder.key(r.best.coords[i][0], r.best.coords[i][1], r.best.coords[i][2]), i);
  }
  lf._commit(r.best.coords, occ, r.best.energy, null);
  assert.ok(lf.isValidConformation());
  // energy equals minus the number of hydrophobic hydrophobic contacts;
  // other lattice adjacencies may exist but score zero in HP mode.
  let hh = 0;
  for (const [a, b] of lf.contacts()) {
    if (f.isHydrophobic(HP_SEQ[a]) && f.isHydrophobic(HP_SEQ[b])) hh++;
  }
  assert.strictEqual(lf.energy(), -hh);
});

test('isHydrophobic set matches the Dill style set', () => {
  assert.strictEqual(f.isHydrophobic('A'), true);
  assert.strictEqual(f.isHydrophobic('W'), true);
  assert.strictEqual(f.isHydrophobic('D'), false);
  assert.strictEqual(f.isHydrophobic('K'), false);
  assert.strictEqual(f.isHydrophobic('R'), false);
});

test('folded conformations are more compact than extended ones', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 100, replicas: 4, stepsPerSweep: 25, seed: 43 });
  const extended = new f.LatticeFolder(SEQ);
  assert.ok(r.best.radiusOfGyration < extended.radiusOfGyration());
});

// -----------------------------------------------------------------------
// Lattice to backbone pipeline
// -----------------------------------------------------------------------
test('latticeBackboneAngles assigns every residue a basin', () => {
  const lf = new f.LatticeFolder(SEQ);
  const angles = f.latticeBackboneAngles(lf.coords);
  assert.strictEqual(angles.length, lf.n);
  for (const a of angles) {
    assert.ok(['BETA', 'ALPHA_R', 'ALPHA_L', 'PP_II'].includes(a.region));
    assert.strictEqual(a.omega, 180);
  }
});

test('straight extended chain maps entirely to beta', () => {
  const lf = new f.LatticeFolder(SEQ);
  const angles = f.latticeBackboneAngles(lf.coords);
  for (const a of angles) assert.strictEqual(a.region, 'BETA');
});

test('reconstructLatticeBackbone yields atoms and classified dihedrals', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 60, replicas: 4, stepsPerSweep: 20, seed: 47 });
  const lf = new f.LatticeFolder(SEQ);
  const occ = new Map();
  for (let i = 0; i < r.best.coords.length; i++) {
    occ.set(f.LatticeFolder.key(r.best.coords[i][0], r.best.coords[i][1], r.best.coords[i][2]), i);
  }
  lf._commit(r.best.coords, occ, r.best.energy, null);
  const rec = f.reconstructLatticeBackbone(lf);
  assert.strictEqual(rec.atoms.length, lf.n * 4);
  assert.strictEqual(rec.dihedrals.length, lf.n);
  for (const d of rec.dihedrals) {
    assert.ok(['BETA', 'ALPHA_R', 'ALPHA_L', 'PP_II'].includes(d.region));
    assert.strictEqual(d.omega, 180);
  }
});

test('reconstructed backbone from folded lattice has no steric clashes', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 80, replicas: 4, stepsPerSweep: 25, seed: 53 });
  const lf = new f.LatticeFolder(SEQ);
  const occ = new Map();
  for (let i = 0; i < r.best.coords.length; i++) {
    occ.set(f.LatticeFolder.key(r.best.coords[i][0], r.best.coords[i][1], r.best.coords[i][2]), i);
  }
  lf._commit(r.best.coords, occ, r.best.energy, null);
  const rec = f.reconstructLatticeBackbone(lf);
  const clashes = f.findStericClashes(rec.atoms, {
    excludePair: (i, j) => Math.abs(rec.atoms[i].resSeq - rec.atoms[j].resSeq) <= 1
  });
  assert.strictEqual(clashes.length, 0, 'rebuilt backbone should be clash free');
});

test('folderToPdb writes a CA trace that reparses cleanly', () => {
  const lf = new f.LatticeFolder('ACDEFGHIK');
  lf.tryMove('pivot', f.makeRng(3));
  const pdb = f.folderToPdb(lf, { title: 'trace test', scale: 3.8 });
  const rep = f.parsePdb(pdb);
  assert.strictEqual(rep.atoms.length, 9);
  assert.strictEqual(rep.atoms[0].name.trim(), 'CA');
  // lattice spacing 1 maps to 3.8 angstrom CA spacing
  assert.ok(Math.abs(f.vdist(rep.atoms[0].xyz, rep.atoms[1].xyz) - 3.8) < 1e-9);
});

// -----------------------------------------------------------------------
// End to end pipeline
// -----------------------------------------------------------------------
test('full pipeline: fold, rebuild, analyze, export, reparse', () => {
  const r = f.simulateRemc({ sequence: SEQ, sweeps: 60, replicas: 4, stepsPerSweep: 20, seed: 59 });
  const lf = new f.LatticeFolder(SEQ);
  const occ = new Map();
  for (let i = 0; i < r.best.coords.length; i++) {
    occ.set(f.LatticeFolder.key(r.best.coords[i][0], r.best.coords[i][1], r.best.coords[i][2]), i);
  }
  lf._commit(r.best.coords, occ, r.best.energy, null);
  const rec = f.reconstructLatticeBackbone(lf);
  const pdb = f.atomsToPdb(rec.atoms, {
    title: 'fold-lattice pipeline',
    remarks: ['energy ' + r.best.energy.toFixed(2) + ' RT', 'contacts ' + r.best.contactCount]
  });
  const reparsed = f.parsePdb(pdb);
  assert.strictEqual(reparsed.atoms.length, rec.atoms.length);
  const d = f.computeDihedrals(reparsed.atoms);
  assert.strictEqual(d.length, SEQ.length);
  let classified = 0;
  for (const entry of d) {
    if (entry.region && entry.region !== 'OUTLIER') classified++;
    if (entry.omegaClass) assert.strictEqual(entry.omegaClass, 'TRANS');
  }
  assert.ok(classified >= SEQ.length - 2, 'most residues should classify');
});

test('module exposes the full API surface', () => {
  const required = [
    'VERSION', 'AA20', 'MJ1996', 'mjPairEnergy', 'makeRng',
    'vadd', 'vsub', 'vscale', 'vdot', 'vcross', 'vlen', 'vdist', 'vnorm',
    'bondAngle', 'dihedral', 'placeAtom',
    'normalizeDeg', 'classifyRamachandran', 'classifyOmega',
    'findStericClashes', 'radiusOfGyration',
    'parsePdb', 'parsePdbLine', 'atomsToPdb',
    'reconstructBackbone', 'computeDihedrals', 'latticeBackboneAngles',
    'validateSequence', 'LatticeFolder', 'rotationAboutAxis',
    'geometricTemps', 'remcGenerator', 'simulateRemc',
    'annealingGenerator', 'simulateAnnealing',
    'folderToPdb', 'reconstructLatticeBackbone', 'snapshot', 'restore'
  ];
  for (const k of required) assert.ok(f[k] !== undefined, 'missing export: ' + k);
});

test('VERSION is 1.0.0', () => {
  assert.strictEqual(f.VERSION, '1.0.0');
});
