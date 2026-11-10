'use strict';

/*
 * Terminal demo for fold-lattice. Folds a few sequences with replica
 * exchange, prints the results, rebuilds an all-atom backbone and prints
 * the first lines of a PDB export.
 *
 * Run: node demo/commands.js
 */

const foldl = require('../src/fold.js');

const SEQUENCES = [
  { name: 'demo 24mer', seq: 'GAGSLYQGSSLAQSGFLSSHHSGF', mode: 'mj' },
  { name: 'block HP 21mer', seq: 'SSSSAAAASSSSAAAASSSAA', mode: 'hp' },
  { name: 'tryptophan zipper like', seq: 'SWTWENGKWTWK', mode: 'mj' }
];

function bar(value, min, max, width) {
  const frac = Math.max(0, Math.min(1, (value - min) / (max - min)));
  const filled = Math.round(frac * width);
  return '[' + '#'.repeat(filled) + '.'.repeat(width - filled) + ']';
}

console.log('fold-lattice ' + foldl.VERSION + ' terminal demo');
console.log('nodes ' + process.version + ', zero dependencies');
console.log('');

for (const job of SEQUENCES) {
  console.log('--- ' + job.name + ' (' + job.seq.length + ' residues, ' + job.mode.toUpperCase() + ' model) ---');
  const t0 = Date.now();
  const result = foldl.simulateRemc({
    sequence: job.seq,
    mode: job.mode,
    replicas: 6,
    sweeps: 150,
    stepsPerSweep: 30,
    seed: 11
  });
  const ms = Date.now() - t0;
  const best = result.best;
  console.log('  best energy        ' + best.energy.toFixed(2) + ' RT');
  console.log('  contacts           ' + best.contactCount);
  console.log('  radius of gyration ' + best.radiusOfGyration.toFixed(2) + ' lattice units');
  console.log('  replica exchanges  ' + result.exchanges + ' of ' + result.exchangeTries + ' accepted');
  console.log('  mc acceptance      ' + (100 * result.accepted / result.proposed).toFixed(1) + '%');
  console.log('  wall time          ' + ms + ' ms for ' + (result.sweeps * 6) + ' sweeps');
  console.log('  energy trace       ' + bar(best.energy, -40, 0, 30) + ' ' + result.trace[0].toFixed(1) + ' to ' + best.energy.toFixed(1));

  // Rebuild an all-atom backbone and classify Ramachandran regions.
  const folder = new foldl.LatticeFolder(job.seq, { mode: job.mode });
  const occ = new Map();
  for (let i = 0; i < best.coords.length; i++) {
    occ.set(foldl.LatticeFolder.key(best.coords[i][0], best.coords[i][1], best.coords[i][2]), i);
  }
  folder._commit(best.coords, occ, best.energy, null);
  const rec = foldl.reconstructLatticeBackbone(folder);
  const counts = {};
  for (const d of rec.dihedrals) counts[d.region] = (counts[d.region] || 0) + 1;
  const clashes = foldl.findStericClashes(rec.atoms, {
    excludePair: (i, j) => Math.abs(rec.atoms[i].resSeq - rec.atoms[j].resSeq) <= 1
  });
  console.log('  rebuilt atoms      ' + rec.atoms.length + ' (N, CA, C, O per residue)');
  console.log('  ramachandran       ' + JSON.stringify(counts));
  console.log('  steric clashes     ' + clashes.length);
  console.log('');
}

// PDB export sample from the first job.
const job = SEQUENCES[0];
const result = foldl.simulateRemc({ sequence: job.seq, replicas: 6, sweeps: 150, stepsPerSweep: 30, seed: 11 });
const folder = new foldl.LatticeFolder(job.seq);
const occ = new Map();
for (let i = 0; i < result.best.coords.length; i++) {
  occ.set(foldl.LatticeFolder.key(result.best.coords[i][0], result.best.coords[i][1], result.best.coords[i][2]), i);
}
folder._commit(result.best.coords, occ, result.best.energy, null);
const rec = foldl.reconstructLatticeBackbone(folder);
const pdb = foldl.atomsToPdb(rec.atoms, {
  title: 'fold-lattice demo ' + job.name,
  remarks: ['energy ' + result.best.energy.toFixed(2) + ' RT']
});
console.log('--- pdb export sample (first 8 lines of ' + pdb.split('\n').length + ') ---');
console.log(pdb.split('\n').slice(0, 8).join('\n'));
console.log('');
console.log('done. open index.html for the interactive 3D viewer.');
