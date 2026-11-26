# fold-lattice

A zero dependency toolkit for coarse grained protein folding, replica-exchange
Monte Carlo, Ramachandran analysis and PDB export. Pure JavaScript, runs in
Node.js and in the browser with no packages to install.

Version 2 adds the HP workbench: a Dill-model algorithm comparator on 2D
square and 3D cubic lattices, exact ground state enumeration for short
chains, a folding funnel landscape, mutation stress scans, an LxL residue
contact map, SHA-256 checksummed audit exports, URL hash bookmarks and an
in page self test runner.

## What it does

**Lattice folding engine.** Protein chains live on the 3D cubic lattice as
self avoiding walks. Energy is scored with either the Dill HP model (H-H
contacts score -1) or the Miyazawa-Jernigan 1996 20x20 contact potential in
RT units. Moves are pull moves (after Lesh, Mitzenmacher and Whitesides),
pivot rotations, corner flips, two-bead crankshafts and endpoint pulls.

**Search algorithms.** Metropolis Monte Carlo simulated annealing with
geometric or linear cooling, and parallel tempering replica-exchange Monte
Carlo (REMC) with a geometric temperature ladder and pairwise Metropolis
exchanges between adjacent replicas.

**All-atom backbone rebuild.** NeRF (Natural Extension Reference Frame)
placement reconstructs N, CA, C, O atoms from per-residue (phi, psi, omega)
dihedrals with Engh and Huber bond geometry. The rebuild reproduces the
input torsions to machine precision and produces clash-free structures for
lattice conformations.

**Ramachandran classification.** Every residue is classified into the alpha
right, alpha left, beta, polyproline II, gamma turn or outlier regions of
the Ramachandran map, plus cis and trans peptide bond checks via omega.

**Steric analysis.** Radius of gyration and a van der Waals steric clash
detector with element-aware radii, bonded pair exclusion and custom pair
filters.

**PDB input and output.** Parse ATOM and HETATM records (first model,
alternate locations A or blank), and export atoms to standard PDB text with
correct column layout, TER records and element columns.

**HP workbench (`src/lab.js` + `lab-app.js`).** Five seeded algorithms run
side by side on the same HP sequence: random SAW sampling, greedy local
pivot, fixed temperature Metropolis Monte Carlo, simulated annealing with a
geometric schedule, and exact exhaustive enumeration over symmetry reduced
self avoiding walks (up to length 16 in 2D and 10 in 3D). Every algorithm
reports best energy, contact count, H-H contacts, radius of gyration,
end to end distance, wall time and the optimality gap against the exact
ground state. The workbench also draws the energy vs compactness funnel
cloud, an interactive LxL residue contact matrix, a per residue mutation
stability scan (H to P and P to H delta-E), and exports a canonicalized
JSON / CSV / Markdown audit report stamped with a pure JavaScript FIPS
180-4 SHA-256 digest (unkeyed tamper detection, not a signature). Runs
bookmark into the URL hash as `#fold=SEQ:DIMd:SEED` and any compared
conformation can be shown on the original 3D stage.

> The Dill HP lattice model is a simplified statistical physics toy. It is
> not an atomistic force field and not an AlphaFold style structure
> predictor; it demonstrates folding thermodynamics and search dynamics on
> a lattice.

## Quick start

```bash
npm test                        # run the 447 unit tests
node demo/commands.js           # run a folding demo in the terminal
node tools/browser-audit.js     # headless Chromium UI audit, needs Chrome and Node 22+
```

Or open `index.html` for the interactive 3D viewer. Nothing to build, no
dependencies to install.

### Fold a protein in 10 lines

```js
const foldl = require('./src/fold.js');

const result = foldl.simulateRemc({
  sequence: 'GAGSLYQGSSLAQSGFLSSHHSGF',
  mode: 'mj',        // Miyazawa-Jernigan energies (or 'hp')
  replicas: 6,       // temperature ladder size
  sweeps: 200,
  stepsPerSweep: 30,
  seed: 11
});

console.log('best energy', result.best.energy);          // RT units
console.log('contacts', result.best.contactCount);
console.log('radius of gyration', result.best.radiusOfGyration);
```

### Rebuild an all-atom backbone and export PDB

```js
const foldl = require('./src/fold.js');

const folder = new foldl.LatticeFolder('GAGSLYQGSSLAQSGFLSSHHSGF');
// ... run moves, or seed from result.best.coords ...
const rec = foldl.reconstructLatticeBackbone(folder);
const pdb = foldl.atomsToPdb(rec.atoms, { title: 'my fold' });

for (const d of rec.dihedrals) {
  console.log(d.resIndex, d.region, 'phi', d.phi, 'psi', d.psi);
}
```

## The Miyazawa-Jernigan potential

The 20x20 contact energy matrix embedded in `src/fold.js` was transcribed
from Table 3 of:

> Miyazawa, S. and Jernigan, R.L. (1996). Residue-residue potentials with a
> favorable contact pair term and an unfavorable high packing density term,
> for simulation and threading. Journal of Molecular Biology 256, 623-644.

Row and column order is C, M, F, I, L, V, W, Y, A, G, T, S, N, Q, D, E, H,
R, K, P. Energies are in RT units. Contacts are counted between residues
that are lattice adjacent but not sequential neighbors in the chain. The
strongest diagonal entry is L-L at -7.37 and the strongest off diagonal
pair is F-L at -7.28, matching the published values.

## Move sets

- **Pull moves**: the complete, reversible local move of Lesh, Mitzenmacher
  and Whitesides, generalized to the cubic lattice. A residue slides into a
  free cell adjacent to its chain neighbor, the corner of the unit square is
  checked, and the chain is pulled two residues at a time until connected.
- **Pivot moves**: rotate the tail beyond a bond by one to three quarter
  turns about the bond axis, using exact integer rotation matrices so every
  result stays on the lattice.
- **Corner flips**: an interior residue jumps to the fourth corner of the
  unit square spanned by its neighbors.
- **Crankshafts**: the two middle residues of a square U turn rotate 180
  degrees about the axis through the ends.
- **Endpoint pulls**: terminal residues walk out into a free path of two
  cells and the chain follows.

## Replica exchange

Parallel tempering runs N replicas at geometrically spaced temperatures
between `minTemp` and `maxTemp` (in RT units). Each sweep performs a fixed
number of Metropolis steps per replica, then attempts pairwise exchanges
between adjacent temperature slots with the criterion
`exp((1/Ta - 1/Tb)(Ea - Eb))`, alternating even and odd pairs. The best
conformation ever seen by any replica is tracked throughout.

## Tests

447 unit tests cover the RNG, vector math, the MJ matrix (symmetry, anchor
values, strongest pairs), sequence validation, bond angles and dihedrals,
NeRF placement round trips, rotation matrices, Ramachandran regions, steric
clashes, radius of gyration, PDB parsing and export (column layout, TER
records, chain breaks, alternate locations, multi-model files), every move
type, occupancy-map integrity, Metropolis edge cases, option validation,
simulated annealing, REMC exchange statistics, HP model behavior, the
lattice to backbone pipeline, viewer security and accessibility markup, and
a full fold to rebuild to analyze to export to reparse end to end run. The
v2 lab suite adds 226 more tests: canonical SAW counts, self avoidance and
unit step invariants under every lab move, energy translation, rotation,
reversal and bound properties, exact enumeration against an independent
brute force DFS, hand verified ground states, comparator determinism and
optimality gaps, mutation delta-E recomputation, contact matrix symmetry
and coding, SHA-256 against FIPS 180-4 vectors plus a platform crypto
oracle at every padding boundary, canonical JSON key ordering, checksum
tamper detection, CSV formula injection escaping, bookmark edge cases, and
static checks on the shipped files (strict mode, no innerHTML sinks, no
external loads, accessibility roles, CSP intact).

`tools/browser-audit.js` additionally drives real headless Chromium over
the DevTools protocol at 1440x900 and 390x844 (42 checks): page load, zero
console and page errors, zero external requests, zero horizontal overflow,
the comparator table, funnel and contact map SVGs, the mutation scan,
checksummed audit verify and tamper flows, the URL hash bookmark, the self
test runner and the view on stage bridge.

## Browser viewer

`index.html` is a self contained interactive 3D molecular viewer with a
hand written canvas renderer. It folds chains in your browser with a worker
driven progress bar, lets you scrub through the replica trace, inspect each
residue's Ramachandran class, switch between HP and MJ potentials, and
export PDB files from folded conformations. Below it, the HP workbench
runs the algorithm comparator, folding funnel, contact map, mutation scan,
checksummed audit exports, bookmarking and self tests described above, all
offline with zero network requests.

## License

MIT
