/*
 * fold-lattice interactive viewer. Hand written canvas 3D renderer, worker
 * driven replica-exchange folding, Ramachandran map, PDB import and export.
 * No dependencies, no build step.
 */
(function () {
  'use strict';

  var F = window.foldl;

  // ---------------------------------------------------------------------
  // dom helpers
  // ---------------------------------------------------------------------
  function $(id) { return document.getElementById(id); }

  var stage = $('stage');
  var legend = $('legend');
  var ribbon = $('ribbon');
  var residueCard = $('residueCard');
  var resBody = $('resBody');
  var resClose = $('resClose');
  var seqInput = $('seqInput');
  var presetSel = $('presetSel');
  var modelSel = $('modelSel');
  var replicasR = $('replicasR');
  var replicasV = $('replicasV');
  var sweepsR = $('sweepsR');
  var sweepsV = $('sweepsV');
  var colorSel = $('colorSel');
  var showAll = $('showAll');
  var foldBtn = $('foldBtn');
  var stopBtn = $('stopBtn');
  var pdbBtn = $('pdbBtn');
  var errMsg = $('errMsg');
  var progressBar = $('progressBar');
  var progressText = $('progressText');
  var stEnergy = $('stEnergy');
  var stContacts = $('stContacts');
  var stRg = $('stRg');
  var stExch = $('stExch');
  var stAcc = $('stAcc');
  var stSweeps = $('stSweeps');
  var pdbFile = $('pdbFile');
  var pdbMsg = $('pdbMsg');
  var traceChart = $('traceChart');
  var ramaChart = $('ramaChart');
  var verTag = $('verTag');

  verTag.textContent = 'v' + F.VERSION + ', ' + F.RAMACHANDRAN_REGIONS.length +
    ' Ramachandran regions, 20x20 MJ matrix';

  // ---------------------------------------------------------------------
  // state
  // ---------------------------------------------------------------------
  var state = {
    source: 'fold',          // 'fold' or 'pdb'
    view: 'lattice',         // 'lattice' or 'atoms'
    points: [],              // { p: [x,y,z], r: worldRadius, resIndex, rgb, element }
    bonds: [],               // [i, j] index pairs into points
    residues: [],            // { letter, num, chain, resName }
    dihedrals: null,
    rebuilt: null,           // cached all atom rebuild of the final snapshot
    snapshots: [],           // { sweep, energy, coords, exchanges, accepted, proposed }
    scrub: 0,
    trace: null,
    result: null,
    selected: null,
    running: false,
    userScrubbed: false,
    seq: null,
    mode: 'mj',
    pdbAtoms: null
  };

  var cam = { yaw: 0.65, pitch: -0.35, zoom: 1 };
  var sceneCenter = [0, 0, 0];
  var sceneR = 10;
  var lastPX = [], lastPY = [], lastPR = [];

  // ---------------------------------------------------------------------
  // colors
  // ---------------------------------------------------------------------
  var REGION_COLORS = {};
  for (var ri = 0; ri < F.RAMACHANDRAN_REGIONS.length; ri++) {
    REGION_COLORS[F.RAMACHANDRAN_REGIONS[ri].name] = F.RAMACHANDRAN_REGIONS[ri].color;
  }
  var OUTLIER_COLOR = '#8fa0bd';
  var HYDRO_COLOR = '#e0913e';
  var POLAR_COLOR = '#4ecdc4';
  var ELEMENT_COLORS = { C: '#9aa5b1', N: '#4e7de0', O: '#e0564e', S: '#e0b13e', H: '#d8dee9', P: '#e0913e' };

  function hexToRgb(hex) {
    var h = hex.replace('#', '');
    return [
      parseInt(h.slice(0, 2), 16),
      parseInt(h.slice(2, 4), 16),
      parseInt(h.slice(4, 6), 16)
    ];
  }

  function hslHex(h, s, l) {
    h = (((h % 360) + 360) % 360) / 360;
    var a = s * Math.min(l, 1 - l);
    function f(n) {
      var k = (n + h * 12) % 12;
      var c = l - a * Math.max(-1, Math.min(k - 3, Math.min(9 - k, 1)));
      return Math.round(255 * c).toString(16).padStart(2, '0');
    }
    return '#' + f(0) + f(8) + f(4);
  }

  function regionColor(name) {
    return REGION_COLORS[name] || OUTLIER_COLOR;
  }

  // ---------------------------------------------------------------------
  // geometry builders
  // ---------------------------------------------------------------------
  function setSceneBounds(pts) {
    var n = pts.length;
    if (!n) { sceneR = 10; sceneCenter = [0, 0, 0]; return; }
    var c = [0, 0, 0];
    for (var i = 0; i < n; i++) {
      c[0] += pts[i].p[0]; c[1] += pts[i].p[1]; c[2] += pts[i].p[2];
    }
    c[0] /= n; c[1] /= n; c[2] /= n;
    var r = 0;
    for (var j = 0; j < n; j++) {
      var dx = pts[j].p[0] - c[0], dy = pts[j].p[1] - c[1], dz = pts[j].p[2] - c[2];
      var d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (d > r) { r = d; }
    }
    sceneCenter = c;
    sceneR = Math.max(r, 1.5);
  }

  function groupResidues(atoms) {
    var groups = [];
    var cur = null;
    for (var i = 0; i < atoms.length; i++) {
      var a = atoms[i];
      if (!cur || cur.chainID !== a.chainID || cur.resSeq !== a.resSeq) {
        cur = { chainID: a.chainID, resSeq: a.resSeq, resName: a.resName, start: i, end: i };
        groups.push(cur);
      }
      cur.end = i;
    }
    return groups;
  }

  function resIndexOfGroup(groups, atomIndex) {
    for (var g = 0; g < groups.length; g++) {
      if (atomIndex >= groups[g].start && atomIndex <= groups[g].end) { return g; }
    }
    return groups.length - 1;
  }

  function setFromLattice(coords) {
    var pts = [];
    var bonds = [];
    for (var i = 0; i < coords.length; i++) {
      pts.push({ p: [coords[i][0], coords[i][1], coords[i][2]], r: 0.30, resIndex: i, element: null });
      if (i > 0) { bonds.push([i - 1, i]); }
    }
    var residues = [];
    var seqStr = state.seq || '';
    for (var k = 0; k < coords.length; k++) {
      var letter = seqStr.charAt(k) || '?';
      residues.push({ letter: letter, num: k + 1, chain: 'A', resName: F.AA1_TO_3[letter] || '?' });
    }
    state.view = 'lattice';
    state.points = pts;
    state.bonds = bonds;
    state.residues = residues;
    setSceneBounds(pts);
  }

  function setFromAtoms(atoms) {
    var groups = groupResidues(atoms);
    var pts = [];
    var bonds = [];
    for (var i = 0; i < atoms.length; i++) {
      var a = atoms[i];
      pts.push({
        p: [a.xyz[0], a.xyz[1], a.xyz[2]],
        r: Math.max(F.atomRadius(a) * 0.45, 0.5),
        resIndex: resIndexOfGroup(groups, i),
        element: (a.element || 'C').toUpperCase()
      });
    }
    for (var j = 0; j < atoms.length; j++) {
      for (var k2 = j + 1; k2 < atoms.length; k2++) {
        var b = atoms[k2];
        var a2 = atoms[j];
        if (b.chainID !== a2.chainID) { break; }
        if (Math.abs(b.resSeq - a2.resSeq) > 1) { break; }
        var dx = b.xyz[0] - a2.xyz[0], dy = b.xyz[1] - a2.xyz[1], dz = b.xyz[2] - a2.xyz[2];
        var d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > 0 && d2 < 4.0) { bonds.push([j, k2]); }
      }
    }
    var residues = [];
    for (var g = 0; g < groups.length; g++) {
      residues.push({
        letter: F.AA3_TO_1[groups[g].resName] || '?',
        num: groups[g].resSeq,
        chain: groups[g].chainID,
        resName: groups[g].resName
      });
    }
    state.view = 'atoms';
    state.points = pts;
    state.bonds = bonds;
    state.residues = residues;
    setSceneBounds(pts);
  }

  // ---------------------------------------------------------------------
  // colors applied to the current point set
  // ---------------------------------------------------------------------
  function applyColors() {
    var cm = colorSel.value;
    var n = state.residues.length;
    var residueRgb = [];
    for (var i = 0; i < n; i++) {
      var hex = null;
      if (cm === 'hydrophobic') {
        hex = F.isHydrophobic(state.residues[i].letter) ? HYDRO_COLOR : POLAR_COLOR;
      } else if (cm === 'ramachandran') {
        var d = state.dihedrals ? state.dihedrals[i] : null;
        hex = d && d.region ? regionColor(d.region) : OUTLIER_COLOR;
      } else if (cm === 'element') {
        hex = null; // resolved per atom below
      } else {
        hex = n > 1 ? hslHex(300 * i / (n - 1), 0.72, 0.62) : hslHex(180, 0.72, 0.62);
      }
      residueRgb.push(hex ? hexToRgb(hex) : null);
    }
    for (var p = 0; p < state.points.length; p++) {
      var pt = state.points[p];
      var rgb;
      if (cm === 'element' && state.view === 'atoms') {
        rgb = hexToRgb(ELEMENT_COLORS[pt.element] || ELEMENT_COLORS.C);
      } else if (cm === 'element' && state.view === 'lattice' && state.seq) {
        rgb = hexToRgb(F.isHydrophobic(state.seq.charAt(pt.resIndex)) ? HYDRO_COLOR : POLAR_COLOR);
      } else {
        rgb = residueRgb[pt.resIndex] || hexToRgb(OUTLIER_COLOR);
      }
      pt.rgb = rgb;
    }
    buildLegend();
    buildRibbon();
  }

  // ---------------------------------------------------------------------
  // renderer
  // ---------------------------------------------------------------------
  function fitCanvas(cv) {
    var dpr = window.devicePixelRatio || 1;
    var w = cv.clientWidth, h = cv.clientHeight;
    var W = Math.round(w * dpr), H = Math.round(h * dpr);
    if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
    return { w: w, h: h, dpr: dpr };
  }

  function rotatePoint(p) {
    var x = p[0] - sceneCenter[0];
    var y = p[1] - sceneCenter[1];
    var z = p[2] - sceneCenter[2];
    var cy = Math.cos(cam.yaw), sy = Math.sin(cam.yaw);
    var x1 = cy * x + sy * z;
    var z1 = -sy * x + cy * z;
    var cp = Math.cos(cam.pitch), sp = Math.sin(cam.pitch);
    var y2 = cp * y - sp * z1;
    var z2 = sp * y + cp * z1;
    return [x1, y2, z2];
  }

  function drawScene() {
    var m = fitCanvas(stage);
    var ctx = stage.getContext('2d');
    ctx.setTransform(m.dpr, 0, 0, m.dpr, 0, 0);
    ctx.clearRect(0, 0, m.w, m.h);
    var pts = state.points;
    lastPX = []; lastPY = []; lastPR = [];
    if (!pts.length) { return; }

    var cx = m.w / 2, cy = m.h / 2;
    var D = Math.max(sceneR * 3, 10);
    var focal = 0.42 * Math.min(m.w, m.h) * (D / sceneR) * cam.zoom;

    var px = [], py = [], pd = [], ps = [];
    for (var i = 0; i < pts.length; i++) {
      var rp = rotatePoint(pts[i].p);
      var depth = D + rp[2];
      if (depth < 0.5) { depth = 0.5; }
      var s = focal / depth;
      px.push(cx + rp[0] * s);
      py.push(cy - rp[1] * s);
      pd.push(depth);
      ps.push(s);
    }

    // bonds first, alpha and width by depth
    for (var b = 0; b < state.bonds.length; b++) {
      var i1 = state.bonds[b][0], i2 = state.bonds[b][1];
      var t = 1 - ((pd[i1] + pd[i2]) / 2 - (D - sceneR)) / (2 * sceneR);
      if (t < 0) { t = 0; } else if (t > 1) { t = 1; }
      var c1 = pts[i1].rgb, c2 = pts[i2].rgb;
      var r = Math.round((c1[0] + c2[0]) / 2);
      var g = Math.round((c1[1] + c2[1]) / 2);
      var bl = Math.round((c1[2] + c2[2]) / 2);
      ctx.strokeStyle = 'rgba(' + r + ',' + g + ',' + bl + ',' + (0.30 + 0.65 * t).toFixed(3) + ')';
      ctx.lineWidth = 0.8 + 1.6 * t;
      ctx.beginPath();
      ctx.moveTo(px[i1], py[i1]);
      ctx.lineTo(px[i2], py[i2]);
      ctx.stroke();
    }

    // beads far to near
    var order = [];
    for (var o = 0; o < pts.length; o++) { order.push(o); }
    order.sort(function (a, b2) { return pd[b2] - pd[a]; });
    for (var q = 0; q < order.length; q++) {
      var idx = order[q];
      var t2 = 1 - (pd[idx] - (D - sceneR)) / (2 * sceneR);
      if (t2 < 0) { t2 = 0; } else if (t2 > 1) { t2 = 1; }
      var rpx = Math.max(pts[idx].r * ps[idx], 1.4);
      var c = pts[idx].rgb;
      ctx.fillStyle = 'rgba(' + c[0] + ',' + c[1] + ',' + c[2] + ',' + (0.35 + 0.65 * t2).toFixed(3) + ')';
      ctx.beginPath();
      ctx.arc(px[idx], py[idx], rpx, 0, 6.2832);
      ctx.fill();
      if (state.selected !== null && pts[idx].resIndex === state.selected) {
        ctx.strokeStyle = 'rgba(255,255,255,0.9)';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.arc(px[idx], py[idx], rpx + 2.5, 0, 6.2832);
        ctx.stroke();
      }
      lastPX[idx] = px[idx]; lastPY[idx] = py[idx]; lastPR[idx] = rpx;
    }

    if (state.view === 'lattice') {
      ctx.fillStyle = 'rgba(223,230,243,0.85)';
      ctx.font = '11px monospace';
      ctx.fillText('N', px[0] - 4, py[0] - (pts[0].r * ps[0] + 6));
      var lastIdx = pts.length - 1;
      ctx.fillText('C', px[lastIdx] - 4, py[lastIdx] - (pts[lastIdx].r * ps[lastIdx] + 6));
    }
  }

  // ---------------------------------------------------------------------
  // camera interaction
  // ---------------------------------------------------------------------
  var drag = null;
  stage.addEventListener('pointerdown', function (e) {
    drag = { x: e.clientX, y: e.clientY, moved: 0 };
    stage.setPointerCapture(e.pointerId);
    stage.classList.add('dragging');
  });
  stage.addEventListener('pointermove', function (e) {
    if (!drag) { return; }
    var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    drag.moved += Math.abs(dx) + Math.abs(dy);
    cam.yaw += dx * 0.008;
    cam.pitch = Math.max(-1.45, Math.min(1.45, cam.pitch + dy * 0.008));
    drag.x = e.clientX; drag.y = e.clientY;
  });
  stage.addEventListener('pointerup', function (e) {
    stage.classList.remove('dragging');
    if (drag && drag.moved < 5) { pickBead(e); }
    drag = null;
  });
  stage.addEventListener('pointercancel', function () { drag = null; stage.classList.remove('dragging'); });
  stage.addEventListener('wheel', function (e) {
    e.preventDefault();
    cam.zoom = Math.max(0.25, Math.min(6, cam.zoom * Math.exp(-e.deltaY * 0.0012)));
  }, { passive: false });

  function pickBead(e) {
    var rect = stage.getBoundingClientRect();
    var mx = e.clientX - rect.left, my = e.clientY - rect.top;
    var bestIdx = -1, bestD = 18;
    for (var i = 0; i < state.points.length; i++) {
      if (lastPX[i] === undefined) { continue; }
      var dx = lastPX[i] - mx, dy = lastPY[i] - my;
      var d = Math.sqrt(dx * dx + dy * dy);
      var reach = Math.max(10, lastPR[i] + 4);
      if (d < reach && d < bestD) { bestIdx = i; bestD = d; }
    }
    if (bestIdx >= 0) { selectResidue(state.points[bestIdx].resIndex); }
    else { selectResidue(null); }
  }

  function selectResidue(idx) {
    state.selected = idx;
    var chips = ribbon.querySelectorAll('.res');
    for (var i = 0; i < chips.length; i++) {
      chips[i].classList.toggle('selected', parseInt(chips[i].getAttribute('data-i'), 10) === idx);
    }
    if (idx === null) { residueCard.style.display = 'none'; return; }
    var res = state.residues[idx];
    if (!res) { residueCard.style.display = 'none'; return; }
    var html = '<b>' + res.resName + ' ' + res.num + '</b> chain ' + res.chain;
    var d = state.dihedrals ? state.dihedrals[idx] : null;
    if (d && d.phi !== null && d.phi !== undefined) {
      html += '<br>region: ' + (d.region || 'outlier');
      html += '<br>phi ' + fmtDeg(d.phi) + ', psi ' + fmtDeg(d.psi);
      if (d.omega !== null && d.omega !== undefined) {
        html += '<br>omega ' + fmtDeg(d.omega) + ' (' + (d.omegaClass || '?') + ')';
      }
    } else if (state.view === 'lattice') {
      var pt = null;
      for (var i2 = 0; i2 < state.points.length; i2++) {
        if (state.points[i2].resIndex === idx) { pt = state.points[i2]; break; }
      }
      if (pt) {
        html += '<br>cell (' + pt.p[0] + ', ' + pt.p[1] + ', ' + pt.p[2] + ')';
        html += '<br>' + (F.isHydrophobic(res.letter) ? 'hydrophobic' : 'polar') + ' residue';
        html += '<br>switch the display to all atom for angles';
      }
    }
    resBody.innerHTML = html;
    residueCard.style.display = 'block';
  }

  function fmtDeg(v) {
    return (Math.round(v * 10) / 10).toFixed(1);
  }

  resClose.addEventListener('click', function () { selectResidue(null); });

  // ---------------------------------------------------------------------
  // legend and ribbon
  // ---------------------------------------------------------------------
  function sw(hex, label) {
    return '<div class="row"><span class="swatch" style="background:' + hex + '"></span>' + label + '</div>';
  }

  function buildLegend() {
    var cm = colorSel.value;
    var html = '';
    if (!state.points.length) { legend.style.display = 'none'; return; }
    if (cm === 'sequence') {
      html = '<div class="row"><span class="swatch" style="width:60px;background:linear-gradient(90deg,' +
        hslHex(0, 0.72, 0.62) + ',' + hslHex(75, 0.72, 0.62) + ',' + hslHex(150, 0.72, 0.62) + ',' +
        hslHex(225, 0.72, 0.62) + ',' + hslHex(300, 0.72, 0.62) + ')"></span>N terminus to C terminus</div>';
    } else if (cm === 'hydrophobic') {
      html = sw(HYDRO_COLOR, 'hydrophobic (ACFILMVWY)') + sw(POLAR_COLOR, 'polar');
    } else if (cm === 'ramachandran') {
      html = sw(regionColor('ALPHA_R'), 'alpha right') + sw(regionColor('ALPHA_L'), 'alpha left') +
        sw(regionColor('PP_II'), 'polyproline II') + sw(regionColor('BETA'), 'beta sheet') +
        sw(regionColor('EPSILON'), 'beta bridge') + sw(regionColor('GAMMA'), 'gamma turn') +
        sw(OUTLIER_COLOR, 'outlier');
    } else if (cm === 'element' && state.view === 'atoms') {
      html = sw(ELEMENT_COLORS.C, 'carbon') + sw(ELEMENT_COLORS.N, 'nitrogen') +
        sw(ELEMENT_COLORS.O, 'oxygen') + sw(ELEMENT_COLORS.S, 'sulfur');
    } else if (cm === 'element' && state.view === 'lattice') {
      html = sw(HYDRO_COLOR, 'hydrophobic') + sw(POLAR_COLOR, 'polar');
    }
    legend.innerHTML = html;
    legend.style.display = html ? 'block' : 'none';
  }

  function buildRibbon() {
    var html = '';
    for (var i = 0; i < state.residues.length; i++) {
      var res = state.residues[i];
      var col = '#333a4d';
      var pt = null;
      for (var p = 0; p < state.points.length; p++) {
        if (state.points[p].resIndex === i) { pt = state.points[p]; break; }
      }
      if (pt && pt.rgb) {
        col = 'rgb(' + pt.rgb[0] + ',' + pt.rgb[1] + ',' + pt.rgb[2] + ')';
      }
      html += '<div class="res' + (state.selected === i ? ' selected' : '') + '" data-i="' + i +
        '" style="background:' + col + '" title="' + res.resName + ' ' + res.num + '">' +
        res.letter + '</div>';
    }
    ribbon.innerHTML = html;
  }

  ribbon.addEventListener('click', function (e) {
    var chip = e.target.closest('.res');
    if (!chip) { return; }
    var idx = parseInt(chip.getAttribute('data-i'), 10);
    selectResidue(state.selected === idx ? null : idx);
  });

  // ---------------------------------------------------------------------
  // contact counting for live snapshots
  // ---------------------------------------------------------------------
  function countContacts(coords) {
    var n = coords.length, count = 0;
    for (var i = 0; i < n; i++) {
      for (var j = i + 2; j < n; j++) {
        var d = Math.abs(coords[i][0] - coords[j][0]) +
          Math.abs(coords[i][1] - coords[j][1]) +
          Math.abs(coords[i][2] - coords[j][2]);
        if (d === 1) { count++; }
      }
    }
    return count;
  }

  // ---------------------------------------------------------------------
  // charts
  // ---------------------------------------------------------------------
  function drawTrace() {
    var m = fitCanvas(traceChart);
    var ctx = traceChart.getContext('2d');
    ctx.setTransform(m.dpr, 0, 0, m.dpr, 0, 0);
    ctx.clearRect(0, 0, m.w, m.h);
    var data = null;
    if (state.result && state.result.trace) { data = state.result.trace; }
    else if (state.snapshots.length > 1) {
      data = state.snapshots.map(function (s) { return s.energy; });
    }
    if (!data || data.length < 2) {
      ctx.fillStyle = '#8fa0bd';
      ctx.font = '12px monospace';
      ctx.fillText(state.source === 'pdb' ? 'no energy trace for an imported PDB' : 'run a fold to see the trace',
        14, 24);
      return;
    }
    var lo = data[0], hi = data[0];
    for (var i = 0; i < data.length; i++) {
      if (data[i] < lo) { lo = data[i]; }
      if (data[i] > hi) { hi = data[i]; }
    }
    if (hi - lo < 1e-9) { hi = lo + 1; }
    var padL = 46, padR = 14, padT = 16, padB = 24;
    var W = m.w - padL - padR, H = m.h - padT - padB;
    function X(i) { return padL + W * i / (data.length - 1); }
    function Y(v) { return padT + H * (1 - (v - lo) / (hi - lo)); }
    ctx.strokeStyle = '#26304a';
    ctx.lineWidth = 1;
    ctx.strokeRect(padL, padT, W, H);
    if (lo < 0 && hi > 0) {
      ctx.strokeStyle = 'rgba(143,160,189,0.4)';
      ctx.beginPath();
      ctx.moveTo(padL, Y(0));
      ctx.lineTo(padL + W, Y(0));
      ctx.stroke();
    }
    ctx.strokeStyle = '#4ecdc4';
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    for (var k = 0; k < data.length; k++) {
      if (k === 0) { ctx.moveTo(X(k), Y(data[k])); } else { ctx.lineTo(X(k), Y(data[k])); }
    }
    ctx.stroke();
    if (state.source === 'fold' && state.snapshots.length) {
      var snap = state.snapshots[Math.min(state.scrub, state.snapshots.length - 1)];
      var xi = Math.min(snap.sweep - 1, data.length - 1);
      if (xi >= 0) {
        ctx.fillStyle = '#e0b13e';
        ctx.beginPath();
        ctx.arc(X(xi), Y(snap.energy), 4, 0, 6.2832);
        ctx.fill();
      }
    }
    ctx.fillStyle = '#8fa0bd';
    ctx.font = '10px monospace';
    ctx.fillText(hi.toFixed(1), 6, padT + 8);
    ctx.fillText(lo.toFixed(1), 6, padT + H);
    ctx.fillText('sweep 0', padL, m.h - 6);
    ctx.fillText('sweep ' + (data.length - 1), m.w - padR - 56, m.h - 6);
  }

  function drawRama() {
    var m = fitCanvas(ramaChart);
    var ctx = ramaChart.getContext('2d');
    ctx.setTransform(m.dpr, 0, 0, m.dpr, 0, 0);
    ctx.clearRect(0, 0, m.w, m.h);
    var padL = 40, padR = 12, padT = 12, padB = 26;
    var W = m.w - padL - padR, H = m.h - padT - padB;
    function X(phi) { return padL + W * (phi + 180) / 360; }
    function Y(psi) { return padT + H * (1 - (psi + 180) / 360); }
    for (var r = 0; r < F.RAMACHANDRAN_REGIONS.length; r++) {
      var reg = F.RAMACHANDRAN_REGIONS[r];
      var rgb = hexToRgb(reg.color);
      var boxes = [[reg.phi, reg.psi]];
      if (reg.psiAlt) { boxes.push([reg.phi, reg.psiAlt]); }
      for (var b = 0; b < boxes.length; b++) {
        var x0 = X(boxes[b][0][0]), x1 = X(boxes[b][0][1]);
        var y0 = Y(boxes[b][1][1]), y1 = Y(boxes[b][1][0]);
        ctx.fillStyle = 'rgba(' + rgb[0] + ',' + rgb[1] + ',' + rgb[2] + ',0.16)';
        ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
        ctx.strokeStyle = 'rgba(' + rgb[0] + ',' + rgb[1] + ',' + rgb[2] + ',0.5)';
        ctx.lineWidth = 1;
        ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
      }
    }
    ctx.strokeStyle = 'rgba(38,48,74,0.9)';
    ctx.lineWidth = 1;
    ctx.strokeRect(padL, padT, W, H);
    ctx.fillStyle = '#8fa0bd';
    ctx.font = '9px monospace';
    for (var a = -90; a <= 90; a += 90) {
      ctx.fillText(String(a), X(a) - 8, m.h - 8);
      ctx.fillText(String(a), 8, Y(a) + 3);
      ctx.strokeStyle = 'rgba(38,48,74,0.45)';
      ctx.beginPath();
      ctx.moveTo(X(a), padT); ctx.lineTo(X(a), padT + H);
      ctx.moveTo(padL, Y(a)); ctx.lineTo(padL + W, Y(a));
      ctx.stroke();
    }
    ctx.fillText('phi', m.w - 26, m.h - 8);
    ctx.fillText('psi', 6, padT + 8);
    var ds = state.dihedrals;
    if (!ds || !ds.length) {
      ctx.fillStyle = '#8fa0bd';
      ctx.font = '12px monospace';
      ctx.fillText('fold a chain to rebuild its backbone and plot the map', 16, 30);
      return;
    }
    for (var i = 0; i < ds.length; i++) {
      var d = ds[i];
      if (d.phi === null || d.psi === null || d.phi === undefined || d.psi === undefined) { continue; }
      var col = regionColor(d.region);
      var rgb2 = hexToRgb(col);
      ctx.fillStyle = 'rgba(' + rgb2[0] + ',' + rgb2[1] + ',' + rgb2[2] + ',0.95)';
      ctx.beginPath();
      ctx.arc(X(d.phi), Y(d.psi), 3.5, 0, 6.2832);
      ctx.fill();
      if (state.selected === i) {
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.4;
        ctx.beginPath();
        ctx.arc(X(d.phi), Y(d.psi), 6, 0, 6.2832);
        ctx.stroke();
      }
    }
  }

  // ---------------------------------------------------------------------
  // stats panel
  // ---------------------------------------------------------------------
  function setStat(el, value, unit) {
    el.innerHTML = value + (unit ? ' <small>' + unit + '</small>' : '');
  }

  function updateStats(snap) {
    if (!snap) { return; }
    setStat(stEnergy, snap.energy.toFixed(2), 'RT');
    if (state.source === 'fold' && snap.coords) {
      setStat(stContacts, String(countContacts(snap.coords)));
      setStat(stRg, F.radiusOfGyration(snap.coords).toFixed(2));
    }
    setStat(stExch, String(snap.exchanges));
    var acc = snap.proposed ? Math.round(100 * snap.accepted / snap.proposed) + '%' : '0%';
    setStat(stAcc, acc);
    setStat(stSweeps, String(snap.sweep));
  }

  function clearStats() {
    setStat(stEnergy, 'n/a');
    setStat(stContacts, 'n/a');
    setStat(stRg, 'n/a');
    setStat(stExch, 'n/a');
    setStat(stAcc, 'n/a');
    setStat(stSweeps, 'n/a');
  }

  // ---------------------------------------------------------------------
  // worker driven folding (with a chunked main thread fallback)
  // ---------------------------------------------------------------------
  var WORKER_SRC = [
    'var foldl = self.foldl;',
    'self.onmessage = function (e) {',
    '  var msg = e.data;',
    '  if (!msg || msg.cmd !== \'fold\') { return; }',
    '  try {',
    '    var it = foldl.remcGenerator({',
    '      sequence: msg.sequence, mode: msg.mode,',
    '      replicas: msg.replicas, sweeps: msg.sweeps,',
    '      stepsPerSweep: msg.stepsPerSweep, seed: msg.seed,',
    '      yieldEvery: msg.yieldEvery',
    '    });',
    '    var r = it.next();',
    '    while (!r.done) {',
    '      self.postMessage({ type: \'progress\', data: r.value });',
    '      r = it.next();',
    '    }',
    '    self.postMessage({ type: \'done\', data: r.value });',
    '  } catch (err) {',
    '    self.postMessage({ type: \'error\', message: String(err && err.message ? err.message : err) });',
    '  }',
    '};'
  ].join('\n');

  var worker = null;

  function getEngineSource() {
    return fetch('src/fold.js').then(function (res) {
      if (!res.ok) { throw new Error('could not load src/fold.js (HTTP ' + res.status + ')'); }
      return res.text();
    });
  }

  function ensureWorker() {
    if (worker) { return Promise.resolve(worker); }
    return getEngineSource().then(function (src) {
      // The engine exports to window.foldl in the browser path; give the
      // worker scope a window binding before the engine text runs.
      var blob = new Blob(['var window = self;\n' + src + '\n' + WORKER_SRC],
        { type: 'text/javascript' });
      worker = new Worker(URL.createObjectURL(blob));
      worker.onmessage = onWorkerMessage;
      worker.onerror = function (e) {
        showError('worker error: ' + (e.message || 'unknown'));
        finishRun();
      };
      return worker;
    });
  }

  function showError(text) {
    errMsg.textContent = text;
  }

  function clearError() {
    errMsg.textContent = '';
  }

  function startRun() {
    clearError();
    var raw = seqInput.value.replace(/[^A-Za-z]/g, '').toUpperCase();
    var seq;
    try {
      seq = F.validateSequence(raw);
    } catch (err) {
      showError('invalid sequence: ' + err.message);
      return;
    }
    var mode = modelSel.value === 'hp' ? 'hp' : 'mj';
    var replicas = parseInt(replicasR.value, 10);
    var sweeps = parseInt(sweepsR.value, 10);
    var seed = Math.floor(Math.random() * 900000) + 1000;
    var yieldEvery = Math.max(5, Math.round(sweeps / 40));
    var job = {
      sequence: seq, mode: mode, replicas: replicas, sweeps: sweeps,
      stepsPerSweep: Math.max(seq.length, 10), seed: seed, yieldEvery: yieldEvery
    };

    state.source = 'fold';
    state.seq = seq;
    state.mode = mode;
    state.snapshots = [];
    state.scrub = 0;
    state.userScrubbed = false;
    state.trace = null;
    state.result = null;
    state.selected = null;
    state.dihedrals = null;
    state.rebuilt = null;
    state.running = true;
    state.pdbAtoms = null;
    pdbMsg.textContent = '';
    selectResidue(null);
    pdbBtn.disabled = true;
    foldBtn.disabled = true;
    stopBtn.disabled = false;
    progressText.textContent = 'starting, seed ' + seed;
    progressBar.style.width = '0%';

    // snapshot 0: the straight extended chain
    var straight = [];
    for (var i = 0; i < seq.length; i++) { straight.push([i, 0, 0]); }
    state.snapshots.push({
      sweep: 0, energy: 0, coords: straight,
      exchanges: 0, accepted: 0, proposed: 0
    });
    refreshDisplay();
    drawTrace();
    drawRama();
    updateStats(state.snapshots[0]);

    ensureWorker().then(function (w) {
      w.postMessage({ cmd: 'fold', sequence: job.sequence, mode: job.mode,
        replicas: job.replicas, sweeps: job.sweeps,
        stepsPerSweep: job.stepsPerSweep, seed: job.seed, yieldEvery: job.yieldEvery });
    }).catch(function () {
      // no worker or no fetch (for example file:// open): run in chunks here
      try {
        startRunChunked(F.remcGenerator(job));
      } catch (err) {
        showError('fold failed: ' + (err.message || err));
        finishRun();
      }
    });
  }

  function startRunChunked(it) {
    function step() {
      if (!state.running) { return; }
      var t0 = Date.now();
      var r = null;
      try {
        do { r = it.next(); } while (!r.done && Date.now() - t0 < 24);
      } catch (err) {
        showError('fold failed: ' + (err.message || err));
        finishRun();
        return;
      }
      if (r.done) {
        onEngineMessage({ data: { type: 'done', data: r.value } });
      } else {
        onEngineMessage({ data: { type: 'progress', data: r.value } });
        setTimeout(step, 0);
      }
    }
    setTimeout(step, 0);
  }

  function onEngineMessage(msg) {
    if (!msg) { return; }
    if (msg.type === 'progress') {
      var d = msg.data;
      state.snapshots.push({
        sweep: d.sweep, energy: d.bestEnergy, coords: d.bestCoords,
        exchanges: d.exchanges, accepted: d.accepted, proposed: d.proposed
      });
      if (!state.userScrubbed) { state.scrub = state.snapshots.length - 1; }
      progressBar.style.width = Math.round(100 * d.sweep / d.sweeps) + '%';
      progressText.textContent = 'sweep ' + d.sweep + ' of ' + d.sweeps +
        ', best energy ' + d.bestEnergy.toFixed(2) + ' RT';
      updateStats(state.snapshots[state.scrub]);
      refreshDisplay();
      drawTrace();
    } else if (msg.type === 'done') {
      state.result = msg.data;
      state.trace = msg.data.trace;
      cacheRebuild();
      finishRun();
      updateStats(finalSnapshot());
      refreshDisplay();
      drawTrace();
      drawRama();
      progressText.textContent = 'finished in ' + msg.data.elapsedMs + ' ms, best energy ' +
        msg.data.best.energy.toFixed(2) + ' RT, ' + msg.data.best.contactCount + ' contacts';
      progressBar.style.width = '100%';
    } else if (msg.type === 'error') {
      showError('fold failed: ' + msg.message);
      finishRun();
    }
  }

  function onWorkerMessage(e) {
    onEngineMessage(e.data);
  }

  function finalSnapshot() {
    return state.snapshots.length ? state.snapshots[state.snapshots.length - 1] : null;
  }

  function finishRun() {
    state.running = false;
    foldBtn.disabled = false;
    stopBtn.disabled = true;
    if (state.snapshots.length > 1 || state.pdbAtoms) { pdbBtn.disabled = false; }
  }

  function stopRun() {
    if (worker) { worker.terminate(); worker = null; }
    finishRun();
    progressText.textContent = 'stopped at sweep ' +
      (finalSnapshot() ? finalSnapshot().sweep : 0);
  }

  function folderFromCoords(coords) {
    var lf = new F.LatticeFolder(state.seq, { mode: state.mode });
    var occ = new Map();
    for (var i = 0; i < coords.length; i++) {
      occ.set(F.LatticeFolder.key(coords[i][0], coords[i][1], coords[i][2]), i);
    }
    lf._commit(coords, occ, 0, null);
    return lf;
  }

  function cacheRebuild() {
    var snap = finalSnapshot();
    if (!snap || !state.seq) { return; }
    try {
      var rec = F.reconstructLatticeBackbone(folderFromCoords(snap.coords));
      state.rebuilt = rec;
      state.dihedrals = rec.dihedrals;
    } catch (err) {
      state.rebuilt = null;
      state.dihedrals = null;
    }
  }

  function refreshDisplay() {
    try {
      if (state.source === 'pdb') { applyColors(); return; }
      var snap = state.snapshots[Math.min(state.scrub, Math.max(state.snapshots.length - 1, 0))];
      if (!snap) { return; }
      if (showAll.value === 'yes') {
        var rec;
        if (state.rebuilt && snap === finalSnapshot()) { rec = state.rebuilt; }
        else { rec = F.reconstructLatticeBackbone(folderFromCoords(snap.coords)); }
        if (snap === finalSnapshot()) { state.dihedrals = rec.dihedrals; }
        setFromAtoms(rec.atoms);
      } else {
        setFromLattice(snap.coords);
      }
      applyColors();
    } catch (err) {
      showError('display error: ' + (err.message || err));
    } finally {
      updateScrubber();
    }
  }

  // ---------------------------------------------------------------------
  // scrub slider
  // ---------------------------------------------------------------------
  var scrubR = document.createElement('input');
  scrubR.type = 'range';
  scrubR.min = '0';
  scrubR.max = '0';
  scrubR.value = '0';
  scrubR.style.cssText = 'width:100%;accent-color:#4ecdc4;display:none;';
  var scrubLabel = document.createElement('div');
  scrubLabel.style.cssText = 'font-size:11px;color:#8fa0bd;margin-bottom:6px;display:none;';
  scrubLabel.textContent = 'trace scrubber: replay the best conformation over sweeps';
  var stageWrap = document.getElementById('stageWrap');
  stageWrap.insertBefore(scrubLabel, ribbon);
  stageWrap.insertBefore(scrubR, ribbon);

  scrubR.addEventListener('input', function () {
    if (!state.snapshots.length) { return; }
    state.userScrubbed = true;
    state.scrub = parseInt(scrubR.value, 10);
    updateStats(state.snapshots[state.scrub]);
    refreshDisplay();
    drawTrace();
  });

  function updateScrubber() {
    var n = state.snapshots.length;
    if (n > 1) {
      scrubR.style.display = 'block';
      scrubLabel.style.display = 'block';
      scrubR.max = String(n - 1);
      scrubR.value = String(Math.min(state.scrub, n - 1));
    }
  }

  // ---------------------------------------------------------------------
  // pdb export
  // ---------------------------------------------------------------------
  function exportPdb() {
    try {
      var text;
      if (state.source === 'pdb' && state.pdbAtoms) {
        text = F.atomsToPdb(state.pdbAtoms, { title: 'fold-lattice round trip export' });
      } else {
        var snap = finalSnapshot();
        if (!snap) { return; }
        var lf = folderFromCoords(snap.coords);
        var rec = F.reconstructLatticeBackbone(lf);
        text = F.atomsToPdb(rec.atoms, {
          title: 'fold-lattice ' + (state.seq || '').slice(0, 30),
          remarks: [
            'mode ' + state.mode,
            'energy ' + snap.energy.toFixed(2) + ' RT',
            'sweep ' + snap.sweep
          ]
        });
      }
      var blob = new Blob([text], { type: 'chemical/x-pdb' });
      var a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = 'fold-lattice.pdb';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      setTimeout(function () { URL.revokeObjectURL(a.href); }, 2000);
    } catch (err) {
      showError('export failed: ' + (err.message || err));
    }
  }

  // ---------------------------------------------------------------------
  // pdb import
  // ---------------------------------------------------------------------
  pdbFile.addEventListener('change', function () {
    pdbMsg.textContent = '';
    var file = pdbFile.files && pdbFile.files[0];
    if (!file) { return; }
    file.text().then(function (text) {
      var parsed = F.parsePdb(text);
      var dihedrals = F.computeDihedrals(parsed.atoms);
      state.source = 'pdb';
      state.pdbAtoms = parsed.atoms;
      state.dihedrals = dihedrals;
      state.result = null;
      state.trace = null;
      state.seq = null;
      state.snapshots = [];
      state.selected = null;
      state.running = false;
      selectResidue(null);
      setFromAtoms(parsed.atoms);
      applyColors();
      var clashes = F.findStericClashes(parsed.atoms, {
        excludePair: function (i, j) {
          var a = parsed.atoms[i], b = parsed.atoms[j];
          return a.chainID === b.chainID && Math.abs(a.resSeq - b.resSeq) <= 1;
        }
      });
      pdbMsg.textContent = 'parsed ' + parsed.atoms.length + ' atoms, ' +
        state.residues.length + ' residues, ' + clashes.length +
        ' non bonded steric clashes, ' + parsed.numModelsRead + ' model(s) read';
      pdbBtn.disabled = false;
      foldBtn.disabled = false;
      stopBtn.disabled = true;
      clearStats();
      progressText.textContent = 'showing imported PDB, ' + parsed.atoms.length + ' atoms';
      progressBar.style.width = '0%';
      drawTrace();
      drawRama();
      showAll.value = 'yes';
    }).catch(function (err) {
      pdbMsg.textContent = 'parse failed: ' + (err.message || err);
    });
  });

  // ---------------------------------------------------------------------
  // controls wiring
  // ---------------------------------------------------------------------
  foldBtn.addEventListener('click', startRun);
  stopBtn.addEventListener('click', stopRun);
  pdbBtn.addEventListener('click', exportPdb);

  replicasR.addEventListener('input', function () { replicasV.textContent = replicasR.value; });
  sweepsR.addEventListener('input', function () { sweepsV.textContent = sweepsR.value; });

  presetSel.addEventListener('change', function () {
    var v = presetSel.value;
    if (!v) { return; }
    seqInput.value = v;
    modelSel.value = (v === 'SSSSAAAASSSSAAAASSSAA') ? 'hp' : 'mj';
  });

  colorSel.addEventListener('change', function () {
    applyColors();
    drawRama();
  });

  showAll.addEventListener('change', function () {
    if (state.source === 'pdb') { return; }
    refreshDisplay();
    drawRama();
  });

  window.addEventListener('resize', function () {
    drawTrace();
    drawRama();
  });

  // ---------------------------------------------------------------------
  // render loop
  // ---------------------------------------------------------------------
  function loop() {
    drawScene();
    requestAnimationFrame(loop);
  }

  // initial idle display: straight chain so the stage is not empty
  (function initIdle() {
    state.seq = 'GAGSLYQGSSLAQSGFLSSHHSGF';
    state.mode = 'mj';
    var straight = [];
    for (var i = 0; i < state.seq.length; i++) { straight.push([i, 0, 0]); }
    setFromLattice(straight);
    applyColors();
  })();

  loop();
})();
