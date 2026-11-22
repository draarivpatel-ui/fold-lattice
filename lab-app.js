/*
 * fold-lattice lab workbench UI. Drives the HP comparator, folding funnel
 * landscape, contact map, mutation stress scan, checksummed audit export,
 * URL hash bookmarks and the in page self test runner. All rendering is
 * hand written SVG or DOM text nodes; no dependencies, no network.
 */
(function () {
  'use strict';

  var F = window.foldl;
  var L = window.foldlLab;
  if (!F || !L) { return; }

  function $(id) { return document.getElementById(id); }

  var hpSeqInput = $('hpSeqInput');
  var hpPresetSel = $('hpPresetSel');
  var dimSel = $('dimSel');
  var seedInput = $('seedInput');
  var samplesInput = $('samplesInput');
  var mcStepsInput = $('mcStepsInput');
  var tempInput = $('tempInput');
  var saSweepsInput = $('saSweepsInput');
  var t0Input = $('t0Input');
  var t1Input = $('t1Input');
  var compareBtn = $('compareBtn');
  var labErr = $('labErr');
  var labTableWrap = $('labTableWrap');
  var labStatus = $('labStatus');
  var funnelSvg = $('funnelSvg');
  var funnelAxisSel = $('funnelAxisSel');
  var funnelInfo = $('funnelInfo');
  var cmapSvg = $('cmapSvg');
  var cmapInfo = $('cmapInfo');
  var mutBtn = $('mutBtn');
  var mutSvg = $('mutSvg');
  var mutInfo = $('mutInfo');
  var auditJsonBtn = $('auditJsonBtn');
  var auditCsvBtn = $('auditCsvBtn');
  var auditMdBtn = $('auditMdBtn');
  var bookmarkBtn = $('bookmarkBtn');
  var verifyBtn = $('verifyBtn');
  var tamperBtn = $('tamperBtn');
  var auditSum = $('auditSum');
  var auditVerdict = $('auditVerdict');
  var auditPreview = $('auditPreview');
  var selfTestBtn = $('selfTestBtn');
  var selfTestList = $('selfTestList');
  var selfTestInfo = $('selfTestInfo');

  var SVGNS = 'http://www.w3.org/2000/svg';

  var ALGO_COLORS = {
    'random-walk': '#8fa0bd',
    'greedy-pivot': '#e0b13e',
    'metropolis': '#4e7de0',
    'annealing': '#4ecdc4',
    'exact': '#e0564e'
  };

  var labState = {
    comparison: null,
    selected: -1,
    funnelCloud: null,
    mutationScan: null,
    auditReport: null,
    auditJson: ''
  };

  // ---------------------------------------------------------------------
  // svg helpers
  // ---------------------------------------------------------------------
  function el(name, attrs) {
    var e = document.createElementNS(SVGNS, name);
    if (attrs) {
      for (var k in attrs) { e.setAttribute(k, attrs[k]); }
    }
    return e;
  }

  function clear(node) {
    while (node.firstChild) { node.removeChild(node.firstChild); }
  }

  function setErr(text) { labErr.textContent = text; }

  // ---------------------------------------------------------------------
  // comparator
  // ---------------------------------------------------------------------
  // The HTML min/max attributes are hints only; every bound is enforced
  // here so pasted values, bookmarks and typed input all stay sane.
  function intIn(el, name, lo, hi) {
    var v = parseInt(el.value, 10);
    if (!Number.isInteger(v) || v < lo || v > hi) {
      throw new Error(name + ' must be an integer between ' + lo + ' and ' + hi);
    }
    return v;
  }
  function numIn(el, name, lo, hi) {
    var v = parseFloat(el.value);
    if (!Number.isFinite(v) || v < lo || v > hi) {
      throw new Error(name + ' must be a number between ' + lo + ' and ' + hi);
    }
    return v;
  }

  function readControls() {
    var seq = L.validateHpSequence(hpSeqInput.value);
    if (seq.length > L.BOOKMARK_MAX_SEQ) {
      throw new Error('sequence is capped at ' + L.BOOKMARK_MAX_SEQ + ' residues for the workbench');
    }
    var dim = parseInt(dimSel.value, 10);
    var seed = intIn(seedInput, 'seed', 0, 4294967295);
    var samples = intIn(samplesInput, 'rw samples', 10, 5000);
    var mcSteps = intIn(mcStepsInput, 'mc steps', 100, 200000);
    var temp = numIn(tempInput, 'mc temp', 0, 10);
    var saSweeps = intIn(saSweepsInput, 'sa sweeps', 5, 2000);
    var T0 = numIn(t0Input, 'sa T0', 0.05, 50);
    var T1 = numIn(t1Input, 'sa T1', 0, 10);
    return {
      sequence: seq, dim: dim, seed: seed,
      samples: samples, mcSteps: mcSteps, temperature: temp,
      saSweeps: saSweeps, T0: T0, T1: T1
    };
  }

  function runComparator() {
    setErr('');
    var opts;
    try {
      opts = readControls();
    } catch (err) {
      setErr(err.message || String(err));
      return;
    }
    compareBtn.disabled = true;
    labStatus.textContent = 'running comparator...';
    // Let the status paint before the (fast but nonzero) synchronous bench.
    setTimeout(function () {
      try {
        var t0 = Date.now();
        var cmp = L.compareAlgorithms(opts);
        var cmpMs = Date.now() - t0;
        labState.comparison = cmp;
        labState.auditReport = null;
        labState.auditJson = '';
        auditPreview.textContent = '';
        auditSum.textContent = '';
        auditVerdict.textContent = '';
        // Default selection: exact ground state when present, else best.
        var sel = 0;
        if (cmp.exact) { sel = cmp.results.indexOf(cmp.exact); }
        else {
          for (var i = 1; i < cmp.results.length; i++) {
            if (cmp.results[i].best.energy < cmp.results[sel].best.energy) { sel = i; }
          }
        }
        labState.funnelCloud = L.sampleFunnel(cmp.sequence, cmp.dim, 260, cmp.seed + 99);
        renderTable();
        selectResult(sel);
        var note = cmp.exactSkipped
          ? 'exact enumeration skipped: chain length ' + cmp.sequence.length +
            ' exceeds the cap (' + L.EXACT_MAX_LEN[cmp.dim] + ' in ' + cmp.dim + 'D)'
          : 'exact ground state ' + cmp.exact.best.energy + ' over ' + cmp.exact.states + ' states';
        labStatus.textContent = 'compared ' + cmp.results.length + ' algorithms in ' +
          cmpMs + ' ms, ' + note;
      } catch (err) {
        setErr('compare failed: ' + (err.message || err));
      } finally {
        compareBtn.disabled = false;
      }
    }, 30);
  }

  function fmt(v, d) {
    if (v === null || v === undefined) { return 'n/a'; }
    return Number(v).toFixed(d === undefined ? 2 : d);
  }

  function renderTable() {
    var cmp = labState.comparison;
    clear(labTableWrap);
    if (!cmp) { return; }
    var table = document.createElement('table');
    table.className = 'labtable';
    var head = document.createElement('tr');
    ['algorithm', 'best E', 'contacts', 'H-H', 'R_g', 'R_ee', 'gap', 'ms', ''].forEach(function (h) {
      var th = document.createElement('th');
      th.textContent = h;
      head.appendChild(th);
    });
    table.appendChild(head);
    cmp.results.forEach(function (r, idx) {
      var tr = document.createElement('tr');
      tr.setAttribute('data-i', String(idx));
      tr.className = idx === labState.selected ? 'sel' : '';
      var cells = [
        r.label + (r.exact ? ' (ground truth)' : ''),
        String(r.best.energy),
        String(r.metrics.contactCount),
        String(r.metrics.hhContacts),
        fmt(r.metrics.rg),
        fmt(r.metrics.ree),
        r.optimalityGap === null ? 'n/a' : '+' + r.optimalityGap,
        String(r.ms)
      ];
      cells.forEach(function (c) {
        var td = document.createElement('td');
        td.textContent = c;
        tr.appendChild(td);
      });
      var tdBtn = document.createElement('td');
      var btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'mini';
      btn.textContent = 'view';
      btn.setAttribute('aria-label', 'show ' + r.label + ' conformation on the 3D stage');
      btn.addEventListener('click', function (ev) {
        ev.stopPropagation();
        selectResult(idx);
        viewOnStage(r);
      });
      tdBtn.appendChild(btn);
      tr.appendChild(tdBtn);
      tr.addEventListener('click', function () { selectResult(idx); });
      table.appendChild(tr);
    });
    labTableWrap.appendChild(table);
  }

  function selectResult(idx) {
    labState.selected = idx;
    var cmp = labState.comparison;
    if (!cmp) { return; }
    var rows = labTableWrap.querySelectorAll('tr[data-i]');
    for (var i = 0; i < rows.length; i++) {
      rows[i].className = parseInt(rows[i].getAttribute('data-i'), 10) === idx ? 'sel' : '';
    }
    var r = cmp.results[idx];
    renderContactMap(r);
    renderFunnel();
    labState.mutationScan = null;
    clear(mutSvg);
    mutInfo.textContent = 'run a mutation scan on the selected conformation (' + r.label + ')';
    buildAudit();
  }

  function viewOnStage(r) {
    if (!window.foldlStage || !labState.comparison) { return; }
    window.foldlStage.showLatticeFold(
      r.best.coords.map(function (c) { return c.slice(); }),
      L.hpToAa(labState.comparison.sequence),
      r.best.energy
    );
    var wrap = document.getElementById('stageWrap');
    if (wrap && wrap.scrollIntoView) { wrap.scrollIntoView(); }
  }

  // ---------------------------------------------------------------------
  // folding funnel chart (SVG)
  // ---------------------------------------------------------------------
  function drawAxes(svg, xLabel, yLabel, xTicks, yTicks, W, H, pad) {
    svg.appendChild(el('line', {
      x1: pad.l, y1: pad.t, x2: pad.l, y2: H - pad.b,
      stroke: '#26304a', 'stroke-width': 1
    }));
    svg.appendChild(el('line', {
      x1: pad.l, y1: H - pad.b, x2: W - pad.r, y2: H - pad.b,
      stroke: '#26304a', 'stroke-width': 1
    }));
    xTicks.forEach(function (t) {
      var tx = el('text', { x: t.px, y: H - pad.b + 14, 'text-anchor': 'middle', 'font-size': 9, fill: '#8fa0bd' });
      tx.textContent = t.label;
      svg.appendChild(tx);
    });
    yTicks.forEach(function (t) {
      var ty = el('text', { x: pad.l - 5, y: t.py + 3, 'text-anchor': 'end', 'font-size': 9, fill: '#8fa0bd' });
      ty.textContent = t.label;
      svg.appendChild(ty);
    });
    var xl = el('text', { x: W - pad.r - 8, y: H - 6, 'text-anchor': 'end', 'font-size': 10, fill: '#8fa0bd' });
    xl.textContent = xLabel;
    svg.appendChild(xl);
    var yl = el('text', { x: pad.l - 5, y: pad.t - 6, 'text-anchor': 'end', 'font-size': 10, fill: '#8fa0bd' });
    yl.textContent = yLabel;
    svg.appendChild(yl);
  }

  function renderFunnel() {
    var cmp = labState.comparison;
    clear(funnelSvg);
    if (!cmp) { return; }
    var xMode = funnelAxisSel.value; // 'rg' or 'step'
    var W = 560, H = 300;
    var pad = { l: 44, r: 12, t: 18, b: 30 };
    funnelSvg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);

    // Collect every point to fix the axis ranges. The random cloud is a
    // compactness backdrop and only makes sense on the R_g axis.
    var xs = [], ys = [];
    var cloud = (labState.funnelCloud && xMode === 'rg') ? labState.funnelCloud.points : [];
    cloud.forEach(function (p) {
      xs.push(p.Rg);
      ys.push(p.E);
    });
    cmp.results.forEach(function (r) {
      (r.points || []).forEach(function (p) {
        xs.push(xMode === 'rg' ? p.Rg : (p.step === undefined ? 0 : p.step));
        ys.push(p.E);
      });
    });
    if (!xs.length) { xs = [0, 1]; ys = [0, 1]; }
    var xLo = Math.min.apply(null, xs), xHi = Math.max.apply(null, xs);
    var yLo = Math.min.apply(null, ys), yHi = Math.max.apply(null, ys);
    if (xMode === 'step') { xLo = 0; }
    if (xHi - xLo < 1e-9) { xHi = xLo + 1; }
    if (yHi - yLo < 1e-9) { yHi = yLo + 1; }
    yLo -= 0.4; yHi += 0.2;

    function X(v) { return pad.l + (W - pad.l - pad.r) * (v - xLo) / (xHi - xLo); }
    function Y(v) { return pad.t + (H - pad.t - pad.b) * (1 - (v - yLo) / (yHi - yLo)); }

    var xTicks = [], yTicks = [];
    for (var i = 0; i <= 4; i++) {
      var xv = xLo + (xHi - xLo) * i / 4;
      var yv = yLo + (yHi - yLo) * i / 4;
      xTicks.push({ px: X(xv), label: xMode === 'rg' ? xv.toFixed(1) : String(Math.round(xv)) });
      yTicks.push({ py: Y(yv), label: yv.toFixed(1) });
    }
    drawAxes(funnelSvg, xMode === 'rg' ? 'R_g' : 'step', 'E', xTicks, yTicks, W, H, pad);

    // Random SAW cloud behind the trajectories (compactness view only).
    cloud.forEach(function (p) {
      funnelSvg.appendChild(el('circle', {
        cx: X(p.Rg), cy: Y(p.E), r: 2,
        fill: 'rgba(143,160,189,0.28)'
      }));
    });

    // Exact ground state marker line.
    if (cmp.exact && !cmp.exact.truncated) {
      var ge = cmp.exact.best.energy;
      var line = el('line', {
        x1: pad.l, y1: Y(ge), x2: W - pad.r, y2: Y(ge),
        stroke: ALGO_COLORS.exact, 'stroke-width': 1, 'stroke-dasharray': '5 4', opacity: 0.7
      });
      funnelSvg.appendChild(line);
      var gt = el('text', { x: W - pad.r - 4, y: Y(ge) - 5, 'text-anchor': 'end', 'font-size': 10, fill: ALGO_COLORS.exact });
      gt.textContent = 'exact ground state ' + ge;
      funnelSvg.appendChild(gt);
    }

    // Trajectory polylines and points per algorithm.
    cmp.results.forEach(function (r) {
      var col = ALGO_COLORS[r.algorithm] || '#dfe6f3';
      var pts = (r.points || []).map(function (p) {
        var xv = xMode === 'rg' ? p.Rg : (p.step === undefined ? 0 : p.step);
        return [X(xv), Y(p.E)];
      });
      if (pts.length > 1) {
        var path = pts.map(function (p, i) { return (i ? 'L' : 'M') + p[0].toFixed(1) + ' ' + p[1].toFixed(1); }).join(' ');
        funnelSvg.appendChild(el('path', {
          d: path, fill: 'none', stroke: col, 'stroke-width': 1.4, opacity: 0.75
        }));
      }
      pts.forEach(function (p) {
        funnelSvg.appendChild(el('circle', { cx: p[0], cy: p[1], r: 2.1, fill: col, opacity: 0.85 }));
      });
      // Best marker.
      var bm = r.metrics;
      var bx = xMode === 'rg' ? bm.rg : (r.steps || (r.trajectory ? r.trajectory.length : 0));
      var bc = el('circle', {
        cx: X(bx), cy: Y(bm.energy), r: 5,
        fill: 'none', stroke: col, 'stroke-width': 2
      });
      funnelSvg.appendChild(bc);
      var bt = el('title');
      bt.textContent = r.label + ' best E ' + bm.energy;
      bc.appendChild(bt);
    });
    funnelInfo.textContent = cmp.results.length + ' algorithms, cloud of ' +
      cloud.length + ' random walks; ringed markers are run bests';
  }

  // ---------------------------------------------------------------------
  // contact map
  // ---------------------------------------------------------------------
  var CMAP_COLORS = {
    0: '#141a26',
    1: '#333a4d',
    2: '#4ecdc4',
    3: '#5c6b8a',
    4: '#1c2333'
  };

  function renderContactMap(result) {
    var cmp = labState.comparison;
    clear(cmapSvg);
    if (!cmp || !result) { return; }
    var seq = cmp.sequence;
    var m = L.contactMatrix(seq, result.best.coords, cmp.dim);
    var n = m.n;
    var cell = Math.max(14, Math.min(30, Math.floor(420 / n)));
    var padL = 22, padT = 18;
    var W = padL + n * cell + 8, H = padT + n * cell + 24;
    cmapSvg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    for (var a = 0; a < n; a++) {
      for (var b = 0; b < n; b++) {
        var code = m.cells[a * n + b];
        var rect = el('rect', {
          x: padL + b * cell, y: padT + a * cell,
          width: cell - 1, height: cell - 1,
          fill: CMAP_COLORS[code], rx: 2,
          'data-i': a, 'data-j': b
        });
        var t = el('title');
        t.textContent = 'residues ' + (a + 1) + ' (' + seq[a] + ') and ' + (b + 1) + ' (' + seq[b] + ')';
        rect.appendChild(t);
        cmapSvg.appendChild(rect);
      }
      var lbl = el('text', { x: padL - 4, y: padT + a * cell + cell * 0.7, 'text-anchor': 'end', 'font-size': 9, fill: '#8fa0bd' });
      lbl.textContent = seq[a] + (a + 1);
      cmapSvg.appendChild(lbl);
      var lbl2 = el('text', { x: padL + a * cell + cell * 0.3, y: padT - 6, 'font-size': 9, fill: '#8fa0bd' });
      lbl2.textContent = seq[a];
      cmapSvg.appendChild(lbl2);
    }
    cmapInfo.textContent = seq.length + ' residues, ' + m.hhCount + ' non bonded H-H contacts (teal), ' +
      (m.contactCount - m.hhCount) + ' other contacts, energy ' + (-m.hhCount) +
      '. click a cell for pair details.';
    cmapSvg.onclick = function (e) {
      var target = e.target;
      if (!target || !target.getAttribute) { return; }
      var i = target.getAttribute('data-i');
      var j = target.getAttribute('data-j');
      if (i === null || j === null) { return; }
      i = parseInt(i, 10); j = parseInt(j, 10);
      var code = m.cells[i * m.n + j];
      var kinds = {
        0: 'not adjacent', 1: 'chain bond', 2: 'non bonded H-H contact',
        3: 'non bonded contact (not H-H)', 4: 'same residue'
      };
      cmapInfo.textContent = 'residues ' + (i + 1) + ' (' + seq[i] + ') and ' + (j + 1) +
        ' (' + seq[j] + '): ' + kinds[code];
      var marks = cmapSvg.querySelectorAll('rect.sel');
      for (var q = 0; q < marks.length; q++) { marks[q].classList.remove('sel'); }
      target.classList.add('sel');
      renderMutationHighlight(i, j);
    };
  }

  // ---------------------------------------------------------------------
  // mutation stress scan
  // ---------------------------------------------------------------------
  function runMutationScan() {
    var cmp = labState.comparison;
    if (!cmp || labState.selected < 0) { return; }
    var r = cmp.results[labState.selected];
    labState.mutationScan = L.mutationScan(cmp.sequence, r.best.coords, cmp.dim);
    labState.mutationScan.algorithm = r.algorithm;
    renderMutationChart();
    buildAudit();
  }

  function renderMutationChart() {
    var scan = labState.mutationScan;
    var cmp = labState.comparison;
    clear(mutSvg);
    if (!scan || !cmp) { return; }
    var n = scan.perResidue.length;
    var W = 560, H = 200;
    var pad = { l: 34, r: 10, t: 16, b: 30 };
    mutSvg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    var lo = 0, hi = 0;
    scan.perResidue.forEach(function (e) {
      if (e.dE < lo) { lo = e.dE; }
      if (e.dE > hi) { hi = e.dE; }
    });
    if (hi - lo < 1e-9) { hi = lo + 1; }
    var barW = (W - pad.l - pad.r) / n;
    function Y(v) { return pad.t + (H - pad.t - pad.b) * (1 - (v - lo) / (hi - lo)); }
    mutSvg.appendChild(el('line', {
      x1: pad.l, y1: Y(0), x2: W - pad.r, y2: Y(0), stroke: '#8fa0bd', 'stroke-width': 1
    }));
    var yl = el('text', { x: pad.l - 5, y: pad.t + 4, 'text-anchor': 'end', 'font-size': 9, fill: '#8fa0bd' });
    yl.textContent = '+' + hi;
    mutSvg.appendChild(yl);
    var yl2 = el('text', { x: pad.l - 5, y: Y(lo) + 4, 'text-anchor': 'end', 'font-size': 9, fill: '#8fa0bd' });
    yl2.textContent = String(lo);
    mutSvg.appendChild(yl2);
    scan.perResidue.forEach(function (e) {
      var x = pad.l + e.index * barW;
      var y0 = Y(0), y1 = Y(e.dE);
      var rect = el('rect', {
        x: x + 1.5, y: Math.min(y0, y1),
        width: Math.max(1, barW - 3), height: Math.max(1.5, Math.abs(y1 - y0)),
        fill: e.dE > 0 ? '#e0564e' : (e.dE < 0 ? '#4ecdc4' : '#333a4d'), rx: 1.5
      });
      var t = el('title');
      t.textContent = 'residue ' + (e.index + 1) + ' ' + e.from + ' to ' + e.to + ', deltaE ' + e.dE;
      rect.appendChild(t);
      mutSvg.appendChild(rect);
      var lab = el('text', {
        x: x + barW / 2, y: H - pad.b + 14, 'text-anchor': 'middle',
        'font-size': Math.min(9, barW * 0.55), fill: '#8fa0bd'
      });
      lab.textContent = cmp.sequence[e.index] + (e.index + 1);
      mutSvg.appendChild(lab);
    });
    var w = scan.worst;
    mutInfo.textContent = 'E0 ' + scan.e0 + ': ' + scan.destabilizing + ' destabilizing, ' +
      scan.stabilizing + ' stabilizing, ' + scan.neutral + ' neutral' +
      (w ? ', worst residue ' + (w.index + 1) + ' ' + w.from + ' to ' + w.to + ' (deltaE +' + w.dE + ')' : '');
  }

  function renderMutationHighlight(i, j) {
    var scan = labState.mutationScan;
    if (!scan) { return; }
    var a = scan.perResidue[i], b = scan.perResidue[j];
    if (a && b) {
      mutInfo.textContent = 'pair ' + (i + 1) + ',' + (j + 1) + ': mutating residue ' + (i + 1) +
        ' gives deltaE ' + a.dE + ', residue ' + (j + 1) + ' gives deltaE ' + b.dE;
    }
  }

  // ---------------------------------------------------------------------
  // audit report
  // ---------------------------------------------------------------------
  function buildAudit() {
    var cmp = labState.comparison;
    if (!cmp) { return; }
    var payload = {
      sequence: cmp.sequence,
      dimension: cmp.dim,
      seed: cmp.seed,
      groundEnergy: cmp.exact && !cmp.exact.truncated ? cmp.exact.best.energy : null,
      results: cmp.results.map(function (r) {
        return {
          algorithm: r.algorithm,
          label: r.label,
          energy: r.best.energy,
          contactCount: r.metrics.contactCount,
          hhContacts: r.metrics.hhContacts,
          rg: Math.round(r.metrics.rg * 10000) / 10000,
          ree: Math.round(r.metrics.ree * 10000) / 10000,
          optimalityGap: r.optimalityGap,
          states: r.states === undefined ? null : r.states,
          ms: r.ms
        };
      }),
      mutationScan: labState.mutationScan ? {
        algorithm: labState.mutationScan.algorithm,
        e0: labState.mutationScan.e0,
        destabilizing: labState.mutationScan.destabilizing,
        stabilizing: labState.mutationScan.stabilizing,
        neutral: labState.mutationScan.neutral,
        worst: labState.mutationScan.worst,
        perResidue: labState.mutationScan.perResidue
      } : null,
      modelDisclosure: 'Dill HP lattice model: simplified statistical physics toy, not atomistic prediction'
    };
    labState.auditReport = L.buildAuditReport(payload);
    labState.auditJson = L.auditToJson(labState.auditReport);
    auditPreview.textContent = labState.auditJson;
    auditSum.textContent = 'sha256 ' + labState.auditReport.checksum;
    auditVerdict.textContent = '';
  }

  function download(name, text, type) {
    var blob = new Blob([text], { type: type });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 2000);
  }

  function verifyNow() {
    if (!labState.auditJson) { auditVerdict.textContent = 'run the comparator first'; return; }
    var v = L.verifyAuditText(labState.auditJson);
    auditVerdict.textContent = v.ok ? 'checksum verified: payload is intact' : 'FAILED: ' + v.reason;
    auditVerdict.className = v.ok ? 'verdict pass' : 'verdict fail';
  }

  function tamperNow() {
    if (!labState.auditJson) { auditVerdict.textContent = 'run the comparator first'; return; }
    var tampered = L.tamperAuditText(labState.auditJson);
    var v = L.verifyAuditText(tampered);
    auditVerdict.textContent = v.ok
      ? 'tamper check FAILED: altered report still verified'
      : 'tamper detected: ' + v.reason;
    auditVerdict.className = v.ok ? 'verdict fail' : 'verdict pass';
    auditPreview.textContent = tampered;
  }

  function bookmarkNow() {
    var cmp = labState.comparison;
    if (!cmp) { auditVerdict.textContent = 'run the comparator first'; return; }
    var mark = L.encodeBookmark({ sequence: cmp.sequence, dim: cmp.dim, seed: cmp.seed });
    try {
      history.replaceState(null, '', '#' + mark);
    } catch (err) {
      /* file:// or odd embedding: still show the bookmark */
    }
    auditVerdict.textContent = 'bookmark set: ' + location.pathname.split('/').pop() + '#' + mark;
    auditVerdict.className = 'verdict pass';
  }

  // ---------------------------------------------------------------------
  // self tests
  // ---------------------------------------------------------------------
  function runSelfTestsUI() {
    selfTestBtn.disabled = true;
    selfTestInfo.textContent = 'running self tests...';
    setTimeout(function () {
      var out = L.runSelfTests();
      clear(selfTestList);
      out.results.forEach(function (r) {
        var li = document.createElement('li');
        li.className = r.pass ? 'pass' : 'fail';
        li.textContent = (r.pass ? 'PASS ' : 'FAIL ') + r.name + ' - ' + r.detail;
        selfTestList.appendChild(li);
      });
      selfTestInfo.textContent = out.passed + '/' + out.total + ' self checks passed';
      selfTestBtn.disabled = false;
    }, 20);
  }

  // ---------------------------------------------------------------------
  // bookmark restore + wiring
  // ---------------------------------------------------------------------
  function restoreBookmark() {
    var b = L.decodeBookmark(location.hash);
    if (!b) { return false; }
    hpSeqInput.value = b.sequence;
    dimSel.value = String(b.dim);
    seedInput.value = String(b.seed);
    return true;
  }

  compareBtn.addEventListener('click', runComparator);
  mutBtn.addEventListener('click', runMutationScan);
  funnelAxisSel.addEventListener('change', renderFunnel);
  auditJsonBtn.addEventListener('click', function () {
    if (labState.auditJson) { download('fold-lattice-audit.json', labState.auditJson, 'application/json'); }
  });
  auditCsvBtn.addEventListener('click', function () {
    if (labState.auditReport) { download('fold-lattice-audit.csv', L.auditToCsv(labState.auditReport), 'text/csv'); }
  });
  auditMdBtn.addEventListener('click', function () {
    if (labState.auditReport) { download('fold-lattice-audit.md', L.auditToMarkdown(labState.auditReport), 'text/markdown'); }
  });
  bookmarkBtn.addEventListener('click', bookmarkNow);
  verifyBtn.addEventListener('click', verifyNow);
  tamperBtn.addEventListener('click', tamperNow);
  selfTestBtn.addEventListener('click', runSelfTestsUI);

  hpPresetSel.addEventListener('change', function () {
    var v = hpPresetSel.value;
    if (!v) { return; }
    var parts = v.split(':');
    hpSeqInput.value = parts[0];
    if (parts[1]) { dimSel.value = parts[1]; }
  });

  var restored = restoreBookmark();
  if (restored) {
    labStatus.textContent = 'restored from bookmark, running comparator';
    runComparator();
  }
  runSelfTestsUI();
})();
