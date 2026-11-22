'use strict';

/*
 * Zero dependency headless Chromium audit for fold-lattice.
 *
 * Drives a real Chrome over the DevTools protocol (Node >= 22 ships a
 * global WebSocket, so no npm packages are needed). For each viewport it
 * loads index.html over file://, records console errors, uncaught page
 * errors, external network requests and horizontal overflow, then
 * exercises the HP workbench: comparator, folding funnel, contact map,
 * mutation scan, checksummed audit export, tamper check, URL bookmark,
 * in page self tests and the view on stage bridge.
 *
 * Usage: node tools/browser-audit.js
 *        CHROME_PATH=/path/to/chrome node tools/browser-audit.js
 *
 * Requires Node 22+ (global WebSocket) and a Chrome or Chromium binary.
 */

const { spawn } = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const ROOT = path.join(__dirname, '..');
const PAGE_URL = 'file://' + path.join(ROOT, 'index.html');
const VIEWPORTS = [
  { name: 'desktop 1440x900', width: 1440, height: 900, mobile: false },
  { name: 'mobile 390x844', width: 390, height: 844, mobile: true }
];

const CHROME_CANDIDATES = [
  process.env.CHROME_PATH,
  '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  '/Applications/Chromium.app/Contents/MacOS/Chromium',
  '/usr/bin/google-chrome',
  '/usr/bin/chromium',
  '/usr/bin/chromium-browser'
].filter(Boolean);

function findChrome() {
  for (const p of CHROME_CANDIDATES) {
    try { fs.accessSync(p, fs.constants.X_OK); return p; } catch (e) { /* next */ }
  }
  throw new Error('no Chrome or Chromium binary found; set CHROME_PATH');
}

function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }

/* ------------------------------ CDP client ----------------------------- */

class Cdp {
  constructor(wsUrl) {
    this.wsUrl = wsUrl;
    this.nextId = 1;
    this.pending = new Map();
    this.listeners = new Map();
    this.ready = new Promise((resolve, reject) => {
      this.ws = new WebSocket(wsUrl);
      this.ws.onopen = resolve;
      this.ws.onerror = (e) => reject(new Error('CDP websocket failed: ' + (e.message || 'error')));
      this.ws.onmessage = (ev) => this._onMessage(ev.data);
    });
  }

  _onMessage(data) {
    let msg;
    try { msg = JSON.parse(data); } catch (e) { return; }
    if (msg.id !== undefined) {
      const p = this.pending.get(msg.id);
      if (p) {
        this.pending.delete(msg.id);
        if (msg.error) { p.reject(new Error(msg.method + ': ' + JSON.stringify(msg.error))); }
        else { p.resolve(msg.result); }
      }
      return;
    }
    const ls = this.listeners.get(msg.method);
    if (ls) { for (const fn of ls) { fn(msg.params || {}); } }
  }

  on(method, fn) {
    if (!this.listeners.has(method)) { this.listeners.set(method, []); }
    this.listeners.get(method).push(fn);
  }

  async send(method, params) {
    await this.ready;
    const id = this.nextId++;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject, method });
      this.ws.send(JSON.stringify({ id, method, params: params || {} }));
    });
  }

  async evaluate(expression, awaitPromise) {
    const r = await this.send('Runtime.evaluate', {
      expression: expression,
      awaitPromise: awaitPromise !== false,
      returnByValue: true
    });
    if (r.exceptionDetails) {
      throw new Error('page evaluation failed: ' + JSON.stringify(r.exceptionDetails.exception && r.exceptionDetails.exception.description || r.exceptionDetails.text));
    }
    return r.result ? r.result.value : undefined;
  }

  close() {
    try { this.ws.close(); } catch (e) { /* ignore */ }
  }
}

/* ------------------------------ audit core ----------------------------- */

async function waitFor(cdp, expression, timeoutMs, label) {
  const t0 = Date.now();
  for (;;) {
    const v = await cdp.evaluate(expression, true);
    if (v) { return v; }
    if (Date.now() - t0 > timeoutMs) { throw new Error('timed out waiting for ' + label); }
    await sleep(60);
  }
}

async function auditViewport(chromePort, vp, checks) {
  // Fresh tab per viewport so counters are cleanly scoped.
  const created = await (await fetch(
    'http://127.0.0.1:' + chromePort + '/json/new?about:blank', { method: 'PUT' })).json();
  const cdp = new Cdp(created.webSocketDebuggerUrl);
  const state = {
    consoleErrors: [],
    pageErrors: [],
    external: []
  };
  try {
    await cdp.send('Page.enable');
    await cdp.send('Runtime.enable');
    await cdp.send('Network.enable');
    cdp.on('Runtime.consoleAPICalled', (p) => {
      if (p.type === 'error') {
        state.consoleErrors.push((p.args || []).map(a => a.value || a.description || '').join(' '));
      }
    });
    cdp.on('Runtime.exceptionThrown', (p) => {
      const d = p.exceptionDetails || {};
      state.pageErrors.push(d.text + ' ' + (d.exception && d.exception.description || ''));
    });
    cdp.on('Network.requestWillBeSent', (p) => {
      const u = p.request && p.request.url || '';
      if (/^https?:|^wss?:|^ftp:/i.test(u)) { state.external.push(u); }
    });
    await cdp.send('Emulation.setDeviceMetricsOverride', {
      width: vp.width, height: vp.height, deviceScaleFactor: 1, mobile: vp.mobile
    });

    const loaded = new Promise((resolve) => { cdp.on('Page.loadEventFired', resolve); });
    await cdp.send('Page.navigate', { url: PAGE_URL });
    await Promise.race([loaded, sleep(15000)]);
    await waitFor(cdp, 'document.readyState === "complete"', 10000, 'document complete');
    checks.push(['page loads over file://', true]);

    // --- workbench presence
    const controls = await cdp.evaluate(
      '!!document.getElementById("labCard") && !!document.getElementById("compareBtn")');
    checks.push(['workbench controls present', controls === true]);

    // --- comparator run
    await cdp.evaluate(
      'document.getElementById("hpSeqInput").value = "HPPHPHPHPHH";' +
      'document.getElementById("dimSel").value = "2";' +
      'document.getElementById("seedInput").value = "7";' +
      'document.getElementById("compareBtn").click(); true', true);
    await waitFor(cdp,
      'document.querySelectorAll("#labTableWrap tr[data-i]").length >= 5', 30000, 'comparator table');
    checks.push(['comparator table lists five algorithms', true]);

    const exactRow = await cdp.evaluate(
      'Array.from(document.querySelectorAll("#labTableWrap td"))' +
      '.some(td => /ground truth/.test(td.textContent))', true);
    checks.push(['exact enumeration row marked as ground truth', exactRow === true]);

    const status = await cdp.evaluate(
      'document.getElementById("labStatus").textContent', true);
    checks.push(['status reports exact ground state', /exact ground state/.test(status)]);

    const funnelPts = await cdp.evaluate(
      'document.getElementById("funnelSvg").querySelectorAll("circle").length', true);
    checks.push(['folding funnel SVG renders points (' + funnelPts + ')', funnelPts > 40]);

    const cmapCells = await cdp.evaluate(
      'document.getElementById("cmapSvg").querySelectorAll("rect").length', true);
    checks.push(['contact map SVG renders LxL cells (' + cmapCells + ')', cmapCells >= 121]);

    // --- view on stage bridge
    await cdp.evaluate(
      'document.querySelector("#labTableWrap tr[data-i] button").click(); true', true);
    await waitFor(cdp,
      '/HP fold from the comparator/.test(document.getElementById("progressText").textContent)',
      8000, 'stage bridge progress text');
    checks.push(['view on stage bridge updates the 3D viewer', true]);

    // --- mutation scan
    await cdp.evaluate('document.getElementById("mutBtn").click(); true', true);
    await waitFor(cdp,
      'document.getElementById("mutSvg").querySelectorAll("rect").length > 5', 15000, 'mutation bars');
    checks.push(['mutation scan renders per residue bars', true]);

    const mutInfo = await cdp.evaluate(
      'document.getElementById("mutInfo").textContent', true);
    checks.push(['mutation summary reports destabilizing count', /destabilizing/.test(mutInfo)]);

    // --- audit: checksum, verify, tamper, bookmark
    const audit = await cdp.evaluate(
      'JSON.parse(document.getElementById("auditPreview").textContent)', true);
    checks.push(['audit JSON present with sha256 checksum',
      !!audit && /^[0-9a-f]{64}$/.test(audit.checksum || '')]);

    await cdp.evaluate('document.getElementById("verifyBtn").click(); true', true);
    const verdict = await cdp.evaluate(
      'document.getElementById("auditVerdict").textContent', true);
    checks.push(['checksum verification reports intact payload', /checksum verified/.test(verdict)]);

    await cdp.evaluate('document.getElementById("tamperBtn").click(); true', true);
    const tamper = await cdp.evaluate(
      'document.getElementById("auditVerdict").textContent', true);
    checks.push(['tamper check detects the altered digest', /tamper detected/.test(tamper)]);

    await cdp.evaluate('document.getElementById("bookmarkBtn").click(); true', true);
    const hash = await cdp.evaluate('location.hash', true);
    checks.push(['bookmark writes a #fold= URL hash', /^#fold=HPPHPHPHPHH:2d:7/.test(hash || '')]);

    const hashOk = await cdp.evaluate(
      '(function(){var b=window.foldlLab.decodeBookmark(location.hash);' +
      'return b && b.sequence==="HPPHPHPHPHH" && b.dim===2 && b.seed===7;})()', true);
    checks.push(['bookmark decodes back to the run parameters', hashOk === true]);

    // --- self tests
    await cdp.evaluate('document.getElementById("selfTestBtn").click(); true', true);
    const stText = await waitFor(cdp,
      'var t=document.getElementById("selfTestInfo").textContent;' +
      '/self checks passed/.test(t) ? t : false; var t;', 15000, 'self test summary');
    const stMatch = /(\d+)\/(\d+) self checks passed/.exec(stText);
    checks.push(['in page self tests all pass (' + stText + ')',
      !!stMatch && stMatch[1] === stMatch[2] && Number(stMatch[2]) >= 12]);

    const failLis = await cdp.evaluate(
      'document.querySelectorAll("#selfTestList li.fail").length', true);
    checks.push(['self test list shows zero failures', failLis === 0]);

    // --- layout + network hygiene (checked last so late errors count too)
    const overflow = await cdp.evaluate(
      'document.documentElement.scrollWidth - document.documentElement.clientWidth', true);
    checks.push(['no horizontal overflow (scrollWidth surplus ' + overflow + 'px)', overflow <= 1]);

    checks.push(['zero console errors', state.consoleErrors.length === 0,
      state.consoleErrors.slice(0, 3).join(' | ')]);
    checks.push(['zero uncaught page errors', state.pageErrors.length === 0,
      state.pageErrors.slice(0, 3).join(' | ')]);
    checks.push(['zero external network requests', state.external.length === 0,
      state.external.slice(0, 3).join(' | ')]);
  } finally {
    cdp.close();
    try {
      await fetch('http://127.0.0.1:' + chromePort + '/json/close/' + created.id);
    } catch (e) { /* tab may already be gone */ }
  }
}

async function main() {
  if (typeof WebSocket === 'undefined') {
    throw new Error('global WebSocket is missing: run the audit on Node 22 or newer');
  }
  const chrome = findChrome();
  const profile = fs.mkdtempSync(path.join(os.tmpdir(), 'foldl-audit-'));
  const args = [
    '--headless=new',
    '--remote-debugging-port=0',
    '--user-data-dir=' + profile,
    '--no-first-run',
    '--no-default-browser-check',
    '--disable-extensions',
    '--disable-component-update',
    '--mute-audio',
    '--force-device-scale-factor=1',
    'about:blank'
  ];
  const proc = spawn(chrome, args, { stdio: ['ignore', 'ignore', 'pipe'] });
  let stderrBuf = '';
  proc.stderr.on('data', (d) => { stderrBuf += d; });

  try {
    // Chrome with port 0 picks a free port and records it in DevToolsActivePort.
    const portFile = path.join(profile, 'DevToolsActivePort');
    let port = 0;
    for (let i = 0; i < 200 && !port; i++) {
      if (fs.existsSync(portFile)) {
        port = parseInt(fs.readFileSync(portFile, 'utf8').split('\n')[0], 10);
      }
      await sleep(50);
    }
    if (!port) { throw new Error('Chrome did not publish DevToolsActivePort\n' + stderrBuf.slice(-800)); }

    let pass = 0, fail = 0;
    for (const vp of VIEWPORTS) {
      const checks = [];
      try {
        await auditViewport(port, vp, checks);
      } catch (err) {
        checks.push(['viewport run completed', false, String(err.message || err)]);
      }
      console.log('\n== ' + vp.name + ' ==');
      for (const [label, ok, extra] of checks) {
        if (ok) { pass++; console.log('  PASS ' + label); }
        else { fail++; console.log('  FAIL ' + label + (extra ? ' :: ' + extra : '')); }
      }
    }
    console.log('\nheadless Chromium checks: ' + pass + ' passed, ' + fail + ' failed, ' +
      (pass + fail) + ' total across ' + VIEWPORTS.length + ' viewports');
    if (fail > 0) { process.exitCode = 1; }
  } finally {
    try { proc.kill('SIGKILL'); } catch (e) { /* already dead */ }
    try { fs.rmSync(profile, { recursive: true, force: true }); } catch (e) { /* ignore */ }
  }
}

if (require.main === module) {
  main().catch((err) => {
    console.error('browser audit crashed: ' + (err && err.stack || err));
    process.exitCode = 2;
  });
}
