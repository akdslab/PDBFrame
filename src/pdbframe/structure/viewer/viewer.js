// Browser side of pf.view(). widget.py wraps this file in a function, after
// `var D = {...}` - the data it prepared (structures, regions, options).
//
// Residues are keyed 'model:chain:index': the structure's position in D.models
// (0 = the leading frame, then the compared ones), the internal one-letter chain
// ID, and the residue's position in its chain (also its residue number in the PDB
// text that 3Dmol.js reads). Chains are keyed 'model:chain'.

var root = document.getElementById('pfv-' + D.uid);
root.style.setProperty('--highlight', D.highlight);
root.style.setProperty('--pick', D.pick);

// ---- load 3Dmol.js once per page, ignoring require.js, as py3Dmol does
if (typeof $3Dmolpromise === 'undefined' || $3Dmolpromise === null) {
  $3Dmolpromise = new Promise(function(resolve) {
    var savedexports, savedmodule;
    if (typeof exports !== 'undefined') savedexports = exports; else exports = {};
    if (typeof module !== 'undefined') savedmodule = module; else module = {};
    var tag = document.createElement('script');
    tag.src = D.js_url;
    tag.onload = function() { exports = savedexports; module = savedmodule; resolve(); };
    document.head.appendChild(tag);
  });
}

// ---- small helpers
function el(tag, cls, text) {
  var e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}
// black or white text, whichever reads better on a background color
function textOn(color) {
  var probe = el('div'); probe.style.color = color;
  document.body.appendChild(probe);
  var rgb = getComputedStyle(probe).color.match(/\d+/g).map(Number);
  probe.remove();
  return (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) / 255 < 0.5 ? 'white' : 'black';
}
function isAdd(e) { return e && (e.ctrlKey || e.metaKey); }   // Ctrl / Cmd adds a pick
var MANY = D.models.length > 1;
if (MANY) root.classList.add('many');

// ---- residue lookups
var info = {}, cells = {}, order = {}, rows = {}, chainLabels = {}, chainOf = {};
var regionsOf = {}, baseColor = {}, ownColor = {}, valueColor = {}, links = D.links || {};
D.regions.forEach(function(rg) {
  rg.keys.forEach(function(k) { (regionsOf[k] = regionsOf[k] || []).push(rg); });
});
D.models.forEach(function(M, m) {
  M.chains.forEach(function(ch) {
    var c = m + ':' + ch.id;
    chainOf[c] = {model: m, ch: ch};
    order[c] = [];
    ch.residues.forEach(function(r, n) {
      var k = c + ':' + r[0];
      info[k] = {model: m, chainKey: c, chain: ch.label, pos: n, number: r[1], name: r[2],
                 aa: r[3], value: r[5]};
      order[c].push(k);
      valueColor[k] = r[4];
      ownColor[k] = r[4] || ch.color;      // value color, else chain color
    });
  });
});

// ---- state
var viewer = null, models3d = [], startView = null;
var current = null;       // hovered residue
var selected = [];        // selected range (keys, one chain)
var drag = null;          // range being dragged
var picks = [];           // picked items: {type: 'res', key} or {type: 'chain', key}
var pairs = [];           // distances between picks
var hidden = {};          // chains hidden in 3D
var offModels = {};       // structures hidden in 3D
var offRegions = {};      // regions switched off in the legend
D.hidden_regions.forEach(function(name) { offRegions[name] = true; });
var surfaceOn = D.surface;

// a residue's color: its first switched-on region, else its value, else its chain
function recolor() {
  Object.keys(info).forEach(function(k) {
    var on = (regionsOf[k] || []).filter(function(g) { return !offRegions[g.name]; });
    baseColor[k] = on.length ? on[0].color : ownColor[k];
    var cell = cells[k];
    if (!cell) return;
    // the sequence shades regions and values only; chain colors stay on the labels
    var bg = on.length ? on[0].color : valueColor[k];
    cell.style.background = bg || ''; cell.style.color = bg ? textOn(bg) : '';
  });
}

// ---- toolbar
var bar = root.querySelector('.bar');
function chainName(c) {
  var x = chainOf[c];
  return (MANY ? D.models[x.model].name + ' ' : '') + x.ch.label;
}
// a dropdown checklist with "All"; items are [key, label, group heading or null]
function checklist(title, items, isOn, apply, describe) {
  var det = el('details', 'dropdown'), sum = el('summary'), box = el('div', 'menu');
  det.appendChild(sum); det.appendChild(box);
  var checks = {}, lastGroup = null;
  function check(text, onChange) {
    var lab = el('label'), b = el('input');
    b.type = 'checkbox';
    b.addEventListener('change', function() { onChange(b.checked); refresh(); });
    lab.appendChild(b); lab.appendChild(document.createTextNode(' ' + text));
    box.appendChild(lab);
    return b;
  }
  function refresh() {
    var keys = Object.keys(checks), on = keys.filter(function(k) { return checks[k].checked; });
    all.checked = on.length === keys.length;
    all.indeterminate = on.length > 0 && on.length < keys.length;
    sum.textContent = title + ': ' + (on.length === keys.length ? 'all' : on.length === 0 ? 'none' :
      on.length <= 3 ? on.map(describe).join(', ') : on.length + ' of ' + keys.length) + ' ▾';
  }
  var all = check('All', function(v) {
    Object.keys(checks).forEach(function(k) { checks[k].checked = v; apply(k, v); });
  });
  items.forEach(function(it) {
    if (it[2] && it[2] !== lastGroup) { box.appendChild(el('div', 'menu-head', it[2])); lastGroup = it[2]; }
    checks[it[0]] = check(it[1], function(v) { apply(it[0], v); });
    checks[it[0]].checked = isOn(it[0]);
  });
  refresh();
  bar.appendChild(det);
  return {set: function(k, v) { checks[k].checked = v; refresh(); }};
}
function button(text, title, onClick, parent) {
  var b = el('button', null, text); b.title = title;
  b.addEventListener('click', onClick); (parent || bar).appendChild(b);
  return b;
}
var chainItems = [];
D.models.forEach(function(M, m) {
  M.chains.forEach(function(ch) { chainItems.push([m + ':' + ch.id, 'Chain ' + ch.label, MANY ? M.name : null]); });
});

checklist('Sequences', chainItems, function(c) { return chainOf[c].ch.shown; }, function(c, on) {
  showRow(c, on);
  updatePanel();
}, chainName);
checklist('3D chains', chainItems, function() { return true; }, function(c, on) {
  if (on) delete hidden[c]; else hidden[c] = true;
  redrawSurface(); draw();
}, chainName);

// structures: show / hide each, or flip through them one at a time
var structMenu = null, oneAtATime = false, currentModel = 0;
if (MANY) {
  structMenu = checklist('Structures', D.models.map(function(M, m) { return [String(m), M.name, null]; }),
    function() { return true; },
    function(m, on) { if (on) delete offModels[m]; else offModels[m] = true; redrawSurface(); draw(); },
    function(m) { return D.models[+m].name; });
  var flip = el('span', 'flip');
  var oneLab = el('label', 'check'), oneBox = el('input');
  oneBox.type = 'checkbox';
  oneBox.addEventListener('change', function() { setOneAtATime(oneBox.checked); });
  oneLab.appendChild(oneBox); oneLab.appendChild(document.createTextNode(' One at a time'));
  oneLab.title = 'Show one structure; step with ◀ ▶ or the arrow keys';
  flip.appendChild(oneLab);
  button('◀', 'Previous structure (left arrow)', function() { step(-1); }, flip);
  button('▶', 'Next structure (right arrow)', function() { step(1); }, flip);
  bar.appendChild(flip);
  root.tabIndex = 0;         // so the arrow keys reach the widget once it is clicked
  root.addEventListener('keydown', function(e) {
    if (!oneAtATime || (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight')) return;
    e.preventDefault(); step(e.key === 'ArrowRight' ? 1 : -1);
  });
}
function showOnly(m) {
  currentModel = m;
  D.models.forEach(function(M, i) {
    var on = i === m;
    if (on) delete offModels[i]; else offModels[i] = true;
    structMenu.set(String(i), on);
  });
  root.querySelector('.info').textContent = 'Showing ' + D.models[m].name +
    ' (' + (m + 1) + ' of ' + D.models.length + ')';
  redrawSurface(); draw();
}
function setOneAtATime(on) {
  oneAtATime = on;
  if (on) { showOnly(currentModel); root.focus(); }
  else {
    D.models.forEach(function(M, i) { delete offModels[i]; structMenu.set(String(i), true); });
    redrawSurface(); draw();
  }
}
function step(d) {
  if (!oneAtATime) { oneBox.checked = true; oneAtATime = true; }
  showOnly((currentModel + d + D.models.length) % D.models.length);
}

var surfLab = el('label', 'check'), surfBox = el('input');
surfBox.type = 'checkbox'; surfBox.checked = surfaceOn;
surfBox.addEventListener('change', function() { surfaceOn = surfBox.checked; redrawSurface(); draw(); });
surfLab.appendChild(surfBox); surfLab.appendChild(document.createTextNode(' Surface'));
bar.appendChild(surfLab);

// arrows and spheres: one checkbox per group, shared, so edges and their nodes switch together
var shapesOff = {}, shapeInfo = null;
(D.arrows || []).concat(D.spheres || []).map(function(a) { return a.group; })
  .filter(function(g, i, all) { return all.indexOf(g) === i; })
  .forEach(function(g) {
    var lab = el('label', 'check'), box = el('input');
    box.type = 'checkbox'; box.checked = true;
    box.addEventListener('change', function() { shapesOff[g] = !box.checked; draw(); });
    lab.appendChild(box); lab.appendChild(document.createTextNode(' ' + g));
    bar.appendChild(lab);
  });
var arrowsOff = shapesOff;   // the old name, for anything reading it
function point(v) { return {x: v[0], y: v[1], z: v[2]}; }
function middle(a) {
  return {x: (a.start[0] + a.end[0]) / 2, y: (a.start[1] + a.end[1]) / 2,
          z: (a.start[2] + a.end[2]) / 2};
}
// clicking a shape that carries `info` writes it in the info line; clicking it again clears
function clicker(shape, where) {
  if (!shape.info) return {};
  return {clickable: true, callback: function() {
    shapeInfo = (shapeInfo && shapeInfo.text === shape.info) ? null
              : {text: shape.info, color: shape.color, at: where};
    draw();
  }};
}
function drawArrows() {
  (D.arrows || []).forEach(function(a) {
    if (shapesOff[a.group]) return;
    var spec = {start: point(a.start), end: point(a.end), radius: a.radius, color: a.color};
    if (a.head === false) {
      viewer.addCylinder(Object.assign(spec, clicker(a, middle(a))));
    } else {
      viewer.addArrow(Object.assign(spec, {radiusRatio: 2.2, mid: 0.85},
                                    clicker(a, middle(a))));
    }
    if (a.label) {
      viewer.addLabel(a.label, {position: point(a.end), fontSize: 12, fontColor: 'white',
                                backgroundColor: a.color, backgroundOpacity: 0.8});
    }
  });
  (D.spheres || []).forEach(function(sp) {
    if (shapesOff[sp.group]) return;
    viewer.addSphere(Object.assign({center: point(sp.center), radius: sp.radius,
                                    color: sp.color}, clicker(sp, point(sp.center))));
    if (sp.label) {
      viewer.addLabel(sp.label, {position: point(sp.center), fontSize: 11, fontColor: 'white',
                                 backgroundColor: sp.color, backgroundOpacity: 0.8});
    }
  });
  if (shapeInfo) {
    viewer.addLabel(shapeInfo.text, {position: shapeInfo.at, fontSize: 12, fontColor: 'white',
                                     backgroundColor: 'black', backgroundOpacity: 0.75});
  }
}
button('Reset view', 'Back to the starting orientation and zoom', function() {
  if (viewer && startView) { viewer.setView(startView); viewer.render(); }
});

// saved orientations: turn the camera so the given frame is seen down -z with y up
function quaternion(m) {   // rotation matrix -> (x, y, z, w), for v' = m v
  var t = m[0][0] + m[1][1] + m[2][2], s;
  if (t > 0) {
    s = Math.sqrt(t + 1) * 2;
    return [(m[2][1] - m[1][2]) / s, (m[0][2] - m[2][0]) / s, (m[1][0] - m[0][1]) / s, s / 4];
  }
  if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
    s = Math.sqrt(1 + m[0][0] - m[1][1] - m[2][2]) * 2;
    return [s / 4, (m[0][1] + m[1][0]) / s, (m[0][2] + m[2][0]) / s, (m[2][1] - m[1][2]) / s];
  }
  if (m[1][1] > m[2][2]) {
    s = Math.sqrt(1 + m[1][1] - m[0][0] - m[2][2]) * 2;
    return [(m[0][1] + m[1][0]) / s, s / 4, (m[1][2] + m[2][1]) / s, (m[0][2] - m[2][0]) / s];
  }
  s = Math.sqrt(1 + m[2][2] - m[0][0] - m[1][1]) * 2;
  return [(m[0][2] + m[2][0]) / s, (m[1][2] + m[2][1]) / s, s / 4, (m[1][0] - m[0][1]) / s];
}
function orient(R, t) {
  if (!viewer) return;
  var v = viewer.getView(), q = quaternion(R);
  v[4] = q[0]; v[5] = q[1]; v[6] = q[2]; v[7] = q[3];
  if (D.orientation_center === 'frame') {
    // the point the transform maps to the origin, -R^T t, goes to the middle
    for (var i = 0; i < 3; i++) v[i] = R[0][i] * t[0] + R[1][i] * t[1] + R[2][i] * t[2];
  } else if (startView) {
    v[0] = startView[0]; v[1] = startView[1]; v[2] = startView[2];
  }
  viewer.setView(v); viewer.render();
}
Object.keys(D.orientations).forEach(function(name) {
  var o = D.orientations[name];
  button(name, 'Turn the view to this orientation', function() { orient(o[0], o[1]); });
});

button('Copy as code', 'Copy the selection or the picks as a pf.select(...) call', copyAsCode);
button('Save PNG', 'Download the 3D view as an image', function() {
  if (!viewer) return;
  var a = el('a'); a.href = viewer.pngURI(); a.download = (D.models[0].name || 'structure') + '.png';
  document.body.appendChild(a); a.click(); a.remove();
});
button('Clear', 'Clear the selection and the picks', function() {
  picks.slice().forEach(togglePick);
  setSelection([]);
});

// help
var help = el('details', 'dropdown help'), helpSum = el('summary', null, '?');
helpSum.title = 'How to use this viewer';
var helpBox = el('div', 'menu helpbox');
helpBox.innerHTML =
  '<b>Sequence panel</b> (choose chains with <i>Sequences</i>)<ul>' +
  '<li>Hover a residue to light it up in 3D (and in the 3D view to find it here).</li>' +
  '<li>Drag across residues to select a range. Click a chain label to select the whole chain. ' +
  'Click a selected residue again to clear.</li>' +
  '<li>Regions: click a name in the legend to switch it off or on; <kbd>Shift</kbd>-click to show ' +
  'only it, and again to show all.</li></ul>' +
  '<b>Measuring</b><ul>' +
  '<li><kbd>Ctrl</kbd>-click (<kbd>⌘</kbd>-click on a Mac) residues or chain labels, in either panel, to pick them.</li>' +
  '<li>Every pair of picks gets its closest atom-atom distance (and Cα–Cα for two residues), ' +
  'listed below and drawn as a dashed line.</li></ul>' +
  '<b>Moving in 3D</b><ul>' +
  '<li>Rotate: drag.</li>' +
  '<li>Move (pan): <kbd>Ctrl</kbd>+drag, or drag with the middle mouse button.</li>' +
  '<li>Zoom: scroll, or <kbd>Shift</kbd>+drag.</li>' +
  '<li><i>Reset view</i> goes back to the start; named buttons turn to saved orientations.</li></ul>' +
  (MANY ? '<b>Several structures</b><ul>' +
  '<li>Hovering or selecting a residue lights up its partners in every structure; their sequences ' +
  'are aligned under the first structure’s.</li>' +
  '<li><i>Structures</i>: show or hide each. <i>One at a time</i>: flip through them with ' +
  '◀ ▶ or the arrow keys (click the viewer first).</li></ul>' : '') +
  '<b>Other controls</b><ul>' +
  '<li><i>3D chains</i>: hide chains in the 3D view. <i>Surface</i>: molecular surface.</li>' +
  '<li><i>Copy as code</i>: the selection or picks as a <code>pf.select(...)</code> call.</li>' +
  '<li><i>Save PNG</i>: an image of the 3D view. <i>Clear</i>: drop the selection and picks.</li></ul>';
help.appendChild(helpSum); help.appendChild(helpBox); bar.appendChild(help);

document.addEventListener('click', function(e) {
  root.querySelectorAll('details.dropdown').forEach(function(d) { if (!d.contains(e.target)) d.open = false; });
});

// ---- sequence panel: the lead's chains, each with its partners aligned underneath.
// One structure: each chain's row wraps. Several: blocks of columns, every row of a
// group stacked in each block, so paired residues line up as in an alignment.
var tracks = root.querySelector('.tracks');
var groupBoxes = [];
var BLOCK = Math.max(20, Math.floor((D.width - 150) / 12));   // columns per block
function chainLabel(c) {
  var x = chainOf[c], ch = x.ch, m = x.model;
  var lab = el('div', 'chain', (MANY && m > 0 ? D.models[m].name + ' ' : '') + 'Chain ' + ch.label);
  if (MANY) lab.classList.add('wide');
  lab.style.background = ch.color; lab.style.color = textOn(ch.color);
  lab.title = (MANY ? D.models[m].name + ', chain ' + ch.label + '. ' : '') +
    'Click: select the chain. Ctrl/Cmd-click: pick it.';
  lab.addEventListener('click', function(e) {
    if (isAdd(e)) togglePick({type: 'chain', key: c});
    else setSelection(order[c].slice());
  });
  (chainLabels[c] = chainLabels[c] || []).push(lab);
  return lab;
}
function residueCell(c, idx, counter) {
  var cell = el('div', 'res');
  if (idx === null) {             // a gap, so paired residues line up
    cell.appendChild(el('div', 'num', ''));
    cell.appendChild(el('div', 'aa gap', '-'));
    return cell;
  }
  var k = c + ':' + idx, r = info[k];
  cell.appendChild(el('div', 'num', counter.n % 10 === 0 ? r.number : ''));
  counter.n++;
  var aa = el('div', 'aa', r.aa);
  aa.addEventListener('mousedown', function(e) {
    e.preventDefault();
    if (isAdd(e)) togglePick({type: 'res', key: k}); else startDrag(k);
  });
  aa.addEventListener('mouseenter', function() { if (drag) extendDrag(k); hover(k, false); });
  aa.addEventListener('mouseleave', function() { hover(null); });
  cell.appendChild(aa);
  cells[k] = aa;
  return cell;
}
D.layout.forEach(function(group) {
  var box = el('div', 'group'); tracks.appendChild(box);
  groupBoxes.push({box: box, chains: group.rows.map(function(r) { return r.chain; })});
  var counters = group.rows.map(function() { return {n: 0}; });
  var width = MANY ? BLOCK : Infinity;
  var ncols = Math.max.apply(null, group.rows.map(function(r) { return r.cells.length; }));
  for (var start = 0; start < ncols; start += width) {
    var block = el('div', 'block'); box.appendChild(block);
    group.rows.forEach(function(row, i) {
      var c = row.chain, part = row.cells.slice(start, start + width);
      if (!part.length) return;
      var line = el('div', 'row');
      line.appendChild(chainLabel(c));
      var seq = el('div', 'seq');
      part.forEach(function(idx) { seq.appendChild(residueCell(c, idx, counters[i])); });
      line.appendChild(seq); block.appendChild(line);
      (rows[c] = rows[c] || []).push(line);
    });
  }
  group.rows.forEach(function(row) { showRow(row.chain, chainOf[row.chain].ch.shown); });
});
function showRow(c, on) { rows[c].forEach(function(line) { line.style.display = on ? '' : 'none'; }); }
function rowShown(c) { return rows[c] && rows[c][0].style.display !== 'none'; }
function updatePanel() {
  // a group's box (and its spacing) shows only while one of its rows does
  groupBoxes.forEach(function(g) { g.box.style.display = g.chains.some(rowShown) ? '' : 'none'; });
}
updatePanel();

// ---- legend: region chips switch regions on and off
var legend = root.querySelector('.legend'), chips = {};
function regionToggle(names, e) {
  var all = D.regions.map(function(g) { return g.name; });
  if (e.shiftKey) {
    // show only these; if that is already the case, show all
    var only = all.every(function(x) { return (names.indexOf(x) >= 0) !== !!offRegions[x]; });
    all.forEach(function(x) { if (only || names.indexOf(x) >= 0) delete offRegions[x]; else offRegions[x] = true; });
  } else {
    var anyOn = names.some(function(x) { return !offRegions[x]; });
    names.forEach(function(x) { if (anyOn) offRegions[x] = true; else delete offRegions[x]; });
  }
  refreshChips(); recolor(); draw();
}
function refreshChips() {
  Object.keys(chips).forEach(function(key) {
    var names = chips[key].names;
    chips[key].el.classList.toggle('off', names.every(function(x) { return offRegions[x]; }));
  });
}
function chip(key, text, color, names) {
  var s = el('span', 'chip toggle', text);
  s.style.background = color; s.style.color = textOn(color);
  s.title = 'Click: switch on / off. Shift-click: show only this, again to show all.';
  s.addEventListener('click', function(e) { regionToggle(names, e); });
  legend.appendChild(s);
  chips[key] = {el: s, names: names};
}
Object.keys(D.region_groups).forEach(function(g) {
  chip('group:' + g, g + ' ▣', '#666666', D.region_groups[g]);
});
D.regions.forEach(function(rg) { chip('region:' + rg.name, rg.name, rg.color, [rg.name]); });
D.value_scales.forEach(function(vs) {
  var scale = el('span', 'scale');
  scale.appendChild(el('span', null, vs.name + ' ' + vs.lo));
  var grad = el('span', 'gradient');
  grad.style.background = 'linear-gradient(to right, ' + vs.colors.join(', ') + ')';
  scale.appendChild(grad); scale.appendChild(el('span', null, String(vs.hi)));
  legend.appendChild(scale);
});
D.models.forEach(function(M, m) {
  if (m === 0 || !M.chains.length) return;
  var c = el('span', 'chip', M.name), color = M.chains[0].color;
  c.style.background = color; c.style.color = textOn(color);
  legend.appendChild(c);
});
refreshChips(); recolor();

// ---- selecting a range: drag across the sequence (one chain)
function group(k) {
  // a residue plus its partners in every structure, through the leading structure
  var lead = info[k].model === 0 ? [k] : (links[k] || []);
  var out = [k];
  lead.forEach(function(j) {
    if (out.indexOf(j) < 0) out.push(j);
    (links[j] || []).forEach(function(p) { if (out.indexOf(p) < 0) out.push(p); });
  });
  return out;
}
function withLinks(keys) {
  var seen = {}, out = [];
  keys.forEach(function(k) { group(k).forEach(function(j) { if (!seen[j]) { seen[j] = 1; out.push(j); } }); });
  return out;
}
var selShown = {};   // selected residues plus their partners in the other structures
function setSelection(keys) {
  Object.keys(selShown).forEach(function(k) { if (cells[k]) cells[k].classList.remove('sel'); });
  selected = keys; selShown = {};
  withLinks(selected).forEach(function(k) { selShown[k] = true; if (cells[k]) cells[k].classList.add('sel'); });
  describe(); draw();
}
function startDrag(k) {
  // clicking the only selected residue clears it; otherwise start a new range
  var again = selected.length === 1 && selected[0] === k;
  setSelection(again ? [] : [k]);
  drag = again ? null : {start: k};
}
function extendDrag(k) {
  var a = info[drag.start], b = info[k];
  if (a.chainKey !== b.chainKey) return;
  var lo = Math.min(a.pos, b.pos), hi = Math.max(a.pos, b.pos);
  setSelection(order[a.chainKey].slice(lo, hi + 1));
}
document.addEventListener('mouseup', function() { drag = null; });

// ---- picking residues / chains (Ctrl/Cmd-click) and measuring between them
function pickId(item) { return item.type + item.key; }
function togglePick(item) {
  var i = picks.map(pickId).indexOf(pickId(item));
  var els = item.type === 'res' ? [cells[item.key]] : chainLabels[item.key];
  var on = i < 0;
  if (on) picks.push(item); else picks.splice(i, 1);
  (els || []).forEach(function(e) { if (e) e.classList.toggle('picked', on); });
  measure(); draw();
}
function sel(k) {
  var p = k.split(':');
  return {model: models3d[+p[0]], chain: p[1], resi: parseInt(p[2])};
}
function chainSel(c) { var p = c.split(':'); return {model: models3d[+p[0]], chain: p[1]}; }
function atomsOf(item) { return viewer.selectedAtoms(item.type === 'res' ? sel(item.key) : chainSel(item.key)); }
function label(k) { var r = info[k]; return r.name + ' ' + r.number; }
function itemName(item) {
  if (item.type === 'chain') return 'Chain ' + chainName(item.key);
  var r = info[item.key];
  return (MANY ? D.models[r.model].name + ' ' : '') + r.chain + ' · ' + label(item.key);
}
function atomName(a) { return label(a.model + ':' + a.chain + ':' + a.resi) + ' ' + a.atom; }
function closest(A, B) {
  var best = Infinity, ba = null, bb = null;
  for (var i = 0; i < A.length; i++) {
    var a = A[i];
    for (var j = 0; j < B.length; j++) {
      var b = B[j], dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z, d = dx * dx + dy * dy + dz * dz;
      if (d < best) { best = d; ba = a; bb = b; }
    }
  }
  return {d: Math.sqrt(best), a: ba, b: bb};
}
function dist(a, b) { return Math.sqrt(Math.pow(a.x - b.x, 2) + Math.pow(a.y - b.y, 2) + Math.pow(a.z - b.z, 2)); }
function regionNames(keys) {
  // every region a residue is in, switched on or not
  var regs = {};
  keys.forEach(function(k) { (regionsOf[k] || []).forEach(function(g) { regs[g.name] = 1; }); });
  return Object.keys(regs).join(', ') || '–';
}
function measure() {
  var box = root.querySelector('.picks');
  pairs = [];
  if (!picks.length) { box.innerHTML = ''; return; }
  if (!viewer) { box.textContent = D.sequence_only ? 'Show 3D to measure.' : 'Loading...'; return; }
  var atoms = picks.map(atomsOf);
  var html = '<table><tr><th>Picked</th><th>Regions</th><th>Atoms</th></tr>';
  picks.forEach(function(item, i) {
    var keys = item.type === 'res' ? [item.key] : order[item.key];
    html += '<tr><td>' + itemName(item) + '</td><td>' + regionNames(keys) + '</td><td>' +
            atoms[i].length + '</td></tr>';
  });
  html += '</table>';
  if (picks.length > 1) {
    html += '<table><tr><th>Between</th><th>Closest atoms</th><th>Cα–Cα</th></tr>';
    for (var i = 0; i < picks.length; i++) {
      for (var j = i + 1; j < picks.length; j++) {
        var c = closest(atoms[i], atoms[j]), caText = '–';
        if (picks[i].type === 'res' && picks[j].type === 'res') {
          var ca1 = atoms[i].filter(function(a) { return a.atom === 'CA'; })[0];
          var ca2 = atoms[j].filter(function(a) { return a.atom === 'CA'; })[0];
          if (ca1 && ca2) caText = dist(ca1, ca2).toFixed(1) + ' Å';
        }
        pairs.push(c);
        html += '<tr><td>' + itemName(picks[i]) + ' ↔ ' + itemName(picks[j]) + '</td><td>' +
                c.d.toFixed(1) + ' Å (' + atomName(c.a) + ' – ' + atomName(c.b) + ')</td><td>' +
                caText + '</td></tr>';
      }
    }
    html += '</table>';
  }
  box.innerHTML = html;
}

// ---- copy the selection / picks as a pf.select(...) call
function copyAsCode() {
  var byChain = {}, whole = {};
  if (picks.length) {
    picks.forEach(function(item) {
      if (item.type === 'chain') { whole[item.key] = true; return; }
      var c = info[item.key].chainKey;
      (byChain[c] = byChain[c] || []).push(item.key);
    });
  } else if (selected.length) {
    byChain[info[selected[0]].chainKey] = selected.slice();
  }
  var parts = {};
  function add(m, args) { (parts[m] = parts[m] || []).push(args); }
  Object.keys(whole).forEach(function(c) { add(chainOf[c].model, "chains='" + chainOf[c].ch.label + "'"); });
  Object.keys(byChain).forEach(function(c) {
    if (whole[c]) return;
    var nums = byChain[c].map(function(k) { return parseInt(info[k].number); });
    var range = !picks.length && nums.length > 1;   // a dragged range is contiguous
    add(chainOf[c].model, "chains='" + chainOf[c].ch.label + "', residues=" +
        (range ? '(' + nums[0] + ', ' + nums[nums.length - 1] + ')' : '[' + nums.join(', ') + ']'));
  });
  // the lead is `pf`; with one compared structure it is `other`, with several `others[i]`
  function varName(m) {
    m = +m;
    if (m === 0) return 'pf';
    return D.models.length === 2 ? 'other' : 'others[' + (m - 1) + ']';
  }
  var code = Object.keys(parts).map(function(m) {
    return parts[m].map(function(args) { return varName(m) + '.select(' + args + ')'; }).join(' + ');
  }).join('\n');
  var infoBox = root.querySelector('.info');
  if (!code) { infoBox.textContent = 'Nothing selected or picked to copy.'; return; }
  function show(prefix) { infoBox.textContent = prefix + code; }
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(code).then(function() { show('Copied: '); },
                                             function() { show('Copy this: '); });
  } else { show('Copy this: '); }
}

// ---- hover and the info line
function describe() {
  var box = root.querySelector('.info');
  if (current) {
    var r = info[current], regs = (regionsOf[current] || []).map(function(g) { return g.name; });
    var val = (r.value !== null && r.value !== undefined && D.value_scales.length) ?
      ' · ' + (D.value_scales.length === 1 ? D.value_scales[0].name : 'value') + ' ' + r.value : '';
    box.textContent = (MANY ? D.models[r.model].name + ' · ' : '') + 'Chain ' + r.chain +
      ' · ' + label(current) + val + (regs.length ? ' · ' + regs.join(', ') : '');
  } else if (selected.length) {
    var first = selected[0], last = selected[selected.length - 1];
    box.textContent = 'Selected: chain ' + info[first].chain + ' · ' + label(first) +
      (selected.length > 1 ? ' – ' + label(last) + ' (' + selected.length + ' residues)' : '');
  } else {
    box.textContent = 'Hover a residue. Drag across the sequence to select a range. ' +
      'Ctrl/Cmd-click residues or chain labels to measure. See ? for help.';
  }
}
var hovered = [];
function hover(k, fromViewer) {
  if (k === current) return;
  hovered.forEach(function(j) { if (cells[j]) cells[j].classList.remove('hov'); });
  current = k;
  hovered = k ? withLinks([k]) : [];
  hovered.forEach(function(j) { if (cells[j]) cells[j].classList.add('hov'); });
  if (k && fromViewer && cells[k] && rowShown(info[k].chainKey)) {
    cells[k].scrollIntoView({block: 'nearest'});
  }
  describe(); draw();
}

// ---- 3D drawing
function selMany(keys) {
  var groups = {};
  keys.forEach(function(k) {
    var p = k.split(':'), g = p[0] + ':' + p[1];
    (groups[g] = groups[g] || []).push(parseInt(p[2]));
  });
  return Object.keys(groups).map(function(g) {
    var p = g.split(':');
    return {model: models3d[+p[0]], chain: p[1], resi: groups[g], _m: +p[0]};
  });
}
function strip(s) { return {model: s.model, chain: s.chain, resi: s.resi}; }
function visible(m) { return !offModels[m]; }
function draw() {
  if (!viewer) return;
  var pickedChains = {};
  picks.forEach(function(item) { if (item.type === 'chain') pickedChains[item.key] = true; });
  models3d.forEach(function(mod, m) {
    if (!visible(m)) { viewer.setStyle({model: mod}, {}); return; }
    var spec = {};
    // one color function per structure: selection, then picked chain, then region / value / chain
    spec[D.models[m].style || D.style] = {colorfunc: function(a) {
      var k = m + ':' + a.chain + ':' + a.resi;
      if (selShown[k]) return D.highlight;
      if (pickedChains[m + ':' + a.chain]) return D.pick;
      return baseColor[k] || '#999999';
    }, opacity: D.opacity === undefined ? 1 : D.opacity};
    viewer.setStyle({model: mod}, spec);
    if (D.ligands) viewer.addStyle({model: mod, hetflag: true}, {stick: {radius: 0.25}});
  });
  D.regions.forEach(function(rg) {
    if (!rg.extra || offRegions[rg.name]) return;
    selMany(rg.keys).forEach(function(s) {
      if (!visible(s._m) || rg.extra === (D.models[s._m].style || D.style)) return;
      var o = {}; o[rg.extra] = {color: rg.color}; viewer.addStyle(strip(s), o);
    });
  });
  selMany(Object.keys(selShown)).forEach(function(s) {
    if (visible(s._m)) viewer.addStyle(strip(s), {stick: {colorscheme: 'greenCarbon', radius: 0.2}});
  });
  picks.forEach(function(item) {
    if (item.type === 'res' && visible(info[item.key].model)) {
      viewer.addStyle(sel(item.key), {stick: {color: D.pick, radius: 0.25}});
    }
  });
  Object.keys(hidden).forEach(function(c) { viewer.setStyle(chainSel(c), {}); });
  viewer.removeAllLabels();
  viewer.removeAllShapes();
  drawArrows();
  pairs.forEach(function(c) {
    viewer.addCylinder({start: {x: c.a.x, y: c.a.y, z: c.a.z}, end: {x: c.b.x, y: c.b.y, z: c.b.z},
                        radius: 0.07, dashed: true, color: D.pick});
    viewer.addLabel(c.d.toFixed(1) + ' Å',
      {position: {x: (c.a.x + c.b.x) / 2, y: (c.a.y + c.b.y) / 2, z: (c.a.z + c.b.z) / 2},
       fontSize: 11, fontColor: 'black', backgroundColor: D.pick, backgroundOpacity: 0.85});
  });
  if (current) {
    hovered.forEach(function(j) {
      if (!hidden[info[j].chainKey] && visible(info[j].model)) {
        viewer.addStyle(sel(j), {stick: {color: D.highlight, radius: 0.3}});
      }
    });
    var r = info[current];
    if (visible(r.model)) {
      viewer.addLabel((MANY ? D.models[r.model].name + ' ' : '') + 'Chain ' + r.chain + ' ' + label(current),
        {fontSize: 12, fontColor: 'white', backgroundColor: 'black', backgroundOpacity: 0.7},
        sel(current));
    }
  }
  viewer.render();
}
function redrawSurface() {
  if (!viewer) return;
  viewer.removeAllSurfaces();
  if (!surfaceOn) return;
  models3d.forEach(function(mod, m) {
    if (!visible(m)) return;
    var shown = D.models[m].chains.map(function(ch) { return ch.id; })
      .filter(function(id) { return !hidden[m + ':' + id]; });
    if (shown.length) {
      viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.55, color: 'white'},
                        {model: mod, chain: shown, hetflag: false});
    }
  });
}

// ---- x / y / z arrows in the corner, turning with the main view
// The letters are HTML placed at the projected arrow tips: 3Dmol.js labels come out
// tiny when an axis points straight at the camera.
function makeAxes() {
  var box = root.querySelector('.axes');
  var g = $3Dmol.createViewer(box, {backgroundColor: 'white', backgroundAlpha: 0});
  var axes = [['X', {x: 1, y: 0, z: 0}, '#D62728'], ['Y', {x: 0, y: 1, z: 0}, '#2CA02C'],
              ['Z', {x: 0, y: 0, z: 1}, '#1F77B4']];
  var letters = axes.map(function(ax) {
    g.addArrow({start: {x: 0, y: 0, z: 0}, end: ax[1], radius: 0.06, radiusRatio: 2.2, mid: 0.72,
                color: ax[2]});
    var s = el('span', 'axis-letter', ax[0]); s.style.color = ax[2]; box.appendChild(s);
    return s;
  });
  // zoomTo sets the clipping planes; then zoom in and put the origin in the middle
  g.zoomTo(); g.zoom(4.2);
  var start = g.getView(); start[0] = start[1] = start[2] = 0; g.setView(start);
  function sync() {
    var v = viewer.getView(), w = g.getView();
    w[4] = v[4]; w[5] = v[5]; w[6] = v[6]; w[7] = v[7];   // copy the rotation only
    g.setView(w); g.render();
    // modelToScreen gives page coordinates; place letters relative to the origin,
    // which sits in the middle of the box
    var o = g.modelToScreen({x: 0, y: 0, z: 0});
    axes.forEach(function(ax, i) {
      var p = g.modelToScreen({x: ax[1].x * 1.3, y: ax[1].y * 1.3, z: ax[1].z * 1.3});
      letters[i].style.left = (box.clientWidth / 2 + p.x - o.x) + 'px';
      letters[i].style.top = (box.clientHeight / 2 + p.y - o.y) + 'px';
    });
  }
  viewer.setViewChangeCallback(sync);
  sync();
  root.pfv.axes = g;
}

describe();
function start3d() {
  var box = root.querySelector('.view');
  box.style.width = D.width + 'px'; box.style.height = D.height + 'px';
  viewer = $3Dmol.createViewer(box, {backgroundColor: 'white'});
  root.pfv = {viewer: viewer, orient: orient};   // for scripting from the page
  D.models.forEach(function(M) { models3d.push(viewer.addModel(M.pdb, 'pdb')); });
  viewer.setHoverable({}, true,
    function(atom) {
      var k = atom.model + ':' + atom.chain + ':' + atom.resi;
      if (info[k] && visible(info[k].model)) hover(k, true);
    },
    function() { hover(null); });
  viewer.setHoverDuration(40);
  viewer.setClickable({}, true, function(atom, v, event) {
    var k = atom.model + ':' + atom.chain + ':' + atom.resi;
    if (!info[k] || !visible(info[k].model)) return;
    if (isAdd(event)) togglePick({type: 'res', key: k});
    else setSelection(selected.length === 1 && selected[0] === k ? [] : [k]);
  });
  draw();
  viewer.zoomTo();
  startView = viewer.getView();
  if (D.orientation) { orient(D.orientation[0], D.orientation[1]); startView = viewer.getView(); }
  if (D.axes) makeAxes(); else root.querySelector('.axes').remove();
  redrawSurface();
  measure(); draw();
}
if (!D.sequence_only) $3Dmolpromise.then(start3d);
else {
  // sequences only: the 3D view starts when asked for
  root.querySelector('.view').style.display = 'none';
  var show3d = button('Show 3D', 'Open the 3D structure under the sequences', function() {
    root.querySelector('.view').style.display = '';
    show3d.remove();
    $3Dmolpromise.then(start3d);
  });
}
