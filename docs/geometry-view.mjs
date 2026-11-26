import {structure, defaults, pose} from './geometry.mjs';

const el = id => document.getElementById(id);
const specs = [
  ['height','TNF origin above plane',0,30,0.1,'nm'],
  ['tilt','Complex tilt',-85,85,1,'°'],
  ['offset','Domain center from Fv end',0,12,0.1,'nm'],
  ['radius','Domain radius',0.1,8,0.1,'nm'],
  ['bearing','Domain bearing at zero tilt',-180,180,1,'°'],
];
el('geometry-sliders').innerHTML = specs.map(([key,label,min,max,step]) =>
  `<label class="geometry-slider" for="geom-${key}"><span id="geom-${key}-label">${label}</span><output id="geom-${key}-value" for="geom-${key}"></output></label><input aria-labelledby="geom-${key}-label" id="geom-${key}" type="range" min="${min}" max="${max}" step="${step}">`).join('');

function scene(result, radius, showDomain) {
  const {rotate, center} = result;
  const traces = Object.entries(structure.traces).map(([id,rows]) => [id,rows.map(([res,...p]) => [res,...rotate(p)])]);
  const anchor = rotate(structure.end);
  const all = traces.flatMap(([,rows]) => rows.map(([,x,y]) => [x,y]));
  if (showDomain) all.push([center[0]-radius,center[1]-radius],[center[0]+radius,center[1]+radius]);
  all.push([-8,-4],[12,16]);
  const xs = all.map(p=>p[0]), ys=all.map(p=>p[1]);
  const lowX = Math.min(...xs), highX = Math.max(...xs), lowY = Math.min(...ys), highY = Math.max(...ys);
  const scale = Math.min(49, 360/(highX-lowX+1), 230/(highY-lowY+1));
  const x = v => 224+(v-(lowX+highX)/2)*scale;
  const y = v => 143-(v-(lowY+highY)/2)*scale;
  const point = p => `${x(p[0])},${y(p[1])}`;
  const stroke = {A:'#365f82',B:'#6b91ae',C:'#91aabd',F:'#087563',I:'#35987d'};
  const plane = y(0);
  const parts = [`<svg viewBox="0 0 448 286" role="img" aria-label="5YOY side projection at the selected tilt, relative to a hypothetical flat membrane">`,
    `<rect x="0" y="${plane}" width="448" height="${286-plane}" fill="#e5e9dd"/>`,
    `<path d="M0 ${plane}H448" stroke="#566858" stroke-dasharray="5 4"/>`,
    `<text x="12" y="${plane+18}" class="geom-svg-label">Membrane plane</text>`,
    `<path d="M382 268h${2*scale}" stroke="#203c31" stroke-width="2" transform="translate(${-2*scale+46} 0)"/>`,
    '<text x="428" y="258" text-anchor="end" class="geom-svg-label">2 nm</text>'];
  if (showDomain) parts.push(
    `<line x1="${x(anchor[0])}" y1="${y(anchor[1])}" x2="${x(center[0])}" y2="${y(center[1])}" stroke="#774376" stroke-dasharray="3 3"/>`,
    `<circle cx="${x(center[0])}" cy="${y(center[1])}" r="${radius*scale}" fill="${result.envelopeClearance < 0 ? '#9b4832':'#774376'}" fill-opacity=".13" stroke="#774376" stroke-width="1.5" stroke-dasharray="5 3"/>`);
  for (const [id,rows] of traces) {
    let d='', prev=null;
    for (const [,a,b,z] of rows) {
      const p=[a,b,z];
      const contiguous=prev && Math.hypot(...p.map((v,i)=>v-prev[i]))<0.55;
      d+=`${contiguous?'L':'M'}${point(p)} `;
      prev=p;
    }
    parts.push(`<path d="${d}" fill="none" stroke="${stroke[id]}" stroke-width="1.5" stroke-linejoin="round" stroke-linecap="round"><title>Chain ${id}</title></path>`);
  }
  for (const [res,...p] of structure.footprint) {
    const [a,b]=rotate(p);
    parts.push(`<circle cx="${x(a)}" cy="${y(b)}" r="3.1" fill="#e7ae42" stroke="#795315" stroke-width=".65"><title>TNF chain A residue ${res}</title></circle>`);
  }
  parts.push(`<circle cx="${x(anchor[0])}" cy="${y(anchor[1])}" r="3" fill="#774376"><title>Fv C-terminal reference point</title></circle>`);
  return parts.join('')+'</svg>';
}

function renderGeometry() {
  const values = Object.fromEntries(specs.map(([key,,,,,unit]) => {
    const value = el(`geom-${key}`).valueAsNumber;
    el(`geom-${key}-value`).value = `${value}${unit==='°'?'°':' nm'}`;
    return [key,value];
  }));
  const show = el('geometry-shell').checked;
  const result = pose(values);
  el('geometry-view').innerHTML = scene(result,values.radius,show);
  const entries = [['TNF to plane',result.tnfClearance],['Fv to plane',result.fvClearance]];
  if (show) entries.push(['Domain to plane',result.envelopeClearance]);
  const overlaps = v => v < -1e-8;
  const minimum = show ? result.minimumOriginHeight : Math.max(0,values.height-Math.min(result.tnfClearance,result.fvClearance));
  const nm = v => `${(Math.abs(v)<0.005 ? 0 : v).toFixed(2)} nm`;
  el('geometry-results').innerHTML = entries.map(([label,v]) =>
    `<div><dt>${label}</dt><dd class="${overlaps(v)?'geom-overlap':''}">${nm(v)}<small>${overlaps(v)?'Overlaps plane':'No overlap'}</small></dd></div>`).join('')
    +`<div><dt>Minimum origin height</dt><dd>${nm(minimum)}<small>At this tilt</small></dd></div>`;
  el('geometry-summary').textContent = overlaps(result.tnfClearance)
    ? 'The TNF trimer itself crosses the plane. This is not a feasible target pose in this model.'
    : overlaps(result.fvClearance) ? 'TNF clears the plane, but the deposited Fv does not at this assumed orientation.'
    : show && overlaps(result.envelopeClearance) ? 'The resolved TNF and Fv clear the plane, but the hypothetical domain crosses it.'
    : 'All displayed components clear the plane. That rules out this one type of overlap, not other obstacles to binding.';
}

function reset() {
  for (const [key] of specs) el(`geom-${key}`).value=defaults[key];
  el('geometry-shell').checked=true;
  renderGeometry();
}
el('geometry-sliders').addEventListener('input',renderGeometry);
el('geometry-shell').addEventListener('change',renderGeometry);
el('geometry-reset').addEventListener('click',reset);
reset();
