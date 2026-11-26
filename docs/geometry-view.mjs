import {structure,fcrn,camera} from './geometry.mjs';

const el = id => document.getElementById(id);
const fmt = n => n.toFixed(2);
const cx=260,cy=156;
const palette={A:'#284f73',B:'#567c9b',C:'#7292a7',D:'#2f7f69',E:'#2f7f69',F:'#08614d',G:'#4e947e',H:'#4e947e',I:'#198266'};
const fcrnColors={A:'#326a91',B:'#388274',D:'#8770a4',E:'#a36222'};

function draw() {
  const yaw=el('geom-yaw').valueAsNumber, pitch=el('geom-pitch').valueAsNumber;
  el('geom-yaw-value').value=`${yaw}°`; el('geom-pitch-value').value=`${pitch}°`;
  const isTnf=el('geom-entry').value==='5YOY';
  const data=isTnf?structure:fcrn;
  const span=Math.max(...Object.values(data.traces).flat().map(([,x,y,z])=>Math.hypot(x,y,z)));
  const scale=141/(span+0.35);
  const xp=x=>fmt(cx+x*scale),yp=y=>fmt(cy-y*scale);
  for(const id of ['geom-tnf-options','geom-tnf-caption','geom-tnf-key']) el(id).hidden=!isTnf;
  for(const id of ['geom-fcrn-caption','geom-fcrn-key']) el(id).hidden=isTnf;
  const rotate=camera(yaw,pitch);
  const selected = !isTnf?'ABDE':el('geom-all').checked ? 'ABCDEFGHI' : `ABC${structure.pairings[0].light}${structure.pairings[0].heavy}`;
  const ink=[];
  for (const id of selected) {
    let prior=null;
    for (const [,x,y,z] of data.traces[id]) {
      const point=rotate([x,y,z]);
      if (prior && Math.hypot(...point.map((v,i)=>v-prior[i]))<0.52) {
        ink.push([(point[2]+prior[2])/2,`<path d="M${xp(prior[0])} ${yp(prior[1])}L${xp(point[0])} ${yp(point[1])}" fill="none" stroke="${(isTnf?palette:fcrnColors)[id]}" stroke-width="${isTnf&&'ABC'.includes(id)?'2.4':'1.8'}" stroke-linecap="round"/>`]);
      }
      prior=point;
    }
  }
  if (isTnf && el('geom-footprint').checked) for (const [res,...p] of structure.footprint) {
    const [x,y,z]=rotate(p);
    ink.push([z,`<circle cx="${xp(x)}" cy="${yp(y)}" r="3.25" fill="#efb64f" stroke="#684511" stroke-width=".8"><title>TNF chain A residue ${res}</title></circle>`]);
  }
  const label=isTnf?`deposited TNF trimer and ${selected.length===9?'three golimumab Fvs':'one golimumab Fv'}`:'deposited human FcRn, beta-2 microglobulin, albumin and one chain of Fc-YTE';
  el('geometry-view').innerHTML=`<svg viewBox="0 0 520 312" role="img" aria-label="Rotated orthographic C-alpha traces of ${label}">${ink.sort((a,b)=>a[0]-b[0]).map(x=>x[1]).join('')}<path d="M${fmt(492-2*scale)} 290H492" stroke="#203c31" stroke-width="2"/><text x="492" y="280" text-anchor="end" class="geom-svg-label">2 nm</text></svg>`;
}
function reset() {
  el('geom-yaw').value=24;
  el('geom-pitch').value=28;
  el('geom-all').checked=true;
  el('geom-footprint').checked=true;
  draw();
}
el('geometry-controls').addEventListener('input',draw);
el('geometry-reset').addEventListener('click',reset);
reset();
