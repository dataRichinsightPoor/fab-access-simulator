import {test} from 'node:test';
import assert from 'node:assert/strict';
import {structure, defaults, pose, atomRadius} from '../docs/geometry.mjs';

const near = (actual,expected) => assert.ok(Math.abs(actual-expected)<1e-9,`${actual} != ${expected}`);
test('coordinate subset covers the deposited TNF trimer, its selected Fv, and documented epitope',() => {
  assert.deepEqual(Object.keys(structure.traces),['A','B','C','F','I']);
  assert.equal(structure.fvAtoms.length,1801);
  assert.ok(structure.tnfAtoms.length>3000);
  assert.deepEqual(structure.footprint.map(([res])=>+res),[24,65,67,71,104,105,106,107,108,111,138,140,141]);
  const origins = ['A','B','C'].map(c => structure.traces[c][0].slice(1));
  for (const i of [0,1,2]) assert.ok(Math.abs(origins.reduce((s,p)=>s+p[i],0)/3)<0.0001);
});

test('rigid tilt and lift preserve the TNF–Fv interface distances', () => {
  const [, ...a]=structure.footprint[0];
  const b=structure.fvAtoms[0];
  const distance = (p,q) => Math.hypot(...p.map((v,i)=>v-q[i]));
  const {rotate}=pose({...defaults,tilt:70,height:12});
  near(distance(a,b),distance(rotate(a),rotate(b)));
});

test('known pose, atom padding, and domain tangency are analytic',() => {
  const data={end:[2,1,0],fvAtoms:[[0,0,0],[2,0,0]],tnfAtoms:[[0,1,0]]};
  const value=pose({height:1,tilt:0,offset:1,radius:3,bearing:90},data);
  near(value.tnfClearance,2-atomRadius);
  near(value.fvClearance,1-atomRadius);
  near(value.envelopeClearance,0);
  near(value.minimumOriginHeight,1);
  const tilted=pose({height:1,tilt:-90,offset:0,radius:1,bearing:0},data);
  near(tilted.fvClearance,-1-atomRadius);
  near(tilted.envelopeClearance,-2);
  near(tilted.minimumOriginHeight,3);
});

test('worked example separates Fv clearance from a hypothetical domain overlap', () => {
  const result=pose(defaults);
  assert.ok(result.tnfClearance>0);
  assert.ok(result.fvClearance>0);
  assert.ok(result.envelopeClearance<0);
  const cleared=pose({...defaults,height:result.minimumOriginHeight});
  near(Math.min(cleared.tnfClearance,cleared.fvClearance,cleared.envelopeClearance),0);
});
