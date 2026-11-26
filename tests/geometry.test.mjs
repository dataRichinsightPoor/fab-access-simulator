import {test} from 'node:test';
import assert from 'node:assert/strict';
import {structure,fcrn,camera} from '../docs/geometry.mjs';

const near = (a,b) => assert.ok(Math.abs(a-b)<1e-9,`${a} != ${b}`);
test('5YOY contains three separate bound Fvs in biological assembly 1',() => {
  assert.equal(structure.id,'5YOY');
  assert.deepEqual(Object.keys(structure.traces),[...'ABCDEFGHI']);
  assert.deepEqual(structure.pairings.map(p=>p.tnf+p.light+p.heavy),['AFI','BDG','CEH']);
  assert.deepEqual(structure.footprint.map(([res])=>+res),[24,65,67,71,104,105,106,107,108,111,138,140,141]);
  assert.deepEqual(Object.values(structure.traces).map(r=>r.length),[154,154,154,108,108,108,126,126,125]);
});
test('4N0U is the separate deposited asymmetric unit with the four expected protein chains',() => {
  assert.equal(fcrn.id,'4N0U');
  assert.deepEqual(Object.keys(fcrn.traces),[...'ABDE']);
  assert.ok(fcrn.traces.D.length>500);
  assert.ok(fcrn.traces.E.length>180);
  assert.equal(fcrn.footprint.length,0);
});
test('camera retains interface distances and closes full turn without changing the structure',() => {
  const [, ...a]=structure.footprint[0], [, ...b]=structure.traces.F[0];
  const distance=(p,q)=>Math.hypot(...p.map((v,i)=>v-q[i]));
  const r=camera(61,-25);
  near(distance(a,b),distance(r(a),r(b)));
  for (const [i,x] of camera(360,0)(a).entries()) near(x,a[i]);
  near(camera(90,0)([1,0,0])[2],1);
});
