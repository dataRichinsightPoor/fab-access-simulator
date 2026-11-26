import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {spawnSync} from 'node:child_process';
import {exactStep, keys, preset, readProfile, run, summarize, sweep} from '../docs/model.mjs';

const close=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);

test('exact constant exposure with occupied-site loss',()=>{
  const b={kon_per_nM_h:0.2,koff_per_h:0.1,bound_loss_per_h:0.3,accessible_fraction:0.6};
  close(exactStep(0,8,b,2),0.6*1.6/2*(1-Math.exp(-4)));
  close(exactStep(0.6,0,b,3),0.6*Math.exp(-1.2));
  assert.equal(exactStep(0.4,0,{...b,koff_per_h:0,bound_loss_per_h:0},4),0.4);
});

test('baseline agrees with the independent Python implementation at every sample',()=>{
  const cfg=JSON.parse(readFileSync(new URL('../example.json',import.meta.url)));
  const py=spawnSync(process.env.PYTHON || 'python',['-c','import json; from fab_access import compare; print(json.dumps(compare(json.load(open("example.json")))))'],{encoding:'utf8',cwd:new URL('../',import.meta.url)});
  assert.equal(py.status,0,py.stderr);
  const expected=JSON.parse(py.stdout);
  run(cfg).forEach((s,i)=>{
    const expectedRows=expected[keys[i]];
    assert.equal(s.rows.length,expectedRows.length);
    s.rows.forEach((row,j)=>row.forEach((v,k)=>close(v,expectedRows[j][k])));
    close(s.mean,[0.4255465313904402,0.26418387884663125,0.6699417784266697,0.442162954528588][i]);
    close(s.above,[29.501994224463267,12.410294094704675,62.36108774774942,27.228376618095943][i]);
  });
});

test('piecewise profile interpolation and grid honors abrupt exposure changes',()=>{
  const c=preset('slow');
  c.dt_h=0.13;
  const rows=run(c)[0].rows;
  assert.ok(rows.some(([t,C])=>t===12.02&&C===0));
  assert.ok(rows.some(([t,C])=>t===12&&C===4));
  assert.ok(rows.filter(([t])=>t>12.02).every(([,C])=>C===0));
  assert.ok(rows.at(-1)[2] < rows.find(([t])=>t===12.02)[2]);
});

test('the same KD can yield slower loading and longer washout',()=>{
  const [ref,,,variant]=run(preset('slow'));
  assert.ok(ref.rows.find(r=>r[0]===3)[2]>variant.rows.find(r=>r[0]===3)[2]);
  assert.ok(ref.end<variant.end);
});

test('the accessible-fraction ceiling holds under saturating exposure',()=>{
  const c=preset('ceiling');
  c.variant_exposure.initial_nM=1e7;
  const obs=run(c).at(-1);
  assert.ok(obs.rows.every(([, , B])=>B>=0&&B<=0.4));
  assert.equal(obs.above,0);
});

test('turnover, no drug, and a threshold crossing behave as expected',()=>{
  const c=preset('turnover');
  const [r,,,v]=run(c);
  assert.ok(v.mean<r.mean);
  c.constant_concentration_nM=0;
  assert.ok(run(c).every(s=>s.peak===0));
  close(summarize([[0,1,0],[2,1,1],[4,1,0]],0.5).above,2);
});

test('strict input handling rejects malformed, negative, duplicate and short profiles',()=>{
  assert.deepEqual(readProfile('# header\r\ntime_h,free_tissue_nM\r\n0,3\r\n24,0',24),[[0,3],[24,0]]);
  for (const csv of ['0,3\n0,2','1,2\n24,1','0,2\n20,1','0,2\n24,-1','0,2\n24,','0,2\n24,NaN','0,2\n24,5,9']) {
    assert.throws(()=>readProfile(csv,24));
  }
  assert.throws(()=>run({...preset(),dt_h:1e-6}));
  assert.throws(()=>run({...preset(),occupancy_threshold:0}));
});

test('sweep preserves direction, thresholds, and expected reference comparison',()=>{
  const cells=sweep(preset());
  assert.equal(cells.length,42);
  assert.equal(cells.find(r=>r.access===0.5&&r.multiplier===2).hours,0);
  close(cells.find(r=>r.access===0.8&&r.multiplier===2).hours,45.13093837269157);
});
