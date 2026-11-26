// Usage: node scripts/build-5yoy.mjs 5YOY.pdb > docs/5yoy-data.mjs
import {createHash} from 'node:crypto';
import {readFileSync} from 'node:fs';

const pdb = readFileSync(process.argv[2], 'utf8');
if (/^HEADER.{56}4N0U/m.test(pdb)) {
  const traces = Object.fromEntries([...'ABDE'].map(c=>[c,[]]));
  for (const l of pdb.split('\n')) {
    if (!l.startsWith('ATOM  ') || !traces[l[21]] || l.slice(12,16).trim()!=='CA' || ![' ','A'].includes(l[16])) continue;
    traces[l[21]].push([l.slice(22,27).trim(),...[30,38,46].map(i=>Number(l.slice(i,i+8))/10)]);
  }
  const all=Object.values(traces).flat();
  const center=[1,2,3].map(i=>all.reduce((s,p)=>s+p[i],0)/all.length);
  for (const rows of Object.values(traces)) for(const p of rows) for(let i=1;i<4;i++) p[i]=Number((p[i]-center[i-1]).toFixed(4));
  console.log('// Generated from public PDB 4N0U, deposited asymmetric unit.');
  console.log(`export default ${JSON.stringify({id:'4N0U',sourceSHA256:createHash('sha256').update(pdb).digest('hex'),traces,footprint:[]})};`);
  process.exit(0);
}
if (!/^HEADER.{56}5YOY/m.test(pdb)) throw Error('Expected the RCSB 5YOY PDB file.');
const chains = Object.fromEntries([...`ABCDEFGHI`].map(c => [c, []]));
for (const line of pdb.split('\n')) {
  if (!line.startsWith('ATOM  ') || !chains[line[21]] || ![' ', 'A'].includes(line[16])) continue;
  const element = line.slice(76,78).trim();
  if (['H','D'].includes(element)) continue;
  chains[line[21]].push({
    name:line.slice(12,16).trim(), res:line.slice(22,27).trim(),
    p:[30,38,46].map(i => Number(line.slice(i,i+8))/10),
  });
}
const dist2 = (a,b) => a.reduce((s,x,i) => s+(x-b[i])**2,0);
const paired = (target,group) => [...group].map(c => [
  c, chains[c].reduce((n,b) => n+chains[target].filter(a => dist2(a.p,b.p)<=0.45**2).length,0),
]).sort((a,b) => b[1]-a[1])[0][0];
const pairings = [...'ABC'].map(c => ({tnf:c, light:paired(c,'DEF'), heavy:paired(c,'GHI')}));
const ca = c => chains[c].filter(a => a.name==='CA');
const mean = points => [0,1,2].map(i => points.reduce((s,p) => s+p[i],0)/points.length);
const sub = (a,b) => a.map((v,i) => v-b[i]);
const norm = a => Math.hypot(...a);
const unit = a => a.map(x => x/norm(a));
const dot = (a,b) => a.reduce((s,x,i) => s+x*b[i],0);
const cross = (a,b) => [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const origin = mean([...`ABC`].flatMap(c => ca(c).map(a => a.p)));
const up = unit(sub(origin,mean([...`ABC`].map(c => ca(c)[0].p))));
const depth = unit(cross(sub(mean(ca('A').map(a=>a.p)),origin),up));
const right = cross(up,depth);
const frame = p => [right,up,depth].map(axis => Number(dot(sub(p,origin),axis).toFixed(4)));
const footprint = new Set([24,65,67,71,104,105,106,107,108,111,138,140,141]);
const result = {
  id:'5YOY', pairings,
  sourceSHA256:createHash('sha256').update(pdb).digest('hex'),
  traces:Object.fromEntries([...'ABCDEFGHI'].map(c => [c,ca(c).map(a => [a.res,...frame(a.p)])])),
  footprint:ca('A').filter(a => footprint.has(Number(a.res))).map(a => [a.res,...frame(a.p)]),
};
console.log('// Generated from public PDB 5YOY. See scripts/build-5yoy.mjs and README.');
console.log(`export default ${JSON.stringify(result)};`);
