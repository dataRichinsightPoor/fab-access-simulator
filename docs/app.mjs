import {preset, run, sweep, readProfile} from './model.mjs';

const $=id=>document.getElementById(id);
const names=['Reference','Variant, same exposure','Exposure alone','Variant, own exposure'];
const colors=['#244d78','#9b5025','#5f676b','#007b65'];
const dashes=['','9 5','2 5',''];
const formats=['reference','variant'];
let activeConfig, activeResults, lastSweep;
const pct = v => `${(v*100).toFixed(1)}%`;
const fixed = v => v.toFixed(1);
const signed = v => `${v>=0?'+':'−'}${Math.abs(v).toFixed(1)}`;

$('formats').innerHTML=formats.map((s,i)=>`<fieldset>
  <legend>${i?'Redesigned variant':'Reference antibody'}</legend>
  <label>Exposure profile<select id="${s}-mode"><option value="decay">Exponential decay</option><option value="points">Custom time/concentration</option></select></label>
  <div class="fields" id="${s}-decay">
    <label>Initial tissue level (nM)<input id="${s}-initial" type="number" min="0" step="any" required></label>
    <label>Tissue half-life (h)<input id="${s}-half" type="number" min="0.000001" step="any" required></label>
  </div>
  <div id="${s}-points" hidden>
    <p class="profile-help">One <code>time_h,free_tissue_nM</code> pair per line. Start at 0 h and cover the full run. Points are joined by straight lines.</p>
    <label>Profile CSV<textarea id="${s}-csv" spellcheck="false"></textarea></label>
  </div>
  <div class="fields">
    <label>Association (nM⁻¹ h⁻¹)<input id="${s}-kon" type="number" min="0" step="any" required></label>
    <label>Dissociation (h⁻¹)<input id="${s}-koff" type="number" min="0" step="any" required></label>
    <label>Accessible sites (%)<input id="${s}-access" type="number" min="0" max="100" step="any" required></label>
    <label>Occupied-site loss (h⁻¹)<input id="${s}-loss" type="number" min="0" step="any" required></label>
  </div>
</fieldset>`).join('');
$('legend').innerHTML=names.map((n,i)=>`<span><i class="swatch ${['','dashed','dotted',''][i]}" style="--swatch:${colors[i]}" aria-hidden="true"></i>${n}</span>`).join('');

function mode(s) {
  const custom=$(`${s}-mode`).value==='points';
  $(`${s}-points`).hidden=!custom;
  $(`${s}-decay`).hidden=custom;
}

function fill(c){
  $('scenario-name').value=typeof c.label==='string'?c.label.slice(0,90):'Custom scenario';
  $('duration').value=c.duration_h; $('step').value=c.dt_h; $('threshold').value=c.occupancy_threshold*100;
  for(const s of formats){
    const b=c[`${s}_binding`], source=c[`${s}_exposure`];
    const e=c.constant_concentration_nM==null?source:{points:[[0,c.constant_concentration_nM],[c.duration_h,c.constant_concentration_nM]]};
    $(`${s}-mode`).value=e.points?'points':'decay';
    $(`${s}-initial`).value=e.initial_nM ?? 15; $(`${s}-half`).value=e.half_life_h ?? 10;
    $(`${s}-csv`).value=`time_h,free_tissue_nM\n${(e.points ?? [[0,15],[c.duration_h,1]]).map(p=>p.join(',')).join('\n')}`;
    $(`${s}-kon`).value=b.kon_per_nM_h; $(`${s}-koff`).value=b.koff_per_h;
    $(`${s}-access`).value=b.accessible_fraction*100; $(`${s}-loss`).value=b.bound_loss_per_h ?? 0;
    mode(s);
  }
}

function value(id){
  const v=$(id).valueAsNumber;
  if (!Number.isFinite(v)) throw Error(`Enter a valid number for ${$(id).parentElement.firstChild.textContent.trim()}.`);
  return v;
}
function getConfig(){
  const c={label:$('scenario-name').value.trim() || 'Custom scenario',duration_h:value('duration'),dt_h:value('step'),occupancy_threshold:value('threshold')/100};
  for(const s of formats){
    c[`${s}_binding`]={kon_per_nM_h:value(`${s}-kon`),koff_per_h:value(`${s}-koff`),accessible_fraction:value(`${s}-access`)/100,bound_loss_per_h:value(`${s}-loss`)};
    c[`${s}_exposure`]=$(`${s}-mode`).value==='points'
      ? {points:readProfile($(`${s}-csv`).value,c.duration_h)}
      : {initial_nM:value(`${s}-initial`),half_life_h:value(`${s}-half`)};
  }
  return c;
}

function chart(rows,c,occupancy){
  const w=540,h=159, l=45,t=11,pw=481,ph=113;
  const col=occupancy?2:1;
  const maximum=occupancy?1:rows.reduce((m,s)=>s.rows.reduce((max,p)=>Math.max(max,p[1]),m),0.001);
  const x=hour=>l+hour/c.duration_h*pw;
  const y=v=>t+ph*(1-v/maximum);
  const pieces=[`<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${occupancy?'Calculated fraction of target occupied over time':'Prescribed free tissue antibody concentration over time'}">`];
  for(let i=0;i<=4;i++){
    const v=maximum*i/4, yy=y(v), tt=c.duration_h*i/4;
    pieces.push(`<line x1="${l}" x2="${l+pw}" y1="${yy}" y2="${yy}" stroke="#d4dfd5"/><text x="${l-7}" y="${yy+4}" text-anchor="end">${occupancy?v*100:Number(v.toPrecision(3))}</text><text x="${x(tt)}" y="143" text-anchor="middle">${Number(tt.toPrecision(3))}</text>`);
  }
  pieces.push('<text x="284" y="157" text-anchor="middle">Hours</text>');
  if(occupancy) pieces.push(`<line x1="${l}" x2="${l+pw}" y1="${y(c.occupancy_threshold)}" y2="${y(c.occupancy_threshold)}" stroke="#787c68" stroke-dasharray="4 4"/>`);
  else pieces.push('<text x="524" y="24" text-anchor="end">Two paired exposure profiles</text>');
  for(const [i,s] of rows.entries()){
    if(!occupancy && ![0,3].includes(i)) continue;
    const stride=Math.max(1,Math.ceil(s.rows.length/1300));
    const shape=s.rows.filter((_,k)=>k%stride===0||k===s.rows.length-1);
    pieces.push(`<polyline vector-effect="non-scaling-stroke" fill="none" stroke="${colors[i]}" stroke-width="2.2" stroke-dasharray="${dashes[i]}" points="${shape.map(r=>`${x(r[0]).toFixed(2)},${y(r[col]).toFixed(2)}`).join(' ')}"><title>${names[i]}</title></polyline>`);
  }
  return pieces.join('')+'</svg>';
}

function render(c, rows){
  activeConfig=c; activeResults=rows; lastSweep=undefined;
  $('active-scenario').textContent=`${c.label} · ${c.duration_h} h run · Threshold ${pct(c.occupancy_threshold)}`;
  $('occupancy').innerHTML=chart(rows,c,true); $('exposure').innerHTML=chart(rows,c,false);
  $('threshold-head').textContent=`Time ≥ ${pct(c.occupancy_threshold)}`;
  $('metrics').innerHTML=rows.map((r,i)=>`<tr><th scope="row">${names[i]}</th><td>${pct(r.mean)}</td><td>${fixed(r.above)} h</td><td>${pct(r.end)}</td></tr>`).join('');
  const [ref,,,mod]=rows, dMean=(mod.mean-ref.mean)*100, dTime=mod.above-ref.above;
  $('interpretation').textContent=`Variant at its own exposure versus reference: ${signed(dMean)} percentage points in mean occupancy; ${signed(dTime)} h at or above the chosen threshold. A threshold is an assumption, not an efficacy cutoff.`;
  $('sweep-result').replaceChildren();
  $('save-sweep').hidden=true;
  $('run-sweep').disabled=false;
  $('dirty').textContent='';
  $('error').hidden=true;
}

function download(content, mime, filename){
  const url=URL.createObjectURL(new Blob([content],{type:mime}));
  const a=document.createElement('a');
  a.href=url; a.download=filename; document.body.append(a); a.click(); a.remove();
  setTimeout(()=>URL.revokeObjectURL(url),2000);
}
function error(e){
  $('error').textContent=e.message || String(e);
  $('error').hidden=false;
}
function selectPreset(name){
  document.querySelectorAll('[data-preset]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.preset===name)));
}
document.querySelectorAll('[data-preset]').forEach(b=>b.addEventListener('click',()=>{
  selectPreset(b.dataset.preset);
  const c=preset(b.dataset.preset);
  fill(c); render(c,run(c));
}));
$('parameters').addEventListener('input',()=>{
  $('error').hidden=true;
  $('dirty').textContent='Inputs changed. Run comparison to update the results and sensitivity sweep.';
  $('run-sweep').disabled=true;
  selectPreset();
});
formats.forEach(s=>$(`${s}-mode`).addEventListener('change',()=>mode(s)));
$('parameters').addEventListener('submit',e=>{
  e.preventDefault();
  try{const c=getConfig(); render(c,run(c));}catch(err){error(err);}
});
$('save-json').addEventListener('click',()=>{
  try{
    const c=getConfig(); run(c); $('error').hidden=true;
    download(JSON.stringify(c,null,2)+'\n','application/json','fab-access-scenario.json');
  }catch(err){error(err);}
});
$('load-json').addEventListener('change',async e=>{
  const file=e.target.files[0];
  if(!file) return;
  try{
    if(file.size>200000) throw Error('Keep a scenario file under 200 KB.');
    const c=JSON.parse(await file.text());
    if(typeof c.label!=='string') c.label='Imported scenario';
    c.label=c.label.slice(0,90);
    const rows=run(c); selectPreset(); fill(c); render(c,rows);
  }catch(err){error(err);}
  finally{e.target.value='';}
});
$('save-csv').addEventListener('click',()=>{
  const csv=['condition,time_h,free_tissue_nM,fraction_total_target_occupied',
    ...activeResults.flatMap(s=>s.rows.map(r=>[s.id,...r].join(',')))].join('\n')+'\n';
  download(csv,'text/csv','fab-access-timecourses.csv');
});
$('run-sweep').addEventListener('click',()=>{
  if(!activeConfig) return;
  lastSweep=sweep(activeConfig);
  const cols=[...new Set(lastSweep.map(r=>r.multiplier))];
  const caps=[...new Set(lastSweep.map(r=>r.access))];
  const heads=cols.map(x=>`<th scope="col">${x}×</th>`).join('');
  const body=caps.map(f=>`<tr><th scope="row">${pct(f)}</th>${lastSweep.filter(r=>r.access===f).map(r=>`<td class="${r.delta_hours>0?'positive':'negative'}">${signed(r.delta_hours)} h</td>`).join('')}</tr>`).join('');
  $('sweep-result').innerHTML=`<table><caption>Variant association multiplier at each accessible fraction; relative to reference time at or above threshold</caption><thead><tr><th scope="col">Accessible sites</th>${heads}</tr></thead><tbody>${body}</tbody></table>`;
  $('save-sweep').hidden=false;
});
$('save-sweep').addEventListener('click',()=>{
  if(!lastSweep) return;
  download([Object.keys(lastSweep[0]).join(','),...lastSweep.map(r=>Object.values(r).join(','))].join('\n')+'\n','text/csv','fab-access-sweep.csv');
});
const c=preset();
fill(c); selectPreset('tradeoff'); render(c,run(c));
