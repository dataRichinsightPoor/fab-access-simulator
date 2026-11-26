export const version = '1.0.0';
export const keys = ['reference', 'variant_same_exposure', 'exposure_only', 'variant_own_exposure'];
const baseline = {
  label: 'Exposure versus access', duration_h: 72, dt_h: 0.05, occupancy_threshold: 0.5,
  reference_exposure: {initial_nM: 15, half_life_h: 10},
  variant_exposure: {initial_nM: 15, half_life_h: 22},
  reference_binding: {kon_per_nM_h: 0.08, koff_per_h: 0.18, accessible_fraction: 1, bound_loss_per_h: 0},
  variant_binding: {kon_per_nM_h: 0.04, koff_per_h: 0.18, accessible_fraction: 0.85, bound_loss_per_h: 0}
};

export function preset(name='tradeoff') {
  const c = structuredClone(baseline);
  if (name === 'slow') {
    c.label = 'Same equilibrium affinity, different rates';
    c.duration_h = 48;
    c.reference_exposure = c.variant_exposure = {points: [[0,4],[12,4],[12.02,0],[48,0]]};
    c.variant_binding = {...c.reference_binding, kon_per_nM_h: 0.02, koff_per_h: 0.045};
  }
  if (name === 'ceiling') {
    c.label = 'A hard access ceiling';
    c.variant_binding.accessible_fraction = 0.4;
    c.variant_exposure.initial_nM = 60;
  }
  if (name === 'turnover') {
    c.label = 'Loss from the occupied surface pool';
    c.variant_exposure = {...c.reference_exposure};
    c.variant_binding = {...c.reference_binding, bound_loss_per_h: 0.35};
  }
  return c;
}

function number(x, label, low, high) {
  if (typeof x !== 'number' || !Number.isFinite(x) || x < low || x > high) {
    throw Error(`${label} must be a number from ${low} to ${high}.`);
  }
  return x;
}

export function readProfile(text, duration) {
  if (text.length > 100000) throw Error('Keep a profile under 100 KB.');
  const lines = text.trim().split(/\r?\n/).filter(x=> x.trim() && !x.trim().startsWith('#'));
  if (/^\s*time_h\s*,\s*free_tissue_nM\s*$/i.test(lines[0] || '')) lines.shift();
  return validatePoints(lines.map((l,i) => {
    const cells = l.split(',').map(s=>s.trim());
    if (cells.length !== 2 || cells.some(v=>v === '')) throw Error(`Profile row ${i+1}: use time, concentration.`);
    return cells.map(Number);
  }), duration);
}

function validatePoints(points, duration) {
  if (!Array.isArray(points) || points.length < 2 || points.length > 2000) {
    throw Error('Provide 2 to 2,000 time/concentration pairs.');
  }
  let prior = -1;
  for (const [i,row] of points.entries()) {
    if (!Array.isArray(row) || row.length !== 2) throw Error(`Profile row ${i+1} needs two values.`);
    const t = number(row[0], `Time at row ${i+1}`, 0, 1e7);
    number(row[1], `Concentration at row ${i+1}`, 0, 1e7);
    if (t <= prior) throw Error('Profile times must strictly increase.');
    prior = t;
  }
  if (points[0][0] !== 0 || prior < duration) throw Error('A profile must start at time 0 and cover the full simulated duration.');
  return points;
}

function exposureFunction(e, duration, name) {
  if (!e || typeof e !== 'object') throw Error(`Missing ${name} exposure.`);
  if (e.points !== undefined) {
    const points = validatePoints(e.points, duration);
    return {knots: points.map(p=>p[0]), at(t) {
      let lo = 0, hi = points.length-1;
      while (hi-lo > 1) {
        const m = (lo+hi) >> 1;
        if (points[m][0] <= t) lo=m;
        else hi=m;
      }
      const [a,y] = points[lo], [b,z] = points[hi];
      return y+(z-y)*(Math.min(b,t)-a)/(b-a);
    }};
  }
  const initial = number(e.initial_nM, `${name} initial tissue concentration`, 0, 1e7);
  const half = number(e.half_life_h, `${name} half-life`, 0.000001, 1e7);
  return {knots: [], at: t => initial * 2 ** (-t / half)};
}

function binding(b,name) {
  if (!b || typeof b !== 'object') throw Error(`Missing ${name} binding parameters.`);
  return {
    kon_per_nM_h: number(b.kon_per_nM_h, `${name} association rate`, 0, 1e6),
    koff_per_h: number(b.koff_per_h, `${name} dissociation rate`, 0, 1e6),
    accessible_fraction: number(b.accessible_fraction, `${name} accessible fraction`, 0, 1),
    bound_loss_per_h: number(b.bound_loss_per_h ?? 0, `${name} occupied-site loss`, 0, 1e6)
  };
}

function config(c) {
  const duration = number(c.duration_h, 'Duration', 0.001, 1e5);
  const h = number(c.dt_h, 'Maximum step', 0.000001, 1e5);
  const threshold = number(c.occupancy_threshold, 'Occupancy threshold', 0.000001, 1);
  if (duration/h > 20000) throw Error('Use a larger step or shorter duration (maximum 20,000 intervals before adding profile knots).');
  const exposures = [
    exposureFunction(c.reference_exposure,duration,'Reference'),
    exposureFunction(c.variant_exposure,duration,'Variant')
  ];
  const binders = [binding(c.reference_binding,'Reference'),binding(c.variant_binding,'Variant')];
  const clamp = c.constant_concentration_nM;
  if (clamp !== undefined && clamp !== null) {
    number(clamp, 'Constant concentration', 0, 1e7);
    exposures.fill({knots: [], at: () => clamp});
  }
  const n = Math.ceil(duration/h);
  const regular = Array.from({length:n+1}, (_,i) => i === n ? duration : i*(duration/n));
  const times = [...new Set([...regular, ...exposures.flatMap(e=>e.knots).filter(t=>t<duration)])].sort((a,b)=>a-b);
  return {threshold, exposures, binders, times};
}

export function exactStep(B, C, b, h) {
  const a = b.kon_per_nM_h * C;
  const r = a + b.koff_per_h + (b.bound_loss_per_h ?? 0);
  if (!r) return B;
  const steady = b.accessible_fraction * (a/r);
  return B + (steady-B)*(-Math.expm1(-r*h));
}

function timecourse(times,e,b) {
  let B=0;
  return times.map((t,i)=>{
    if (i > 0) {
      const last = times[i-1];
      B = exactStep(B,e.at((last+t)/2),b,t-last);
    }
    return [t,e.at(t),B];
  });
}

export function summarize(rows,threshold) {
  let exposure=0, occupancy=0, above=0, peak=0;
  for (let i=1; i<rows.length; i++) {
    const [t,C,B]=rows[i], [s,D,A]=rows[i-1], h=t-s;
    exposure+=(C+D)*h/2;
    occupancy+=(B+A)*h/2;
    peak=Math.max(peak,B);
    if (Math.min(A,B)>=threshold) above+=h;
    else if (Math.max(A,B)>threshold) above+=h*(Math.max(A,B)-threshold)/Math.abs(B-A);
  }
  return {exposure, mean:occupancy/rows.at(-1)[0], above, peak, end:rows.at(-1)[2]};
}

export function run(c) {
  const {threshold, exposures:[ref,mod], binders:[normal,variant], times} = config(c);
  return [[ref,normal],[ref,variant],[mod,normal],[mod,variant]].map(([e,b],i)=>{
    const rows=timecourse(times,e,b);
    return {id:keys[i],rows,...summarize(rows,threshold)};
  });
}

export function sweep(c) {
  const {threshold, exposures:[ref,mod], binders:[normal,variant], times} = config(c);
  const baselineMetrics = summarize(timecourse(times,ref,normal),threshold);
  return [0.25,0.5,0.65,0.8,0.9,1].flatMap(f => [0.1,0.25,0.5,0.75,1,1.5,2].map(mult=>{
    const b = {...variant,accessible_fraction:f,kon_per_nM_h:variant.kon_per_nM_h*mult};
    const m = summarize(timecourse(times,mod,b),threshold);
    return {access:f,multiplier:mult,mean:m.mean,hours:m.above,delta_mean:m.mean-baselineMetrics.mean,delta_hours:m.above-baselineMetrics.above};
  }));
}
