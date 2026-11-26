import structure from './5yoy-data.mjs';

export {structure};
export const atomRadius = 0.17; // Uniform heavy-atom padding, not a fitted surface.
export const defaults = {height:3, tilt:-45, radius:2, offset:3, bearing:-45};
const radians = degrees => degrees*Math.PI/180;

export function pose({height,tilt,radius,offset,bearing}, data=structure) {
  const theta = radians(tilt);
  const cos = Math.cos(theta), sin = Math.sin(theta);
  const rotate = ([x,y,z]) => [x*cos-y*sin, x*sin+y*cos+height, z];
  const beta = radians(bearing);
  const center = rotate([data.end[0]+offset*Math.cos(beta),data.end[1]+offset*Math.sin(beta),data.end[2]]);
  const clearance = atoms => atoms.reduce((min,[x,y]) => Math.min(min,x*sin+y*cos+height-atomRadius),Infinity);
  const fvClearance = clearance(data.fvAtoms);
  const tnfClearance = clearance(data.tnfAtoms);
  const envelopeClearance = center[1]-radius;
  return {
    rotate, center, tnfClearance, fvClearance, envelopeClearance,
    minimumOriginHeight:Math.max(0,height-Math.min(tnfClearance,fvClearance,envelopeClearance)),
  };
}
