import structure from './5yoy-data.mjs';
import fcrn from './4n0u-data.mjs';
export {structure,fcrn};
const rad = degrees => degrees*Math.PI/180;
export function camera(yaw,pitch) {
  const cy=Math.cos(rad(yaw)), sy=Math.sin(rad(yaw));
  const cp=Math.cos(rad(pitch)), sp=Math.sin(rad(pitch));
  return ([x,y,z]) => {
    const d=x*sy+z*cy;
    return [x*cy-z*sy,y*cp-d*sp,y*sp+d*cp];
  };
}
