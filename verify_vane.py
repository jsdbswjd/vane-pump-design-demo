from dataclasses import replace
from vane_model import *
import json
from pathlib import Path
base=Design();z=calculate(base)
# Independent geometry conditions at both tangent junctions.
aa,rs,g=quarter(base);t=g['tangent'];f=g['foot'];v=g['direction']
assert abs(np.dot(f,v))<1e-10
assert np.linalg.norm(t-f-np.dot(t,v)*v)<1e-10
assert abs(rs[0]-g['rmin'])<1e-10 and abs(rs[-1]-g['rmax'])<1e-10
# Symmetry and endpoints.
p=np.linspace(0,2*np.pi,91)
assert np.max(abs(radius(p,base)-radius(p+np.pi,base)))<1e-10
# Width scales displacement, RPM scales flow but not displacement.
assert abs(calculate(replace(base,width_mm=7))['displacement']/z['displacement']-.5)<1e-10
assert abs(calculate(replace(base,rpm=1500))['ideal_flow']/z['ideal_flow']-.5)<1e-10
assert calculate(replace(base,k=3))['valid'] is False
# Total volume of all adjacent chambers equals integrated annulus volume.
N=1440;step=N//12
integrand=.5*(z['r']**2-19**2)*14
annulus=np.sum((integrand[:-1]+integrand[1:])/2*np.diff(z['phi']))/1000
assert abs(sum(z['volume'][::step])-annulus)<1e-9
# Convergence in angular integration; outline sampling independently refined.
fine=calculate(base,5760,1600)
relative=abs(z['displacement']/fine['displacement']-1)
assert relative<.001
cases=[]
for k in [1.5,2,2.5,3]:
 for a in [4,10,15]:
  out=calculate(Design(k,a),720)
  cases.append(dict(k=k,a=a,disp=out['displacement'],rmin=out['rmin'],pa=out['max_pressure_angle'],valid=out['valid']))
Path('verification_cases.json').write_text(json.dumps(cases))
print('PASS: tangent conditions, symmetry, volume conservation, width/RPM scaling, infeasibility')
print('Displacement refinement relative change:',relative)
print('Baseline cc/rev:',z['displacement'])
# Exhaustive grid selection under explicitly assumed geometric constraints.
best=None
for k in np.linspace(1.5,3,31):
 for a in np.linspace(4,15,23):
  s=calculate(Design(float(k),float(a),1),720)
  if not s['valid']:continue
  b=max(8,13/s['displacement'])
  if b<=18 and (best is None or b<best['width_mm']):best=dict(k=float(k),alpha=float(a),width_mm=b)
print('Grid minimum width for 13 cc/rev:',best)
