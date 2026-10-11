"""python optimize_vane.py --target 13 ; writes grid_results.csv."""
import argparse,csv
import numpy as np
from vane_model import Design,calculate
p=argparse.ArgumentParser();p.add_argument('--target',type=float,default=13);a=p.parse_args()
if a.target<=0:p.error('Target must be positive')
rows=[];best=None
for k in np.linspace(1.5,3,31):
 for alpha in np.linspace(4,15,23):
  d=Design(float(k),float(alpha),1);z=calculate(d,720)
  width=max(8,a.target/z['displacement'])
  valid=z['valid'] and width<=18
  row=dict(K=float(k),alpha_deg=float(alpha),required_width_mm=width,feasible=valid,
   rmin_mm=z['rmin'],max_pressure_angle_deg=z['max_pressure_angle'])
  rows.append(row)
  if valid and (best is None or width<best['required_width_mm']):best=row
with open('grid_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
if best:
 # Re-evaluate chosen candidate at finer geometry/integration resolution.
 fine=calculate(Design(best['K'],best['alpha_deg'],1),5760,1600)
 refined_width=max(8,a.target/fine['displacement'])
 print('Grid candidate:',best)
 print(f'Refined required width = {refined_width:.6f} mm for {a.target} cc/rev')
else:print('No feasible candidate in this grid')
print('Saved grid_results.csv; grid optimum is not a global-optimum guarantee.')
