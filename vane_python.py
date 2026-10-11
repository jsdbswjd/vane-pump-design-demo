"""python vane_python.py (sliders); python vane_python.py --save (static plot)."""
import sys
if '--save' in sys.argv:
 import matplotlib;matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from vane_model import Design,calculate
fig,axs=plt.subplots(1,3,figsize=(15,7));fig.subplots_adjust(bottom=.36,top=.79,wspace=.3)
fig.suptitle('VANE PUMP | patent-informed geometry exploration',fontsize=18,y=.97)
fig.text(.5,.92,'KR101646052B1 · C1 tangent reconstruction · Assumed dimensions / ideal displacement',ha='center',color='#956114')
status=fig.text(.06,.845,'',fontsize=11)
sliders=[]
for i,(name,lo,hi,v,step) in enumerate([('Profile K',1.5,3,2,.01),('Tangent angle [deg]',4,15,10,.1),('Width [mm]',8,18,14,.5),('Speed [rpm]',500,5000,3000,100)]):
 sliders.append(Slider(fig.add_axes([.23,.25-i*.052,.65,.025]),name,lo,hi,valinit=v,valstep=step))
def update(_=None):
 d=Design(*[s.val for s in sliders]);z=calculate(d)
 for ax in axs:ax.clear();ax.grid(alpha=.15)
 x=z['r']*np.cos(z['phi']);y=z['r']*np.sin(z['phi'])
 axs[0].plot(x,y,color='#177c88',lw=2,label='Cam inner profile')
 axs[0].plot(26.75*np.cos(z['phi']),26.75*np.sin(z['phi']),color='#9ba5b0',label='Outer boundary')
 axs[0].plot(19*np.cos(z['phi']),19*np.sin(z['phi']),color='#253a58',label='Assumed rotor')
 for a in np.arange(12)*2*np.pi/12:
  rr=np.interp(a,z['phi'],z['r']);axs[0].plot([19*np.cos(a),rr*np.cos(a)],[19*np.sin(a),rr*np.sin(a)],color='#d68c2c')
 axs[0].set(aspect='equal',xlabel='x [mm]',ylabel='y [mm]',title='1. Reconstructed cross-section');axs[0].legend(fontsize=8)
 axs[1].plot(np.rad2deg(z['phi'][:-1]),z['volume'],color='#177c88');axs[1].set(xlabel='Rotor angle [deg]',ylabel='One chamber volume [cc]',title='2. Ideal chamber expansion')
 ks=np.linspace(1.5,3,61);ys=[]
 for k in ks:
  v=calculate(Design(k,d.alpha_deg,d.width_mm,d.rpm));ys.append(v['displacement'] if v['valid'] else np.nan)
 axs[2].plot(ks,ys,color='#177c88');axs[2].scatter([d.k],[z['displacement']],color='#253a58');axs[2].set(xlabel='Profile K [-]',ylabel='Ideal displacement [cc/rev]',title='3. Feasible K sweep')
 status.set_text(('FEASIBLE' if z['valid'] else 'INVALID: '+', '.join(z['reasons']))+f" | displacement {z['displacement']:.3f} cc/rev | ideal flow {z['ideal_flow']:.2f} L/min")
 fig.canvas.draw_idle()
for s in sliders:s.on_changed(update)
fig.text(.06,.025,'No leakage / finite vane thickness / port timing / wear prediction. Values are not the patent published values.',fontsize=10)
update()
if '--save' in sys.argv:fig.savefig('vane_preview.png',dpi=140)
else:plt.show()
