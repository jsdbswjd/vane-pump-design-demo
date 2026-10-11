"""Patent-informed C1 tangent reconstruction, not exact proprietary geometry."""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Design:
    k: float=2.0
    alpha_deg: float=10.0
    width_mm: float=14.0
    rpm: float=3000.0
    outer_diameter_mm: float=53.5
    wall_mm: float=2.0
    rotor_radius_mm: float=19.0
    vanes: int=12
    max_extension_mm: float=6.0

def quarter(d, samples=400):
    if not (1.5<=d.k<=3 and 4<=d.alpha_deg<=15 and d.width_mm>0):
        raise ValueError('Outside parameter limits')
    a=np.deg2rad(d.alpha_deg)
    rmax=d.outer_diameter_mm/2-d.wall_mm
    R=rmax/(d.k+2)
    theta0=2*a
    t=np.array([R*(theta0-np.sin(theta0))-np.pi*R,R*(1-np.cos(theta0)+d.k)])
    direction=np.array([np.sin(a),np.cos(a)])
    length=np.dot(t,direction)
    foot=t-length*direction
    rmin=np.linalg.norm(foot)
    beta=np.arctan2(foot[1],-foot[0])
    th=np.linspace(theta0,np.pi,samples)
    cycl=np.column_stack([R*(th-np.sin(th))-np.pi*R,R*(1-np.cos(th)+d.k)])
    arc_angles=np.linspace(0,beta,100)
    arc=np.column_stack([-rmin*np.cos(arc_angles),rmin*np.sin(arc_angles)])
    line=np.linspace(foot,t,100)
    points=np.vstack([arc,line[1:],cycl[1:]])
    ang=np.arctan2(points[:,1],-points[:,0])
    radii=np.linalg.norm(points,axis=1)
    if length<=0 or beta<0 or beta>np.pi/2 or np.any(np.diff(ang)<-1e-10) or rmin>=rmax:
        raise ValueError('Tangent geometry invalid')
    return ang,radii,{'R':R,'rmin':rmin,'rmax':rmax,'tangent':t,'foot':foot,'direction':direction}

def radius(phi,d,samples=400):
    angles,rs,_=quarter(d,samples)
    # quarter angle beta measured from either x-axis towards y-axis
    beta=np.arctan2(np.abs(np.sin(phi)),np.abs(np.cos(phi)))
    return np.interp(beta,angles,rs)

def calculate(d,steps=1440,samples=400):
    if steps%d.vanes: raise ValueError('Angular step count must be a multiple of vane count')
    phi=np.linspace(0,2*np.pi,steps+1)
    r=radius(phi,d,samples)
    _,_,geo=quarter(d,samples)
    extension=geo['rmax']-d.rotor_radius_mm
    reasons=[]
    if geo['rmin']-d.rotor_radius_mm<.5: reasons.append('Rotor clearance below 0.5 mm')
    if extension>d.max_extension_mm: reasons.append('Vane extension above limit')
    if d.wall_mm<1: reasons.append('Outer wall below 1 mm')
    # Ideal zero-thickness vanes; finite chamber sector integrated in mm^3.
    area=.5*(r*r-d.rotor_radius_mm**2)
    seg=(area[:-1]+area[1:])/2*(2*np.pi/steps)*d.width_mm
    cumulative=np.r_[0,np.cumsum(np.tile(seg,2))]
    shift=steps//d.vanes
    volumes=(cumulative[np.arange(steps)+shift]-cumulative[np.arange(steps)])/1000
    displacement=2*d.vanes*(volumes.max()-volumes.min())
    ext=r[:-1]-d.rotor_radius_mm
    vel=np.gradient(np.r_[ext[-1],ext,ext[0]],2*np.pi/steps)[1:-1]*(d.rpm*2*np.pi/60)
    # Pressure angle = angle between radial direction and contour normal.
    deriv=np.gradient(r,2*np.pi/steps)
    pressure_angle=np.rad2deg(np.arctan(np.abs(deriv/r)))
    if pressure_angle.max()>30: reasons.append('Pressure angle above 30 deg')
    return dict(phi=phi,r=r,volume=volumes,displacement=displacement,
        ideal_flow=displacement*d.rpm/1000,extension=extension,velocity=vel,
        max_pressure_angle=float(pressure_angle.max()),valid=not reasons,reasons=reasons,**geo)
