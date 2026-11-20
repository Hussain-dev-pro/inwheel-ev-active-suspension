import sys; sys.path.insert(0,'.')
import numpy as np, json
from core import *; from sim import *
ev=params(True); ice=params(False); D=json.load(open('out/lqr_designs.json'))
t,zr,w=road_random(60,seed=0); sk=5000
tb,zb,wb=road_bump(30,0.075,3.7,2.5)
# semi-active tuning on a 20 s road: (a) pure skyhook best comfort; (b) hybrid best comfort with tyre <= passive EV
_,_,w20=road_random(20,seed=3)
evS=summarize(simulate2(ev,('passive',22000,1500),w20),sk)
sky=[];hyb=[]
for csk in (3000,4500,6000):
    for cmin in (300,600):
        s=summarize(simulate2(ev,('hyb',1.0,csk,0,cmin,4500),w20),sk); sky.append((s['aw'],csk,cmin))
for al in (0.2,0.35,0.5,0.65):
  for csk in (3000,5000):
    for cg in (1500,3000):
      for cmin in (300,600,900):
        s=summarize(simulate2(ev,('hyb',al,csk,cg,cmin,4500),w20),sk)
        if s['tyre']<=evS['tyre']: hyb.append((s['aw'],al,csk,cg,cmin))
sk1=min(sky); hb=min(hyb)
Kc,_=lqr_lag(ev,1.0,1/0.03**2,10**D['lag_comfort']['lqt'],10**D['lag_comfort']['lr'])
Kg,_=lqr_lag(ev,1.0,1/0.03**2,10**D['lag_grip']['lqt'],10**D['lag_grip']['lr'])
L,Cm,Dm=kalman_design(ev)
pc=D['passive_comfort']
cases=[("Normal car, passive stock",ice,('passive',22000,1500)),
       ("EV, passive stock",ev,('passive',22000,1500)),
       ("EV, passive re-tuned",ev,('passive',pc['ks'],pc['cs'])),
       ("EV, semi-active skyhook",ev,('hyb',1.0,sk1[1],0,sk1[2],4500)),
       ("EV, semi-active hybrid",ev,('hyb',*hb[1:],4500)),
       ("EV, active LQR comfort-first",ev,('llqr',Kc)),
       ("EV, active LQG comfort-first",ev,('llqg',Kc,L,Cm,Dm)),
       ("EV, active LQR grip-first",ev,('llqr',Kg))]
R=dict(design=dict(sky=list(map(float,sk1[1:])),hyb=list(map(float,hb[1:])),
       eig_lqr_c=float(max(np.real(cl_eigs_lag(ev,Kc)))),eig_lqg_c=float(max(np.real(cl_eigs_lag(ev,Kc,L,Cm,Dm)))),
       eig_lqr_g=float(max(np.real(cl_eigs_lag(ev,Kg))))),rows=[])
B={}
for name,p,c in cases:
    Y=simulate3(p,c,w); s=summarize(Y,sk); Yb=simulate3(p,c,wb)
    row=dict(name=name,aw=s['aw'],label=iso_label(s['aw']),tyre=s['tyre'],trav3=3000*s['trav'],umax=s['umax'],urms=s['urms'],
             sat=100*s['satfrac'] if c[0] in('llqr','llqg') else 0.0,Pcons=float(np.mean(np.maximum(Y[sk:,4],0))) if c[0] in('llqr','llqg') else 0.0,
             bump_acc=float(np.abs(Yb[:,0]).max()),bump_tyre=float(Yb[:,2].max()),bump_trav=float(np.abs(Yb[:,1]).max()*1000))
    R['rows'].append(row); B[name]=Yb
    print("%-32s aw %.3f (%s) tyre %.4f 3s %.1fmm umax %.0f urms %.0f sat %.2f%% P %.0fW | hump acc %.2f tyre %.2f trav %.1fmm"%(
        name,row['aw'],row['label'],row['tyre'],row['trav3'],row['umax'],row['urms'],row['sat'],row['Pcons'],row['bump_acc'],row['bump_tyre'],row['bump_trav']))
print(R['design'])
json.dump(R,open('out/phase4.json','w'),indent=1)
np.savez('out/phase4_bump.npz',t=tb,zr=zb,names=np.array([c[0] for c in cases]),**{f"c{i}":B[c[0]] for i,c in enumerate(cases)})
