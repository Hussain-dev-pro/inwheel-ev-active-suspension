import sys; sys.path.insert(0,'.')
import numpy as np, json
from core import *; from sim import *
ev=params(True); ice=params(False); D=json.load(open('out/lqr_designs.json')); P4=json.load(open('out/phase4.json'))
Kc,_=lqr_lag(ev,1.0,1/0.03**2,10**D['lag_comfort']['lqt'],10**D['lag_comfort']['lr'])
L,Cm,Dm=kalman_design(ev)          # designed on NOMINAL car, kept fixed
hb=P4['design']['hyb']
_,_,w=road_random(20,seed=7); sk=3000
N=100; rng=np.random.default_rng(42)
S=dict(ms=350*rng.uniform(0.9,1.3,N), mu=70*rng.uniform(0.85,1.15,N), kt=200000*rng.uniform(0.8,1.2,N), cs=1500*rng.uniform(0.8,1.0,N))
ctrls=[("Passive stock",lambda p:('passive',22000,p['cs'])),
       ("Semi-active hybrid",lambda p:('hyb',*hb,4500)),
       ("LQR comfort-first",lambda p:('llqr',Kc)),
       ("LQG comfort-first",lambda p:('llqg',Kc,L,Cm,Dm))]
out={c[0]:dict(aw=[],tyre=[],trav3=[],sat=[],stable=[]) for c in ctrls}
for i in range(N):
    p=dict(ev); p.update(ms=S['ms'][i],mu=S['mu'][i],kt=S['kt'][i],cs=S['cs'][i])
    for name,mk in ctrls:
        c=mk(p); Y=simulate3(p,c,w); s=summarize(Y,sk); o=out[name]
        o['aw'].append(s['aw']); o['tyre'].append(s['tyre']); o['trav3'].append(3000*s['trav']); o['sat'].append(100*s['satfrac'] if c[0] in('llqr','llqg') else 0)
        if c[0]=='llqr': st=max(np.real(cl_eigs_lag(p,Kc)))<0
        elif c[0]=='llqg': st=max(np.real(cl_eigs_lag(p,Kc,L,Cm,Dm)))<0
        else: st=True
        o['stable'].append(bool(st))
    if i%10==9: print(i+1,flush=True)
# nominal references on same road
ref={}
for name,p,c in (("Normal car nominal",ice,('passive',22000,1500)),("EV passive nominal",ev,('passive',22000,1500))):
    s=summarize(simulate3(p,c,w),sk); ref[name]=dict(aw=s['aw'],tyre=s['tyre'])
json.dump(dict(samples={k:list(map(float,v)) for k,v in S.items()},out=out,ref=ref),open('out/phase5.json','w'))
for name in out:
    o=out[name]; a=np.array(o['aw']); t=np.array(o['tyre'])
    print("%-20s aw median %.3f [p5 %.3f p95 %.3f max %.3f] tyre median %.4f [p95 %.4f max %.4f] stable %d/%d maxsat %.2f%%"%(
        name,np.median(a),np.percentile(a,5),np.percentile(a,95),a.max(),np.median(t),np.percentile(t,95),t.max(),sum(o['stable']),N,max(o['sat'])))
print(ref)
