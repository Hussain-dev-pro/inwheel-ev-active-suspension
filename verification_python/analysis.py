import numpy as np, json, sys
sys.path.insert(0,'.')
from core import *
from scipy.linalg import expm
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BLUE, ORANGE = "#2a78d6", "#eb6834"; INK="#0b0b0b"; INK2="#52514e"; SURF="#fcfcfb"; GRIDC="#e4e3df"
plt.rcParams.update({"font.family":"DejaVu Sans","figure.facecolor":SURF,"axes.facecolor":SURF,"axes.edgecolor":"#b8b7b0",
  "axes.labelcolor":INK2,"xtick.color":INK2,"ytick.color":INK2,"text.color":INK,"axes.grid":True,"grid.color":GRIDC,
  "grid.linewidth":0.8,"axes.spines.top":False,"axes.spines.right":False,"axes.titlesize":11,"axes.titleweight":"bold",
  "axes.titlelocation":"left","font.size":9.5})
R = {}
ice, ev = params(False), params(True)

# ---------------- PHASE 1 ----------------
fi, zi = freqs(ice); fe, ze = freqs(ev)
R['modes'] = dict(ICE=dict(f=fi.tolist(), zeta=zi.tolist()), EV=dict(f=fe.tolist(), zeta=ze.tolist()))
def step(p, h=0.05, T=2.5, n=1500):
    A,Bw,C = build(p); t=np.linspace(0,T,n); x0=np.array([0,0,-h,0.])
    X = np.array([expm(A*tk)@x0 for tk in t]).T
    zs = X[0]+X[2]+h; acc = (A[1,:]@X); tyre=X[2]
    return t, zs, acc, tyre
out={}
for n,p in (("ICE",ice),("EV",ev)):
    t,zs,acc,ty = step(p); h=0.05
    over = (zs.max()-h)/h*100
    band = np.abs(zs-h)>0.02*h; ts = t[np.where(band)[0].max()+1] if band.any() else 0
    out[n]=dict(overshoot_pct=over, settle2pct_s=ts, peak_acc=float(np.abs(acc).max()), peak_tyre_defl_mm=float(np.abs(ty).max()*1000))
R['step']=out
# Fig 1
fig,ax=plt.subplots(1,2,figsize=(11,4.1))
def wheel(p,h=0.05,T=0.6,n=900):
    A,Bw,C=build(p); t=np.linspace(0,T,n); x0=np.array([0,0,-h,0.])
    X=np.array([expm(A*tk)@x0 for tk in t]).T
    return t,(X[2]+h)*1000
w1=wheel(ice); w2=wheel(ev)
ax[0].plot(w1[0],w1[1],color=BLUE,lw=2); ax[0].plot(w2[0],w2[1],color=ORANGE,lw=2)
ax[0].axhline(50,color=INK2,lw=0.8,ls=(0,(4,3)))
ax[0].set_title("Wheel response to a 50 mm road step"); ax[0].set_xlabel("Time (s)"); ax[0].set_ylabel("Wheel displacement (mm)")
ax[0].text(0.30,66,"Normal car: rings at 11.7 Hz",color=BLUE,fontsize=9.5,fontweight="bold")
ax[0].text(0.30,61,"In-wheel EV: rings at 8.9 Hz, slower to settle",color=ORANGE,fontsize=9.5,fontweight="bold")
ax[0].text(0.30,45,"Body overshoot is the same for both (55%)",color=INK2,fontsize=8.5)
ax[0].set_xlim(0,0.6)
def Hacc(p):
    A,Bw,C=build(p); w=2*np.pi*F
    return np.array([abs((C[0]@np.linalg.solve(1j*wi*np.eye(4)-A,Bw)).item()*1j*wi) for wi in w])
ax[1].loglog(F,Hacc(ice),color=BLUE,lw=2); ax[1].loglog(F,Hacc(ev),color=ORANGE,lw=2)
ax[1].set_title("Body acceleration per metre of road height"); ax[1].set_xlabel("Frequency (Hz)"); ax[1].set_ylabel("|body accel / road height|  (m/s² per m)")
for f,lab,c,dy in ((fi[0],"body bounce\n%.2f Hz"%fi[0],INK2,0),(fi[1],"wheel hop\n%.1f Hz (normal)"%fi[1],BLUE,0),(fe[1],"wheel hop\n%.1f Hz (EV)"%fe[1],ORANGE,0)):
    ax[1].axvline(f,color=c,lw=0.8,ls=(0,(4,3)))
ax[1].text(fi[0]*1.06,3,"body bounce\n%.2f Hz"%fi[0],color=INK2,fontsize=8.5,va="bottom")
ax[1].text(fe[1]*0.93,6,"EV hop\n%.1f Hz"%fe[1],color=ORANGE,fontsize=8.5,ha="right",fontweight="bold")
ax[1].text(fi[1]*1.06,6,"normal car hop\n%.1f Hz"%fi[1],color=BLUE,fontsize=8.5,fontweight="bold")
ax[1].set_xlim(0.5,80); ax[1].set_ylim(2,3e4)
fig.tight_layout(); fig.savefig("out/phase1_ice_vs_ev.png",dpi=200); plt.close(fig)

# ---------------- PHASE 3 ----------------
KS=np.linspace(15000,43000,57); CS=np.linspace(400,5000,116)
G={}
for n,p in (("ICE",ice),("EV",ev)):
    d={k:np.zeros((KS.size,CS.size)) for k in ('aw','trav','tyre')}
    for i,ks in enumerate(KS):
        for j,cs in enumerate(CS):
            m=metrics(p,ks,cs)
            for k in d: d[k][i,j]=m[k]
    G[n]=d
base={n:metrics(p) for n,p in (("ICE",ice),("EV",ev))}
R['stock']={n:{k:float(v) for k,v in base[n].items()} for n in base}
trav_lim=0.08/3
def pareto(d):
    feas=d['trav']<=trav_lim; pts=np.column_stack([d['aw'][feas],d['tyre'][feas]]); pts=pts[np.argsort(pts[:,0])]
    fr=[];b=np.inf
    for a,t in pts:
        if t<b-1e-12: fr.append((a,t)); b=t
    return np.array(fr)
FR={n:pareto(G[n]) for n in G}
lo=np.array([15000,400.]); hi=np.array([43000,5000.])
def solve(p,L):
    def f(z):
        z=np.clip(z,0,1); x=lo+z*(hi-lo); m=metrics(p,x[0],x[1])
        return m['aw']+1e3*max(0,m['tyre']/L-1)**2+1e3*max(0,3*m['trav']/0.08-1)**2
    best=None
    for z0 in ([0.1,0.2],[0.3,0.3],[0.5,0.5],[0.2,0.6],[0.05,0.1]):
        r=minimize(f,z0,method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-8,maxiter=600))
        if best is None or r.fun<best.fun: best=r
    x=lo+np.clip(best.x,0,1)*(hi-lo); m=metrics(p,x[0],x[1]); return x,m
opt={}
for tol in (0.05,0.10,0.20):
    L=base['ICE']['tyre']*(1+tol)
    for n,p in (("ICE",ice),("EV",ev)):
        x,m=solve(p,L)
        opt[f"{n}_tol{int(tol*100)}"]=dict(ks=float(x[0]),cs=float(x[1]),aw=float(m['aw']),tyre=float(m['tyre']),trav3s_mm=float(m['trav']*3000),
              f_body=float(freqs(p,x[0],x[1])[0][0]),zeta=float(freqs(p,x[0],x[1])[1][0]))
R['opt']=opt
# lowest achievable tyre (grip-only) per vehicle
for n in G:
    ok=G[n]['trav']<=trav_lim; t=np.where(ok,G[n]['tyre'],np.inf); i,j=np.unravel_index(t.argmin(),t.shape)
    R[f'gripbest_{n}']=dict(ks=float(KS[i]),cs=float(CS[j]),tyre=float(G[n]['tyre'][i,j]),aw=float(G[n]['aw'][i,j]))
    # best grip at stock-ICE comfort level
    ok2=ok&(G[n]['aw']<=base['ICE']['aw']); t=np.where(ok2,G[n]['tyre'],np.inf); i,j=np.unravel_index(t.argmin(),t.shape)
    R[f'grip_at_iceComfort_{n}']=dict(ks=float(KS[i]),cs=float(CS[j]),tyre=float(G[n]['tyre'][i,j]),aw=float(G[n]['aw'][i,j]))
json.dump(R,open("out/results.json","w"),indent=1)

# Fig 2 carpet plots
cmap=LinearSegmentedColormap.from_list("bl",["#eaf1fb","#9cc0ee","#2a78d6","#12376b"])
fig,axs=plt.subplots(2,2,figsize=(11,7.6),sharex=True,sharey=True)
for _a in axs.ravel(): _a.grid(False)
KK,CC=np.meshgrid(CS,KS/1000)   # x=cs, y=ks
for r,(key,ttl,unit,lv) in enumerate((("aw","Ride harshness (ISO 2631 weighted accel)","m/s²",np.arange(0.6,2.0001,0.1)),("tyre","Road-holding loss (dynamic tyre load ratio, RMS)","–",np.arange(0.20,0.4601,0.02)))):
    for c,n in enumerate(("ICE","EV")):
        ax=axs[r,c]; Z=G[n][key]
        cf=ax.contourf(CS,KS/1000,Z,levels=lv,cmap=cmap,extend="both")
        cs_=ax.contour(CS,KS/1000,Z,levels=lv[::2],colors="white",linewidths=0.6)
        ax.clabel(cs_,fmt="%.2f",fontsize=7,inline=True,inline_spacing=3)
        ax.scatter([ice['cs']],[ice['ks']/1000],s=55,facecolor=SURF,edgecolor=INK,linewidths=1.8,zorder=5)
        ax.set_title(("Normal car" if n=="ICE" else "In-wheel EV (+30 kg/wheel)")+" – "+ttl.split(" (")[0],fontsize=10)
        if r==1: ax.set_xlabel("Damping c_s (N·s/m)")
        if c==0: ax.set_ylabel("Spring rate k_s (kN/m)")
    fig.colorbar(cf,ax=axs[r,:],shrink=0.9,pad=0.02,label=ttl.split("(")[1].rstrip(")")+(" ["+unit+"]"))
axs[0,0].text(0.02,0.04,"○ stock tune",transform=axs[0,0].transAxes,fontsize=8.5,color=INK)
fig.savefig("out/phase3_carpet_plots.png",dpi=200,bbox_inches="tight"); plt.close(fig)

# Fig 3 Pareto
fig,ax=plt.subplots(figsize=(8.6,5.4))
for n,c,lab in (("ICE",BLUE,"Normal car"),("EV",ORANGE,"In-wheel EV")):
    F_=FR[n]; ax.plot(F_[:,0],F_[:,1],color=c,lw=2.2); ax.plot(base[n]['aw'],base[n]['tyre'],"o",ms=9,mfc=SURF,mec=c,mew=2.2,zorder=5)
ax.set_xlabel("Ride harshness: ISO 2631 weighted RMS body acceleration (m/s²)   → worse"); ax.set_ylabel("Road-holding loss: RMS dynamic tyre load ratio   → worse")
ax.set_title("Best achievable trade-off, passive suspension (ISO 8608 Class C road, 80 km/h)")
ax.text(0.80,0.40,"Normal car",color=BLUE,fontweight="bold",fontsize=10); ax.text(1.55,0.40,"In-wheel EV",color=ORANGE,fontweight="bold",fontsize=10)
ax.annotate("stock tune",(base['ICE']['aw'],base['ICE']['tyre']),xytext=(0.78,0.212),fontsize=8.5,color=INK2,arrowprops=dict(arrowstyle="-",color=INK2,lw=0.8))
ax.annotate("stock tune",(base['EV']['aw'],base['EV']['tyre']),xytext=(1.45,0.285),fontsize=8.5,color=INK2,arrowprops=dict(arrowstyle="-",color=INK2,lw=0.8))
ax.axhline(base['ICE']['tyre'],color=INK2,lw=0.8,ls=(0,(4,3)))
ax.text(0.72,base['ICE']['tyre']+0.003,"normal car stock road-holding",fontsize=8,color=INK2)
ax.set_xlim(0.7,1.9); ax.set_ylim(0.2,0.42)
fig.tight_layout(); fig.savefig("out/phase3_pareto.png",dpi=200); plt.close(fig)
print(json.dumps(R,indent=1))
