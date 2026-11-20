import sys; sys.path.insert(0,'.')
import numpy as np, json
from core import *; from sim import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Polygon
from PIL import Image
ev=params(True); D=json.load(open('out/lqr_designs.json'))
Kc,_=lqr_lag(ev,1,1/0.03**2,10**D['lag_comfort']['lqt'],10**D['lag_comfort']['lr']); L,Cm,Dm=kalman_design(ev)
vk,h,Lh,T=30,0.075,3.7,2.6; v=vk/3.6
t,zr,w=road_bump(vk,h,Lh,T)
res=[]
for c in (('passive',22000,1500),('llqg',Kc,L,Cm,Dm)):
    Y=simulate3(ev,c,w); x3=Y[:,2]*(ev['ms']+ev['mu'])*g/ev['kt']; zu=zr+x3; zs=zu+Y[:,1]; res.append((zs,zu,Y[:,0]))
ORANGE,AQUA,INK2="#eb6834","#1baf7a","#52514e"; names=["EV, normal suspension","EV, smart (active) suspension"]; cols=[ORANGE,AQUA]
k=2; Rw=0.30; gap=0.45; bw=1.3; bh=0.32
road=lambda p: h*np.sin(np.pi*np.clip(p,0,Lh)/Lh)**2
def frame(n, fig=None):
    fig=plt.figure(figsize=(11,6.2),dpi=90,facecolor="#fcfcfb")
    fig.suptitle("75 mm speed hump at 30 km/h   ·   slow motion 0.5×   ·   vertical motion drawn ×2",fontsize=10.5,color=INK2,y=0.98)
    s=(t[n]-0.2)*v; px=np.linspace(-2.2,2.2,200)
    for j in range(2):
        ax=fig.add_axes([0.04+0.48*j,0.47,0.44,0.44]); ax.set_xlim(-2.2,2.2); ax.set_ylim(-0.05,1.55); ax.set_aspect('equal'); ax.axis('off')
        zs,zu,acc=res[j]; yu=Rw+k*zu[n]; yb=Rw+gap+k*zs[n]
        ax.fill(np.r_[px,2.2,-2.2],np.r_[k*road(px+s),-0.05,-0.05],color="#8c8c8c")
        ax.add_patch(Circle((0,yu),Rw,color="#262626")); ax.add_patch(Circle((0,yu),0.08,color="#b3b3b3"))
        yy=np.linspace(yu+0.08,yb,17); xx=-0.25+np.r_[0,0.09*(-1)**np.arange(1,16),0]; ax.plot(xx,yy,color="#1f3a5f",lw=2)
        mid=(yu+yb)/2; ax.add_patch(Polygon([[0.19,yu+0.08],[0.31,yu+0.08],[0.31,mid+0.04],[0.19,mid+0.04]],fc="#c8d9ee",ec="#1f3a5f"))
        ax.plot([0.25,0.25],[mid-0.02,yb],color="#1f3a5f",lw=2.5)
        ax.add_patch(Rectangle((-bw/2,yb),bw,bh,fc=cols[j],ec="k"))
        ax.plot([-0.9,0.9],[Rw+gap+bh/2]*2,':',color="#999")
        ax.set_title(names[j],color=cols[j],fontsize=12.5,fontweight="bold")
        ax.text(-2.15,1.42,"shake: %+.1f m/s²"%acc[n],color=cols[j],fontsize=10.5,fontweight="bold")
    axb=fig.add_axes([0.08,0.15,0.88,0.27]); axb.set_facecolor("#fcfcfb")
    for j in range(2): axb.plot(t[:n+1],res[j][2][:n+1],color=cols[j],lw=2,label=names[j])
    axb.set_xlim(0,T); axb.set_ylim(-7,5); axb.grid(color="#e4e3df"); axb.set_xlabel("Time (s)"); axb.set_ylabel("Body acceleration (m/s²)")
    axb.legend(loc="lower right",frameon=False,fontsize=9); axb.set_title("How hard the passengers get shaken (less = more comfortable)",fontsize=10.5,loc="left",color=INK2)
    for sp in ("top","right"): axb.spines[sp].set_visible(False)
    return fig
frames=[]
for n in range(0,len(t),20):
    f=frame(n); f.canvas.draw(); im=Image.frombuffer('RGBA',f.canvas.get_width_height(),f.canvas.buffer_rgba()).convert('RGB'); frames.append(im); plt.close(f)
frames[0].save('out/hump_animation.gif',save_all=True,append_images=frames[1:],duration=40,loop=0,optimize=True)
# filmstrip: 4 key moments
picks=[(0.30,"a) hitting the hump"),(0.58,"b) coming down"),(0.92,"c) after the hump"),(1.50,"d) settling")]
strip=[]
for tt,lab in picks:
    n=int(round(tt/DT)); f=frame(n); f.text(0.5,0.005,lab+"  (t = %.2f s)"%tt,ha="center",fontsize=11,fontweight="bold"); f.canvas.draw()
    strip.append(Image.frombuffer('RGBA',f.canvas.get_width_height(),f.canvas.buffer_rgba()).convert('RGB')); plt.close(f)
W,H=strip[0].size; m=Image.new('RGB',(W*2,H*2),'white')
for i,im in enumerate(strip): m.paste(im,((i%2)*W,(i//2)*H))
m.save('out/animation_frames.png'); print(len(frames),"frames; peak", [round(float(np.abs(r[2]).max()),2) for r in res])
