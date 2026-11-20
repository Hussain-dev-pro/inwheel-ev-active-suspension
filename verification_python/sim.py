import numpy as np
from scipy.linalg import expm, solve_continuous_are
from core import params, build, wk_mag, road_psd, GD0, g
DT = 1e-3
TAU = 1/(2*np.pi*30)      # actuator lag (30 Hz bandwidth)
UMAX = 2500.0             # actuator force limit [N]

def plant(p, cs=None):
    A,Bw,_ = build(p, None, cs)
    ms,mu,kt = p['ms'],p['mu'],p['kt']
    Bu = np.array([[0],[1/ms],[0],[-1/mu]])
    Cz = np.vstack([A[1,:],[1,0,0,0],[0,0,kt/((ms+mu)*g),0]]); Dz = np.array([[1/ms],[0],[0]])
    return A,Bu,Bw,Cz,Dz

def zoh(A,B,dt=DT):
    n,m = B.shape; M=np.zeros((n+m,n+m)); M[:n,:n]=A*dt; M[:n,n:]=B*dt; E=expm(M)
    return E[:n,:n], E[:n,n:]

def care_ham(A,B,Q,R,N):
    """CARE with cross term via Hamiltonian eigenvectors (same method as MATLAB port)."""
    Ri=np.linalg.inv(R); At=A-B@Ri@N.T; Qt=Q-N@Ri@N.T
    H=np.block([[At,-B@Ri@B.T],[-Qt,-At.T]]); w,V=np.linalg.eig(H)
    Vs=V[:,np.real(w)<0]; n=A.shape[0]; X=np.real(Vs[n:]@np.linalg.inv(Vs[:n])); X=(X+X.T)/2
    return X

def lqr_design(p, qa, qd, qt, r):
    A,Bu,Bw,Cz,Dz = plant(p)
    Qz=np.diag([qa,qd,qt]); Q=Cz.T@Qz@Cz; N=Cz.T@Qz@Dz; R=Dz.T@Qz@Dz+r
    X=care_ham(A,Bu,Q,R,N); K=np.linalg.solve(R,Bu.T@X+N.T)
    Xs=solve_continuous_are(A,Bu,Q,R,s=N); Ks=np.linalg.solve(R,Bu.T@Xs+N.T)
    return K, np.max(np.abs(K-Ks))/np.max(np.abs(Ks))

def kalman_design(p, v=80/3.6, cls='C', sig_a=0.05, sig_d=0.5e-3, scale=100.0):
    """Continuous steady-state Kalman. Process noise: road velocity (white for ISO 8608, w=2)."""
    A,Bu,Bw,Cz,Dz = plant(p)
    q = 2*np.pi**2*GD0[cls]*0.1**2*v          # two-sided intensity of road velocity
    Cm = Cz[:2]; 
    Rn = np.diag([sig_a**2*DT, sig_d**2*DT])*scale  # sampled noise intensity x conservatism factor (limits noise feed-through)
    Y = care_ham(A.T, Cm.T, Bw@Bw.T*q, Rn, np.zeros((4,2)))
    L = Y@Cm.T@np.linalg.inv(Rn)
    return L, Cm, Dz[:2]

def road_random(T, v=80/3.6, cls='C', seed=0, fmin=0.5, fmax=80, nk=3000):
    rng=np.random.default_rng(seed); fk=np.linspace(fmin,fmax,nk); df=fk[1]-fk[0]
    amp=np.sqrt(2*road_psd(fk,v,cls)*df); ph=rng.uniform(0,2*np.pi,nk)
    t=np.arange(0,T,DT); zr=np.zeros_like(t); w=np.zeros_like(t)
    for i in range(0,nk,300):
        s=slice(i,i+300); arg=2*np.pi*np.outer(fk[s],t)+ph[s,None]
        zr+=amp[s]@np.sin(arg); w+=(amp[s]*2*np.pi*fk[s])@np.cos(arg)
    return t,zr,w

def road_bump(v_kmh=20, h=0.08, L=0.6, T=2.0, t0=0.2):
    v=v_kmh/3.6; t=np.arange(0,T,DT); x=(t-t0)*v; on=(x>=0)&(x<=L)
    zr=np.where(on,h*np.sin(np.pi*x/L)**2,0.0)              # smooth (raised-cosine) hump
    w=np.where(on,h*np.pi*v/L*np.sin(2*np.pi*x/L),0.0)
    return t,zr,w

def simulate(p, ctrl, w, seed=1, lag=True, sat=True):
    """ctrl: ('passive',) | ('sky',csky,cmin,cmax) | ('lqr',K) | ('lqg',K,L,Cm,Dm)"""
    kind=ctrl[0]
    cs_pass = 0.0 if kind=='sky' else p['cs']
    A,Bu,Bw,Cz,Dz = plant(p, cs=cs_pass)
    # augmented with actuator lag state
    Aa=np.zeros((5,5)); Aa[:4,:4]=A; Aa[:4,4:5]=Bu; Aa[4,4]=-1/TAU if lag else 0
    Ba=np.zeros((5,2)); Ba[4,0]=1/TAU if lag else 0; Ba[:4,1:2]=Bw
    if not lag: Aa[4,4]=0
    Ad,Bd=zoh(Aa,Ba)
    if kind=='lqg':
        K,L,Cm,Dm=ctrl[1:]
        Ae=A-L@Cm; Be=np.hstack([Bu-L@Dm, L]); Aed,Bed=zoh(Ae,Be)
        xh=np.zeros((4,1)); rng=np.random.default_rng(seed); sig=np.array([[0.05],[0.5e-3]])
    x=np.zeros((5,1)); n=len(w); Y=np.zeros((n,5))   # acc, defl, tyre, force, power
    for k in range(n):
        xs=x[:4]; ua=x[4,0]
        if kind=='passive': uc=0.0
        elif kind=='sky':
            csky,cmin,cmax=ctrl[1:]; vs=xs[1,0]; vr=xs[1,0]-xs[3,0]
            c = np.clip(csky*vs/vr,cmin,cmax) if (vs*vr>0 and abs(vr)>1e-9) else cmin
            uc=-c*vr; ua=uc            # ideal variable damper: no lag
        elif kind=='lqr': uc=(-(ctrl[1]@xs)).item()
        else:
            uc=(-(K@xh)).item()
        if sat and kind!='sky': uc=np.clip(uc,-UMAX,UMAX)
        if kind=='sky' or not lag: ua=uc; x[4,0]=uc
        acc=(Cz[0]@xs).item()+ua/p["ms"]
        Y[k]=[acc,xs[0,0],(Cz[2]@xs).item(),ua,ua*(xs[1,0]-xs[3,0])]
        if kind=='lqg':
            ym=np.array([[acc],[xs[0,0]]])+sig*rng.standard_normal((2,1))
            xh=Aed@xh+Bed@np.vstack([[[ua]],ym])
        x=Ad@x+Bd@np.array([[uc],[w[k]]])
        if kind=='sky' or not lag: x[4,0]=0.0
    return Y

def wrms(acc, skip=0):
    a=acc[skip:]; n=len(a); F=np.fft.rfftfreq(n,DT); S=np.fft.rfft(a)
    band=(F>=0.5)&(F<=80); W=np.zeros_like(F); W[band]=wk_mag(F[band])
    # Parseval on one-sided spectrum
    p2=np.abs(S*W)**2; p2[1:-1]*=2
    return np.sqrt(p2.sum()/n**2)

def summarize(Y, skip):
    a,d,t,u,P=[Y[skip:,i] for i in range(5)]
    return dict(aw=wrms(Y[:,0],skip), trav=d.std(), tyre=t.std(), umax=np.abs(u).max(), urms=u.std(),
                Pavg=P.mean(), Pabs=np.abs(P).mean(), satfrac=float(np.mean(np.abs(u)>=UMAX-1e-6)))

LABELS=[(0.315,"not uncomfortable"),(0.63,"a little uncomfortable"),(1.0,"fairly uncomfortable"),
        (1.6,"uncomfortable"),(2.5,"very uncomfortable"),(np.inf,"extremely uncomfortable")]
def iso_label(aw):
    for th,l in LABELS:
        if aw<th: return l

# ---------- frequency-shaped LQR (ISO Wk filter inside the design) ----------
def biquad(b2,b1,b0,a1,a0):
    A=np.array([[0,1],[-a0,-a1]]); B=np.array([[0],[1]]); C=np.array([[b0-b2*a0, b1-b2*a1]]); D=np.array([[b2]])
    return A,B,C,D
def series(s1,s2):
    A1,B1,C1,D1=s1; A2,B2,C2,D2=s2
    A=np.block([[A1,np.zeros((A1.shape[0],A2.shape[0]))],[B2@C1,A2]]); B=np.vstack([B1,B2@D1])
    C=np.hstack([D2@C1,C2]); D=D2@D1; return A,B,C,D
def wk_ss():
    w1,w3,w4,w5,w6=[2*np.pi*f for f in (0.4,12.5,12.5,2.37,3.35)]
    Q1=1/np.sqrt(2);Q4=0.63;Q5=Q6=0.91
    Hh=biquad(1,0,0,w1/Q1,w1**2)
    Ht=biquad(0,w4**2/w3,w4**2,w4/Q4,w4**2)
    Hs=biquad(1,w5/Q5,w5**2,w6/Q6,w6**2)
    return series(series(Hh,Ht),Hs)

def shaped_plant(p):
    A,Bu,Bw,Cz,Dz=plant(p); Af,Bf,Cf,Df=wk_ss(); nf=Af.shape[0]
    Aa=np.block([[A,np.zeros((4,nf))],[Bf@Cz[:1],Af]]); Bua=np.vstack([Bu,Bf*Dz[0,0]]); Bwa=np.vstack([Bw,np.zeros((nf,1))])
    Za=np.vstack([np.hstack([Df@Cz[:1],Cf]), np.hstack([Cz[1:],np.zeros((2,nf))])]); Dza=np.vstack([Df*Dz[0,0],np.zeros((2,1))])
    return Aa,Bua,Bwa,Za,Dza

def slqr_design(p,qa,qd,qt,r):
    Aa,Bua,Bwa,Za,Dza=shaped_plant(p)
    Qz=np.diag([qa,qd,qt]); Q=Za.T@Qz@Za; N=Za.T@Qz@Dza; R=Dza.T@Qz@Dza+r
    X=care_ham(Aa,Bua,Q,R,N); K=np.linalg.solve(R,Bua.T@X+N.T)
    Xs=solve_continuous_are(Aa,Bua,Q,R,s=N); Ks=np.linalg.solve(R,Bua.T@Xs+N.T)
    return K, np.max(np.abs(K-Ks))/np.max(np.abs(Ks))

def cl_metrics(p, K, v=80/3.6, cls='C', pdes=None):
    """Exact freq-domain metrics of plant p under shaped-LQR gain K (designed on pdes; no saturation/lag)."""
    Aa,Bua,Bwa,Za,Dza=shaped_plant(p)
    Acl=Aa-Bua@K
    _,_,_,Cz,Dz=plant(p)
    Cacc=np.hstack([Cz[0:1],np.zeros((1,Aa.shape[0]-4))])-Dz[0,0]*K
    Cd=np.hstack([Cz[1:3],np.zeros((2,Aa.shape[0]-4))]); Cu=-K
    from core import F
    S=road_psd(F,v,cls); W=wk_mag(F); I=np.eye(Aa.shape[0])
    H=np.array([ (np.vstack([Cacc,Cd,Cu])@np.linalg.solve(1j*2*np.pi*f*I-Acl,Bwa)).ravel()*1j*2*np.pi*f for f in F]).T
    r=lambda m2: np.sqrt(np.trapezoid(m2*S,F))
    return dict(aw=r(np.abs(H[0])**2*W**2), trav=r(np.abs(H[1])**2), tyre=r(np.abs(H[2])**2), urms=r(np.abs(H[3])**2),
                stable=bool(np.all(np.real(np.linalg.eigvals(Acl))<0)))

def simulate2(p, ctrl, w, seed=1, sig=(0.05,0.5e-3)):
    """ctrl: ('passive',ks,cs) | ('hyb',alpha,csky,cgnd,cmin,cmax) | ('slqr',K) | ('slqg',K,L,Cm,Dm)
       Active: saturation +-UMAX, first-order lag TAU. Semi-active: ideal variable damper (passive spring kept)."""
    kind=ctrl[0]; pp=dict(p)
    if kind=='passive': pp['ks'],pp['cs']=ctrl[1],ctrl[2]
    if kind=='hyb': pp['cs']=0.0
    A,Bu,Bw,Cz,Dz=plant(pp)
    Aa=np.zeros((5,5)); Aa[:4,:4]=A; Aa[:4,4:5]=Bu; Aa[4,4]=-1/TAU
    Ba=np.zeros((5,2)); Ba[4,0]=1/TAU; Ba[:4,1:2]=Bw; Ad,Bd=zoh(Aa,Ba)
    Ap,Bp=zoh(np.hstack([A,np.zeros((4,0))]),np.hstack([Bu,Bw]))   # no-lag plant for semi-active/passive
    Af,Bf,Cf,Df=wk_ss(); Afd,Bfd=zoh(Af,Bf); xf=np.zeros((6,1))
    if kind=='slqg':
        K,L,Cm,Dm=ctrl[1:]; Ae=A-L@Cm; Be=np.hstack([Bu-L@Dm,L]); Aed,Bed=zoh(Ae,Be); xh=np.zeros((4,1))
        rng=np.random.default_rng(seed); sg=np.array([[sig[0]],[sig[1]]])
    x=np.zeros((5,1)); n=len(w); Y=np.zeros((n,5))
    for k in range(n):
        xs=x[:4]
        if kind in('passive','hyb'):
            if kind=='passive': ua=0.0
            else:
                al,csk,cg,cmin,cmax=ctrl[1:]; vs=xs[1,0]; vu=xs[3,0]; vr=vs-vu
                if abs(vr)<1e-9: c=cmin
                else:
                    c_s=csk*vs/vr if vs*vr>0 else 0.0
                    c_g=cg*(-vu)/vr if (-vu)*vr>0 else 0.0
                    c=np.clip(al*c_s+(1-al)*c_g,cmin,cmax)
                ua=-c*vr
            acc=(Cz[0]@xs).item()+ua/pp['ms']
            Y[k]=[acc,xs[0,0],(Cz[2]@xs).item(),ua,ua*(xs[1,0]-xs[3,0])]
            xn=Ap@xs+Bp@np.array([[ua],[w[k]]]); x[:4]=xn; continue
        ua=x[4,0]; acc=(Cz[0]@xs).item()+ua/pp['ms']
        if kind=='slqr':
            K=ctrl[1]; uc=(-(K@np.vstack([xs,xf]))).item()
            xf=Afd@xf+Bfd*acc
        else:
            uc=(-(K@np.vstack([xh,xf]))).item()
            acc_hat=(Cz[0]@xh).item()+ua/pp['ms']; xf=Afd@xf+Bfd*acc_hat
            ym=np.array([[acc],[xs[0,0]]])+sg*rng.standard_normal((2,1))
            xh=Aed@xh+Bed@np.vstack([[[ua]],ym])
        uc=float(np.clip(uc,-UMAX,UMAX))
        Y[k]=[acc,xs[0,0],(Cz[2]@xs).item(),ua,ua*(xs[1,0]-xs[3,0])]
        x=Ad@x+Bd@np.array([[uc],[w[k]]])
    return Y

def cl_eigs(p, K, L=None, Cm=None, Dm=None):
    """Closed-loop eigenvalues incl. actuator lag (+ estimator for LQG). Linear (no saturation)."""
    A,Bu,Bw,Cz,Dz=plant(p); Af,Bf,Cf,Df=wk_ss(); Kx,Kf=K[:,:4],K[:,4:]
    # states: x(4), ua(1), xf(6) [, xh(4)]
    if L is None:
        n=11; M=np.zeros((n,n))
        M[:4,:4]=A; M[:4,4:5]=Bu
        M[4,:4]=-(Kx)/TAU; M[4,4]=-1/TAU; M[4,5:11]=-(Kf)/TAU
        M[5:11,:4]=Bf@Cz[:1]; M[5:11,4:5]=Bf*Dz[0,0]; M[5:11,5:11]=Af
    else:
        n=15; M=np.zeros((n,n))
        M[:4,:4]=A; M[:4,4:5]=Bu
        M[4,4]=-1/TAU; M[4,5:11]=-Kf/TAU; M[4,11:15]=-Kx/TAU
        M[5:11,11:15]=Bf@Cz[:1]; M[5:11,4:5]=Bf*Dz[0,0]; M[5:11,5:11]=Af
        # xh' = A xh + Bu ua + L(Cm x + Dm ua - Cm xh - Dm ua)
        M[11:15,:4]=L@Cm; M[11:15,4:5]=Bu; M[11:15,11:15]=A-L@Cm
    return np.linalg.eigvals(M)

# ---------- lag-aware shaped LQR: states [x(4); ua; xf(6)], input = force command ----------
def lag_plant(p):
    A,Bu,Bw,Cz,Dz=plant(p); Af,Bf,Cf,Df=wk_ss(); nf=6; n=11
    M=np.zeros((n,n)); M[:4,:4]=A; M[:4,4:5]=Bu; M[4,4]=-1/TAU
    M[5:,:4]=Bf@Cz[:1]; M[5:,4:5]=Bf*Dz[0,0]; M[5:,5:]=Af
    Bc=np.zeros((n,1)); Bc[4,0]=1/TAU; Bwa=np.zeros((n,1)); Bwa[:4]=Bw
    Za=np.zeros((3,n)); Za[0,:4]=(Df@Cz[:1]).ravel(); Za[0,4]=Df[0,0]*Dz[0,0]; Za[0,5:]=Cf.ravel()
    Za[1,0]=1; Za[2,:4]=Cz[2]
    Cacc=np.zeros((1,n)); Cacc[0,:4]=Cz[0]; Cacc[0,4]=Dz[0,0]
    return M,Bc,Bwa,Za,Cacc

def lqr_lag(p,qa,qd,qt,r):
    M,Bc,Bwa,Za,_=lag_plant(p); Qz=np.diag([qa,qd,qt]); Q=Za.T@Qz@Za+1e-12*np.eye(11)
    X=care_ham(M,Bc,Q,np.array([[r]]),np.zeros((11,1))); K=Bc.T@X/r
    Xs=solve_continuous_are(M,Bc,Q,np.array([[r]])); Ks=Bc.T@Xs/r
    return K, float(np.max(np.abs(K-Ks))/np.max(np.abs(Ks)))

def cl_lag(p,K,v=80/3.6,cls='C'):
    M,Bc,Bwa,Za,Cacc=lag_plant(p); Acl=M-Bc@K
    from core import F
    C=np.vstack([Cacc,Za[1:3],np.eye(11)[4:5]])       # acc, defl, tyre, actual force
    S=road_psd(F,v,cls); W=wk_mag(F); I=np.eye(11)
    H=np.array([(C@np.linalg.solve(1j*2*np.pi*f*I-Acl,Bwa)).ravel()*1j*2*np.pi*f for f in F]).T
    r=lambda m2: np.sqrt(np.trapezoid(m2*S,F))
    return dict(aw=r(np.abs(H[0])**2*W**2),trav=r(np.abs(H[1])**2),tyre=r(np.abs(H[2])**2),urms=r(np.abs(H[3])**2),
                stable=bool(np.all(np.real(np.linalg.eigvals(Acl))<0)))

def cl_eigs_lag(p,K,L=None,Cm=None,Dm=None):
    M,Bc,Bwa,Za,Cacc=lag_plant(p)
    if L is None: return np.linalg.eigvals(M-Bc@K)
    A,Bu,Bw,Cz,Dz=plant(p); Kx=K[:,:4]; Kr=K[:,4:]
    n=15; G=np.zeros((n,n)); G[:11,:11]=M; G[:11,:11]-=0
    # command uses xh for plant states, true ua and xf (xf driven by estimated accel)
    G[:11,4:11]-=Bc@Kr; G[:11,11:15]-=Bc@Kx
    G[5:11,:4]=0; G[5:11,11:15]+=wk_ss()[1]@Cz[:1]
    G[11:15,:4]=L@Cm; G[11:15,4:5]=Bu; G[11:15,11:15]=A-L@Cm
    return np.linalg.eigvals(G)

def simulate3(p, ctrl, w, seed=1, sig=(0.05,0.5e-3)):
    """lag-aware LQR/LQG: ('llqr',K) | ('llqg',K,L,Cm,Dm); others delegated to simulate2"""
    kind=ctrl[0]
    if kind not in('llqr','llqg'): return simulate2(p,ctrl,w,seed,sig)
    A,Bu,Bw,Cz,Dz=plant(p)
    Aa=np.zeros((5,5)); Aa[:4,:4]=A; Aa[:4,4:5]=Bu; Aa[4,4]=-1/TAU
    Ba=np.zeros((5,2)); Ba[4,0]=1/TAU; Ba[:4,1:2]=Bw; Ad,Bd=zoh(Aa,Ba)
    Af,Bf,Cf,Df=wk_ss(); Afd,Bfd=zoh(Af,Bf); xf=np.zeros((6,1)); K=ctrl[1]
    if kind=='llqg':
        _,_,L,Cm,Dm=ctrl; Ae=A-L@Cm; Be=np.hstack([Bu-L@Dm,L]); Aed,Bed=zoh(Ae,Be); xh=np.zeros((4,1))
        rng=np.random.default_rng(seed); sg=np.array([[sig[0]],[sig[1]]])
    x=np.zeros((5,1)); n=len(w); Y=np.zeros((n,5))
    for k in range(n):
        xs=x[:4]; ua=x[4,0]; acc=(Cz[0]@xs).item()+ua/p['ms']
        xc = xs if kind=='llqr' else xh
        uc=float(np.clip((-(K@np.vstack([xc,[[ua]],xf]))).item(),-UMAX,UMAX))
        Y[k]=[acc,xs[0,0],(Cz[2]@xs).item(),ua,ua*(xs[1,0]-xs[3,0])]
        acc_c = acc if kind=='llqr' else (Cz[0]@xh).item()+ua/p['ms']
        xf=Afd@xf+Bfd*acc_c
        if kind=='llqg':
            ym=np.array([[acc],[xs[0,0]]])+sg*rng.standard_normal((2,1)); xh=Aed@xh+Bed@np.vstack([[[ua]],ym])
        x=Ad@x+Bd@np.array([[uc],[w[k]]])
    return Y
