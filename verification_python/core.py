import numpy as np
g = 9.81

def params(ev=False, dm_ev=30.0):
    p = dict(ms=350.0, mu=40.0, ks=22000.0, cs=1500.0, kt=200000.0)
    if ev: p['mu'] += dm_ev
    return p

def build(p, ks=None, cs=None):
    ms, mu, kt = p['ms'], p['mu'], p['kt']
    ks = p['ks'] if ks is None else ks
    cs = p['cs'] if cs is None else cs
    A = np.array([[0,1,0,-1],
                  [-ks/ms,-cs/ms,0,cs/ms],
                  [0,0,0,1],
                  [ks/mu,cs/mu,-kt/mu,-cs/mu]])
    Bw = np.array([[0],[0],[-1],[0]])
    C = np.vstack([A[1,:], [1,0,0,0], [0,0,kt/((ms+mu)*g),0]])
    return A, Bw, C

def wk_mag(f):
    s = 1j*2*np.pi*np.asarray(f, float)
    f1,f2,f3,f4,f5,f6 = 0.4,100.0,12.5,12.5,2.37,3.35
    w1,w2,w3,w4,w5,w6=[2*np.pi*x for x in (f1,f2,f3,f4,f5,f6)]
    Q1=Q2=1/np.sqrt(2);Q4=0.63;Q5=Q6=0.91
    Hh = s**2/(s**2+w1/Q1*s+w1**2)
    Hl = w2**2/(s**2+w2/Q2*s+w2**2)
    Ht = (w4**2/w3)*(s+w3)/(s**2+w4/Q4*s+w4**2)
    Hs = (s**2+w5/Q5*s+w5**2)/(s**2+w6/Q6*s+w6**2)
    return np.abs(Hh*Hl*Ht*Hs)

GD0 = {'A':16e-6,'B':64e-6,'C':256e-6,'D':1024e-6}   # ISO 8608 Gd(n0=0.1 c/m) [m^3]
def road_psd(f, v, cls='C'):
    """one-sided displacement PSD [m^2/Hz] vs temporal freq f [Hz]; Gd(n)=Gd0 (n/n0)^-2"""
    n0 = 0.1
    return GD0[cls]*n0**2*v/np.asarray(f)**2

F = np.logspace(np.log10(0.5), np.log10(80), 400)   # ISO 2631-1 band
def metrics(p, ks=None, cs=None, v=80/3.6, cls='C'):
    A, Bw, C = build(p, ks, cs)
    w = 2*np.pi*F
    I = np.eye(4)
    H = np.empty((3, F.size), complex)
    for i, wi in enumerate(w):
        x = np.linalg.solve(1j*wi*I - A, Bw)        # state per unit road VELOCITY
        H[:, i] = (C @ x).ravel() * (1j*wi)          # per unit road DISPLACEMENT
    S = road_psd(F, v, cls)
    def rms(mag2): return np.sqrt(np.trapezoid(mag2*S, F))
    a_w  = rms(np.abs(H[0])**2 * wk_mag(F)**2)   # ISO 2631 weighted body accel [m/s^2]
    a_u  = rms(np.abs(H[0])**2)                  # unweighted
    trav = rms(np.abs(H[1])**2)                  # rms suspension deflection [m]
    tyre = rms(np.abs(H[2])**2)                  # rms dynamic tyre load ratio [-]
    return dict(aw=a_w, au=a_u, trav=trav, tyre=tyre)

def freqs(p, ks=None, cs=None):
    A,_,_ = build(p, ks, cs)
    ev = np.linalg.eigvals(A)
    ev = ev[np.imag(ev) > 0]
    wn = np.abs(ev); zeta = -np.real(ev)/wn
    idx = np.argsort(wn)
    return wn[idx]/(2*np.pi), zeta[idx]
