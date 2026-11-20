function Y = simulate_car(p, ctrl, w, opt)
%SIMULATE_CAR  Time-domain simulation, 1 kHz, exact ZOH discretisation.
%  ctrl.type: 'passive' (ctrl.ks, ctrl.cs) | 'hybrid' (alpha,csky,cgnd,cmin,cmax) |
%             'lqr' (K) | 'lqg' (K,L,Cm,Dm)
%  Active: force limit opt.umax, actuator lag opt.tau.  LQG: sensor noise opt.sig.
%  Y columns: body accel, suspension defl, tyre-load ratio, actuator force, actuator power
dt=opt.dt; n=numel(w); Y=zeros(n,5); pp=p;
switch ctrl.type
    case 'passive', pp.ks=ctrl.ks; pp.cs=ctrl.cs;
    case 'hybrid',  pp.cs=0;
end
[A,Bu,Bw,C,Du] = build_model(pp);
if any(strcmp(ctrl.type,{'passive','hybrid'}))
    [Ad,Bd]=zoh_disc(A,[Bu Bw],dt); x=zeros(4,1);
    for k=1:n
        ua=0;
        if strcmp(ctrl.type,'hybrid')
            vs=x(2); vu=x(4); vr=vs-vu; c=ctrl.cmin;
            if abs(vr)>1e-9
                cs_=0; cg=0;
                if vs*vr>0,    cs_=ctrl.csky*vs/vr;  end
                if -vu*vr>0,   cg=ctrl.cgnd*(-vu)/vr; end
                c=min(max(ctrl.alpha*cs_+(1-ctrl.alpha)*cg, ctrl.cmin), ctrl.cmax);
            end
            ua=-c*vr;
        end
        Y(k,:)=[C(1,:)*x+Du(1)*ua, x(1), C(3,:)*x, ua, ua*(x(2)-x(4))];
        x=Ad*x+Bd*[ua; w(k)];
    end
    return
end
% active: plant + actuator lag
Aa=zeros(5); Aa(1:4,1:4)=A; Aa(1:4,5)=Bu; Aa(5,5)=-1/opt.tau;
Ba=zeros(5,2); Ba(5,1)=1/opt.tau; Ba(1:4,2)=Bw; [Ad,Bd]=zoh_disc(Aa,Ba,dt);
[Af,Bf] = wk_ss(); [Afd,Bfd]=zoh_disc(Af,Bf,dt); xf=zeros(6,1); K=ctrl.K;
isg = strcmp(ctrl.type,'lqg');
if isg
    [Aed,Bed]=zoh_disc(A-ctrl.L*ctrl.Cm, [Bu-ctrl.L*ctrl.Dm, ctrl.L], dt); xh=zeros(4,1);
    rs=RandStream('mt19937ar','Seed',opt.seed);
end
x=zeros(5,1);
for k=1:n
    xs=x(1:4); ua=x(5); acc=C(1,:)*xs+Du(1)*ua;
    if isg, xc=xh; else, xc=xs; end
    uc=min(max(-K*[xc; ua; xf], -opt.umax), opt.umax);
    Y(k,:)=[acc, xs(1), C(3,:)*xs, ua, ua*(xs(2)-xs(4))];
    if isg
        xf=Afd*xf+Bfd*(C(1,:)*xh+Du(1)*ua);
        ym=[acc; xs(1)] + opt.sig(:).*randn(rs,2,1);
        xh=Aed*xh+Bed*[ua; ym];
    else
        xf=Afd*xf+Bfd*acc;
    end
    x=Ad*x+Bd*[uc; w(k)];
end
end
