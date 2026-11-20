function m = cl_lag_metrics(p, K, tau, v, cls)
%CL_LAG_METRICS  Exact frequency-domain metrics with the lag-aware LQR in the loop (no saturation).
[M,Bc,Bw,Z,Cacc] = lag_plant(p, tau); Acl=M-Bc*K;
f=logspace(log10(0.5),log10(80),400); Cy=[Cacc; Z(2:3,:); [0 0 0 0 1 zeros(1,6)]];
H=zeros(4,numel(f));
for i=1:numel(f), s=1i*2*pi*f(i); H(:,i)=(Cy*((s*eye(11)-Acl)\Bw))*s; end
S=iso8608_psd(f,v,cls); W=wk_weight(f); r=@(m2) sqrt(trapz(f,m2.*S));
m.aw=r(abs(H(1,:)).^2.*W.^2); m.trav=r(abs(H(2,:)).^2); m.tyre=r(abs(H(3,:)).^2); m.urms=r(abs(H(4,:)).^2);
m.stable = all(real(eig(Acl))<0);
end
