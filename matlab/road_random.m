function [t,zr,w] = road_random(T, v, cls, seed, dt)
%ROAD_RANDOM  ISO 8608 random road as a sum of sinusoids (0.5-80 Hz). Returns height zr and velocity w.
rng(seed); nk=3000; fk=linspace(0.5,80,nk); df=fk(2)-fk(1);
amp=sqrt(2*iso8608_psd(fk,v,cls)*df); ph=2*pi*rand(1,nk);
t=0:dt:T-dt; zr=zeros(size(t)); w=zeros(size(t));
for i=1:300:nk
    j=i:min(i+299,nk); arg=2*pi*fk(j)'*t+ph(j)';
    zr=zr+amp(j)*sin(arg); w=w+(amp(j).*2*pi.*fk(j))*cos(arg);
end
end
