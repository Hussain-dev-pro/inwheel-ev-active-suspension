function [t,zr,w] = road_bump(v_kmh, h, L, T, dt)
%ROAD_BUMP  Smooth (raised-cosine) speed hump of height h [m], length L [m], at v_kmh.
v=v_kmh/3.6; t=0:dt:T-dt; x=(t-0.2)*v; on=(x>=0)&(x<=L);
zr=zeros(size(t)); w=zeros(size(t));
zr(on)=h*sin(pi*x(on)/L).^2; w(on)=h*pi*v/L*sin(2*pi*x(on)/L);
end
