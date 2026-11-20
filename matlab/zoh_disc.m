function [Ad,Bd] = zoh_disc(A,B,dt)
%ZOH_DISC  Exact zero-order-hold discretisation (core MATLAB expm).
n=size(A,1); m=size(B,2); E=expm([A B; zeros(m,n+m)]*dt); Ad=E(1:n,1:n); Bd=E(1:n,n+1:end);
end
