function [M,Bc,Bw,Z,Cacc] = lag_plant(p, tau)
%LAG_PLANT  Design model for active control: quarter car + actuator lag + Wk filter.
%  States [z_s-z_u; dz_s; z_u-z_r; dz_u; actuator force; 6 Wk filter states]
%  Input  force COMMAND (actuator follows with first-order lag tau)
%  Z = [Wk-weighted body accel; suspension deflection; tyre-load ratio]
[A,Bu,Bwp,C,Du] = build_model(p);
[Af,Bf,Cf,Df] = wk_ss();
M = zeros(11); M(1:4,1:4)=A; M(1:4,5)=Bu; M(5,5)=-1/tau;
M(6:11,1:4)=Bf*C(1,:); M(6:11,5)=Bf*Du(1); M(6:11,6:11)=Af;
Bc = zeros(11,1); Bc(5)=1/tau;  Bw = zeros(11,1); Bw(1:4)=Bwp;
Z = zeros(3,11); Z(1,1:4)=Df*C(1,:); Z(1,5)=Df*Du(1); Z(1,6:11)=Cf; Z(2,1)=1; Z(3,1:4)=C(3,:);
Cacc = zeros(1,11); Cacc(1:4)=C(1,:); Cacc(5)=Du(1);
end
