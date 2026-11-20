function [L,Cm,Dm] = kalman_design(p, v, cls, sig, dt, scale)
%KALMAN_DESIGN  Steady-state Kalman gain. Sensors: body accelerometer + suspension-travel sensor.
%  Process noise = road velocity, which is WHITE for an ISO 8608 road (PSD ~ 1/n^2).
%  scale (default 100) makes the filter trust the sensors less, so sensor noise does not
%  leak into the actuator force (checked: 17x less >80 Hz jitter, same comfort).
if nargin < 6, scale = 100; end
[A,Bu,Bw,C,Du] = build_model(p);
Gd0 = struct('A',16e-6,'B',64e-6,'C',256e-6,'D',1024e-6);
q  = 2*pi^2*Gd0.(cls)*0.1^2*v;                 % two-sided intensity of road velocity
Cm = C(1:2,:); Dm = Du(1:2);
Rn = diag(sig.^2*dt)*scale;                      % sampled sensor noise intensity x conservatism
Y  = care_ham(A', Cm', Bw*Bw'*q, Rn, zeros(4,2));
L  = Y*Cm'/Rn;
end
