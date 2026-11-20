function [A,Bu,Bw,C,Du] = build_model(p, ks, cs)
%BUILD_MODEL  State-space quarter-car model.
%   States  x = [ z_s-z_u ; dz_s ; z_u-z_r ; dz_u ]   (rattle space, body vel, tyre defl, wheel vel)
%   Inputs  u = actuator force (active, Phase 4);  w = road velocity dz_r/dt
%   Outputs y = [ body accel (m/s^2) ; suspension deflection (m) ; dynamic tyre load ratio (-) ]
%   dx/dt = A x + Bu u + Bw w ;   y = C x + Du u
if nargin < 2 || isempty(ks), ks = p.ks; end
if nargin < 3 || isempty(cs), cs = p.cs; end
ms = p.ms;  mu = p.mu;  kt = p.kt;
A  = [ 0      1      0     -1;
      -ks/ms -cs/ms  0      cs/ms;
       0      0      0      1;
       ks/mu  cs/mu -kt/mu -cs/mu ];
Bu = [0; 1/ms; 0; -1/mu];
Bw = [0; 0; -1; 0];
C  = [ A(2,:);
       1 0 0 0;
       0 0 kt/((ms+mu)*p.g) 0 ];
Du = [1/ms; 0; 0];
end
