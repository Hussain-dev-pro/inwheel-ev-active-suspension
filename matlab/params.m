function p = params(variant)
%PARAMS  Quarter-car parameters: generic C-segment compact car, FRONT axle corner.
%   p = params('ice')  normal combustion/EV-skateboard car with inboard motor
%   p = params('ev')   same car with an in-wheel motor (+30 kg unsprung per wheel)
%
%   Sources / basis (see README):
%     ms  : measured front-corner weight ~398 kg (NHTSA test, 2020 C-segment sedan,
%           front axle 798 kg) minus ~40 kg unsprung  ->  ~350 kg
%     mu  : engineering estimate: wheel+tyre ~21, hub/knuckle ~9, brake ~10, share of arms ~3
%     ks  : wheel rate sized for ~1.2 Hz body bounce (typical front 1.0-1.5 Hz)
%     kt  : typical passenger-car tyre, ASSUMED 200 kN/m (literature range ~135-250 kN/m)
%     cs  : gives damping ratio ~0.23-0.26 (typical 0.2-0.4)
%     EV  : +30 kg/wheel (Lotus Focus test; Protean Pd18 adds 36 kg/corner)
if nargin < 1, variant = 'ice'; end
p.ms = 350;      % sprung mass per corner          [kg]
p.mu = 40;       % unsprung mass per corner        [kg]
p.ks = 22000;    % suspension spring (wheel rate)  [N/m]
p.cs = 1500;     % suspension damper               [N s/m]
p.kt = 200000;   % tyre vertical stiffness         [N/m]
p.g  = 9.81;     % gravity                         [m/s^2]
switch lower(variant)
    case 'ice'
    case 'ev'
        p.mu = p.mu + 30;
    otherwise
        error('params:variant', 'variant must be ''ice'' or ''ev''');
end
p.variant = lower(variant);
end
