function m = ride_metrics(p, ks, cs, v, cls)
%RIDE_METRICS  Frequency-domain ride/handling metrics on an ISO 8608 random road.
%   m = ride_metrics(p)                       stock tune, Class C road, 80 km/h
%   m = ride_metrics(p, ks, cs, v, cls)       any spring/damper, speed v [m/s], class 'A'..'D'
%   m.aw   ISO 2631 Wk-weighted RMS body acceleration   [m/s^2]  (ride harshness)
%   m.au   unweighted RMS body acceleration             [m/s^2]
%   m.trav RMS suspension deflection                    [m]      (rattle space use)
%   m.tyre RMS dynamic tyre load / static load          [-]      (road-holding loss)
%   Closed-form transfer functions from road height (equal to the state-space model;
%   verified in selfcheck.m). Exact, fast, and needs no toolbox.
if nargin < 2 || isempty(ks), ks = p.ks; end
if nargin < 3 || isempty(cs), cs = p.cs; end
if nargin < 4 || isempty(v),  v  = 80/3.6; end
if nargin < 5 || isempty(cls), cls = 'C'; end
ms = p.ms;  mu = p.mu;  kt = p.kt;
f  = logspace(log10(0.5), log10(80), 400);          % ISO 2631-1 band
s  = 1i*2*pi*f;
k   = cs*s + ks;
den = k*ms.*s.^2 + (k + ms*s.^2).*(mu*s.^2 + kt);
D   = -kt*ms*s.^2 ./ den;                           % suspension deflection / road height
acc = -k.*D/ms;                                     % body acceleration  / road height
Zu  = -D.*(k + ms*s.^2) ./ (ms*s.^2);               % wheel displacement / road height
tyr = kt/((ms+mu)*p.g) * (Zu - 1);                  % dyn. tyre load ratio / road height
S   = iso8608_psd(f, v, cls);
W   = wk_weight(f);
rmsf = @(mag2) sqrt(trapz(f, mag2.*S));
m.aw   = rmsf(abs(acc).^2 .* W.^2);
m.au   = rmsf(abs(acc).^2);
m.trav = rmsf(abs(D).^2);
m.tyre = rmsf(abs(tyr).^2);
end
