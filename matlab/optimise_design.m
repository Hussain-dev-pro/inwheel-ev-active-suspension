function [x, m] = optimise_design(p, tyreLimit)
%OPTIMISE_DESIGN  Best passive (ks, cs): minimise weighted body accel subject to
%   RMS tyre-load ratio <= tyreLimit  and  3-sigma suspension travel <= 80 mm.
%   Uses fminsearch (core MATLAB) from several starts, with box bounds via a clip + penalty.
%   x = [ks cs],  bounds: ks 15-43 kN/m (ride frequency ~1.0-1.6 Hz), cs 400-5000 N s/m.
lo = [15000 400];  hi = [43000 5000];
starts = [0.1 0.2; 0.3 0.3; 0.5 0.5; 0.2 0.6; 0.05 0.1];
opts = optimset('TolX',1e-5, 'TolFun',1e-8, 'MaxIter',600, 'Display','off');
best = inf;  zb = starts(1,:);
for k = 1:size(starts,1)
    [z, fv] = fminsearch(@(z) cost(z, p, tyreLimit, lo, hi), starts(k,:), opts);
    if fv < best, best = fv; zb = z; end
end
x = lo + min(max(zb,0),1).*(hi - lo);
m = ride_metrics(p, x(1), x(2));
end

function J = cost(z, p, L, lo, hi)
x = lo + min(max(z,0),1).*(hi - lo);
m = ride_metrics(p, x(1), x(2));
J = m.aw + 1e3*max(0, m.tyre/L - 1)^2 + 1e3*max(0, 3*m.trav/0.08 - 1)^2;
end
