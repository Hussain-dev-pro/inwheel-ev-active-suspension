%% SELFCHECK  Run this first. Verifies the building blocks before trusting any result.
%  1) Wk weighting vs ISO 2631-1 table      2) closed-form vs state-space model
%  3) stock-car metrics vs independently validated reference values
clear; clc;

% --- 1) Wk weighting -----------------------------------------------------------
fT = [1 1.25 1.6 2 2.5 3.15 4 5 6.3 8 10 12.5 16 20 25 31.5 40 50];
wT = [0.482 0.484 0.494 0.531 0.631 0.804 0.967 1.039 1.054 1.036 0.988 0.902 ...
      0.768 0.636 0.513 0.405 0.314 0.246];
err = max(abs(wk_weight(fT)./wT - 1))*100;
fprintf('1) Wk weighting: max error vs ISO 2631-1 table = %.2f %%   %s\n', err, pf(err < 1));

% --- 2) closed form vs state-space ---------------------------------------------
for v = {'ice','ev'}
    p = params(v{1});
    a = ride_metrics(p);
    b = metrics_ss(p);
    d = max(abs([a.aw a.trav a.tyre]./[b.aw b.trav b.tyre] - 1));
    fprintf('2) %-3s closed-form vs state-space: max rel. diff = %.1e   %s\n', upper(v{1}), d, pf(d < 1e-6));
end

% --- 3) reference values (validated against a time-domain random-road simulation) --
ref = struct('ice', [1.1475 0.012035 0.2217], 'ev', [1.2203 0.012493 0.2539]);
for v = {'ice','ev'}
    m = ride_metrics(params(v{1}));
    d = max(abs([m.aw m.trav m.tyre]./ref.(v{1}) - 1));
    fprintf('3) %-3s stock: aw=%.4f  trav=%.5f m  tyre=%.4f   diff vs reference = %.2f %%   %s\n', ...
        upper(v{1}), m.aw, m.trav, m.tyre, d*100, pf(d < 0.005));
end

function s = pf(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end

function m = metrics_ss(p)   % independent state-space evaluation (matrix solve per frequency)
[A,~,Bw,C] = build_model(p);
f = logspace(log10(0.5), log10(80), 400);
H = zeros(3, numel(f));
for i = 1:numel(f)
    w = 2*pi*f(i);
    H(:,i) = (C*((1i*w*eye(4) - A)\Bw))*(1i*w);
end
S = iso8608_psd(f, 80/3.6, 'C');  W = wk_weight(f);
r = @(mag2) sqrt(trapz(f, mag2.*S));
m.aw = r(abs(H(1,:)).^2.*W.^2);  m.trav = r(abs(H(2,:)).^2);  m.tyre = r(abs(H(3,:)).^2);
end
