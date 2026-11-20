%% PHASE 1 - Quarter-car model: normal car vs in-wheel-motor EV
%  Builds the 2-DOF model, extracts natural frequencies / damping ratios,
%  and compares the step response and frequency response of both vehicles.
clear; clc; close all;
if ~exist('media','dir'), mkdir('media'); end
variants = {'ice','ev'};  label = {'Normal car','In-wheel EV'};  col = [0.165 0.471 0.839; 0.922 0.408 0.204];

fprintf('\nPHASE 1 - natural frequencies (stock suspension)\n');
fprintf('%-12s %6s %6s | %12s %10s | %12s %10s\n','vehicle','ms','mu','body bounce','zeta','wheel hop','zeta');
for k = 1:2
    p = params(variants{k});  A = build_model(p);
    ev = eig(A);  ev = ev(imag(ev) > 0);  [wn,ix] = sort(abs(ev));  ev = ev(ix);
    fn = wn/(2*pi);  zeta = -real(ev)./wn;
    fprintf('%-12s %6.0f %6.0f | %9.2f Hz %10.3f | %9.2f Hz %10.3f\n', label{k}, p.ms, p.mu, fn(1), zeta(1), fn(2), zeta(2));
    R(k).fn = fn;  R(k).zeta = zeta;  R(k).p = p;  R(k).A = A;
end

%% Step response: 50 mm road step (free response from x0; no toolbox needed)
h = 0.05;  x0 = [0; 0; -h; 0];
tB = linspace(0, 2.5, 1500);  tW = linspace(0, 0.6, 900);
fprintf('\nStep response to a %.0f mm road step\n', h*1000);
for k = 1:2
    XB = zeros(4,numel(tB));  for i = 1:numel(tB), XB(:,i) = expm(R(k).A*tB(i))*x0; end
    zs = XB(1,:) + XB(3,:) + h;                         % body displacement
    over = (max(zs) - h)/h*100;
    out = find(abs(zs - h) > 0.02*h, 1, 'last');  ts = tB(min(out+1, numel(tB)));
    XW = zeros(4,numel(tW));  for i = 1:numel(tW), XW(:,i) = expm(R(k).A*tW(i))*x0; end
    R(k).wheel = (XW(3,:) + h)*1000;                    % wheel displacement [mm]
    fprintf('%-12s body overshoot %.1f %%   2%% settling time %.2f s\n', label{k}, over, ts);
end

%% Frequency response: body acceleration per metre of road height
f = logspace(log10(0.5), log10(80), 400);  s = 1i*2*pi*f;
for k = 1:2
    p = R(k).p;  ks = p.ks;  cs = p.cs;  ms = p.ms;  mu = p.mu;  kt = p.kt;
    kk = cs*s + ks;  den = kk*ms.*s.^2 + (kk + ms*s.^2).*(mu*s.^2 + kt);
    R(k).Hacc = abs(-kk.*(-kt*ms*s.^2./den)/ms);
end

fig = figure('Position',[100 100 1100 420]);
subplot(1,2,1); hold on; grid on;
for k = 1:2, plot(tW, R(k).wheel, 'Color', col(k,:), 'LineWidth', 2); end
yline(50, '--k'); xlabel('Time (s)'); ylabel('Wheel displacement (mm)');
title('Wheel response to a 50 mm road step'); legend(label, 'Location','northeast');
subplot(1,2,2);
loglog(f, R(1).Hacc, 'Color', col(1,:), 'LineWidth', 2); hold on; grid on;
loglog(f, R(2).Hacc, 'Color', col(2,:), 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('|body accel / road height|  (m/s^2 per m)');
title('Body acceleration per metre of road height'); legend(label, 'Location','southwest');
try exportgraphics(fig, 'media/phase1_ice_vs_ev.png', 'Resolution', 200); catch, saveas(fig, 'media/phase1_ice_vs_ev.png'); end
save('phase1_results.mat', 'R');
fprintf('\nSaved media/phase1_ice_vs_ev.png and phase1_results.mat\n');
