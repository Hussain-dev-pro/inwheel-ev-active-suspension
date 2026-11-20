%% PHASE 4 - Control design: semi-active and active suspension for the in-wheel-motor EV
%  Controllers: passive (stock / re-tuned), semi-active skyhook, semi-active hybrid sky-groundhook,
%  frequency-shaped lag-aware LQR (comfort-first, grip-first), LQG (2 noisy sensors + Kalman).
%  Real-world limits: force +-2500 N, 30 Hz actuator lag, sensor noise, actuator power.
clear; clc; close all; if ~exist('media','dir'), mkdir('media'); end
opt.dt=1e-3; opt.tau=1/(2*pi*30); opt.umax=2500; opt.sig=[0.05 0.5e-3]; opt.seed=1;
ev=params('ev'); ice=params('ice'); v=80/3.6;

%% 1) Invariant point: no force between body and wheel changes body accel at sqrt(kt/mu)
f0 = sqrt(ev.kt/ev.mu)/(2*pi);
fprintf('Tyre-only frequency of the EV wheel: %.2f Hz  (inside the 4-8 Hz human-sensitive band edge)\n', f0);

%% 2) Designs (weights found by optimisation; see README). qd = 1/0.03^2 on travel, qa = 1.
Kc = lqr_lag(ev, 1, 1/0.03^2, 10^1.30, 1e-9, opt.tau);   % comfort-first: tyre-load <= EV stock
Kg = lqr_lag(ev, 1, 1/0.03^2, 10^1.60, 1e-9, opt.tau);   % grip-first:    tyre-load <= normal-car stock
[L,Cm,Dm] = kalman_design(ev, v, 'C', opt.sig, opt.dt);
mc = cl_lag_metrics(ev,Kc,opt.tau,v,'C');  mg = cl_lag_metrics(ev,Kg,opt.tau,v,'C');
fprintf('\nExact design check (frequency domain, no saturation)\n');
fprintf('LQR comfort-first: aw %.3f  tyre %.4f  force rms %.0f N  stable %d\n', mc.aw, mc.tyre, mc.urms, mc.stable);
fprintf('LQR grip-first   : aw %.3f  tyre %.4f  force rms %.0f N  stable %d\n', mg.aw, mg.tyre, mg.urms, mg.stable);

%% 3) Time-domain bench: same 60 s ISO 8608 Class C road at 80 km/h + 75 mm x 3.7 m speed hump at 30 km/h
[~,~,w]  = road_random(60, v, 'C', 0, opt.dt);  skip=5000;
[tb,zb,wb] = road_bump(30, 0.075, 3.7, 2.5, opt.dt);
C = {};
C(end+1,:) = {'Normal car, passive stock',     ice, struct('type','passive','ks',22000,'cs',1500)};
C(end+1,:) = {'EV, passive stock',             ev,  struct('type','passive','ks',22000,'cs',1500)};
C(end+1,:) = {'EV, passive re-tuned',          ev,  struct('type','passive','ks',15000,'cs',1565)};
C(end+1,:) = {'EV, semi-active skyhook',       ev,  struct('type','hybrid','alpha',1,'csky',6000,'cgnd',0,'cmin',300,'cmax',4500)};
C(end+1,:) = {'EV, semi-active hybrid',        ev,  struct('type','hybrid','alpha',0.5,'csky',3000,'cgnd',3000,'cmin',300,'cmax',4500)};
C(end+1,:) = {'EV, active LQR comfort-first',  ev,  struct('type','lqr','K',Kc)};
C(end+1,:) = {'EV, active LQG comfort-first',  ev,  struct('type','lqg','K',Kc,'L',L,'Cm',Cm,'Dm',Dm)};
C(end+1,:) = {'EV, active LQR grip-first',     ev,  struct('type','lqr','K',Kg)};
fprintf('\n%-30s %7s %-22s %7s %8s %7s %7s %6s | %8s %6s %7s\n','case','aw','ISO 2631 label','tyre','3s trav','F rms','F max','P(W)','hump acc','tyre','travel');
for i=1:size(C,1)
    Y = simulate_car(C{i,2}, C{i,3}, w, opt);  s = td_metrics(Y, skip, opt.dt, opt.umax);
    Yb= simulate_car(C{i,2}, C{i,3}, wb, opt); B{i}=Yb; %#ok<SAGROW>
    T(i)=s; %#ok<SAGROW>
    fprintf('%-30s %7.3f %-22s %7.4f %6.1fmm %7.0f %7.0f %6.0f | %8.2f %6.2f %5.1fmm\n', C{i,1}, s.aw, iso_label(s.aw), s.tyre, 3000*s.trav, ...
        s.urms, s.umax, s.P, max(abs(Yb(:,1))), max(Yb(:,3)), 1000*max(abs(Yb(:,2))));
end

%% 4) Figures
col=[0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478];
fig=figure('Position',[100 100 1100 420]); idx=[1 2 7];
subplot(1,2,1); hold on; grid on;
for j=1:3, plot(tb, B{idx(j)}(:,1), 'Color', col(j,:), 'LineWidth', 1.8); end
xlabel('Time (s)'); ylabel('Body acceleration (m/s^2)'); title('75 mm speed hump at 30 km/h'); legend(C(idx,1),'Location','southeast');
subplot(1,2,2); hold on; grid on;
for j=1:3, plot(tb, B{idx(j)}(:,3), 'Color', col(j,:), 'LineWidth', 1.8); end
xlabel('Time (s)'); ylabel('Dynamic tyre-load ratio (-)'); title('Tyre load (1 = wheel lifts off)');
try exportgraphics(fig,'media/phase4_hump.png','Resolution',200); catch, saveas(fig,'media/phase4_hump.png'); end
save('phase4_results.mat','C','T','Kc','Kg','L','Cm','Dm','opt');
fprintf('\nSaved media/phase4_hump.png and phase4_results.mat\n');
