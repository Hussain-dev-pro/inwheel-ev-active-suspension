%% PHASE 3 - Passive suspension: sweep, Pareto front, optimisation (normal car vs in-wheel EV)
%  Road: ISO 8608 Class C, 80 km/h.  Comfort: ISO 2631-1 Wk-weighted RMS acceleration.
%  Road holding: RMS dynamic tyre-load ratio.  Constraint: 3-sigma suspension travel <= 80 mm.
clear; clc; close all;
if ~exist('media','dir'), mkdir('media'); end
variants = {'ice','ev'};  label = {'Normal car','In-wheel EV (+30 kg/wheel)'};
blue = [0.165 0.471 0.839];  orange = [0.922 0.408 0.204];
cm = interp1([0 .33 .66 1]', [0.918 0.945 0.984; 0.612 0.753 0.933; 0.165 0.471 0.839; 0.071 0.216 0.420], linspace(0,1,64)');

%% 1) Sweep spring rate x damping  (ride frequency window ~1.0-1.6 Hz)
KS = linspace(15000, 43000, 57);  CS = linspace(400, 5000, 116);
for k = 1:2
    p = params(variants{k});
    G(k).aw = zeros(numel(KS),numel(CS));  G(k).trav = G(k).aw;  G(k).tyre = G(k).aw;
    for i = 1:numel(KS)
        for j = 1:numel(CS)
            m = ride_metrics(p, KS(i), CS(j));
            G(k).aw(i,j) = m.aw;  G(k).trav(i,j) = m.trav;  G(k).tyre(i,j) = m.tyre;
        end
    end
    stock(k) = ride_metrics(p);
end
fprintf('\nSTOCK TUNE (ks=22 kN/m, cs=1500 N s/m)\n');
for k = 1:2
    fprintf('%-28s weighted accel %.3f m/s^2   tyre-load %.4f   3-sigma travel %.1f mm\n', label{k}, stock(k).aw, stock(k).tyre, 3000*stock(k).trav);
end
fprintf('EV vs normal car: comfort %+.1f %%, road holding %+.1f %%\n', ...
    (stock(2).aw/stock(1).aw-1)*100, (stock(2).tyre/stock(1).tyre-1)*100);

%% 2) Pareto fronts (comfort vs road holding, travel-feasible designs)
travLim = 0.08/3;
for k = 1:2
    feas = G(k).trav <= travLim;
    pts = [G(k).aw(feas), G(k).tyre(feas)];
    [~,ix] = sort(pts(:,1));  pts = pts(ix,:);
    front = zeros(0,2);  best = inf;
    for i = 1:size(pts,1)
        if pts(i,2) < best - 1e-12, front(end+1,:) = pts(i,:); best = pts(i,2); end %#ok<SAGROW>
    end
    PF{k} = front;
    ok = feas;  t = G(k).tyre;  t(~ok) = inf;  [tmin,ii] = min(t(:));
    [r,c] = ind2sub(size(t), ii);
    fprintf('%-28s best possible road holding %.4f  (ks=%.0f, cs=%.0f, weighted accel %.3f)\n', ...
        label{k}, tmin, KS(r), CS(c), G(k).aw(r,c));
end

%% 3) fminsearch optimisation at three road-holding tolerances vs the normal car's stock value
fprintf('\nOPTIMUM PASSIVE TUNE: minimise weighted accel, tyre-load <= (1+tol) x normal-car stock\n');
fprintf('%-6s %-28s %8s %8s %8s %8s %9s\n','tol','vehicle','ks','cs','accel','tyre','3s-travel');
for tol = [0.05 0.10 0.20]
    for k = 1:2
        [x,m] = optimise_design(params(variants{k}), stock(1).tyre*(1+tol));
        fprintf('+%-4.0f%% %-28s %8.0f %8.0f %8.3f %8.4f %7.1f mm\n', tol*100, label{k}, x(1), x(2), m.aw, m.tyre, 3000*m.trav);
        OPT.(sprintf('%s_tol%d', variants{k}, round(tol*100))) = struct('x',x,'m',m);
    end
end

%% 4) Figures
fig1 = figure('Position',[100 100 1100 760]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
fld = {'aw','tyre'};  lvl = {0.6:0.1:2.0, 0.20:0.02:0.46};
ttl = {'Ride harshness: ISO 2631 weighted accel (m/s^2)','Road-holding loss: RMS tyre-load ratio'};
for r = 1:2
    for k = 1:2
        nexttile; contourf(CS, KS/1000, G(k).(fld{r}), lvl{r}, 'LineColor','none'); hold on;
        [c,hh] = contour(CS, KS/1000, G(k).(fld{r}), lvl{r}(1:2:end), 'w'); clabel(c,hh,'FontSize',7,'Color','w');
        plot(1500, 22, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerEdgeColor','k');
        colormap(gca, cm); caxis(lvl{r}([1 end])); colorbar;
        title([label{k} ' - ' ttl{r}]); if r == 2, xlabel('Damping c_s (N s/m)'); end
        if k == 1, ylabel('Spring rate k_s (kN/m)'); end
    end
end
try exportgraphics(fig1, 'media/phase3_carpet_plots.png', 'Resolution', 200); catch, saveas(fig1, 'media/phase3_carpet_plots.png'); end

fig2 = figure('Position',[100 100 860 540]); hold on; grid on;
plot(PF{1}(:,1), PF{1}(:,2), 'Color', blue,   'LineWidth', 2.2);
plot(PF{2}(:,1), PF{2}(:,2), 'Color', orange, 'LineWidth', 2.2);
plot(stock(1).aw, stock(1).tyre, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerEdgeColor', blue);
plot(stock(2).aw, stock(2).tyre, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerEdgeColor', orange);
yline(stock(1).tyre, '--k');  xlim([0.7 1.9]);  ylim([0.2 0.42]);
xlabel('Ride harshness: ISO 2631 weighted RMS body acceleration (m/s^2)  -> worse');
ylabel('Road-holding loss: RMS dynamic tyre-load ratio  -> worse');
title('Best achievable trade-off, passive suspension (ISO 8608 Class C, 80 km/h)');
legend({'Normal car','In-wheel EV','Normal car stock tune','EV stock tune'}, 'Location','northeast');
try exportgraphics(fig2, 'media/phase3_pareto.png', 'Resolution', 200); catch, saveas(fig2, 'media/phase3_pareto.png'); end

save('phase3_results.mat', 'KS','CS','G','PF','stock','OPT');
fprintf('\nSaved media/phase3_carpet_plots.png, media/phase3_pareto.png, phase3_results.mat\n');
