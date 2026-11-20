%% ANIMATE_HUMP - Side-by-side animation over a speed hump:
%  in-wheel-motor EV with normal suspension  vs  same EV with smart (LQG active) suspension.
%  Saves media/hump_animation.gif (and a video file if this MATLAB supports it).
clear; clc; close all;
if ~exist('media','dir'), mkdir('media'); end
opt.dt=1e-3; opt.tau=1/(2*pi*30); opt.umax=2500; opt.sig=[0.05 0.5e-3]; opt.seed=1;
ev = params('ev');
Kc = lqr_lag(ev, 1, 1/0.03^2, 10^1.30, 1e-9, opt.tau);
[L,Cm,Dm] = kalman_design(ev, 80/3.6, 'C', opt.sig, opt.dt);

% speed hump: 75 mm high, 3.7 m long, at 30 km/h
vkmh = 30; h = 0.075; Lh = 3.7; T = 2.6; v = vkmh/3.6;
[t,zr,w] = road_bump(vkmh, h, Lh, T, opt.dt);  t = t(:); zr = zr(:);
cars  = {struct('type','passive','ks',22000,'cs',1500), ...
         struct('type','lqg','K',Kc,'L',L,'Cm',Cm,'Dm',Dm)};
names = {'EV - normal suspension', 'EV - smart (active) suspension'};
col   = [0.922 0.408 0.204; 0.106 0.686 0.478];
for j = 1:2
    Y = simulate_car(ev, cars{j}, w, opt);
    x3 = Y(:,3)*(ev.ms+ev.mu)*ev.g/ev.kt;      % tyre deflection z_u - z_r
    zu{j} = zr + x3;  zs{j} = zu{j} + Y(:,2);  acc{j} = Y(:,1);  F{j} = Y(:,4); %#ok<SAGROW>
end

%% drawing set-up (motion exaggerated x2 so it is easy to see)
k = 2;  Rw = 0.30;  gap = 0.45;  bw = 1.3;  bh = 0.32;
road = @(p) h*sin(pi*min(max(p,0),Lh)/Lh).^2;
fig = figure('Color','w','Position',[80 60 1250 700]);
for j = 1:2
    ax(j) = subplot(2,2,j); hold(ax(j),'on'); axis(ax(j),'equal'); axis(ax(j),[-2.2 2.2 -0.05 1.55]); %#ok<SAGROW>
    set(ax(j),'XTick',[],'YTick',[],'Box','on'); title(ax(j), names{j}, 'Color', col(j,:), 'FontSize', 13);
    G(j).road  = fill(ax(j), nan, nan, [0.55 0.55 0.55], 'EdgeColor', 'none');
    G(j).wheel = rectangle(ax(j),'Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[0.15 0.15 0.15],'EdgeColor','k');
    G(j).hub   = rectangle(ax(j),'Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
    G(j).spring= plot(ax(j), nan, nan, 'Color', [0.12 0.23 0.37], 'LineWidth', 2);
    G(j).dcyl  = fill(ax(j), nan, nan, [0.75 0.85 0.95], 'EdgeColor', [0.12 0.23 0.37], 'LineWidth', 1.2);
    G(j).drod  = plot(ax(j), nan, nan, 'Color', [0.12 0.23 0.37], 'LineWidth', 2.5);
    G(j).body  = fill(ax(j), nan, nan, col(j,:), 'EdgeColor', 'k', 'LineWidth', 1.2);
    G(j).txt   = text(ax(j), -2.1, 1.45, '', 'FontSize', 11, 'FontWeight', 'bold');
    G(j).ref   = plot(ax(j), [-0.9 0.9], (Rw+gap+bh/2)*[1 1], ':', 'Color', [0.4 0.4 0.4]);
end
axb = subplot(2,2,[3 4]); hold(axb,'on'); grid(axb,'on');
for j = 1:2, G(j).trace = plot(axb, nan, nan, 'Color', col(j,:), 'LineWidth', 2); end
xlim(axb,[0 T]); ylim(axb,[-7 5]); xlabel(axb,'Time (s)'); ylabel(axb,'Body acceleration (m/s^2)');
legend(axb, names, 'Location', 'southeast'); title(axb, 'How hard the passengers get shaken  (less = more comfortable)');
annotation(fig,'textbox',[0.01 0.955 0.98 0.04],'String', ...
    sprintf('75 mm speed hump at %d km/h  |  slow motion 0.5x  |  vertical motion drawn x%d', vkmh, k), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',11);

%% video outputs
gifFile = fullfile('media','hump_animation.gif');  vw = [];
try
    vw = VideoWriter(fullfile('media','hump_animation.mp4'),'MPEG-4');
catch
    try, vw = VideoWriter(fullfile('media','hump_animation.avi'),'Motion JPEG AVI'); catch, vw = []; end
end
if ~isempty(vw), vw.FrameRate = 25; open(vw); end

%% animate: one frame every 20 ms of simulated time, played at 25 fps (0.5x)
idx = 1:20:numel(t);  first = true;
for n = idx
    s  = (t(n) - 0.2)*v;                        % distance travelled onto the hump
    px = linspace(-2.2, 2.2, 200);
    for j = 1:2
        yu = Rw + k*zu{j}(n);  yb = Rw + gap + k*zs{j}(n);
        set(G(j).road, 'XData', [px 2.2 -2.2], 'YData', [k*road(px+s) -0.05 -0.05]);
        set(G(j).wheel,'Position', [-Rw yu-Rw 2*Rw 2*Rw]);
        set(G(j).hub,  'Position', [-0.08 yu-0.08 0.16 0.16]);
        [sx,sy] = zigzag(-0.25, yu+0.08, yb, 7, 0.09);  set(G(j).spring, 'XData', sx, 'YData', sy);
        mid = (yu+yb)/2;
        set(G(j).dcyl, 'XData', 0.25+[-0.06 0.06 0.06 -0.06], 'YData', [yu+0.08 yu+0.08 mid+0.04 mid+0.04]);
        set(G(j).drod, 'XData', [0.25 0.25], 'YData', [mid-0.02 yb]);
        set(G(j).body, 'XData', [-bw bw bw -bw]/2, 'YData', [yb yb yb+bh yb+bh]);
        set(G(j).trace,'XData', t(1:n), 'YData', acc{j}(1:n));
        set(G(j).txt,  'String', sprintf('shake: %+.1f m/s^2', acc{j}(n)), 'Color', col(j,:));
    end
    drawnow;
    fr = getframe(fig);
    if ~isempty(vw), writeVideo(vw, fr); end
    [im,cm] = rgb2ind(frame2im(fr), 128);
    if first, imwrite(im,cm,gifFile,'gif','LoopCount',inf,'DelayTime',0.04); first = false;
    else,     imwrite(im,cm,gifFile,'gif','WriteMode','append','DelayTime',0.04); end
end
if ~isempty(vw), close(vw); end
fprintf('Peak shake  - normal: %.2f m/s^2   smart: %.2f m/s^2  (%.0f%% less)\n', ...
    max(abs(acc{1})), max(abs(acc{2})), 100*(1-max(abs(acc{2}))/max(abs(acc{1}))));
fprintf('Saved %s\n', gifFile);

function [x,y] = zigzag(x0, y0, y1, n, wdt)
yy = linspace(y0, y1, 2*n+3);
xx = x0 + [0, wdt*(-1).^(1:2*n+1), 0];
x = xx;  y = yy;
end
