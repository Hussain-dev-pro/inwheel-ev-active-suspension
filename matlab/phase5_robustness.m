%% PHASE 5 - Robustness: Monte Carlo over load, unsprung mass, tyre pressure and damper wear.
%  Controller gains stay FIXED at their nominal design (as on a real car).
clear; clc; close all; if ~exist('media','dir'), mkdir('media'); end
load('phase4_results.mat','Kc','L','Cm','Dm','opt');
ev=params('ev'); v=80/3.6; N=100;                 % reduce N to save MATLAB Online hours
[~,~,w]=road_random(20, v, 'C', 7, opt.dt); skip=3000;
rng(42);
ms=350*(0.9+0.4*rand(N,1)); mu=70*(0.85+0.3*rand(N,1)); kt=2e5*(0.8+0.4*rand(N,1)); cs=1500*(0.8+0.2*rand(N,1));
names={'Passive stock','Semi-active hybrid','LQR comfort-first','LQG comfort-first'};
AW=zeros(N,4); TY=AW; ST=true(N,4);
for i=1:N
    p=ev; p.ms=ms(i); p.mu=mu(i); p.kt=kt(i); p.cs=cs(i);
    ctr={struct('type','passive','ks',22000,'cs',cs(i)), ...
         struct('type','hybrid','alpha',0.5,'csky',3000,'cgnd',3000,'cmin',300,'cmax',4500), ...
         struct('type','lqr','K',Kc), struct('type','lqg','K',Kc,'L',L,'Cm',Cm,'Dm',Dm)};
    for j=1:4
        s=td_metrics(simulate_car(p,ctr{j},w,opt),skip,opt.dt,opt.umax); AW(i,j)=s.aw; TY(i,j)=s.tyre;
    end
    [M,Bc]=lag_plant(p,opt.tau); ST(i,3)=all(real(eig(M-Bc*Kc))<0);
    ST(i,4)=lqg_stable(p,Kc,L,Cm,opt.tau);
    if mod(i,10)==0, fprintf('%d/%d\n',i,N); end
end
fprintf('\n%-20s %10s %10s %10s | %10s %10s | %s\n','controller','aw median','aw p95','aw max','tyre med','tyre p95','stable');
for j=1:4
    fprintf('%-20s %10.3f %10.3f %10.3f | %10.4f %10.4f | %d/%d\n',names{j},median(AW(:,j)),prctile_(AW(:,j),95),max(AW(:,j)),median(TY(:,j)),prctile_(TY(:,j),95),sum(ST(:,j)),N);
end
fig=figure('Position',[100 100 1000 420]);
subplot(1,2,1); boxplot_(AW,names); ylabel('ISO 2631 weighted accel (m/s^2)'); title('Ride harshness across 100 car variants');
subplot(1,2,2); boxplot_(TY,names); ylabel('RMS tyre-load ratio'); title('Road-holding loss across 100 car variants');
try exportgraphics(fig,'media/phase5_robustness.png','Resolution',200); catch, saveas(fig,'media/phase5_robustness.png'); end
save('phase5_results.mat','AW','TY','ST','ms','mu','kt','cs');

function ok = lqg_stable(p,K,L,Cm,tau)
[A,Bu,~,C,Du]=build_model(p); [M,Bc]=lag_plant(p,tau); [~,Bf]=wk_ss();
G=zeros(15); G(1:11,1:11)=M; G(1:11,5:11)=G(1:11,5:11)-Bc*K(5:11); G(1:11,12:15)=-Bc*K(1:4);
G(6:11,1:4)=0; G(6:11,12:15)=Bf*C(1,:); G(12:15,1:4)=L*Cm; G(12:15,5)=Bu; G(12:15,12:15)=A-L*Cm;
ok=all(real(eig(G))<0);
end
function q = prctile_(x,p), x=sort(x); q=interp1(linspace(0,100,numel(x)),x,p); end
function boxplot_(X,names)   % toolbox-free box plot
hold on; grid on;
for j=1:size(X,2)
    x=sort(X(:,j)); q=interp1(linspace(0,100,numel(x)),x,[5 25 50 75 95]);
    plot([j j],q([1 2]),'k'); plot([j j],q([4 5]),'k');
    rectangle('Position',[j-0.25 q(2) 0.5 q(4)-q(2)],'FaceColor',[0.612 0.753 0.933],'EdgeColor',[0.165 0.471 0.839]);
    plot([j-0.25 j+0.25],[q(3) q(3)],'Color',[0.071 0.216 0.42],'LineWidth',2);
end
set(gca,'XTick',1:size(X,2),'XTickLabel',names); xtickangle(20); xlim([0.4 size(X,2)+0.6]);
end
