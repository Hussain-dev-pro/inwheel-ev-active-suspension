function s = td_metrics(Y, skip, dt, umax)
%TD_METRICS  Metrics from a time-domain run (first 'skip' samples discarded).
%  ISO 2631 weighting applied via FFT (Parseval) over 0.5-80 Hz.
a=Y(skip+1:end,1); n=numel(a); F=(0:floor(n/2))'/(n*dt); S=fft(a); S=S(1:numel(F));
W=zeros(size(F)); b=F>=0.5 & F<=80; W(b)=wk_weight(F(b)).';
p2=abs(S.*W).^2; p2(2:end-1)=2*p2(2:end-1);
s.aw   = sqrt(sum(p2)/n^2);
s.trav = std(Y(skip+1:end,2));
s.tyre = std(Y(skip+1:end,3));
s.umax = max(abs(Y(skip+1:end,4)));
s.urms = std(Y(skip+1:end,4));
s.sat  = 100*mean(abs(Y(skip+1:end,4)) >= umax-1e-6);
s.P    = mean(max(Y(skip+1:end,5),0));
end
