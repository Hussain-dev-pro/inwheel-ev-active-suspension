function K = lqr_lag(p, qa, qd, qt, r, tau)
%LQR_LAG  Frequency-shaped, lag-aware LQR gain. u_cmd = -K*[x; force; xf]
%  qa: weight on ISO-weighted body accel, qd: suspension travel, qt: tyre-load ratio, r: force
[M,Bc,~,Z] = lag_plant(p, tau);
Q = Z'*diag([qa qd qt])*Z + 1e-12*eye(11);
X = care_ham(M, Bc, Q, r, zeros(11,1));
K = (Bc'*X)/r;
end
