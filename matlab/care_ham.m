function X = care_ham(A,B,Q,R,N)
%CARE_HAM  Continuous algebraic Riccati equation (with cross term N) via the
%  Hamiltonian eigenvector method. Core MATLAB, no Control System Toolbox.
if nargin < 5, N = zeros(size(B)); end
Ri = inv(R); At = A - B*Ri*N'; Qt = Q - N*Ri*N';
H = [At, -B*Ri*B'; -Qt, -At'];
[V,E] = eig(H); ev = diag(E); Vs = V(:, real(ev) < 0);
n = size(A,1); X = real(Vs(n+1:end,:)/Vs(1:n,:)); X = (X+X')/2;
end
