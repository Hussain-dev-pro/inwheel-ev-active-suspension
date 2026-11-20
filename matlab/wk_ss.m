function [Af,Bf,Cf,Df] = wk_ss()
%WK_SS  State-space ISO 2631 Wk weighting (high-pass x transition x step; 6 states).
%  Used INSIDE the controller design so LQR minimises the comfort people feel.
w1=2*pi*0.4; w3=2*pi*12.5; w4=2*pi*12.5; w5=2*pi*2.37; w6=2*pi*3.35;
Q1=1/sqrt(2); Q4=0.63; Q5=0.91; Q6=0.91;
S1 = biquad(1,0,0,w1/Q1,w1^2);
S2 = biquad(0,w4^2/w3,w4^2,w4/Q4,w4^2);
S3 = biquad(1,w5/Q5,w5^2,w6/Q6,w6^2);
[Af,Bf,Cf,Df] = series2(series2(S1,S2),S3);
end
function S = biquad(b2,b1,b0,a1,a0)
S.A=[0 1;-a0 -a1]; S.B=[0;1]; S.C=[b0-b2*a0, b1-b2*a1]; S.D=b2;
end
function varargout = series2(S1,S2)
A=[S1.A zeros(size(S1.A,1),size(S2.A,1)); S2.B*S1.C S2.A]; B=[S1.B; S2.B*S1.D];
C=[S2.D*S1.C S2.C]; D=S2.D*S1.D;
if nargout<=1, varargout{1}=struct('A',A,'B',B,'C',C,'D',D); else, varargout={A,B,C,D}; end
end
