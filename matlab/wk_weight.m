function W = wk_weight(f)
%WK_WEIGHT  ISO 2631-1 Wk frequency weighting (vertical, z axis): magnitude at f [Hz].
%   Band-limiting (0.4 / 100 Hz) x accel-velocity transition x upward step.
%   Checked against the ISO table in selfcheck.m (error < 0.2 %).
s  = 1i*2*pi*f(:).';
f1 = 0.4;  f2 = 100;  f3 = 12.5;  f4 = 12.5;  f5 = 2.37;  f6 = 3.35;
Q1 = 1/sqrt(2);  Q2 = 1/sqrt(2);  Q4 = 0.63;  Q5 = 0.91;  Q6 = 0.91;
w1 = 2*pi*f1;  w2 = 2*pi*f2;  w3 = 2*pi*f3;  w4 = 2*pi*f4;  w5 = 2*pi*f5;  w6 = 2*pi*f6;
Hh = s.^2 ./ (s.^2 + (w1/Q1)*s + w1^2);
Hl = w2^2 ./ (s.^2 + (w2/Q2)*s + w2^2);
Ht = (w4^2/w3) * (s + w3) ./ (s.^2 + (w4/Q4)*s + w4^2);
Hs = (s.^2 + (w5/Q5)*s + w5^2) ./ (s.^2 + (w6/Q6)*s + w6^2);
W  = abs(Hh .* Hl .* Ht .* Hs);
end
