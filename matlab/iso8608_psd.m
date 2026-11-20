function S = iso8608_psd(f, v, cls)
%ISO8608_PSD  One-sided road-height PSD [m^2/Hz] at temporal frequency f [Hz], speed v [m/s].
%   ISO 8608:  Gd(n) = Gd(n0) (n/n0)^-2,  n0 = 0.1 cycle/m,  n = f/v.
%   Class A 16e-6 | B 64e-6 | C 256e-6 | D 1024e-6  m^3  (smooth ... rough).
Gd0 = struct('A',16e-6, 'B',64e-6, 'C',256e-6, 'D',1024e-6);
n0  = 0.1;
S   = Gd0.(upper(cls)) * n0^2 * v ./ f.^2;
end
