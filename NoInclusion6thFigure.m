%% no inclusion - WORKING!!
clear;clc;
n1 = 2.0; ka1 = 2.5;

M = 12;

Tm = double(T_matrix_L_mag(0,M,ka1, n1, 0,0));

fin = linspace(0, 2*pi, 30);

fb = 0;
ff = 0;
for n = -M:M
   for np = -M:M
      fb = fb + (-1j)^n .* (1j) ^ np .* exp(1j*(n * (fin + pi) - np * fin)) * Tm(n + M + 1, np + M + 1);
      ff = ff + (-1j)^n .* (1j) ^ np .* exp(1j*(n * fin - np * fin)) * Tm(n + M + 1, np + M + 1);
   end    
end

sigb = 2*(abs((1-1j)*fb).^2);
sigf = 2*(abs((1-1j)*ff).^2);
plot(fin, sigb)
hold on
plot(fin, sigf)
xlim([0,2*pi])
ylim([0, 140])

