clear;clc;
ka1 = 0.5; ka2 = 2.5; n1 = 1e308; n2 = 2.0; kd1 = 1.5; th1 = 0;


fin = linspace(0, 2*pi, 50);
M = 12;
Tm = double(T_matrix_L_mag(1, M, [ka1 ka2], [n1 n2], kd1, th1));


fb = 0;
ff = 0;

for n = -M:M
    for np = -M:M
       fb = fb + (-1j)^n * (1j)^np .* exp(1j*(n * (fin + pi) - np * fin)) * Tm(n + M + 1, np + M + 1);
       ff = ff + (-1j)^n * (1j)^np .* exp(1j*(n * fin - np * fin)) * Tm(n + M + 1, np + M + 1);
    end
end

sigb = 2 * (abs((1-1j)*fb).^2);
sigf = 2 * (abs((1-1j)*ff).^2);

hold on
plot(fin, sigb,'r')
plot(fin, sigf,'b')
xlim([0, 2*pi]);
ylim([0, 140]);
grid on


