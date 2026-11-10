% Fig. 12: differential scattering cross-section k0*sigma_d vs. observation angle,
% air hole, dielectric (n = 4) and perfect conductor, TM and TE, incidence phi_in = 270 deg (wave travelling in -y).
clc; clear;
M = 12;
ka1 = 0.3; ka2 = 0.3; ka3 = 0.3;
kd1 = 0.35; kd2 = 1; kd3 = 1.65; t1 = pi/2; t2 = deg2rad(210); t3 = deg2rad(330);
ka4 = 2.5; n4 = 2;
n1 = 1; n2 = 4; n3 = Inf;         % FIX: exact perfect conductor (was 1e10)
fin = 3*pi/2;                       % FIX: 270 deg (was 2*pi/3 = 120 deg)
fi0 = linspace(0, 2*pi, 361);       % FIX: finer grid (was 35 points)
L = 3;

TM = T_matrix_L_mag(L, M, [ka1, ka2, ka3, ka4], [n1, n2, n3, n4], [kd1, kd2, kd3], [t1, t2, t3]);
TE = T_matrix_L_El(L, M, [ka1, ka2, ka3, ka4], [n1, n2, n3, n4], [kd1, kd2, kd3], [t1, t2, t3]);
sigTM = k0sigma(TM, fin, fi0);
sigTE = k0sigma(TE, fin, fi0);

figure; hold on; grid on;
plot(rad2deg(fi0), sigTM, 'b-',  'LineWidth', 1.5);   % FIX: plot against phi0 directly (no flip)
plot(rad2deg(fi0), sigTE, 'r--', 'LineWidth', 1.5);
legend('TM', 'TE');
xlim([0 360]); xticks(0:90:360); ylim([0 140]);
xlabel('\phi_0 (deg)'); ylabel('k_0\sigma_d');
title('Fig. 12: air hole, dielectric (n = 4) and perfect conductor');
