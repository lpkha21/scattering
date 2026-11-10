% Fig. 8: differential scattering cross-section k0*sigma_d vs. observation angle,
% homogeneous cylinder (no inclusions), TM and TE, incidence phi_in = 270 deg (wave travelling in -y).
clc; clear;
M = 12;
ka1 = 2.5; n1 = 2;
fin = 3*pi/2;                       % FIX: 270 deg (was 2*pi/3 = 120 deg)
fi0 = linspace(0, 2*pi, 361);       % FIX: finer grid (was 35 points)
L = 0;

TM = T_matrix_L_mag(L, M, ka1, n1, [], []);
TE = T_matrix_L_El(L, M, ka1, n1, [], []);
sigTM = k0sigma(TM, fin, fi0);
sigTE = k0sigma(TE, fin, fi0);

figure; hold on; grid on;
plot(rad2deg(fi0), sigTM, 'b-',  'LineWidth', 1.5);   % FIX: plot against phi0 directly (no flip)
plot(rad2deg(fi0), sigTE, 'r--', 'LineWidth', 1.5);
legend('TM', 'TE');
xlim([0 360]); xticks(0:90:360); ylim([0 140]);
xlabel('\phi_0 (deg)'); ylabel('k_0\sigma_d');
title('Fig. 8: homogeneous cylinder (no inclusions)');
