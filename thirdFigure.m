% Fig. 3: TM backscattering cross-section of a host cylinder with three
% inclusions vs. the radial position k0*d1 of cylinder 1.
clc; clear;
M = 12;
ka1 = 0.25; ka2 = 0.25; ka3 = 0.25; n1 = 2; n2 = 2; n3 = 2;
kd2 = 1.5; kd3 = 1.5; t1 = pi/2; t2 = deg2rad(210); t3 = deg2rad(330);
ka4 = 3; n4 = 1.5;
kd1 = linspace(0, 2.75, 56);        % FIX: finer grid (was 25 points)
fins = [0, pi/2, 3*pi/2];           % (1) 0 deg, (2) 90 deg, (3) 270 deg
L = 3;

sigb = zeros(3, length(kd1));
parfor v = 1:length(kd1)
    T = T_matrix_L_mag(L, M, [ka1, ka2, ka3, ka4], [n1, n2, n3, n4], [kd1(v), kd2, kd3], [t1, t2, t3]);
    s = zeros(3, 1);
    for c = 1:3
        s(c) = k0sigma(T, fins(c), fins(c) + pi);
    end
    sigb(:, v) = s;
end

% accuracy check (paper Eq. 32): reciprocity for one configuration
T = T_matrix_L_mag(L, M, [ka1, ka2, ka3, ka4], [n1, n2, n3, n4], [2.25, kd2, kd3], [t1, t2, t3]);
fprintf('reciprocity error: %.2e\n', abs(k0sigma(T, 0.3, 2.1) - k0sigma(T, 2.1 + pi, 0.3 + pi)));

figure; hold on; grid on;
plot(kd1, sigb(1,:), 'k-',  'LineWidth', 1.5);
plot(kd1, sigb(2,:), 'b:',  'LineWidth', 1.5);
plot(kd1, sigb(3,:), 'r--', 'LineWidth', 1.5);
legend('(1) \phi^{in} = 0', '(2) \phi^{in} = 90', '(3) \phi^{in} = 270');
xlim([0 2.75]); ylim([6 12]);
xlabel('k_0d_1'); ylabel('k_0\sigma_b');
title('Fig. 3: three inclusions, TM');
