% Fig. 2: backscattering cross-section k0*sigma_b of a host cylinder with two
% identical inclusions vs. the index difference dn = n1 - n3 (n3 = host index).
clc; clear;
M = 12;
ka1 = 0.25; ka2 = 0.25; n1 = 2; n2 = 2; kd2 = 1.5; t1 = 0; t2 = pi;
ka3 = 3;
fin = 0; fi0 = fin + pi;            % incidence along +x, observe straight back
dn = linspace(0, 1, 101);           % FIX: uniform grid in dn (was 35 points uniform in epsilon)
n3 = n1 - dn;
L = 2;

kd1s = [1.5, 2.25];                 % symmetric (1),(3) and asymmetric (2),(4)
sigTM = zeros(2, length(dn));
sigTE = zeros(2, length(dn));
for c = 1:2
    kd1 = kd1s(c);
    tm = zeros(1, length(dn));
    te = zeros(1, length(dn));
    parfor v = 1:length(dn)
        TM = T_matrix_L_mag(L, M, [ka1, ka2, ka3], [n1, n2, n3(v)], [kd1, kd2], [t1, t2]);
        TE = T_matrix_L_El(L, M, [ka1, ka2, ka3], [n1, n2, n3(v)], [kd1, kd2], [t1, t2]);
        tm(v) = k0sigma(TM, fin, fi0);
        te(v) = k0sigma(TE, fin, fi0);
    end
    sigTM(c, :) = tm;
    sigTE(c, :) = te;
end

figure; hold on; grid on;
plot(dn, sigTE(1,:), 'r--', 'LineWidth', 1.5);
plot(dn, sigTE(2,:), 'r-',  'LineWidth', 1.5);
plot(dn, sigTM(1,:), 'b--', 'LineWidth', 1.5);
plot(dn, sigTM(2,:), 'b-',  'LineWidth', 1.5);
legend('(1) TE, k_0d_1 = 1.5', '(2) TE, k_0d_1 = 2.25', ...
       '(3) TM, k_0d_1 = 1.5', '(4) TM, k_0d_1 = 2.25');
xlim([0 1]); ylim([0 25]);
xlabel('\Delta n = n_1 - n_3'); ylabel('k_0\sigma_b');
title('Fig. 2: two inclusions, \phi^{in} = 0');
