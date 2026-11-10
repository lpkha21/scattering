% Fig. 6: backscattering k0*sigma_b and forward scattering k0*sigma_f vs. incidence
% angle, TM, perfectly conducting inclusions.
clc; clear;
M = 12;
ka1 = 0.5; ka2 = 0.5; kd1 = 1.5; kd2 = 1; t1 = 0; t2 = pi;
ka3 = 2.5; n3 = 2;
nin = Inf;                     % index of the inclusions (Inf = perfect conductor)
Tfun = @T_matrix_L_mag;      % TM
fin = linspace(0, 2*pi, 361);         % FIX: finer grid (was 50 points)

Ts = cell(1, 3);
Ts{1} = Tfun(1, M, [ka1, ka3], [nin, n3], kd1, t1);                          % (1) single inclusion
Ts{2} = Tfun(2, M, [ka1, ka2, ka3], [nin, nin, n3], [kd1, kd2], [t1, t2]);   % (2) two inclusions
Ts{3} = Tfun(0, M, ka3, n3, [], []);                                          % (3) no inclusions (FIX: was missing)

sb = zeros(3, length(fin));
sf = zeros(3, length(fin));
for c = 1:3
    for v = 1:length(fin)
        sb(c, v) = k0sigma(Ts{c}, fin(v), fin(v) + pi);   % backscattering: look back at the source
        sf(c, v) = k0sigma(Ts{c}, fin(v), fin(v));        % FIX: forward = sigma_d at phi0 = phi_in
    end                                                    % (was sigma_t - sigma_b, which is wrong)
end

cols  = [0.85 0.33 0.10; 0 0.45 0.74; 0 0 0];
names = {'(1) one inclusion', '(2) two inclusions', '(3) no inclusions'};
deg = rad2deg(fin);
figure; hold on; grid on;
for c = 1:3
    plot(deg, sf(c,:), '-',  'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', [names{c} ', forward']);
    plot(deg, sb(c,:), '--', 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', [names{c} ', back']);
end
legend('show', 'Location', 'bestoutside');
xlim([0 360]); xticks(0:90:360); ylim([0 140]);
xlabel('\phi^{in} (deg)'); ylabel('k_0\sigma_b, k_0\sigma_f');
title('Fig. 6: TM, perfectly conducting inclusions');
