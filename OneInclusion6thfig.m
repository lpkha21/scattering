clear;clc;
ka1 = 0.5; ka2 = 2.5; n1 = 1e308; n2 = 2.0; kd1 = 1.5; th1 = 0;


fin = linspace(0, 2*pi, 50);
M = 12;
Tm = double(T_matrix_L_mag(1, M, [ka1 ka2], [n1 n2], kd1, th1));

f0 = linspace(0,2*pi, 1000);


for v = 1:length(fin)
    fp1 = 0;
    fp2 = 0;
    
    for n = -M:M
        for np = -M:M
            fp1 = fp1 + (-1j)^n * (1j)^np .*  exp(1j*(n * f0 - np * fin(v)))* Tm(n + M + 1, np + M + 1);
            fp2 = fp2 + (-1j)^n * (1j)^np .*  exp(1j*(n * fin(v) - np * fin(v)))* Tm(n + M + 1, np + M + 1);
            
        end
    end
    
    sige = -2 * real((1+1j)*(1-1j)*fp2);
    sigs = trapz(f0, abs((1-1j)/sqrt(pi)*fp1).^2);
    
    abs(sige - sigs)
end





% hold on
% plot(fin, sigb,'r')
% plot(fin, sigf,'b')
% xlim([0, 2*pi]);
% ylim([0, 140]);
% grid on


