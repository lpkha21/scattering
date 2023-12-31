clc;
N = 12; M = 12;
fin = 0; ka1 = 0.25; ka2 = 0.25; n1 = 2; n2 = 2; kd2 = 1.5; t1 = 0; t2 = pi; 
ka3 = 3;
ep0 = 8.854187 * 1e-12;
mi0 = 1.256637 * 1e-6;
ep3 = linspace(ep0,4*ep0,100); % MAIN independent argument
n3 = sqrt(ep3/ep0); % MAIN independent argument
sigmab = zeros(1,length(n3));

% --------- (1) and (3) graphs

kd1 = 1.5;
fi0 = fin + pi;
L = 2;
parfor v = 1:length(n3)
    
    TL1 = T_matrix_L_mag(L, M, [ka1, ka2, ka3], [n1, n2, n3(v)], [kd1, kd2], [0 pi]); 
    f = 0;
    for n = -N:N
        for np = -N:N
            f = f + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0 - np*fin)) .* TL1(n + N + 1, np + N + 1);
        end
    end
    f = (1-1j)/(sqrt(pi)) * f;
    sigmab(v) = 2*pi*((abs(f)).^2); % function to plot
end
plot(n1-n3,sigmab);
% hold on

% kd1 = 2.25;
% fi0 = fin + pi;
% for v = 1:length(n3)
%     
%     TL1 = T_matrix_L_mag(2, M, [ka1, ka2, ka3], [n1, n2, n3(v)], [kd1, kd2], [0 pi]); 
%     f = 0;
% 
%     for n = -N:N
%         for np = -N:N
%             f = f + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0 - np*fin)) .* TL1(n + N + 1, np + N + 1);
%         end
%     end
%     f = (1-1j)/(sqrt(pi)) * f;
%     sigmab(v) = 2*pi*((abs(f)).^2); % function to plot
% end
% 
% plot(n1-n3,sigmab);