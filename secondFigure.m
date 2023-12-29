clc;
N = 12; M = 12;
fin = 0; ka1 = 0.25; ka2 = 0.25; n1 = 2; n2 = 2; kd2 = 1.5; t1 = 0; t2 = pi; 
ka3 = 3;
ep0 = 8.854187 * 1e-12;
mi0 = 1.256637 * 1e-6;
ep3 = linspace(ep0,4*ep0); % MAIN independent argument
n3 = sqrt(ep3/ep0);

% --------- (1) and (3) 
kd1 = 1.5;
fi0 = fin + pi;
TL1 = T_matrix_L(M,0,0); % TODO T(L+1) aggregate matrix

f = 0;

for n = -N:N
    for np = -N:N
        f = f + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0 - np*fin)) .* TL1(n + N + 1, np + N + 1);
    end
end
f = (1-1j)/(sqrt(pi)) * f;
sigmab = 2*pi*abs(f); % function to plot
plot(ep3,sigmab,'x');