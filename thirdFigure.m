clc;
N = 12; M = 12;
ka1 = 0.25; ka2 = 0.25; ka3 = 0.25; n1 = 2; n2 = 2; n3 = 2; kd2 = 1.5; kd3 = 1.5; t1 = pi/2; t2 = deg2rad(210); t3 = deg2rad(330); 
ka4 = 3; n4 = 1.5;
ep0 = 8.854187 * 1e-12;
mi0 = 1.256637 * 1e-6;
kd1 = linspace(0,3,25);
sigmabMag1 = zeros(1,length(kd1));
sigmabMag2 = zeros(1,length(kd1));
sigmabMag3 = zeros(1,length(kd1));
L = 3;


for v = 1:length(kd1)

    TL1 = T_matrix_L_mag(L, M, [ka1, ka2, ka3, ka4], [n1, n2, n3, n4], [kd1(v), kd2, kd3], [t1, t2, t3]); 
    fL1 = 0;
    fL2 = 0;
    fL3 = 0;
    for n = -N:N
        for np = -N:N
            fin = 0;
            fi0 = fin + pi;
            fL1 = fL1 + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0 - np*fin)) .* TL1(n + N + 1, np + N + 1);
            
            fin = pi/2;
            fi0 = fin + pi;
            fL2 = fL2 + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0 - np*fin)) .* TL1(n + N + 1, np + N + 1);
            
            fin = 3*pi/2;
            fi0 = fin + pi;
            fL3 = fL3 + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0 - np*fin)) .* TL1(n + N + 1, np + N + 1);
        end
    end
    fL1 = (1-1j)/(sqrt(pi)) * fL1;
    sigmabMag1(v) = 2*pi*((abs(fL1)).^2); % function to plot TL    

    fL2 = (1-1j)/(sqrt(pi)) * fL2;
    sigmabMag2(v) = 2*pi*((abs(fL2)).^2); % function to plot TL   
   
    fL3 = (1-1j)/(sqrt(pi)) * fL3;
    sigmabMag3(v) = 2*pi*((abs(fL3)).^2); % function to plot TL   
end
plot(kd1,sigmabMag1);
hold on;
plot(kd1,sigmabMag2,'.');
plot(kd1,sigmabMag3,'--');
ylim([6 12]);
