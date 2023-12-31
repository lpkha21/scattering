clc;
N = 12; M = 12;
ka1 = 2.5; n1 = 2; 

ep0 = 8.854187 * 1e-12;
mi0 = 1.256637 * 1e-6;

fi0 = linspace(0,2*pi,35);
fin = 2*pi/3;
sigmadMag = (zeros(1,length(fi0)));
sigmadEl =  (zeros(1,length(fi0)));


L = 0;
TL1 = T_matrix_L_mag(L, M, [ka1], [n1], [], []); 
for v = 1:length(fi0)
    
    fL = 0;
    for n = -N:N
        for np = -N:N
            fL = fL + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0(v) - np*fin)) .* TL1(n + N + 1, np + N + 1);
        end
    end
    fL = (1-1j)/(sqrt(pi)) * fL;
    sigmadMag(v) = 2*pi*((abs(fL)).^2); % function to plot TL
    
    
    TE1 = T_matrix_L_El(L, M, [ka1], [n1], [], []); 
    fE = 0;
    for n = -N:N
        for np = -N:N
            fE = fE + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0(v) - np*fin)) .* TE1(n + N + 1, np + N + 1);
        end
    end
    fE = (1-1j)/(sqrt(pi)) * fE;
    sigmadEl(v) = 2*pi*((abs(fE)).^2); % function to plot TE
    
end
plot(-fi0+2*pi,sigmadMag);
hold on
plot(-fi0+2*pi,sigmadEl,'--');