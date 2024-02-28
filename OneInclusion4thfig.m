clc; clear;
ka1 = 0.5; ka2 = 0.5; ka3 = 2.5; n1 = 4; n2 = 4; n3 = 2; kd1 = 1.5; kd2 = 1; th1 = 0; th2 = pi;

N = 12;
Tm = double(T_matrix_L_mag(0,N,[ka3, ka3], [ n3, n3], [kd1], [th1]));


fin = linspace(0, 2*pi, 40);

sigmab = zeros(size(fin));
sigmaf = zeros(size(fin));

for v = 1:length(fin)
    fL1 = 0;
    fL2 = 0;
    f1t = 0;
    f2t = 0;
    
    
    
    fi = linspace(0, 2*pi, 10000);
    f1opt = 0;
    f2opt = 0;
    
    
    for n = -N:N

        for np = -N:N
            
            f1opt = f1opt + ((-1j)^n) * ((1j)^np) * exp(1j*(n*fi - fin(v)*np)) * Tm(n + N + 1, np + N + 1);
            f2opt = f2opt + ((-1j)^n) * ((1j)^np) * exp(1j*(n*fin(v) - fin(v)*np)) * Tm(n + N + 1, np + N + 1);            
            
            f1t = f1t + (-1j)^n * (1j)^np * exp(1j*((0)*n - np*fin(v))) * Tm(n + N + 1, np + N + 1);
            f2t = f2t + (-1j)^n * (1j)^np * exp(1j*((fin(v) + pi)*n - np*pi)) * Tm(n + N + 1, np + N + 1);
            
            
            fi0 = fin(v) + pi;
            fL1 = fL1 + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fi0 - np*fin(v))) .* Tm(n + N + 1, np + N + 1);
            fL2 = fL2 + (-1j).^n .* (1j) .^ np .* exp(1j *(n*fin(v) - np*fin(v))) .* Tm(n + N + 1, np + N + 1);
        end
    end
    f1opt = f1opt*(1-1j)/sqrt(pi);
    f2opt = f2opt*(1-1j)/sqrt(pi);
    
    sigs = trapz(fi, (abs(f1opt)).^2);
    sige = -2*sqrt(pi)*real((1+1j)*f2opt);
    f1t = (1-1j)*f1t;
    f2t = (1-1j)*f2t;
    
    sigd1 = 2*(abs(f1t))^2;
    sigd2 = 2*(abs(f2t))^2;
    
    abs(sigd2 - sigd1)
    abs(sigs - sige)
    
    
    fL1 = (1-1j)/(sqrt(pi)) * fL1;
    fL2 = (1-1j)/(sqrt(pi)) * fL2;
    sigmab(v) = 2*pi*((abs(fL1)).^2); % function to plot TL    
    sigmaf(v) = 2*pi*((abs(fL2)).^2);
end

plot(fin, sigmab)
hold on
plot(fin, sigmaf)
xlim([0,2*pi])
ylim([0, 120])
