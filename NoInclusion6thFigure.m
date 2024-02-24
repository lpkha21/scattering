clc; clear;
ka1 = 0.5; ka2 = 2.5; n1 = 2e40; n2 = 2; th1 = 0; kd1 = 1.5;


fin = linspace(0, 2*pi, 40);

sigb = zeros(size(fin));
sigf = zeros(size(fin));
sigt = zeros(size(fin));
Tm = double(T_matrix_L_mag(1, 12, [ka1, ka2], [n1, n2], kd1, th1));



for v = 1:length(fin)
    fb = 0;
    ft = 0;
    
    f1t = 0;
    f2t = 0;
    
    
    fi = linspace(0, 2*pi, 10000);
    f1opt = 0;
    f2opt = 0;
    for n = -12:12
        finside = 0;
        for np = -12:12
            
            f1opt = f1opt + ((-1j)^n) * ((1j)^np) * exp(1j*(n*fi - fin(v)*np)) * Tm(n + 12 + 1, np + 12 + 1);
            f2opt = f2opt + ((-1j)^n) * ((1j)^np) * exp(1j*(n*fin(v) - fin(v)*np)) * Tm(n + 12 + 1, np + 12 + 1);
            
            f1t = f1t + (-1j)^n * (1j)^np * exp(1j*((0)*n - np*fin(v))) * Tm(n + 13, np + 13);
            f2t = f2t + (-1j)^n * (1j)^np * exp(1j*((fin(v) + pi)*n - np*pi)) * Tm(n + 13, np + 13);
            
            fb = fb + (-1j)^n * (1j)^np * exp(1j*((fin(v)+2*pi)*n - np*fin(v))) * Tm(n + 13, np + 13);
            finside = finside + (1j)^np * exp(-1j * np * fin(v)) * Tm(n + 13, np + 13);
        end
        
        ft = ft + ((abs(finside))^2);
    end
    % testing theorems -------
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
    % ----------------------
    
    fb = (1-1j)*fb;
    sigb(v) = 2 * (abs(ft))^2;
    sigt(v) = 4 * ft;
    
end

plot(fin, sigb)
hold on
% plot(fin, sigt - sigb,'--')
legend('sigb', 'sigf');
xlim([0, 2*pi])

