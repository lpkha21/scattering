clc; clear; close all;

ka = 2.5; n3 = 2;
fi0 = 0;

Tm = double(T_matrix_L_mag(0 ,12,[2.5], [2],[1.5, 1],[0, pi]));
N = 12;
sigb = zeros(1,100);
in = 1;
for fin = linspace(0,2*pi,100)
    
    % reciprocity ----------------------
    f1 = 0;
    f2 = 0;
    fPlot = 0;
    for n = -N:N
       for np = -N:N
%           f1 = f1 + ((-1j)^n) * ((1j)^np) * exp(1j*(n*fi0 - fin*np)) * Tm(n + N + 1, np + N + 1);
%           f2 = f2 + ((-1j)^n) * ((1j)^np) * exp(1j*(n*(fin + pi) - (fi0 + pi)*np)) * Tm(n + N + 1, np + N + 1);
          fPlot = fPlot + ((-1j)^n) * ((1j)^np) * exp(1j*(n*(fin + pi) - (fin)*np)) * Tm(n + N + 1, np + N + 1);
       end
    end
    fPlot = (1-1j)/sqrt(pi) * fPlot;
    sigb(in) = 2*pi*((abs(fPlot))^2);
    in = in + 1;
%     f1 = (abs((1-1j)/sqrt(pi)*f1))^2;
%     f2 = (abs((1-1j)/sqrt(pi)*f2))^2;
%     
%     formatted_number = sprintf('%.5e\n', f2 - f1);
%     result = ['reciprocity ', formatted_number];
%     disp(result);
    
%     fi = 0:0.000001:2*pi;
%     f1 = 0;
%     f2 = 0;
%     for n = -N:N
%        for np = -N:N
%           f1 = f1 + ((-1j)^n) * ((1j)^np) * exp(1j*(n*fi - fin*np)) * Tm(n + N + 1, np + N + 1);
%           f2 = f2 + ((-1j)^n) * ((1j)^np) * exp(1j*(n*fin - fin*np)) * Tm(n + N + 1, np + N + 1);
%       end
%     end
%     f1 = f1*(1-1j)/sqrt(pi);
%     f2 = f2*(1-1j);
%     
%     sigs = trapz(fi, ((abs(f1)).^2));
%     sige = -2*real((1+1j)*f2);
%     
%     formatted_number = sprintf('%.5e\n', sige - sigs);
%     result = ['optical ', formatted_number];
%     disp(result);
    
end



plot(linspace(0,2*pi,100), sigb)