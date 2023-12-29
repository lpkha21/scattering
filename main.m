clc
% TODO - 
% T matrix

M = 12;
N = 12;
eps0 = 8.8541878128 * 10e-12;
mi0 = 1.25663706212 * 10e-6;

ephsi = @(ni) ni * ni * eps0;
niu = @(eps, miu) sqrt(eps / miu);


% permiability - everywhere miu0
% L = 2 
% k0a1 = k0a2 = 0.25;  n1=n2=2; k0d2 = 1.5; th1 = 0; th2 = 180; 
% k0a3 = 3; fin = 0;
% TM wave 3) k0d1 = 1.5 4) k0d1 = 2.25
n1 = 2; n2 = 2;
n3 = linspace(n1 - 1, n1, 1000); % deltaN = n1-n3

% k0-s ar vuwer radgan k0*sigmab unda avago sadac k0 ikveceba.
TT = tMatMag(2, [ephsi(n1), ephsi(n2)], [0.25, 0.25], [1.5, 1.5], 2, [0, 180], n3);

% Use linear indices to access elements in TT
t = TT;

[n, np] = ndgrid(-N:1:N, -N:1:N);

% Reshape the matrix to match the size of n and np


f = @(fi0, fi1) (1 - 1j) / sqrt(pi) .* sum(sum((-1j).^n .* (-1j).^np .* exp(1j .* (n .* fi0 - np .* fi1)) .* t));

sigmab = @(Fin) 2 .* pi .* (abs(f(Fin + pi, Fin))) .^ 2;

plot(n3, sigmab(0), 'x');