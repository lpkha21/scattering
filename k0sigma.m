function s = k0sigma(T, fin, fi0)
% Normalised differential scattering cross-section k0*sigma_d (paper Eqs. 28-30).
%   T   : aggregate T-matrix T_{L+1}, size (2N+1) x (2N+1)
%   fin : incidence angle (direction the plane wave travels), radians
%   fi0 : observation angle(s), radians (scalar or array)
% Backscattering: fi0 = fin + pi.  Forward scattering: fi0 = fin.
% k0*sigma_d = 2*pi*k0*|f|^2 = 4*|sum_n (-i)^n e^{i n fi0} a_n|^2,  a = T*A_in.
N = (size(T,1) - 1)/2;
n = (-N:N).';
a = T * ((1j).^n .* exp(-1j*n*fin));            % scattered coefficients
S = ((-1j).^n .* a).' * exp(1j*n*fi0(:).');      % one value per observation angle
s = reshape(4*abs(S).^2, size(fi0));
end
