function T = T_matrix_L_mag(L,M,k0ai,ni,k0di,Thi)
% TM aggregate T-matrix T_{L+1} of the host cylinder with L inclusions
% (paper Eqs. 7, 24, 26-27).
%   k0ai = [k0*a_1 ... k0*a_L, k0*a_host]
%   ni   = [n_1 ... n_L, n_host]       (Inf = perfect conductor)
%   k0di = [k0*d_1 ... k0*d_L],  Thi = [theta_1 ... theta_L] in radians

% FIX: waves are translated between cylinders inside the host, where the
% wavenumber is k_host = n_host*k0, so alpha/beta need k_host*d, not k0*d.
k0di = ni(L+1) * k0di;

% ------ T_L: aggregate T-matrix of the L inclusions seen from the host centre, Eq. (7)
Tm = zeros(2*M + 1);
for k = 1:L
    Tm = Tm + beta(M,k0di,Thi,0,k) * TmLMag(L,L,M,k0ai,k0di,Thi,k,ni);
end
% ------ outer boundary of the host
f  = besJ(M,L+1,L+1,k0ai,ni)  + besH(M,L+1,L+1,k0ai,ni)  * Tm;    % Eq. (26)
df = dbesJ(M,L+1,L+1,k0ai,ni) + dbesH(M,L+1,L+1,k0ai,ni) * Tm;    % Eq. (27)
G = df / f;

first  = dbesH(M,0,L+1,k0ai,ni) - ni(L+1)*G*besH(M,0,L+1,k0ai,ni);
second = dbesJ(M,0,L+1,k0ai,ni) - ni(L+1)*G*besJ(M,0,L+1,k0ai,ni);

T = -(first \ second);                                            % Eq. (24)
end
