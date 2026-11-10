function t = T_matrix_base_case_mag(L,M,k0ai,ni,si)
% TM (E_z) T-matrix of inclusion si alone in the host medium, paper Eq. (12).
% With mu = mu0 everywhere eta_i/eta_0 = n_i, so refractive indices are used directly.
% ni(si) = Inf means a perfect conductor: T = -J/H (exact limit).
JL  = besJ(M,L+1,si,k0ai,ni);   HL  = besH(M,L+1,si,k0ai,ni);    % host side, k_host*a_si
dJL = dbesJ(M,L+1,si,k0ai,ni);  dHL = dbesH(M,L+1,si,k0ai,ni);
if isinf(ni(si))
    t = -(HL \ JL);
    return;
end
R = dbesJ(M,si,si,k0ai,ni) / besJ(M,si,si,k0ai,ni);             % J'(k_i a) / J(k_i a)
t1 = ni(L+1)*dHL - ni(si)*R*HL;
t2 = ni(L+1)*dJL - ni(si)*R*JL;
t = -(t1 \ t2);
end
