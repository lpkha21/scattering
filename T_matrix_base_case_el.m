function t = T_matrix_base_case_el(L,M,k0ai,ni,si)
% TE (H_z) T-matrix of inclusion si alone in the host medium, paper Eq. (13).
% Same as the TM case with the two refractive-index prefactors swapped.
% ni(si) = Inf means a perfect conductor: T = -J'/H' (exact limit).
JL  = besJ(M,L+1,si,k0ai,ni);   HL  = besH(M,L+1,si,k0ai,ni);    % host side, k_host*a_si
dJL = dbesJ(M,L+1,si,k0ai,ni);  dHL = dbesH(M,L+1,si,k0ai,ni);
if isinf(ni(si))
    t = -(dHL \ dJL);
    return;
end
R = dbesJ(M,si,si,k0ai,ni) / besJ(M,si,si,k0ai,ni);             % J'(k_i a) / J(k_i a)
t1 = ni(si)*dHL - ni(L+1)*R*HL;
t2 = ni(si)*dJL - ni(L+1)*R*JL;
t = -(t1 \ t2);
end
