function t = T_matrix_base_case_mag(L,M,k0ai,ni,si)
ep0 = 8.854 * 1e-12;
mi0 = 1.256637 * 1e-6;
niu0 = sqrt(ep0/mi0);
niuL1 = ni(L+1) * niu0;
niusi = ni(si) * niu0;

t1 = niuL1* dbesH(M,L+1,si,k0ai,ni) - niusi* (dbesJ(M,si,si,k0ai,ni) / besJ(M,si,si,k0ai,ni)) *besH(M,L + 1,si,k0ai,ni);
t2 = niuL1* dbesJ(M,L+1,si,k0ai,ni) - niusi* (dbesJ(M,si,si,k0ai,ni) / besJ(M,si,si,k0ai,ni)) *besJ(M,L + 1,si,k0ai,ni);

t = - (t1\t2);
end