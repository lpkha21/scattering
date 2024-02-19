function t = T_matrix_base_case_mag(L,M,k0ai,ni,si)
ep0 = 8.8541878128 * 1e-12;
mi0 = 1.25663706212 * 1e-6;
niu0 = vpa(sqrt(ep0/mi0),100);
niuL1 = vpa(ni(L+1) * niu0,100);
niusi = vpa(ni(si) * niu0,100);
  
t1 = niuL1* vpa(dbesH(M,L+1,si,k0ai,ni),50) - niusi* double(dbesJ(M,si,si,k0ai,ni)) / vpa((besJ(M,si,si,k0ai,ni)),50) *vpa(besH(M,L + 1,si,k0ai,ni),50);
t2 = niuL1* vpa(dbesJ(M,L+1,si,k0ai,ni),50) - niusi* double(dbesJ(M,si,si,k0ai,ni)) / vpa((besJ(M,si,si,k0ai,ni)),50) *vpa(besJ(M,L + 1,si,k0ai,ni),50);

t = - ((vpa(t1,100)\vpa(t2,100)));
end