function T = T_matrix_L_mag(L,M,k0ai,ni,k0di,Thi)
ep0 = 8.8541878128 * 1e-12;
mi0 = 1.25663706212 * 1e-6;
niu0 = vpa(sqrt(ep0/mi0),100);
niuL = vpa(ni(L+1)*niu0,100);
% ------ calculating TL aggregate T-matrix for the L cylinders
Tm = 0;
for k = 1:L
    Tm = Tm + vpa(beta(M,k0di,Thi,0,k),100)*vpa(TmLMag(L,L, M, k0ai, k0di, Thi, k, ni),100);
end
% ------
f = vpa(besJ(M,L+1,L+1,k0ai,ni),100) + vpa(besH(M,L+1,L+1,k0ai,ni),100) * Tm;
df = vpa(dbesJ(M,L+1,L+1,k0ai,ni),100) + vpa(dbesH(M,L+1,L+1,k0ai,ni),100) * Tm;

first = niu0*dbesH(M,0,L+1,k0ai,ni) - niuL*((vpa(df,100) / vpa(f,100)) * besH(M,0,L+1,k0ai,ni));
second = niu0*dbesJ(M,0,L+1,k0ai,ni) - niuL*((vpa(df,100) / vpa(f,100)) * besJ(M,0,L+1,k0ai,ni));

T = -(vpa(first,100) \ vpa(second,100));
end