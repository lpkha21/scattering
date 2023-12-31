function T = T_matrix_L_mag(L,M,k0ai,ni,k0di,Thi)
ep0 = 8.854187 * 1e-12;
mi0 = 1.256637 * 1e-6;
niu0 = sqrt(ep0/mi0);
niuL = ni(L+1)*niu0;
% ------ calculating TL aggregate T-matrix for the L cylinders
Tm = zeros(2*M + 1);
for k = 1:L
    q =  beta(M,k0di,Thi,0,k)*vpa(TmLMag(L,L, M, k0ai, k0di, Thi, k, ni));
    Tm = Tm + q;
end
% ------
f = besJ(M,L+1,L+1,k0ai,ni) + besH(M,L+1,L+1,k0ai,ni) * Tm;
df = dbesJ(M,L+1,L+1,k0ai,ni) + dbesH(M,L+1,L+1,k0ai,ni) * Tm;
    
first = niu0*dbesH(M,0,L+1,k0ai,ni) - niuL*((vpa(df) / vpa(f)) * besH(M,0,L+1,k0ai,ni));
second = niu0*dbesJ(M,0,L+1,k0ai,ni) - niuL*((vpa(df) / vpa(f)) * besJ(M,0,L+1,k0ai,ni));

T = -vpa(vpa(first) \ vpa(second));
end

