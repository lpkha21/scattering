function T = T_matrix_L_el(L,M,k0ai,ni,k0di)
ep0 = 8.854187 * 1e-12;
mi0 = 1.256637 * 1e-6;
niu0 = sqrt(ep0/mi0);niuL = ni(L+1)*niu0;
first = niuL*dbesH(M,L+1,1,k0ai,ni,1) - niu0*(dFL(L + 1).* inv(FL(L+1))).*besH(M,L+1,1,k0ai,ni,1);
second = niuL*dbesJ(M,L+1,1,k0ai,ni,1) - niu0*(dFL(L + 1).*inv(FL(L+1))).*besJ(M,L+1,1,k0ai,ni,1);
T = -inv(first) * second;


    function f = FL(L)
        f = besJ(M,L,L,k0ai,ni,0) + besH(M,L,L,k0ai,ni,0) .* TmL(L-1);
    end
    function f = dFL(L)
        f = dbesJ(M,L,L,k0ai,ni,0) + dbesH(M,L,L,k0ai,ni,0) .* TmL(L-1);
    end

end