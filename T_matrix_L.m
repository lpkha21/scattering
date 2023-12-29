function T = T_matrix_L(L,M,k0ai,ni,k0di)
first = ni0*dH0(L+1,1,1) - niL*(dFL(L + 1)\FL(L+1))*H0(L+1,1,1);
second = ni0*dJ0(L+1,1,1) - niL*(dFL(L + 1)\FL(L+1))*J0(L+1,1,1);
T = -inv(first) * second;

% function definitions
% -------------------
% besJi - J(i,i)
% besjL1 - J(L+1,i)
% besHL1 - H(L+1,i)
% dbesJi - dJ(i,i)
% dbesJL1 - dJ(L+1,i)
% dbesHL1 - dH(L+1,i)
% H0 - H(0, L+1)
% J0 - J(0, L + 1)
% dH0 - dH(0, L + 1)
% dJ0 - dJ(0, L + 1)
    



    
    function f = FL(L)
        f = J0(L,ni(L)) + H0(L,ni(L)) .* TL(L-1);
    end
    function f = dFL(L)
        f = dJ0(L,ni(L)) + dH0(L,ni(L)) .* TL(L-1);
    end
    
end






