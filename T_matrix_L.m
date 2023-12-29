function T = T_matrix_L(L,M,ka,kLa)
first = ni0*dH0(L+1) - niL*(dF(L + 1)\F(L+1))*H0(L+1);
second = ni0*dJ0(L+1) - niL*(dF(L + 1)\F(L+1))*J0(L+1);
T = -inv(first) * second;

% function definitions
% -------------------
% besJi - J(i,i)
% besjL1 - J(L+1,i)
% besHL1 - H(L+1,i)
% dbesJi - dJ(i,i)
% dbesJL1 - dJ(L+1,i)
% dbesHL1 - dH(L+1,i)


end