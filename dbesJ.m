function dJ = dbesJ(M,I,J,k0ai,ni)
if(I == 0)
    n = 1;
else
    n = ni(I);
end
dJ = zeros(2*M +1);
z = n*k0ai(J);
for m = -M:1:M
    dJ(m + M + 1,m + M + 1) = m*besselj(m,z)/z - besselj(m+1,z);
end
end