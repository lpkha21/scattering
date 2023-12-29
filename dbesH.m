function dH = dbesH(I,J,k0ai,ni,z)
n = ni(J);
if(z == 1)
    n = 1;
end
dH = zeros(2*M +1);
for m = -M:1:M
    dH(m + M + 1,m + M + 1) = m*besselh(m,n*k0ai(I))/(n*k0ai(I)) - besselh(m+1,n*k0ai(I));
end
end