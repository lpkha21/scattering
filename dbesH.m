function dH = dbesH(M,I,J,k0ai,ni)
if(I == 0)
    n = 1;
else
    n = ni(I);
end
dH = zeros(2*M +1);
z = n*k0ai(J);
for m = -M:M
    
    dH(m + M + 1,m + M + 1) = m*besselh(m,z)/z - besselh(m+1,z);

end
end