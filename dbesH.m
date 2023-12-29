function dH = dbesH(M,I,J,k0ai,ni)
if(I == 0)
    n = 1;
else
    n = ni(I);
end
dH = zeros(2*M +1);
for m = -M:1:M
    if m < 0
        dH(m + M + 1,m + M + 1) = (-1)^(-m) * (-m*besselh(-m,n*k0ai(J))/(n*k0ai(J)) - besselh(-m+1,n*k0ai(J)));
    else
        dH(m + M + 1,m + M + 1) = m*besselh(m,n*k0ai(J))/(n*k0ai(J)) - besselh(m+1,n*k0ai(J));
    end
end
end