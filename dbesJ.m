function dJ = dbesJ(M,I,J,k0ai,ni)
if(I == 0)
    n = 1;
else
    n = ni(I);
end
dJ = zeros(2*M +1);
for m = -M:1:M
    if(m < 0)
        dJ(m + M + 1,m + M + 1) = 0.5 * (-1)^(-m) * (besselj(-m-1,n*k0ai(J)) - besselj(-m+1,n*k0ai(J)));
    else
        dJ(m + M + 1,m + M + 1) = 0.5 * (besselj(m-1,n*k0ai(J)) - besselj(m+1,n*k0ai(J)));
    end
end
end