function H = besH(M,I,J,k0ai,ni)
if(I == 0)
    n = 1;
else
    n = ni(I);
end
H = zeros(2*M +1);
for m = -M:M
    H(m + M + 1,m + M + 1) = besselh(m,n*k0ai(J));
end
end