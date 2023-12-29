function H = besH(I,J,k0ai,ni,z)
n = ni(J);
if(z == 1)
    n = 1;
end

H = zeros(2*M +1);
for m = -M:1:M
    H(m + M + 1,m + M + 1) = besselh(m,n*k0ai(I));
end
end