function bes = besJ(I,J,k0ai,ni,z)
n = ni(J);
if(z == 1)
    n = 1;
end
bes = zeros(2*M +1);
for m = -M:1:M
    bes(m + M + 1,m + M + 1) = besselj(m,n*k0ai(I));
end
end