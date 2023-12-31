function b = besJ(M,I,J,k0ai,ni)
if(I == 0)
    n = 1;
else
    n = ni(I);
end
b = zeros(2*M +1);
for m = -M:1:M
    b(m + M + 1,m + M + 1) = besselj(m,n*k0ai(J));
end

end