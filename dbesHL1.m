function H = dbesHL1(I,M,kLa)
H = zeros(2*M + 1);
for m = -M:1:M
    for n = -M:1:M
        H(m + M + 1,n + M + 1) = m*besselh(m,kLa(I))/kLa(I) - besselh(m+1,kLa(I));
    end
end
end