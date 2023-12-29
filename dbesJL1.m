function J = dbesJL1(I,M,kLa)
J = zeros(2*M + 1);
for m = -M:1:M
    for n = -M:1:M
        J(m + M + 1,n + M + 1) = m*besselj(m,kLa(I))/kLa(I) - besselj(m+1,kLa(I));
    end
end

end