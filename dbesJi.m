function J = dbesJi(I,M,ka)
J = zeros(2*M + 1);
for m = -M:1:M
    for n = -M:1:M
        J(m + M + 1,n + M + 1) = m*besselj(m,ka(I))/ka(I) - besselj(m+1,ka(I));
    end
end

end