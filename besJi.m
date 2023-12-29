function J = besJi(I,M,ka)
J = zeros(2*M + 1);
for m = -M:1:M
    for n = -M:1:M
        J(m + M + 1,n + M + 1) = besselj(m,ka(I));
    end
end

end