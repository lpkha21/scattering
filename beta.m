function b = beta(M,k0di,Thi,J,I)
% beta_{J,I}: translation of regular waves between cylinder I and cylinder J
% (index 0 = host centre). k0di must already be k_host*d.
if J == I
    b = eye(2*M + 1);
    return;
end
b = zeros(2*M + 1);

if J < I
    if (J == 0)
        dist = k0di(I);
        ang = Thi(I);
    else
        dist = dij(k0di(J), k0di(I), Thi(J), Thi(I));
        ang = fij(k0di(J), k0di(I), Thi(J), Thi(I));
    end
    for m = -M:M
        for n = -M:M
            b(n + M + 1,m + M + 1) = besselj(n-m,dist) * exp(1j*(m-n)*ang);
        end
    end
    return;
end
% I < J: inverse translation
c = beta(M, k0di, Thi, I, J);
[L, U, P] = lu(c);

b = U \ (L \ P);

end
