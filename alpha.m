function a = alpha(M,k0di,Thi,J,I)
% alpha_{J,I}: outgoing waves of cylinder I re-expanded as regular waves about
% cylinder J (index 0 = host centre). k0di must already be k_host*d
% (T_matrix_L_El / T_matrix_L_mag do this scaling).
a = zeros(2*M + 1);

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
            a(n + M + 1,m + M + 1) = besselh(n-m,1,dist) * exp(1j*(m-n)*ang);
        end
    end
    return;
end
% I < J
if (I == 0)
    dist = k0di(J);
    ang = Thi(J);
else
    dist = dij(k0di(I), k0di(J), Thi(I), Thi(J));
    ang = fij(k0di(I), k0di(J), Thi(I), Thi(J));
end
for m = -M:M
    for n = -M:M
        a(n + M + 1,m + M + 1) = (-1)^(n-m)*besselh(n-m,1,dist) * exp(1j*(m-n)*ang);
    end
end

end
