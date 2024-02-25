function a = alpha(M,k0di,Thi,J,I)
digits(100)
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
            a(n + M + 1,m + M + 1) = vpa(besselh(n-m,dist) * exp(1j*(m-n)*ang));
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
        a(n + M + 1,m + M + 1) = vpa((-1)^(n-m)*besselh(n-m,dist) * exp(1j*(m-n)*ang));
    end
end

end