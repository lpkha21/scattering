function a = alpha(M,k0di,Thi,J,I)
a = zeros(2*M + 1);

if J < I
    if (J == 0)
        dist = k0di(I);
        ang = Thi(I);
    else
        dist = vpa(dij(k0di(J), k0di(I), Thi(J), Thi(I)));
        ang = vpa(fij(k0di(J), k0di(I), Thi(J), Thi(I)));
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
    dist = vpa(dij(k0di(I), k0di(J), Thi(I), Thi(J)));
    ang = vpa(fij(k0di(I), k0di(J), Thi(I), Thi(J)));
end
for m = -M:M
    for n = -M:M
        a(n + M + 1,m + M + 1) = (-1)^(n-m)*besselh(n-m,1,dist) * exp(1j*(m-n)*ang);
    end
end

end