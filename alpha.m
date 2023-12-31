function a = alpha(M,k0di,Thi,J,I)
a = zeros(2*M + 1);

if (I ~= 0 && J ~= 0)
    dist = dij(k0di(J),k0di(I),Thi(J),Thi(I));
    ang = fij(k0di(J),k0di(I),Thi(J),Thi(I));
end
if(J == 0)
    dist = k0di(I);
    ang = Thi(I);
end
if(I == 0)
    dist = k0di(J);
    ang = mod((pi + Thi(J)),2*pi);
end
for m = -M:M
    for n = -M:M
        a(m+M + 1,n + M + 1) = besselh(n - m, dist) * exp(-1j*(n - m) * ang);

    end
end
end