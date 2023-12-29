function Tl = TmL(L , M, k0ai,k0di, Thi, si, ni)
% recursive function counting T(si,(L)) * beta(si,0)

% base case
if(L == 1)
    Tl = T_matrix_base_case_mag(L,M,k0ai,ni,si);
    Tl = Tl + eye(size(Tl));
    return;
end

% si -> summation index [1;L]

% recursive step
Tlminus1 = TmL(L - 1,M,k0ai,k0di,Thi,si,ni);

% --------------------------------
first_sum = 0;
for I = 1:L
    first_sum = first_sum + alpha(M,k0di,Thi,L,si) * Tlminus1 * beta(M,k0di,Thi,0,si)*alpha(M,k0di,Thi,si,L);
end
second_sum = 0;
for I = 1:L
    second_sum = second_sum + alpha(M,k0di,Thi,L,si) * Tlminus1;
end


first_part = eye(2*M + 1) - T_matrix_base_case_mag(L,M,k0ai,ni,L)*first_sum;
second_part = T_matrix_base_case_mag(L,M,k0ai,ni,L)* ( beta(M,k0di,Thi,L,0) + second_sum);

TLB = first_part \ second_part;
Tl = Tlminus1  + Tlminus1*beta(M,k0di,Thi,0,si)*alpha(M,k0di,Thi,si,L)*TLB;

end

