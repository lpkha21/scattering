function Tl = TmLMag(L ,l, M, k0ai,k0di, Thi, si, ni)
% TM: returns T_{si(l)} * beta_{si,0}, the T-matrix of cylinder si in the presence
% of the first l cylinders (recursive algorithm, paper Eqs. 9-10).
% k0di must already be k_host*d.

% base case
if(l == 1)
    Tl = T_matrix_base_case_mag(L,M,k0ai,ni,si)*beta(M,k0di,Thi,si,0);
    return;
end

% recursive step: compute the (l-1)-cylinder results once and reuse them
Tprev = cell(1, l-1);
for I = 1:(l-1)
    Tprev{I} = TmLMag(L,l-1,M,k0ai,k0di,Thi,I,ni);
end
first_sum = zeros(2*M + 1);
second_sum = zeros(2*M + 1);
for I = 1:(l-1)
    first_sum  = first_sum  + alpha(M,k0di,Thi,l,I) * Tprev{I} * beta(M,k0di,Thi,0,I) * alpha(M,k0di,Thi,I,l);
    second_sum = second_sum + alpha(M,k0di,Thi,l,I) * Tprev{I};
end
Tl1 = T_matrix_base_case_mag(L,M,k0ai,ni,l);
TLB = (eye(2*M + 1) - Tl1*first_sum) \ (Tl1*(beta(M,k0di,Thi,l,0) + second_sum));   % Eq. (9)

if(si == l)
   Tl = TLB;
   return;
end
Tl = Tprev{si} + Tprev{si}*beta(M,k0di,Thi,0,si)*alpha(M,k0di,Thi,si,l)*TLB;   % Eq. (10)
end
