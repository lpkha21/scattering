function Tl = TmLMag(L, l, M, k0ai, k0di, Thi, si, ni)
digits(50)
% Base case
if l == 1
    Tl = T_matrix_base_case_mag(L, M, k0ai, ni, si) * beta(M, k0di, Thi, si, 0);
    return;
end

% Recursive step
first_sum = zeros(2*M + 1); % Initialize for efficiency
second_sum = zeros(2*M + 1);

for I = 1:(l-1)
    alpha_term = vpa(alpha(M, k0di, Thi, l, I));
    T_prev = TmLMag(L, l-1, M, k0ai, k0di, Thi, I, ni);  
    first_sum = first_sum + alpha_term * T_prev * beta(M, k0di, Thi, 0, I) * alpha_term;
    second_sum = second_sum + alpha_term * T_prev;
end

% Use LU decomposition for more robust inversion
[L_matrix, U_matrix, P_matrix] = lu(eye(2*M + 1) - T_matrix_base_case_mag(L, M, k0ai, ni, l) * first_sum);
TLB = U_matrix \ (L_matrix \ (P_matrix * (T_matrix_base_case_mag(L, M, k0ai, ni, l) * (beta(M, k0di, Thi, l, 0) + second_sum)))); % Use backslash '\' instead

% Check if calculation is for the final layer
if si == l
    Tl = vpa(TLB);
    return; 
end

Tlminus1 = TmLMag(L, l - 1, M, k0ai, k0di, Thi, si, ni);
Tl = vpa(Tlminus1 + Tlminus1 * beta(M, k0di, Thi, 0, si) * alpha(M, k0di, Thi, si, l) * TLB);

end
