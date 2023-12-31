function Tl = TmLMag(L ,l, M, k0ai,k0di, Thi, si, ni)
% % Tlminus1 = T_matrix_base_case_mag(L,M,k0ai,ni,si)*beta(M,k0di,Thi,si,0);
% % Tl = Tlminus1;
% % for l = 2:L
% %     first_sum = 0;
% %     for I = 1:(L-1)
% %         first_sum = first_sum + alpha(M,k0di,Thi,L,I) * TmL(L-1,M,k0ai,k0di,Thi,I,ni) * beta(M,k0di,Thi,0,I)*alpha(M,k0di,Thi,I,L);
% %     end
% %     second_sum = 0;
% %     for I = 1:(L-1)
% %         second_sum = second_sum + alpha(M,k0di,Thi,L,I) * TmL(L-1,M,k0ai,k0di,Thi,I,ni);
% %     end
% %     
% %     first_part = eye(2*M + 1) - T_matrix_base_case_mag(L,M,k0ai,ni,L)*first_sum;
% %     second_part = T_matrix_base_case_mag(L,M,k0ai,ni,L)* ( beta(M,k0di,Thi,L,0) + second_sum);
% %     
% %     TLB = first_part \ second_part;
% %     Tl = Tlminus1  + Tlminus1*beta(M,k0di,Thi,0,si)*alpha(M,k0di,Thi,si,L)*TLB;
% %     Tlminus1 = Tl;
% % end

% base case
if(l == 1)
    Tl = T_matrix_base_case_mag(L,M,k0ai,ni,si)*beta(M,k0di,Thi,si,0);
    return;
end

% recursive step


% --------------------------------
first_sum = 0;
parfor I = 1:(l-1)
    first_sum = first_sum + alpha(M,k0di,Thi,l,I) * TmLMag(L,l-1,M,k0ai,k0di,Thi,I,ni) * beta(M,k0di,Thi,0,I)*alpha(M,k0di,Thi,I,l);
end
second_sum = 0;
parfor I = 1:(l-1)
    second_sum = second_sum + alpha(M,k0di,Thi,l,I) * TmLMag(L,l-1,M,k0ai,k0di,Thi,I,ni);
end

first_part = eye(2*M + 1) - T_matrix_base_case_mag(L,M,k0ai,ni,l)*first_sum;

second_part = T_matrix_base_case_mag(L,M,k0ai,ni,l)* ( beta(M,k0di,Thi,l,0) + second_sum);

TLB = vpa(first_part )\ vpa(second_part );
if(si == l)
   Tl = TLB;
   return;
end
Tlminus1 = TmLMag(L,l - 1,M,k0ai,k0di,Thi,si,ni);

Tl = Tlminus1  + Tlminus1*beta(M,k0di,Thi,0,si)*alpha(M,k0di,Thi,si,l)*TLB;
% Tl
end