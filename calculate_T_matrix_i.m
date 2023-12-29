function T = calculate_T_matrix_i(i, ka, L, th, M)
    T = zeros(length(i)); % Initialize T as a square matrix

    ni0 = sqrt(8.8541878128 * 10e-12 / (1.25663706212 * 10e-6));

    % Define n and m
    n = 0:M;
    m = 0:M;

    for p = 1:length(i)
        % Define beta and T1 functions outside the loop
        
        T1 = @(ii, l) calculate_T_matrix_i(ii, ka, l, th, M);

        % Initialize matrices
        J = zeros((M + 1) * L, (M + 1) * L);
        dJ = zeros((M + 1) * L, (M + 1) * L);
        H = zeros((M + 1) * L, (M + 1) * L);
        dH = zeros((M + 1) * L, (M + 1) * L);

        for ii = 1:L
            % Calculate J, dJ, H, dH using element-wise operations
            J((ii - 1) * (M + 1) + 1:ii * (M + 1), (ii - 1) * (M + 1) + 1:ii * (M + 1)) = besselj(m, ka(ii)) .* besselj(m', ka(ii));
            dJ((ii - 1) * (M + 1) + 1:ii * (M + 1), (ii - 1) * (M + 1) + 1:ii * (M + 1)) = 0.5 * (besselj(m - 1, ka(ii)) - besselj(m + 1, ka(ii))) .* besselj(m', ka(ii));
            H((ii - 1) * (M + 1) + 1:ii * (M + 1), (ii - 1) * (M + 1) + 1:ii * (M + 1)) = besselh1(m, ka(ii)) .* besselj(m', ka(ii));
            dH((ii - 1) * (M + 1) + 1:ii * (M + 1), (ii - 1) * (M + 1) + 1:ii * (M + 1)) = 0.5 * (besselh1(m - 1, ka(ii)) - besselh1(m + 1, ka(ii))) .* besselj(m', ka(ii));
        end

        % Calculate T-matrix element using the given formula
        T(p) = - (ni0 * dH(1, i(p)) - ni0 * dH(i(p)*(L-1)+1, i(p)) * (H(i(p)*(L-1)+1:i(p)*L, i(p)*(L-1)+1:i(p)*L) \ H(1:i(p), i(p)))) \ ...
            (ni0 * dJ(1, i(p)) - ni0 * dH(i(p)*(L-1)+1, i(p)) * (H(i(p)*(L-1)+1:i(p)*L, i(p)*(L-1)+1:i(p)*L) \ J(1:i(p), i(p))));
    end

    function result = besselh1(m, x)
        % Bessel function of the third kind (Hankel function) of order 1
        result = besselh(m, 1, x);
    end
end
