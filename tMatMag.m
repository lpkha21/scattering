function [t] = tMatMag(l, Eph, ka, kd, L, th,n3)
    % Eph - epsilon vector
    % ka - k(i) * a(i) vector - (a - radius)
    % kd - k(i) * d(i)
    % L - total N of inclusion
    % l - curr # of inclusion
    % th - angle vector

    ni0 = sqrt(8.8541878128 * 10e-12 / (1.25663706212 * 10e-6));
    niL = sqrt(Eph(l) / (1.25663706212 * 10e-6));
    M = 12;
    % Define n and m
    n = 0:M;
    m = 0:M;

    % T1 - Ti(l)
    % Jn-m(k0*dij) -> besselj(n-m,k0*dij);
    beta = @(x, y) besselj(n - m, kd(L - l + 1)) .* exp(-1j * (n - m) * th(L - l + 1));

    % Calculate J, H, dJ, dH using Kronecker delta
    J = zeros((M+1)*L, (M+1)*L);
    dJ = zeros((M+1)*L, (M+1)*L);
    H = zeros((M+1)*L, (M+1)*L);
    dH = zeros((M+1)*L, (M+1)*L);

    for ii = 1:L
        m = 0:M;

        % Update J, dJ, H, dH using Kronecker delta
        rows = (ii-1)*(M+1)+1:ii*(M+1);
        cols = (ii-1)*(M+1)+1:ii*(M+1);
        J(rows, cols) = besselj(m, ka(ii)) .* besselj(m', ka(ii));
        dJ(rows, cols) = besselj_prime(m, ka(ii)) .* besselj(m', ka(ii));
        H((ii-1)*(M+1)+1:ii*(M+1), (ii-1)*(M+1)+1:ii*(M+1)) = besselh1(m, ka(ii)) .* besselj(m', ka(ii));
        dH((ii-1)*(M+1)+1:ii*(M+1), (ii-1)*(M+1)+1:ii*(M+1)) = besselh1_prime(m, ka(ii)) .* besselj(m', ka(ii));
    end
    i = 1:M;  % i is a vector representing the range from 1 to M

    T1 = @(ii, l) calculate_T_matrix_i(ii, ka, l, th, M);

    T_l = sum(beta(0, i) .* (T1(i, l).' .* beta(i, 0)));

    T_l_replicated = repmat(T_l, size(H,1) / size(T_l,1) , size(H,2) / size(T_l,2));

    F = J + H .* T_l_replicated;
    dF = dJ + dH .* T_l_replicated;

    epsilon = 1e-5;  % Small regularization parameter
    F = F + epsilon * eye(size(F));
    dF = dF + epsilon * eye(size(dF));
    % Calculate T-matrix element using the given formula
    q1 = (ni0 * dH(1, l) - niL * dF * (F \ H));
    q2 = (ni0 * dJ(1, l) - niL * dF * (F \ J));

    t = -  q1 .\ q2;

    function result = besselj_prime(m, x)
        % Derivative of Bessel function of the first kind with respect to its argument
        result = 0.5 * (besselj(m - 1, x) - besselj(m + 1, x));
    end

    function result = besselh1(m, x)
        % Bessel function of the third kind (Hankel function) of order 1
        result = besselh(m, 1, x);
    end

    function result = besselh1_prime(m, x)
        % Derivative of Bessel function of the third kind (Hankel function) of order 1
        result = 0.5 * (besselh(m - 1, 1, x) - besselh(m + 1, 1, x));
    end
end
