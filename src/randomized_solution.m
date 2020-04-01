function [vector] = randomized_solution(matrix)
    % Function:
    %   - generate a randomized complex vector for homogeneous QCQP with two double-sided constraints
    %
    % InputArg(s):
    %   - matrix [Q] (N * N): a Hermitian matrix
    %
    % OutputArg(s):
    %   - vector [v] (N * 1): a vector that satisfies (1) trace(\Lambda * v * v') = trace(\Lambda); (2) trace(Q * v * v') = trace(Q); (3) v' * v = N
    %
    % Comment(s):
    %   - solve the complex version of the QCQP with two double-sided constraints
    %   - \Lambda can be any real diagonal matrix and is set to I
    %
    % Reference(s):
    %   - Y. Huang and D. P. Palomar, "Randomized Algorithms for Optimal Solutions of Double-Sided QCQP With Applications in Signal Processing," IEEE Transactions on Signal Processing, vol. 62, no. 5, pp. 1093–1108, 2014.
    %
    % Author & Date: Yang (i@snowztail.com) - 26 Mar 20



    Q1 = real(matrix);
    Q2 = imag(matrix);
    delta2 = trace(matrix);
    isValid = false;
    while ~isValid
        xi = 1 - 2 * randi([0, 1], [length(matrix), 1]);
        eta = 1 - 2 * randi([0, 1], [length(matrix), 1]);
        isValid = (xi.' * Q1 * xi - delta2) * (eta.' * Q1 * eta - delta2) <= 0;
    end
    if xi.' * Q1 * xi == delta2
        vector = xi;
    elseif eta.' * Q1 * eta == delta2
        vector = eta;
    else
        gamma0 = (xi.' * Q2 * eta + sqrt((eta.' * Q2 * xi) ^ 2 - (eta.' * Q1 * eta - delta2) * (xi.' * Q1 * xi - delta2))) / (xi.' * Q1 * xi - delta2);
        vector = gamma0 / sqrt(1 + gamma0 ^ 2) * xi + 1i / sqrt(1 + gamma0 ^ 2) * eta;
    end

end
