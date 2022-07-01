function [Q, R] = classic_gram(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(m, n);
    for j = 1 : n
        Q(:, j) = A(:, j);
        for i = 1 : j - 1
            R(i, j) = Q(:, i)' * A(:, j);
            Q(:, j) = Q(:, j) - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(Q(:, j));
        Q(:, j) = Q(:, j) / R(j, j);
    end
end