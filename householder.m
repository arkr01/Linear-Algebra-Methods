function [Q, R] = householder(A)
    [m, n] = size(A);
    R = A;
    Q = eye(m);
    for i = 1 : n
        A(i : m, i) = R(i : m, i);
        U = A(i : m, i) + sign(A(i, i)) * norm(A(i : m, i)) * ...
            eye(size(A(i : m, i)));
        U = U / norm(U);
        reflector = (eye(m - i + 1)) - 2 * (U * U.');
        R(i : m, i : n) = reflector * R(i : m, i : n);
        Q = [eye(i - 1), zeros(i - 1, m - i + 1); ...
            zeros(m - i + 1, i - 1), reflector] * Q;
    end
    Q = Q.';
end