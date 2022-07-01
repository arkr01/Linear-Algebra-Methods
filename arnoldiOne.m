function [h_next, Q] = arnoldiOne(A, Q, k)
    z = A*Q(:, k);
    h_next = zeros(k + 1, 1);
    for i = 1:k
        h_next(i) = ((Q(:, i)).') * z;
        z = z - h_next(i) * Q(:, i);
    end
    h_next(k + 1) = norm(z);
    if h_next(k + 1) == 0
        return
    end
    Q(:, k + 1) = z / h_next(k + 1);
end