function [H, Q] = arnoldi(A, Q, k)
    H = zeros(k + 1, k);
    for j = 1:k
        z = A*Q(:, j);
        for i = 1:j
            H(i, j) = ((Q(:, i)).') * z;
            z = z - H(i, j) * Q(:, i);
        end
        H(j + 1, j) = norm(z);
        if H(j + 1, j) == 0
            return
        end
        Q(:, j + 1) = z / H(j + 1, j);
    end
end