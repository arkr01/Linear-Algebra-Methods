function [x, rel_res, num] = my_jacobi(A, b, x, max_iter)
    d = diag(A);
    epsilon = 10^(-8);
    num = 0;
    rel_res = zeros(1, max_iter);
    while true
        num = num + 1;
        r = b - A * x;
        x = x + r ./ d;
        rel_res(num) = norm(r) / norm(b);
        if rel_res(num) <= epsilon || num == max_iter
            break;
        end
    end
end