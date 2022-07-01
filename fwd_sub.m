function x = fwd_sub(A, b)
    x = zeros(size(b));
    n = length(b);
    for i = 1 : n
        x(i) = (b(i) - sum(A(i, 1 : i - 1) * x(1 : i - 1))) / (A(i, i));
    end
end