function C = mykmeans(X, k, tol)
    n = size(X, 1);
    C_old = X(randperm(n, k), :);
    C = zeros(size(C_old));
    init = 0;

    while diag(pdist2(C, C_old)) > tol
        if init ~= 0
            C_old = C;
        end
        init = 1;
        D = pdist2(X, C_old);
        [~, ownerships] = min(D, [], 2);
        
        for j = 1 : k
            C(j, :) = mean(X(ownerships == j, :));
        end
    end
end