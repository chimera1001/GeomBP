function [x_bounded,fval] = boundedKS(W, P, D, c)
    % W: weights (n x 1)
    % P: profits (n x 1)
    % D: demand vector (n x 1)
    % c: knapsack capacity (scalar)
    % Returns: x_bounded (n x 1), the bounded knapsack solution

    W = W(:); P = P(:); D = D(:);
    n = length(W);

    % Expand into binary items
    W_exp = []; P_exp = []; map = [];
    for i = 1:n
        W_exp = [W_exp; repmat(W(i), D(i), 1)];
        P_exp = [P_exp; repmat(P(i), D(i), 1)];
        map = [map; repmat(i, D(i), 1)];
    end

    % Call combo on binary knapsack problem
    m = length(W_exp);
    [x_bin, ~] = combo(int64(W_exp), int64(P_exp), int64(c), ...
                       int64(0), int64(0), int32(1), int32(0), int32(m));

    % Reconstruct bounded solution
    x_bounded = zeros(n, 1);
    for j = 1:m
        if x_bin(j)
            x_bounded(map(j)) = x_bounded(map(j)) + 1;
        end
    end

    fval = P'*x_bounded;
end
