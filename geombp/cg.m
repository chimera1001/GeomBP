function [obj_val,optimal_rhs,optimal_bins,UB,infeas] = cg(n,W,D,c,used_bins,forb_bins,PARAMS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLUMN GENERATION FOR BIN PACKING PROBLEM (BPP)
%
% Solves the following BPP instance:
%
%   Given:
%     - n ∈ ℕ : number of item types
%     - W ∈ ℕ^n : item weights
%     - D ∈ ℕ^n : item demands
%     - c ∈ ℕ : bin capacity
%
%   Objective:
%     Minimize the number of bins used to satisfy D,
%     using packing patterns x ∈ ℕ^n such that Wᵀx ≤ c.
%
%   Input parameters:
%     used_bins : fixed (must-use) patterns from branching
%     forb_bins : forbidden patterns from branching
%     PARAMS    : algorithmic parameters (tolerances, scaling)
%
%   Output:
%     obj_val       : optimal number of bins
%     optimal_rhs   : right-hand side of the final LP solution
%     optimal_bins  : pattern matrix associated to solution
%     UB            : upper bound from heuristics
%     infeas        : feasibility flag (true if demand infeasible)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Preprocessing

infeas = false;
if exist('used_bins', 'var') && ~isempty(used_bins)
    % Subtract total usage from demands
    usage = sum(used_bins, 2);
    D = D - usage;

    % Infeasibility check
    if any(D < 0)
        obj_val = NaN;
        optimal_rhs = [];
        optimal_bins = [];
        UB = NaN;
        infeas = true;
        return;
    end

    % Filter out forbidden bins involving now-zero-demand items
    if exist('forb_bins', 'var') && ~isempty(forb_bins)
        valid_idx = [];
        for k = 1:size(forb_bins, 2)
            if all(D >= forb_bins(:,k))
                valid_idx(end+1) = k;  % keep bin if still valid
            end
        end
        forb_bins = forb_bins(:, valid_idx);
    end

    % Keep only items with positive demand
    active = D > 0;
    D = D(active);
    W = W(active);
    if exist('forb_bins', 'var') && ~isempty(forb_bins)
        forb_bins = forb_bins(active, :);
    end

    n = length(D);
end


%% Run Heuristic UB1 (P = W)

UB1_bins = {};
tmpD = D;
P = W;
while true
    [x,~] = boundedKS(W, P, tmpD, c);
    % [x, ~] = bouknap(int64(W), int64(P), int32(tmpD), int64(c), int32(n));
    UB1_bins{end+1} = x;
    tmpD = tmpD - x;
    if sum(tmpD) == 0
        break
    end
end


%% Run Heuristic UB2 (P = W.^2)

UB2_bins = {};
tmpD = D;
P = W .^ 2;
while true
    [x,~] = boundedKS(W, P, tmpD, c);
    % [x, ~] = bouknap(int64(W), int64(P), int32(tmpD), int64(c), int32(n));
    UB2_bins{end+1} = x;
    tmpD = tmpD - x;
    if sum(tmpD) == 0
        break
    end
end


%% Run Heuristic UB3 (Gaussian-weighted W)

UB3_bins = {};
tmpD = D;
P = W .* exp((W - mean(W)).^2 / (2 * var(W)));
while true
    [x,~] = boundedKS(W, P, tmpD, c);
    % [x, ~] = bouknap(int64(W), int64(P), int32(tmpD), int64(c), int32(n));
    UB3_bins{end+1} = x;
    tmpD = tmpD - x;
    if sum(tmpD) == 0
        break
    end
end


%% Run Heuristic UB4 (P = W.^4)
UB4_bins = {};
tmpD = D;
P = W .^ 4;
while true
    [x,~] = boundedKS(W, P, tmpD, c);
    % [x, ~] = bouknap(int64(W), int64(P), int32(tmpD), int64(c), int32(n));
    UB4_bins{end+1} = x;
    tmpD = tmpD - x;
    if sum(tmpD) == 0
        break
    end
end


%% Compute Upper Bound

len1 = length(UB1_bins);
len2 = length(UB2_bins);
len3 = length(UB3_bins);
len4 = length(UB4_bins);

% Find best upper bound and its source
[UB, selector] = min([len1, len2, len3,len4]);

% Retrieve best cutting plan
switch selector
    case 1
        UB_patterns = UB1_bins;
    case 2
        UB_patterns = UB2_bins;
    case 3
        UB_patterns = UB3_bins;
    case 4
        UB_patterns = UB4_bins;
end


%% Aggregate All Patterns and Construct Matrix A

allPatterns = [UB1_bins, UB2_bins, UB3_bins,UB4_bins];
patternMatrix = cell2mat(cellfun(@(x) x(:), allPatterns, 'UniformOutput', false));
patternMatrix = reshape(patternMatrix, n, [])';

% Filter out singleton patterns (only one non-zero entry)
nonSingletonRows = sum(patternMatrix > 0, 2) > 1;
patternMatrix = patternMatrix(nonSingletonRows, :);

% Remove forbidden bins if provided
if exist('forb_bins', 'var') && ~isempty(forb_bins)
    is_forbidden = false(size(patternMatrix, 1), 1);
    for k = 1:size(patternMatrix, 1)
        for j = 1:size(forb_bins, 2)
            if norm(patternMatrix(k,:)' - forb_bins(:,j), Inf) < 1e-10
                is_forbidden(k) = true;
                break;
            end
        end
    end
    patternMatrix = patternMatrix(~is_forbidden, :);
end

% Extract unique non-singleton, non-forbidden patterns
[~, idx] = unique(patternMatrix, 'rows');
A = patternMatrix(idx, :)';

% Handle degenerate case: no non-singleton patterns found
if isempty(A)
    obj_val = n;
    optimal_bins = eye(n);
    optimal_rhs = D;

    if exist('active', 'var')
        obj_val = obj_val + size(used_bins, 2);  % include must-use bins in objective
        n0 = length(active);  % original number of items
        full_bins = zeros(n0, size(optimal_bins, 2));
        full_bins(active, :) = optimal_bins;
        optimal_bins = full_bins;
    end

    return;
end


%% Initialize Restricted Master Problem

invB = eye(n);
A = [eye(n), A];
% cCG = [ones(1, n), ones(1, size(A, 2) - n)];

idxB = 1:n;
idxNB = n+1:size(A, 2);

itRS = 0;
MAX_IT_RS = 10000;

while true
    pi = sum(invB,1);
    redcost_vec = pi * A(:, idxNB) - 1;
    [redCOST, idxMaxRed] = max(redcost_vec);

    if redCOST <= PARAMS.TOL_TERMINATION || itRS > MAX_IT_RS
        break;
    end

    entCOL = invB * A(:, idxNB(idxMaxRed));
    rhs = invB * D;

    best_key = [Inf, Inf];
    idxEXIT = 0;
    for ctr = 1:n
        if entCOL(ctr) > PARAMS.TOL_RATIO
            ratiocand = rhs(ctr) / entCOL(ctr);
            key = [ratiocand, ctr];
            if ratiocand < best_key(1) || (abs(ratiocand - best_key(1)) < 1e-10 && ctr < best_key(2))
                best_key = key;
                idxEXIT = ctr;
            end
        end
    end

    pivot = entCOL(idxEXIT);
    entCOL = -entCOL / pivot;
    entCOL(idxEXIT) = 1 / pivot;

    r1 = invB(idxEXIT,:);
    invB(idxEXIT,:) = 0;
    invB = invB + entCOL * r1;

    % Update basis
    entering = idxNB(idxMaxRed);
    leaving = idxB(idxEXIT);
    idxB(idxEXIT) = entering;
    idxNB(idxMaxRed) = leaving;

    itRS = itRS + 1;
end

if itRS > MAX_IT_RS
    warning("Initial simplex failed to converge within iteration limit.");
end

A = A(:, idxB);


%% Column Generation

MAX_IT = 500000;

itCG = 0;

while true
    % Compute and perturb duals
    pi = sum(invB, 1);

    % Normalize and scale duals
    norm_pi = norm(pi(pi > 0));

    P = floor(PARAMS.MULT_OBJ * pi / norm_pi);

    % Solve pricing problem
    [x, z] = combo(int64(W), int64(P), int64(c), int64(0), int64(0), int32(1), int32(0), int32(n));

    % Handle forbidden bins
    if ~isempty(forb_bins)
        is_forbidden = true;
        while is_forbidden
            is_forbidden = false;
            for k = 1:size(forb_bins, 2)
                if norm(forb_bins(:,k) - x, Inf) < 1e-10
                    idxP = find(P > 0);
                    z = z - min(P(idxP));
                    [x_sub, z] = twodks(int32(P(idxP)), int32(W(idxP)), int32(P(idxP)), c, z);
                    x = zeros(n,1);
                    x(idxP) = x_sub;
                    is_forbidden = true;
                    break;
                end
            end
        end
    end

    % Reduced cost check
    redCOST = pi * x - 1;
    if redCOST <= PARAMS.TOL_TERMINATION
        break;
    end

    % Entering column and pivot
    rhs = invB * D;
    entCOL = invB * x;

    % Lexicographic ratio test
    best_key = [Inf, Inf];
    idxEXIT = 0;
    for ctr = 1:n
        if entCOL(ctr) > PARAMS.TOL_RATIO
            ratiocand = rhs(ctr) / entCOL(ctr);
            key = [ratiocand, ctr];
            if ratiocand < best_key(1) || (abs(ratiocand - best_key(1)) < 1e-10 && ctr < best_key(2))
                best_key = key;
                idxEXIT = ctr;
            end
        end
    end

    if idxEXIT == 0
        warning('No valid pivot found. Terminating.');
        break;
    end

    % Pivot
    pivot = entCOL(idxEXIT);
    entCOL = -entCOL / pivot;
    entCOL(idxEXIT) = 1 / pivot;
    r1 = invB(idxEXIT,:);
    invB(idxEXIT,:) = 0;
    invB = invB + entCOL * r1;

    % Replace column
    A(:,idxEXIT) = x;

    % Iteration counter
    itCG = itCG + 1;
    if itCG > MAX_IT
        warning("Column generation stopped after max iterations.");
        break;
    end
end


% Final extract
obj_val = round(sum(invB,1) * D, PARAMS.DEC_ROUND_OBJ);
rhs = round(invB * D, PARAMS.DEC_ROUND_X);
selected_patterns = find(rhs > 0);
optimal_bins = A(:, selected_patterns);
optimal_rhs = rhs(selected_patterns);

if exist('active', 'var')
    obj_val = obj_val + size(used_bins, 2);
    n0 = length(active);
    full_bins = zeros(n0, size(optimal_bins, 2));
    full_bins(active, :) = optimal_bins;
    optimal_bins = full_bins;
end


end
