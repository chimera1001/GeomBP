% =========================================================================
%   GeomBP.m,    Masoud Ataei                                 July 2025
% =========================================================================
%
%   This script implements the Geometrical Branch-and-Price (GeomBP) algorithm
%   for solving the Bin Packing Problem (BPP), as introduced in:
%
%     M. Ataei and S. Chen,
%     "A Geometrical Branch-and-Price Algorithm for the Bin Packing Problem",
%     University of Toronto & York University, 2025.
%
%   This implementation is designed for academic use only.
%
%   Correspondence:
%     Dr. Masoud Ataei
%     Department of Mathematical and Computational Sciences
%     University of Toronto Mississauga
%     Email: masoud.ataei@utoronto.ca
%
%   License:
%     This code is provided for academic and non-commercial purposes only.
%
% =========================================================================


%% Load data

clc; clear;

data = readmatrix("../Data/Falkenauer_CSP/Falkenauer_T/Falkenauer_t60_00.txt", 'Range', 'A1');

n = data(1,1);
c = data(2,1);
W = data(3:end,1);
D = data(3:end,2);

clear data;


%% Parameters

PARAMS.MULT_OBJ = 1e8;
PARAMS.TOL_RATIO = 1e-8;
PARAMS.TOL_TERMINATION = 1e-9;
PARAMS.DEC_ROUND_OBJ = 5;
PARAMS.DEC_ROUND_X = 4;
PARAMS.BIG_M = 1e6;
BATCH_LIMIT = 100;


%% Branch and Price Algorithm

% Root node
root.used = []; root.forb = [];
[root.obj, root.rhs, root.bins, UB, infeas] = cg(n,W,D,c,root.used,root.forb,PARAMS);
if infeas, error('Root infeasible'); end

best_UB = UB; best_bins = root.bins; best_rhs = root.rhs;
nodes = struct('node',root,'score',root.obj,'label','R');
nodeCount = 0;

while ~isempty(nodes)
    scores = [nodes.score];
    labels = {nodes.label};
    [~, order] = sortrows([scores(:), label_to_priority(labels)']);
    current_struct = nodes(order(1)); nodes(order(1)) = [];
    current = current_struct.node; current_label = current_struct.label;
    nodeCount = nodeCount + 1;
    fprintf('Node%3d (%s): LB=%.4f, UB=%d\n', nodeCount, current_label, current.obj, best_UB);

    if ceil(current.obj) >= best_UB
        fprintf(' → pruned\n'); continue;
    end

    if all(abs(current.rhs - round(current.rhs)) < 1e-8)
        if round(current.obj) < best_UB
            best_UB = round(current.obj);
            best_bins = [current.used, current.bins];
            best_rhs = [ones(1,size(current.used,2)), round(current.rhs)'];
            fprintf(' → New UB from integrality = %d\n', best_UB);
            nodes = nodes(arrayfun(@(x) ceil(x.score)<best_UB, nodes));
        end
        continue;
    end

    if ~strcmp(current_label,'D')
        fprintf(' → randomized diving...\n');
        idxs = randperm(size(current.bins,2)); used_now=current.used;
        for k = idxs
            if isempty(used_now)
                used_now = current.bins(:,k);
            else
                tmp = sum(used_now,2)+ current.bins(:,k);
                if all(tmp <= D + 1e-8)
                    used_now = [used_now, current.bins(:,k)];
                end
            end
        end
        [child.obj,child.rhs,child.bins,~,infeas] = cg(n,W,D,c,used_now,current.forb,PARAMS);
        if ~infeas && ceil(child.obj) < best_UB
            child.used = used_now;
            child.forb = current.forb;
            nodes(end+1)=struct('node',child,'score',child.obj,'label','D');
        end
    end

    fr = find(abs(current.rhs - round(current.rhs))>1e-8);
    [~,jpos] = max(current.rhs(fr)); j = fr(jpos);
    xj = current.bins(:,j);

    % Left: enforce
    nodeL.used = [current.used, xj];
    nodeL.forb = current.forb;
    [nodeL.obj,nodeL.rhs,nodeL.bins,~,infeasL] = cg(n,W,D,c,nodeL.used,nodeL.forb,PARAMS);
    if ~infeasL && ceil(nodeL.obj)<best_UB
        nodeL.used = [current.used, xj];
        nodeL.forb = current.forb;
        nodes(end+1) = struct('node',nodeL,'score',nodeL.obj,'label','E');
    end

    % Right: forbid
    isf = false;
    if ~isempty(current.forb)
        isf = any(all(current.forb == xj, 1));
    end
    if ~isf
        nodeR.used = current.used;
        nodeR.forb = [current.forb, xj];
        [nodeR.obj,nodeR.rhs,nodeR.bins,~,infeasR] = cg(n,W,D,c,nodeR.used,nodeR.forb,PARAMS);
        if ~infeasR && ceil(nodeR.obj)<best_UB
            nodeR.used = current.used;
            nodeR.forb = [current.forb, xj];
            nodes(end+1)=struct('node',nodeR,'score',nodeR.obj,'label','F');
        end
    end
end

fprintf('\nDONE: best integer = %d bins\n', best_UB);

function p = label_to_priority(labels)
    map = struct('D',1,'E',2,'F',3,'R',4);
    p = zeros(size(labels));
    for i=1:numel(labels), p(i)=map.(labels{i}); end
end