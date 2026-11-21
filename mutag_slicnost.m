1;

function Z = blondelZ(A, B, Z0, maxIter, tol)
% Blondel iteration with a given starting matrix Z0 (size(B,1) x size(A,1)).
% Even number of steps, normalized by the Frobenius norm each step.
    Z = Z0 / norm(Z0, 'fro');
    for it = 1:maxIter
        Zold = Z;
        for step = 1:2
            Z = B*Z*A' + B'*Z*A;
            Z = Z / norm(Z, 'fro');
        end
        if norm(Z - Zold, 'fro') < tol
            break;
        end
    end
end

function s = graphScore(A, B, la, lb, useLabels, maxIter, tol)
% One number describing how well G_A and G_B match (not yet normalized).
    if useLabels
        Z0 = double(lb(:) == la(:)');           % 1 where atom types match
    else
        Z0 = ones(size(B,1), size(A,1));
    end
    if norm(Z0, 'fro') == 0
        s = 0; return;                          % no common atom types at all
    end
    Z  = blondelZ(A, B, Z0, maxIter, tol);
    Zs = Z * sqrt(numel(Z));                    % remove size effect
    s  = (mean(max(Zs, [], 1)) + mean(max(Zs, [], 2))) / 2;   % symmetric best match
end

%% ===== Load MUTAG (TU Dortmund format) =====
useLabels = true;
maxIter   = 100;
tol       = 1e-9;

edges = dlmread('MUTAG_A.txt', ',');            % global node ids, both directions
gi    = load('MUTAG_graph_indicator.txt');      % graph id of each node
nl    = load('MUTAG_node_labels.txt');          % atom type of each node
nG    = max(gi);
Aall  = sparse(edges(:,1), edges(:,2), 1, numel(gi), numel(gi));

graphs = cell(nG, 1);
labs   = cell(nG, 1);
for g = 1:nG
    nodes     = find(gi == g);
    graphs{g} = full(Aall(nodes, nodes));
    labs{g}   = nl(nodes);
end
fprintf('Loaded %d molecules.\n', nG);

%% ===== Self scores (for normalization) =====
selfScore = zeros(nG, 1);
for g = 1:nG
    selfScore(g) = graphScore(graphs{g}, graphs{g}, labs{g}, labs{g}, useLabels, maxIter, tol);
end

%% ===== Pairwise similarity matrix between molecules =====
Sim = eye(nG);
for g = 1:nG
    for h = g+1:nG
        s = graphScore(graphs{g}, graphs{h}, labs{g}, labs{h}, useLabels, maxIter, tol);
        s = s / sqrt(selfScore(g) * selfScore(h));
        Sim(g, h) = s;
        Sim(h, g) = s;
    end
end

%% ===== Usage =====
g = 1; h = 2;
fprintf('Similarity of molecule %d and %d: %.4f\n', g, h, Sim(g, h));

T = triu(Sim, 1);
[vals, order] = sort(T(:), 'descend');
[I, J] = ind2sub([nG nG], order(1:10));
fprintf('\nTop 10 most similar molecule pairs\n');
for k = 1:10
    fprintf('%5d %5d %8.4f\n', I(k), J(k), vals(k));
end
