1;

function s = graphScore(A, B, la, lb, useLabels, tol)

    if useLabels
        Z0 = double(lb(:) == la(:)');
    else
        Z0 = ones(size(B,1), size(A,1));
    end
    if norm(Z0, 'fro') == 0
        s = 0; return;
    end
    Z  = blondel_slicnost(A, B, Z0, 1, tol);
    Zs = Z * numel(Z);
    s  = (mean(max(Zs, [], 1)) + mean(max(Zs, [], 2))) / 2;
end


useLabels = true;
maxIter   = 100;
tol       = 1e-9;


edges = dlmread('MUTAG_A.txt', ',');
gi    = load('MUTAG_graph_indicator.txt');
nl    = load('MUTAG_node_labels.txt');
gl    = load('MUTAG_graph_labels.txt');
nG    = max(gi);
Aall  = sparse(edges(:,1), edges(:,2), 1, numel(gi), numel(gi));

graphs = cell(nG, 1);
labs   = cell(nG, 1);
for g = 1:nG
    nodes     = find(gi == g);
    graphs{g} = full(Aall(nodes, nodes));
    labs{g}   = nl(nodes);
end
fprintf('Učitao %d molekula.\n', nG);


selfScore = zeros(nG, 1);
for g = 1:nG
    selfScore(g) = graphScore(graphs{g}, graphs{g}, labs{g}, labs{g}, useLabels, tol);
end


Sim = eye(nG);
for g = 1:nG
    for h = g+1:nG
        s = graphScore(graphs{g}, graphs{h}, labs{g}, labs{h}, useLabels, tol);
        s = s / sqrt(selfScore(g) * selfScore(h));
        Sim(g, h) = s;
        Sim(h, g) = s;
    end
end


g = 1; h = 2;
fprintf('Sličnost molekula %d and %d: %.4f\n', g, h, Sim(g, h));

T = triu(Sim, 1);
[vals, order] = sort(T(:), 'descend');
[I, J] = ind2sub([nG nG], order(1:10));

className = {'non-mutagen', 'mutagen'};

fprintf('\nTop 10 najsličnijih parova molekula\n');
fprintf('%5s %5s %8s   %-12s %-12s %s\n', 'mol1', 'mol2', 'sim', 'klasa1', 'klasa2', 'ista klasa?');
nSame = 0;
for k = 1:10
    c1 = gl(I(k));
    c2 = gl(J(k));
    if c1 == c2
        nSame = nSame + 1;
        verdict = ['da (obje ' className{(c1 == 1) + 1} ')'];
    else
        verdict = 'ne';
    end
    fprintf('%5d %5d %8.4f   %-12s %-12s %s\n', I(k), J(k), vals(k), ...
            className{(c1 == 1) + 1}, className{(c2 == 1) + 1}, verdict);
end
fprintf('\n%d od 10 najsličnijih parova su iste klase.\n', nSame);
