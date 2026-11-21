1;

function Z = blondelZ(A, B, Z0, max_Iteracija, granica_Konvergencije)
% Blondel iteration with a graf_Ideviven starting matrix Z0 (size(B,1) x size(A,1)).
% Even number of steps, normalized by the Frobenius norm each step.
    Z = Z0 / norm(Z0, 'fro');
    for it = 1:max_Iteracija
        Zold = Z;
        for step = 1:2
            Z = B*Z*A' + B'*Z*A;
            Z = Z / norm(Z, 'fro');
        end
        if norm(Z - Zold, 'fro') < granica_Konvergencije
            break;
        end
    end
end

function s = graphScore(A, B, la, lb, ratlikuj_Atome, max_Iteracija, granica_Konvergencije)
% One number describing how well G_A and G_B match (not yet normalized).
    if ratlikuj_Atome
        Z0 = double(lb(:) == la(:)');           % 1 where atom types match
    else
        Z0 = ones(size(B,1), size(A,1));
    end
    if norm(Z0, 'fro') == 0
        s = 0; return;                          % no common atom types at all
    end
    Z  = blondelZ(A, B, Z0, max_Iteracija, granica_Konvergencije);
    Zs = Z * sqrt(numel(Z));                    % remove size effect
    s  = (mean(max(Zs, [], 1)) + mean(max(Zs, [], 2))) / 2;   % symmetric best match
end

%% ===== Load MUTAG (TU Dortmund format) =====
ratlikuj_Atome = true;
max_Iteracija   = 100;
granica_Konvergencije       = 1e-9;

bridovi = dlmread('MUTAG_A.txt', ',');            % global node ids, both directions
graf_Idevi    = load('MUTAG_graph_indicator.txt');      % graph id of each node
tipovi_Atoma    = load('MUTAG_node_labels.txt');          % atom type of each node
broj_Grafova    = max(graf_Idevi);
A_Puna_Matrica  = sparse(bridovi(:,1), bridovi(:,2), 1, numel(graf_Idevi), numel(graf_Idevi));

grafovi = cell(broj_Grafova, 1);
tipovi_Atoma_Grafova   = cell(broj_Grafova, 1);
for g = 1:broj_Grafova
    atomi     = find(graf_Idevi == g);
    grafovi{g} = full(A_Puna_Matrica(atomi, atomi));
    tipovi_Atoma_Grafova{g}   = tipovi_Atoma(atomi);
end
fprintf('Loaded %d molecules.\n', broj_Grafova);

%% ===== Self scores (for normalization) =====
self_Slicnost = zeros(broj_Grafova, 1);
for g = 1:broj_Grafova
    self_Slicnost(g) = graphScore(grafovi{g}, grafovi{g}, tipovi_Atoma_Grafova{g}, tipovi_Atoma_Grafova{g}, ratlikuj_Atome, max_Iteracija, granica_Konvergencije);
end

%% ===== Pairwise similarity matrix between molecules =====
matrica_Slicnosti = eye(broj_Grafova);
for g = 1:broj_Grafova
    for h = g+1:broj_Grafova
        s = graphScore(grafovi{g}, grafovi{h}, tipovi_Atoma_Grafova{g}, tipovi_Atoma_Grafova{h}, ratlikuj_Atome, max_Iteracija, granica_Konvergencije);
        s = s / sqrt(self_Slicnost(g) * self_Slicnost(h));
        matrica_Slicnosti(g, h) = s;
        matrica_Slicnosti(h, g) = s;
    end
end

%% ===== Usage =====
g = 1; h = 2;
fprintf('Matrica slicnosti molekule %d i %d: %.4f\n', g, h, matrica_Slicnosti(g, h));

gornje_Trokut_Matrica_Slicnosti = triu(matrica_Slicnosti, 1);
[vrijednost_Najslicnijih_Vrhova, redni_Broj_Najslicnijih_Vrhova] = sort(gornje_Trokut_Matrica_Slicnosti(:), 'descend');
[I, J] = ind2sub([broj_Grafova broj_Grafova], redni_Broj_Najslicnijih_Vrhova(1:10));
fprintf('\nTop 10 most similar molecule pairs\n');
for k = 1:10
    fprintf('%5d %5d %8.4f\n', I(k), J(k), vrijednost_Najslicnijih_Vrhova(k));
end
