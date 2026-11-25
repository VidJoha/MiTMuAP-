1;
function [S, common, d] = cosineSimilarity(B)
% Cosine similarity between the rows of a binary sparse matrix B.
%   S(i,j)      = cosine of rows i and j
%   common(i,j) = number of shared nonzero columns
%   d(i)        = number of nonzeros in row i
    N      = size(B, 1);
    d      = full(sum(B, 2));
    Bn     = spdiags(1 ./ sqrt(d), 0, N, N) * B;   % rows scaled to unit length
    S      = full(Bn * Bn');
    common = full(B * B');
end

clear; clc;

minimalan_Broj_Glasova  = 10;      % keep voters who cast at least this many votes
max_Broj_Glasaca      = 500;     % memory limit: at most this many voters (most active first)
broj_Najslicnijih_Parova      = 10;      % how many most similar pairs to print
max_Iteracija   = 100;     % max number of double steps for the iterative method
granica_Konvergencije       = 1e-9;

txtfile = 'wiki-Vote.txt';
gzfile  = 'wiki-Vote.txt.gz'; % Pripazite na velika/mala slova na Linuxu/Macu

if ~exist(txtfile, 'file')
    if exist(gzfile, 'file')
        gunzip(gzfile, '.');
    elseif exist('wiki-Vote.txt.gz', 'file')
        gunzip('wiki-Vote.txt.gz', '.');
    else
        error('Datoteka %s niti njena .gz arhiva nisu pronađene!', txtfile);
    end
end

fid = fopen(txtfile, 'r');
if fid == -1
    error('Nije moguće otvoriti datoteku %s', txtfile);
end
C = textscan(fid, '%f %f', 'CommentStyle', '#');
fclose(fid);

raw = unique([C{1}, C{2}], 'rows');

[originalni_indexi, ~, poravnati_indexi] = unique(raw(:));
E = reshape(poravnati_indexi, size(raw));
broj_Glasova = numel(originalni_indexi);
A = sparse(E(:,1), E(:,2), 1, broj_Glasova, broj_Glasova);
A = A - spdiags(diag(A), 0, broj_Glasova, broj_Glasova);
A = spones(A);



Broj_Glasova_Po_Osobi = full(sum(A, 2));
kandidati   = find(Broj_Glasova_Po_Osobi >= minimalan_Broj_Glasova);
fprintf('Glasači sa barem %d glasova: %d\n', minimalan_Broj_Glasova, numel(kandidati));
if numel(kandidati) > max_Broj_Glasaca
    [~, o] = sort(Broj_Glasova_Po_Osobi(kandidati), 'descend');
    kandidati = kandidati(o(1:max_Broj_Glasaca));
end
broj_Kandidata = numel(kandidati);


id_Kandidata   = originalni_indexi(kandidati);
gbroj_Glasova_Kandidata = Broj_Glasova_Po_Osobi(kandidati);
fid = fopen('selected_voters.txt', 'w');
fprintf(fid, '# OriginalId\tVotes\n');
fprintf(fid, '%d\t%d\n', [id_Kandidata, gbroj_Glasova_Kandidata]');
fclose(fid);



B = A(kandidati, :);
[S, common, d] = cosineSimilarity(B);

T = triu(S, 1);
[vals, order] = sort(T(:), 'descend');
[I, J] = ind2sub([broj_Kandidata broj_Kandidata], order(1:broj_Najslicnijih_Parova));

fprintf('\nTop %d najsličnijih parova glasača (kosinus)\n', broj_Najslicnijih_Parova);
fprintf('%10s %10s %8s %8s %8s %8s\n', 'id_1', 'id_2', 'kosinus', 'votes_1', 'votes_2', 'common');
for k = 1:broj_Najslicnijih_Parova
    fprintf('%10d %10d %8.3f %8d %8d %8d\n', originalni_indexi(kandidati(I(k))), originalni_indexi(kandidati(J(k))), ...
            vals(k), d(I(k)), d(J(k)), common(I(k), J(k)));
end

a = originalni_indexi(kandidati(I(1)));
b = originalni_indexi(kandidati(J(1)));
va = unique(raw(raw(:,1) == a, 2));
vb = unique(raw(raw(:,1) == b, 2));
fprintf('\nNajbolji par:\n');
fprintf('  Glasač %d (%d): %d glasova\n', a, kandidati(I(1)), numel(va));
fprintf('  Glasač %d (%d): %d glasova\n', b, kandidati(J(1)), numel(vb));
fprintf('  Zajednički glasovi: %d\n', numel(intersect(va, vb)));


As = A(kandidati, kandidati);
Z0 = ones(size(As,1),size(As,1));

Sb= blondel_slicnost(As, As,Z0, max_Iteracija, granica_Konvergencije);

Tb = triu(Sb, 1);
[valsB, orderB] = sort(Tb(:), 'descend');
[Ib, Jb] = ind2sub([broj_Kandidata broj_Kandidata], orderB(1:broj_Najslicnijih_Parova));

fprintf('\nTop %d najsličnijih parova glasača (iterativna metoda)\n', broj_Najslicnijih_Parova);
fprintf('%10s %10s %10s\n', 'id_1', 'id_2', 'sličnost');
for k = 1:broj_Najslicnijih_Parova
    fprintf('%10d %10d %10.5f\n', originalni_indexi(kandidati(Ib(k))), originalni_indexi(kandidati(Jb(k))), valsB(k));
end


[V, D] = eigs(S, 2, 'lm');
[~, eig_order] = sort(diag(D), 'descend');
V = V(:, eig_order);

[~, o] = sort(V(:, 2));
figure;
imagesc(S(o, o)); colorbar; axis square;
title(sprintf('Cosine similarity of %d voters with at least %d votes', broj_Kandidata, minimalan_Broj_Glasova));
xlabel('voters (reordered)'); ylabel('voters (reordered)');

[Vb, Db] = eigs(Sb, 2, 'lm');
[~, eo] = sort(abs(diag(Db)), 'descend');
Vb = Vb(:, eo);
[~, ob] = sort(Vb(:, 2));
figure;
imagesc(Sb(ob, ob)); colorbar; axis square;
title(sprintf('Iterative (Blondel) similarity of %d voters', broj_Kandidata));
xlabel('voters (reordered)'); ylabel('voters (reordered)');

print(figure(1), 'kosinus_metoda.png', '-dpng', '-r200');
print(figure(2), 'iterativna_metoda.png', '-dpng', '-r200');

