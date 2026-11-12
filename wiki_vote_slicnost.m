clear; clc;

minVotes  = 10;      % keep voters who cast at least this many votes
maxN      = 500;     % memory limit: at most this many voters (most active first)
topK      = 10;      % how many most similar pairs to print

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

raw = unique([C{1}, C{2}], 'rows');          % uklanja eventualne duplikate
fprintf('Unique rows in file: %d\n', size(raw, 1));

[ids, ~, idx] = unique(raw(:));              % ids(k) = original id of node k
E = reshape(idx, size(raw));
n = numel(ids);
A = sparse(E(:,1), E(:,2), 1, n, n);         % A(i,j) = 1: i voted on j
A = A - spdiags(diag(A), 0, n, n);           % remove self-votes
A = spones(A);                               % osigurava da su sve vrijednosti 0 ili 1
fprintf('Nodes: %d, edges: %d\n', n, nnz(A));


outdeg = full(sum(A, 2));
cand   = find(outdeg >= minVotes);
fprintf('Voters with at least %d votes: %d\n', minVotes, numel(cand));
if numel(cand) > maxN
    [~, o] = sort(outdeg(cand), 'descend');
    cand = cand(o(1:maxN));
    fprintf('Using the %d most active of them.\n', maxN);
end
N = numel(cand);

% Original ids of the selected voters, saved to a text file
voterIds   = ids(cand);
voterVotes = outdeg(cand);
fid = fopen('selected_voters.txt', 'w');
fprintf(fid, '# OriginalId\tVotes\n');
fprintf(fid, '%d\t%d\n', [voterIds, voterVotes]');
fclose(fid);
fprintf('Original ids of all selected voters saved to selected_voters.txt\n');


B  = A(cand, :);                             % rows = selected voters, columns = candidates
d  = full(sum(B, 2));                        % votes per voter
Bn = spdiags(1 ./ sqrt(d), 0, N, N) * B;     % rows scaled to unit length
S  = full(Bn * Bn');                         % cosine similarity, N x N
common = full(B * B');                       % number of shared candidates


T = triu(S, 1);                              % each pair once, no diagonal
[vals, order] = sort(T(:), 'descend');
[I, J] = ind2sub([N N], order(1:topK));

fprintf('\nTop %d most similar voter pairs (cosine)\n', topK);
fprintf('%10s %10s %8s %8s %8s %8s\n', 'id_1', 'id_2', 'cosine', 'votes_1', 'votes_2', 'common');
for k = 1:topK
    fprintf('%10d %10d %8.3f %8d %8d %8d\n', ids(cand(I(k))), ids(cand(J(k))), ...
            vals(k), d(I(k)), d(J(k)), common(I(k), J(k)));
end


a = ids(cand(I(1)));
b = ids(cand(J(1)));
va = unique(raw(raw(:,1) == a, 2));
vb = unique(raw(raw(:,1) == b, 2));
fprintf('\nBest pair, checked in the raw file:\n');
fprintf('  voter %d (%d): %d votes\n', a, cand(I(1)), numel(va));
fprintf('  voter %d (%d): %d votes\n', b, cand(J(1)), numel(vb));
fprintf('  shared candidates: %d\n', numel(intersect(va, vb)));
fprintf('  (search in Wiki-Vote.txt for lines starting with %d or %d)\n', a, b);


% Octave eigs ponekad ne vraća svojstvene vrijednosti sortirane silazno, pa ih ručno sortiramo
[V, D] = eigs(S, 2, 'lm');
[~, eig_order] = sort(diag(D), 'descend');
V = V(:, eig_order);

[~, o] = sort(V(:, 2));                      % order voters along 2nd eigenvector
figure;
imagesc(S(o, o)); colorbar; axis square;
title(sprintf('Cosine similarity of %d voters with at least %d votes', N, minVotes));
xlabel('voters (reordered)'); ylabel('voters (reordered)');
