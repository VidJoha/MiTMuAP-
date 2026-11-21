% Load edge list (comma-separated, 1-indexed node ids)
E = dlmread('MUTAG_A.txt', ',');

N = max(E(:));                       % number of nodes (3371 for MUTAG)
A = sparse(E(:,1), E(:,2), 1, N, N); % sparse adjacency matrix

fprintf('Nodes: %d, nonzeros: %d, symmetric: %d\n', N, nnz(A), isequal(A, A'));

figure;
spy(A);                              % shows the nonzero pattern
title('MUTAG adjacency matrix');

gi = dlmread('MUTAG_graph_indicator.txt');   % graph id of each node
k  = 1;                                      % which graph to show (1..188)

idx = find(gi == k);                         % nodes of graph k
Ak  = full(A(idx, idx));                     % its adjacency matrix

figure;
imagesc(Ak);
colormap(flipud(gray));                      % edges dark on white
axis square;
colorbar;
title(sprintf('Adjacency matrix of graph %d (%d nodes)', k, numel(idx)));
